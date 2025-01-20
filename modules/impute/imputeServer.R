imputeServer <- function(input, output, session) {
  # 添加日志函数
  logMessage <- function(msg, type = "INFO") {
    tryCatch(
      {
        log_dir <- "/srv/shiny-server/logs"
        log_file <- file.path(log_dir, "impute.log")

        # 检查日志目录是否存在，如果不存在则尝试创建
        if (!dir.exists(log_dir)) {
          dir.create(log_dir, recursive = TRUE, mode = "0775")
        }

        # 检查日志文件是否可写
        if (!file.exists(log_file)) {
          file.create(log_file)
          Sys.chmod(log_file, mode = "0664")
        }

        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        message <- sprintf("[%s] [%s] %s\n", timestamp, type, msg)
        cat(message, file = log_file, append = TRUE)
      },
      error = function(e) {
        # 如果写入日志失败，将错误信息打印到控制台
        warning(paste("写入日志失败:", e$message))
      }
    )
  }

  # 设置beagle路径
  BEAGLE_JAR <- "/srv/shiny-server/beagle/beagle.06Aug24.a91.jar"

  # 验证beagle文件存在
  if (!file.exists(BEAGLE_JAR)) {
    stop("Beagle JAR file not found at: ", BEAGLE_JAR)
  }

  # 添加文件管理函数
  create_task_dirs <- function(task_id) {
    tryCatch(
      {
        # 获取当前年月
        current_month <- format(Sys.Date(), "%Y%m")

        # 基础目录
        base_dir <- "/srv/shiny-server/data"

        # 检查基础目录权限
        if (!file.access(base_dir, mode = 2) == 0) {
          stop(paste("No write permission for base directory:", base_dir))
        }

        # 创建任务相关目录
        dirs <- list(
          upload_dir = file.path(base_dir, "uploads", current_month, task_id),
          result_dir = file.path(base_dir, "results", current_month, task_id),
          temp_dir = file.path(base_dir, "temp", task_id)
        )

        # 创建目录
        for (dir in dirs) {
          if (!dir.exists(dir)) {
            # 使用系统命令创建目录并设置权限
            cmd <- sprintf("mkdir -p '%s' && chmod 775 '%s'", dir, dir)
            system(cmd)

            # 验证目录是否创建成功
            if (!dir.exists(dir)) {
              stop(paste("Failed to create directory:", dir))
            }

            # 验证目录权限
            if (!file.access(dir, mode = 2) == 0) {
              stop(paste("No write permission for directory:", dir))
            }
          }
        }

        return(dirs)
      },
      error = function(e) {
        msg <- paste("Directory creation failed:", e$message)
        warning(msg)
        stop(msg)
      }
    )
  }

  # 清理过期文件
  cleanup_old_files <- function() {
    # 清理30天前的文件
    days_to_keep <- 30
    base_dir <- "/srv/shiny-server/data"

    # 获取所有月份目录
    months <- list.dirs(file.path(base_dir, "uploads"), recursive = FALSE)

    for (month_dir in months) {
      # 检查目录的修改时间
      if (difftime(Sys.time(), file.info(month_dir)$mtime, units = "days") > days_to_keep) {
        unlink(month_dir, recursive = TRUE)
        # 同时删除对应的结果目录
        result_dir <- sub("uploads", "results", month_dir)
        unlink(result_dir, recursive = TRUE)
      }
    }

    # 清理临时文件
    temp_files <- list.files(file.path(base_dir, "temp"), full.names = TRUE)
    for (temp_file in temp_files) {
      if (difftime(Sys.time(), file.info(temp_file)$mtime, units = "days") > 1) {
        unlink(temp_file, recursive = TRUE)
      }
    }
  }

  # 初始化任务ID和目录
  task <- sample(1000:9999, 1)
  dirs <- create_task_dirs(paste0(task, "_impute"))

  # 定期清理旧文件
  observe({
    invalidateLater(24 * 60 * 60 * 1000) # 每24小时执行一次
    cleanup_old_files()
  })

  output$demo_impute_fam <- DT::renderDataTable(DT::datatable(demo_gs_fam_dat,
    colnames = "",
    rownames = F
  ))

  ## gs Data upload------------------------------------------------------------------

  impute_density <- reactive({
    req(input$density)
    waiter <- Waiter$new(id = "density_summary", html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    return(input$density)
  })

  ### Genotype
  impute_gen_dat <- reactive({
    req(input$impute_gen_upload)
    waiter <- Waiter$new(id = "impute_gen_summary", html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())

    tryCatch(
      {
        file_path <- input$impute_gen_upload$datapath
        file_name <- input$impute_gen_upload$name

        # 改进文件类型检查
        if (endsWith(file_name, ".vcf.gz")) {
          ext <- "vcf.gz"
        } else {
          ext <- tools::file_ext(file_name)
        }

        if (!ext %in% c("vcf", "vcf.gz")) {
          stop("Invalid file type. Please upload .vcf or .vcf.gz file")
        }

        # 确保目标目录存在并可写
        if (!dir.exists(dirs$upload_dir) || !file.access(dirs$upload_dir, mode = 2) == 0) {
          cmd <- sprintf("mkdir -p '%s' && chmod 775 '%s'", dirs$upload_dir, dirs$upload_dir)
          system(cmd)
        }

        # 移动文件到上传目录
        new_path <- file.path(dirs$upload_dir, file_name)

        # 使用系统命令复制文件并设置权限
        cmd <- sprintf("cp '%s' '%s' && chmod 664 '%s'", file_path, new_path, new_path)
        if (system(cmd) != 0) {
          stop(paste("Failed to copy file to:", new_path))
        }

        # 删除临时文件
        unlink(file_path)

        # 验证文件
        if (!file.exists(new_path)) {
          stop(paste("File not found after copying:", new_path))
        }

        # 读取VCF文件并验证格式
        tryCatch(
          {
            if (ext == "vcf.gz") {
              con <- gzfile(new_path, "r")
            } else {
              con <- file(new_path, "r")
            }

            if (!is.null(attr(con, "status"))) {
              close(con)
              stop(paste("Cannot open file:", new_path))
            }

            header_found <- FALSE
            variant_count <- 0
            sample_count <- 0

            while (length(line <- readLines(con, n = 1)) > 0) {
              if (startsWith(line, "#CHROM")) {
                header_found <- TRUE
                fields <- strsplit(line, "\t")[[1]]
                sample_count <- length(fields) - 9 # VCF fixed fields
                break
              }
            }

            if (!header_found) {
              stop("Invalid VCF format: header not found")
            }

            # Count variants
            variant_count <- 0
            for (i in 1:100) {
              line <- readLines(con, n = 1)
              if (length(line) == 0) break
              if (!startsWith(line, "#")) variant_count <- variant_count + 1
            }

            close(con)

            res <- list()
            res$summary <- data.frame(
              Files = file_name,
              No_ind = sample_count,
              No_var = variant_count
            )
            return(res)
          },
          error = function(e) {
            if (exists("con") && isOpen(con)) {
              close(con)
            }
            stop(paste("Failed to read VCF file:", e$message))
          }
        )
      },
      error = function(e) {
        validate(paste("[ERROR:]", e$message))
      }
    )
  })

  ## data check
  output$density_summary <- renderUI({
    req(input$density)
    div(paste0("The ", input$density, " has been selected."))
  })

  ### Genotype
  output$impute_gen_summary <- renderUI({
    req(input$impute_gen_upload)
    div(paste0(
      "Found ", impute_gen_dat()$summary$No_var, "+ variants and ",
      impute_gen_dat()$summary$No_ind, " samples in the VCF file"
    ))
  })

  # 提交确认对话框
  submit_confirm <- modalDialog(
    div(
      p("确认提交任务？"),
      p("请确保："),
      tags$ul(
        tags$li("已选择填充密度"),
        tags$li("已上传基因型文件")
      )
    ),
    title = "确认提交",
    footer = tagList(
      actionButton("impute_ok", "确认", class = "btn btn-primary"),
      actionButton("cancel", "取消", class = "btn btn-default")
    )
  )

  # 密度检查对话框
  density_check <- modalDialog(
    div(
      p("请选择填充密度")
    ),
    title = "选择填充密度",
    footer = tagList(
      actionButton("cancel", "确定", class = "btn btn-warning")
    )
  )

  # Submit按钮事件处理
  observeEvent(input$impute_submit_bttn, {
    logMessage("Submit button clicked")

    # 检查密度选择
    if (is.null(input$density) || input$density == "") {
      logMessage("Density not selected", "ERROR")
      showModal(density_check)
      return()
    }

    # 检查文件上传
    if (is.null(input$impute_gen_upload)) {
      logMessage("No file uploaded", "ERROR")
      showNotification("请先上传基因型文件", type = "error")
      return()
    }

    # 显示确认对话框
    showModal(submit_confirm)
  })

  observeEvent(input$impute_ok, {
    logMessage("Confirmation dialog accepted")
    removeModal() # 立即关闭确认对话框

    tryCatch(
      {
        res <- list()
        res$dir <- dirs$result_dir
        res$upload_dir <- dirs$upload_dir
        res$temp_dir <- dirs$temp_dir
        res$gen_file <- input$impute_gen_upload$name
        res$task <- paste0(task, "_impute")
        res$density <- input$density

        logMessage(sprintf("Task info - ID: %s, File: %s", res$task, res$gen_file))

        # 保存参数
        par_file <- file.path(dirs$temp_dir, "par.RData")
        save(res, file = par_file)
        logMessage(sprintf("Parameters saved to %s", par_file))

        # 准备填补命令
        vcf_file <- file.path(dirs$upload_dir, res$gen_file)
        ref_dir <- "/srv/shiny-server/data_base/imputation"
        output_prefix <- file.path(dirs$result_dir, "imputed")

        cmd <- sprintf(
          "nohup python3 /srv/shiny-server/data_base/imputation/run_beagle_imputation.py --vcf '%s' --out '%s' --ref-dir '%s' --threads 4 > '%s/imputation.log' 2>&1 &",
          vcf_file,
          output_prefix,
          ref_dir,
          dirs$temp_dir
        )

        logMessage(sprintf("Prepared command: %s", cmd))

        # 写入执行脚本
        script_path <- file.path(dirs$temp_dir, "run_imputation.sh")
        writeLines(c("#!/bin/bash", cmd), script_path)
        system(paste("chmod +x", script_path))
        logMessage(sprintf("Script created at %s", script_path))

        # 在后台执行填补脚本
        system(script_path, wait = FALSE)
        logMessage("Imputation script executed in background")

        # 显示提交成功消息
        showNotification("填补任务已提交到后台运行", type = "message")
        showNotification(sprintf("任务ID: %s", res$task), type = "message")
        showNotification(sprintf("结果将保存在: %s", dirs$result_dir), type = "message")
        showNotification(sprintf("日志文件: %s/imputation.log", dirs$temp_dir), type = "message")

        # 更新任务ID和目录
        task <<- sample(1000:9999, 1)
        dirs <<- create_task_dirs(paste0(task, "_impute"))
      },
      error = function(e) {
        logMessage(sprintf("Error in submit handler: %s", e$message), "ERROR")
        showNotification(paste("提交失败：", e$message), type = "error")
      }
    )
  })

  # 取消按钮事件处理
  observeEvent(input$cancel, {
    logMessage("Cancel button clicked")
    removeModal()
  })
}
