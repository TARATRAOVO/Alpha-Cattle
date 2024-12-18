pdataServer <- function(input, output, session) {
  flag <- 0
  
  # 获取用户选择的表型列名，如果为空，则使用数据的列名
  pdata_name <- reactive({
    req(pdata_dat())
    if (is.null(input$pdata_colname)) {
      names(pdata_dat())
    } else {
      input$pdata_colname
    }
  })

  # 当用户上传数据文件时更新表型列选择器
  observeEvent(input$pdata_upload, {
    updatePickerInput(session = session, inputId = "pdata_colname", choices = pdata_name())
  })

  # 处理上传的表型数据
  pdata_dat <- reactive({
    req(input$pdata_upload)
    
    ext <- tools::file_ext(input$pdata_upload$name)
    dat <- switch(ext,
                  csv = data.table::fread(input$pdata_upload$datapath, encoding = "UTF-8"),
                  txt = data.table::fread(input$pdata_upload$datapath, encoding = "UTF-8"),
                  validate("[ERROR:] Invalid phenotype file; Please upload a .csv or .txt file")
    )
    dir.create(path_phe, showWarnings = FALSE)
    system(paste0("cp ", input$pdata_upload$datapath, " ", path_phe, "/", input$pdata_upload$name))

    # 处理 0 转换为 NA，过滤 NA 和标准差筛选逻辑
    if (!is.null(input$pdata_colname)) {
      if (input$pdata_trans0 == "Yes") {
        dat[dat == 0] <- NA
      }
      if ("ID" %in% input$pdata_colname) {
        name_filter <- input$pdata_colname[-which(input$pdata_colname == "ID")]
      } else {
        name_filter <- input$pdata_colname
      }
      test <- dat %>%
        summarise(across(
          .cols = name_filter,
          .fns = list(mean = mean, sd = sd), na.rm = TRUE
        ))
      for (i in name_filter) {
        a <- select(test, paste0(i, "_mean")) - input$pdata_sd * select(test, paste0(i, "_sd"))
        b <- select(test, paste0(i, "_mean")) + input$pdata_sd * select(test, paste0(i, "_sd"))
        dat[[i]][(select(dat, i) < a[[1]] | select(dat, i) > b[[1]])] <- NA
      }

      if (input$pdata_FilterNA == "Yes") {
        dat <- dat %>% filter(!across(name_filter, .fns = is.na))
      }
    }

    return(dat)
  })

  # 数据表渲染
  output$pdata_dat <- DT::renderDataTable({
    select(pdata_dat(), pdata_name())
  })

  # 摘要表渲染
  output$pdata_summary <- DT::renderDataTable({
    pdata_summary_num <- map_dbl(pdata_name(), function(x) sum(!is.na(pdata_dat()[[x]])))
    pdata_summary_mean <- map_dbl(pdata_name(), function(x) mean(pdata_dat()[[x]], na.rm = TRUE))
    pdata_summary_sd <- map_dbl(pdata_name(), function(x) sd(pdata_dat()[[x]], na.rm = TRUE))

    pdata <- data.frame(
      Trait = pdata_name(),
      Size = pdata_summary_num,
      Mean = pdata_summary_mean,
      SD = pdata_summary_sd
    )
    pdata
  })

  # 直方图渲染
  output$pdata_plot <- renderCombineWidgets({
    plot_dat <- pdata_dat() %>%
      select(pdata_name()) %>%
      select_if(is.numeric)
    plist <- lapply(names(plot_dat), function(x) {
      ggplotly(ggplot(plot_dat, aes_string(x)) +
                 geom_histogram(fill = "#69b3a2", bins = 30,
                                color = "#e9ecef", alpha = 0.9) +
                 theme_bw() +
                 labs(x = x, y = "Counts"))
    })

    # 使用 combineWidgets 组合多个图表
    p <- combineWidgets(list = plist, nrow = ceiling(length(plist) / 3))
    p
  })

  # 文件下载处理
  output$downloadPData <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$pdata_upload$name), "_filter.", tools::file_ext(input$pdata_upload$name))
    },
    content = function(file) {
      switch(tools::file_ext(input$pdata_upload$name),
             csv = write.csv(pdata_dat(), file, row.names = FALSE, quote = FALSE, na = ""),
             txt = write.table(pdata_dat(), file, row.names = FALSE, quote = FALSE, na = "")
      )
    }
  )
}