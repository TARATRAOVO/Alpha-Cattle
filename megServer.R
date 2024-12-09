# 在脚本开头添加
source("vcfTools.R") # 替换为实际路径

library(shiny)
library(vcfR)
library(biomaRt)
library(dplyr)
library(plotly)

# 在文件开头添加必要的系统命令检查函数
check_system_requirements <- function() {
  # 检查bgzip是否安装
  bgzip_check <- system("which bgzip", ignore.stderr = TRUE)
  if (bgzip_check != 0) {
    stop("未找到bgzip。请安装htslib工具包。")
  }
  
  # 检查tabix是否安装
  tabix_check <- system("which tabix", ignore.stderr = TRUE)
  if (tabix_check != 0) {
    stop("未找到tabix。请安装htslib工具包。")
  }
}

# 添加处理上传文件的函数
prepare_vcf_file <- function(input_path, is_gz = FALSE) {
  # 创建临时目录
  temp_dir <- tempdir()
  
  # 构造输出文件路径
  output_gz <- file.path(temp_dir, "uploaded_vcf.vcf.gz")
  
  tryCatch({
    if (!is_gz) {
      # 如果输入文件不是gz格式，使用bgzip压缩
      system2("bgzip", 
             args = c("-c", input_path), 
             stdout = output_gz)
    } else {
      # 如果已经是gz格式，直接复制
      file.copy(input_path, output_gz, overwrite = TRUE)
    }
    
    # 创建tabix索引
    system2("tabix", 
           args = c("-p", "vcf", output_gz))
    
    return(output_gz)
  }, error = function(e) {
    stop("处理上传文件时发生错误：", e$message)
  })
}

magServer <- function(input, output, session) {
  # 读取 SNP 信息文件并确保列名正确处理
  snp_info <- read.csv("data_base/meg_snp_info/snp_info.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) # 添加 check.names = FALSE

  # 提取唯一的表型并确保不为空
  unique_traits <- unique(snp_info$`Variant Phenotype`)

  # 将表型传递给 UI，添加验证
  output$trait_choices <- renderUI({
    if (length(unique_traits) > 0) {
      selectInput("selected_trait",
        "选择分析表型:",
        choices = unique_traits,
        selected = unique_traits[1]
      )
    } else {
      div(
        style = "color: red;",
        "错误：未能加载表型数据。请检查数据文件。"
      )
    }
  })

  # 根据数据来源动态显示上传或选择文件的UI
  output$upload_ui <- renderUI({
    if (input$data_source == "上传数据") {
      fileInput("vcf_file", "选择VCF文件:", accept = c(".vcf", ".vcf.gz"))
    }
  })

  # 添加：当选择表型时立即显示相关SNP信息
  observeEvent(input$selected_trait, {
    req(input$selected_trait)

    # 获取选中表型的SNP信息
    selected_snps <- dplyr::filter(snp_info, `Variant Phenotype` == input$selected_trait)

    # 创建用于显示的数据框
    snp_info_display <- data.frame(
      表型 = selected_snps$`Variant Phenotype`,
      染色体 = selected_snps$Chr.,
      物理位置 = selected_snps$start,
      基因 = selected_snps$Gene,
      变异类型 = selected_snps$`Type of Variant`,
      变异描述 = selected_snps$variant_description,
      突变 = selected_snps$mutation
    )

    # 显示SNP信息表格
    output$selected_snp_info <- renderTable({
      snp_info_display
    })
  })

  # 读取VCF文件并显示
  observeEvent(input$load_vcf, {
    req(input$data_source, input$selected_trait)
    
    # 检查系统要求
    check_system_requirements()

    withProgress(message = "正在处理数据...", value = 0, {
      incProgress(0.1, detail = "准备加载VCF文件...")

      # 判断数据来源并处理VCF文件
      tryCatch({
        if (input$data_source == "上传数据" && !is.null(input$vcf_file)) {
          # 处理上传的文件
          showNotification("正在处理上传的VCF文件...", type = "message")
          
          # 检查文件是否为gz格式
          is_gz <- grepl("\\.gz$", input$vcf_file$name)
          
          # 准备VCF文件（压缩并建立索引）
          vcf_path <- prepare_vcf_file(input$vcf_file$datapath, is_gz)
          
          showNotification("VCF文件处理完成", type = "message")
        } else if (input$data_source == "公共数据库数据") {
          # 使用公共数据库的VCF文件
          selected_snps <- dplyr::filter(snp_info, `Variant Phenotype` == input$selected_trait)
          chr_number <- unique(selected_snps$Chr.)
          
          vcf_path <- file.path(
            "/srv/shiny-server/data_base/genomic_data",
            sprintf("samples-Chr%s-Run8-TAUIND-public_miss05_snponly_nomes_phased.vcf.gz", chr_number)
          )
          
          if (!file.exists(vcf_path)) {
            stop(sprintf("找不到染色体%s的VCF文件。", chr_number))
          }
        } else {
          stop("请选择有效的数据来源和文件。")
        }
        
        incProgress(0.4, detail = "筛选相关SNP...")
        
        # 根据选择的表型筛选相关SNP - 修改列名以匹配新格式
        selected_trait <- input$selected_trait
        relevant_snps <- dplyr::filter(snp_info, `Variant Phenotype` == selected_trait) %>%
          dplyr::transmute(
            SNP_Name = paste(Chr., start, mutation, sep = "_"),
            Chromosome = Chr.,
            Position = start,
            Mutation = mutation,
            Gene = Gene
          )

        incProgress(0.6, detail = "计算REF比例...")

        # 获取VCF文件路径
        if (input$data_source == "上传数据" && !is.null(input$vcf_file)) {
          vcf_path <- input$vcf_file$datapath
        } else {
          # 使用公共数据库的VCF文件路径
          chr_number <- unique(selected_snps$Chr.)
          vcf_path <- file.path(
            "/srv/shiny-server/data_base/genomic_data",
            sprintf("samples-Chr%s-Run8-TAUIND-public_miss05_snponly_nomes_phased.vcf.gz", chr_number)
          )
        }

        # 对每个SNP位点计算REF比例
        ref_percent <- sapply(1:nrow(relevant_snps), function(i) {
          chr <- relevant_snps$Chromosome[i]
          pos <- relevant_snps$Position[i]
          
          result <- query_ref_percent(vcf_path, chr, pos)
          
          showNotification(
            sprintf("位点 %s:%s 查询结果: %s",
                   chr, pos,
                   if(is.na(result)) "未找到" else sprintf("%.2f%%", result)),
            type = if(is.na(result)) "warning" else "message"
          )
          
          return(result)
        })

        incProgress(0.8, detail = "准备显示结果...")

        # 在计算完ref_percent后，生成饼图
        output$allele_plots <- renderUI({
          plot_output_list <- lapply(1:nrow(relevant_snps), function(i) {
            plot_id <- paste0("plot_", i)
            
            # 提取位点信息
            chr <- relevant_snps$Chromosome[i]
            pos <- relevant_snps$Position[i]
            gene <- relevant_snps$Gene[i]
            ref_pct <- ref_percent[i]
            
            if (!is.na(ref_pct)) {
              alt_pct <- 100 - ref_pct
              
              # 创建饼图数据
              plot_data <- data.frame(
                type = c("REF", "ALT"),
                value = c(ref_pct, alt_pct)
              )
              
              # 使用plotly创建交互式饼图
              output[[plot_id]] <- renderPlotly({
                plot_ly(plot_data, labels = ~type, values = ~value, type = 'pie',
                       textinfo = 'label+percent',
                       hoverinfo = 'text',
                       text = ~paste(type, sprintf("%.2f%%", value)),
                       marker = list(colors = c('#1f77b4', '#ff7f0e'))) %>%
                  layout(
                    title = list(
                      text = sprintf("位点 %s:%s (%s) 等位基因分布",
                                   chr, pos, gene),
                      y = 0.9  # 将标题向下移动
                    ),
                    margin = list(t = 100),  # 增加顶部边距
                    showlegend = TRUE
                  )
              })
              
              # 返回plotly输出容器
              div(
                style = "margin-bottom: 20px;",
                plotlyOutput(plot_id, height = "400px")
              )
            } else {
              # 如果没有数据，显示提示信息
              div(
                style = "margin-bottom: 20px; color: red;",
                sprintf("位点 %s:%s (%s) 未找到数据", chr, pos, gene)
              )
            }
          })
          
          # 将所有图表组合在一起
          do.call(tagList, plot_output_list)
        })

        # 创建结果数据框 - 修改以包含更多信息
        snp_ref_df <- data.frame(
          表型 = selected_snps$`Variant Phenotype`,
          染色体 = relevant_snps$Chromosome,
          物理位置 = relevant_snps$Position,
          基因 = relevant_snps$Gene,
          变异类型 = selected_snps$`Type of Variant`,
          突变 = relevant_snps$Mutation,
          REF等位基因比例 = sapply(ref_percent, function(x) {
            if (is.na(x)) {
              "未找到相关SNP"
            } else if (x == 0) {
              "0.00%"
            } else {
              sprintf("%.2f%%", x)
            }
          })
        )

        # 显示结果
        output$snp_ref_display <- renderTable({
          snp_ref_df
        })

        incProgress(1, detail = "完成！")

        showNotification("据处理完成。", type = "message")
      }, error = function(e) {
        showNotification(paste("错误：", e$message), type = "error")
        return(NULL)
      }, finally = {
        # 如果是上传的文件，清理临时文件
        if (input$data_source == "上传数据") {
          temp_dir <- tempdir()
          unlink(file.path(temp_dir, "uploaded_vcf.vcf.gz*"))
        }
      })
    })
  })

  # 其他已有的 server 逻辑...
}
