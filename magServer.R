# 在脚本开头添加
source("vcfTools.R")  # 替换为实际路径

library(shiny)
library(vcfR)
library(biomaRt)
library(dplyr)

magServer <- function(input, output, session) {
  
  # 读取 SNP 信息文件
  snp_info <- read.table("data_base/meg_snp_info/snp_info.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
  
  # 提取唯一的性状
  unique_traits <- unique(snp_info$Trait)
  
  # 将性状传递给 UI
  output$trait_choices <- renderUI({
    selectInput("selected_trait", "选择分析性状:", choices = unique_traits, selected = unique_traits[1])
  })
  
  # 根据数据来源动态显示上传或选择文件的UI
  output$upload_ui <- renderUI({
    if (input$data_source == "上传数据") {
      fileInput("vcf_file", "选择VCF文件:", accept = c(".vcf", ".vcf.gz"))
    }
  })
  
  output$server_ui <- renderUI({
    if (input$data_source == "公共数据库数据") {
      selectInput("server_vcf", "从服务器选择VCF文件:", 
                  choices = list.files("/srv/shiny-server/data_base/genomic_data", pattern = "\\.vcf$", full.names = TRUE))
    }
  })
  
  # 显示数据处理进度
  output$processing_progress <- renderUI({
    # 这里可以根据需要自定义进度条或日志显示
    # 例如，使用进度条组件
    h4("处理进度将在加载数据时显示。")
  })
  
  # 读取VCF文件并显示
  observeEvent(input$load_vcf, {
    req(input$data_source)
    
    # 使用进度条显示数据处理过程
    withProgress(message = '正在处理数据...', value = 0, {
      
      incProgress(0.1, detail = "加载VCF文件中...")
      
      # 判断数据来源
      if (input$data_source == "上传数据" && !is.null(input$vcf_file)) {
        vcf_path <- input$vcf_file$datapath
      } else if (input$data_source == "公共数据库数据" && input$server_vcf != "") {
        vcf_path <- input$server_vcf
      } else {
        showNotification("请选择有效的数据来源和文件。", type = "error")
        return(NULL)
      }
      
      incProgress(0.2, detail = "读取基因型矩阵...")
      
      # 读取基因型矩阵
      genotype_matrix <- read_vcf_to_matrix(vcf_path)
      
      incProgress(0.3, detail = "显示基因型矩阵...")
      
      # 显示基因型矩阵的前几行
      output$vcf_display <- renderTable({
        head(genotype_matrix)
      })
      
      incProgress(0.4, detail = "筛选相关SNP...")
      
      # 根据选择的性状筛选相关SNP
      selected_trait <- input$selected_trait
      relevant_snps <- dplyr::select(dplyr::filter(snp_info, Trait == selected_trait), SNP_Name, Chromosome, Position)
      
      incProgress(0.6, detail = "计算REF比例...")
      
      # 计算REF比例
      ref_percent <- sapply(1:nrow(relevant_snps), function(i){
        return(query_ref_percent(genotype_matrix, relevant_snps$Chromosome[i], relevant_snps$Position[i]))
      })
      
      incProgress(0.8, detail = "准备显示结果...")
      
      
      # 创建结果数据框
      # 如果ref_percent为NA，则显示"未找到相关SNP"
      snp_ref_df <- data.frame(
        SNP_Name = relevant_snps$SNP_Name,
        REF_Percentage = ifelse(is.na(ref_percent), "未找到相关SNP", ref_percent)
      )
      
      # 显示结果
      output$snp_ref_display <- renderTable({
        snp_ref_df
      })
      
      incProgress(1, detail = "完成！")
      
      showNotification("数据处理完成。", type = "message")
      
    })
    
  })
  
  # 其他已有的 server 逻辑...
}
