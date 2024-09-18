library(shiny)
library(vcfR)

magServer <- function(input, output, session) {
  observeEvent(input$load_vcf, {
    req(input$vcf_file, input$server_vcf)
    
    # 判断是从本地上传还是从服务器加载
    if (!is.null(input$vcf_file)) {
      vcf_data <- read_vcf_to_matrix(input$vcf_file$datapath)
    } else if (input$server_vcf != "") {
      vcf_data <- read_vcf_to_matrix(input$server_vcf)
    }
    
    # 处理vcf_data，例如显示或进一步分析
  })
}
