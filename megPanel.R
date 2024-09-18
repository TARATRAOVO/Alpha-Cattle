megPanel <- tabPanel(
  div(class = "inlay", style = "height:70px;width:100%;background-color: white;"),
  title = "Major-Effect Gene",
  value = "MajorEffectGenePage",
  icon = icon("biohazard"),
  
  # 添加选择VCF文件的功能
  sidebarLayout(
    sidebarPanel(
      fileInput("vcf_file", "选择VCF文件:", accept = c(".vcf", ".vcf.gz")),
      selectInput("server_vcf", "从服务器选择VCF文件:", 
                  choices = list.files("/srv/shiny-server/data_base/genomic_data", pattern = "\\.vcf$", full.names = TRUE)),
      actionButton("load_vcf", "从数据库中加载")
    ),
    mainPanel(
      # 这里可以添加显示VCF文件内容的UI组件
    )
  )
)