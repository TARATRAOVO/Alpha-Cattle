megPanel <- tabPanel(
  div(class = "inlay", style = "height:70px;width:100%;background-color: white;"),
  title = "Major-Effect Gene",
  value = "MajorEffectGenePage",
  icon = icon("biohazard"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("data_source", "选择数据来源:", 
                  choices = c("上传数据", "公共数据库数据")),
      
      # 上传数据部分
      uiOutput("upload_ui"),
      
      # 公共数据库数据部分
      uiOutput("server_ui"),
      
      actionButton("load_vcf", "加载基因组数据"),
      
      # 添加选择性状的功能
      uiOutput("trait_choices")
    ),
    mainPanel(
      # 显示VCF文件内容的UI组件
      tableOutput("vcf_display"),
      
      # 新增：显示SNP名称和REF比例的表格
      h3("SNP REF比例查询结果"),
      tableOutput("snp_ref_display"),
      
      # 新增：显示数据处理进度的UI组件
      h3("数据处理进度"),
      uiOutput("processing_progress")
      
      # 其他显示分析结果的UI组件
    )
  )
)