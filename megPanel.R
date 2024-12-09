megPanel <- tabPanel(
  title = "Major-Effect Gene",
  value = "MajorEffectGenePage",
  icon = icon("biohazard"),

  # 添加全局样式
  tags$head(
    tags$style(HTML("
      .meg-section {
        padding: 30px 0;
        background: #f8f9fa;
      }

      .meg-card {
        background: white;
        border-radius: 12px;
        padding: 25px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        margin-bottom: 20px;
      }

      .meg-title {
        font-size: 2.5rem;
        font-weight: 700;
        color: #2c3e50;
        margin-bottom: 1.5rem;
        text-align: center;
      }

      .meg-subtitle {
        font-size: 1.8rem;
        color: #2c3e50;
        margin: 20px 0;
        padding-bottom: 10px;
        border-bottom: 2px solid #eee;
      }

      .control-panel {
        background: white;
        padding: 20px;
        border-radius: 12px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
      }

      .action-button {
        width: 100%;
        margin-top: 20px;
        padding: 10px;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        border-radius: 30px;
        transition: all 0.3s ease;
      }

      .action-button:hover {
        transform: translateY(-2px);
        box-shadow: 0 5px 15px rgba(0,0,0,0.2);
      }

      .result-section {
        margin-top: 20px;
      }
    "))
  ),

  # 顶部装饰条
  div(class = "inlay", style = "height:120px;width:100%;background-color: white;"),

  # 主标题
  div(
    class = "container-fluid",
    h1(class = "meg-title", "Major-Effect Gene Analysis")
  ),

  # 主要内容区域
  div(
    class = "meg-section",
    div(
      class = "container-fluid",
      sidebarLayout(
        # 侧边控制面板
        sidebarPanel(
          width = 3,
          div(
            class = "control-panel",
            selectInput("data_source", "选择数据来源:",
              choices = c("上传数据", "公共数据库数据"),
              width = "100%"
            ),

            # 上传数据部分
            uiOutput("upload_ui"),

            # 添加选择性状的功能
            uiOutput("trait_choices"),
            actionButton("load_vcf", "开始分析",
              class = "btn-primary action-button"
            )
          )
        ),

        # 主要显示区域
        mainPanel(
          width = 9,
          # SNP信息卡片
          div(
            class = "meg-card",
            h3(class = "meg-subtitle", "选中表型的SNP信息"),
            tableOutput("selected_snp_info")
          ),

          # VCF预览卡片
          div(
            class = "meg-card",
            h3(class = "meg-subtitle", "VCF文件内容预览"),
            tableOutput("vcf_display")
          ),

          # SNP REF比例卡片
          div(
            class = "meg-card",
            h3(class = "meg-subtitle", "SNP REF比例查询结果"),
            tableOutput("snp_ref_display")
          ),

          # 等位基因分布卡片
          div(
            class = "meg-card",
            h3(class = "meg-subtitle", "等位基因比例分布"),
            div(
              class = "result-section",
              uiOutput("allele_plots")
            )
          )
        )
      )
    )
  )
)
