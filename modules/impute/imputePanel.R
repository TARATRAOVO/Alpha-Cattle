imputePanel <- tabPanel(
  title = "Genomic Chip Imputation",
  value = "ImputationPage",

  # 添加全局样式
  tags$head(
    tags$style(HTML("
      .impute-section {
        padding: 20px 50px;
        background: white;
        width: 100%;
        max-width: none;
      }

      .impute-card {
        background: white;
        border-radius: 12px;
        padding: 25px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        margin-bottom: 20px;
        width: 100%;
      }

      .impute-title {
        font-size: 2.5rem;
        font-weight: 700;
        color: #2c3e50;
        margin-bottom: 1.5rem;
        text-align: center;
      }

      .impute-subtitle {
        font-size: 1.8rem;
        color: #2c3e50;
        margin: 5px 0;
        padding-bottom: 10px;
        border-bottom: 2px solid #eee;
      }

      .control-panel {
        background: white;
        padding: 20px;
        border-radius: 12px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        height: 635px;
        overflow-y: auto;
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

      .tutorial {
        font-size: 16px;
        line-height: 1.6;
      }

      .tutorial li {
        margin-bottom: 15px;
      }

      /* 添加表格相关样式 */
      .table-container {
        background: #fff;
        border-radius: 8px;
        padding: 15px;
        margin: 10px 0;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }

      .dataTables_wrapper {
        padding: 0;
        margin: 0;
        font-size: 14px;
      }

      .dataTable {
        width: 100% !important;
        margin: 0 !important;
        border-collapse: collapse;
      }

      .dataTable thead th {
        background-color: #f8f9fa;
        color: #2c3e50;
        font-weight: 600;
        padding: 12px 8px;
        border-bottom: 2px solid #dee2e6;
      }

      .dataTable tbody td {
        padding: 8px;
        border: 1px solid #dee2e6;
      }

      .dataTable.cell-border tbody tr:hover {
        background-color: #f8f9fa;
      }

      .dataTables_scrollBody {
        border: 1px solid #dee2e6;
        border-radius: 4px;
      }

      /* 美化上传组件样式 */
      .upload-section {
        margin-bottom: 30px;
        padding-bottom: 20px;
        border-bottom: 1px solid #eee;
      }

      .upload-section:last-child {
        border-bottom: none;
      }

      .upload-section .shiny-input-container {
        margin-bottom: 10px;
      }

      .section-summary {
        color: #666;
        font-size: 0.9em;
        margin-top: 8px;
      }

      /* 美化文件上传框 */
      .upload-section .form-group {
        margin-bottom: 0;
      }

      .upload-section .form-control {
        border: 2px dashed #ddd;
        background: #fafafa;
        height: auto;
        padding: 12px 15px;
        transition: all 0.3s ease;
      }

      .upload-section .form-control:hover {
        border-color: #2980b9;
        background: #f8f9fa;
      }

      .upload-section .btn-file {
        background: #fff;
        border: 1px solid #e0e0e0;
        color: #2c3e50;
        padding: 8px 12px;
        border-radius: 4px;
        font-weight: 500;
        transition: all 0.3s ease;
      }

      .upload-section .btn-file:hover {
        background: #f8f9fa;
        border-color: #2980b9;
        color: #2980b9;
      }

      /* 美化选择框 */
      .upload-section .selectize-input {
        border: 1px solid #e0e0e0;
        border-radius: 4px;
        padding: 8px 12px;
        box-shadow: none;
        transition: all 0.3s ease;
      }

      .upload-section .selectize-input.focus {
        border-color: #2980b9;
        box-shadow: 0 0 0 1px rgba(41, 128, 185, 0.2);
      }

      .upload-section .selectize-dropdown {
        border: 1px solid #e0e0e0;
        border-radius: 4px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.1);
      }

      .upload-section .selectize-dropdown-content {
        padding: 5px 0;
      }

      .upload-section .selectize-dropdown-content .option {
        padding: 8px 12px;
        color: #2c3e50;
      }

      .upload-section .selectize-dropdown-content .option:hover {
        background: #f8f9fa;
      }

      .upload-section .selectize-dropdown-content .active {
        background: #2980b9;
        color: white;
      }
    "))
  ),


  # 主标题
  div(
    class = "container-fluid",
    style = "margin-top: 60px;",
    h1(class = "impute-title", "Genomic Chip Imputation")
  ),

  # 主要内容区域
  div(
    class = "impute-section",
    div(
      class = "container-fluid",

      # 文件上传部分
      div(
        class = "impute-card",
        h3(class = "impute-subtitle", "Upload Files"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            div(
              class = "control-panel",
              div(
                id = "density_BS",
                class = "upload-section",
                selectInput(
                  inputId = "density",
                  label = "Choose the imputation density",
                  choice = c("", "Whole_Genome"),
                  selected = F
                ),
                div(
                  class = "section-summary",
                  uiOutput("density_summary")
                )
              ),
              bsPopover(
                id = "density_BS",
                title = "Density:",
                content = "You need to determine the target density of your chip filling process.",
                placement = "right",
                trigger = "hover",
                options = list(container = "body")
              ),
              div(
                id = "impute_gen_upload_BS",
                class = "upload-section",
                fileInput(
                  "impute_gen_upload",
                  "Upload genotype file (.vcf or .vcf.gz)",
                  multiple = F,
                  accept = c(
                    ".vcf",
                    ".vcf.gz",
                    "application/x-gzip",
                    "application/gzip"
                  )
                ),
                div(
                  class = "section-summary",
                  uiOutput("impute_gen_summary")
                )
              ),
              bsPopover(
                id = "impute_gen_upload_BS",
                title = "Genotype file:",
                content = "Please upload your genome sequencing data in VCF format (.vcf or .vcf.gz)",
                placement = "right",
                trigger = "hover",
                options = list(container = "body")
              )
            )
          ),
          mainPanel(
            width = 9,
            div(
              class = "impute-card",
              style = "overflow-y:scroll;height:635px;",
              h3(class = "impute-subtitle", "Tutorial"),
              tags$div(
                class = "document tutorial",
                div(
                  h4("Work flow", style = "color: #2c3e50; margin-bottom: 25px; font-weight: 600;"),
                  tags$ol(
                    tags$li(
                      style = "margin-bottom: 25px;",
                      strong("Determine Target Density:"),
                      tags$ul(
                        style = "margin-top: 10px;",
                        tags$li(
                          tags$span(style = "color: #2980b9; font-weight: 500;", "Whole Genome:"),
                          " Complete genomic information for imputation"
                        )
                      )
                    ),
                    tags$li(
                      style = "margin-bottom: 25px;",
                      strong("File Upload Requirements:"),
                      tags$ul(
                        style = "margin-top: 10px;",
                        tags$li(
                          "Upload genome chip sequencing data in Plink format (.bed, .bim, .fam)"
                        ),
                        tags$li(
                          "Files can be separated by chromosomes"
                        ),
                        tags$li(
                          "Each chromosome requires .bed, .bim, and .fam files"
                        )
                      )
                    ),
                    tags$li(
                      style = "margin-bottom: 25px;",
                      strong("File Format Specifications:"),
                      tags$ul(
                        style = "margin-top: 10px;",
                        tags$li(
                          "FAM file requirements:",
                          div(
                            style = "margin-top: 10px;",
                            div(
                              class = "table-container",
                              DT::dataTableOutput("demo_impute_fam")
                            )
                          )
                        ),
                        tags$li("No header row allowed in FAM files"),
                        tags$li("First and second columns must both contain animal ID")
                      )
                    )
                  )
                ),
                div(
                  style = "margin-top: 40px;",
                  h4("Reference Article", style = "color: #2c3e50; margin-bottom: 25px; font-weight: 600;"),
                  p(
                    style = "line-height: 1.6;",
                    "For information about Holstein cow chip imputation accuracy and the impact of chip density on GBLUP genome prediction accuracy, please refer to:",
                    tags$br(),
                    tags$br(),
                    tags$strong("Nazzicari N, Biscarini F. Stacked kinship CNN vs. GBLUP for genomic predictions of additive and complex continuous phenotypes. Sci Rep. 2022 Nov 18;12(1):19889."),
                    tags$br(),
                    tags$br(),
                    "Consider the conclusions of this article along with your specific needs when choosing between 150K_SNP or whole genome imputation."
                  )
                )
              )
            )
          )
        )
      ),

      # 提交部分
      div(
        class = "impute-card",
        h3(class = "impute-subtitle", "Submit Job"),
        div(
          class = "row",
          div(
            class = "col-md-6",
            textInput(
              inputId = "impute_email",
              label = "E-mail",
              value = "Type your email address...",
              width = "100%"
            ),
            bsPopover(
              "impute_email",
              title = "Your e-mail address",
              content = "We do not do anything with it...apart from sending your imputation results. :)",
              placement = "right",
              trigger = "hover",
              options = list(container = "body")
            ),
            h5("Once analysis is done, you'll receive results by email.")
          ),
          div(
            class = "col-md-6",
            style = "display: flex; align-items: flex-end;",
            actionButton(
              inputId = "impute_submit_bttn",
              label = "Submit!",
              class = "btn-info action-button",
              icon = icon("paper-plane")
            )
          )
        )
      )
    )
  )
)
