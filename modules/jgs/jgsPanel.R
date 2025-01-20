jgsPanel <- tabPanel(
  title = "Joint Genomic Selection",
  value = "JointGenomicSelectionPage",

  # 添加全局样式
  tags$head(
    tags$style(HTML("
      .jgs-section {
        padding: 20px 50px;
        background: white;
        width: 100%;
        max-width: none;
      }

      .jgs-card {
        background: white;
        border-radius: 12px;
        padding: 25px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        margin-bottom: 20px;
        width: 100%;
      }

      .jgs-title {
        font-size: 2.5rem;
        font-weight: 700;
        color: #2c3e50;
        margin-bottom: 1.5rem;
        text-align: center;
      }

      .jgs-subtitle {
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

      .model-check {

        padding: 15px;
        border-radius: 8px;
        border-left: 4px solid #2c3e50;
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

      #model_box {
        width: 100%;
        max-width: 1200px;
        margin: 0 auto;
        padding: 0 15px;
        display: flex;
        flex-wrap: nowrap;
        justify-content: space-between;
        gap: 4px;
      }

      #model_box .col-sm-2 {
        flex: 1;
        min-width: 200px;
        padding: 0 2px;
        margin-bottom: 0;
      }

      /* 移除container-fluid的限制 */
      .container-fluid {
        max-width: none;
        padding-left: 15px;
        padding-right: 15px;
        width: 100%;
      }

      #model_box .shiny-input-container {
        width: 100%;
        margin: 0;
      }

      #model_box .multi-input-container {
        width: 100%;
        margin: 0;
      }

      #model_box .multi-input-container .non-selected,
      #model_box .multi-input-container .selected {
        width: 100%;
        margin: 0;
      }

      /* 修改sidebarPanel样式 */
      .well {
        background-color: white !important;
        border: none !important;
        box-shadow: none !important;
        padding: 0 !important;
        margin: 0 !important;
      }

      /* 添加上传组件样式 */
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

      /* 文件名显示样式 */
      .upload-section .name-display {
        color: #2c3e50;
        font-size: 0.9em;
        margin-top: 5px;
        padding: 5px 0;
      }

      /* 上传按钮样式 */
      .upload-section .btn-file {
        position: relative;
        overflow: hidden;
      }

      .upload-section .btn-file input[type=file] {
        position: absolute;
        top: 0;
        right: 0;
        min-width: 100%;
        min-height: 100%;
        font-size: 100px;
        text-align: right;
        filter: alpha(opacity=0);
        opacity: 0;
        outline: none;
        background: white;
        cursor: pointer;
        display: block;
      }

      /* 文件上传状态提示 */
      .upload-section .file-input-label {
        font-weight: 500;
        color: #2c3e50;
        margin-bottom: 8px;
      }

      .upload-section .file-input-hint {
        font-size: 0.85em;
        color: #7f8c8d;
        margin-top: 4px;
      }
    "))
  ),

  # 主标题
  div(
    class = "container-fluid",
    style = "margin-top: 60px;",
    h1(class = "jgs-title", "Joint Genomic Selection Analysis")
  ),

  # 主要内容区域
  div(
    class = "jgs-section",
    div(
      class = "container-fluid",

      # 文件上传部分
      div(
        class = "jgs-card",
        h3(class = "jgs-subtitle", "Upload Files"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            div(
              class = "control-panel",
              div(
                id = "breed",
                class = "upload-section",
                selectInput(
                  inputId = "breed", label = "Choose the breed",
                  choice = c("", "Holstein"), selected = F
                ),
                div(
                  class = "section-summary",
                  uiOutput("breed_summary")
                )
              ),
              bsPopover(
                id = "breed",
                title = "Breed:",
                content = "You have to choose a breed to combine your data with the existing data for analysis.",
                placement = "right", trigger = "hover", options = list(container = "body")
              ),
              div(
                id = "type",
                class = "upload-section",
                selectInput(
                  inputId = "type", label = "Choose the Phenotype Type",
                  choice = c(
                    "", "Milk_Production_Traits_inTDM", "Milk_Production_Traits_inAVE",
                    "Reproductive_Traits", "Body_Type_Traits"
                  ), selected = F
                ),
                div(
                  class = "section-summary",
                  uiOutput("type_summary")
                )
              ),
              bsPopover(
                id = "type",
                title = "Type:",
                content = "You have to choose a type to determine whether the Test Day Model will be used or not.",
                placement = "right", trigger = "hover", options = list(container = "body")
              ),
              div(
                id = "jgs_phe_upload_BS",
                class = "upload-section",
                fileInput("jgs_phe_upload", "Upload phenotypic file (.txt or .csv)"),
                div(
                  class = "section-summary",
                  uiOutput("jgs_phe_summary")
                )
              ),
              bsPopover(
                id = "jgs_phe_upload_BS",
                title = "Phenotype file:",
                content = "This file is required, and the first column must be Animal ID",
                placement = "right", trigger = "hover", options = list(container = "body")
              ),
              div(
                id = "jgs_ped_upload_BS",
                class = "upload-section",
                fileInput("jgs_ped_upload", "Upload pedigree file (.txt or .csv)"),
                div(
                  class = "section-summary",
                  uiOutput("jgs_ped_summary"),
                  uiOutput("jgs_phe_ped_summary")
                )
              ),
              bsPopover(
                id = "jgs_ped_upload_BS",
                title = "Pedigree file:",
                content = "This file is required for ssGBLUP. The first column is Animal ID; second is sire ID; third is dam ID; fourth is order",
                placement = "right", trigger = "hover", options = list(container = "body")
              ),
              div(
                id = "jgs_gen_upload_BS",
                class = "upload-section",
                fileInput("jgs_gen_upload", "Upload genotype files (.bed, .bim, .fam)", multiple = T),
                div(
                  class = "section-summary",
                  uiOutput("jgs_gen_summary"),
                  uiOutput("jgs_phe_gen_summary"),
                  uiOutput("jgs_ped_gen_summary")
                )
              ),
              bsPopover(
                id = "jgs_gen_upload_BS",
                title = "Genotype file:",
                content = "These files are required for GBLUP or ssGBLUP. The format is in Plink bed format and can be separated in different chromosomes and each chromosome contains .bed, .bim and .fam files",
                placement = "right", trigger = "hover", options = list(container = "body")
              ),
            )
          ),
          mainPanel(
            width = 9,
            div(
              class = "jgs-card",
              style = "overflow-y:scroll;height:635px;",
              h3(class = "jgs-subtitle", "Tutorial"),
              tags$div(
                class = "document tutorial-section",
                style = "padding: 20px;",
                div(
                  h4("Work flow", style = "color: #2c3e50; margin-bottom: 25px; font-weight: 600;"),
                  tags$ol(
                    class = "tutorial",
                    tags$li(
                      style = "margin-bottom: 25px;",
                      strong("Model Selection Based on Phenotypic Types:"),
                      tags$ul(
                        style = "margin-top: 10px;",
                        tags$li(
                          tags$span(style = "color: #2980b9; font-weight: 500;", "Milk Production Traits:"),
                          " Calculating breeding values using the Test Day Model"
                        ),
                        tags$li(
                          tags$span(style = "color: #2980b9; font-weight: 500;", "Reproductive Traits:"),
                          " Calculating breeding values using the Repeatability Model (include permanent environmental effects)"
                        ),
                        tags$li(
                          tags$span(style = "color: #2980b9; font-weight: 500;", "Body Type Traits:"),
                          " Calculating breeding values using the common Animal Model"
                        )
                      )
                    ),
                    tags$li(
                      style = "margin-bottom: 25px;",
                      strong("File Combinations & Analysis Methods:"),
                      tags$ul(
                        style = "margin-top: 10px;",
                        tags$li(
                          tags$span(style = "color: #2980b9; font-weight: 500;", "ssGBLUP Analysis:"),
                          " Upload phenotype file + pedigree file + genotype file"
                        ),
                        tags$li(
                          tags$span(style = "color: #2980b9; font-weight: 500;", "GBLUP Analysis:"),
                          " Upload phenotype file + genotype file"
                        )
                      )
                    ),
                    tags$li(
                      style = "margin-bottom: 25px;",
                      strong("Model Specification:"),
                      p(
                        style = "margin-top: 10px;",
                        "Specify model components based on column names in uploaded phenotype file. ",
                        tags$span(
                          style = "color: #e74c3c;",
                          "Note: Currently only single-trait analysis is supported. When multiple traits are specified, they will be analyzed one by one using the same model."
                        )
                      )
                    ),
                    tags$li(
                      style = "margin-bottom: 25px;",
                      strong("Model Verification:"),
                      " Review your specified model in the bottom left section"
                    ),
                    tags$li(
                      style = "margin-bottom: 25px;",
                      strong("Job Submission:"),
                      " Provide your email address and submit. Results will be sent to you once the analysis is complete"
                    )
                  )
                ),
                div(
                  style = "margin-top: 40px;",
                  h4("File Format Requirements", style = "color: #2c3e50; margin-bottom: 25px; font-weight: 600;"),
                  tags$ol(
                    class = "tutorial",
                    tags$li(
                      style = "margin-bottom: 25px;",
                      strong("Phenotype File Format:"),
                      tags$ul(
                        style = "margin-top: 10px;",
                        tags$li("File can be comma or white space separated"),
                        tags$li("Must include column headers"),
                        tags$li("First column must be animal ID"),
                        tags$li("Other columns can contain model variables"),
                        tags$li("Missing values should be marked as NA"),
                        tags$li(
                          "Example format:",
                          div(
                            style = "margin-top: 10px;",
                            div(
                              class = "table-container",
                              DT::dataTableOutput("demo_jgs_phe")
                            )
                          )
                        )
                      )
                    ),
                    tags$li(
                      style = "margin-bottom: 25px;",
                      strong("Pedigree File Format:"),
                      tags$ul(
                        style = "margin-top: 10px;",
                        tags$li("File can be comma or white space separated"),
                        tags$li("Must include column headers"),
                        tags$li(
                          "Required columns in order:",
                          tags$ol(
                            style = "margin-top: 5px;",
                            tags$li("Animal ID"),
                            tags$li("Sire ID"),
                            tags$li("Dam ID"),
                            tags$li("Order variable (e.g., birth date or numeric sorting value)")
                          )
                        ),
                        tags$li("Missing parents should be marked as NA"),
                        tags$li(
                          "Example format:",
                          div(
                            style = "margin-top: 10px;",
                            div(
                              class = "table-container",
                              DT::dataTableOutput("demo_jgs_ped")
                            )
                          )
                        )
                      )
                    ),
                    tags$li(
                      style = "margin-bottom: 25px;",
                      strong("Genotype File Format:"),
                      p(
                        style = "margin-top: 10px;",
                        "Files must be in Plink bed format. See ",
                        tags$a(
                          href = "http://zzz.bwh.harvard.edu/plink/binary.shtml",
                          "Plink documentation",
                          style = "color: #3498db; text-decoration: underline;"
                        ),
                        " for details."
                      ),
                      p(
                        style = "margin-top: 10px; color: #e74c3c;",
                        strong("Important FAM file requirements:")
                      ),
                      tags$ul(
                        style = "margin-top: 10px;",
                        tags$li("No header row allowed"),
                        tags$li("First and second columns must both contain animal ID"),
                        tags$li(
                          "Example format:",
                          div(
                            style = "margin-top: 10px;",
                            div(
                              class = "table-container",
                              DT::dataTableOutput("demo_jgs_fam")
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      ),

      # 模型组件部分
      div(
        class = "jgs-card",
        h3(class = "jgs-subtitle", "Model Components"),
        div(
          class = "container-fluid",
          fluidRow(
            column(
              width = 12,
              id = "model_box",
              style = "max-width: 1200px; margin: 0 auto;",
              mod_comp_box(
                id = "jgs_trt",
                label = "Traits to be calculated :",
                title = "Traits",
                content = "You can specify multiple traits, which will be analyzed one by one"
              ),
              mod_comp_box(
                id = "jgs_fct",
                label = "Fixed effect (factor) :",
                title = "Fixed effects",
                content = "You can specify multiple fixed factors"
              ),
              mod_comp_box(
                id = "jgs_cov",
                label = "Fixed effect (covariate) :",
                title = "Covariate",
                content = "You can specify multiple covariate varibles which are also considered fixed effects"
              ),
              mod_comp_box(
                id = "jgs_gen",
                label = "Random effect (genetic) :",
                title = "Genetic random effects",
                content = "This is for the random effect that the pedigree-based or genotype-based or their combined relationship matrix is exerted on. Please note that HEGS now only support unique random genetic effect. You may supply more varibles but only the first one is considered"
              ),
              mod_comp_box(
                id = "jgs_rnd",
                label = "Random effect (other) :",
                title = "Non-genetic random effects",
                content = "This is for the non-genetic random effect, such as permenent environmental effect, dam effect, etc.",
                placement = "left"
              )
            )
          )
        )
      ),

      # 模型检查和提交部分
      div(
        class = "jgs-card",
        div(
          class = "row",
          style = "display: flex; margin: 0 -15px;",
          div(
            class = "col-md-6",
            style = "padding: 0 15px;",
            h3(class = "jgs-subtitle", "Check Model"),
            div(
              style = "min-height: 100px;",
              htmlOutput("jgs_model_check")
            )
          ),
          div(
            class = "col-md-6",
            style = "padding: 0 15px;",
            h3(class = "jgs-subtitle", "Submit Job"),
            div(
              class = "control-panel",
              style = "height: auto; min-height: 100px;",
              textInput(
                inputId = "jgs_email", label = "E-mail",
                value = "Type your email address...", width = "100%"
              ),
              bsPopover("jgs_email",
                title = "Your e-mail address",
                content = "We do not do anything with it...apart from sending your GS results. :)",
                placement = "right", trigger = "hover",
                options = list(container = "body")
              ),
              h5("Once analysis is done, you'll receive results by email."),
              actionButton(
                inputId = "jgs_submit_bttn",
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
)
