gsPanel <- tabPanel(
  title = "Regular Genomic Selection",
  value = "RegularGenomicSelectionPage",

  # 添加全局样式
  tags$head(
    tags$style(HTML("
      .gs-section {
        padding: 20px 50px;
        background: white;
        width: 100%;
        max-width: none;
      }

      .gs-card {
        background: white;
        border-radius: 12px;
        padding: 25px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        margin-bottom: 20px;
        width: 100%;
      }

      .gs-title {
        font-size: 2.5rem;
        font-weight: 700;
        color: #2c3e50;
        margin-bottom: 1.5rem;
        text-align: center;
      }

      .gs-subtitle {
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

      .table-container {
        background: #fff;
        border-radius: 8px;
        padding: 15px;
        margin: 10px 0;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
    ")),
  ),

  # 主标题
  div(
    class = "container-fluid",
    style = "margin-top: 60px;",
    h1(class = "gs-title", "Regular Genomic Selection Analysis")
  ),

  # 主要内容区域
  div(
    class = "gs-section",
    div(
      class = "container-fluid",

      # 文件上传部分
      div(
        class = "gs-card",
        h3(class = "gs-subtitle", "Upload Files"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            div(
              class = "control-panel",
              div(
                id = "typegs",
                class = "upload-section",
                selectInput(
                  inputId = "typegs",
                  label = "Choose the Phenotype Type",
                  choice = c("", "Milk_Production_Traits_inTDM", "Milk_Production_Traits_inAVE", "Reproductive_Traits", "Body_Type_Traits"),
                  selected = F
                ),
                div(
                  class = "section-summary",
                  uiOutput("typegs_summary")
                )
              ),
              bsPopover(
                id = "typegs",
                title = "Type:",
                content = "You have to choose a type to determine whether the Test Day Model will be used or not.",
                placement = "right",
                trigger = "hover",
                options = list(container = "body")
              ),
              div(
                id = "gs_phe_upload_BS",
                class = "upload-section",
                fileInput("gs_phe_upload", "Upload phenotypic file (.txt or .csv)"),
                div(
                  class = "section-summary",
                  uiOutput("gs_phe_summary")
                )
              ),
              bsPopover(
                id = "gs_phe_upload_BS",
                title = "Phenotype file:",
                content = "This file is required, and the first column must be Animal ID",
                placement = "right",
                trigger = "hover",
                options = list(container = "body")
              ),
              div(
                id = "gs_ped_upload_BS",
                class = "upload-section",
                fileInput("gs_ped_upload", "Upload pedigree file (.txt or .csv)"),
                div(
                  class = "section-summary",
                  uiOutput("gs_ped_summary"),
                  uiOutput("gs_phe_ped_summary")
                )
              ),
              bsPopover(
                id = "gs_ped_upload_BS",
                title = "Pedigree file:",
                content = "This file is required for BLUP or ssGBLUP. The first column is Animal ID; second is sire ID; third is dam ID; fourth is order",
                placement = "right",
                trigger = "hover",
                options = list(container = "body")
              ),
              div(
                id = "gs_gen_upload_BS",
                class = "upload-section",
                fileInput("gs_gen_upload", "Upload genotype files (.bed, .bim, .fam)", multiple = T),
                div(
                  class = "section-summary",
                  uiOutput("gs_gen_summary"),
                  uiOutput("gs_phe_gen_summary"),
                  uiOutput("gs_ped_gen_summary")
                )
              ),
              bsPopover(
                id = "gs_gen_upload_BS",
                title = "Genotype file:",
                content = "These files are required for GBLUP or ssGBLUP. The format is in Plink bed format and can be separated in different chromosomes and each chromosome contains .bed, .bim and .fam files",
                placement = "right",
                trigger = "hover",
                options = list(container = "body")
              ),
              div(
                id = "gs_var_upload_BS",
                class = "upload-section",
                fileInput("gs_var_upload", "Upload variance component files (.txt or .csv)"),
                div(
                  class = "section-summary",
                  uiOutput("gs_var_summary")
                )
              ),
              bsPopover(
                id = "gs_var_upload_BS",
                title = "Variance component file:",
                content = "This file is required when you want to skip the variance component estimation step and use the components you provide to predict EBV directly.",
                placement = "right",
                trigger = "hover",
                options = list(container = "body")
              )
            )
          ),
          mainPanel(
            width = 9,
            div(
              class = "control-panel",
              h3(class = "gs-subtitle", "Tutorial"),
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
                          tags$span(style = "color: #2980b9; font-weight: 500;", "BLUP Analysis:"),
                          " Upload phenotype file + pedigree file"
                        ),
                        tags$li(
                          tags$span(style = "color: #2980b9; font-weight: 500;", "GBLUP Analysis:"),
                          " Upload phenotype file + genotype file"
                        ),
                        tags$li(
                          tags$span(style = "color: #2980b9; font-weight: 500;", "Direct EBV Prediction:"),
                          " Upload variance component file + Any of the file combinations above"
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
                              DT::dataTableOutput("demo_gs_phe")
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
                              DT::dataTableOutput("demo_gs_ped")
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
                              DT::dataTableOutput("demo_gs_fam")
                            )
                          )
                        )
                      )
                    ),
                    tags$li(
                      style = "margin-bottom: 25px;",
                      strong("Variance Component File Format:"),
                      tags$ul(
                        style = "margin-top: 10px;",
                        tags$li("No header row allowed"),
                        tags$li("Each row represents variance components for each trait"),
                        tags$li("Each column represents variance components for each random effect"),
                        tags$li(
                          "Example format:",
                          div(
                            style = "margin-top: 10px;",
                            div(
                              class = "table-container",
                              DT::dataTableOutput("demo_gs_var")
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
        class = "gs-card",
        h3(class = "gs-subtitle", "Model Components"),
        div(
          class = "container-fluid",
          fluidRow(
            column(
              width = 12,
              id = "model_box",
              style = "max-width: 1200px; margin: 0 auto;",
              mod_comp_box(
                id = "gs_trt",
                label = "Traits to be calculated :",
                title = "Traits",
                content = "You can specify multiple traits, which will be analyzed one by one"
              ),
              mod_comp_box(
                id = "gs_fct",
                label = "Fixed effect (factor) :",
                title = "Fixed effects",
                content = "You can specify multiple fixed factors"
              ),
              mod_comp_box(
                id = "gs_cov",
                label = "Fixed effect (covariate) :",
                title = "Covariate",
                content = "You can specify multiple covariate varibles which are also considered fixed effects"
              ),
              mod_comp_box(
                id = "gs_gen",
                label = "Random effect (genetic) :",
                title = "Genetic random effects",
                content = "This is for the random effect that the pedigree-based or genotype-based or their combined relationship matrix is exerted on. Please note that HEGS now only support unique random genetic effect. You may supply more varibles but only the first one is considered"
              ),
              mod_comp_box(
                id = "gs_rnd",
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
        class = "gs-card",
        div(
          class = "row",
          style = "display: flex; margin: 0 -15px;",
          div(
            class = "col-md-6",
            style = "padding: 0 15px;",
            h3(class = "gs-subtitle", "Check Model"),
            div(
              style = "min-height: 100px;",
              htmlOutput("gs_model_check")
            )
          ),
          div(
            class = "col-md-6",
            style = "padding: 0 15px;",
            h3(class = "gs-subtitle", "Submit Job"),
            div(
              class = "control-panel",
              style = "height: auto; min-height: 100px;",
              textInput(
                inputId = "gs_email",
                label = "E-mail",
                value = "Type your email address...",
                width = "100%"
              ),
              bsPopover(
                "gs_email",
                title = "Your e-mail address",
                content = "We do not do anything with it...apart from sending your GS results. :)",
                placement = "right",
                trigger = "hover",
                options = list(container = "body")
              ),
              h5("Once analysis is done, you'll receive results by email."),
              actionButton(
                inputId = "gs_submit_bttn",
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
