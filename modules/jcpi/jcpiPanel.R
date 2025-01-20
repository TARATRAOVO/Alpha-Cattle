jcpiPanel <- tabPanel(
  title = "Joint GCPI Calculation",
  value = "JointGCPICalculationPage",

  # 添加全局样式
  tags$head(
    tags$style(HTML("
      .jcpi-section {
        padding: 20px 50px;
        background: white;
        width: 100%;
        max-width: none;
      }

      .jcpi-card {
        background: white;
        border-radius: 12px;
        padding: 25px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        margin-bottom: 20px;
        width: 100%;
      }

      .jcpi-title {
        font-size: 2.5rem;
        font-weight: 700;
        color: #2c3e50;
        margin-bottom: 1.5rem;
        text-align: center;
      }

      .jcpi-subtitle {
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
    "))
  ),

  # 主标题
  div(
    class = "container-fluid",
    style = "margin-top: 60px;",
    h1(class = "jcpi-title", "Joint GCPI Calculation")
  ),

  # 主要内容区域
  div(
    class = "jcpi-section",
    div(
      class = "container-fluid",

      # 文件上传部分
      div(
        class = "jcpi-card",
        h3(class = "jcpi-subtitle", "Upload Files"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            div(
              class = "control-panel",
              div(
                id = "jcpi_phe_upload_BS",
                class = "upload-section",
                fileInput("jcpi_phe_upload", "Upload GCPI related phenotypic file (.txt or .csv)"),
                div(
                  class = "section-summary",
                  uiOutput("jcpi_phe_summary")
                )
              ),
              bsPopover(
                id = "jcpi_phe_upload_BS",
                title = "GCPI Related Phenotype file:",
                content = "This file is required for GBLUP to calculate GEBV for GCPI, and the first column must be Animal ID. The specific format can refer from the example table on the right of the page",
                placement = "right", trigger = "hover", options = list(container = "body")
              ),
              div(
                id = "jcpi_gen_upload_BS",
                class = "upload-section",
                fileInput("jcpi_gen_upload", "Upload genotype files (.bed, .bim, .fam)", multiple = T),
                div(
                  class = "section-summary",
                  uiOutput("jcpi_gen_summary"),
                  uiOutput("jcpi_phe_gen_summary")
                )
              ),
              bsPopover(
                id = "jcpi_gen_upload_BS",
                title = "Genotype file:",
                content = "These files are required for GBLUP to calculate GEBV for GCPI. The format is in Plink bed format and can be separated in different chromosomes and each chromosome contains .bed, .bim and .fam files",
                placement = "right", trigger = "hover", options = list(container = "body")
              )
            )
          ),
          mainPanel(
            width = 9,
            div(
              class = "control-panel",
              tags$div(
                class = "document",
                div(tags$b("Tutorial", style = "font-size:40px;")),
                br(),
                div(tags$b("Work flow", style = "font-size:25px;")),
                br(),
                tags$ol(
                  class = "tutorial",
                  tags$li(
                    "This page is only for GCPI calculation. You should upload all 3 files:",
                    tags$ul(
                      tags$li("GCPI--You should upload both 2 files at the left, the file formats are shown follow"),
                      tags$li("This page will combine the file you upload with our reference individuals, which will make the result more accurate"),
                      tags$li("Your Phenotype file should contain: Milk, Fatpct, Pr opct, SCS, Type, MS, FL, we will first calculate the GEBV for each trait with the file you upload"),
                      tags$li("Then we will use the GEBV to calculate GCPI for you"),
                      tags$li("GCPI is applicable to domestic dairy cows and progeny determination bulls with both production performance and body type identification results and genome information"),
                    )
                  ),
                  div(class = "table-container", DT::dataTableOutput("jtraits_symbol_table")),
                  tags$li(
                    "GCPI is calculated as follow, it is according to the 2016th standards of China Dairy Cattle Association:",
                    tags$ul(
                      tags$li("GCPI=20*[30*GEBV-Milk/800+15*GEBV-Fatpct/0.3+25*GEBV-Pr opct/0.12+5*GEBV-Type/5+10*GEBV-MS/5+5*GEBV-FL/5-10*(GEBV-SCS-3)/0.46]+80"),
                    )
                  ),
                  tags$li("Please upload your Phenotype file with 7 traits in the Symbol Box above in this page, which contains  Milk, Fatpct, Pr opct, SCS, Type, MS, FL.
                               See the file format introduction for details."),
                  tags$li("Submit your email address and wait the results.")
                ),
                br(),
                div(tags$b("File format instruction", style = "font-size:25px;")),
                br(),
                tags$ol(
                  class = "tutorial",
                  tags$li(
                    "Phenotype file",
                    tags$ul(
                      tags$li("Phenotype file can be comma or white space separated"),
                      tags$li("Phenotype file MUST have column names (header)"),
                      tags$li("The first column of Phenotype file MUST be animal ID"),
                      tags$li("Other columns can be other variables that may be needed in the model"),
                      tags$li("The missing value is also allowed and can be indicated as NA"),
                      tags$li("The Last 7 columns MUST be the record of 7 traits: Milk, Fatpct, Propct, SCS, Type, MS and FL"),
                      tags$li("If you do not have the record of some traits, fill the column of the traits in the phenotype file with NA"),
                      tags$li("An example of phenotype file is below:"),
                      div(class = "table-container", DT::dataTableOutput("demo_jcpi_phe"))
                    )
                  ),
                  tags$li(
                    "Genotype file:",
                    div(
                      "The genotype files were in Plink bed format,
                                                  and the standard format can be found ",
                      tags$a("here", href = "http://zzz.bwh.harvard.edu/plink/binary.shtml"),
                      ".Besides going by basic rules of Plink,
                                           the fam file is required to be prepared using the following rules:"
                    ),
                    tags$ul(
                      tags$li("No header is allowed"),
                      tags$li("The first and second columns of fam file MUST both be animal ID"),
                      tags$li("An example of fam file is below:"),
                      div(class = "table-container", DT::dataTableOutput("demo_jcpi_fam"))
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
        class = "jcpi-card",
        h3(class = "jcpi-subtitle", "Model Components in Milk Production Traits"),
        div(
          class = "container-fluid",
          fluidRow(
            column(
              width = 12,
              id = "model_box",
              style = "max-width: 1200px; margin: 0 auto;",
              mod_comp_box1(
                id = "jcpi_trt1",
                label = "Traits to be calculated :",
                title = "Traits",
                content = "You MUST choose traits the SAME as the traits in your trait file. No more and No less, which means you should choose Milk, Fatpct, Pr opct, SCS, Type, MS and FL and upload them in your Phenotype file!"
              ),
              mod_comp_box1(
                id = "jcpi_fct1",
                label = "Fixed effect (factor) :",
                title = "Fixed effects",
                content = "You can specify multiple fixed factors"
              ),
              mod_comp_box1(
                id = "jcpi_cov1",
                label = "Fixed effect (covariate) :",
                title = "Covariate",
                content = "You can specify multiple covariate varibles which are also considered fixed effects"
              ),
              mod_comp_box1(
                id = "jcpi_gen1",
                label = "Random effect (genetic) :",
                title = "Genetic random effects",
                content = "This is for the random effect that the pedigree-based or genotype-based or their combined relationship matrix is exerted on. Please note that Alpha-Cattle now only support unique random genetic effect. You may supply more varibles but only the first one is considered"
              ),
              mod_comp_box1(
                id = "jcpi_rnd1",
                label = "Random effect (other) :",
                title = "Non-genetic random effects",
                content = "This is for the non-genetic random effect, such as permenent environmental effect, dam effect, etc.",
                placement = "left"
              )
            )
          )
        )
      ),

      # 模型组件部分 2
      div(
        class = "jcpi-card",
        h3(class = "jcpi-subtitle", "Model Components in Body Shape Traits"),
        div(
          class = "container-fluid",
          fluidRow(
            column(
              width = 12,
              id = "model_box",
              style = "max-width: 1200px; margin: 0 auto;",
              mod_comp_box2(
                id = "jcpi_trt2",
                label = "Traits to be calculated :",
                title = "Traits",
                content = "You MUST choose traits the SAME as the traits in your trait file. No more and No less, which means you should choose Milk, Fatpct, Pr opct, SCS, Type, MS and FL and upload them in your Phenotype file!"
              ),
              mod_comp_box2(
                id = "jcpi_fct2",
                label = "Fixed effect (factor) :",
                title = "Fixed effects",
                content = "You can specify multiple fixed factors"
              ),
              mod_comp_box2(
                id = "jcpi_cov2",
                label = "Fixed effect (covariate) :",
                title = "Covariate",
                content = "You can specify multiple covariate varibles which are also considered fixed effects"
              ),
              mod_comp_box2(
                id = "jcpi_gen2",
                label = "Random effect (genetic) :",
                title = "Genetic random effects",
                content = "This is for the random effect that the pedigree-based or genotype-based or their combined relationship matrix is exerted on. Please note that Alpha-Cattle now only support unique random genetic effect. You may supply more varibles but only the first one is considered"
              ),
              mod_comp_box2(
                id = "jcpi_rnd2",
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
        class = "jcpi-card",
        div(
          class = "row",
          style = "display: flex; margin: 0 -15px;",
          div(
            class = "col-md-6",
            style = "padding: 0 15px;",
            h3(class = "jcpi-subtitle", "Check Model"),
            div(
              style = "min-height: 100px;",
              htmlOutput("jcpi_model_check1"),
              htmlOutput("jcpi_model_check2"),
              h4("GBLUP Model Need At Least:"),
              h5("Milk Production Traits=Mean(fix)+Herd(fix)+YearSeason(fix)+Parity(fix)+ID(gen)+pe(ID(rnd))+e"),
              h5("Body Shape Traits=Mean(fix)+Herd(fix)+YearSeason(fix)+Age(fix)+sex(fix,if both sire and dam included)+ID(gen)+pe(ID(rnd))+e")
            )
          ),
          div(
            class = "col-md-6",
            style = "padding: 0 15px;",
            h3(class = "jcpi-subtitle", "Submit Job"),
            div(
              class = "control-panel",
              style = "height: auto; min-height: 100px;",
              textInput(
                inputId = "jcpi_email",
                label = "E-mail",
                value = "Type your email address...",
                width = "100%"
              ),
              bsPopover(
                "jcpi_email",
                title = "Your e-mail address",
                content = "We do not do anything with it...apart from sending your GCPI results. :)",
                placement = "right",
                trigger = "hover",
                options = list(container = "body")
              ),
              h5("Once analysis is done, you'll receive results by email."),
              actionButton(
                inputId = "jcpi_submit_bttn",
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
