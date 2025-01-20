cpiPanel <- tabPanel(
  title = "CPI or GCPI Calculation",

  # 添加全局样式
  tags$head(
    tags$style(HTML("
      .cpi-section {
        padding: 20px 50px;
        background: white;
        width: 100%;
        max-width: none;
      }

      .cpi-card {
        background: white;
        border-radius: 12px;
        padding: 25px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        margin-bottom: 20px;
        width: 100%;
      }

      .cpi-title {
        font-size: 2.5rem;
        font-weight: 700;
        color: #2c3e50;
        margin-bottom: 1.5rem;
        text-align: center;
      }

      .cpi-subtitle {
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
    "))
  ),

  # 主标题
  div(
    class = "container-fluid",
    style = "margin-top: 60px;",
    h1(class = "cpi-title", "CPI or GCPI Calculation")
  ),

  # 主要内容区域
  div(
    class = "cpi-section",
    div(
      class = "container-fluid",

      # 文件上传部分
      div(
        class = "cpi-card",
        h3(class = "cpi-subtitle", "Upload Files"),
        sidebarLayout(
          sidebarPanel(
            width = 3,
            div(
              class = "control-panel",
              div(
                id = "typecpi",
                class = "upload-section",
                selectInput(
                  inputId = "typecpi",
                  label = "Choose the method used in Milk Production Traits EBV calculation",
                  choice = c("", "Test_Day_Model", "Average_305Days"),
                  selected = F
                ),
                div(
                  class = "section-summary",
                  uiOutput("typecpi_summary")
                )
              ),
              bsPopover(
                id = "typecpi",
                title = "Method:",
                content = "You have to choose a Method to use in Milk Production Traits EBV calculation process.",
                placement = "right",
                trigger = "hover",
                options = list(container = "body")
              ),
              div(
                id = "cpi_phe_upload_BS",
                class = "upload-section",
                fileInput("cpi_phe_upload", "Upload CPI related phenotypic file (.txt or .csv)"),
                div(
                  class = "section-summary",
                  uiOutput("cpi_phe_summary")
                )
              ),
              bsPopover(
                id = "cpi_phe_upload_BS",
                title = "CPI Related Phenotype file:",
                content = "This file is required for all CPI calculation process.and the first column must be Animal ID. The specific format can refer from the example table on the right of the page",
                placement = "right",
                trigger = "hover",
                options = list(container = "body")
              ),
              div(
                id = "cpi_ped_upload_BS",
                class = "upload-section",
                fileInput("cpi_ped_upload", "Upload CPI related pedigree file (.txt or .csv)"),
                div(
                  class = "section-summary",
                  uiOutput("cpi_ped_summary"),
                  uiOutput("cpi_phe_ped_summary")
                )
              ),
              bsPopover(
                id = "cpi_ped_upload_BS",
                title = "CPI Related Pedigree file:",
                content = "This file is required for BLUP to calculate EBV for CPI, and the first column must be Animal ID; second is sire ID; third is dam ID; fourth is order. The specific format can refer from the example table on the right of the page",
                placement = "right",
                trigger = "hover",
                options = list(container = "body")
              ),
              div(
                id = "cpi_gen_upload_BS",
                class = "upload-section",
                fileInput("cpi_gen_upload", "Upload genotype files (.bed, .bim, .fam)", multiple = T),
                div(
                  class = "section-summary",
                  uiOutput("cpi_gen_summary"),
                  uiOutput("cpi_phe_gen_summary"),
                  uiOutput("cpi_ped_gen_summary")
                )
              ),
              bsPopover(
                id = "cpi_gen_upload_BS",
                title = "Genotype file:",
                content = "These files are required for GBLUP to calculate GEBV for GCPI. The format is in Plink bed format and can be separated in different chromosomes and each chromosome contains .bed, .bim and .fam files",
                placement = "right",
                trigger = "hover",
                options = list(container = "body")
              )
            )
          ),
          mainPanel(
            width = 9,
            div(
              class = "cpi-card",
              style = "overflow-y:scroll;height:635px;",
              h3(class = "cpi-subtitle", "Tutorial"),
              tags$div(
                class = "document tutorial-section",
                style = "padding: 20px;",
                div(tags$b("Work flow", style = "font-size:25px;")),
                br(),
                tags$ol(
                  class = "tutorial",
                  tags$li(
                    "Different kinds of CPI calculations need different Upload files. If you upload Phenotype and Ped file:",
                    tags$ul(
                      tags$li("CPI1/CPI3--Your Phenotype file should contain: Milk, Fatpct, Propct, SCS, Type, MS and FL"),
                      tags$li("CPI2--Your Phenotype file should contain: Milk, Fatpct, Propct and SCS"),
                      tags$li("CPI1 is applicable to domestic cows with both production performance and body type identification results, as well as progeny determination bulls"),
                      tags$li("CPI2 is applicable to domestic cows with only production performance data and progeny determination bulls"),
                      tags$li("CPI3 is applicable to the verification of progeny test results of IMPORTED bulls")
                    )
                  ),
                  div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                  tags$li(
                    "Different kinds of CPI calculations need different Upload files. If you upload Phenotype and Genotype file:",
                    tags$ul(
                      tags$li("GCPI--You should upload both Phenotype file and Genotype file, the file formats are shown as follow"),
                      tags$li("NOTE that if you want to calculate GEBV only through the individual genome uploaded by yourself, you need to upload the phenotype file. Else if you want to combine our reference genome, please go to Joint Genomic Selection page"),
                      tags$li("Your Phenotype file should contain: Milk, Fatpct, Propct, SCS, Type, MS, FL, we will first calculate the GEBV for each trait with the file you upload"),
                      tags$li("Then we will use the GEBV to calculate GCPI for you"),
                      tags$li("GCPI is applicable to domestic dairy cows and progeny determination bulls with both production performance and body type identification results and genome information")
                    )
                  ),
                  div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                  tags$li(
                    "The symbols corresponding to various traits and the standard deviation of domestic verified bulls are shown in the table below:",
                    tags$ul(
                      DT::dataTableOutput("traits_symbol_table", width = "50%")
                    )
                  ),
                  div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                  tags$li(
                    "Different kinds of CPI are calculated as follow, they are set according to the 2016th standards of China Dairy Cattle Association:",
                    tags$ul(
                      tags$li("CPI1=20*[30*EBV-Milk/459+15*EBV-Fatpct/0.16+25*EBV-Propct/0.08+5*EBV-Type/5+10*EBV-MS/5+5*EBV-FL/5-10*(EBV-SCS-3)/0.16]"),
                      tags$li("CPI2=20*[30*EBV-Milk/459+15*EBV-Fatpct/0.16+25*EBV-Propct/0.08-10*(EBV-SCS-3)/0.16]"),
                      tags$li("CPI3=20*[30*EBV-Milk/800+10*EBV-Fatpct/0.3+20*EBV-Propct/0.12+5*EBV-Type/5+15*EBV-MS/5+10*EBV-FL/5-10*(EBV-SCS-3)/0.46]"),
                      tags$li("GCPI=20*[30*GEBV-Milk/800+15*GEBV-Fatpct/0.3+25*GEBV-Propct/0.12+5*GEBV-Type/5+10*GEBV-MS/5+5*GEBV-FL/5-10*(GEBV-SCS-3)/0.46]+80")
                    )
                  ),
                  div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                  tags$li("If you both upload Phenotype file and Pedigree file, we will send you the EBV for each trait as well as CPI1/CPI2/CPI3.
                               If your Phenotype Trait file contains both production and bodyshape Traits: Milk, Fatpct, Propct, SCS, Type, MS and FL, we will send you both CPI1 and CPI3.
                               Else if you upload Trait file only contains production Traits: Milk, Fatpct, Propct, SCS, we will send you only CPI2.
                               If you both upload Genotype file and Phenotype file, we will send you the GEBV for each trait as well as GCPI. And your Phenotype file should also contain Traits: Milk, Fatpct, Propct, SCS, Type, MS and FL.
                               Check that Phenotype file in to calculate CPI1/CPI3 contains 7 Traits (the other columns are other effects) while CPI2 Phenotype file contains 5 Traits."),
                  div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                  tags$li("Submit your email address and wait the results.")
                )
              )
            )
          )
        )
      ),

      # 模型组件部分
      div(
        class = "cpi-card",
        h3(class = "cpi-subtitle", "Model Components in Milk Production Traits"),
        div(
          class = "container-fluid",
          fluidRow(
            column(
              width = 12,
              id = "model_box",
              style = "max-width: 1200px; margin: 0 auto;",
              mod_comp_box1(
                id = "cpi_trt1",
                label = "Traits to be calculated :",
                title = "Traits",
                content = "You MUST choose traits the SAME as the traits in your trait file. No more and No less, which means you should choose Milk, Fatpct, Propct, SCS, Type, MS and FL and upload them in your Phenotype file!"
              ),
              mod_comp_box1(
                id = "cpi_fct1",
                label = "Fixed effect (factor) :",
                title = "Fixed effects",
                content = "You can specify multiple fixed factors"
              ),
              mod_comp_box1(
                id = "cpi_cov1",
                label = "Fixed effect (covariate) :",
                title = "Covariate",
                content = "You can specify multiple covariate varibles which are also considered fixed effects"
              ),
              mod_comp_box1(
                id = "cpi_gen1",
                label = "Random effect (genetic) :",
                title = "Genetic random effects",
                content = "This is for the random effect that the pedigree-based or genotype-based or their combined relationship matrix is exerted on"
              ),
              mod_comp_box1(
                id = "cpi_rnd1",
                label = "Random effect (other) :",
                title = "Non-genetic random effects",
                content = "This is for the non-genetic random effect, such as permenent environmental effect, dam effect, etc.",
                placement = "left"
              )
            )
          )
        )
      ),

      # 第二个模型组件部分
      div(
        class = "cpi-card",
        h3(class = "cpi-subtitle", "Model Components in Body Shape Traits"),
        div(
          class = "container-fluid",
          fluidRow(
            column(
              width = 12,
              id = "model_box",
              style = "max-width: 1200px; margin: 0 auto;",
              mod_comp_box2(
                id = "cpi_trt2",
                label = "Traits to be calculated :",
                title = "Traits",
                content = "You MUST choose traits the SAME as the traits in your trait file"
              ),
              mod_comp_box2(
                id = "cpi_fct2",
                label = "Fixed effect (factor) :",
                title = "Fixed effects",
                content = "You can specify multiple fixed factors"
              ),
              mod_comp_box2(
                id = "cpi_cov2",
                label = "Fixed effect (covariate) :",
                title = "Covariate",
                content = "You can specify multiple covariate varibles which are also considered fixed effects"
              ),
              mod_comp_box2(
                id = "cpi_gen2",
                label = "Random effect (genetic) :",
                title = "Genetic random effects",
                content = "This is for the random effect that the pedigree-based or genotype-based or their combined relationship matrix is exerted on"
              ),
              mod_comp_box2(
                id = "cpi_rnd2",
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
        class = "cpi-card",
        div(
          class = "row",
          style = "display: flex; margin: 0 -15px;",
          div(
            class = "col-md-6",
            style = "padding: 0 15px;",
            h3(class = "cpi-subtitle", "Check Model"),
            div(
              style = "min-height: 100px;",
              htmlOutput("cpi_model_check1"),
              htmlOutput("cpi_model_check2"),
              h4("GBLUP Model Need At Least:"),
              h5("Milk Production Traits=Mean(fix)+Herd(fix)+YearSeason(fix)+Parity(fix)+ID(gen)+pe(ID(rnd))+e"),
              h5("Body Shape Traits=Mean(fix)+Herd(fix)+YearSeason(fix)+Age(fix)+sex(fix,if both sire and dam included)+ID(gen)+pe(ID(rnd))+e")
            )
          ),
          div(
            class = "col-md-6",
            style = "padding: 0 15px;",
            h3(class = "cpi-subtitle", "Submit Job"),
            div(
              class = "control-panel",
              style = "height: auto; min-height: 100px;",
              textInput(
                inputId = "cpi_email",
                label = "E-mail",
                value = "Type your email address...",
                width = "100%"
              ),
              bsPopover(
                "cpi_email",
                title = "Your e-mail address",
                content = "We do not do anything with it...apart from sending your CPI results. :)",
                placement = "right",
                trigger = "hover",
                options = list(container = "body")
              ),
              h5("Once analysis is done, you'll receive results by email."),
              actionButton(
                inputId = "cpi_submit_bttn",
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
