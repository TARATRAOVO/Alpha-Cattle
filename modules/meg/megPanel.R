megPanel <- tabPanel(
  title = "Major-Effect Gene",
  value = "MajorEffectGenePage",
  icon = icon("biohazard"),

  # Add global styles
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

  # Main title
  div(
    class = "container-fluid",
    style = "margin-top: 60px;",
    h1(class = "meg-title", "Major-Effect Gene Analysis")
  ),

  # Main content area
  div(
    class = "meg-section",
    div(
      class = "container-fluid",
      sidebarLayout(
        # Sidebar control panel
        sidebarPanel(
          width = 3,
          div(
            class = "control-panel",
            selectInput("data_source", "Select Data Source:",
              choices = c("Upload Data", "Public Database Data"),
              width = "100%"
            ),

            # Upload data section
            uiOutput("upload_ui"),

            # Add trait selection functionality
            uiOutput("trait_choices"),
            actionButton("load_vcf", "Start Analysis",
              class = "btn-primary action-button"
            )
          )
        ),

        # Main display area
        mainPanel(
          width = 9,
          # SNP information card
          div(
            class = "meg-card",
            h3(class = "meg-subtitle", "SNP Information for Selected Phenotype"),
            tableOutput("selected_snp_info")
          ),

          # VCF preview card
          div(
            class = "meg-card",
            h3(class = "meg-subtitle", "VCF File Content Preview"),
            tableOutput("vcf_display")
          ),

          # SNP REF ratio card
          div(
            class = "meg-card",
            h3(class = "meg-subtitle", "SNP REF Ratio Query Results"),
            tableOutput("snp_ref_display")
          ),

          # Allele distribution card
          div(
            class = "meg-card",
            h3(class = "meg-subtitle", "Allele Ratio Distribution"),
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
