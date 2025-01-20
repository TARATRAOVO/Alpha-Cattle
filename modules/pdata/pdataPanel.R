pdataPanel <- tabPanel(
  tags$head(
    tags$style(HTML("
      .pdata-section {
        padding: 30px 0;
        background: #f8f9fa;
      }

      .pdata-card {
        background: white;
        border-radius: 12px;
        padding: 25px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        margin-bottom: 20px;
      }

      .pdata-title {
        font-size: 2.5rem;
        font-weight: 700;
        color: #2c3e50;
        margin-bottom: 1.5rem;
        text-align: center;
      }

      .pdata-subtitle {
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

      .filter-option-inner-inner {color: white;}
      .btn-file {
        background-color: #2c3e50;
        border-color: #2c3e50;
      }
      .bs-placeholder {
        background-color: #2c3e50;
      }
      .nav-pills > li.active > a,
      .nav-pills > li.active > a:focus,
      .nav-pills > li.active > a:hover {
        background-color: #2c3e50;
      }
      .btn-primary {
        background-color: #2c3e50;
        border-color: #2c3e50;
      }
      .btn-primary:hover {
        background-color: #34495e;
        border-color: #34495e;
      }
      .well {
        background-color: #f8f9fa;
        border: none;
        box-shadow: 0 0 10px rgba(0,0,0,0.1);
        border-radius: 8px;
      }
      .control-label {
        color: #2c3e50;
        font-weight: 600;
      }
      .bootstrap-select .dropdown-menu {
        border-radius: 8px;
      }
      .slider-handle {
        background-color: #2c3e50;
        background-image: none;
      }
      .slider-selection {
        background-color: #95a5a6;
        background-image: none;
      }
    ")),
  ),
  value = "PublicDataPage",
  title = "Explore Public Data",
  icon = icon("database"),
  div(
    class = "container-fluid",
    style = "margin-top: 60px;",
    h1(class = "pdata-title", "Public Data: Phenotypes")
  ),
  div(
    class = "pdata-section",
    div(
      class = "container-fluid",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          div(
            class = "control-panel",
            div(
              id = "pdata_upload_BS",
              style = "margin-bottom: 30px;",
              pickerInput(
                "pdata_server_file",
                "Select Database File:",
                choices = list.files("./data_base/pheno_data"),
                multiple = FALSE,
                options = list(`live-search` = TRUE)
              )
            ),
            bsPopover(
              id = "pdata_upload_BS",
              title = "Select Database File:",
              content = "Select the data you want to view from the database",
              placement = "right",
              trigger = "hover"
            ),
            div(
              id = "pdata_colname_BS",
              style = "margin-bottom: 30px;",
              pickerInput(
                inputId = "pdata_colname",
                label = "Select Traits to Analyze:",
                multiple = TRUE,
                choices = NULL,
                options = list(`live-search` = TRUE)
              )
            ),
            bsPopover(
              id = "pdata_colname_BS",
              title = "Select Traits to Filter:",
              content = "You can select one or more traits, and we will filter each trait",
              placement = "right",
              trigger = "hover"
            ),
            div(
              id = "pdata_trans0_BS",
              style = "margin-bottom: 30px;",
              radioButtons(
                "pdata_trans0",
                "Transform 0 to NA:",
                choices = c("Yes", "No"),
                selected = "Yes",
                inline = TRUE
              )
            ),
            div(
              id = "pdata_FilterNA_BS",
              style = "margin-bottom: 30px;",
              radioButtons(
                "pdata_FilterNA",
                "Filter NA:",
                choices = c("Yes", "No"),
                selected = "Yes",
                inline = TRUE
              )
            ),
            div(
              id = "pdata_sd_BS",
              style = "margin-bottom: 30px;",
              sliderInput(
                "pdata_sd",
                "Filter Range by Standard Deviation",
                value = 3,
                min = 0,
                max = 10,
                step = 0.1
              )
            ),
            div(
              id = "plot_type_BS",
              style = "margin-bottom: 30px;",
              selectInput(
                "plot_type",
                "Select Plot Type:",
                choices = c(
                  "Histogram" = "histogram",
                  "Box Plot" = "boxplot",
                  "Scatter Plot" = "scatter",
                  "Correlation Heatmap" = "correlation",
                  "Density Plot" = "density",
                  "Violin Plot" = "violin"
                ),
                selected = "histogram"
              )
            ),
            conditionalPanel(
              condition = "input.plot_type == 'scatter'",
              div(
                style = "margin-bottom: 30px;",
                pickerInput(
                  "scatter_x",
                  "X-axis Variable:",
                  choices = NULL,
                  options = list(`live-search` = TRUE)
                ),
                pickerInput(
                  "scatter_y",
                  "Y-axis Variable:",
                  choices = NULL,
                  options = list(`live-search` = TRUE)
                )
              )
            )
          )
        ),
        mainPanel(
          width = 9,
          div(
            class = "nav nav-pills",
            style = "margin-bottom: 20px;",
            tabsetPanel(
              tabPanel(
                h4("Plot"),
                div(
                  class = "pdata-card",
                  uiOutput("plot_ui")
                )
              ),
              tabPanel(
                h4("Summary"),
                div(
                  class = "pdata-card",
                  DT::dataTableOutput("pdata_summary")
                )
              ),
              tabPanel(
                h4("Table"),
                div(
                  class = "pdata-card",
                  DT::dataTableOutput("pdata_dat")
                )
              )
            )
          )
        )
      )
    )
  )
)
