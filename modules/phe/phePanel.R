phePanel <- tabPanel(
  title = "Phenotype processing",
  value = "PhenotypeProcessingPage",

  # 添加全局样式
  tags$head(
    tags$style(HTML("
      .phe-section {
        padding: 30px 0;
        background: #f8f9fa;
      }

      .phe-card {
        background: white;
        border-radius: 12px;
        padding: 25px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        margin-bottom: 20px;
      }

      .phe-title {
        font-size: 2.5rem;
        font-weight: 700;
        color: #2c3e50;
        margin-bottom: 1.5rem;
        text-align: center;
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

      /* 表格样式 */
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

      /* 上传组件样式 */
      .upload-section {
        margin-bottom: 30px;
        padding-bottom: 20px;
        border-bottom: 1px solid #eee;
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

      /* 选择框样式 */
      .selectize-input {
        border: 1px solid #e0e0e0;
        border-radius: 4px;
        padding: 8px 12px;
        box-shadow: none;
        transition: all 0.3s ease;
      }

      .selectize-input.focus {
        border-color: #2980b9;
        box-shadow: 0 0 0 1px rgba(41, 128, 185, 0.2);
      }

      /* 按钮样式 */
      .btn-primary {
        background-color: #2c3e50;
        border-color: #2c3e50;
        transition: all 0.3s ease;
      }

      .btn-primary:hover {
        background-color: #34495e;
        border-color: #34495e;
        transform: translateY(-2px);
        box-shadow: 0 5px 15px rgba(0,0,0,0.2);
      }

      /* 滑块样式 */
      .irs--shiny .irs-bar {
        background: #2c3e50;
      }

      .irs--shiny .irs-handle {
        border-color: #2c3e50;
      }
    "))
  ),
  div(
    class = "container-fluid",
    style = "margin-top: 60px;",
    h1(class = "phe-title", "Phenotype Processing")
  ),
  div(
    class = "phe-section",
    div(
      class = "container-fluid",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          div(
            class = "control-panel",
            div(
              id = "phe_upload_BS",
              class = "upload-section",
              fileInput("phe_upload", "Upload phenotype file (.txt or .csv)")
            ),
            bsPopover(
              id = "phe_upload_BS",
              title = "Phenotype file:",
              content = "This file is required, and the first column must be Animal ID",
              placement = "right",
              trigger = "hover",
              options = list(container = "body")
            ),
            div(
              id = "phe_colname_BS",
              class = "upload-section",
              pickerInput(
                inputId = "phe_colname",
                label = "Select the traits to be analyzed:",
                multiple = TRUE,
                choices = NULL,
                options = list(`live-search` = TRUE)
              )
            ),
            bsPopover(
              id = "phe_colname_BS",
              title = "Select the traits to be filtered:",
              content = "You can select one or more traits and we will filter each trait",
              placement = "right",
              trigger = "hover",
              options = list(container = "body")
            ),
            div(
              id = "trans0_BS",
              class = "upload-section",
              radioButtons(
                "trans0",
                "Transform 0 to NA:",
                choices = c("Yes", "No"),
                selected = "Yes",
                inline = TRUE
              )
            ),
            bsPopover(
              id = "trans0_BS",
              title = "Transform 0 values to NA:",
              content = "If Yes is selected, we will replace 0 with NA, as typically 0 values in phenotype data are meaningless and may affect subsequent analysis",
              placement = "right",
              trigger = "hover",
              options = list(container = "body")
            ),
            div(
              id = "FilterNA_BS",
              class = "upload-section",
              radioButtons(
                "FilterNA",
                "Filter NA values:",
                choices = c("Yes", "No"),
                selected = "Yes",
                inline = TRUE
              )
            ),
            bsPopover(
              id = "FilterNA_BS",
              title = "Filter NA values:",
              content = "In descriptive statistics, NA will be ignored. If Yes: When all selected traits in a record are NA, this record will be deleted.",
              placement = "right",
              trigger = "hover",
              options = list(container = "body")
            ),
            div(
              id = "sd_BS",
              class = "upload-section",
              sliderInput(
                "sd",
                "Filter range by standard deviation:",
                value = 3,
                min = 0,
                max = 10,
                step = 0.1
              )
            ),
            bsPopover(
              id = "sd_BS",
              title = "Standard deviation filter:",
              content = "Filter data outside n standard deviations from the mean. For multiple traits, out-of-range values are converted to NA.",
              placement = "right",
              trigger = "hover",
              options = list(container = "body")
            ),
            div(
              id = "plot_type_BS",
              class = "upload-section",
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
            ),
            div(
              style = "display: flex; gap: 10px; margin-top: 20px;",
              downloadButton("downloadData", "Download", class = "btn btn-primary"),
              actionButton("start", "Start!", class = "btn btn-primary")
            )
          )
        ),
        mainPanel(
          width = 9,
          div(
            class = "nav nav-pills",
            tabsetPanel(
              tabPanel(
                h4("Plot"),
                div(
                  class = "phe-card",
                  uiOutput("plot_ui")
                )
              ),
              tabPanel(
                h4("Summary"),
                div(
                  class = "phe-card",
                  DT::dataTableOutput("phe_summary")
                )
              ),
              tabPanel(
                h4("Table"),
                div(
                  class = "phe-card",
                  DT::dataTableOutput("phe_dat")
                )
              )
            )
          )
        )
      )
    )
  )
)
