pdataPanel <- tabPanel(
  tags$style(HTML(".filter-option-inner-inner {color:white;}
                   .btn-file{
                   background-color: #2c3e50;
                   }
                   .jumbotron h1{
                   margin-top: 90px;
                   }
                   .bs-placeholder{
                   background-color: #2c3e50;
                  }")),
  value = "PublicDataPage",
  title = "Explore Public Data",
  icon = icon("database"),
  div(class = "inlay", style = "height:120px;width:100%;background-color: white;"),
  column(
    width = 12,
    h1("Public Data: Phenotypes")
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(
        id = "pdata_upload_BS", style = "margin-top:30px",
        pickerInput("pdata_server_file", "选择数据库文件:", choices = list.files("./data_base"), multiple = FALSE)
      ),
      bsPopover(
        id = "pdata_upload_BS",
        title = "选择数据库文件:",
        content = "从数据库中选择想要查看的数据",
        placement = "right", trigger = "hover", options = list(container = "body")
      ),
      div(id = "pdata_colname_BS", style = "margin-top:50px", pickerInput(
        inputId = "pdata_colname",
        label = "Select the traits to be analyzed:",
        multiple = TRUE,
        choices = NULL
      )),
      bsPopover(
        id = "pdata_colname_BS",
        title = "选择要筛选的性状:",
        content = "您可以选择一个或多个性状，并过滤每个性状",
        placement = "right", trigger = "hover", options = list(container = "body")
      ),
      div(id = "pdata_trans0_BS", style = "margin-top:50px", radioButtons("pdata_trans0", "Trans 0 to NA:", choices = c("Yes", "No"), selected = "Yes")),
      bsPopover(
        id = "pdata_trans0_BS",
        title = "将0转换为NA:",
        content = "选择Yes，将0替换为NA，因为表型中的0通常没有意义。",
        placement = "right", trigger = "hover", options = list(container = "body")
      ),
      div(id = "pdata_FilterNA_BS", style = "margin-top:50px", radioButtons("pdata_FilterNA", "Filter NA:", choices = c("Yes", "No"), selected = "Yes")),
      bsPopover(
        id = "pdata_FilterNA_BS",
        title = "过滤NA值:",
        content = "在描述性统计中，NA值会被忽略。选择Yes时，记录中所有性状为NA的记录将被删除。",
        placement = "right", trigger = "hover", options = list(container = "body")
      ),
      div(id = "pdata_sd_BS", style = "margin-top:50px", sliderInput("pdata_sd", "根据标准差过滤范围", value = 3, min = 0, max = 10)),
      bsPopover(
        id = "pdata_sd_BS",
        title = "选择过滤标准:",
        content = "根据与平均值相差n个标准差的范围过滤数据，超出范围的值将转换为NA。",
        placement = "right", trigger = "hover", options = list(container = "body")
      ),

      #   div(id = "downloadPData_BS", style = "margin-top:50px;text-align:center", downloadButton("downloadPData", "Download", class = "btn btn-primary")),
      #   div(id = "pdata_start_BS", style = "margin-top:50px;float:right", actionButton("pdata_start", "Start!"))
    ),
    mainPanel(
      width = 9,
      div(class = "nav nav-pills", tabsetPanel(
        tabPanel(h4("Plot"), combineWidgetsOutput("pdata_plot", height = "600px")),
        tabPanel(h4("Summary"), DT::dataTableOutput("pdata_summary")),
        tabPanel(h4("Table"), DT::dataTableOutput("pdata_dat"))
      ))
    )
  )
)
