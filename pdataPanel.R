pdataPanel <- tabPanel(
  # 设置页面样式
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
  value = "PublicDataPage",   # 设置页面值，用于导航
  title = "Explore Public Data",  # 设置页面标题
  icon = icon("database"),  # 页面图标，数据库图标
  div(class = "inlay", style = "height:70px;width:100%;background-color: white;"),  # 页面顶部间距
  column(
    width = 12,
    h1("Phenotype processing")  # 页面标题
  ),
  div(class = "inlay", style = "height:90px;width:100%;background-color: white;"),  # 页面底部间距
  
  # 主页面布局，侧边栏和主面板
  sidebarLayout(
    sidebarPanel(
      width = 3,
      # 列名选择器，用于选择需要分析的性状
      div(id = "phe_colname_BS", style = "margin-top:50px", pickerInput(
        inputId = "phe_colname",
        label = "Select the traits to be analyzed:", 
        multiple = TRUE,
        choices = c("Trait1", "Trait2", "Trait3"),
      )),
      bsPopover(id = "phe_colname_BS", 
                title = "Select the traits to be filtered:", 
                content = "您可以选择一个或多个性状，系统将根据选择的性状进行过滤",
                placement = "right", trigger = "hover", options = list(container = "body")),
      
      # 是否将0转换为NA的选项
      div(id = "trans0_BS", style = "margin-top:50px", radioButtons("trans0", "Trans 0 to NA:", choices = c("Yes", "No"), selected = "Yes")),
      bsPopover(id = "trans0_BS", 
                title = "Verify that the 0 value is converted to NA:", 
                content = "如果选择 Yes，0 值将被替换为 NA，因为通常表型中的 0 值是无意义的，会影响后续分析",
                placement = "right", trigger = "hover", options = list(container = "body")),
      
      # 是否过滤掉NA值的选项
      div(id = "FilterNA_BS", style = "margin-top:50px", radioButtons("FilterNA", "Filter NA:", choices = c("Yes", "No"), selected = "Yes")),
      bsPopover(id = "FilterNA_BS", 
                title = "Verify that filter NA:", 
                content = "在描述性统计中，NA 会被忽略。选择 Yes 时，当所有选择的性状在某条记录中均为 NA 时，该记录将被删除。",
                placement = "right", trigger = "hover", options = list(container = "body")),
      
      # 滑块输入，设置标准差范围，用于筛选数据
      div(id = "sd_BS", style = "margin-top:50px", sliderInput("sd", "Filter range according to standard deviation", value = 3, min = 0, max = 10)),
      bsPopover(id = "sd_BS", 
                title = "Select filter criteria:", 
                content = "根据均值的 n 个标准差范围筛选数据。当选择多个性状时，超出范围的值将被转换为 NA。",
                placement = "right", trigger = "hover", options = list(container = "body")),
      
      # 下载按钮，用户可以下载筛选后的数据
      div(id = "download_BS", style = "margin-top:50px;text-align:center", downloadButton("downloadData", "Download", class = "btn btn-primary")),
      
      # 开始按钮，用于触发数据分析
      div(id = "start_BS", style = "margin-top:50px;float:right", actionButton("start", "Start!"))
    ),
    
    # 主面板，显示绘图、摘要和表格
    mainPanel(
      width = 9,
      
      # 设置选项卡面板，包含绘图、摘要和表格选项
      div(class = "nav nav-pills", tabsetPanel(
        
        # 绘图选项卡，显示性状的直方图
        tabPanel(h4("Plot"), combineWidgetsOutput("plot", height = "600px")),
        
        # 摘要选项卡，显示性状的统计摘要信息
        tabPanel(h4("Summary"), DT::dataTableOutput("phe_summary")),
        
        # 表格选项卡，显示数据表
        tabPanel(h4("Table"), DT::dataTableOutput("phe_dat"))
      ))
    )
  )
)   