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
  div(class = "inlay", style = "height:70px;width:100%;background-color: white;"),  
  column(
    width = 12,
    h1("Public Data: Phenotypes")  
  ),
  div(class = "inlay", style = "height:90px;width:100%;background-color: white;"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(id = "pdata_upload_BS", style = "margin-top:30px", fileInput("pdata_upload", "Upload phenotype file (.txt or .csv):")),
      bsPopover(id = "pdata_upload_BS", 
                title = "Phenotype file:", 
                content = "This file is required, and the first column must be Animal ID",
                placement = "right", trigger = "hover", options = list(container = "body")),
      
      div(id = "pdata_colname_BS", style = "margin-top:50px", pickerInput(
        inputId = "pdata_colname",
        label = "Select the traits to be analyzed:", 
        multiple = TRUE,
        choices = NULL
      )),
      bsPopover(id = "pdata_colname_BS", 
                title = "Select the traits to be filtered:", 
                content = "You can select one or more traits, and we will filter each trait",
                placement = "right", trigger = "hover", options = list(container = "body")),
      
      div(id = "pdata_trans0_BS", style = "margin-top:50px", radioButtons("pdata_trans0", "Trans 0 to NA:", choices = c("Yes", "No"), selected = "Yes")),
      bsPopover(id = "pdata_trans0_BS", 
                title = "Verify that the 0 value is converted to NA:", 
                content = "If Yes is selected, we will replace the 0 with NA because the value of 0 in the phenotype is usually meaningless.",
                placement = "right", trigger = "hover", options = list(container = "body")),
      
      div(id = "pdata_FilterNA_BS", style = "margin-top:50px", radioButtons("pdata_FilterNA", "Filter NA:", choices = c("Yes", "No"), selected = "Yes")),
      bsPopover(id = "pdata_FilterNA_BS", 
                title = "Verify that filter NA:", 
                content = "In descriptive statistics, NA will be ignored. Yes: When all traits in a record are NA, the record will be deleted.",
                placement = "right", trigger = "hover", options = list(container = "body")),
      
      div(id = "pdata_sd_BS", style = "margin-top:50px", sliderInput("pdata_sd", "Filter range according to standard deviation", value = 3, min = 0, max = 10)),
      bsPopover(id = "pdata_sd_BS", 
                title = "Select filter criteria:", 
                content = "Filter the data based on n standard deviations from the mean. Out-of-range values are converted to NA when multiple traits are selected.",
                placement = "right", trigger = "hover", options = list(container = "body")),
      
      div(id = "downloadPData_BS", style = "margin-top:50px;text-align:center", downloadButton("downloadPData", "Download", class = "btn btn-primary")),
      div(id = "pdata_start_BS", style = "margin-top:50px;float:right", actionButton("pdata_start", "Start!"))
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