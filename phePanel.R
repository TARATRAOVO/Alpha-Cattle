 phePanel<-tabPanel(
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
   title = "Phenotype processing",
   value = "PhenotypeProcessingPage",
   div(class = "inlay", style = "height:70px;width:100%;background-color: white;"),
   column(
     width=12,
     h1("Phenotype processing")
   ),div(class = "inlay", style = "height:90px;width:100%;background-color: white;"),
 sidebarLayout(
    sidebarPanel(
         width=3,
          div(id="phe_upload_BS",style="margin-top:30px",fileInput("phe_upload","Upload phenotype file(.txt or .csv):")),
          bsPopover(id="phe_upload_BS", 
                        title="Phenotype file:", 
                        content="This file is required, and the first column must be Animal ID",
                        placement="right",trigger="hover", options = list(container = "body")),
         
          div(id="phe_colname_BS",style="margin-top:50px",pickerInput(
               inputId = "phe_colname",
               label = "Select the traits to be analyzed:", 
               multiple = TRUE,
               choices = colnames(demo_gs_phe_dat)
          )),
          bsPopover(id="phe_colname_BS", 
                    title="Select the traits to be filtered :", 
                    content="You can select one or more traits and we will filter each trait",
                    placement="right",trigger="hover", options = list(container = "body")),
         ###��0ת��ΪNA������ʱ�Ƿ����NA
         div(id="trans0_BS", style="margin-top:50px",radioButtons("trans0","Trans 0 to NA:",choices = c("Yes","No"),selected = "Yes")),
         bsPopover(id="trans0_BS", 
                   title="Verify that the 0 value is converted to NA :", 
                   content="If Yes is selected, we will replace the 0  with NA, because normally the value of 0 in the phenotype is meaningless and will affect the subsequent analysis",
                   placement="right",trigger="hover", options = list(container = "body")),

         div(id="FilterNA_BS", style="margin-top:50px",radioButtons("FilterNA","Filter NA:",choices = c("Yes","No"),selected = "Yes")),
         bsPopover(id="FilterNA_BS", 
                   title="Verify that filter NA :", 
                   content="In descriptive statistics, NA will be ignored. Yes: When all traits selected in a recode are NA, this recode will be deleted.",
                   placement="right",trigger="hover", options = list(container = "body")),
         div(id="sd_BS",style="margin-top:50px",sliderInput("sd","Filter range according to standard deviation",value = 3,min=0,max=10)),
         bsPopover(id="sd_BS", 
                   title="Select filter criteria :", 
                   content="Filter the data out of n standard deviation of the mean.When multiple traits are selected, out-of-range values are converted to NA.",
                   placement="right",trigger="hover", options = list(container = "body")),
         div(id="download_BS",style="margin-top:50px;text-align:center",downloadButton("downloadData", "Download",class = "btn btn-primary")),
         div(id="start_BS",style="margin-top:50px;float:right",actionButton("start", "Start!"))
         
    )
         ,
         mainPanel(
            width=9,
             
              div(class="nav nav-pills",tabsetPanel(
                   
                   tabPanel(h4("Plot"), combineWidgetsOutput("plot",height = "600px")),
                   
                   tabPanel(h4("Summary"), DT::dataTableOutput("phe_summary")),
                   
                   tabPanel(h4("Table"),  DT::dataTableOutput("phe_dat"))
              )
                            
                )
         
                     
             )
        )
   )

 