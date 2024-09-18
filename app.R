library(shinydashboard)
library(shiny)
library(rintrojs)
library(shinyBS)
library(plotly)
library(waiter)
#library(ggiraph)
library(shinycssloaders)
library(hchinamap)
library(ggplot2)
library(ggpubr)
library(shinyjs)
library(stringr)
library(shinyWidgets)
library(tidyverse)
library(shinyalert) 
# library(leaflet)
# library(leafletCN)
# library(vroom)
# library(tools)
library(shinydashboardPlus)
# library(shinythemes)
library(dygraphs)
library(manipulateWidget)
library(shinyLP)


options(shiny.maxRequestSize = 10 * 1024^3)
dir.create("./temp",showWarnings = F)
task<<-sample(1000:9999,1)
path_gs <<- paste0("temp/",task,"_gs")
path_phe <<- paste0("temp/",task,"_phe")
path_impute <<- paste0("temp/",task,"_impute")
path_jgs <<- paste0("temp/",task,"_jgs")
path_cpi <<- paste0("temp/",task,"_cpi")
path_jcpi <<- paste0("temp/",task,"_jcpi")


demo_gs_phe_dat<<-data.frame(
  ID=paste0("IND",1:6),
  Mean=rep(1,6),
  Herd=c(rep(1,3),rep(2,3)),
  YearSeason=rep(1:3,2),
  Sex=sample(c("F","M"),6,replace = T),
  Parity=sample(1:3,6,replace = T),
  Trait1=abs(round(rnorm(6,0,1),2)),
  Trait2=abs(round(rnorm(6,10,100),2))
)

demo_gs_ped_dat<<-data.frame(ID=paste0("IND",1:6),
                            Sire=c("<NA>","<NA>","IND1","IND1","IND3","IND3"),
                            Dam=c("<NA>","<NA>","<NA>","IND2","IND4","IND4"),
                            Order=1:6
)

demo_gs_fam_dat<<-data.frame(ID1=paste0("IND",1:6),
                            ID2=paste0("IND",1:6),
                            V3=0,V4=0,V5=0,V6=-9)

demo_gs_var_dat<<-data.frame(V1=c(10,20),V2=c(5,3),V3=c(40,80))

traits_symbol_table_dat<<-data.frame(Trait=c("产奶量305天","乳脂率","乳蛋白率","体型总分","泌乳系统","肢蹄","体细胞评分"),
                             Symbol=c("Milk","Fatpct","Propct","Type","MS","FL","SCS"),
                             Type=c("Milk Production Trait","Milk Production Trait","Milk Production Trait","Body Shape Trait","Body Shape Trait","Body Shape Trait","Milk Production Trait"))

demo_trait_type1_dat<<-data.frame(
                             ID=paste0("IND",1:6),
                             Mean=rep(1,6),
                             Herd=c(rep(1,3),rep(2,3)),
                             YearSeason=rep(1:3,2),
                             Parity=sample(1:3,6,replace = T),
                             Age=c(27,15,18,25,35,32),
                             Milk=c(670,690,720,750,780,800),
                             Fatpct=c(0.24,0.27,0.30,0.33,0.36,0.39),
                             Propct=c(2.6,2.8,3.0,3.2,3.4,3.6),
                             Type=c(5,6,7,8,9,10),
                             MS=c(5,6,7,8,9,10),
                             FL=c(5,6,7,8,9,10),
                             SCS=c(5,6,7,8,9,10))

demo_trait_type2_dat<<-data.frame(
                             ID=paste0("IND",1:6),
                             Mean=rep(1,6),
                             Herd=c(rep(1,3),rep(2,3)),
                             YearSeason=rep(1:3,2),
                             Parity=sample(1:3,6,replace = T),
                             Age=c(27,15,18,25,35,32),
                             Milk=c(670,690,720,750,780,800),
                             Fatpct=c(0.24,0.27,0.30,0.33,0.36,0.39),
                             Propct=c(2.6,2.8,3.0,3.2,3.4,3.6),
                             SCS=c(5,6,7,8,9,10))



mod_comp_box<<-function(id,label,title,content, placement="right",choices=colnames(demo_gs_phe_dat),single=FALSE){
  if(!single){
    box<-column(
      width=2,
      
      div(id=paste0(id,"_BS"),
         

            multiInput(
              inputId = id,
              label = span(id=paste0(id,"_BS"),label),
              choices = choices,
              options = list(enable_search = TRUE, non_selected_header = "Choose between:",
                             selected_header = "You have selected:"),
              width = "100%"
            )
          
      ),
      
      bsPopover(id=paste0(id,"_BS"), 
                title=title, 
                content=content,
                placement=placement,trigger="hover", options = list(container = "body"))
      
    )
  }else{
    box<-column(
      width=2,
      
      div(id=paste0(id,"_BS"),
          box(
            width=12,
            title = "",
            status = "danger",
            solidHeader = F,
            multiInput(
              inputId = id,
              label = span(id=paste0(id,"_BS"),label),
              choices = choices,
              options = list(enable_search = TRUE, non_selected_header = "Choose between:",
                             selected_header = "You have selected:",limit=1),
              width = "100%"
            )
          )
      ),
      
      bsPopover(id=paste0(id,"_BS"), 
                title=title, 
                content=content,
                placement=placement,trigger="hover", options = list(container = "body"))
      
    )
  }
  return(box)
}

mod_comp_box1<<-function(id,label,title,content, placement="right",choices=colnames(demo_trait_type1_dat)[c(1,2,3,4,5,7,8,9,13)],single=FALSE){
  if(!single){
    box<-column(
      width=2,
      
      div(id=paste0(id,"_BS"),
         

            multiInput(
              inputId = id,
              label = span(id=paste0(id,"_BS"),label),
              choices = choices,
              options = list(enable_search = TRUE, non_selected_header = "Choose between:",
                             selected_header = "You have selected:"),
              width = "100%"
            )
          
      ),
      
      bsPopover(id=paste0(id,"_BS"), 
                title=title, 
                content=content,
                placement=placement,trigger="hover", options = list(container = "body"))
      
    )
  }else{
    box<-column(
      width=2,
      
      div(id=paste0(id,"_BS"),
          box(
            width=12,
            title = "",
            status = "danger",
            solidHeader = F,
            multiInput(
              inputId = id,
              label = span(id=paste0(id,"_BS"),label),
              choices = choices,
              options = list(enable_search = TRUE, non_selected_header = "Choose between:",
                             selected_header = "You have selected:",limit=1),
              width = "100%"
            )
          )
      ),
      
      bsPopover(id=paste0(id,"_BS"), 
                title=title, 
                content=content,
                placement=placement,trigger="hover", options = list(container = "body"))
      
    )
  }
  return(box)
}


mod_comp_box2<<-function(id,label,title,content, placement="right",choices=colnames(demo_trait_type1_dat)[c(1,2,3,4,6,10,11,12)],single=FALSE){
  if(!single){
    box<-column(
      width=2,
      
      div(id=paste0(id,"_BS"),
         

            multiInput(
              inputId = id,
              label = span(id=paste0(id,"_BS"),label),
              choices = choices,
              options = list(enable_search = TRUE, non_selected_header = "Choose between:",
                             selected_header = "You have selected:"),
              width = "100%"
            )
          
      ),
      
      bsPopover(id=paste0(id,"_BS"), 
                title=title, 
                content=content,
                placement=placement,trigger="hover", options = list(container = "body"))
      
    )
  }else{
    box<-column(
      width=2,
      
      div(id=paste0(id,"_BS"),
          box(
            width=12,
            title = "",
            status = "danger",
            solidHeader = F,
            multiInput(
              inputId = id,
              label = span(id=paste0(id,"_BS"),label),
              choices = choices,
              options = list(enable_search = TRUE, non_selected_header = "Choose between:",
                             selected_header = "You have selected:",limit=1),
              width = "100%"
            )
          )
      ),
      
      bsPopover(id=paste0(id,"_BS"), 
                title=title, 
                content=content,
                placement=placement,trigger="hover", options = list(container = "body"))
      
    )
  }
  return(box)
}

isValidEmail <<- function(x) {
  grepl("\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>", as.character(x), ignore.case=TRUE)
}

submit_confirm <<- function(id){
  modalDialog(
    "Are you sure to submit (please check your data and model once more)? If you click \"OK\", the page will be refreshed",
    title = "Submit?",
    footer = tagList(
      actionButton("cancel", "Cancel"),
      actionButton(paste0(id,"_ok"), "OK", class = "btn btn-danger")
    )
    
  )

}
email_check <<- modalDialog(
    "Your email address is invalid. Please check",
    title = "Email address error",
    footer = tagList(
      actionButton("cancel", "Cancel", class = "btn btn-warning"),
    )
  )


gs_method_check <<- modalDialog(
  div(
    p("Phenotype+pedigree+genotype files => ssGBLUP;"),
    p("Phenotype+genotype files => GBLUP;"),
    p("Phenotype+pedigree files => BLUP;"),
    p("Any combination above + variance componet file => corresponding calculation without variance component estimation;")
  ),
  title = "Missing some files?",
  footer = tagList(
    actionButton("cancel", "Cancel", class = "btn btn-warning"),
  )
)

jgs_method_check <<- modalDialog(
  div(
    p("Phenotype+pedigree+genotype files => ssGBLUP;"),
    p("Phenotype+genotype files => GBLUP;")
  ),
  title = "Missing some files?",
  footer = tagList(
    actionButton("cancel", "Cancel", class = "btn btn-warning"),
  )
)

cpi_method_check <<- modalDialog(
  div(
    p("Phenotype+genotype files => GBLUP;")
  ),
  title = "Missing some files?",
  footer = tagList(
    actionButton("cancel", "Cancel", class = "btn btn-warning"),
  )
)


source("homePanel.R")
source("imputePanel.R")
source("phePanel.R")
source("gsPanel.R")
source("cpiPanel.R")


gsMenu<-navbarMenu(
  title = "Genomic Selection",
  icon = icon("dna"),
  imputePanel,
  phePanel,
  gsPanel,
  cpiPanel
)

source("jgsPanel.R")
source("jcpiPanel.R")
source("pdataPanel.R")
source("docPanel.R")
source("megPanel.R")


jgsMenu<-navbarMenu(
  title = "Joint Genomic Selection",
  icon = icon("mixcloud"),
  jgsPanel,
  jcpiPanel
)

ui<-navbarPage(
  id="tabs",
  title=span(
    "Alpha-Cattle",
    style='font-size:30px;margin-right:50px'
  ),
  
  position = 'fixed-top',
  header = tagList(
    useShinydashboard(),
    useShinydashboardPlus(),
    use_waiter(),
    useShinyjs(),
    #extendShinyjs(text = "shinyjs.refresh = function() { location.reload(); }", functions = c("pageCol")),
    #setBackgroundColor(color = c("ghostwhite")),
    tags$head(
      tags$style( HTML(
        " .navbar-fixed-top{
        display: flex;
        }   
        .navbar-brand {
    FONT-VARIANT: JIS04;
    line-height: 0.6;
}"
      )),
    tags$link(rel="stylesheet", type="text/css", href="styles.css"),
    tags$style("@import url(https://use.fontawesome.com/releases/v6.1.1/css/all.css);")
    )
              
  ),
  theme = shinythemes::shinytheme(theme = "flatly"),

  homePanel,
  gsMenu,
  jgsMenu,
  pdataPanel,
  megPanel,
  docPanel
  
  
  #footer = div(class="footer",
   #             includeHTML("footer.html")
   #)
  
)

server<-function(input,output,session){
  
  source("jgsServer.R")
  source("gsServer.R")
  source("pheServer.R")
  source("imputeServer.R")
  source("cpiServer.R")
  source("jcpiServer.R")
  source("pdataServer.R")
  source("magServer.R")
 
  server_gs(input,output,session)
  server_jgs(input,output,session)
  pheServer(input,output,session)
  pdataServer(input,output,session)
  imputeServer(input,output,session)
  server_cpi(input,output,session)
  server_jcpi(input,output,session)
  magServer(input, output, session)
  output$pic<-renderImage({
    list(src="19.jpg",
         width="70%")
  },deleteFile = FALSE)
  output$pic_pca<-renderImage({
    list(src="pca_accuracy.jpg",
         width="90%")
  },deleteFile = FALSE)
  output$pic_cor<-renderImage({
    list(src="cor_he.jpg",
         width="90%")
  },deleteFile = FALSE)
  output$pic_joint<-renderImage({
    list(src="joint.jpg",
         width="90%")
  },deleteFile = FALSE)
  observeEvent(input$tabBut, {
    #browser()
    updateNavbarPage(session = session, inputId = "tabs", selected = "About HCGSP")
  })
  observeEvent(input$btn_phenotype, {
  updateNavbarPage(session = session, inputId = "tabs", selected = "PhenotypeProcessingPage")
})

observeEvent(input$btn_genomic_selection, {
  updateNavbarPage(session = session, inputId = "tabs", selected = "RegularGenomicSelectionPage")
})

observeEvent(input$btn_joint_selection, {
  updateNavbarPage(session = session, inputId = "tabs", selected = "JointGenomicSelectionPage")
})

observeEvent(input$btn_cpi_gcpi, {
  updateNavbarPage(session = session, inputId = "tabs", selected = "JointGCPICalculationPage")
})
  
}


shinyApp(ui, server)

