gsPanel<-tabPanel(
  title = "Regular Genomic Selection",
  value = "RegularGenomicSelectionPage",
  div(class = "inlay", style = "height:70px;width:100%;background-color: white;"),
  column(
    width=12,
    h1("Upload files"),
  ),div(class = "inlay", style = "height:90px;width:100%;background-color: white;"),
  
    sidebarLayout(
      sidebarPanel(
        width=3,
        div(style = "overflow-y:scroll;height:600px;",

                                                  div(id="typegs",style="margin-top:30px",selectInput(inputId="typegs",label="Chose the Phenotype Type",choice=c("","Milk_Production_Traits_inTDM","Milk_Production_Traits_inAVE","Reproductive_Traits","Body_Type_Traits"),selected=F)),
                                                  bsPopover(id="typegs",
                                                            title="Type:",
                                                            content="You have to choose a type to determine whether the Test Day Model will be used or not.",
                                                            placement="right",trigger="hover", options = list(container = "body")),
                                                  uiOutput('typegs_summary'),

            div(id="gs_phe_upload_BS",style="margin-top:30px",fileInput("gs_phe_upload", "Upload phenotypic file (.txt or .csv)")),
            bsPopover(id="gs_phe_upload_BS", 
                      title="Phenotype file:", 
                      content="This file is required, and the first column must be Animal ID",
                      placement="right",trigger="hover", options = list(container = "body")),
            uiOutput('gs_phe_summary'),
           
            div(id="gs_ped_upload_BS",style="margin-top:50px",fileInput("gs_ped_upload", "Upload pedigree file (.txt or .csv)")),
            bsPopover(id="gs_ped_upload_BS", 
                      title="Pedigree file:", 
                      content="This file is required for BLUP or ssGBLUP. The first column is Animal ID; second is sire ID; third is dam ID; fourth is order",
                      placement="right",trigger="hover", options = list(container = "body")),
            uiOutput('gs_ped_summary'),
            uiOutput('gs_phe_ped_summary'),
       
            
            div(id="gs_gen_upload_BS",style="margin-top:50px",fileInput("gs_gen_upload", "Upload genotype files (.bed, .bim, .fam)", multiple=T)),
            bsPopover(id="gs_gen_upload_BS", 
                      title="Genotype file:", 
                      content="These files are required for GBLUP or ssGBLUP. The format is in Plink bed format and can be separated in different chromosomes and each chromosome contains .bed, .bim and .fam files",
                      placement="right",trigger="hover", options = list(container = "body")),
            uiOutput('gs_gen_summary'),
            uiOutput('gs_phe_gen_summary'),
            uiOutput('gs_ped_gen_summary'),
   
            div(id="gs_var_upload_BS",style="margin-top:50px",fileInput("gs_var_upload", "Upload variance component files (.txt or .csv)")),
            bsPopover(id="gs_var_upload_BS",
                      title="Variance component file:", 
                      content="This file is required when you want to skip the variance component estimation step and use the components you provide to predict EBV directly.",
                      placement="right",trigger="hover", options = list(container = "body")),
            uiOutput('gs_var_summary')
       
        )
      ),

     
      mainPanel(
        width=9,
        div(style = 'overflow-y:scroll;height:635px;',
            tags$div(
              class = "document",
              div(tags$b("Tutorial", style='font-size:40px;')),
              br(),
              div(tags$b("Work flow", style='font-size:25px;')),
              br(),
              tags$ol(class="tutorial",
                       tags$li("We will use different computational models for different phenotypic types:",
                              tags$ul(
                                tags$li("Select Milk Production Traits: Calculating breeding values using the Test Day Model"),
                                tags$li("Select Reproductive Traits: Calculating breeding values using the Repeatability Model. Please remember add permanent environmental effects"),
                                tags$li("Select Body Type Traits: Calculating breeding values using the common Animal Model."),
                              )),
                      div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                      tags$li("Upload files and different combination of files means different methods:",
                              tags$ul(
                                tags$li("Upload phenotype file + pedigree file + genotype file: sinle-step genomic best linear unbiased prediction (ssGBLUP)"),
                                tags$li("Upload phenotype file + pedigree file : best linear unbiased prediction (BLUP)"),
                                tags$li("Upload phenotype file + genotype file : genomic best linear unbiased prediction (GBLUP)"),
                                tags$li("Upload variance component file + Any of the file combinations above: predict estimated breeding values without variance components estimation (fastest)"),
                              )),
                      div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                      tags$li("Specify model components based on column names in uploaded phenotype file. 
                        Please note that only single-trait analysis is supported at present version, 
                        although multiple traits can be specified. When multiple traits are specified, 
                        the traits will be analysed one by one using the same model."),
                      div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                      tags$li("Check the model you specified at the left bottom."),
                      div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                      tags$li("Submit your email address and wait the results.")
              ),
              br(),
              div(tags$b("File format instruction", style='font-size:25px;')),
              br(),
              tags$ol(class="tutorial",
                      tags$li("Phenotype file:",
                              tags$ul(
                                tags$li("Phenotype file can be comma or white space separated"),
                                tags$li("Phenotype file MUST have column names (header)"),
                                tags$li("The first column of phenotype file MUST be animal ID"),
                                tags$li("Other columns can be other variables that may be needed in the model"),
                                tags$li("The missing value is also allowed and can be indicated as NA"),
                                tags$li("An example of phenotype file is below:"),
                                div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                                DT::dataTableOutput("demo_gs_phe",width="90%")
                              )),
                      div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                      tags$li("Pedigree file:",
                              tags$ul(
                                tags$li("Pedigree file can be comma or white space separated"),
                                tags$li("Pedigree file MUST have column names (header)"),
                                tags$li("The first column of phenotype file MUST be animal ID; 
                                        the second column MUST be the sire of individuals in the first column;
                                        the third column MUST be the dam of individuals in the first column;
                                        the fourth column MUST be the order variables, such as birth date or
                                        other numberic values that can be used to sort the pedigree"),
                                tags$li("The missing parent is allowed and can be indicated as NA"),
                                tags$li("An example of pedigree file is below:"),
                                div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                                DT::dataTableOutput("demo_gs_ped",width="70%")
                              )),
                      div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                      tags$li("Genotype file:",
                              div("The genotype files were in Plink bed format, 
                                        and the standard format can be found ",
                                  tags$a("here", href = "http://zzz.bwh.harvard.edu/plink/binary.shtml"),
                                  ". Besides going by basic rules of Plink,
                                  the fam file is required to be prepared using the following rules:"),
                              tags$ul(
                                tags$li("No header is allowed"),
                                tags$li("The first and second columns of fam file MUST both be animal ID"),
                                tags$li("An example of fam file is below:"),
                                div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                                DT::dataTableOutput("demo_gs_fam",width="70%")
                              )),
                      div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                      tags$li("Variance component file:",
                              tags$ul(
                                tags$li("No header is allowed"),
                                tags$li("Each row represents the variance comonents for each trait, 
                                        and each column represents the variance components 
                                        for each random effect across different traits"),
                                tags$li("An example of variance component file is below, corresponding to the model containing 
                                        three random effect in total, animal effect, dam effect, and random residual effect 
                                        for two traits:"),
                                div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                                DT::dataTableOutput("demo_gs_var",width="50%")
                              ))
              )
            ))
        )
      
  
  ),
  br(),
  
  column(
    width=12,
    h1("Specify model components")
  ),div(class = "inlay", style = "height:100px;width:100%;background-color: white;"),
  
  fluidRow(
    width=12,
    column(
      width=12,
      id="model_box",
      
      mod_comp_box(id="gs_trt",label="Traits to be calculated :",title="Traits",
	  content="You can specify multiple traits, which will be analyzed one by one"),
      mod_comp_box(id="gs_fct",label="Fixed effect (factor) :", title="Fixed effects",
                   content="You can specify multiple fixed factors"),
      mod_comp_box(id="gs_cov",label="Fixed effect (covariate) :", title="Covariate",
                   content="You can specify multiple covariate varibles which are also considered fixed effects"),
      mod_comp_box(id="gs_gen",label="Random effect (genetic) :", title="Genetic random effects",
                   content="This is for the  random effect that the pedigree-based or genotype-based or their combined relationship matrix is exerted on. Please note that HEGS now only support unique random genetic effect. You may supply more varibles but only the first one is considered"),
      mod_comp_box(id="gs_rnd",label="Random effect (other) :", title="Non-genetic random effects",
                   content="This is for the non-genetic random effect, such as permenent environmental effect, dam effect, etc.",
                   placement="left")
    )
    
  ),
  
  br(),
  
  fluidRow(
    column(
      width=12,
      column(
        width=6,
        h1("Check model"),
        htmlOutput("gs_model_check")
      ),
      column(
        width=6,
        h1("Submit job"),
        textInput(inputId = "gs_email", label = "E-mail", value = "Type your email address...", width = "400px"),
        bsPopover("gs_email", title='Your e-mail address', content='We do not do anything with it...apart from sending your GS results. :)', placement="right",trigger="hover", options = list(container = "body")),
        h5("Once analysis is done, you'll receive results by email."),
        # maybe also good to put a submit button -- to avoid refreshes that are annoying and slow down application
        #submitButton(text = "Submit!", icon("paper-plane"), width="250px")
        actionButton(
          inputId = "gs_submit_bttn",
          label = "Submit!",
          #width = "100%",
          class = "btn-info",
          icon=icon("paper-plane"),
          style='font-size:150%;width:280px'
          #style="color: #fff; background-color: #337ab7; border-color: light-blue"
        )
      )
    )
  )
)