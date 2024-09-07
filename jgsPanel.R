jgsPanel<-tabPanel(
  title = "Joint Genomic Selection",
  value = "JointGenomicSelectionPage",
  div(class = "inlay", style = "height:70px;width:100%;background-color: white;"),
  column(
    width=12,
    h1("Upload files")
  ),div(class = "inlay", style = "height:90px;width:100%;background-color: white;"),
  
  sidebarLayout(
    sidebarPanel(
      width=3,
        div(style="overflow-y:scroll;height:600px;",

			  div(id="breed",style="margin-top:30px",selectInput(inputId="breed",label="Chose the breed",choice=c("","Holstein"),selected=F)),
			  bsPopover(id="breed", 
			            title="Breed:", 
			            content="You have to choose a breed to combine your data with the existing data for analysis.",
			            placement="right",trigger="hover", options = list(container = "body")),
			  uiOutput('breed_summary'),

                                                  div(id="type",style="margin-top:30px",selectInput(inputId="type",label="Chose the Phenotype Type",choice=c("","Milk_Production_Traits_inTDM","Milk_Production_Traits_inAVE","Reproductive_Traits","Body_Type_Traits"),selected=F)),
                                                  bsPopover(id="type",
                                                            title="Type:",
                                                            content="You have to choose a type to determine whether the Test Day Model will be used or not.",
                                                            placement="right",trigger="hover", options = list(container = "body")),
                                                  uiOutput('type_summary'),
		
            div(id="jgs_phe_upload_BS",style="margin-top:50px",fileInput("jgs_phe_upload", "Upload phenotypic file (.txt or .csv)")),
            bsPopover(id="jgs_phe_upload_BS", 
                      title="Phenotype file:", 
                      content="This file is required, and the first column must be Animal ID",
                      placement="right",trigger="hover", options = list(container = "body")),
            uiOutput('jgs_phe_summary'),

            div(id="jgs_ped_upload_BS",style="margin-top:50px",fileInput("jgs_ped_upload", "Upload pedigree file (.txt or .csv)")),
            bsPopover(id="jgs_ped_upload_BS", 
                      title="Pedigree file:", 
                      content="This file is required for ssGBLUP. The first column is Animal ID; second is sire ID; third is dam ID; fourth is order",
                      placement="right",trigger="hover", options = list(container = "body")),
            uiOutput('jgs_ped_summary'),
            uiOutput('jgs_phe_ped_summary'),

            div(id="jgs_gen_upload_BS",style="margin-top:50px",fileInput("jgs_gen_upload", "Upload genotype files (.bed, .bim, .fam)", multiple=T)),
            bsPopover(id="jgs_gen_upload_BS", 
                      title="Genotype file:", 
                      content="These files are required for GBLUP or ssGBLUP. The format is in Plink bed format and can be separated in different chromosomes and each chromosome contains .bed, .bim and .fam files",
                      placement="right",trigger="hover", options = list(container = "body")),
            uiOutput('jgs_gen_summary'),
            uiOutput('jgs_phe_gen_summary'),
            uiOutput('jgs_ped_gen_summary')
            
        
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
                       tags$li("Upload files and different combination of files means different methods using combined population:",
                              tags$ul(
                                tags$li("Upload phenotype file + pedigree file + genotype file: sinle-step genomic best linear unbiased prediction (ssGBLUP)"),
                                tags$li("Upload phenotype file + genotype file : genomic best linear unbiased prediction (GBLUP)"),
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
                                DT::dataTableOutput("demo_jgs_phe",width="90%")
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
                                DT::dataTableOutput("demo_jgs_ped",width="70%")
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
                                DT::dataTableOutput("demo_jgs_fam",width="70%")
                              ))
              )        
            )
        )
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
      
      mod_comp_box(id="jgs_trt",label="Traits to be calculated :",title="Traits",content=
                     "You can specify multiple traits, which will be analyzed one by one"),
      mod_comp_box(id="jgs_fct",label="Fixed effect (factor) :", title="Fixed effects",
                   content="You can specify multiple fixed factors"),
      mod_comp_box(id="jgs_cov",label="Fixed effect (covariate) :", title="Covariate",
                   content="You can specify multiple covariate varibles which are also considered fixed effects"),
      mod_comp_box(id="jgs_gen",label="Random effect (genetic) :", title="Genetic random effects",
                   content="This is for the  random effect that the pedigree-based or genotype-based or their combined relationship matrix is exerted on. Please note that HEGS now only support unique random genetic effect. You may supply more varibles but only the first one is considered"),
      mod_comp_box(id="jgs_rnd",label="Random effect (other) :", title="Non-genetic random effects",
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
        htmlOutput("jgs_model_check")
      ),
      column(
        width=6,
        h1("Submit job"),
        textInput(inputId = "jgs_email", label = "E-mail", value = "Type your email address...", width = "400px"),
        bsPopover("jgs_email", title='Your e-mail address', content='We do not do anything with it...apart from sending your GS results. :)', placement="right",trigger="hover", options = list(container = "body")),
        h5("Once analysis is done, you'll receive results by email."),
        # maybe also good to put a submit button -- to avoid refreshes that are annoying and slow down application
        #submitButton(text = "Submit!", icon("paper-plane"), width="250px")
        actionButton(
          inputId = "jgs_submit_bttn",
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