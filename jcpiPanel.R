jcpiPanel<-tabPanel(
  title = "Joint GCPI Calculation",
  value = "JointGCPICalculationPage",
  div(class = "inlay", style = "height:70px;width:100%;background-color: white;"),
  column(
    width=12,
    h1("Upload files"),
  ),div(class= "inlay", style = "height:90px;width:100%;background-color: white;"),

    sidebarLayout(
      sidebarPanel(
         width=3,
         div(style = 'overflow-y:scroll;height:600px;',
             div(id="jcpi_phe_upload_BS",style="margin-top:30px",fileInput("jcpi_phe_upload", "Upload GCPI related phenotypic file (.txt or .csv)")),
             bsPopover(id="jcpi_phe_upload_BS",
                       title="GCPI Related Phenotype file:",
                       content="This file is required for GBLUP to calculate GEBV for GCPI, and the first column must be Animal ID. The specific format can refer from the example table on the right of the page",
                       placement="right",trigger="hover" , options = list(container = "body")),
              uiOutput('jcpi_phe_summary'),

              div(id="jcpi_gen_upload_BS",style="margin-top:50px",fileInput("jcpi_gen_upload", "Upload genotype files (.bed, .bim, .fam)", multiple=T)),
              bsPopover(id="jcpi_gen_upload_BS",
                        title="Genotype file:",
                        content="These files are required for GBLUP to calculate GEBV for GCPI. The format is in Plink bed format and can be separated in different chromosomes and each chromosome contains .bed, .bim and .fam files",
                        placement="right",trigger="hover", options =list(container = "body")),
               uiOutput('jcpi_gen_summary'),
               uiOutput('jcpi_phe_gen_summary')

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
                           tags$li("This page is only for GCPI calculation. You should upload all 3 files:",
                                   tags$ul(
                                     tags$li("GCPI--You should upload both 2 files at the left, the file formats are shown follow"),
                                     tags$li("This page will combine the file you upload with our reference individuals, which will make the result more accurate"),
                                     tags$li("Your Phenotype file should contain: Milk, Fatpct, Pr opct, SCS, Type, MS, FL, we will first calculate the GEBV for each trait with the file you upload"),
                                     tags$li("Then we will use the GEBV to calculate GCPI for you"),
                                     tags$li("GCPI is applicable to domestic dairy cows and progeny determination bulls with both production performance and body type identification results and genome information"), 
                                   )),
                            div(class= "inlay", style = "height:20px;width:100%;background-color: white;"),
                            tags$li("The symbols corresponding to various traits and the standard deviation of domestic verified bulls are shown in the table below:",
                                    tags$ul(
                                    DT::dataTableOutput("jtraits_symbol_table",width="50%")
                                   )),
                            div(class= "inlay", style = "height:20px;width:100%;background-color: white;"),
                            tags$li("GCPI is calculated as follow, it is according to the 2016th standards of China Dairy Cattle Association:",
                                    tags$ul(
                                      tags$li("GCPI=20*[30*GEBV-Milk/800+15*GEBV-Fatpct/0.3+25*GEBV-Pr opct/0.12+5*GEBV-Type/5+10*GEBV-MS/5+5*GEBV-FL/5-10*(GEBV-SCS-3)/0.46]+80"),
                                    )),
                            div(class= "inlay", style = "height:20px;width:100%;background-color: white;"),
                            tags$li("Please upload your Phenotype file with 7 traits in the Symbol Box above in this page, which contains  Milk, Fatpct, Pr opct, SCS, Type, MS, FL.
                               See the file format introduction for details."),
                            div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                            tags$li("Submit your email address and wait the results.")
                     ),
                     br(),
                     div(tags$b("File format instruction", style='font-size:25px;')),
                     br(),
                     tags$ol(class="tutorial",
                              tags$li("Phenotype file",
                                       tags$ul(
                                            tags$li("Phenotype file can be comma or white space separated"),
                                            tags$li("Phenotype file MUST have column names (header)"),
                                            tags$li("The first column of Phenotype file MUST be animal ID"),
                                            tags$li("Other columns can be other variables that may be needed in the model"),
                                            tags$li("The missing value is also allowed and can be indicated as NA"),
                                            tags$li("The Last 7 columns MUST be the record of 7 traits: Milk, Fatpct, Propct, SCS, Type, MS and FL"),
                                            tags$li("If you do not have the record of some traits, fill the column of the traits in the phenotype file with NA"),
                                            tags$li("An example of phenotype file is below:"),
                                            div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                                            DT::dataTableOutput("demo_jcpi_phe",width="90%")
                                        )),
                               div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                               tags$li("Genotype file:",
                                       div("The genotype files were in Plink bed format,
                                                  and the standard format can be found ",
                                           tags$a("here",href = "http://zzz.bwh.harvard.edu/plink/binary.shtml"),
                                           ".Besides going by basic rules of Plink,
                                           the fam file is required to be prepared using the following rules:"),
                                        tags$ul(
                                           tags$li("No header is allowed"),
                                           tags$li("The first and second columns of fam file MUST both be animal ID"),
                                           tags$li("An example of fam file is below:"),
                                           div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                                           DT::dataTableOutput("demo_jcpi_fam",width="70%")
                                         ))
                 )
               ))
           )


  ),
  br(),

  column(
    width= 12,
    h1("Specify model components in Milk Production Traits")
  ),div(class = "inlay", style = "height:100px;width:100%;background-color: white;"),

  fluidRow(
    width=12,
    column(
      width=12,
      id="model_box",

      mod_comp_box1(id="jcpi_trt1",label="Traits to be calculated :",title="Traits",
                  content="You MUST choose traits the SAME as the traits in your trait file. No more and No less, which means you should choose Milk, Fatpct, Pr opct, SCS, Type, MS and FL and upload them in your Phenotype file!"),
      mod_comp_box1(id="jcpi_fct1",label="Fixed effect (factor) :", title="Fixed effects",
                  content="You can specify multiple fixed factors"),
      mod_comp_box1(id="jcpi_cov1",label="Fixed effect (covariate) :", title="Covariate",
                  content="You can specify multiple covariate varibles which are also considered fixed effects"),
      mod_comp_box1(id="jcpi_gen1",label="Random effect (genetic) :", title="Genetic random effects",
                  content="This is for the  random effect that the pedigree-based or genotype-based or their combined relationship matrix is exerted on. Please note that Alpha-Cattle now only support unique random genetic effect. You may supply more varibles but only the first one is considered"),
      mod_comp_box1(id="jcpi_rnd1",label="Random effect (other) :", title="Non-genetic random effects",
                  content="This is for the non-genetic random effect, such as permenent environmental effect, dam effect, etc.",
                  placement="left")
   )

  ),

  br(),


  column(
    width= 12,
    h1("Specify model components in Body Shape Traits")
  ),div(class = "inlay", style = "height:100px;width:100%;background-color: white;"),

  fluidRow(
    width=12,
    column(
      width=12,
      id="model_box",

      mod_comp_box2(id="jcpi_trt2",label="Traits to be calculated :",title="Traits",
                  content="You MUST choose traits the SAME as the traits in your trait file. No more and No less, which means you should choose Milk, Fatpct, Pr opct, SCS, Type, MS and FL and upload them in your Phenotype file!"),
      mod_comp_box2(id="jcpi_fct2",label="Fixed effect (factor) :", title="Fixed effects",
                  content="You can specify multiple fixed factors"),
      mod_comp_box2(id="jcpi_cov2",label="Fixed effect (covariate) :", title="Covariate",
                  content="You can specify multiple covariate varibles which are also considered fixed effects"),
      mod_comp_box2(id="jcpi_gen2",label="Random effect (genetic) :", title="Genetic random effects",
                  content="This is for the  random effect that the pedigree-based or genotype-based or their combined relationship matrix is exerted on. Please note that Alpha-Cattle now only support unique random genetic effect. You may supply more varibles but only the first one is considered"),
      mod_comp_box2(id="jcpi_rnd2",label="Random effect (other) :", title="Non-genetic random effects",
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
         htmlOutput("jcpi_model_check1"),
         htmlOutput("jcpi_model_check2"),
         h4("GBLUP Model Need At Least:"),
         h5("Milk Production Traits=Mean(fix)+Herd(fix)+YearSeason(fix)+Parity(fix)+ID(gen)+pe(ID(rnd))+e"),
         h5("Body Shape Traits=Mean(fix)+Herd(fix)+YearSeason(fix)+Age(fix)+sex(fix,if both sire and dam included)+ID(gen)+pe(ID(rnd))+e")
        ),
        column(
          width=6,
          h1("Submit job"),
          textInput(inputId = "jcpi_email", label = "E-mail", value = "Type your email address...", width = "400px"),
          bsPopover("jcpi_email", title='Your e-mail address', content='We do not do anything with it...apart from sending your GCPI results. :)', placement="right",trigger="hover", options = list(container = "body")),
          h5("Once analysis is done, you'll receive results by email."),
          actionButton(
            inputId = "jcpi_submit_bttn",
            label = "Submit!",
            class = "btn-info",
            icon=icon("paper-plane"),
            style='font-size:150%;width:280px'
           )
         )
       )
     )
   )

