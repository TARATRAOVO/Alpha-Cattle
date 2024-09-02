cpiPanel<-tabPanel(
  title = "CPI or GCPI Calculation",
  div(class = "inlay", style = "height:70px;width:100%;background-color: white;"),
  column(
    width=12,
    h1("Upload files"),
  ),div(class= "inlay", style = "height:90px;width:100%;background-color: white;"),

    sidebarLayout(
      sidebarPanel(
         width=3,
         div(style = 'overflow-y:scroll;height:600px;',
             div(id="cpi_phe_upload_BS",style="margin-top:30px",fileInput("cpi_phe_upload", "Upload CPI related phenotypic file (.txt or .csv)")),
             bsPopover(id="cpi_phe_upload_BS",
                       title="CPI Related Phenotype file:",
                       content="This file is required for all CPI calculation process.and the first column must be Animal ID. The specific format can refer from the example table on the right of the page",
                       placement="right",trigger="hover" , options = list(container = "body")),
             uiOutput('cpi_phe_summary'),

             div(id="cpi_ped_upload_BS",style="margin-top:30px",fileInput("cpi_ped_upload", "Upload CPI related pedigree file (.txt or .csv)")),
             bsPopover(id="cpi_ped_upload_BS",
                       title="CPI Related Pedigree file:",
                       content="This file is required for BLUP to calculate EBV for CPI, and the first column must be Animal ID; second is sire ID; third is dam ID; fourth is order. The specific format can refer from the example table on the right of the page",
                       placement="right",trigger="hover" , options = list(container = "body")),
              uiOutput('cpi_ped_summary'),
              uiOutput('cpi_phe_ped_summary'),

              div(id="cpi_gen_upload_BS",style="margin-top:50px",fileInput("cpi_gen_upload", "Upload genotype files (.bed, .bim, .fam)", multiple=T)),
              bsPopover(id="cpi_gen_upload_BS",
                        title="Genotype file:",
                        content="These files are required for GBLUP to calculate GEBV for GCPI. The format is in Plink bed format and can be separated in different chromosomes and each chromosome contains .bed, .bim and .fam files",
                        placement="right",trigger="hover", options =list(container = "body")),
               uiOutput('cpi_gen_summary'),
               uiOutput('cpi_phe_gen_summary'),
               uiOutput('cpi_ped_gen_summary'),

               div(id="typecpi",style="margin-top:30px",selectInput(inputId="typecpi",label="Chose the method used in Milk Production Traits EBV calculation",choice=c("","Test_Day_Model","Average_305Days"),selected=F)),
               bsPopover(id="typecpi",
                         title="Method:",
                         content="You have to choose a Method to use in Milk Production Traits EBV calculation process.",
                         placement="right",trigger="hover", options = list(container = "body")),
                         uiOutput('typecpi_summary')

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
                          tags$li("Different kinds of CPI calculations need different Upload files. If you upload Phenotype and Ped file:",
                                  tags$ul(
                                    tags$li("CPI1/CPI3--Your Phenotype file should contain: Milk, Fatpct, Propct, SCS, Type, MS and FL"),
                                    tags$li("CPI2--Your Phenotype file should contain: Milk, Fatpct, Propct and SCS"),
                                    tags$li("CPI1 is applicable to domestic cows with both production performance and body type identification results, as well as progeny determination bulls"),
                                    tags$li("CPI2 is applicable to domestic cows with only production performance data and progeny determination bulls"),
                                    tags$li("CPI3 is applicable to the verification of progeny test results of IMPORTED bulls"),
                                   )),
                           div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                           tags$li("Different kinds of CPI calculations need different Upload files. If you upload Phenotype and Genotype file:",
                                   tags$ul(
                                     tags$li("GCPI--You should upload both Phenotype file and Genotype file, the file formats are shown as follow"),
                                     tags$li("NOTE that if you want to calculate GEBV only through the individual genome uploaded by yourself, you need to upload the phenotype file. Else if you want to combine our reference genome, please go to Joint Genomic Selection page"),
                                     tags$li("Your Phenotype file should contain: Milk, Fatpct, Propct, SCS, Type, MS, FL, we will first calculate the GEBV for each trait with the file you upload"),
                                     tags$li("Then we will use the GEBV to calculate GCPI for you"),
                                     tags$li("GCPI is applicable to domestic dairy cows and progeny determination bulls with both production performance and body type identification results and genome information"), 
                                   )),
                            div(class= "inlay", style = "height:20px;width:100%;background-color: white;"),
                            tags$li("The symbols corresponding to various traits and the standard deviation of domestic verified bulls are shown in the table below:",
                                    tags$ul(
                                    DT::dataTableOutput("traits_symbol_table",width="50%")
                                   )),
                            div(class= "inlay", style = "height:20px;width:100%;background-color: white;"),
                            tags$li("Different kinds of CPI are calculated as follow, they are set according to the 2016th standards of China Dairy Cattle Association:",
                                    tags$ul(
                                      tags$li("CPI1=20*[30*EBV-Milk/459+15*EBV-Fatpct/0.16+25*EBV-Propct/0.08+5*EBV-Type/5+10*EBV-MS/5+5*EBV-FL/5-10*(EBV-SCS-3)/0.16]"),
                                      tags$li("CPI2=20*[30*EBV-Milk/459+15*EBV-Fatpct/0.16+25*EBV-Propct/0.08-10*(EBV-SCS-3)/0.16]"),
                                      tags$li("CPI3=20*[30*EBV-Milk/800+10*EBV-Fatpct/0.3+20*EBV-Propct/0.12+5*EBV-Type/5+15*EBV-MS/5+10*EBV-FL/5-10*(EBV-SCS-3)/0.46]"),
                                      tags$li("GCPI=20*[30*GEBV-Milk/800+15*GEBV-Fatpct/0.3+25*GEBV-Propct/0.12+5*GEBV-Type/5+10*GEBV-MS/5+5*GEBV-FL/5-10*(GEBV-SCS-3)/0.46]+80"),
                                    )),
                            div(class= "inlay", style = "height:20px;width:100%;background-color: white;"),
                            tags$li("If you both upload Phenotype file and Pedigree file, we will send you the EBV for each trait as well as CPI1/CPI2/CPI3.
                               If your Phenotype Trait file contains both production and bodyshape Traits: Milk, Fatpct, Propct, SCS, Type, MS and FL, we will send you both CPI1 and CPI3. 
                               Else if you upload Trait file only contains production Traits: Milk, Fatpct, Propct, SCS, we will send you only CPI2.
                               If you both upload Genotype file and Phenotype file, we will send you the GEBV for each trait as well as GCPI. And your Phenotype file should also contain Traits: Milk, Fatpct, Propct, SCS, Type, MS and FL.
                               Check that Phenotype file in to calculate CPI1/CPI3 contains 7 Traits (the other columns are other effects) while CPI2 Phenotype file contains 5 Traits."),
                             div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                             tags$li("Submit your email address and wait the results.")
                     ),
                     br(),
                     div(tags$b("File format instruction", style='font-size:25px;')),
                     br(),
                     tags$ol(class="tutorial",
                             tags$li("Phenotype file in type 1(for CPI1/CPI3 calculation):",
                                     tags$ul(
                                        tags$li("Phenotype file can be comma or white space separated"),
                                        tags$li("Phenotype file MUST have column names (header)"),
                                        tags$li("The first column of Phenotype file MUST be animal ID"),
                                        tags$li("Other columns can be other variables that may be needed in the model"),
                                        tags$li("The missing value is also allowed and can be indicated as NA"),
                                        tags$li("The Last 7 columns MUST be the record of 7 traits: Milk, Fatpct, Propct, SCS, Type, MS and FL"),
                                         tags$li("An example of Type 1 Phenotype file is below:"),
                                         div(class= "inlay", style = "height:20px;width:100%;background-color: white;"),
                                         DT::dataTableOutput("demo_cpi_phe1",width="90%") 
                                      )),
                             div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                             tags$li("Phenotype file in type 2(for CPI2 calculation):",
                                     tags$ul(
                                          tags$li("Phenotype file can be comma or white space separated"),
                                          tags$li("Phenotype file MUST have column names (header)"),
                                          tags$li("The first column of Phenotype file MUST be animal ID"),
                                          tags$li("Other columns can be other variables that may be needed in the model"),
                                          tags$li("The missing value is also allowed and can be indicated as NA"),
                                          tags$li("The Last 5 columns MUST be the record of 5 traits: Milk, Fatpct, Propct and SCS"),
                                           tags$li("An example of Type 2 Phenotype file is below:"),
                                           div(class= "inlay", style = "height:20px;width:100%;background-color: white;"),
                                           DT::dataTableOutput("demo_cpi_phe2",width="90%")
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
                                            tags$li("The missing parent is allowed and can be indicated as NA, but missing order is not allowed. You can try package visPedigree for correct Pedigree format in R"),
                                            tags$li("An example of phenotype file is below:"),
                                            div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                                            DT::dataTableOutput("demo_cpi_ped",width="90%")
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
                                           DT::dataTableOutput("demo_cpi_fam",width="70%")
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

      mod_comp_box1(id="cpi_trt1",label="Traits to be calculated :",title="Traits",
                  content="You MUST choose traits the SAME as the traits in your trait file. No more and No less, which means you should choose Milk, Fatpct, Propct, SCS, Type, MS and FL and upload them in your Phenotype file!"),
      mod_comp_box1(id="cpi_fct1",label="Fixed effect (factor) :", title="Fixed effects",
                  content="You can specify multiple fixed factors"),
      mod_comp_box1(id="cpi_cov1",label="Fixed effect (covariate) :", title="Covariate",
                  content="You can specify multiple covariate varibles which are also considered fixed effects"),
      mod_comp_box1(id="cpi_gen1",label="Random effect (genetic) :", title="Genetic random effects",
                  content="This is for the  random effect that the pedigree-based or genotype-based or their combined relationship matrix is exerted on. Please note that HCGSP now only support unique random genetic effect. You may supply more varibles but only the first one is considered"),
      mod_comp_box1(id="cpi_rnd1",label="Random effect (other) :", title="Non-genetic random effects",
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

      mod_comp_box2(id="cpi_trt2",label="Traits to be calculated :",title="Traits",
                  content="You MUST choose traits the SAME as the traits in your trait file. No more and No less, which means you should choose Milk, Fatpct, Propct, SCS, Type, MS and FL and upload them in your Phenotype file!"),
      mod_comp_box2(id="cpi_fct2",label="Fixed effect (factor) :", title="Fixed effects",
                  content="You can specify multiple fixed factors"),
      mod_comp_box2(id="cpi_cov2",label="Fixed effect (covariate) :", title="Covariate",
                  content="You can specify multiple covariate varibles which are also considered fixed effects"),
      mod_comp_box2(id="cpi_gen2",label="Random effect (genetic) :", title="Genetic random effects",
                  content="This is for the  random effect that the pedigree-based or genotype-based or their combined relationship matrix is exerted on. Please note that HCGSP now only support unique random genetic effect. You may supply more varibles but only the first one is considered"),
      mod_comp_box2(id="cpi_rnd2",label="Random effect (other) :", title="Non-genetic random effects",
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
         htmlOutput("cpi_model_check1"),
         htmlOutput("cpi_model_check2"),
         h4("GBLUP Model Need At Least:"),
         h5("Milk Production Traits=Mean(fix)+Herd(fix)+YearSeason(fix)+Parity(fix)+ID(gen)+pe(ID(rnd))+e"),
         h5("Body Shape Traits=Mean(fix)+Herd(fix)+YearSeason(fix)+Age(fix)+sex(fix,if both sire and dam included)+ID(gen)+pe(ID(rnd))+e")
        ),
        column(
          width=6,
          h1("Submit job"),
          textInput(inputId = "cpi_email", label = "E-mail", value = "Type your email address...", width = "400px"),
          bsPopover("cpi_email", title='Your e-mail address', content='We do not do anything with it...apart from sending your CPI results. :)', placement="right",trigger="hover", options = list(container = "body")),
          h5("Once analysis is done, you'll receive results by email."),
          actionButton(
            inputId = "cpi_submit_bttn",
            label = "Submit!",
            class = "btn-info",
            icon=icon("paper-plane"),
            style='font-size:150%;width:280px'
           )
         )
       )
     )
   )

