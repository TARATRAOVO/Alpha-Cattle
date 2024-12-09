imputePanel <- tabPanel(
  title = "Genomic Chip Imputation",
  div(class = "inlay", style = "height:120px;width:100%;background-color: white;"),
  column(
    width = 12,
    h1("Upload files")
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(
        style = "overflow-y:scroll;height:600px;",
        div(id = "density_BS", style = "margin-top:30px", selectInput(inputId = "density", label = "Chose the imputation density", choice = c("", "150K_SNP", "Whole_Genome"), selected = F)),
        bsPopover(
          id = "density_BS",
          title = "Density:",
          content = "You need to determine the target density of your chip filling process.",
          placement = "right", trigger = "hover", options = list(container = "body")
        ),
        uiOutput("density_summary"),
        div(id = "impute_gen_upload_BS", style = "margin-top:50px", fileInput("impute_gen_upload", "Upload genotype files (.bed, .bim, .fam)", multiple = T)),
        bsPopover(
          id = "impute_gen_upload_BS",
          title = "Genotype file:",
          content = "These files are required for Imputation. Please upload your genome sequencing chip in designative format. The format is in Plink bed format and can be separated in different chromosomes and each chromosome contains .bed, .bim and .fam files",
          placement = "right", trigger = "hover", options = list(container = "body")
        ),
        uiOutput("impute_gen_summary"),
      )
    ),
    mainPanel(
      width = 9,
      div(
        style = "overflow-y:scroll;height:635px;",
        tags$div(
          class = "document",
          div(tags$b("Tutorial", style = "font-size:40px;")),
          br(),
          div(tags$b("Work flow", style = "font-size:25px;")),
          br(),
          tags$ol(
            class = "tutorial",
            tags$li("Determine the target density of your chip filling process. Note that it may take more time to impute into whole genome and 150k density chip is enough for genomic prediction such as GBLUP."),
            div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
            tags$li("Upload your genome chip sequencing data in specific format(.bed, .bim, .fam).
                         You need to convert your chip data to plink format before uploading.
                         Note that the upload file can be separated in different chromosomes and each chromosome contains .bed, .bim and .fam files."),
            div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
            tags$li("Submit your email address and wait the results.")
          ),
          br(),
          div(tags$b("File format instruction", style = "font-size:25px;")),
          br(),
          tags$ol(
            class = "tutorial",
            tags$li(
              "Genotype file:",
              div(
                "The genotype files were in Plink bed format,
                                        and the standard format can be found ",
                tags$a("here", href = "http://zzz.bwh.harvard.edu/plink/binary.shtml"),
                ". Besides going by basic rules of Plink,
                                  the fam file is required to be prepared using the following rules:"
              ),
              tags$ul(
                tags$li("No header is allowed"),
                tags$li("The first and second columns of fam file MUST both be animal ID"),
                tags$li("An example of fam file is below:"),
                div(class = "inlay", style = "height:20px;width:100%;background-color: white;"),
                DT::dataTableOutput("demo_impute_fam", width = "70%")
              )
            )
          )
        )
      )
    )
  ),
  br(),
  fluidRow(
    column(
      width = 12,
      column(
        width = 6,
        h1("Reference Article"),
        h5("The article about the accuracy of Holstein cow chip imputation and the impact of chip density on the accuracy of GBLUP genome prediction can be referred to:"),
        h4("Nazzicari N, Biscarini F. Stacked kinship CNN vs. GBLUP for genomic predictions of additive and complex continuous phenotypes. Sci Rep. 2022 Nov 18;12(1):19889."),
        h5("You can refer to the conclusions of the article and consider your own needs, then choose to impute into 150K_SNP or whole genome.")
      ),
      column(
        width = 6,
        h1("Submit job"),
        textInput(inputId = "impute_email", label = "E-mail", value = "Type your email address...", width = "400px"),
        bsPopover("impute_email", title = "Your e-mail address", content = "We do not do anything with it...apart from sending your GS results. :)", placement = "right", trigger = "hover", options = list(container = "body")),
        h5("Once analysis is done, you'll receive results by email."),
        # maybe also good to put a submit button -- to avoid refreshes that are annoying and slow down application
        # submitButton(text = "Submit!", icon("paper-plane"), width="250px")
        actionButton(
          inputId = "impute_submit_bttn",
          label = "Submit!",
          # width = "100%",
          class = "btn-info",
          icon = icon("paper-plane"),
          style = "font-size:150%;width:280px"
          # style="color: #fff; background-color: #337ab7; border-color: light-blue"
        )
      )
    )
  )
)
