docPanel <- tabPanel(
  div(class = "inlay", style = "height:70px;width:100%;background-color: white;"),
  title = "About Alpha-Cattle",
  value = "About HCGSP",
  icon = icon("dochub"),

  # 页面标题
  fluidRow(
    column(12, 
      h1("Alpha-Cattle: Holistic Cattle Genomic Selection Platform", style = "text-align: center; margin-top: 20px;"),
      p("Welcome to Alpha-Cattle, a comprehensive platform dedicated to enhancing cattle breeding programs through the power of advanced genomic selection. Our goal is to provide farmers, ranchers, and researchers with the tools and insights necessary to optimize genetic improvement in cattle populations.", 
        style = "text-align: center; font-size: 18px; margin-top: 20px;")
    )
  ),

  # 项目概述
  fluidRow(
    column(10, offset = 1,
      h3("Project Overview", style = "margin-top: 30px;"),
      p("Alpha-Cattle integrates state-of-the-art genomic selection techniques with large-scale data from cattle populations to provide accurate predictions of breeding values. Our platform is designed to help users make informed decisions about cattle breeding based on comprehensive genetic and phenotypic data.", 
        style = "font-size: 16px;"),
      tags$ul(
        tags$li("Comprehensive Database: Includes over 3,800 reference populations with phenotype and genome-wide information."),
        tags$li("Accurate Genomic Predictions: Utilizing combined reference populations to enhance prediction accuracy."),
        tags$li("CPI and GCPI Calculations: Tailored trait evaluation for Holstein cattle, including milk yield, body type, and more.")
      )
    )
  ),

  # 项目优势
  fluidRow(
    column(10, offset = 1,
      h3("Why Choose Alpha-Cattle?", style = "margin-top: 30px;"),
      p("Alpha-Cattle offers unparalleled insights into genomic selection, helping users achieve greater genetic gain and optimize the traits that matter most. With advanced analytics and easy-to-use tools, Alpha-Cattle is your partner in modern cattle breeding.",
        style = "font-size: 16px;"),
      tags$ul(
        tags$li("Advanced Genomic Selection Models."),
        tags$li("Easy Upload and Management of Genomic Data."),
        tags$li("Actionable Insights for Breeding Decisions.")
      )
    )
  ),

  # 联系方式
  fluidRow(
    column(12, 
      h4("Get In Touch", style = "text-align: center; margin-top: 40px;"),
      p("For more information, contact us at: info@alpha-cattle.com", 
        style = "text-align: center; font-size: 16px;")
    )
  )
)