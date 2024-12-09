docPanel <- tabPanel(
  div(class = "inlay", style = "height:120px;width:100%;background-color: white;"),
  title = "About Alpha-Cattle",
  value = "About HCGSP",
  icon = icon("dochub"),

  # Hero Section
  div(
    class = "hero-section",
    style = "background: linear-gradient(135deg, #1e88e5 0%, #1565c0 100%); 
           color: white; 
           padding: 60px 0; 
           margin-bottom: 40px;",
    fluidRow(
      column(10, offset = 1,
        h1("Alpha-Cattle", 
           style = "text-align: center; 
                  font-size: 3.5em; 
                  font-weight: 700; 
                  margin-bottom: 20px;"),
        h2("Holistic Cattle Genomic Selection Platform",
           style = "text-align: center; 
                  font-weight: 300; 
                  margin-bottom: 30px;"),
        p("Empowering farmers and researchers with advanced genomic selection tools", 
          style = "text-align: center; 
                 font-size: 1.2em; 
                 opacity: 0.9;")
      )
    )
  ),

  # 主要内容区域
  div(
    style = "max-width: 1200px; margin: 0 auto; padding: 0 20px;",
    
    # 项目概述卡片
    div(
      class = "content-card",
      style = "background: white; 
             border-radius: 8px; 
             box-shadow: 0 2px 4px rgba(0,0,0,0.1); 
             padding: 30px; 
             margin-bottom: 40px;",
      h3("Project Overview", 
         style = "color: #1565c0; 
                border-bottom: 2px solid #1565c0; 
                padding-bottom: 10px;"),
      p("Alpha-Cattle integrates state-of-the-art genomic selection techniques with large-scale data from cattle populations to provide accurate predictions of breeding values.", 
        style = "font-size: 16px; line-height: 1.6;"),
      
      div(
        style = "display: flex; flex-wrap: wrap; gap: 20px; margin-top: 20px;",
        div(
          style = "flex: 1; min-width: 250px; 
                 background: #f5f5f5; 
                 padding: 20px; 
                 border-radius: 6px;",
          icon("database"),
          h4("Comprehensive Database"),
          p("3,800+ reference populations with phenotype and genome-wide information")
        ),
        div(
          style = "flex: 1; min-width: 250px; 
                 background: #f5f5f5; 
                 padding: 20px; 
                 border-radius: 6px;",
          icon("chart-line"),
          h4("Accurate Predictions"),
          p("Enhanced prediction accuracy through combined reference populations")
        ),
        div(
          style = "flex: 1; min-width: 250px; 
                 background: #f5f5f5; 
                 padding: 20px; 
                 border-radius: 6px;",
          icon("calculator"),
          h4("Advanced Analytics"),
          p("CPI and GCPI calculations for comprehensive trait evaluation")
        )
      )
    ),

    # 优势卡片
    div(
      class = "content-card",
      style = "background: white; 
             border-radius: 8px; 
             box-shadow: 0 2px 4px rgba(0,0,0,0.1); 
             padding: 30px; 
             margin-bottom: 40px;",
      h3("Why Choose Alpha-Cattle?", 
         style = "color: #1565c0; 
                border-bottom: 2px solid #1565c0; 
                padding-bottom: 10px;"),
      p("Your partner in modern cattle breeding, offering unparalleled insights into genomic selection.",
        style = "font-size: 16px; line-height: 1.6;"),
      
      div(
        style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin-top: 20px;",
        div(
          class = "feature-item",
          style = "text-align: center; padding: 20px;",
          icon("dna", style = "font-size: 2em; color: #1565c0;"),
          h4("Advanced Models"),
          p("State-of-the-art genomic selection models")
        ),
        div(
          class = "feature-item",
          style = "text-align: center; padding: 20px;",
          icon("upload", style = "font-size: 2em; color: #1565c0;"),
          h4("Easy Management"),
          p("Simplified data upload and management")
        ),
        div(
          class = "feature-item",
          style = "text-align: center; padding: 20px;",
          icon("lightbulb", style = "font-size: 2em; color: #1565c0;"),
          h4("Smart Insights"),
          p("Data-driven breeding decisions")
        )
      )
    ),

    # 联系方式卡片
    div(
      class = "content-card",
      style = "background: white; 
             border-radius: 8px; 
             box-shadow: 0 2px 4px rgba(0,0,0,0.1); 
             padding: 30px; 
             text-align: center;",
      h3("Get In Touch", 
         style = "color: #1565c0; 
                margin-bottom: 20px;"),
      div(
        style = "display: flex; justify-content: center; gap: 20px; flex-wrap: wrap;",
        div(
          icon("envelope", style = "color: #1565c0; margin-right: 10px;"),
          "info@alpha-cattle.com"
        ),
        div(
          icon("phone", style = "color: #1565c0; margin-right: 10px;"),
          "+1 (555) 123-4567"
        )
      )
    )
  )
)