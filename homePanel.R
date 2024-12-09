homePanel <- tabPanel(
  title = "Home",
  icon = icon("home"),

  # 全局样式定义
  tags$head(
    tags$style(HTML("
      :root {
        --primary-color: #2c3e50;
        --secondary-color: #6c757d;
        --bg-light: #f8f9fa;
        --spacing-lg: 80px;
        --spacing-md: 60px;
        --spacing-sm: 30px;
      }

      .container-fluid {
        max-width: 1200px;
        margin: 0 auto;
        padding: 0 20px;
      }

      .section {
        padding: var(--spacing-md) 0;
      }

      .feature-card {
        background: white;
        border-radius: 12px;
        padding: 25px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        transition: all 0.3s ease;
        height: 100%;
      }

      .feature-card:hover {
        transform: translateY(-5px);
        box-shadow: 0 6px 25px rgba(0,0,0,0.12);
      }

      .hero-title {
        font-size: 3.5rem;
        font-weight: 700;
        color: var(--primary-color);
        margin-bottom: 1.5rem;
      }

      .hero-subtitle {
        font-size: 1.25rem;
        color: var(--secondary-color);
        line-height: 1.6;
      }

      .btn-custom {
        padding: 0.75rem 2.5rem;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        border-radius: 30px;
        transition: all 0.3s ease;
      }

      .btn-custom:hover {
        transform: translateY(-2px);
        box-shadow: 0 5px 15px rgba(0,0,0,0.2);
      }

      .content-title {
        font-size: 2.25rem;
        font-weight: 700;
        color: var(--primary-color);
        margin-bottom: 2rem;
      }

      .footer {
        background: var(--primary-color);
        color: white;
        padding: var(--spacing-sm) 0;
      }

      /* 添加新的样式 */
      .image-container {
        position: relative;
        overflow: hidden;
        border-radius: 8px;
        margin: 30px 0;
        box-shadow: 0 4px 12px rgba(0,0,0,0.1);
      }
      
      .section-subtitle {
        font-size: 1.8em;
        color: #2c3e50;
        margin: 40px 0 20px 0;
        padding-bottom: 10px;
        border-bottom: 2px solid #eee;
      }
      
      .content-block {
        margin: 25px 0;
      }
      
      .feature-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
        gap: 30px;
        padding: 20px 0;
      }
      
      .scroll-reveal {
        opacity: 0;
        transform: translateY(20px);
        transition: all 0.6s ease-out;
      }
      
      .scroll-reveal.visible {
        opacity: 1;
        transform: translateY(0);
      }
      
      .bilingual-block {
        background: #f8f9fa;
        padding: 20px;
        border-radius: 8px;
        margin: 15px 0;
      }
      
      .highlight-text {
        background: #e3f2fd;
        padding: 2px 5px;
        border-radius: 3px;
      }
    "))
  ),

  # Hero Section
  div(class = "inlay", style = "height:70px;width:100%;background-color: white;"),
  div(
    class = "section hero-section",
    style = "background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);",
    div(
      class = "container-fluid",
      column(
        width = 8, offset = 2,
        div(
          class = "text-center",
          h1(class = "hero-title", "Alpha-Cattle"),
          p(
            class = "hero-subtitle",
            "Alpha-Cattle provides a platform and software to enable multiple modes of genomic selection and Prediction in Holstein Cattles."
          ),
          p(
            class = "hero-subtitle",
            "Alpha-Cattle 提供支持多模式基因组选择、荷斯坦奶牛遗传评估及性状预测的平台与软件。"
          ),
          div(
            style = "margin-top: 2rem",
            actionButton(
              inputId = "tabBut",
              label = "Learn More",
              class = "btn-custom btn-primary"
            )
          )
        )
      )
    )
  ),

  # Features Section
  div(
    class = "section features-section",
    div(
      class = "container-fluid",
      h2(class = "text-center content-title", "Our Features"),
      fluidRow(
        column(
          3,
          div(
            class = "feature-card",
            onclick = "Shiny.setInputValue('btn_phenotype', Math.random());",
            icon("chart-line", class = "fa-3x mb-4"),
            h4("Phenotype Processing", class = "mb-3"),
            p("Users can quickly filter their data and process phenotypic information.")
          )
        ),
        column(
          3,
          div(
            class = "feature-card",
            onclick = "Shiny.setInputValue('btn_genotype', Math.random());",
            icon("dna", class = "fa-3x mb-4"),
            h4("Genotype Processing", class = "mb-3"),
            p("Efficient processing and analysis of genotype data.")
          )
        ),
        column(
          3,
          div(
            class = "feature-card",
            onclick = "Shiny.setInputValue('btn_prediction', Math.random());",
            icon("brain", class = "fa-3x mb-4"),
            h4("Genomic Prediction", class = "mb-3"),
            p("Advanced algorithms for genomic prediction and selection.")
          )
        ),
        column(
          3,
          div(
            class = "feature-card",
            onclick = "Shiny.setInputValue('btn_analysis', Math.random());",
            icon("chart-bar", class = "fa-3x mb-4"),
            h4("Data Analysis", class = "mb-3"),
            p("Comprehensive tools for statistical analysis and visualization.")
          )
        )
      )
    )
  ),


  # Footer
  div(
    class = "footer",
    div(
      class = "container-fluid",
      column(
        width = 12,
        div(
          class = "text-center",
          h6("Copyright 2022 College of Agriculture and Biology, ShangHaiJiaoTong University.",
            style = "margin-bottom: 0.75rem;"
          ),
          h6("800# Dongchuan Road, Shanghai, China",
            style = "margin-bottom: 0.5rem;"
          ),
          p("E-mail: panyuchun1963@aliyun.com",
            style = "margin-bottom: 0;"
          )
        )
      )
    )
  )
)
