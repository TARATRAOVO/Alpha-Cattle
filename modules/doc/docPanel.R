# 文档页面的UI文本
doc_text <- list(
  cn = list(
    subtitle = "整体性牛只基因组选择平台",
    description = "为农户和研究人员提供先进的基因组选择工具",
    project_overview = "项目概述",
    project_overview_desc = "Alpha-Cattle 整合了最先进的基因组选择技术和大规模牛群数据，为育种价值提供准确预测。",
    database_title = "综合数据库",
    database_desc = "3,800+ 参考群体的表型和全基因组信息",
    predictions_title = "准确预测",
    predictions_desc = "通过组合参考群体提高预测准确性",
    analytics_title = "高级分析",
    analytics_desc = "CPI 和 GCPI 计算，用于全面性状评估",
    why_choose = "为什么选择 Alpha-Cattle？",
    why_choose_desc = "您的现代化牛只育种伙伴，提供无与伦比的基因组选择见解。",
    advanced_models = "先进模型",
    advanced_models_desc = "最先进的基因组选择模型",
    easy_management = "便捷管理",
    easy_management_desc = "简化的数据上传和管理",
    smart_insights = "智能洞察",
    smart_insights_desc = "数据驱动的育种决策",
    contact_title = "联系我们"
  ),
  en = list(
    subtitle = "Holistic Cattle Genomic Selection Platform",
    description = "Empowering farmers and researchers with advanced genomic selection tools",
    project_overview = "Project Overview",
    project_overview_desc = "Alpha-Cattle integrates state-of-the-art genomic selection techniques with large-scale data from cattle populations to provide accurate predictions of breeding values.",
    database_title = "Comprehensive Database",
    database_desc = "3,800+ reference populations with phenotype and genome-wide information",
    predictions_title = "Accurate Predictions",
    predictions_desc = "Enhanced prediction accuracy through combined reference populations",
    analytics_title = "Advanced Analytics",
    analytics_desc = "CPI and GCPI calculations for comprehensive trait evaluation",
    why_choose = "Why Choose Alpha-Cattle?",
    why_choose_desc = "Your partner in modern cattle breeding, offering unparalleled insights into genomic selection.",
    advanced_models = "Advanced Models",
    advanced_models_desc = "State-of-the-art genomic selection models",
    easy_management = "Easy Management",
    easy_management_desc = "Simplified data upload and management",
    smart_insights = "Smart Insights",
    smart_insights_desc = "Data-driven breeding decisions",
    contact_title = "Get In Touch"
  )
)

docPanel <- function(language_state) {
  div(
    class = "doc-panel",
    style = "margin-top: 50px;",

    # Hero Section
    div(
      class = "hero-section",
      style = "background: linear-gradient(135deg, #1e88e5 0%, #1565c0 100%);
             color: white;
             padding: 60px 0;
             margin-bottom: 40px;",
      fluidRow(
        column(10,
          offset = 1,
          h1("Alpha-Cattle",
            style = "text-align: center;
                    font-size: 3.5em;
                    font-weight: 700;
                    margin-bottom: 20px;"
          ),
          h2(textOutput("docSubtitle"),
            style = "text-align: center;
                    font-weight: 300;
                    margin-bottom: 30px;"
          ),
          p(textOutput("docDescription"),
            style = "text-align: center;
                   font-size: 1.2em;
                   opacity: 0.9;"
          )
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
        h3(textOutput("projectOverviewTitle"),
          style = "color: #1565c0;
                  border-bottom: 2px solid #1565c0;
                  padding-bottom: 10px;"
        ),
        p(textOutput("projectOverviewDesc"),
          style = "font-size: 16px; line-height: 1.6;"
        ),
        div(
          style = "display: flex; flex-wrap: wrap; gap: 20px; margin-top: 20px;",
          div(
            style = "flex: 1; min-width: 250px;
                   background: #f5f5f5;
                   padding: 20px;
                   border-radius: 6px;",
            icon("database"),
            h4(textOutput("databaseTitle")),
            p(textOutput("databaseDesc"))
          ),
          div(
            style = "flex: 1; min-width: 250px;
                   background: #f5f5f5;
                   padding: 20px;
                   border-radius: 6px;",
            icon("chart-line"),
            h4(textOutput("predictionsTitle")),
            p(textOutput("predictionsDesc"))
          ),
          div(
            style = "flex: 1; min-width: 250px;
                   background: #f5f5f5;
                   padding: 20px;
                   border-radius: 6px;",
            icon("calculator"),
            h4(textOutput("analyticsTitle")),
            p(textOutput("analyticsDesc"))
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
        h3(textOutput("whyChooseTitle"),
          style = "color: #1565c0;
                  border-bottom: 2px solid #1565c0;
                  padding-bottom: 10px;"
        ),
        p(textOutput("whyChooseDesc"),
          style = "font-size: 16px; line-height: 1.6;"
        ),
        div(
          style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin-top: 20px;",
          div(
            class = "feature-item",
            style = "text-align: center; padding: 20px;",
            icon("dna", style = "font-size: 2em; color: #1565c0;"),
            h4(textOutput("advancedModelsTitle")),
            p(textOutput("advancedModelsDesc"))
          ),
          div(
            class = "feature-item",
            style = "text-align: center; padding: 20px;",
            icon("upload", style = "font-size: 2em; color: #1565c0;"),
            h4(textOutput("easyManagementTitle")),
            p(textOutput("easyManagementDesc"))
          ),
          div(
            class = "feature-item",
            style = "text-align: center; padding: 20px;",
            icon("lightbulb", style = "font-size: 2em; color: #1565c0;"),
            h4(textOutput("smartInsightsTitle")),
            p(textOutput("smartInsightsDesc"))
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
        h3(textOutput("contactTitle"),
          style = "color: #1565c0;
                  margin-bottom: 20px;"
        ),
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
}

# 文档页面的服务器逻辑
docServer <- function(input, output, session, language_state) {
  # 文档页面的输出
  output$docSubtitle <- renderText({
    doc_text[[language_state()]]$subtitle
  })

  output$docDescription <- renderText({
    doc_text[[language_state()]]$description
  })

  output$projectOverviewTitle <- renderText({
    doc_text[[language_state()]]$project_overview
  })

  output$projectOverviewDesc <- renderText({
    doc_text[[language_state()]]$project_overview_desc
  })

  output$databaseTitle <- renderText({
    doc_text[[language_state()]]$database_title
  })

  output$databaseDesc <- renderText({
    doc_text[[language_state()]]$database_desc
  })

  output$predictionsTitle <- renderText({
    doc_text[[language_state()]]$predictions_title
  })

  output$predictionsDesc <- renderText({
    doc_text[[language_state()]]$predictions_desc
  })

  output$analyticsTitle <- renderText({
    doc_text[[language_state()]]$analytics_title
  })

  output$analyticsDesc <- renderText({
    doc_text[[language_state()]]$analytics_desc
  })

  output$whyChooseTitle <- renderText({
    doc_text[[language_state()]]$why_choose
  })

  output$whyChooseDesc <- renderText({
    doc_text[[language_state()]]$why_choose_desc
  })

  output$advancedModelsTitle <- renderText({
    doc_text[[language_state()]]$advanced_models
  })

  output$advancedModelsDesc <- renderText({
    doc_text[[language_state()]]$advanced_models_desc
  })

  output$easyManagementTitle <- renderText({
    doc_text[[language_state()]]$easy_management
  })

  output$easyManagementDesc <- renderText({
    doc_text[[language_state()]]$easy_management_desc
  })

  output$smartInsightsTitle <- renderText({
    doc_text[[language_state()]]$smart_insights
  })

  output$smartInsightsDesc <- renderText({
    doc_text[[language_state()]]$smart_insights_desc
  })

  output$contactTitle <- renderText({
    doc_text[[language_state()]]$contact_title
  })
}
