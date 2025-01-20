library(shinydashboard)
library(shiny)
library(rintrojs)
library(shinyBS)
library(plotly)
library(waiter)
# library(ggiraph)
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

# Source utility functions
source("modules/utils.R")

# Source module files
source("modules/home/homePanel.R")
source("modules/impute/imputePanel.R")
source("modules/phe/phePanel.R")
source("modules/gs/gsPanel.R")
source("modules/cpi/cpiPanel.R")
source("modules/jgs/jgsPanel.R")
source("modules/jcpi/jcpiPanel.R")
source("modules/pdata/pdataPanel.R")
source("modules/doc/docPanel.R")
source("modules/meg/megPanel.R")

options(shiny.maxRequestSize = 10 * 1024^3)
dir.create("./temp", showWarnings = F)
task <<- sample(1000:9999, 1)
path_gs <<- paste0("temp/", task, "_gs")
path_phe <<- paste0("temp/", task, "_phe")
path_impute <<- paste0("temp/", task, "_impute")
path_jgs <<- paste0("temp/", task, "_jgs")
path_cpi <<- paste0("temp/", task, "_cpi")
path_jcpi <<- paste0("temp/", task, "_jcpi")

# UI text definitions
ui_text <- list(
  cn = list(
    features = "功能特点",
    phenotype = "表型处理",
    phenotype_desc = "用户可以快速过滤数据并处理表型信息。",
    genotype = "基因型处理",
    genotype_desc = "高效的基因型数据处理与分析。",
    prediction = "基因组预测",
    prediction_desc = "先进的基因组预测和选择算法。",
    analysis = "数据分析",
    analysis_desc = "全面的统计分析和可视化工具。",
    learn_more = "了解更多",
    home = "首页",
    genomic_selection = "基因组选择",
    joint_genomic_selection = "联合基因组选择",
    phenotype_data = "表型数据",
    documentation = "文档",
    major_effect_gene = "主效基因"
  ),
  en = list(
    features = "Our Features",
    phenotype = "Phenotype Processing",
    phenotype_desc = "Users can quickly filter their data and process phenotypic information.",
    genotype = "Genotype Processing",
    genotype_desc = "Efficient processing and analysis of genotype data.",
    prediction = "Genomic Prediction",
    prediction_desc = "Advanced algorithms for genomic prediction and selection.",
    analysis = "Data Analysis",
    analysis_desc = "Comprehensive tools for statistical analysis and visualization.",
    learn_more = "Learn More",
    home = "Home",
    genomic_selection = "Genomic Selection",
    joint_genomic_selection = "Joint Genomic Selection",
    phenotype_data = "Phenotype Data",
    documentation = "Documentation",
    major_effect_gene = "Major Effect Gene"
  )
)

mod_comp_box <<- function(id, label, title, content, placement = "right", choices = colnames(demo_gs_phe_dat), single = FALSE) {
  if (!single) {
    box <- column(
      width = 2,
      div(
        id = paste0(id, "_BS"),
        multiInput(
          inputId = id,
          label = span(id = paste0(id, "_BS"), label),
          choices = choices,
          options = list(
            enable_search = TRUE, non_selected_header = "Choose between:",
            selected_header = "You have selected:"
          ),
          width = "100%"
        )
      ),
      bsPopover(
        id = paste0(id, "_BS"),
        title = title,
        content = content,
        placement = placement, trigger = "hover", options = list(container = "body")
      )
    )
  } else {
    box <- column(
      width = 2,
      div(
        id = paste0(id, "_BS"),
        box(
          width = 12,
          title = "",
          status = "danger",
          solidHeader = F,
          multiInput(
            inputId = id,
            label = span(id = paste0(id, "_BS"), label),
            choices = choices,
            options = list(
              enable_search = TRUE, non_selected_header = "Choose between:",
              selected_header = "You have selected:", limit = 1
            ),
            width = "100%"
          )
        )
      ),
      bsPopover(
        id = paste0(id, "_BS"),
        title = title,
        content = content,
        placement = placement, trigger = "hover", options = list(container = "body")
      )
    )
  }
  return(box)
}

mod_comp_box1 <<- function(id, label, title, content, placement = "right", choices = colnames(demo_trait_type1_dat)[c(1, 2, 3, 4, 5, 7, 8, 9, 13)], single = FALSE) {
  if (!single) {
    box <- column(
      width = 2,
      div(
        id = paste0(id, "_BS"),
        multiInput(
          inputId = id,
          label = span(id = paste0(id, "_BS"), label),
          choices = choices,
          options = list(
            enable_search = TRUE, non_selected_header = "Choose between:",
            selected_header = "You have selected:"
          ),
          width = "100%"
        )
      ),
      bsPopover(
        id = paste0(id, "_BS"),
        title = title,
        content = content,
        placement = placement, trigger = "hover", options = list(container = "body")
      )
    )
  } else {
    box <- column(
      width = 2,
      div(
        id = paste0(id, "_BS"),
        box(
          width = 12,
          title = "",
          status = "danger",
          solidHeader = F,
          multiInput(
            inputId = id,
            label = span(id = paste0(id, "_BS"), label),
            choices = choices,
            options = list(
              enable_search = TRUE, non_selected_header = "Choose between:",
              selected_header = "You have selected:", limit = 1
            ),
            width = "100%"
          )
        )
      ),
      bsPopover(
        id = paste0(id, "_BS"),
        title = title,
        content = content,
        placement = placement, trigger = "hover", options = list(container = "body")
      )
    )
  }
  return(box)
}


mod_comp_box2 <<- function(id, label, title, content, placement = "right", choices = colnames(demo_trait_type1_dat)[c(1, 2, 3, 4, 6, 10, 11, 12)], single = FALSE) {
  if (!single) {
    box <- column(
      width = 2,
      div(
        id = paste0(id, "_BS"),
        multiInput(
          inputId = id,
          label = span(id = paste0(id, "_BS"), label),
          choices = choices,
          options = list(
            enable_search = TRUE, non_selected_header = "Choose between:",
            selected_header = "You have selected:"
          ),
          width = "100%"
        )
      ),
      bsPopover(
        id = paste0(id, "_BS"),
        title = title,
        content = content,
        placement = placement, trigger = "hover", options = list(container = "body")
      )
    )
  } else {
    box <- column(
      width = 2,
      div(
        id = paste0(id, "_BS"),
        box(
          width = 12,
          title = "",
          status = "danger",
          solidHeader = F,
          multiInput(
            inputId = id,
            label = span(id = paste0(id, "_BS"), label),
            choices = choices,
            options = list(
              enable_search = TRUE, non_selected_header = "Choose between:",
              selected_header = "You have selected:", limit = 1
            ),
            width = "100%"
          )
        )
      ),
      bsPopover(
        id = paste0(id, "_BS"),
        title = title,
        content = content,
        placement = placement, trigger = "hover", options = list(container = "body")
      )
    )
  }
  return(box)
}

isValidEmail <<- function(x) {
  grepl("\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>", as.character(x), ignore.case = TRUE)
}

submit_confirm <<- function(id) {
  modalDialog(
    "Are you sure to submit (please check your data and model once more)? If you click \"OK\", the page will be refreshed",
    title = "Submit?",
    footer = tagList(
      actionButton("cancel", "Cancel"),
      actionButton(paste0(id, "_ok"), "OK", class = "btn btn-danger")
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


gsMenu <- navbarMenu(
  title = uiOutput("gsMenuTitle"),
  icon = icon("dna"),
  imputePanel,
  phePanel,
  gsPanel,
  cpiPanel
)

# Source module files with correct paths
source("modules/jgs/jgsPanel.R")
source("modules/jcpi/jcpiPanel.R")
source("modules/pdata/pdataPanel.R")
source("modules/doc/docPanel.R")
source("modules/meg/megPanel.R")

jgsMenu <- navbarMenu(
  title = uiOutput("jgsMenuTitle"),
  icon = icon("mixcloud"),
  jgsPanel,
  jcpiPanel
)

# 添加全局变量用于存储语言状态
language_state <- reactiveVal("cn")

ui <- tagList(
  # Language Switch Button
  div(
    class = "language-switch",
    actionButton(
      inputId = "switchLanguage",
      label = "中文/EN",
      class = "language-btn"
    )
  ),
  navbarPage(
    id = "tabs",
    title = span(
      "Alpha-Cattle",
      style = "font-size:30px;margin-right:50px"
    ),
    position = "fixed-top",
    header = tagList(
      use_waiter(),
      useShinyjs(),
      tags$head(
        tags$style(HTML("
          .navbar-fixed-top{
            display: flex;
          }
          .navbar-brand {
            FONT-VARIANT: JIS04;
            line-height: 0.6;
          }

          /* 调整导航栏标题和图标的样式 */
          .navbar-nav > li > a {
            display: flex !important;
            align-items: center !important;
            gap: 5px !important;
          }

          /* 确保图标和文字在同一行 */
          .navbar-nav .fa {
            margin-right: 5px !important;
            display: inline-block !important;
            vertical-align: middle !important;
          }

          /* 调整下拉菜单的样式 */
          .dropdown-menu {
            margin-top: 0 !important;
          }

          /* 确保标题文字和图标垂直居中 */
          .shiny-html-output {
            display: inline-block !important;
            vertical-align: middle !important;
          }

          /* 调整语言切换按钮的位置 */
          .language-switch {
            position: fixed !important;  /* 改为 fixed */
            top: 60px !important;
            right: 20px;
            z-index: 1001 !important;  /* 确保按钮始终在最上层 */
          }

          /* 美化语言切换按钮 */
          .language-btn {
            background: white;
            border: 1px solid #ddd;
            padding: 5px 15px;
            border-radius: 20px;
            cursor: pointer;
            transition: all 0.3s ease;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            color: #2c3e50;
            font-weight: 500;
          }

          .language-btn:hover {
            background: #f8f9fa;
            box-shadow: 0 3px 8px rgba(0,0,0,0.15);
            transform: translateY(-1px);
            color: #1a252f;
          }
        ")),
        tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
        tags$style("@import url(https://use.fontawesome.com/releases/v6.1.1/css/all.css);")
      )
    ),
    theme = shinythemes::shinytheme(theme = "flatly"),
    tabPanel(
      title = uiOutput("homeTitle"),
      icon = icon("home"),
      homePanel
    ),
    gsMenu,
    jgsMenu,
    tabPanel(
      title = uiOutput("pdataTitle"),
      icon = icon("table"),
      pdataPanel
    ),
    tabPanel(
      title = uiOutput("megTitle"),
      icon = icon("dna"),
      megPanel
    ),
    tabPanel(
      title = uiOutput("docTitle"),
      icon = icon("book"),
      docPanel(language_state)
    )
  )
)

server <- function(input, output, session) {
  # 添加新的 UI 输出
  output$homeTitle <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$home
  })

  output$gsMenuTitle <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$genomic_selection
  })

  output$jgsMenuTitle <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$joint_genomic_selection
  })

  output$pdataTitle <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$phenotype_data
  })

  output$docTitle <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$documentation
  })

  # 添加语言切换处理
  observeEvent(input$switchLanguage, {
    current_lang <- language_state()
    language_state(ifelse(current_lang == "cn", "en", "cn"))
  })

  # 更新UI文本
  output$heroSubtitle <- renderUI({
    lang <- language_state()
    if (lang == "cn") {
      p(
        class = "hero-subtitle",
        "Alpha-Cattle 提供支持多模式基因组选择、荷斯坦奶牛遗传评估及性状预测的平台与软件。"
      )
    } else {
      p(
        class = "hero-subtitle",
        "Alpha-Cattle provides a platform and software to enable multiple modes of genomic selection and Prediction in Holstein Cattles."
      )
    }
  })

  output$learnMoreText <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$learn_more
  })

  output$featuresTitle <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$features
  })

  output$phenotypeTitle <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$phenotype
  })

  output$phenotypeDesc <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$phenotype_desc
  })

  output$genotypeTitle <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$genotype
  })

  output$genotypeDesc <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$genotype_desc
  })

  output$predictionTitle <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$prediction
  })

  output$predictionDesc <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$prediction_desc
  })

  output$analysisTitle <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$analysis
  })

  output$analysisDesc <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$analysis_desc
  })

  # 调用文档页面的服务器逻辑
  docServer(input, output, session, language_state)

  source("modules/jgs/jgsServer.R", local = TRUE)
  source("modules/gs/gsServer.R", local = TRUE)
  source("modules/phe/pheServer.R", local = TRUE)
  source("modules/impute/imputeServer.R", local = TRUE)
  source("modules/cpi/cpiServer.R", local = TRUE)
  source("modules/jcpi/jcpiServer.R", local = TRUE)
  source("modules/pdata/pdataServer.R", local = TRUE)
  source("modules/meg/megServer.R", local = TRUE)

  server_gs(input, output, session)
  server_jgs(input, output, session)
  pheServer(input, output, session)
  pdataServer(input, output, session)
  imputeServer(input, output, session)
  server_cpi(input, output, session)
  server_jcpi(input, output, session)
  magServer(input, output, session)
  output$pic <- renderImage(
    {
      list(
        src = "19.jpg",
        width = "70%"
      )
    },
    deleteFile = FALSE
  )
  output$pic_pca <- renderImage(
    {
      list(
        src = "pca_accuracy.jpg",
        width = "90%"
      )
    },
    deleteFile = FALSE
  )
  output$pic_cor <- renderImage(
    {
      list(
        src = "cor_he.jpg",
        width = "90%"
      )
    },
    deleteFile = FALSE
  )
  output$pic_joint <- renderImage(
    {
      list(
        src = "joint.jpg",
        width = "90%"
      )
    },
    deleteFile = FALSE
  )
  observeEvent(input$tabBut, {
    # browser()
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

  # 添加主效基因标题输出
  output$megTitle <- renderUI({
    lang <- language_state()
    ui_text[[lang]]$major_effect_gene
  })
}


shinyApp(ui, server)
