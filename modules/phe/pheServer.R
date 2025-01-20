pheServer <- function(input, output, session) {
  flag <- 0
  phe_name <- reactive({
    if (is.null(input$phe_colname)) {
      names(phe_dat())
    } else {
      input$phe_colname
    }
  })

  observeEvent(input$phe_upload, {
    updatePickerInput(session = session, inputId = "phe_colname", choices = phe_name())
  })

  phe_dat <- reactive({
    req(input$phe_upload)

    ext <- tools::file_ext(input$phe_upload$name)
    dat <- switch(ext,
      csv = data.table::fread(input$phe_upload$datapath, encoding = "UTF-8"),
      txt = data.table::fread(input$phe_upload$datapath, encoding = "UTF-8"),
      validate("[ERROR:] Invalid phenotype file; Please upload a .csv or .txt file")
    )
    dir.create(path_phe, showWarnings = F)
    system(paste0("cp ", input$phe_upload$datapath, " ", path_phe, "/", input$phe_upload$name))

    if (!is.null(input$phe_colname)) {
      if (input$trans0 == "Yes") {
        dat[dat == 0] <- NA
      }
      if ("ID" %in% input$phe_colname) {
        name_filter <- input$phe_colname[-which(input$phe_colname == "ID")]
      } else {
        name_filter <- input$phe_colname
      }
      test <- dat %>%
        summarise(across(
          .cols = name_filter,
          .fns = list(mean = mean, sd = sd), na.rm = TRUE
        ))
      for (i in name_filter) {
        a <- select(test, paste0(i, "_mean")) - input$sd * select(test, paste0(i, "_sd"))
        b <- select(test, paste0(i, "_mean")) + input$sd * select(test, paste0(i, "_sd"))
        dat[[i]][(select(dat, i) < a[[1]] | select(dat, i) > b[[1]])] <- NA
      }

      if (input$FilterNA == "Yes") {
        dat <- dat %>% filter(!across(name_filter, .fns = is.na))
      }
    }

    return(dat)
  })

  # 更新散点图的X和Y轴选择器
  observe({
    req(phe_dat())
    num_cols <- names(phe_dat())[sapply(phe_dat(), is.numeric)]
    updatePickerInput(session, "scatter_x", choices = num_cols)
    updatePickerInput(session, "scatter_y", choices = num_cols)
  })

  # 动态生成Plot UI
  output$plot_ui <- renderUI({
    req(input$plot_type)

    if (input$plot_type == "correlation") {
      plotlyOutput("correlation_plot", height = "800px")
    } else if (input$plot_type == "scatter") {
      plotlyOutput("scatter_plot", height = "600px")
    } else {
      plotOutput("main_plot", height = "600px")
    }
  })

  # 主绘图输出
  output$main_plot <- renderPlot({
    req(phe_dat(), input$plot_type, input$phe_colname)

    # 确保数据是数据框格式
    plot_data <- as.data.frame(phe_dat())
    # 只选择用户选择的列
    plot_data <- plot_data[, input$phe_colname, drop = FALSE]
    # 只保留数值型列
    numeric_cols <- sapply(plot_data, is.numeric)
    plot_data <- plot_data[, numeric_cols, drop = FALSE]

    if (ncol(plot_data) == 0) {
      return(NULL)
    }

    if (input$plot_type == "histogram") {
      if (ncol(plot_data) > 1) {
        # 多个变量的情况
        plot_list <- lapply(names(plot_data), function(col) {
          ggplot(plot_data, aes_string(x = col)) +
            geom_histogram(fill = "#69b3a2", color = "#e9ecef", bins = 30, alpha = 0.9) +
            theme_minimal() +
            labs(title = col, x = col, y = "Frequency") +
            theme(plot.title = element_text(hjust = 0.5))
        })
        grid.arrange(grobs = plot_list, ncol = min(2, length(plot_list)))
      } else {
        # 单个变量的情况
        ggplot(plot_data, aes_string(x = names(plot_data)[1])) +
          geom_histogram(fill = "#69b3a2", color = "#e9ecef", bins = 30, alpha = 0.9) +
          theme_minimal() +
          labs(title = names(plot_data)[1], x = names(plot_data)[1], y = "Frequency") +
          theme(plot.title = element_text(hjust = 0.5))
      }
    } else if (input$plot_type == "boxplot") {
      # 转换数据为长格式
      plot_data_long <- tidyr::pivot_longer(plot_data,
        cols = everything(),
        names_to = "variable",
        values_to = "value"
      )
      ggplot(plot_data_long, aes(x = variable, y = value)) +
        geom_boxplot(fill = "#69b3a2", alpha = 0.9) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Variable", y = "Value", title = "Box Plot") +
        theme(plot.title = element_text(hjust = 0.5))
    } else if (input$plot_type == "density") {
      if (ncol(plot_data) > 1) {
        # 多个变量的情况
        plot_list <- lapply(names(plot_data), function(col) {
          ggplot(plot_data, aes_string(x = col)) +
            geom_density(fill = "#69b3a2", alpha = 0.9) +
            theme_minimal() +
            labs(title = col, x = col, y = "Density") +
            theme(plot.title = element_text(hjust = 0.5))
        })
        grid.arrange(grobs = plot_list, ncol = min(2, length(plot_list)))
      } else {
        # 单个变量的情况
        ggplot(plot_data, aes_string(x = names(plot_data)[1])) +
          geom_density(fill = "#69b3a2", alpha = 0.9) +
          theme_minimal() +
          labs(title = names(plot_data)[1], x = names(plot_data)[1], y = "Density") +
          theme(plot.title = element_text(hjust = 0.5))
      }
    } else if (input$plot_type == "violin") {
      # 转换数据为长格式
      plot_data_long <- tidyr::pivot_longer(plot_data,
        cols = everything(),
        names_to = "variable",
        values_to = "value"
      )
      ggplot(plot_data_long, aes(x = variable, y = value)) +
        geom_violin(fill = "#69b3a2", alpha = 0.9) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Variable", y = "Value", title = "Violin Plot") +
        theme(plot.title = element_text(hjust = 0.5))
    }
  })

  # 散点图
  output$scatter_plot <- renderPlotly({
    req(phe_dat(), input$scatter_x, input$scatter_y)

    plot_ly(phe_dat(),
      x = as.formula(paste0("~", input$scatter_x)),
      y = as.formula(paste0("~", input$scatter_y)),
      type = "scatter",
      mode = "markers",
      marker = list(color = "#69b3a2", size = 8, opacity = 0.7)
    ) %>%
      layout(
        title = paste(input$scatter_y, "vs", input$scatter_x),
        xaxis = list(title = input$scatter_x),
        yaxis = list(title = input$scatter_y)
      )
  })

  # 相关性热图
  output$correlation_plot <- renderPlotly({
    req(phe_dat(), input$phe_colname)

    # 确保数据是数据框格式
    plot_data <- as.data.frame(phe_dat())

    # 只选择用户选择的列
    if (length(input$phe_colname) > 0) {
      plot_data <- plot_data[, input$phe_colname, drop = FALSE]
    }

    # 只保留数值型列
    numeric_cols <- sapply(plot_data, is.numeric)
    plot_data <- plot_data[, numeric_cols, drop = FALSE]

    # 检查数值型变量的数量
    if (ncol(plot_data) < 2) {
      # 返回一个带有提示信息的空白图
      return(
        plot_ly() %>%
          add_annotations(
            text = "Please select at least two numeric variables to generate correlation heatmap",
            showarrow = FALSE,
            font = list(size = 16)
          )
      )
    }

    # 计算相关系数矩阵
    cor_matrix <- cor(plot_data, use = "pairwise.complete.obs")

    # 创建热图
    plot_ly(
      x = colnames(cor_matrix),
      y = colnames(cor_matrix),
      z = cor_matrix,
      type = "heatmap",
      colors = colorRamp(c("#4575b4", "white", "#d73027")),
      text = round(cor_matrix, 2),
      hoverongaps = FALSE
    ) %>%
      layout(
        title = "Correlation Heatmap",
        xaxis = list(title = ""),
        yaxis = list(title = "")
      )
  })

  output$phe_dat <- DT::renderDataTable({
    select(phe_dat(), phe_name())
  })

  output$phe_summary <- DT::renderDataTable({
    phe_summary_num <- map_dbl(phe_name(), function(x) sum(!is.na(phe_dat()[[x]])))
    phe_summary_mean <- map_dbl(phe_name(), function(x) mean(phe_dat()[[x]], na.rm = T))
    phe_summary_sd <- map_dbl(phe_name(), function(x) sd(phe_dat()[[x]], na.rm = T))

    phe <- data.frame(
      Trait = phe_name(),
      Size = phe_summary_num,
      Mean = phe_summary_mean,
      SD = phe_summary_sd
    )
    phe
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$phe_upload$name), "_filter.", tools::file_ext(input$phe_upload$name))
    },
    content = function(file) {
      switch(tools::file_ext(input$phe_upload$name),
        csv = write.csv(phe_dat(), file, row.names = FALSE, quote = F, na = ""),
        txt = write.table(phe_dat(), file, row.names = FALSE, quote = F, na = "")
      )
    }
  )

  observeEvent(input$start, {
    browser()
  })
}
