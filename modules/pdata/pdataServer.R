library(ggplot2)
library(plotly)
library(gridExtra)
library(dplyr)
library(tidyr)

pdataServer <- function(input, output, session) {
  flag <- 0

  # 获取用户选择的表型列名，如果为空，则使用数据的列名（排除ID列）
  pdata_name <- reactive({
    req(pdata_dat())
    all_names <- names(pdata_dat())
    # 排除ID列
    all_names <- all_names[all_names != "ID"]
    if (is.null(input$pdata_colname)) {
      all_names
    } else {
      input$pdata_colname
    }
  })

  # 当用户选择服务器文件时更新表型列选择器
  observeEvent(input$pdata_server_file, {
    updatePickerInput(session = session, inputId = "pdata_colname", choices = pdata_name())
  })

  # 处理选中的表型数据
  pdata_dat <- reactive({
    req(input$pdata_server_file)

    file_path <- paste0("./data_base/pheno_data/", input$pdata_server_file)
    ext <- tools::file_ext(file_path)
    dat <- switch(ext,
      csv = data.table::fread(file_path, encoding = "UTF-8"),
      txt = data.table::fread(file_path, encoding = "UTF-8"),
      validate("[ERROR:] Invalid phenotype file. Please select a .csv or .txt file")
    )

    # 处理 0 转换为 NA，过滤 NA 和标准差筛选逻辑
    if (!is.null(input$pdata_colname)) {
      if (input$pdata_trans0 == "Yes") {
        dat[dat == 0] <- NA
      }

      # 对选中的列进行标准差筛选
      name_filter <- input$pdata_colname
      test <- dat[, lapply(.SD, function(x) list(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))), .SDcols = name_filter]
      for (i in name_filter) {
        a <- test[[i]]$mean - input$pdata_sd * test[[i]]$sd
        b <- test[[i]]$mean + input$pdata_sd * test[[i]]$sd
        dat[[i]][(dat[[i]] < a | dat[[i]] > b)] <- NA
      }

      if (input$pdata_FilterNA == "Yes") {
        dat <- dat[!Reduce(`|`, lapply(name_filter, function(x) is.na(dat[[x]])))]
      }
    }

    return(dat)
  })

  # 数据表渲染 - Limit to first 50 rows
  output$pdata_dat <- DT::renderDataTable({
    dat <- as.data.frame(pdata_dat()[, pdata_name(), drop = FALSE])
    head(dat, 50)
  })

  # 摘要表渲染
  output$pdata_summary <- DT::renderDataTable({
    pdata_summary_num <- map_dbl(pdata_name(), function(x) sum(!is.na(pdata_dat()[[x]])))
    pdata_summary_mean <- map_dbl(pdata_name(), function(x) mean(pdata_dat()[[x]], na.rm = TRUE))
    pdata_summary_sd <- map_dbl(pdata_name(), function(x) sd(pdata_dat()[[x]], na.rm = TRUE))

    pdata <- data.frame(
      Trait = pdata_name(),
      Size = pdata_summary_num,
      Mean = pdata_summary_mean,
      SD = pdata_summary_sd
    )
    pdata
  })

  # 直方图渲染
  output$pdata_plot <- renderCombineWidgets({
    plot_dat <- as.data.frame(pdata_dat()[, pdata_name(), drop = FALSE]) %>%
      select_if(is.numeric)
    plist <- lapply(names(plot_dat), function(x) {
      ggplotly(ggplot(plot_dat, aes_string(x)) +
        geom_histogram(
          fill = "#69b3a2", bins = 30,
          color = "#e9ecef", alpha = 0.9
        ) +
        theme_bw() +
        labs(x = x, y = "Counts"))
    })

    # 使用 combineWidgets 组合多个图表
    p <- combineWidgets(list = plist, nrow = ceiling(length(plist) / 3))
    p
  })

  # 更新散点图的X和Y轴选择器
  observe({
    req(pdata_dat())
    num_cols <- names(pdata_dat())[sapply(pdata_dat(), is.numeric)]
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
    req(pdata_dat(), input$plot_type, input$pdata_colname)

    # 确保数据是数据框格式
    plot_data <- as.data.frame(pdata_dat())
    # 只选择用户选择的列
    plot_data <- plot_data[, input$pdata_colname, drop = FALSE]
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
    req(pdata_dat(), input$scatter_x, input$scatter_y)

    plot_ly(pdata_dat(),
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
    req(pdata_dat(), input$pdata_colname)

    # 确保数据是数据框格式
    plot_data <- as.data.frame(pdata_dat())

    # 只选择用户选择的列
    if (length(input$pdata_colname) > 0) {
      plot_data <- plot_data[, input$pdata_colname, drop = FALSE]
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
          ) %>%
          layout(
            title = "Note",
            xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
            yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
          )
      )
    }

    # 计算相关性矩阵
    cor_matrix <- try(
      {
        cor(plot_data, use = "pairwise.complete.obs")
      },
      silent = TRUE
    )

    if (inherits(cor_matrix, "try-error")) {
      return(
        plot_ly() %>%
          add_annotations(
            text = "Error calculating correlation. Please check if the data contains sufficient non-missing values",
            showarrow = FALSE,
            font = list(size = 16)
          ) %>%
          layout(
            title = "Error",
            xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
            yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
          )
      )
    }

    # 创建热图
    plot_ly(
      x = colnames(cor_matrix),
      y = colnames(cor_matrix),
      z = as.matrix(cor_matrix),
      type = "heatmap",
      colors = colorRamp(c("#4575b4", "white", "#d73027"))
    ) %>%
      add_annotations(
        text = round(cor_matrix, 2),
        show = TRUE,
        font = list(color = "black"),
        showarrow = FALSE
      ) %>%
      layout(
        title = list(
          text = "Correlation Heatmap",
          y = 0.95
        ),
        xaxis = list(
          title = "",
          tickangle = 45
        ),
        yaxis = list(
          title = ""
        ),
        margin = list(t = 100, b = 100)
      )
  })

  # 文件下载处理
  # output$downloadPData <- downloadHandler(
  #   filename = function() {
  #     paste0(tools::file_path_sans_ext(input$pdata_server_file), "_filter.", tools::file_ext(input$pdata_server_file))
  #   },
  #   content = function(file) {
  #     switch(tools::file_ext(input$pdata_server_file),
  #            csv = write.csv(pdata_dat(), file, row.names = FALSE, quote = FALSE, na = ""),
  #            txt = write.table(pdata_dat(), file, row.names = FALSE, quote = FALSE, na = "")
  #     )
  #   }
  # )
}
