pdataServer <- function(input, output, session) {
  flag <- 0
  
  # 获取输入的列名，如果未指定列名，则使用数据中的列名
  phe_name <- reactive({
    if (is.null(input$phe_colname)) {
      names(phe_dat())
    } else {
      input$phe_colname
    }
  })
  
  # 文件库中数据的列名来更新选择框
  updatePickerInput(session = session, inputId = "phe_colname", choices = names(phe_dat()))

  # 读取上传的表型数据，并根据用户的选择对数据进行筛选和过滤
  phe_dat <- reactive({
    # 读取用于分析的文件
    dat <- data.table::fread("tdm_growth.txt", encoding = "UTF-8")
    # 创建目录并将上传的文件复制到该目录
    dir.create(path_phe, showWarnings = FALSE)
    system(paste0("cp ", input$phe_upload$datapath, " ", path_phe, "/", input$phe_upload$name))
    
    # 如果用户指定了列名，对数据进行过滤处理
    if (!is.null(input$phe_colname)) {
      # 如果选择了将0转换为NA
      if (input$trans0 == "Yes") {
        dat[dat == 0] <- NA
      }
      
      # 如果ID在列名中，将ID列去除
      if ("ID" %in% input$phe_colname) {
        name_filter <- input$phe_colname[-which(input$phe_colname == "ID")]
      } else {
        name_filter <- input$phe_colname
      }
      
      # 计算选择列的均值和标准差
      test <- dat %>%
        summarise(across(
          .cols = name_filter, 
          .fns = list(mean = mean, sd = sd), na.rm = TRUE
        ))
      
      # 根据用户指定的标准差范围过滤数据
      for (i in name_filter) {
        a <- select(test, paste0(i, "_mean")) - input$sd * select(test, paste0(i, "_sd"))
        b <- select(test, paste0(i, "_mean")) + input$sd * select(test, paste0(i, "_sd"))
        dat[[i]][(select(dat, i) < a[[1]] | select(dat, i) > b[[1]])] <- NA
      }
      
      # 如果选择了过滤掉NA值
      if (input$FilterNA == "Yes") {
        dat <- dat %>% filter(!across(name_filter, .fns = is.na))
      }
    }
    
    return(dat)
  })
  
  # 表格显示上传数据
  output$phe_dat <- DT::renderDataTable({
    select(phe_dat(), phe_name())
  })

  # 显示数据的摘要信息（数量、均值、标准差）
  output$phe_summary <- DT::renderDataTable({
    phe_summary_num <- map_dbl(phe_name(), function(x) sum(!is.na(phe_dat()[[x]])))
    phe_summary_mean <- map_dbl(phe_name(), function(x) mean(phe_dat()[[x]], na.rm = TRUE))
    phe_summary_sd <- map_dbl(phe_name(), function(x) sd(phe_dat()[[x]], na.rm = TRUE))
    
    phe <- data.frame(
      Trait = phe_name(),
      Size = phe_summary_num,
      Mean = phe_summary_mean,
      SD = phe_summary_sd
    )
    
    phe
  })
  
  # 绘制选择列的直方图
  output$plot <- renderCombineWidgets({
    plot_dat <- phe_dat() %>%
      select(phe_name()) %>%
      select_if(is.numeric)
    
    # 为每个数值列绘制直方图
    plist <- lapply(names(plot_dat), function(x) {
      ggplotly(ggplot(plot_dat, aes_string(x)) +
        geom_histogram(fill = "#69b3a2", bins = 30, color = "#e9ecef", alpha = 0.9) +
        theme_bw() +
        labs(x = x, y = "Counts"))
    })
    
    # 将多个图表组合在一起显示
    p <- combineWidgets(list = plist, nrow = ceiling(length(plist) / 3))
    p
  })
  
  # 文件下载处理器，生成过滤后的文件
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$phe_upload$name), "_filter.", tools::file_ext(input$phe_upload$name))
    },
    content = function(file) {
      switch(tools::file_ext(input$phe_upload$name),
             csv = write.csv(phe_dat(), file, row.names = FALSE, quote = FALSE, na = ""),
             txt = write.table(phe_dat(), file, row.names = FALSE, quote = FALSE, na = "")
      )
    }
  )

  # 调试模式
  observeEvent(input$start, {
    browser()
  })
}