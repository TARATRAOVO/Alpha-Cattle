pheServer<-function(input,output,session){
  flag<-0
  phe_name<-reactive({
    if(is.null(input$phe_colname)){
    names(phe_dat())}else{
      input$phe_colname
    }
    })
  
  observeEvent(input$phe_upload,{
    updatePickerInput(session = session,inputId="phe_colname",choices=phe_name())
  })


  phe_dat<-reactive({
    req(input$phe_upload)

    ext<-tools::file_ext(input$phe_upload$name)
    dat<-switch(ext,
                csv = data.table::fread(input$phe_upload$datapath,encoding = "UTF-8"),
                txt = data.table::fread(input$phe_upload$datapath,encoding = "UTF-8"),
                validate("[ERROR:] Invalid phenotype file; Please upload a .csv or .txt file")
                
    )
    dir.create(path_phe,showWarnings = F)
    system(paste0("cp ",input$phe_upload$datapath," ",path_phe,"/",input$phe_upload$name))

    if(!is.null(input$phe_colname)){
      if(input$trans0 == "Yes"){
        dat[dat==0]<-NA
      }
      if("ID" %in% input$phe_colname){
       name_filter<-input$phe_colname[-which(input$phe_colname=="ID")]
      }else{
       name_filter<-input$phe_colname
      }
      test<-dat%>%
        summarise(across(
          .cols = name_filter, 
          .fns = list(mean = mean, sd = sd), na.rm = TRUE 
         
        ))
    for(i in name_filter){
     a<-select(test,paste0(i,"_mean")) - input$sd * select(test,paste0(i,"_sd"))
      b<-select(test,paste0(i,"_mean")) + input$sd * select(test,paste0(i,"_sd"))
      dat[[i]][(select(dat,i)<a[[1]] |select(dat,i)>b[[1]] )]<-NA
      
    }

      if(input$FilterNA =="Yes"){
        dat<- dat %>% filter(!across(name_filter, .fns = is.na))
        }
      }
   

    

    return(dat)
  })
  ###�ϴ��ļ�֮����������ͼƬ��table
  ##table

  output$phe_dat<-DT::renderDataTable({
    select(phe_dat(),phe_name())
  })

  ##summary

  output$phe_summary<-DT::renderDataTable({
    
    phe_summary_num<-map_dbl(phe_name(),function(x)sum(!is.na(phe_dat()[[x]])))
    phe_summary_mean<-map_dbl(phe_name(),function(x)mean(phe_dat()[[x]],na.rm=T))
    phe_summary_sd<-map_dbl(phe_name(),function(x)sd(phe_dat()[[x]],na.rm=T))
    
    phe<-data.frame(
      Trait=phe_name(),
      Size=phe_summary_num,
      Mean=phe_summary_mean,
      SD=phe_summary_sd)
    phe
  })
  ###����numeric���л�����

  output$plot<-renderCombineWidgets({

    plot_dat<-phe_dat()%>%
                select(phe_name())%>%
                    select_if(is.numeric)
    plist<-lapply(names(plot_dat),function(x){
      ggplotly(ggplot(plot_dat,aes_string(x))+
          geom_histogram(fill="#69b3a2",bins = 30,
                         color="#e9ecef", alpha=0.9)+
          theme_bw()+
          labs(x=x,y="Counts"))
      })

    #p<-subplot(plist,nrows = ceiling(length(plist)/2))
    p<-combineWidgets(list =plist, nrow = ceiling(length(plist)/3))
    p
  })


      
  output$downloadData <- downloadHandler(
        
        filename = function() {
          paste0(tools::file_path_sans_ext(input$phe_upload$name), "_filter.",tools::file_ext(input$phe_upload$name))
        },
        content = function(file) {
          switch(tools::file_ext(input$phe_upload$name),
                 csv = write.csv(phe_dat(), file, row.names = FALSE,quote = F,na=""),
                 txt = write.table(phe_dat(), file, row.names = FALSE,quote = F,na="")
          )
          
        }
      )
    observeEvent(input$start,{
      browser()
    })
    
    
  }
  
  
  
  