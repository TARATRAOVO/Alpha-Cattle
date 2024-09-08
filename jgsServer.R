server_jgs<-function(input,output,session){
  ## joint
  output$demo_jgs_phe<-DT::renderDataTable(demo_gs_phe_dat, 
                                           options = list(lengthChange = FALSE))
  
  output$demo_jgs_ped<-DT::renderDataTable(demo_gs_ped_dat, 
                                           options = list(lengthChange = FALSE))
  
  output$demo_jgs_fam<-DT::renderDataTable(DT::datatable(demo_gs_fam_dat, colnames="",
                                                         rownames=F))
  
  
  ##joint GS model check
  jgs_phe_name<-reactive(names(jgs_phe_dat()))
  
  observeEvent(input$jgs_phe_upload,{
    updateMultiInput(session = session,inputId="jgs_trt",choices=jgs_phe_name())
    updateMultiInput(session = session,inputId="jgs_fct",choices=jgs_phe_name())
    updateMultiInput(session = session,inputId="jgs_cov",choices=jgs_phe_name())
    updateMultiInput(session = session,inputId="jgs_gen",choices=jgs_phe_name())
    updateMultiInput(session = session,inputId="jgs_rnd",choices=jgs_phe_name())
  })
  
  jgs_model_rhs<-reactive({
    comp<-c(paste0(input$jgs_fct, ifelse(is.null(input$jgs_fct), "", " (fix)")),
            paste0(input$jgs_cov, ifelse(is.null(input$jgs_cov), ""," (cov)")),
            paste0(input$jgs_gen, ifelse(is.null(input$jgs_gen), ""," (gen)")),
            paste0(input$jgs_rnd, ifelse(is.null(input$jgs_rnd), ""," (rnd)"))
    )
    comp<-comp[!comp==""]
    return(paste(comp,collapse = " + "))
  })
  
  output$jgs_model_check <- renderUI({
    if(is.null(input$jgs_trt)){
      div("No model specified",class="model_check")
    }else{
      div(HTML(paste(map(input$jgs_trt,function(x)paste0(x," = ",jgs_model_rhs())),
                     collapse = "<br/>")),class="model_check")
    }
  })
  
  ## gs Data upload------------------------------------------------------------------
  jgs_breed<-reactive({
    req(input$breed)
    waiter<-Waiter$new(id='breed_summary',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    return(input$breed)
    
  })
 
  jgs_type<-reactive({
    req(input$type)
    waiter<-Waiter$new(id='type_summary',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    return(input$type)
    
  })
  
  jgs_phe_dat<-reactive({
    req(input$jgs_phe_upload)
    waiter<-Waiter$new(id='jgs_phe_summary',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    ext<-tools::file_ext(input$jgs_phe_upload$name)
    dat<-switch(ext,
                csv = data.table::fread(input$jgs_phe_upload$datapath),
                txt = data.table::fread(input$jgs_phe_upload$datapath),
                validate("[ERROR:] Invalid phenotype file; Please upload a .csv or .txt file")
                
    )
    dir.create(path_jgs,showWarnings = F)
    system(paste0("cp ",input$jgs_phe_upload$datapath, " ", path_jgs, "/", input$jgs_phe_upload$name))
    return(dat)
  })
  ### Pedigree
  jgs_ped_dat<-reactive({
    req(input$jgs_ped_upload)
    waiter<-Waiter$new(id='jgs_ped_summary',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    ext<-tools::file_ext(input$jgs_ped_upload$name)
    dat<-switch(ext,
                csv = data.table::fread(input$jgs_ped_upload$datapath),
                txt = data.table::fread(input$jgs_ped_upload$datapath),
                validate("[ERROR:] Invalid pedigree file; Please upload a .csv or .txt file")
                
    )
    
    dir.create(path_jgs,showWarnings = F)
    system(paste0("cp ",input$jgs_ped_upload$datapath, " ", path_jgs, "/", input$jgs_ped_upload$name))
    return(dat)
  })
  ### Genotype
  
  jgs_gen_dat<-reactive({
    req(input$jgs_gen_upload)
    waiter<-Waiter$new(id='jgs_gen_summary',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    ext<-tools::file_ext(input$jgs_gen_upload$name)
    path<-unique(dirname(input$jgs_gen_upload$datapath))
    
    
    # Rename files with their ordinary names
    full_name<-paste0(dirname(input$jgs_gen_upload$datapath),"/",input$jgs_gen_upload$name)
    file.rename(input$jgs_gen_upload$datapath,full_name)
    
    
    # file names with path
    name<-unique(tools::file_path_sans_ext(full_name))
    
    # file names without path
    file_name_cln<-basename(input$jgs_gen_upload$name)
    name_cln<-unique(tools::file_path_sans_ext(file_name_cln))
    
    # Supposed to be these files
    supp_files<-unlist(map(name,~paste0(.x,c(".bed",".bim",".fam"))))
    supp_files_cln<-unlist(map(name_cln,~paste0(.x,c(".bed",".bim",".fam"))))
    
    # Check file extension
    file_err_ext<-!ext%in%c("bim","bed","fam")
    if(any(file_err_ext)){
      validate(paste0("[ERROR:] Invalid file type (",
                      file_name_cln[file_err_ext],
                      ") ; Please upload .bed, .bim and .fam files"))
    }
    
    # Check if some files are supplied in supposed file list
    file_err<-!full_name%in%supp_files
    if(any(file_err)){
      validate(paste0("[ERROR:] Invalid files (",
                      file_name_cln[file_err],
                      ") ; Please check file names"))
    }
    
    # Check if some files that are out of supposed file list
    file_not_exist<-!supp_files%in%full_name
    if(any(file_not_exist)){
      validate(paste0("[ERROR:] We cannot find some files (",
                      supp_files_cln[file_not_exist],
                      ") or you provide wrong file names, please check; "))
    }
    
    # Check individuals
    fam<-data.table::fread(paste0(path,"/", name_cln[1], ".fam"))
    
    if(length(name)>1){
      for (i in 2:length(name)){
        fam2<-data.table::fread(paste0(path,"/", name_cln[i], ".fam"))
        if(length(fam2[[1]])!=length(fam[[1]])){
          validate("[ERROR:] Number of individuals must be the same for different plink fam files")
        }else if(!all(fam2[[1]]==fam[[1]]) | !all(fam2[[2]]==fam[[2]])){
          validate("[ERROR:] Individuals in different fam files must be the same ")
        }
      }
    }
    
    # Plink check format
    num_of_variant<-rep(NA,length(name))
    num_of_individual<-rep(NA,length(name))
    
    dir.create(path_jgs,showWarnings = F)
    for (i in 1:length(name)){
      system(paste0("./plink/plink --bfile ",
                    name[i]," --out ",path_jgs,"/",name_cln[i]," --make-bed"))
      
      con  <- file(paste0(path_jgs,"/",name_cln[i],".log"), open = "r")
      
      while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
        
        if(str_detect(oneLine, "^\\d+(?= variant)")){
          num_of_variant[i] <- as.numeric(str_extract(oneLine, "\\d+(?= variant)"))
        }
        
        if(str_detect(oneLine, "^\\d+(?= variant)")){
          num_of_individual[i] <- as.numeric(str_extract(oneLine, "\\d+(?= people)"))
        }
        
        if(str_detect(oneLine,'Error')){
          validate(
            need(!str_detect(oneLine,'Error'), 
                 paste0("[ERROR:] Plink detect an error for ", 
                        name_cln[i], " files. Please check!")
            )
          )
          break
        }
      } 
      
      close(con)
    }
    
    res<-list()
    res$summary<-data.frame(Files=name_cln,No_ind=num_of_individual,No_var=num_of_variant)
    res$ind<-fam[[1]]
    return(res)
  })
  ## data check
  output$breed_summary<-renderUI({
    req(input$breed)
    div(paste0("The ",input$breed, " has been selected."))
  })
  ## type
  output$type_summary<-renderUI({
     req(input$type)
     div(paste0("The ",input$type, " has been selected."))
  })
  ### Phenotype
  output$jgs_phe_summary<-renderUI({
    req(input$jgs_phe_upload)
    div(paste0("There are ", nrow(unique(jgs_phe_dat())), 
               " records for ",length(unique(jgs_phe_dat()[[1]])),
               " individuals in the phenotype file."))
  })
  
  ### Pedigree
  #### Check number of individuals in pedigree
  output$jgs_ped_summary<-renderUI({
    req(input$jgs_ped_upload)
    div(paste0("There are ", nrow(unique(jgs_ped_dat()))), "individuals in the pedigree.")
  })
  
  #### Check number of individuals in pedigree having phenotype
  output$jgs_phe_ped_summary<-renderUI({
    req(input$jgs_ped_upload, input$jgs_phe_upload)
    div(paste0("There are ", sum(as.character(jgs_ped_dat()[[1]])%in%as.character(jgs_phe_dat()[[1]])),
               " individuals in pedigree having phenotype."))
  })
  
  ### Genotype
  output$jgs_gen_summary<-renderUI({
    req(input$jgs_gen_upload)
    div(paste0("There are ", sum(jgs_gen_dat()$summary$No_var)), "variants in ", 
        length(jgs_gen_dat()$ind), " individuals")
  })
  
  output$jgs_phe_gen_summary<-renderUI({
    req(input$jgs_gen_upload, input$jgs_phe_upload)
    div(paste0("There are ", sum(as.character(jgs_gen_dat()$ind)%in%as.character(jgs_phe_dat()[[1]])),
               " genotyped individuals having phenotype."))
  })
  
  output$jgs_ped_gen_summary<-renderUI({
    req(input$jgs_gen_upload, input$jgs_ped_upload)
    div(paste0("There are ", sum(as.character(jgs_gen_dat()$ind)%in%as.character(jgs_ped_dat()[[1]])),
               " genotyped individuals in pedigree"))
  })
  
  
  
  ## GS submit button
  jgs_method_selected<-reactive({
    
    if(!is.null(input$jgs_phe_upload) & !is.null(input$jgs_ped_upload) & !is.null(input$jgs_gen_upload) ){
      method<-"ssGBLUP"
    }else if(!is.null(input$jgs_phe_upload) & is.null(input$jgs_ped_upload) & !is.null(input$jgs_gen_upload) ){
      method<-"GBLUP"
    }else{
      method<-NULL
    }
    return(method)
  })
  breed_check <- modalDialog(
    div(
      p("Please selecte a breed")
      
    ),
    title = "Select a breed",
    footer = tagList(
      actionButton("cancel", "Cancel", class = "btn btn-warning"),
    )
  )

type_check <- modalDialog(
    div(
      p("Please selecte a Phenotype type")
      
    ),
    title = "Select a type",
    footer = tagList(
      actionButton("cancel", "Cancel", class = "btn btn-warning"),
    )
  )

  observeEvent(input$jgs_submit_bttn,{
    if(!isValidEmail(input$jgs_email)){
      showModal(email_check)
    }else if(is.null(jgs_method_selected())){
      showModal(jgs_method_check)
    }else if(is.null(jgs_breed)){
      showModal(breed_check)
    }else if(is.null(jgs_type)){
      showModal(type_check)
    }else{
      showModal(submit_confirm("jgs"))
    }
  })
  
  observeEvent(input$jgs_ok,{
    res<-list()
    res$method<-jgs_method_selected()
    res$model<-list()
    res$model$traits<-input$jgs_trt
    res$model$fix_fct<-input$jgs_fct
    res$model$fix_cov<-input$jgs_cov
    res$model$rnd_gen<-input$jgs_gen
    res$model$rnd_non<-input$jgs_rnd
    res$dir<-path_jgs
    res$phe_file<-input$jgs_phe_upload$name
    res$ped_file<-input$jgs_ped_upload$name
    res$gen_file<-input$jgs_gen_upload$name
    res$var_file<-input$jgs_var_upload$name
    res$email<-input$jgs_email
    res$task<-paste0(task,"_jgs")
    res$breed<-input$breed
    res$type<-input$type
    if(!is.null(input$jgs_gen_upload)){
      res$num_of_snps<-sum(jgs_gen_dat()$summary$No_var)
    }
    save(res,file=paste0(path_jgs,"/par.RData"))
    showNotification("Submitted. We will send e-mail to you, once the job is done. 
                     The page will be refreshed in 5 seconds", type="message")
    removeModal()
    
    Sys.sleep(1)
    showNotification("4 seconds", type="error")
    Sys.sleep(1)
    showNotification("3 seconds", type="error")
    Sys.sleep(1)
    showNotification("2 seconds", type="error")
    Sys.sleep(1)
    showNotification("1 second", type="error")
    Sys.sleep(1)
    refresh()
    ##system(paste0("nohup Rscript joint_GS_job_submit.R /srv/shiny-server/",path_jgs," > log/",task, ".log 2>&1 &"))
    task<<-sample(1000:9999,1)
    path_jgs <<- paste0("temp/",task,"_jgs")
    # #dir.create(path_jgs,showWarnings = F)
  })
  
  observeEvent(input$cancel,{
    removeModal()
  })
  
}
