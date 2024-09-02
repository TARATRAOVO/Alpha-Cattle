server_gs<-function(input,output,session){
  ## GS demo phe table
  output$demo_gs_phe<-DT::renderDataTable(demo_gs_phe_dat, 
                                          options = list(lengthChange = FALSE))
  
  output$demo_gs_ped<-DT::renderDataTable(demo_gs_ped_dat, 
                                          options = list(lengthChange = FALSE))
  
  output$demo_gs_fam<-DT::renderDataTable(DT::datatable(demo_gs_fam_dat, colnames="",
                                                        rownames=F))
  
  output$demo_gs_var<-DT::renderDataTable(DT::datatable(demo_gs_var_dat, colnames="",
                                                        rownames=F))
  ## GS model check
  gs_phe_name<-reactive(names(gs_phe_dat()))
  
  observeEvent(input$gs_phe_upload,{
    updateMultiInput(session = session,inputId="gs_trt",choices=gs_phe_name())
    updateMultiInput(session = session,inputId="gs_fct",choices=gs_phe_name())
    updateMultiInput(session = session,inputId="gs_cov",choices=gs_phe_name())
    updateMultiInput(session = session,inputId="gs_gen",choices=gs_phe_name())
    updateMultiInput(session = session,inputId="gs_rnd",choices=gs_phe_name())
  })
  
  gs_model_rhs<-reactive({
    comp<-c(paste0(input$gs_fct, ifelse(is.null(input$gs_fct), "", " (fix)")),
            paste0(input$gs_cov, ifelse(is.null(input$gs_cov), ""," (cov)")),
            paste0(input$gs_gen, ifelse(is.null(input$gs_gen), ""," (gen)")),
            paste0(input$gs_rnd, ifelse(is.null(input$gs_rnd), ""," (rnd)"))
    )
    comp<-comp[!comp==""]
    return(paste(comp,collapse = " + "))
  })
  output$gs_model_check <- renderUI({
    if(is.null(input$gs_trt)){
      div("No model specified",class="model_check")
    }else{
      div(HTML(paste(map(input$gs_trt,function(x)paste0(x," = ",gs_model_rhs())),
                     collapse = "<br/>")),class="model_check")
    }
  })
  ## gs Data upload------------------------------------------------------------------
 
   gs_typegs<-reactive({
    req(input$typegs)
    waiter<-Waiter$new(id='typegs_summary',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    return(input$typegs)
    
  })

  ### Phenotype
  gs_phe_dat<-reactive({
    req(input$gs_phe_upload)
    waiter<-Waiter$new(id='gs_phe_summary',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    ext<-tools::file_ext(input$gs_phe_upload$name)
    dat<-switch(ext,
                csv = data.table::fread(input$gs_phe_upload$datapath),
                txt = data.table::fread(input$gs_phe_upload$datapath),
                validate("[ERROR:] Invalid phenotype file; Please upload a .csv or .txt file")
                
    )
    dir.create(path_gs,showWarnings = F)
    system(paste0("cp ",input$gs_phe_upload$datapath, " ", path_gs, "/", input$gs_phe_upload$name))
    return(dat)
  })
  ### Pedigree
  gs_ped_dat<-reactive({
    req(input$gs_ped_upload)
    waiter<-Waiter$new(id='gs_ped_summary',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    ext<-tools::file_ext(input$gs_ped_upload$name)
    dat<-switch(ext,
                csv = data.table::fread(input$gs_ped_upload$datapath),
                txt = data.table::fread(input$gs_ped_upload$datapath),
                validate("[ERROR:] Invalid pedigree file; Please upload a .csv or .txt file")
                
    )
    
    dir.create(path_gs,showWarnings = F)
    system(paste0("cp ",input$gs_ped_upload$datapath, " ", path_gs, "/", input$gs_ped_upload$name))
    return(dat)
  })
  
  ### Genotype
  
  gs_gen_dat<-reactive({
    req(input$gs_gen_upload)
    waiter<-Waiter$new(id='gs_gen_summary',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    ext<-tools::file_ext(input$gs_gen_upload$name)
    path<-unique(dirname(input$gs_gen_upload$datapath))
    
    
    # Rename files with their ordinary names
    full_name<-paste0(dirname(input$gs_gen_upload$datapath),"/",input$gs_gen_upload$name)
    file.rename(input$gs_gen_upload$datapath,full_name)
    
    
    # file names with path
    name<-unique(tools::file_path_sans_ext(full_name))
    
    # file names without path
    file_name_cln<-basename(input$gs_gen_upload$name)
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
    
    dir.create(path_gs,showWarnings = F)
    for (i in 1:length(name)){
      system(paste0("/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink --bfile ",
                    name[i]," --out ",path_gs,"/", name_cln[i], " --make-bed"))
      
      con  <- file(paste0(path_gs,"/",name_cln[i],".log"), open = "r")
      
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
  
  ###Variance component file
  gs_var_dat<-reactive({
    req(input$gs_var_upload)
    ext<-tools::file_ext(input$gs_var_upload$name)
    dat<-switch(ext,
                csv = data.table::fread(input$gs_var_upload$datapath),
                txt = data.table::fread(input$gs_var_upload$datapath),
                validate("[ERROR:] Invalid variance component file; Please upload a .csv or .txt file")
                
    )
    
    dir.create(path_gs,showWarnings = F)
    system(paste0("cp ",input$gs_var_upload$datapath, " ", path_gs, "/", input$gs_var_upload$name))
    return(dat)
  })
  ## gs data check
  ## type
  output$typegs_summary<-renderUI({
     req(input$typegs)
     div(paste0("The ",input$typegs, " has been selected."))
  })
  ### Phenotype
  output$gs_phe_summary<-renderUI({
    req(input$gs_phe_upload)
    div(paste0("There are ", nrow(unique(gs_phe_dat())), 
               " records for ",length(unique(gs_phe_dat()[[1]])),
               " individuals in the phenotype file."))
  })
  
  ### Pedigree
  #### Check number of individuals in pedigree
  output$gs_ped_summary<-renderUI({
    req(input$gs_ped_upload)
    div(paste0("There are ", nrow(unique(gs_ped_dat()))), "individuals in the pedigree.")
  })
  
  #### Check number of individuals in pedigree having phenotype
  output$gs_phe_ped_summary<-renderUI({
    req(input$gs_ped_upload, input$gs_phe_upload)
    div(paste0("There are ", sum(as.character(gs_ped_dat()[[1]])%in%as.character(gs_phe_dat()[[1]])),
               " individuals in pedigree having phenotype."))
  })
  
  ### Genotype
  output$gs_gen_summary<-renderUI({
    req(input$gs_gen_upload)
    div(paste0("There are ", sum(gs_gen_dat()$summary$No_var)), "variants in ", 
        length(gs_gen_dat()$ind), " individuals")
  })
  
  output$gs_phe_gen_summary<-renderUI({
    req(input$gs_gen_upload, input$gs_phe_upload)
    div(paste0("There are ", sum(as.character(gs_gen_dat()$ind)%in%as.character(gs_phe_dat()[[1]])),
               " genotyped individuals having phenotype."))
  })
  
  output$gs_ped_gen_summary<-renderUI({
    req(input$gs_gen_upload, input$gs_ped_upload)
    div(paste0("There are ", sum(as.character(gs_gen_dat()$ind)%in%as.character(gs_ped_dat()[[1]])),
               " genotyped individuals in pedigree"))
  })
  
  ### Variance component
  output$gs_var_summary<-renderUI({
    req(input$gs_var_upload)
    div(paste0(ncol(gs_var_dat()), " variance components for ", nrow(gs_var_dat())," traits are loaded."))
  })
  
  ## GS submit button
  gs_method_selected<-reactive({
    if(!is.null(input$gs_phe_upload) & !is.null(input$gs_ped_upload) & is.null(input$gs_gen_upload) & is.null(input$gs_var_upload)){
      method<-"BLUP"
    }else if(!is.null(input$gs_phe_upload) & !is.null(input$gs_ped_upload) & !is.null(input$gs_gen_upload) & is.null(input$gs_var_upload)){
      method<-"ssGBLUP"
    }else if(!is.null(input$gs_phe_upload) & is.null(input$gs_ped_upload) & !is.null(input$gs_gen_upload) & is.null(input$gs_var_upload)){
      method<-"GBLUP"
    }else if(!is.null(input$gs_phe_upload) & !is.null(input$gs_ped_upload) & is.null(input$gs_gen_upload) & !is.null(input$gs_var_upload)){
      method<-"BLUP_VAR"
    }else if(!is.null(input$gs_phe_upload) & !is.null(input$gs_ped_upload) & !is.null(input$gs_gen_upload) & !is.null(input$gs_var_upload)){
      method<-"ssGBLUP_VAR"
    }else if(!is.null(input$gs_phe_upload) & is.null(input$gs_ped_upload) & !is.null(input$gs_gen_upload) & !is.null(input$gs_var_upload)){
      method<-"GBLUP_VAR"
    }else{
      method<-NULL
    }
    return(method)
  })

typegs_check <- modalDialog(
    div(
      p("Please selecte a Phenotype type")
      
    ),
    title = "Select a type",
    footer = tagList(
      actionButton("cancel", "Cancel", class = "btn btn-warning"),
    )
  )
  
  observeEvent(input$gs_submit_bttn,{
    if(!isValidEmail(input$gs_email)){
      showModal(email_check)
    }else if(is.null(gs_method_selected())){
      showModal(gs_method_check)
    }else if(is.null(gs_typegs)){
      showModal(typegs_check)
    }else{
      showModal(submit_confirm("gs"))
    }
  })
  
  observeEvent(input$gs_ok,{
    res<-list()
    res$method<-gs_method_selected()
    res$model<-list()
    res$model$traits<-input$gs_trt
    res$model$fix_fct<-input$gs_fct
    res$model$fix_cov<-input$gs_cov
    res$model$rnd_gen<-input$gs_gen
    res$model$rnd_non<-input$gs_rnd
    res$dir<-path_gs
    res$phe_file<-input$gs_phe_upload$name
    res$ped_file<-input$gs_ped_upload$name
    res$gen_file<-input$gs_gen_upload$name
    res$var_file<-input$gs_var_upload$name
    res$email<-input$gs_email
    res$task<-paste0(task,"_gs")
    res$typegs<-input$typegs
    if(!is.null(input$gs_gen_upload)){
      res$num_of_snps<-sum(gs_gen_dat()$summary$No_var)
    }
    save(res,file=paste0(path_gs,"/par.RData"))
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
    system(paste0("nohup Rscript GS_job_submit.R /srv/shiny-server/",path_gs," > log/",task,".log 2>&1 &"))
    task<<-sample(1000:9999,1)
    path_gs <<- paste0("temp/",task,"_gs")
    #dir.create(path_gs,showWarnings = F)
  })
  
  observeEvent(input$cancel,{
    removeModal()
  })
}
