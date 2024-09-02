imputeServer<-function(input,output,session){

  output$demo_impute_fam<-DT::renderDataTable(DT::datatable(demo_gs_fam_dat, colnames="",
                                                         rownames=F))
  
  ## gs Data upload------------------------------------------------------------------
    
    impute_density<-reactive({
    req(input$density)
    waiter<-Waiter$new(id='density_summary',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    return(input$density)
    
  })

  ### Genotype
  
  impute_gen_dat<-reactive({
    req(input$impute_gen_upload)
    waiter<-Waiter$new(id='impute_gen_summary',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    ext<-tools::file_ext(input$impute_gen_upload$name)
    path<-unique(dirname(input$impute_gen_upload$datapath))
    
    
    # Rename files with their ordinary names
    full_name<-paste0(dirname(input$impute_gen_upload$datapath),"/",input$impute_gen_upload$name)
    file.rename(input$impute_gen_upload$datapath,full_name)
    
    
    # file names with path
    name<-unique(tools::file_path_sans_ext(full_name))
    
    # file names without path
    file_name_cln<-basename(input$impute_gen_upload$name)
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
    
    dir.create(path_impute,showWarnings = F)
    for (i in 1:length(name)){
      system(paste0("/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink --bfile ",
                    name[i]," --out ",path_impute,"/",name_cln[i]," --make-bed"))
      
      con  <- file(paste0(path_impute,"/",name_cln[i],".log"), open = "r")
      
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
  output$density_summary<-renderUI({
    req(input$density)
    div(paste0("The ",input$density, " has been selected."))
  })

  ### Genotype
  output$impute_gen_summary<-renderUI({
    req(input$impute_gen_upload)
    div(paste0("There are ", sum(impute_gen_dat()$summary$No_var)), "variants in ", 
        length(impute_gen_dat()$ind), " individuals")
  }) 
  
  ## GS submit button

  density_check <- modalDialog(
    div(
      p("Please selecte a target density")
      
    ),
    title = "Select a target density",
    footer = tagList(
      actionButton("cancel", "Cancel", class = "btn btn-warning"),
    )
  )
  observeEvent(input$impute_submit_bttn,{
    if(!isValidEmail(input$impute_email)){
      showModal(email_check)
    }else if(is.null(impute_density)){
      showModal(density_check)
    }else{
      showModal(submit_confirm("impute"))
    }
  })
  
  observeEvent(input$impute_ok,{
    res<-list()
    res$dir<-path_impute
    res$gen_file<-input$impute_gen_upload$name
    res$email<-input$impute_email
    res$task<-paste0(task,"_impute")
    res$density<-input$density
    if(!is.null(input$impute_gen_upload)){
      res$num_of_snps<-sum(impute_gen_dat()$summary$No_var)
    }
    save(res,file=paste0(path_impute,"/par.RData"))
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
    system(paste0("nohup Rscript Impute_job_submit.R /srv/shiny-server/",path_impute," > log/",task, ".log 2>&1 &"))
    task<<-sample(1000:9999,1)
    path_impute <<- paste0("temp/",task,"_impute")
    # #dir.create(path_impute,showWarnings = F)
  })
  
  observeEvent(input$cancel,{
    removeModal()
  })
  
}
