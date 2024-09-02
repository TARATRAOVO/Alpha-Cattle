server_cpi<-function(input,output,session){

  output$traits_symbol_table<-DT::renderDataTable(traits_symbol_table_dat,
                                                   options = list(lengthChange = FALSE))
  output$demo_cpi_phe1<-DT::renderDataTable(demo_trait_type1_dat,
                                                   options = list(lengthChange=FALSE))
  output$demo_cpi_phe2<-DT::renderDataTable(demo_trait_type2_dat,
                                                   options = list(lengthChange=FALSE))
  output$demo_cpi_ped<-DT::renderDataTable(demo_gs_ped_dat,
                                                   options = list(lengthChange=FALSE))
  output$demo_cpi_fam<-DT::renderDataTable(DT::datatable(demo_gs_fam_dat, colnames="",
                                                        rownames=F))
  
  ## GS model check
  cpi_phe_name1<-reactive(names(cpi_phe_dat()))

  observeEvent(input$cpi_phe_upload,{
    updateMultiInput(session = session,inputId="cpi_trt1",choices=cpi_phe_name1())
    updateMultiInput(session = session,inputId="cpi_fct1",choices=cpi_phe_name1())
    updateMultiInput(session = session,inputId="cpi_cov1",choices=cpi_phe_name1())
    updateMultiInput(session = session,inputId="cpi_gen1",choices=cpi_phe_name1())
    updateMultiInput(session = session,inputId="cpi_rnd1",choices=cpi_phe_name1())
  })

  cpi_model_rhs1<-reactive({
    comp1<-c(paste0(input$cpi_fct1, ifelse(is.null(input$cpi_fct1), "", " (fix)")),
            paste0(input$cpi_cov1, ifelse(is.null(input$cpi_cov1), ""," (cov)")),
            paste0(input$cpi_gen1, ifelse(is.null(input$cpi_gen1), ""," (gen)")),
            paste0(input$cpi_rnd1, ifelse(is.null(input$cpi_rnd1), ""," (rnd)"))
    )
    comp1<-comp1[!comp1==""]
    return(paste(comp1,collapse = " + "))
   })
   output$cpi_model_check1 <- renderUI({
     if(is.null(input$cpi_trt1)){
        div("No model specified",class="model_check")
     }else{
       div(HTML(paste(map(input$cpi_trt1,function(x)paste0(x," = ",cpi_model_rhs1())),
                      collapse = "<br/>")),class="model_check")

     }
    })

  cpi_phe_name2<-reactive(names(cpi_phe_dat()))

  observeEvent(input$cpi_phe_upload,{
    updateMultiInput(session = session,inputId="cpi_trt2",choices=cpi_phe_name2())
    updateMultiInput(session = session,inputId="cpi_fct2",choices=cpi_phe_name2())
    updateMultiInput(session = session,inputId="cpi_cov2",choices=cpi_phe_name2())
    updateMultiInput(session = session,inputId="cpi_gen2",choices=cpi_phe_name2())
    updateMultiInput(session = session,inputId="cpi_rnd2",choices=cpi_phe_name2())
  })

  cpi_model_rhs2<-reactive({
    comp2<-c(paste0(input$cpi_fct2, ifelse(is.null(input$cpi_fct2), "", " (fix)")),
            paste0(input$cpi_cov2, ifelse(is.null(input$cpi_cov2), ""," (cov)")),
            paste0(input$cpi_gen2, ifelse(is.null(input$cpi_gen2), ""," (gen)")),
            paste0(input$cpi_rnd2, ifelse(is.null(input$cpi_rnd2), ""," (rnd)"))
    )
    comp2<-comp2[!comp2==""]
    return(paste(comp2,collapse = " + "))
   })
   output$cpi_model_check2 <- renderUI({
     if(is.null(input$cpi_trt2)){
        div("No model specified",class="model_check")
     }else{
       div(HTML(paste(map(input$cpi_trt2,function(x)paste0(x," = ",cpi_model_rhs2())),
                      collapse = "<br/>")),class="model_check")

     }
    })

  ##gs Data upload-----------------------------------------------------------------------------------------------------------------

   cpi_typecpi<-reactive({
    req(input$typecpi)
    waiter<-Waiter$new(id='typecpi_summary',html = spin_fading_circles())
    waiter$show()
    on.exit(waiter$hide())
    return(input$typecpi)

  })

  ###Phenotype
  cpi_phe_dat<-reactive({
     req(input$cpi_phe_upload)
     waiter<-Waiter$new(id='cpi_phe_summary',html = spin_fading_circles())
     waiter$show()
     on.exit(waiter$hide())
     ext<-tools::file_ext(input$cpi_phe_upload$name)
     dat<-switch(ext,
                  csv = data.table::fread(input$cpi_phe_upload$datapath),
                  txt = data.table::fread(input$cpi_phe_upload$datapath),
                  validate("[ERROR:] Invalid phenotype file; Please upload a .csv or .txt file")

    )
    dir.create(path_cpi,showWarnings = F)
    system(paste0("cp ",input$cpi_phe_upload$datapath, " ", path_cpi, "/", input$cpi_phe_upload$name))
    return(dat)
   })
  ###Pedigree
  cpi_ped_dat<-reactive({
     req(input$cpi_ped_upload)
     waiter<-Waiter$new(id='cpi_ped_summary',html = spin_fading_circles())
     waiter$show()
     on.exit(waiter$hide())
     ext<-tools::file_ext(input$cpi_ped_upload$name)
     dat<-switch(ext,
                  csv = data.table::fread(input$cpi_ped_upload$datapath),
                  txt = data.table::fread(input$cpi_ped_upload$datapath),
                  validate("[ERROR:] Invalid pedigree file; Please upload a .csv or .txt file")


      )
      dir.create(path_cpi,showWarnings = F)
      system(paste0("cp ",input$cpi_ped_upload$datapath, " ", path_cpi, "/", input$cpi_ped_upload$name))
      return(dat)
     })
     ###Genotype
     cpi_gen_dat<-reactive({
       req(input$cpi_gen_upload)
       waiter<-Waiter$new(id='cpi_gen_summary',html = spin_fading_circles())
       waiter$show()
       on.exit(waiter$hide())
       ext<-tools::file_ext(input$cpi_gen_upload$name)
       path<-unique(dirname(input$cpi_gen_upload$datapath))


       # Rename files with their ordinary names
       full_name<-paste0(dirname(input$cpi_gen_upload$datapath),"/",input$cpi_gen_upload$name)
       file.rename(input$cpi_gen_upload$datapath,full_name)


        # file names with path
        name<-unique(tools::file_path_sans_ext(full_name))

        # file names without path
        file_name_cln<-basename(input$cpi_gen_upload$name)
        name_cln<-unique(tools::file_path_sans_ext(file_name_cln))
 
         # Supposed to be these files
         supp_files<-unlist(map(name,~paste0(.x,c(".bed",".bim",".fam"))))
         supp_files_cln<-unlist(map(name_cln,~paste0(.x,c(".bed",".bim",".fam"))))

         #Check file extension
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

           dir.create(path_cpi,showWarnings =F)
           for (i in 1:length(name)){
             system(paste0("/public/home/chenjc/miniconda3/pkgs/plink-1.90b6.21-h516909a_0/bin/plink --bfile ",
                           name[i]," --out ",path_cpi,"/", name_cln[i], " --make-bed"))

             con <- file(paste0(path_cpi,"/",name_cln[i],".log"), open = "r")

             while (length(oneLine <- readLines(con, n=1, warn=FALSE)) >0) {

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
          ## cpi data check
          ## type
          output$typecpi_summary<-renderUI({
            req(input$typecpi)
            div(paste0("The ",input$typecpi, " has been selected."))
        })

          ### Phenotype
          output$cpi_phe_summary<-renderUI({
            req(input$cpi_phe_upload)
            div(paste0("There are ", nrow(unique(cpi_phe_dat())),
                       " records for ",length(unique(cpi_phe_dat()[[1]])),
                       " individuals in the phenotype file."))
        })

          ### Pedigree
          #### Check number of individuals in pedigree
          output$cpi_ped_summary<-renderUI({
             req(input$cpi_ped_upload)
             div(paste0("There are ", nrow(unique(cpi_ped_dat()))), "individuals in the pedigree.")
        })

            #### Check number of individuals in pedigree having phenotype
            output$cpi_phe_ped_summary<-renderUI({
               req(input$cpi_ped_upload, input$cpi_phe_upload)
               div(paste0("There are ", sum(cpi_ped_dat()[[1]]%in%cpi_phe_dat()[[1]]),
                           " individuals in pedigree having phenotype."))
        })
         
          ### Genotype
          output$cpi_gen_summary<-renderUI({
             req(input$cpi_gen_upload)
             div(paste0("There are ", sum(cpi_gen_dat()$summary$No_var)), "variants in ",
                 length(cpi_gen_dat()$ind), " individuals")
        })

           output$cpi_phe_gen_summary<-renderUI({
              req(input$cpi_gen_upload, input$cpi_phe_upload)
              div(paste0("There are ", sum(cpi_gen_dat()$ind%in%cpi_phe_dat()[[1]]),
                         " genotyped individuals having phenotype."))
        })

            output$cpi_ped_gen_summary<-renderUI({
               req(input$cpi_gen_upload, input$cpi_ped_upload)
               div(paste0("There are ", sum(cpi_gen_dat()$ind%in%cpi_ped_dat()[[1]]),
                        " genotyped individuals in pedigree"))
       })
           
        ## CPI submit button
        cpi_method_selected<-reactive({
          if(!is.null(input$cpi_phe_upload) & is.null(input$cpi_ped_upload) & !is.null(input$cpi_gen_upload)){
            method<-"GBLUP"
          }else if(!is.null(input$cpi_phe_upload) & !is.null(input$cpi_ped_upload) & is.null(input$cpi_gen_upload)){
            method<-"BLUP"
          }else{
             method<-NULL
          }
          return(method)
        })

typecpi_check <- modalDialog(
    div(
      p("Please select a Method used in Milk Production Traits EBV calculation process")

    ),
    title = "Select a method",
    footer = tagList(
      actionButton("cancel", "Cancel", class = "btn btn-warning"),
    )
  )

        observeEvent(input$cpi_submit_bttn,{
          if(!isValidEmail(input$cpi_email)){
            showModal(email_check)
          }else if(is.null(cpi_method_selected())){
            showModal(cpi_method_check)
          }else if(is.null(cpi_typecpi)){
             showModal(typecpi_check)
          }else{
            showModal(submit_confirm("cpi"))
          }
        })

         observeEvent(input$cpi_ok,{
           res<-list()
           res$method<-cpi_method_selected()
           res$model1<-list()
           res$model1$traits<-input$cpi_trt1
           res$model1$fix_fct<-input$cpi_fct1
           res$model1$fix_cov<-input$cpi_cov1
           res$model1$rnd_gen<-input$cpi_gen1
           res$model1$rnd_non<-input$cpi_rnd1
           res$model2<-list()
           res$model2$traits<-input$cpi_trt2
           res$model2$fix_fct<-input$cpi_fct2
           res$model2$fix_cov<-input$cpi_cov2
           res$model2$rnd_gen<-input$cpi_gen2
           res$model2$rnd_non<-input$cpi_rnd2
           res$dir<-path_cpi
           res$phe_file<-input$cpi_phe_upload$name
           res$ped_file<-input$cpi_ped_upload$name
           res$gen_file<-input$cpi_gen_upload$name
           res$email<-input$cpi_email
           res$task<-paste0(task,"_cpi")
           res$typecpi<-input$typecpi
           if(!is.null(input$cpi_gen_upload)){
            res$num_of_snps<-sum(cpi_gen_dat()$summary$No_var)
           }
           save(res,file=paste0(path_cpi,"/par.RData"))
           showNotification("Submitted. We will send email to you, once the job is done.
                              The page will be refreshed in 5 seconds", type="message")
           removeModal()

           Sys.sleep(1)
           showNotification("4 seconds", type="error")
           Sys.sleep(1)
           showNotification("3 seconds", type="error")
           Sys.sleep(1)
           showNotification("2 seconds", type="error")
           Sys.sleep(1)
           showNotification("1 seconds", type="error")
           Sys.sleep(1)
           refresh()
           ##system(paste0("nohup Rscript CPI_job_submit.R /srv/shiny-server/" ,path_cpi," > log/",task,".log 2>&1 &"))
           task<<-sample(1000:9999,1)
           path_cpi<<-paste0("temp/",task,"_cpi")
           #dir.create(path_cpi,showWarnings = F)
        })

        observeEvent(input$cancel,{
          removeModal()
        })
}

             