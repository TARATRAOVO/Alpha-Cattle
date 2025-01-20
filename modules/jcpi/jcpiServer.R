server_jcpi<-function(input,output,session){

  output$jtraits_symbol_table<-DT::renderDataTable(traits_symbol_table_dat,
                                                   options = list(lengthChange = FALSE))
  output$demo_jcpi_phe<-DT::renderDataTable(demo_trait_type1_dat,
                                                   options = list(lengthChange=FALSE))
  output$demo_jcpi_fam<-DT::renderDataTable(DT::datatable(demo_gs_fam_dat, colnames="",
                                                        rownames=F))
  
  ## GS model check
  jcpi_phe_name1<-reactive(names(jcpi_phe_dat()))

  observeEvent(input$jcpi_phe_upload,{
    updateMultiInput(session = session,inputId="jcpi_trt1",choices=jcpi_phe_name1())
    updateMultiInput(session = session,inputId="jcpi_fct1",choices=jcpi_phe_name1())
    updateMultiInput(session = session,inputId="jcpi_cov1",choices=jcpi_phe_name1())
    updateMultiInput(session = session,inputId="jcpi_gen1",choices=jcpi_phe_name1())
    updateMultiInput(session = session,inputId="jcpi_rnd1",choices=jcpi_phe_name1())
  })

  jcpi_model_rhs1<-reactive({
    comp1<-c(paste0(input$jcpi_fct1, ifelse(is.null(input$jcpi_fct1), "", " (fix)")),
            paste0(input$jcpi_cov1, ifelse(is.null(input$jcpi_cov1), ""," (cov)")),
            paste0(input$jcpi_gen1, ifelse(is.null(input$jcpi_gen1), ""," (gen)")),
            paste0(input$jcpi_rnd1, ifelse(is.null(input$jcpi_rnd1), ""," (rnd)"))
    )
    comp1<-comp1[!comp1==""]
    return(paste(comp1,collapse = " + "))
   })
   output$jcpi_model_check1 <- renderUI({
     if(is.null(input$jcpi_trt1)){
        div("No model specified",class="model_check")
     }else{
       div(HTML(paste(map(input$jcpi_trt1,function(x)paste0(x," = ",jcpi_model_rhs1())),
                      collapse = "<br/>")),class="model_check")

     }
    })

  jcpi_phe_name2<-reactive(names(jcpi_phe_dat()))

  observeEvent(input$jcpi_phe_upload,{
    updateMultiInput(session = session,inputId="jcpi_trt2",choices=jcpi_phe_name2())
    updateMultiInput(session = session,inputId="jcpi_fct2",choices=jcpi_phe_name2())
    updateMultiInput(session = session,inputId="jcpi_cov2",choices=jcpi_phe_name2())
    updateMultiInput(session = session,inputId="jcpi_gen2",choices=jcpi_phe_name2())
    updateMultiInput(session = session,inputId="jcpi_rnd2",choices=jcpi_phe_name2())
  })

  jcpi_model_rhs2<-reactive({
    comp2<-c(paste0(input$jcpi_fct2, ifelse(is.null(input$jcpi_fct2), "", " (fix)")),
            paste0(input$jcpi_cov2, ifelse(is.null(input$jcpi_cov2), ""," (cov)")),
            paste0(input$jcpi_gen2, ifelse(is.null(input$jcpi_gen2), ""," (gen)")),
            paste0(input$jcpi_rnd2, ifelse(is.null(input$jcpi_rnd2), ""," (rnd)"))
    )
    comp2<-comp2[!comp2==""]
    return(paste(comp2,collapse = " + "))
   })
   output$jcpi_model_check2 <- renderUI({
     if(is.null(input$jcpi_trt2)){
        div("No model specified",class="model_check")
     }else{
       div(HTML(paste(map(input$jcpi_trt2,function(x)paste0(x," = ",jcpi_model_rhs2())),
                      collapse = "<br/>")),class="model_check")

     }
    })

  ##gs Data upload-----------------------------------------------------------------------------------------------------------------
  ###Phenotype
  jcpi_phe_dat<-reactive({
     req(input$jcpi_phe_upload)
     waiter<-Waiter$new(id='jcpi_phe_summary',html = spin_fading_circles())
     waiter$show()
     on.exit(waiter$hide())
     ext<-tools::file_ext(input$jcpi_phe_upload$name)
     dat<-switch(ext,
                  csv = data.table::fread(input$jcpi_phe_upload$datapath),
                  txt = data.table::fread(input$jcpi_phe_upload$datapath),
                  validate("[ERROR:] Invalid phenotype file; Please upload a .csv or .txt file")


      )
      dir.create(path_jcpi,showWarnings = F)
      system(paste0("cp ",input$jcpi_phe_upload$datapath, " ", path_jcpi, "/", input$jcpi_phe_upload$name))
      return(dat)
     })
     ###Genotype
     jcpi_gen_dat<-reactive({
       req(input$jcpi_gen_upload)
       waiter<-Waiter$new(id='jcpi_gen_summary',html = spin_fading_circles())
       waiter$show()
       on.exit(waiter$hide())
       ext<-tools::file_ext(input$jcpi_gen_upload$name)
       path<-unique(dirname(input$jcpi_gen_upload$datapath))


       # Rename files with their ordinary names
       full_name<-paste0(dirname(input$jcpi_gen_upload$datapath),"/",input$jcpi_gen_upload$name)
       file.rename(input$jcpi_gen_upload$datapath,full_name)


        # file names with path
        name<-unique(tools::file_path_sans_ext(full_name))

        # file names without path
        file_name_cln<-basename(input$jcpi_gen_upload$name)
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

           dir.create(path_jcpi,showWarnings =F)
           for (i in 1:length(name)){
             system(paste0("./plink/plink --bfile ",
                           name[i]," --out ",path_jcpi,"/", name_cln[i], " --make-bed"))

             con <- file(paste0(path_jcpi,"/",name_cln[i],".log"), open = "r")

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
          ## jcpi data check
          
          ### Phenotype
          output$jcpi_phe_summary<-renderUI({
            req(input$jcpi_phe_upload)
            div(paste0("There are ", nrow(unique(jcpi_phe_dat())),
                       " records for ",length(unique(jcpi_phe_dat()[[1]])),
                       " individuals in the phenotype file."))
        })
         
          ### Genotype
          output$jcpi_gen_summary<-renderUI({
             req(input$jcpi_gen_upload)
             div(paste0("There are ", sum(jcpi_gen_dat()$summary$No_var)), "variants in ",
                 length(jcpi_gen_dat()$ind), " individuals")
        })

           output$jcpi_phe_gen_summary<-renderUI({
              req(input$jcpi_gen_upload, input$jcpi_phe_upload)
              div(paste0("There are ", sum(jcpi_gen_dat()$ind%in%jcpi_phe_dat()[[1]]),
                         " genotyped individuals having phenotype."))
        })

           output$jcpi_phe_gen_summary<-renderUI({
              req(input$jcpi_gen_upload, input$jcpi_phe_upload)
              div(paste0("There are ", sum(jcpi_gen_dat()$ind%in%jcpi_phe_dat()[[1]]),
                         " genotyped individuals having phenotype."))
        })
           
        ## GCPI submit button
        jcpi_method_selected<-reactive({
          if(!is.null(input$jcpi_phe_upload) & !is.null(input$jcpi_gen_upload)){
            method<-"GBLUP"
          }else{
             method<-NULL
          }
          return(method)
        })

        observeEvent(input$jcpi_submit_bttn,{
          if(!isValidEmail(input$jcpi_email)){
            showModal(email_check)
          }else if(is.null(jcpi_method_selected())){
            showModal(jcpi_method_check)
          }else{
            showModal(submit_confirm("jcpi"))
          }
        })

         observeEvent(input$jcpi_ok,{
           res<-list()
           res$method<-jcpi_method_selected()
           res$model1<-list()
           res$model1$traits<-input$jcpi_trt1
           res$model1$fix_fct<-input$jcpi_fct1
           res$model1$fix_cov<-input$jcpi_cov1
           res$model1$rnd_gen<-input$jcpi_gen1
           res$model1$rnd_non<-input$jcpi_rnd1
           res$model2<-list()
           res$model2$traits<-input$jcpi_trt2
           res$model2$fix_fct<-input$jcpi_fct2
           res$model2$fix_cov<-input$jcpi_cov2
           res$model2$rnd_gen<-input$jcpi_gen2
           res$model2$rnd_non<-input$jcpi_rnd2
           res$dir<-path_jcpi
           res$phe_file<-input$jcpi_phe_upload$name
           res$gen_file<-input$jcpi_gen_upload$name
           res$email<-input$jcpi_email
           res$task<-paste0(task,"_jcpi")
           if(!is.null(input$jcpi_gen_upload)){
            res$num_of_snps<-sum(jcpi_gen_dat()$summary$No_var)
           }
           save(res,file=paste0(path_jcpi,"/par.RData"))
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
           system(paste0("nohup Rscript GCPI_job_submit.R /srv/shiny-server/" ,path_jcpi," > log/",task,".log 2>&1 &"))
           task<<-sample(1000:9999,1)
           path_jcpi<<-paste0("temp/",task,"_jcpi")
           #dir.create(path_jcpi,showWarnings = F)
        })

        observeEvent(input$cancel,{
          removeModal()
        })
}

             