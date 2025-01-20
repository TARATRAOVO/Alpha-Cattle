# Add at the beginning of the script
source("vcfTools.R") # Replace with actual path

library(shiny)
library(vcfR)
library(biomaRt)
library(dplyr)
library(plotly)

# Add system requirement check function at the beginning
check_system_requirements <- function() {
  # Check if bgzip is installed
  bgzip_check <- system("which bgzip", ignore.stderr = TRUE)
  if (bgzip_check != 0) {
    stop("bgzip not found. Please install htslib package.")
  }

  # Check if tabix is installed
  tabix_check <- system("which tabix", ignore.stderr = TRUE)
  if (tabix_check != 0) {
    stop("tabix not found. Please install htslib package.")
  }
}

# Add function to handle uploaded files
prepare_vcf_file <- function(input_path, is_gz = FALSE) {
  # Create temporary directory
  temp_dir <- tempdir()

  # Construct output file path
  output_gz <- file.path(temp_dir, "uploaded_vcf.vcf.gz")

  tryCatch(
    {
      if (!is_gz) {
        # If input file is not gz format, compress using bgzip
        system2("bgzip",
          args = c("-c", input_path),
          stdout = output_gz
        )
      } else {
        # If already in gz format, just copy
        file.copy(input_path, output_gz, overwrite = TRUE)
      }

      # Create tabix index
      system2("tabix",
        args = c("-p", "vcf", output_gz)
      )

      return(output_gz)
    },
    error = function(e) {
      stop("Error processing uploaded file: ", e$message)
    }
  )
}

magServer <- function(input, output, session) {
  # Read SNP information file and ensure correct column names
  snp_info <- read.csv("data_base/meg_snp_info/snp_info.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Extract unique traits and ensure not empty
  unique_traits <- unique(snp_info$`Variant Phenotype`)

  # Pass traits to UI, add validation
  output$trait_choices <- renderUI({
    if (length(unique_traits) > 0) {
      selectInput("selected_trait",
        "Select Analysis Trait:",
        choices = unique_traits,
        selected = unique_traits[1]
      )
    } else {
      div(
        style = "color: red;",
        "Error: Unable to load trait data. Please check data file."
      )
    }
  })

  # Dynamically display upload or file selection UI based on data source
  output$upload_ui <- renderUI({
    if (input$data_source == "Upload Data") {
      fileInput("vcf_file", "Select VCF File:", accept = c(".vcf", ".vcf.gz"))
    }
  })

  # Add: Display related SNP information immediately when trait is selected
  observeEvent(input$selected_trait, {
    req(input$selected_trait)

    # Get SNP information for selected trait
    selected_snps <- dplyr::filter(snp_info, `Variant Phenotype` == input$selected_trait)

    # Create display dataframe
    snp_info_display <- data.frame(
      Phenotype = selected_snps$`Variant Phenotype`,
      Chromosome = selected_snps$Chr.,
      PhysicalPosition = selected_snps$start,
      Gene = selected_snps$Gene,
      VariantType = selected_snps$`Type of Variant`,
      VariantDescription = selected_snps$variant_description,
      Mutation = selected_snps$mutation
    )

    # Display SNP information table
    output$selected_snp_info <- renderTable({
      snp_info_display
    })
  })

  # Read VCF file and display
  observeEvent(input$load_vcf, {
    req(input$data_source, input$selected_trait)

    # Check system requirements
    check_system_requirements()

    withProgress(message = "Processing data...", value = 0, {
      incProgress(0.1, detail = "Preparing to load VCF file...")

      # Process VCF file based on data source
      tryCatch({
        if (input$data_source == "Upload Data" && !is.null(input$vcf_file)) {
          # Process uploaded file
          showNotification("Processing uploaded VCF file...", type = "message")

          # Check if file is gz format
          is_gz <- grepl("\\.gz$", input$vcf_file$name)

          # Prepare VCF file (compress and index)
          vcf_path <- prepare_vcf_file(input$vcf_file$datapath, is_gz)

          showNotification("VCF file processing complete", type = "message")
        } else if (input$data_source == "Public Database Data") {
          # Use VCF file from public database
          selected_snps <- dplyr::filter(snp_info, `Variant Phenotype` == input$selected_trait)
          chr_number <- unique(selected_snps$Chr.)

          vcf_path <- file.path(
            "/srv/shiny-server/data_base/genomic_data",
            sprintf("samples-Chr%s-Run8-TAUIND-public_miss05_snponly_nomes_phased.vcf.gz", chr_number)
          )

          if (!file.exists(vcf_path)) {
            stop(sprintf("VCF file for chromosome %s not found.", chr_number))
          }
        } else {
          stop("Please select valid data source and file.")
        }

        incProgress(0.4, detail = "Filtering relevant SNPs...")

        # Filter relevant SNPs based on selected trait - modify column names to match new format
        selected_trait <- input$selected_trait
        relevant_snps <- dplyr::filter(snp_info, `Variant Phenotype` == selected_trait) %>%
          dplyr::transmute(
            SNP_Name = paste(Chr., start, mutation, sep = "_"),
            Chromosome = Chr.,
            Position = start,
            Mutation = mutation,
            Gene = Gene
          )

        incProgress(0.6, detail = "Calculating REF ratio...")

        # Get VCF file path
        if (input$data_source == "Upload Data" && !is.null(input$vcf_file)) {
          vcf_path <- input$vcf_file$datapath
        } else {
          # Use public database VCF file path
          chr_number <- unique(selected_snps$Chr.)
          vcf_path <- file.path(
            "/srv/shiny-server/data_base/genomic_data",
            sprintf("samples-Chr%s-Run8-TAUIND-public_miss05_snponly_nomes_phased.vcf.gz", chr_number)
          )
        }

        # Calculate REF ratio for each SNP position
        ref_percent <- sapply(1:nrow(relevant_snps), function(i) {
          chr <- relevant_snps$Chromosome[i]
          pos <- relevant_snps$Position[i]

          result <- query_ref_percent(vcf_path, chr, pos)

          showNotification(
            sprintf(
              "Position %s:%s query result: %s",
              chr, pos,
              if (is.na(result)) "Not found" else sprintf("%.2f%%", result)
            ),
            type = if (is.na(result)) "warning" else "message"
          )

          return(result)
        })

        incProgress(0.8, detail = "Preparing results...")

        # Generate pie charts after calculating ref_percent
        output$allele_plots <- renderUI({
          plot_output_list <- lapply(1:nrow(relevant_snps), function(i) {
            plot_id <- paste0("plot_", i)

            # Extract position information
            chr <- relevant_snps$Chromosome[i]
            pos <- relevant_snps$Position[i]
            gene <- relevant_snps$Gene[i]
            ref_pct <- ref_percent[i]

            if (!is.na(ref_pct)) {
              alt_pct <- 100 - ref_pct

              # Create pie chart data
              plot_data <- data.frame(
                type = c("REF", "ALT"),
                value = c(ref_pct, alt_pct)
              )

              # Create interactive pie chart using plotly
              output[[plot_id]] <- renderPlotly({
                plot_ly(plot_data,
                  labels = ~type, values = ~value, type = "pie",
                  textinfo = "label+percent",
                  hoverinfo = "text",
                  text = ~ paste(type, sprintf("%.2f%%", value)),
                  marker = list(colors = c("#1f77b4", "#ff7f0e"))
                ) %>%
                  layout(
                    title = list(
                      text = sprintf(
                        "Position %s:%s (%s) Allele Distribution",
                        chr, pos, gene
                      ),
                      y = 0.9 # Move title down
                    ),
                    margin = list(t = 100), # Increase top margin
                    showlegend = TRUE
                  )
              })

              # Return plotly output container
              div(
                style = "margin-bottom: 20px;",
                plotlyOutput(plot_id, height = "400px")
              )
            } else {
              # Display message if no data found
              div(
                style = "margin-bottom: 20px; color: red;",
                sprintf("No data found for position %s:%s (%s)", chr, pos, gene)
              )
            }
          })

          # Combine all charts
          do.call(tagList, plot_output_list)
        })

        # Create result dataframe - modify to include more information
        snp_ref_df <- data.frame(
          Phenotype = selected_snps$`Variant Phenotype`,
          Chromosome = relevant_snps$Chromosome,
          PhysicalPosition = relevant_snps$Position,
          Gene = relevant_snps$Gene,
          VariantType = selected_snps$`Type of Variant`,
          Mutation = relevant_snps$Mutation,
          REFAlleleRatio = sapply(ref_percent, function(x) {
            if (is.na(x)) {
              "No related SNP found"
            } else if (x == 0) {
              "0.00%"
            } else {
              sprintf("%.2f%%", x)
            }
          })
        )

        # Display results
        output$snp_ref_display <- renderTable({
          snp_ref_df
        })

        incProgress(1, detail = "Complete!")

        showNotification("Data processing complete.", type = "message")
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
        return(NULL)
      }, finally = {
        # Clean up temporary files if using uploaded file
        if (input$data_source == "Upload Data") {
          temp_dir <- tempdir()
          unlink(file.path(temp_dir, "uploaded_vcf.vcf.gz*"))
        }
      })
    })
  })

  # Other existing server logic...
}
