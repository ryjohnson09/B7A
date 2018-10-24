# B7A Humichip Response Ratio App

library(shiny)
library(tidyverse)

#######################
#### Read in Files ####
#######################

humichip <- suppressWarnings(suppressMessages(read_tsv("Merged_humichip_B7A.tsv")))
metadata <- suppressWarnings(suppressMessages(read_tsv("B7A_metadata.tsv")))


######################
### Set up choices ###
######################

matched_choices <- list("All Samples" = "all_samples",
                        "Matched Samples" = "matched_samples")

visit_choices <- c("A", "B", "C")

diarrhea_choices <- c("Yes", "No", "Both")

compare_choices <- c("Sample", "Moderate-Severe diarrhea Yes/No")

cat_choices <- list("Gene Category" = "geneCategory",
                    "Subcategory 1" = "subcategory1",
                    "Subcategory 2" = "subcategory2",
                    "Phylum" = "Phylum")

phylum_choices <- humichip %>%
  filter(str_detect(lineage, "Bacteria")) %>%
  filter(str_detect(lineage, ";phylum")) %>%
  mutate(Phylum = gsub(x = lineage, # Phylum column
                       pattern = ".*;phylum:(\\w*\\s*[-]*\\w*);.*", 
                       replacement = "\\1")) %>%
  select(Phylum) %>%
  distinct() %>%
  pull(Phylum)





## UI ------------------------------------------------------------------------------------

ui <- fluidPage(
  
  titlePanel("B7A Humichip Response Ratios"),
  sidebarLayout(
    sidebarPanel(
      
      # Matched or All samples
      radioButtons("matched", label = "Matched or All samples",
                   choices = matched_choices, inline = TRUE, selected = "matched_samples"),
      helpText("Matched = Only patients that provided samples for all selected visits"),
      
      ##########################################
      ### Inputs that Select Certain Samples ###
      ##########################################
      fluidRow(
        h3("Sample Selection:"),
        
        column(12, 
               wellPanel(
                 
                 # Filter by Visit
                 checkboxGroupInput('Visit_Letter', 'Visit:', 
                                    choices = visit_choices, 
                                    selected = c("A", "B"), 
                                    inline = TRUE),
                 helpText(code("A:"), " Sample prior to innoculation"),
                 helpText(code("B:"), " Sample at acute symtoms or 5 days post innoculation, prior to antibiotics"),
                 helpText(code("C:"), " Convalescent Sample (~20-30 days post antibiotics)"),
                 
                 # Filter by Diarrhea
                 radioButtons('diarrhea_choice', 'Moderate to Severe Diarrhea?', 
                              choices = diarrhea_choices, 
                              selected = c("Both"), 
                              inline = TRUE),
                 
                 
                 
                 # Comparing Groups
                 radioButtons('compare_groups', 'Groups to Compare:', choices = compare_choices, selected = "Sample", inline = TRUE),
                 helpText("Select patient groups to compare by visit")))),
      
      
      ########################
      ### Y axis categories ##
      ########################
      fluidRow(
        h3("Group Selections"),
        
        column(12, 
               wellPanel(
                 # geneCategory vs subcategory1 vs subcategory2
                 selectInput("cat_choice", "Y-axis Categories", choices = cat_choices, selected = "geneCategory"),
                 helpText("Select Categories to Compare"),
                 
                 # Number of minimum samples that must have Category
                 sliderInput("cat_min_number", "Minimum Observation in Category:",
                             min = 1, max = 20,
                             value = 10)))),
      
      #######################
      ### Probe Selection ###
      #######################
      fluidRow(
        h3("Probe Selection"),
        
        column(12, 
               wellPanel(
                 
                 # Select Phylum probes?
                 h5("Phylum Specific Probes"),
                 checkboxInput("select_phylum", label = "Select Phylum Probes?", value = FALSE),
                 helpText("If selected, only probes from selected bacterial phlya will be included in analysis"),
                 
                 # Phyla output
                 uiOutput("phyla"),
                 helpText("Select probes based on bacterial phyla")))),
      
      
      fluidRow(
        h3("Ordination and Aesthetics"),
        
        column(12, 
               wellPanel(
                 # Only keep significant RR's
                 checkboxInput("keeper", label = "Remove non-significant RR's", value = FALSE)
               )),
        
        
        downloadButton('downloadPlot','Download Plot')),
      
      
      #sidebar width
      width = 4),
    
    # Plot
    mainPanel(
      plotOutput("plot", width = "800px", height = "800px"),
      br()
      
      # Table to see patients (not needed, but useful for troubleshooting)
      #fluidRow(column(12,tableOutput('table')))
    )))



## Server --------------------------------------------------------------------------------

server <- function(input, output){
  
  ##############################################
  ### Filter metadata for samples in Humichip ###
  ##############################################
  
  metadata <- metadata %>%
    filter(glomics_ID %in% colnames(humichip))
  
  ##############################################
  ### Only keep patients with matched visits ###
  ##############################################
  
  humichip_matched <- reactive({
    
    # Get vector of patients that have matched visits
    if(input$matched == "matched_samples"){
      matched_samples <- metadata %>%
        filter(Sample %in% input$Visit_Letter) %>%
        filter(!is.na(`Sample Date`)) %>%
        group_by(`Subject ID`) %>%
        mutate(n_samples = n()) %>%
        filter(n_samples == length(input$Visit_Letter)) %>%
        pull(glomics_ID)
      
    } else {
      matched_samples <- metadata %>%
        filter(Sample %in% input$Visit_Letter) %>%
        filter(!is.na(`Sample Date`)) %>%
        pull(glomics_ID)
    }
    
    # Subset the humichip data
    humichip %>%
      select_if(colnames(.) %in% c("Genbank ID", "gene", "species", 
                                   "geneCategory", "subcategory1",
                                   "subcategory2", "lineage", "annotation", matched_samples))
  })
  
  
  
  
  ##############################
  ### Subset by Disease Type ###
  ##############################

  humichip_disease <- reactive({

    if(input$diarrhea_choice == "Yes"){
      disease_samples <- metadata %>%
        filter(`Moderate-Severe diarrhea Yes/No` == "Yes") %>%
        pull(glomics_ID)
    } else if (input$diarrhea_choice == "No"){
      disease_samples <- metadata %>%
        filter(`Moderate-Severe diarrhea Yes/No` == "No") %>%
        pull(glomics_ID)
    } else if (input$diarrhea_choice == "Both"){
      disease_samples <- metadata %>%
        pull(glomics_ID)
    } else {
      stopApp("Problem when filtering by diarrhea severity")
    }

    # Subset the humichip data
    humichip_matched() %>%
      select_if(colnames(.) %in% c("Genbank ID", "gene", "species",
                                   "geneCategory", "subcategory1",
                                   "subcategory2", "lineage", "annotation",
                                   disease_samples))
  })
  
  ########################
  ### Filter by Phylum ###
  ########################
  
  humichip_phylum <- reactive({
    
    # If filtering by phylum, and phylum column (only bacteria with phylum designation)
    if(input$select_phylum){
      humichip_disease() %>%
        filter(str_detect(lineage, "Bacteria")) %>%
        filter(str_detect(lineage, ";phylum")) %>%
        mutate(Phylum = gsub(x = lineage, # Phylum column
                             pattern = ".*;phylum:(\\w*\\s*[-]*\\w*);.*", 
                             replacement = "\\1")) %>%
        filter(Phylum %in% input$phylum) # Filter step
      
    } else if (!input$select_phylum) {
      humichip_disease() %>%
        filter(str_detect(lineage, "Bacteria")) %>%
        filter(str_detect(lineage, ";phylum")) %>%
        mutate(Phylum = gsub(x = lineage, # Phylum column
                             pattern = ".*;phylum:(\\w*\\s*[-]*\\w*);.*", 
                             replacement = "\\1"))
    } else {
      stopApp()
    }
  })
  
  #########################################
  ### Render list of Phyla if applicable ##
  #########################################
  output$phyla <- renderUI({
    
    if(input$select_phylum){
      checkboxGroupInput("phylum", "Select Phyla:",
                         choices = phylum_choices, inline = TRUE)
    }
  })

  ###############################################
  ### Select for Functional or STR_SPE probes ###
  ###############################################

  humichip_probe <- reactive({

    # If cat_choice == Phylum, filter for STR_SPE
    # and add Phylum category
    if(input$cat_choice == "Phylum"){
      humichip_phylum() %>%
        filter(gene == "STR_SPE")
      
      # If cat_choice == Functional group
    } else if (input$cat_choice != "Phylum"){
      humichip_phylum() %>%
        filter(gene != "STR_SPE")
    } else {
      stopApp("Problem filtering STR_SPE vs Functional Probes")
    }
  })



  ############################################
  ### Filter humichip for certain patients ###
  ### And Make Long ##########################
  ############################################


  humichip_final <- reactive({

    humichip_probe() %>%

      # Make Long
      gather(key = glomics_ID,
             value = Signal,
             -`Genbank ID`, -gene,
             -species, -lineage, -annotation, -geneCategory, -subcategory1,
             -subcategory2, -Phylum) %>%
      # Merge in ID decoder
      left_join(., metadata, by = c("glomics_ID"))
  })


  #################################
  ### Calculate Rel Abundance #####
  #################################
  humichip_rel_abun <- reactive({

    y_axis_cat <- sym(input$cat_choice)
    grouping_cat <- sym(input$compare_groups)

    humichip_final() %>%

      group_by(glomics_ID) %>%
      mutate(Signal_Relative_Abundance = (Signal / sum(Signal, na.rm = TRUE)* 100)) %>%

      # Remove columns not needed
      select(glomics_ID, !!y_axis_cat, Signal_Relative_Abundance) %>%


      # # Remove any rows with NA in the signal category or y_axis_cat
      filter(!is.na(Signal_Relative_Abundance)) %>%
      filter(!is.na(!!y_axis_cat)) %>%


      # Calculate relative abundance
      ungroup() %>%
      group_by(!!y_axis_cat, glomics_ID) %>%
      summarise(cat_relative_abundance = sum(Signal_Relative_Abundance, na.rm = TRUE)) %>%

      # Merge in metadata
      left_join(., metadata, by = c("glomics_ID")) %>%
      select(!!y_axis_cat, glomics_ID, Sample, cat_relative_abundance, !!grouping_cat) %>%
      filter(!is.na(!!grouping_cat)) %>%
      ungroup()
  })



  ################################
  ### Calculate Response Ratio ###
  ################################



  humichip_RR <- reactive({

    y_axis_cat <- sym(input$cat_choice)
    grouping_cat <- sym(input$compare_groups)

    humichip_rel_abun() %>%

      # Calculate mean, sd, and n
      group_by(!!y_axis_cat, !!grouping_cat) %>%
      summarise(mean_signal = mean(cat_relative_abundance),
                sd_signal = sd(cat_relative_abundance),
                n = sum(!is.na(cat_relative_abundance))) %>%


      # Spread the signal mean by visit number
      ungroup() %>%
      group_by(!!y_axis_cat) %>%
      spread(!!grouping_cat, mean_signal) %>%

      # Rename mean columns
      rename(group1_mean = colnames(.)[length(colnames(.)) - 1],
             group2_mean = colnames(.)[length(colnames(.))]) %>%

      # Spread the sd and n columns by visit
      mutate(sd_group1 = ifelse(!is.na(group1_mean), sd_signal, NA)) %>%
      mutate(sd_group2 = ifelse(!is.na(group2_mean), sd_signal, NA)) %>%
      mutate(n_group1 = ifelse(!is.na(group1_mean), n, NA)) %>%
      mutate(n_group2 = ifelse(!is.na(group2_mean), n, NA)) %>%
      select(-sd_signal, -n) %>%

      # Compress NAs
      ungroup() %>%
      group_by(!!y_axis_cat) %>%
      summarise_all(funs(sum(., na.rm = T))) %>%

      # Must have at least __ observations in each subcategory
      filter(n_group1 >= input$cat_min_number) %>%
      filter(n_group2 >= input$cat_min_number) %>%

      # Calculate SEM for each mean
      mutate(SEM_group1 = sd_group1 / sqrt(n_group1)) %>%
      mutate(SEM_group2 = sd_group2 / sqrt(n_group2)) %>%

      # Calculate the Response Ratio (RR)
      mutate(RR = log(group2_mean / group1_mean)) %>%

      # Calculate the Standard error for the RR
      mutate(SE_RR = sqrt((SEM_group1**2 / group1_mean**2) + (SEM_group2**2 / group2_mean**2))) %>%

      # Calcualte the 95% confidence interval for each RR
      mutate(CI95 = abs(1.96 * SE_RR)) %>%

      # Add in keeper column if does not overlap 0
      mutate(keeper = ifelse(0 > (RR - CI95) & 0 < (RR + CI95), "No", "Yes")) %>%

      # Make labels pretty
      mutate(pretty_cat = str_to_title(!!y_axis_cat)) %>%
      mutate(pretty_cat = str_replace_all(pretty_cat, "_"," ")) %>%

      # Factor columns
      mutate(pretty_cat = fct_reorder(pretty_cat, RR)) %>%
      ungroup()
  })


  ##########################################
  ### Filter for RR that don't overlap 0 ###
  ##########################################

  humichip_significant <- reactive({

    if (input$keeper){
      humichip_RR() %>%
        filter(keeper == "Yes")
    } else if (!input$keeper){
      humichip_RR()
    } else {
      stopApp("Failure at 'keeper' column creation")
    }
  })



  ##############################
  #### Create plots/tables #####
  ##############################

  plotInput <- reactive({

    y_axis_cat <- sym(input$cat_choice)
    grouping_cat <- sym(input$compare_groups)

    group_labels <- humichip_final() %>%
      select(!!grouping_cat) %>%
      distinct() %>%
      pull(!!grouping_cat)


    ggplot(data = humichip_significant()) +
      geom_hline(yintercept = 0, linetype = "dashed", size = 1) +

      # points and error bar
      geom_point(aes(x = pretty_cat, y = RR), size = 4) +
      geom_errorbar(aes(ymin = RR - CI95,
                        ymax = RR + CI95,
                        x = pretty_cat),
                    width = 0.25) +

      # Group labels
      annotate(geom = "text", label = group_labels[1], x = Inf, y = -Inf, hjust = 0, vjust = 1,
               size = 5, color = "red", fontface = 2) +
      annotate(geom = "text", label = group_labels[2], x = Inf, y = Inf, hjust = 1, vjust = 1,
               size = 5, color = "red", fontface = 2) +

      # plot labels
      labs(title = "Response Ratio",
           x = "Category",
           y = "Response Ratio") +

      theme_minimal() +
      coord_flip() +
      theme(
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5)
      )

  })



  ####################
  ### Display Plot ###
  ####################
  output$plot <- renderPlot({
    print(plotInput())
  })

  #####################
  ### Download plot ###
  #####################

  output$downloadPlot <- downloadHandler(
    filename = function(){paste("shiny_plot",'.png',sep='')},
    content = function(file){
      ggsave(file, plot=plotInput())
    })
  
  
  
  
  
  
  #output$table <- renderTable({head(humichip_matched(), 100)})
  
}







shinyApp(server = server, ui = ui)