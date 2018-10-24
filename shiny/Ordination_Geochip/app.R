# B7A Geochip Ordination App

library(shiny)
library(tidyverse)
library(vegan)
library(ape)
library(ggExtra)

#######################
#### Read in Files ####
#######################
geochip <- suppressWarnings(suppressMessages(read_tsv("Merged_geochip_B7A.tsv")))
metadata <- suppressWarnings(suppressMessages(read_tsv("B7A_metadata.tsv")))


######################
### Set up choices ###
######################

matched_choices <- list("All Samples" = "all_samples",
                        "Matched Samples" = "matched_samples")

visit_choices <- c("A", "B", "C")

geneCategory_choices <- c("All", unique(na.omit(geochip$Gene_category)))

ordination_choices <- list("PCA" = "PCA",
                           "DCA" = "DCA", 
                           "PCoA" = "PCoA")

disease_choices <- c("All", 
                     "No Disease", 
                     "Moderate / Severe Diarrhea")

phylum_choices <- geochip %>%
  filter(str_detect(Lineage, "Bacteria")) %>%
  filter(str_detect(Lineage, ";phylum")) %>%
  mutate(Phylum = gsub(x = Lineage, # Phylum column
                       pattern = ".*;phylum:(\\w*\\s*[-]*\\w*);.*", 
                       replacement = "\\1")) %>%
  select(Phylum) %>%
  distinct() %>%
  pull(Phylum)


color_choices <- list(
  "None" = "None",
  "Visit" = "Sample",
  "Moderate / Severe Diarrhea?" = "`Moderate-Severe diarrhea Yes/No`")

point_size_choices <- list(
  "None" = "None",
  "Inoculum dose (Target)" = "`Inoculum dose (Target)`",
  "Inoculum dose (Actual)" = "`Inoculum dose (Actual)`",
  "Disease Severity" = "`ETEC disease severity score`")



## UI ------------------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Geochip B7A Ordination"),
  sidebarLayout(
    sidebarPanel(
      
      # Matched or All samples
      radioButtons("matched", label = "Matched or All samples",
                   choices = matched_choices, inline = FALSE, selected = "matched_samples"),
      helpText(code("Matched "), "= Only patients that provided samples for all selected visits"),
      
      
      ################################################
      ### Inputs that select samples from metadata ###
      ################################################
      fluidRow(
        h3("Sample Selection:"),
        
        column(12, 
               wellPanel(
                 
                 # Visit
                 checkboxGroupInput('Visit_Letter', 'Visit:', 
                                    choices = visit_choices, 
                                    selected = c("A", "B", "C"), 
                                    inline = TRUE),
                 helpText(code("A:"), " Sample prior to innoculation"),
                 helpText(code("B:"), " Sample at acute symtoms or 5 days post innoculation, prior to antibiotics"),
                 helpText(code("C:"), " Convalescent Sample (~20-30 days post antibiotics)"),
                 
                 br(),
                 
                 # Disease Severity
                 selectInput("disease_severity", "Disease:", choices = disease_choices, selected = "All"),
                 helpText("Select patients that either did or did not develop moderate to severe diarrhea", br(),
                          code("All:"), " all paitients included in analysis")))),
      
      
      
      ##############
      ### PROBES ###
      ##############
      fluidRow(
        h3("Probe Selection"),
        
        column(12, 
               wellPanel(
                 
                 # Probe Functional Category
                 selectInput("geneCategory", "Probe Functional Category", choices = geneCategory_choices, selected = "All"),
                 helpText("Geochip Functional Category"), br(),
                 
                 # Select Phylum probes?
                 h5("Phylum Specific Probes"),
                 checkboxInput("select_phylum", label = "Select Phylum Probes?", value = FALSE),
                 helpText("If selected, only probes from selected bacterial phlya will be included in analysis"),
                 
                 # Phyla output
                 uiOutput("phyla"),
                 helpText("Select probes based on bacterial phyla")))),
      
      
      #################################
      ### Ordination and Aesthetics ###
      #################################
      fluidRow(
        h3("Ordination and Aesthetics"),
        
        column(12, 
               wellPanel(
                 # Ordination Type
                 radioButtons("ordination", label = h3("Ordination Type:"),
                              choices = ordination_choices, selected = "PCA"),
                 
                 # Color Choices (main choices)
                 radioButtons("color", label = h3("Color Points By:"),
                              choices = color_choices, inline = TRUE, selected = "Sample"),
                 
                 # Point Size
                 radioButtons("point_size", label = h3("Point Size"),
                              choices = point_size_choices, inline = TRUE, selected = "None"),
                 
                 # Connect points?
                 checkboxInput("lines", label = "Connect Points by study id", value = FALSE))),
        
        downloadButton('downloadPlot','Download Plot')),
      
      
      
      
      #sidebar width
      width = 4),
    
    # Plot
    mainPanel(
      plotOutput("plot", width = "800px", height = "800px"),
      
      # Table to see patients (not needed, but useful for troubleshooting)
      fluidRow(column(12,tableOutput('table')))
    )))
      








## Server ---------------------------------------------------------------------------------
server <- function(input, output){
  
  ##############################################
  ### Filter metadata for samples in Geochip ###
  ##############################################
  
  metadata <- metadata %>%
    filter(glomics_ID %in% colnames(geochip))

  ##############################################
  ### Only keep patients with matched visits ###
  ##############################################

  geochip_matched <- reactive({
    if(input$matched == "matched_samples"){
      # Get vector of patients that have matched visits
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
      
    # Subset the geochip data
    geochip %>%
      select_if(colnames(.) %in% c("Genbank ID", "Gene", "Organism", 
                                   "Gene_category", "Subcategory1",
                                   "Subcategory2", "Lineage", matched_samples))
    
  })
  
  ##############################
  ### Subset by Disease Type ###
  ##############################
  
  geochip_disease <- reactive({
    if(input$disease_severity == "Moderate / Severe Diarrhea"){
      disease_samples <- metadata %>%
        filter(`Moderate-Severe diarrhea Yes/No` == "Yes") %>%
        pull(glomics_ID)
    } else if (input$disease_severity == "No Disease"){
      disease_samples <- metadata %>%
        filter(`Moderate-Severe diarrhea Yes/No` == "No") %>%
        pull(glomics_ID)
    } else {
      disease_samples <- metadata %>%
        pull(glomics_ID)
    }
    
    # Subset the geochip data
    geochip_matched() %>%
      select_if(colnames(.) %in% c("Genbank ID", "Gene", "Organism", 
                                   "Gene_category", "Subcategory1",
                                   "Subcategory2", "Lineage", disease_samples))
  })
  
  
  #######################################
  ### Filter by Probe Functional Type ###
  #######################################
  
  geochip_probe <- reactive({
    if(input$geneCategory != "All"){
      geochip_disease() %>%
        filter(Gene_category == input$geneCategory)
    } else if (input$geneCategory == "All"){
      geochip_disease() 
    } else {
      stopApp()
    }
  })
  
  
  
  ########################
  ### Filter by Phylum ###
  ########################
  
  geochip_phylum <- reactive({
      
    if(input$select_phylum){
      geochip_probe() %>%
        filter(!is.na(Lineage)) %>%
        filter(str_detect(Lineage, "Bacteria")) %>%
        filter(str_detect(Lineage, ";phylum")) %>%
        mutate(Phylum = gsub(x = Lineage, # Phylum column
                             pattern = ".*;phylum:(\\w*\\s*\\w*);.*", 
                             replacement = "\\1")) %>%
        filter(Phylum %in% input$phylum)
      
    } else if (!input$select_phylum) {
      geochip_probe()
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
  
  
  
  ##########################################
  ### Convert geochip to numeric matrix ###
  ##########################################
  
  # Filer geochip to only include specified Visit_Number samples
  geochip_matrix <- reactive({
    geochip1 <- geochip_phylum() %>%
      select(starts_with("B7A"))
  
    # Set NA's to 0 and values not NA to original value
    geochip1 <- geochip1 %>%
      mutate_all(funs(ifelse(is.na(.), 0, .)))
    
    # Remove rows that equal 0
    geochip1 <- geochip1[rowSums(geochip1) != 0,]
    
    as.matrix(geochip1)
  })
 
  
  
  #############################
  #### Ordination Analysis ####
  #############################
  
  # Ordination choices
  geochip_ord <- reactive({
    # Perform PCA analysis using vegan
    if(input$ordination == "PCA"){
      vegan::rda(t(geochip_matrix()))
      
      # Perform DCA analysis using vegan
    } else if (input$ordination == "DCA") {
      vegan::decorana(t(geochip_matrix()))
      
      # Perform PCoA using ape package
    } else if (input$ordination == "PCoA"){
      
      geochip_dist <- vegan::vegdist(as.matrix(t(geochip_matrix())))
      ape::pcoa(geochip_dist, correction = "none")
    } else {
      stopApp()
    }
  })
  
  
  geochip_coords <- reactive({
    
    # PCA and DCA coordinates
    if(input$ordination == "PCA" | input$ordination == "DCA"){
      
      # Get the coordinated for the samples
      ord_coords <- scores(geochip_ord(), display = "sites")
      
      # Make tibble
      ord_coords <- as.data.frame(ord_coords) %>%
        rownames_to_column(var = "glomics_ID")
      
      ord_coords
      
      # PCoA Coords
    } else if (input$ordination == "PCoA") {
      
      ord_coords <- geochip_ord()$vectors[,1:2]
      
      # Make tibble
      ord_coords <- as.data.frame(ord_coords) %>%
        rownames_to_column(var = "glomics_ID")
      
      ord_coords
    } else {
      stopApp()
    }
  })
  
  
  geochip_prop_exp <- reactive({
    
    # Extract proportion explained by first couple PC's for PCA
    if (input$ordination == "PCA"){
      summary(eigenvals(geochip_ord()))[2,] * 100
    } else if (input$ordination == "PCoA"){
      geochip_ord()$values$Relative_eig * 100
    } else {
      stopApp()
    }
  })
  
  
  ####################################
  #### Merge coords into Metadata ####
  ####################################
  
  geochip_coords_metadata <- reactive({
    
    geochip_coords() %>%
      left_join(., metadata, by = c("glomics_ID")) %>%
      
      # Dummy column for Null coloring
      mutate(dummy = "None")
      
      # Factor Columns
     # mutate("ETEC disease severity score" = factor(`ETEC disease severity score`))
  }) 
  
  
  
  
  ##############################
  #### Create plots/tables #####
  ##############################
  
  plotInput <- reactive({
    
    # Get color input for points
    my_fill <- ifelse(input$color == "None", "dummy", input$color)

    
    # Get size input for points
    my_size <- ifelse(input$point_size == "None", "4", input$point_size)
    
    
    # Aesthetic sizes
    axis_title_size <- 16
    axis_text_size <- 14
    title_size <- 18
    legend_text_size <- 12
    
    
    # PCA Plot
    if(input$ordination == "PCA"){
      pca_plot <- ggplot(geochip_coords_metadata(),
                         aes_string(x = "PC1", 
                                    y = "PC2", 
                                    color = my_fill,
                                    size = my_size)) +
        
        # Set up proportion explained
        xlab(paste0("PC1(", round(geochip_prop_exp()[[1]], 2), "%)")) +
        ylab(paste0("PC2(", round(geochip_prop_exp()[[2]], 2), "%)")) +
        geom_point(pch = 1, alpha = 1) +
        geom_point(pch = 19, alpha = 0.8) +
        ggtitle("PCA Analysis") +
        theme_minimal() +
        theme(
          axis.title.x = element_text(size = axis_title_size),
          axis.text.x = element_text(size = axis_text_size, hjust = 1),
          axis.text.y = element_text(size = axis_text_size),
          axis.title.y = element_text(size = axis_title_size),
          plot.title = element_text(size = title_size, face = "bold"),
          legend.text = element_text(size = legend_text_size),
          legend.title = element_blank()) +
        guides(fill = guide_legend(override.aes = list(size=7))) # make legend points larger
      
      # Connect lines?
      if (input$lines == TRUE){
        pca_plot <- pca_plot + geom_line(aes(group = `Subject ID`), linetype = 1, color = "black", size = 0.8)
      }
      
      ggMarginal(pca_plot, groupColour = TRUE, groupFill = TRUE)
      
      
    } else if (input$ordination == "DCA"){
      
      
      # DCA plot
      dca_plot <- ggplot( geochip_coords_metadata(),
                         aes_string(x = "DCA1",
                                    y = "DCA2", 
                                    color = my_fill,
                                    size = my_size)) +
        xlab("DCA1") +
        ylab("DCA2") +
        geom_point(pch = 1, alpha = 1) +
        geom_point(pch = 19, alpha = 0.8) +
        ggtitle("DCA Analysis") +
        theme_minimal() +
        theme(
          axis.title.x = element_text(size = axis_title_size),
          axis.text.x = element_text(size = axis_text_size, hjust = 1),
          axis.text.y = element_text(size = axis_text_size),
          axis.title.y = element_text(size = axis_title_size),
          plot.title = element_text(size = title_size, face = "bold"),
          legend.text = element_text(size = legend_text_size),
          legend.title = element_blank()) +
        guides(fill = guide_legend(override.aes = list(size=7))) # make legend points larger
      
      # Connect lines?
      if (input$lines == TRUE){
        dca_plot <- dca_plot + geom_line(aes(group = `Subject ID`), linetype = 1, color = "black", size = 0.8)
      }
      
      ggMarginal(dca_plot, groupColour = TRUE, groupFill = TRUE)
      
      
    } else if (input$ordination == "PCoA"){
      
      # PCoA Plot
      pcoa_plot <- ggplot( geochip_coords_metadata(), 
                          aes_string(x = "Axis.1", 
                                     y = "Axis.2", 
                                     color = my_fill, 
                                     size = my_size)) +
        scale_fill_discrete(name = "Visit Number") +
        xlab(paste0("PCoA1(", round(geochip_prop_exp()[1], 2), "%)")) +
        ylab(paste0("PCoA2(", round(geochip_prop_exp()[2], 2), "%)")) +
        geom_point(pch = 1, alpha = 1) +
        geom_point(pch = 19, alpha = 0.8) + 
        ggtitle("PCoA Analysis") + 
        theme_minimal() +
        theme(
          axis.title.x = element_text(size = axis_title_size), 
          axis.text.x = element_text(size = axis_text_size, hjust = 1),
          axis.text.y = element_text(size = axis_text_size),
          axis.title.y = element_text(size = axis_title_size), 
          plot.title = element_text(size = title_size, face = "bold"),
          legend.text = element_text(size = legend_text_size),
          legend.title = element_blank()) +
        guides(fill = guide_legend(override.aes = list(size=5))) # make legend points larger
      
      # Connect lines?
      if (input$lines == TRUE){
        pcoa_plot <- pcoa_plot + geom_line(aes(group = `Subject ID`), linetype = 1, color = "black", size = 0.8)
      }
      
      ggMarginal(pcoa_plot, groupColour = TRUE, groupFill = TRUE)
    }
  })
  
  
  ####################
  ### Display Plot ###
  ####################
  output$plot <- renderPlot({
    print(plotInput())
  })
  
  output$table <- renderTable({geochip_coords_metadata()})
  
  #####################
  ### Download plot ###
  #####################
  
  output$downloadPlot <- downloadHandler(
    filename = function(){paste("shiny_plot",'.png',sep='')},
    content = function(file){
      ggsave(file, plot=plotInput())
    })
  
  
}







## Run App --------------------------------------------------------------------------------
shinyApp(server = server, ui = ui)