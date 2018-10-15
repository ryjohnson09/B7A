# B7A Humichip Ordination App

library(shiny)
library(tidyverse)
library(vegan)
library(ape)
library(ggExtra)

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

probe_choices <- c("Functional", "Strain/Species", "All")

geneCategory_choices <- c("All", unique(na.omit(humichip$geneCategory)))

ordination_choices <- list("PCA" = "PCA",
                           "DCA" = "DCA", 
                           "PCoA" = "PCoA")

color_choices <- list(
  "None" = "None",
  "Visit" = "Sample",
  "Moderate / Severe Diarrhea?" = "`Moderate-Severe diarrhea Yes/No`",
  "Disease Severity" = "`ETEC disease severity score`")

point_size_choices <- list(
  "None" = "None",
  "Inoculum dose (Target)" = "`Inoculum dose (Target)`",
  "Inoculum dose (Actual)" = "`Inoculum dose (Actual)`")



## UI ------------------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Humichip B7A Ordination"),
  sidebarLayout(
    sidebarPanel(
      
      # Matched or All samples
      radioButtons("matched", label = "Matched or All samples",
                   choices = matched_choices, inline = FALSE, selected = "matched_samples"),
      helpText("Matched = Only patients that provided samples for all selected visits"),
      
      
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
                 helpText("A: Sample prior to innoculation"),
                 helpText("B: Sample at acute symtoms or 5 days post innoculation, prior to antibiotics"),
                 helpText("C: Convalescent Sample (~20-30 days post antibiotics)")))),
      
      
      
      ##############
      ### PROBES ###
      ##############
      fluidRow(
        h3("Probe Selection"),
        
        column(12, 
               wellPanel(
                 # Probe Type
                 selectInput('probe', 'Probe Type:', choices = probe_choices, selected = "All"),
                 helpText(code("Functional:"), " probes specific for functional genes", br(),
                          code("Strain/Species:"), " probes specific for microbial species and strains", br(),
                          code("All:"), " all probes included in analysis"),
                 br(),
                 # Probe Functional Category
                 selectInput("geneCategory", "Probe Functional Category", choices = geneCategory_choices, selected = "All"),
                 helpText("If ", code("Probe Type "), "= ", strong("Functional"), 
                          ", can select by functional group")))),
      
      
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
      plotOutput("plot", width = "800px", height = "800px")
      
      # Table to see patients (not needed, but useful for troubleshooting)
      #fluidRow(column(12,tableOutput('table'))),
    )))
      








## Server ---------------------------------------------------------------------------------
server <- function(input, output){

##############################################
### Only keep patients with matched visits ###
##############################################

  humichip_matched <- reactive({
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
      
    # Subset the humichip data
    humichip %>%
      select_if(colnames(.) %in% c("Genbank ID", "gene", "species", "lineage",
                                   "annotation", "geneCategory", "subcategory1",
                                   "subcategory2", matched_samples))
    
  })
  
  
  #############################
  ### Functional or STR_SPE ###
  #############################
  
  humichip_probe <- reactive({
    if(input$probe == "Functional" & input$geneCategory == "All"){
      humichip_matched() %>%
        filter(gene != "STR_SPE")
    } else if (input$probe == "Functional" & input$geneCategory != "All"){
      humichip_matched() %>%
        filter(gene != "STR_SPE") %>%
        filter(geneCategory == input$geneCategory)
    } else if (input$probe == "Strain/Species"){
      humichip_matched() %>%
        filter(gene == "STR_SPE")
    } else {
      humichip_matched()
    }
  })

  
  ##########################################
  ### Convert Humichip to numeric matrix ###
  ##########################################
  
  # Filer Humichip to only include specified Visit_Number samples
  humichip_matrix <- reactive({
    humichip1 <- humichip_probe() %>%
      select(-`Genbank ID`, -gene, -species, -lineage,
           -annotation, -geneCategory, -subcategory1,
           -subcategory2)
  
    # Set NA's to 0 and values not NA to original value
    humichip1 <- humichip1 %>%
      mutate_all(funs(ifelse(is.na(.), 0, .)))
    
    # Remove rows that equal 0
    humichip1 <- humichip1[rowSums(humichip1) != 0,]
    
    as.matrix(humichip1)
  })
 
  
  
  #############################
  #### Ordination Analysis ####
  #############################
  
  # Ordination choices
  humichip_ord <- reactive({
    # Perform PCA analysis using vegan
    if(input$ordination == "PCA"){
      vegan::rda(t(humichip_matrix()))
      
      # Perform DCA analysis using vegan
    } else if (input$ordination == "DCA") {
      vegan::decorana(t(humichip_matrix()))
      
      # Perform PCoA using ape package
    } else if (input$ordination == "PCoA"){
      
      humichip_dist <- vegan::vegdist(as.matrix(t(humichip_matrix())))
      ape::pcoa(humichip_dist, correction = "none")
    } else {
      stopApp()
    }
  })
  
  
  humichip_coords <- reactive({
    
    # PCA and DCA coordinates
    if(input$ordination == "PCA" | input$ordination == "DCA"){
      
      # Get the coordinated for the samples
      ord_coords <- scores(humichip_ord(), display = "sites")
      
      # Make tibble
      ord_coords <- as.data.frame(ord_coords) %>%
        rownames_to_column(var = "glomics_ID")
      
      ord_coords
      
      # PCoA Coords
    } else if (input$ordination == "PCoA") {
      
      ord_coords <- humichip_ord()$vectors[,1:2]
      
      # Make tibble
      ord_coords <- as.data.frame(ord_coords) %>%
        rownames_to_column(var = "glomics_ID")
      
      ord_coords
    } else {
      stopApp()
    }
  })
  
  
  humichip_prop_exp <- reactive({
    
    # Extract proportion explained by first couple PC's for PCA
    if (input$ordination == "PCA"){
      summary(eigenvals(humichip_ord()))[2,] * 100
    } else if (input$ordination == "PCoA"){
      humichip_ord()$values$Relative_eig * 100
    } else {
      stopApp()
    }
  })
  
  
  ####################################
  #### Merge coords into Metadata ####
  ####################################
  
  humichip_coords_metadata <- reactive({
    
    humichip_coords() %>%
      left_join(., B7A_metadata, by = c("glomics_ID")) %>%
      
      # Factor Columns
      mutate("ETEC disease severity score" = factor(`ETEC disease severity score`))
  }) 
  
  
  
  
  ##############################
  #### Create plots/tables #####
  ##############################
  
  plotInput <- reactive({
    
    # Get color input for points
    my_fill <- ifelse(input$color == "None", "NULL", input$color)

    
    # Get size input for points
    my_size <- ifelse(input$point_size == "None", "4", input$point_size)
    
    
    # Aesthetic sizes
    axis_title_size <- 16
    axis_text_size <- 14
    title_size <- 18
    legend_text_size <- 12
    
    
    # PCA Plot
    if(input$ordination == "PCA"){
      pca_plot <- ggplot(humichip_coords_metadata(),
                         aes_string(x = "PC1", 
                                    y = "PC2", 
                                    color = my_fill,
                                    size = my_size)) +
        
        # Set up proportion explained
        xlab(paste0("PC1(", round(humichip_prop_exp()[[1]], 2), "%)")) +
        ylab(paste0("PC2(", round(humichip_prop_exp()[[2]], 2), "%)")) +
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
      dca_plot <- ggplot( humichip_coords_metadata(),
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
      pcoa_plot <- ggplot( humichip_coords_metadata(), 
                          aes_string(x = "Axis.1", 
                                     y = "Axis.2", 
                                     color = my_fill, 
                                     size = my_size)) +
        scale_fill_discrete(name = "Visit Number") +
        xlab(paste0("PCoA1(", round(humichip_prop_exp()[1], 2), "%)")) +
        ylab(paste0("PCoA2(", round(humichip_prop_exp()[2], 2), "%)")) +
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
  
  #output$table <- renderTable({ humichip_coords_metadata()})
  
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