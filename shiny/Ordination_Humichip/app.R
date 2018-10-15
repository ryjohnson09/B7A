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
                 helpText("B: Sample at acute symtoms or 5 days post innoculation"),
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
      
      # Subset the humichip data
      humichip %>%
        select_if(colnames(.) %in% c("Genbank.ID", "gene", "species", "lineage",
                                     "annotation", "geneCategory", "subcategory1",
                                     "subcategory2", matched_samples))
      
    } else {
      humichip
    }
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

}







## Run App --------------------------------------------------------------------------------
shinyApp(server = server, ui = ui)