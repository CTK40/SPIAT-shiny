# Made it read csv file and enable column selection.
# TODO: add selecting Cell ID 
library(shiny)
library(dplyr)
library(SPIAT)
library(purrr)
source("functions.R")
ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            # Select file format
            selectInput("format", "Select file format",
                        choices = c("general", "inForm", "HALO", "Xenium", 
                                    "cellprofiler", "CODEX", "IMC", "CosMx")),
            uiOutput("format"),
            # checkboxInput("header", "Header", TRUE)
            
        ),
        mainPanel(
            # Show part of the table
            fluidRow(
                # This box is for selecting columns
                column(6, uiOutput('variables1')),
                column(6,  uiOutput('variables2'))),
            # Select columns for genes/markers
            fluidRow(
                # This box is for selecting columns
                column(6, uiOutput('var_gene_select')),
                column(6,  uiOutput('var_gene_ignore'))),
            # Select the range to show rows in the data (fast)
            fluidRow(
                column(3, numericInput("row1", "Show rows from", value = 1, min = 1)),
                column(3,  numericInput("row2", "Show rows to", value = 10, min = 1)),
                # Add a button to save the df object for marker intensity/gene expression
                column(3, actionButton("do", "Save the marker/gene data frame"))),
            tableOutput("image"),
            tableOutput("markerOrGene"),
            fluidRow(
                column(3, uiOutput("phenotype")),
                column(3, uiOutput("coord_x")),
                column(3, uiOutput("coord_y"))
            ),
            tableOutput("metadata")
        )
    )
)

server <- function(input, output, session) {
    # Read in files of different format
    format <- reactive(input$format)
    # The data reading field changes with different formats
    output$format <- renderUI({
        if (format() == "inForm" || format() == "HALO"){
            fileInput("file_markers", "Choose image file",
                      accept = c(".csv", ".gz", ".txt", ".tsv"))
        }else if (format() == "Xenium"){
            textInput("dir", "Type your folder path")
        }else if (format() == "IMC"){
            map(c("file_markers2", "file_meta"),  # how to add different labels here???
                ~ fileInput(.x, NULL, 
                      accept = c(".csv", ".gz", ".txt", ".tsv")))
        }
    })
    
    # read in the markers/genes
    image <- reactive(
        if (format() == "inForm"){
            vroom::vroom(input$file_markers$datapath, delim = "\t")
        }else if(format() == "HALO" || format() == "IMC"){
            vroom::vroom(input$file_markers$datapath)
        }else if (format() == "Xenium"){
            spe <- read_Xenium(samples = input$dir, type = "HDF5", data = "cell")
            data.frame(t(assay((spe))))
        }
    )
    
    markerOrGene <- reactive(
        vroom::vroom(input$file_markers2$datapath)
    )
    
    # read in metadata
    metadata <- reactive(
        if (format() == "IMC"){
            vroom::vroom(input$file_meta$datapath)
        }
    )
    
    # Select columns subject to the file selected
    # select columns 
    output$variables1 <- renderUI(
        varSelectInput("variables1", label = "Variable to select:", data = image(),
                       multiple = TRUE, width = "500px"))
    # ignore columns
    output$variables2 <- renderUI(
        varSelectInput("variables2", label = "Variable to neglect:", data = image(),
                       multiple = TRUE, width = "500px"))
    
    # select columns for markers/genes
    output$var_gene_select <- renderUI(
        varSelectInput("var_gene_select", label = "Variable to select:", data = markerOrGene(),
                       multiple = TRUE, width = "500px"))
    # ignore columns for markers/genes
    output$var_gene_ignore <- renderUI(
        varSelectInput("var_gene_ignore", label = "Variable to neglect:", data = markerOrGene(),
                       multiple = TRUE, width = "500px"))
    
    # Visualise the table
    output$image <- renderTable({
        if (length(input$variables1) == 0 && length(input$variables2) == 0) {
            return(image()[input$row1:input$row2, ])
        }else if (length(input$variables1) != 0 && length(input$variables2) == 0){
            image()[input$row1:input$row2, ] %>% select(!!!input$variables1)
        }else if (length(input$variables1) == 0 && length(input$variables2) != 0){
            image()[input$row1:input$row2, ] %>% select(!(!!!input$variables2))
        }
    }, rownames = TRUE)
    
    # Visualise the table
    output$markerOrGene <- renderTable({
        if (length(input$var_gene_select) == 0 && length(input$var_gene_ignore) == 0) {
            return(markerOrGene()[input$row1:input$row2, ])
        }else if (length(input$var_gene_select) != 0 && length(input$var_gene_ignore) == 0){
            markerOrGene()[input$row1:input$row2, ] %>% select(!!!input$var_gene_select)
        }else if (length(input$var_gene_select) == 0 && length(input$var_gene_ignore) != 0){
            # Couldn't find a way to directly 
            temp <- markerOrGene()[input$row1:input$row2, ] 
            for (i in 1:length(input$var_gene_ignore)){
                temp <- temp %>% select(!(!!input$var_gene_ignore[[i]]))
            }
            temp
        }
    }, rownames = TRUE)
    
    # visualise the metadata
    output$phenotype <- renderUI(
        varSelectInput("phenotype", label = "Variable of phenotype to select:", data = metadata(),
                   multiple = FALSE, width = "500px")
    )
    output$coord_x <- renderUI(
        varSelectInput("coord_x", label = "Variable of x coordinates to select:", 
                     data = metadata(), multiple = FALSE, width = "500px")
    )
    output$coord_y <- renderUI(
    varSelectInput("coord_y", label = "Variable of y coordinates to select:", 
                   data = metadata(), multiple = FALSE, width = "500px")
    )
    output$metadata <- renderTable({
        return(metadata()[1:10, ])
    })
    # If click the button
    observeEvent(input$do, {
       # try save object into a var
        # Note this intensity matrix code is wrong (need to save all data if there is no selected vars)
        if (length(input$var_gene_select == 0) && length(input$var_gene_ignore == 0)) {
            new_df <- markerOrGene()
        }else if (length(input$var_gene_select != 0)){
            new_df <- data.frame(markerOrGene() %>% select(!!!input$var_gene_select))
        }else{
            temp <- markerOrGene()
            for (i in 1:length(input$var_gene_ignore)){
                temp <- temp %>% select(!(!!input$var_gene_ignore[[i]]))
            }
            new_df <- temp
        }
        
        # Note for these metadata, they are still data frames but not vectors, so need to add [, 1] at the end
        phenotypes <- data.frame(metadata() %>% select(!!!input$phenotype))[, 1]
        coord_x <- data.frame(metadata() %>% select(!!!input$coord_x))[, 1]
        coord_y <- data.frame(metadata() %>% select(!!!input$coord_y))[, 1]
        Cell_IDs <- 1:dim(markerOrGene())[1]
        general_format_image <- format_image_to_spe(format = "general", 
                                                    intensity_matrix = new_df,
                                                    Cell_IDs = Cell_IDs,
                                                    phenotypes = phenotypes,
                                                    coord_x = coord_x, coord_y = coord_y)
        save(general_format_image, file = "spe.Rda")
    })
    
}
options(shiny.maxRequestSize=10000*1024^2)
shinyApp(ui, server)

