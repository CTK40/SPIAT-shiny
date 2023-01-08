# Made it read csv file and enable column selection.
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
                        choices = c("general", "inForm", "HALO", "Xenium", "Visium", 
                                     "MERSCOPE", "CosMX", "cellprofiler", "CODEX")),
            uiOutput("format"),
            # checkboxInput("header", "Header", TRUE)
        ),
        mainPanel(
            # Show part of the table
            # Select columns for genes/markers
            fluidRow(
                # This box is for selecting columns
                column(3, uiOutput('cellID_gene')),
                column(4, uiOutput('var_gene_select')),
                column(4,  uiOutput('var_gene_ignore'))),
            # Select the range to show rows in the data (fast)
            fluidRow(
                column(3, numericInput("row1", "Show rows from", value = 1, min = 1)),
                column(3,  numericInput("row2", "Show rows to", value = 10, min = 1)),
                # Add a button to save the df object for marker intensity/gene expression
                column(3, actionButton("do", "Save the marker/gene data frame"))),
            tableOutput("markerOrGene"),
            fluidRow(
                column(3, uiOutput('cellID_metadata')),
                column(3, uiOutput("sample")),
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
        }else if (format() == "Xenium" || format() == "Visium"){
            textInput("dir", "Type your folder path")
        }else if (format() == "MERSCOPE" || format() == "CosMX"){
            fluidRow(
                column(6, fileInput("file_markers2", "Select file for gene expressions:", 
                      accept = c(".csv", ".gz", ".txt", ".tsv"))) , 
                column(6, fileInput("file_meta", "Select file for metadata:", 
                      accept = c(".csv", ".gz", ".txt", ".tsv"))))
        }
    })
    
    # read in the markers/genes
    markerOrGene <- reactive(
        if (format() == "inForm"){
            vroom::vroom(input$file_markers$datapath, delim = "\t")
        }else if(format() == "HALO"){
            vroom::vroom(input$file_markers$datapath)
        }else if (format() == "MERSCOPE" || format() == "CosMX"){
            vroom::vroom(input$file_markers2$datapath)
        }
    )
    
    # read in metadata
    metadata <- reactive(
        if (format() == "MERSCOPE" || format() == "CosMX"){
            vroom::vroom(input$file_meta$datapath)
        }
    )
    
    # Select columns subject to the file selected
    # select column for cell ID
    output$cellID_gene <- renderUI(
        varSelectInput("cellID_gene", label = "Variable for Cell ID:", 
                       data = markerOrGene(), multiple = FALSE))
    # select columns for markers/genes
    output$var_gene_select <- renderUI(
        varSelectInput("var_gene_select", label = "Variable to select:", data = markerOrGene(),
                       multiple = TRUE, width = "500px"))
    # ignore columns for markers/genes
    output$var_gene_ignore <- renderUI(
        varSelectInput("var_gene_ignore", label = "Variable to neglect:", data = markerOrGene(),
                       multiple = TRUE, width = "500px"))
    
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
    output$cellID_metadata <- renderUI(
        varSelectInput("cellID_metadata", label = "Variable for Cell ID:", 
                       data = metadata(), multiple = FALSE))
    output$sample <- renderUI(
        varSelectInput("sample", label = "Varaible of sample ID (fov) to select:", 
                       data = metadata(), multiple = FALSE, width = "500px")
    )
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
        if (format() == "MERSCOPE" || format() == "CosMX"){
            # try save object into a var
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
            Cell_IDs <- data.frame(metadata() %>% select(!!!input$cellID_metadata))[, 1]
            Sample_IDs <- as.character(data.frame(metadata() %>% select(!!!input$sample))[, 1])
            phenotypes <- data.frame(metadata() %>% select(!!!input$phenotype))[, 1]
            coord_x <- data.frame(metadata() %>% select(!!!input$coord_x))[, 1]
            coord_y <- data.frame(metadata() %>% select(!!!input$coord_y))[, 1]
            
            metadata_df <- data.frame(Cell_IDs, Sample_IDs, phenotypes, coord_x, coord_y)
            
            # match the cell IDs with the cell IDs from gene expression matrix
            metadata_df_update <- metadata_df[match(metadata_df$Cell_IDs, 
                                                    data.frame(markerOrGene() %>% select(!!!input$cellID_gene))[, 1]), ]
            
            general_format_image <- format_image_to_spe(format = "general", 
                                                        intensity_matrix = new_df,
                                                        Cell_IDs = metadata_df_update$Cell_IDs,
                                                        Sample_IDs = metadata_df_update$Sample_IDs, 
                                                        phenotypes = metadata_df_update$phenotypes,
                                                        coord_x = metadata_df_update$coord_x, 
                                                        coord_y = metadata_df_update$coord_y)
            save(general_format_image, file = "spe.Rda")
        }else if (format() == "Xenium"){
            Xenium_spe <- read_Xenium(samples = input$dir, type = "HDF5", data = "cell")
            save(Xenium_spe, file = "Xenium_spe.Rda")
        }else if (format() == "Visium"){
            
            Visium_spe <- SpatialExperiment::read10xVisium(
                samples = input$dir, type = "HDF5", data = "raw")
            save(Visium_spe, file = "Visium_spe.Rda")
        }
    })
    
}
options(shiny.maxRequestSize=200000*1024^2)
shinyApp(ui, server)
