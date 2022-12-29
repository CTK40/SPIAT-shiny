# Made it read csv file and enable column selection.
# TODO: Select column for "marker intensities" (gene expressions)
# TODO: Add a button to save df object
library(shiny)
library(dplyr)
library(SPIAT)
source("functions.R")
ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            # Select file format
            selectInput("format", "Select file format",
                        choices = c("general", "inForm", "HALO", "Xenium", 
                                    "cellprofiler", "CODEX", "IMC", "CosMx")),
            textInput("dir", "Type your folder path"),
            # Select file
            fileInput("file1", "Choose image File/Folder",
                      accept = c(".csv", ".gz", ".txt", ".tsv")),
            # checkboxInput("header", "Header", TRUE),
            # Select file for marker intensity/gene expression
            fileInput("file2", "Choose file for assay",
                      accept = c(".csv", ".gz", ".txt", ".tsv")),
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
                column(3,  numericInput("row2", "Show rows to", value = 50, min = 1)),
                # Add a button to save the df object for marker intensity/gene expression
                column(3, actionButton("do", "Save the marker/gene data frame"))),
            tableOutput("image"),
            tableOutput("markerOrGene")
        )
    )
)

server <- function(input, output, session) {
    # Read in files of different format
    format <- reactive(input$format)
    image <- reactive(
        if (format() == "inForm"){
            vroom::vroom(input$file1$datapath, delim = "\t")
        }else if(format() == "HALO"){
            vroom::vroom(input$file1$datapath)
        }else if (format() == "Xenium"){
            spe <- read_Xenium(samples = input$dir, type = "HDF5", data = "cell")
            data.frame(t(assay((spe))))
        }
    )
    
    markerOrGene <- reactive(
        vroom::vroom(input$file2$datapath)
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
            markerOrGene()[input$row1:input$row2, ] %>% select(!(!!!input$var_gene_ignore))
        }
    }, rownames = TRUE)
    # If click the button
    observeEvent(input$do, {
       # try save object into a var
        new_df <- markerOrGene() %>% select(!!!input$var_gene_select)
    })
    
}
options(shiny.maxRequestSize=10000*1024^2)
shinyApp(ui, server)

