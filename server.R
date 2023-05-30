# TODO: be able to format HALO images
library(SPIAT)
library(dplyr)
library(purrr)
library(digest)
library(RColorBrewer)
# Save some colours in the namespace
qual_col_pals <- reactive(brewer.pal.info[brewer.pal.info$category == 'qual',])
col_vector <- reactive(unlist(mapply(brewer.pal, qual_col_pals()$maxcolors, rownames(qual_col_pals()))))

source("functions.R")
# increase the max size of file upload
options(shiny.maxRequestSize=30000*1024^2)

server <- function(input, output, session) {
    # Tab 1 ####
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
    markerOrGene <- reactive({
        read_markerGene(format(), path1 = input$file_markers$datapath,
                        path2 = input$file_markers2$datapath)})
    
    # read in metadata
    metadata <- reactive(
        if (format() == "MERSCOPE" || format() == "CosMX"){
            vroom::vroom(input$file_meta$datapath)
        }
    )
    reactive({
        if (format() != "Visium" && format() != "Xenium"){
            # Select columns subject to the file selected
            # select column for cell ID
            output$cellID_gene <- renderUI(
                varSelectInput("cellID_gene", label = "Variable for Cell ID:", 
                               data = markerOrGene(), multiple = FALSE))
            output$fov_gene <- renderUI(
                if (format() == "CosMX"){
                    varSelectInput("fov_gene", label = "Variable for fov:", 
                                   data = markerOrGene(), multiple = FALSE)
                })
            # The number of textboxes to input markers depends on the number of markers
            marker_names <- reactive(paste0("marker", seq_len(input$n_markers)))
            output$marker <- renderUI({
                map(marker_names(), ~ column(3, textInput(.x, NULL, value = "AMACR")))
            })
            markers <- reactive({
                temp <- c()
                for (i in seq_len(input$n_markers)){
                    temp <- c(temp, eval(parse(text = paste0("input$marker", i))))
                }
                temp
            })
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
            reactive({
                if (!is.null(metadata())){
                    output$cellID_metadata <- renderUI_varSelect("cellID_metadata", "Variable for Cell ID:",
                                                                 data = metadata())
                    output$sample <- renderUI_varSelect("sample", "Varaible of sample ID (fov) to select:", 
                                                        data = metadata(), width = "500px")
                    output$phenotype <- renderUI_varSelect("phenotype", "Variable of phenotype to select:", 
                                                           data = metadata(),width = "500px")
                    output$coord_x <- renderUI_varSelect("coord_x", "Variable of x coordinates to select:", 
                                                         data = metadata(), width = "500px")
                    output$coord_y <- renderUI_varSelect("coord_y", "Variable of y coordinates to select:", 
                                                         data = metadata(),  width = "500px")
                    output$metadata <- renderTable({
                        return(metadata()[1:10, ])
                    })
                }
            })
        }
    })
    
    
   
    # If click the button, save the object
    observeEvent(input$do, {
        format_image(format = format(), var_gene_select = input$var_gene_select, var_gene_ignore = input$var_gene_ignore,
                     df = markerOrGene(), df_meta = metadata(),
                     cellID_metadata = input$cellID_metadata, 
                     sample = input$sample, phenotype = input$phenotype,
                     coord_x = input$coord_x, coord_y = input$coord_y,
                     cellID_gene = input$cellID_gene,
                     fov_gene = input$fov_gene,
                     dir = input$dir, markers = markers(),
                     path1 = input$file_markers$datapath)
    })
    
    # load object and Tab 2 ####
    # load spe object
    observe({
        shinyFileChoose(input, "Btn_GetFile", roots = c(wd='.'), session = session)

    })
    spe <- reactive({
        if(!is.null(input$Btn_GetFile)){
            # browser()
            file_selected<-parseFilePaths(c(wd='.'), input$Btn_GetFile)
            file_path <- as.character(file_selected$datapath)
            if (length(file_path) != 0){
                load(file_path)
                temp_var <- load(file_path)
                eval(parse(text = temp_var))}}
    })
    
    # format spe object
    coldata <- reactive({ get_colData(spe())})
    
    # get the arguments from the ui
    output$feature_colname <- renderUI(
        varSelectInput("feature_colname", "Column name", data = coldata()))
    output$categories <- renderUI(
        selectInput("categories", label = "Cells of interest", 
                     choices = unique(coldata()[[input$feature_colname]]),
                    multiple = TRUE))
    
    output$colour <- renderUI(
        selectInput('colour', label = "Colour vector", choices = col_vector(),
                    selected = col_vector()[1:length(input$categories)],
                    multiple = TRUE))

    output$plot_cells <- renderPlot({
        if (length(spe()) != 0 ){
            g <- plot_cell_categories(spe(), categories_of_interest = input$categories,
                                      feature_colname = input$feature_colname,
                                      colour_vector = input$colour)
            plot(g, width = 10, height = 10)}
        })
    
    # Tab 3 ####
    markers2 <- reactive(rownames(assay(spe())))
    output$marker_level <- renderUI(
        selectInput('marker_level', label = "Select the marker:", 
                    choices = markers2()))
    output$plot_marker <- renderPlot({
        if (length(spe()) != 0 ){
            p <- plot_cell_marker_levels(spe(), marker = input$marker_level)
            plot(p, width = 10, height = 10)}
        })
    #####
}
