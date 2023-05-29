#TODO: split this file into separate files
library(SPIAT)
library(dplyr)
library(purrr)
library(digest)
library(RColorBrewer)
source("functions.R")
# increase the max size of file upload
options(shiny.maxRequestSize=30000*1024^2)
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
    output$fov_gene <- renderUI(
        if (format() == "CosMX"){
            varSelectInput("fov_gene", label = "Variable for fov:", 
                           data = markerOrGene(), multiple = FALSE)
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
        if (format() == "MERSCOPE" ){
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
        }else if (format() == "CosMX") {
            new_df <- markerOrGene()
            newID_df <- data.frame(
                Cell_IDs = as.numeric(data.frame(new_df %>% select(!!!input$cellID_gene))[, 1]),
                Sample_IDs = as.numeric(data.frame(new_df %>% select(!!!input$fov_gene))[, 1] ))
            newID_df$Cell_IDs <- as.character(newID_df$Cell_IDs)
            newID_df$Sample_IDs <- as.character(newID_df$Sample_IDs)
            new_df$new_cellID <- apply(newID_df[, c("Cell_IDs", "Sample_IDs")], 1, digest)
            
            # Note for these metadata, they are still data frames but not vectors, so need to add [, 1] at the end
            Cell_IDs <- data.frame(metadata() %>% select(!!!input$cellID_metadata))[, 1]
            Sample_IDs <- as.character(data.frame(metadata() %>% select(!!!input$sample))[, 1])
            phenotypes <- data.frame(metadata() %>% select(!!!input$phenotype))[, 1]
            coord_x <- data.frame(metadata() %>% select(!!!input$coord_x))[, 1]
            coord_y <- data.frame(metadata() %>% select(!!!input$coord_y))[, 1]
            
            metadata_df <- data.frame(Cell_IDs, Sample_IDs, phenotypes, coord_x, coord_y)
            cellID_df <- metadata_df[, c("Cell_IDs", "Sample_IDs")]
            cellID_df$Cell_IDs <- as.character(cellID_df$Cell_IDs)
            cellID_df$Sample_IDs <- as.character(cellID_df$Sample_IDs)
            metadata_df$new_cellID <- apply(cellID_df, 1, digest)
            
            # match the cell IDs with the cell IDs from gene expression matrix
            metadata_df_update <- metadata_df[match(metadata_df$new_cellID, 
                                                    new_df$new_cellID), ]
            metadata_df_update <- metadata_df_update[complete.cases(metadata_df_update), ]
            
            new_df_update <- new_df[match(new_df$new_cellID, metadata_df_update$new_cellID),]
            new_df_update <- new_df_update[complete.cases(new_df_update), ]
            
            # clean up the gene expression matrix
            if (length(input$var_gene_select == 0) && length(input$var_gene_ignore == 0)) {
                meaningless <- 0
            }else if (length(input$var_gene_select != 0)){
                new_df_update <- data.frame(new_df_update %>% select(!!!input$var_gene_select))
            }else{
                temp <- new_df_update
                for (i in 1:length(input$var_gene_ignore)){
                    temp <- temp %>% select(!(!!input$var_gene_ignore[[i]]))
                }
                new_df_update <- temp
            }
            
            
            general_format_image <- format_image_to_spe(format = "general", 
                                                        intensity_matrix = new_df_update,
                                                        Cell_IDs = metadata_df_update$new_cellID,
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
    spe <- reactive({
        shinyFileChoose(input, "Btn_GetFile", roots = c(wd='.'), session = session)

        if(!is.null(input$Btn_GetFile)){
            # browser()
            file_selected<-parseFilePaths(c(wd='.'), input$Btn_GetFile)
            file_path <- as.character(file_selected$datapath)
            if (length(file_path) != 0){
                load(file_path)
                temp_var <- load(file_path)
                eval(parse(text = temp_var))
            }
        }
    })
    
    coldata <- reactive({
        get_colData(spe())
    })
    output$feature_colname <- renderUI(
        varSelectInput("feature_colname", label = "Column name", 
                       data = coldata(), multiple = FALSE))
    output$categories <- renderUI(
        selectInput("categories", label = "Cells of interest", 
                     choices = unique(coldata()[[input$feature_colname]]),
                    multiple = TRUE)
        )
    qual_col_pals <- reactive(brewer.pal.info[brewer.pal.info$category == 'qual',])
    col_vector <- reactive(unlist(mapply(brewer.pal, qual_col_pals()$maxcolors, rownames(qual_col_pals()))))
    
    output$colour <- renderUI(
        selectInput('colour', label = "Colour vector", choices = col_vector(),
                    selected = col_vector()[1:length(input$categories)],
                    multiple = TRUE)
        )

    output$plot_cells <- renderPlot(
        {
            if (length(spe()) != 0 ){
                g <- plot_cell_categories(spe(), categories_of_interest = input$categories, 
                                     feature_colname = input$feature_colname,
                                     colour_vector = input$colour)
                plot(g, width = 10, height = 10)
            }
        }
    )
    markers <- reactive(rownames(assay(spe())))
    output$marker_level <- renderUI(
        selectInput('marker_level', label = "Select the marker:", 
                    choices = markers())
    )
    output$plot_marker <- renderPlot(
        {
            if (length(spe()) != 0 ){
                p <- plot_cell_marker_levels(spe(), marker = input$marker_level)
                plot(p, width = 10, height = 10)
            }
        }
    )
}
