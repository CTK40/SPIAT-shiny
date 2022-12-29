# Made it read csv file and enable column selection.
# TODO: Select column for "marker intensities" (gene expressions)
# TODO: Add a button to save df object
library(shiny)
library(dplyr)
library(SPIAT)
# define functions
.read_xyz <- function(x) {
    cnms <- c(
        "cell_id", "x_centroid", "y_centroid", "transcript_counts", 
        "control_probe_counts", "control_codeword_counts", "total_counts",
        "cell_area", "cell_area")
    df <- read.csv(x, col.names=cnms)
    
    return(df)
}
read_Xenium <- function (samples = "", 
                         sample_id = paste0("sample", sprintf("%02d", seq_along(samples))), 
                         type = c("HDF5", "sparse"), data = c("filtered", "raw", "cell"), images = "lowres", load = TRUE) 
{
    type <- match.arg(type)
    data <- match.arg(data)
    imgs <- c("lowres", "hires", "detected", "aligned")
    imgs <- match.arg(images, imgs, several.ok = TRUE)
    if (is.null(sids <- names(samples))) {
        if (is.null(sids <- sample_id)) {
            stop("'sample_id' mustn't be NULL when 'samples' are unnamed")
        }
        else if (!is.character(sample_id) && length(unique(sample_id)) != 
                 length(samples)) 
            stop("'sample_id' should contain as many unique values as 'samples'")
    }
    else if (length(unique(sids)) != length(samples)) 
        stop("names of 'samples' should be unique")
    names(samples) <- sids
    i <- basename(samples) != "outs"
    samples[i] <- file.path(samples[i], "outs")
    fns <- paste0(data, "_feature_matrix", switch(type, HDF5 = ".h5",  ""))
    counts <- file.path(samples, fns)
    dir <- file.path(samples) # what is this for? the coord? !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    xyz <- file.path(rep(dir, each = length(sids)), "cells.csv.gz")
    xyz <- xyz[file.exists(xyz)]
    names(xyz) <- sids
    # sfs <- file.path(dir, "scalefactors_json.json")
    # names(xyz) <- names(sfs) <- sids
    # img_fns <- list(lowres = "tissue_lowres_image.png", hires = "tissue_hires_image.png", 
    #                 detected = "detected_tissue_image.jpg", aligned = "aligned_fiducials.jpg")
    # img_fns <- img_fns[imgs]
    # img_fns <- lapply(dir, file.path, img_fns)
    # img_fns <- unlist(img_fns)
    # nan <- !file.exists(img_fns)
    # if (all(nan)) {
    #     stop(sprintf("No matching files found for 'images=c(%s)", 
    #                  paste(dQuote(imgs), collapse = ", ")))
    # }
    # else if (any(nan)) {
    #     message("Skipping missing images\n  ", paste(img_fns[nan], 
    #                                                  collapse = "\n  "))
    #     img_fns <- img_fns[!nan]
    # }
    # img <- readImgData(samples, sids, img_fns, sfs, load)
    spel <- lapply(seq_along(counts), function(i) {
        sce <- read10xCounts(samples = counts[i], sample.names = sids[i], 
                             col.names = TRUE)
        spd <- .read_xyz(xyz[i])
        obs <- intersect(colnames(sce), rownames(spd))
        sce <- sce[, obs]
        spd <- spd[obs, ]
        SpatialExperiment(assays = assays(sce), rowData = DataFrame(symbol = rowData(sce)$Symbol), 
                          sample_id = sids[i], colData = DataFrame(spd), spatialCoordsNames = c("x_centroid", 
                                                                                                "y_centroid"))
    })
    spe <- do.call(cbind, spel)
    # imgData(spe) <- img
    return(spe)
}



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

