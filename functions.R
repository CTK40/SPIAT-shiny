library(DropletUtils)
# TODO: To add selecting sample column

renderUI_varSelect <- function(inputId, label, data, selected = NULL, multiple = FALSE, width = NULL){
    renderUI(varSelectInput(inputId, label, data, selected, multiple, width))
}

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

format_image <- function(format, var_gene_select = input$var_gene_select, var_gene_ignore = input$var_gene_ignore,
                         df = markerOrGene(), df_meta = metadata(),
                         cellID_metadata = input$cellID_metadata, 
                         sample = input$sample, phenotype = input$phenotype,
                         coord_x = input$coord_x, coord_y = input$coord_y,
                         cellID_gene = input$cellID_gene,
                         fov_gene = input$fov_gene,
                         dir = input$dir, markers = NULL,
                         path1 = input$file_markers$datapath){
    if (format == "MERSCOPE" ){
        # try save object into a var
        if (length(var_gene_select) == 0 && length(var_gene_ignore) == 0) {
            new_df <- df
        }else if (length(var_gene_select) != 0){
            new_df <- data.frame(df %>% select(!!!var_gene_select))
        }else{
            temp <- df
            for (i in 1:length(var_gene_ignore)){
                temp <- temp %>% select(!(!!var_gene_ignore[[i]]))
            }
            new_df <- temp
        }
        # Note for these metadata, they are still data frames but not vectors, so need to add [, 1] at the end
        Cell_IDs <- data.frame(df_meta %>% select(!!!cellID_metadata))[, 1]
        Sample_IDs <- as.character(data.frame(df_meta %>% select(!!!sample))[, 1])
        phenotypes <- data.frame(df_meta %>% select(!!!phenotype))[, 1]
        coord_x <- data.frame(df_meta %>% select(!!!coord_x))[, 1]
        coord_y <- data.frame(df_meta %>% select(!!!coord_y))[, 1]
        
        metadata_df <- data.frame(Cell_IDs, Sample_IDs, phenotypes, coord_x, coord_y)
        
        # match the cell IDs with the cell IDs from gene expression matrix
        metadata_df_update <- metadata_df[match(metadata_df$Cell_IDs, 
                                                data.frame(df %>% select(!!!cellID_gene))[, 1]), ]
        
        general_format_image <- format_image_to_spe(format = "general", 
                                                    intensity_matrix = new_df,
                                                    Cell_IDs = metadata_df_update$Cell_IDs,
                                                    Sample_IDs = metadata_df_update$Sample_IDs, 
                                                    phenotypes = metadata_df_update$phenotypes,
                                                    coord_x = metadata_df_update$coord_x, 
                                                    coord_y = metadata_df_update$coord_y)
        save(general_format_image, file = "Objects/spe.Rda")
    }else if (format == "CosMX") {
        new_df <- df
        newID_df <- data.frame(
            Cell_IDs = as.numeric(data.frame(new_df %>% select(!!!cellID_gene))[, 1]),
            Sample_IDs = as.numeric(data.frame(new_df %>% select(!!!fov_gene))[, 1] ))
        newID_df$Cell_IDs <- as.character(newID_df$Cell_IDs)
        newID_df$Sample_IDs <- as.character(newID_df$Sample_IDs)
        new_df$new_cellID <- apply(newID_df[, c("Cell_IDs", "Sample_IDs")], 1, digest)
        
        # Note for these metadata, they are still data frames but not vectors, so need to add [, 1] at the end
        Cell_IDs <- data.frame(df_meta %>% select(!!!cellID_metadata))[, 1]
        Sample_IDs <- as.character(data.frame(df_meta %>% select(!!!sample))[, 1])
        phenotypes <- data.frame(df_meta %>% select(!!!phenotype))[, 1]
        coord_x <- data.frame(df_meta %>% select(!!!coord_x))[, 1]
        coord_y <- data.frame(df_meta %>% select(!!!coord_y))[, 1]
        
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
        if (length(var_gene_select == 0) && length(var_gene_ignore == 0)) {
            meaningless <- 0
        }else if (length(var_gene_select != 0)){
            new_df_update <- data.frame(new_df_update %>% select(!!!var_gene_select))
        }else{
            temp <- new_df_update
            for (i in 1:length(var_gene_ignore)){
                temp <- temp %>% select(!(!!var_gene_ignore[[i]]))
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
        save(general_format_image, file = "Objects/spe.Rda")
    }else if (format == "Xenium"){
        Xenium_spe <- read_Xenium(samples = dir, type = "HDF5", data = "cell")
        save(Xenium_spe, file = "Objects/Xenium_spe.Rda")
    }else if (format == "Visium"){
        Visium_spe <- SpatialExperiment::read10xVisium(
            samples = dir, type = "HDF5", data = "raw")
        spatialCoordsNames(Visium_spe) <- c("Cell.X.Position", "Cell.Y.Position")
        save(Visium_spe, file = "Objects/Visium_spe.Rda")
    }else if (format == "inForm"){
        cols <- substr(as.character(unlist(var_gene_select)),2,stringr::str_length(as.character(unlist(var_gene_select)))-1)
        print(cols)
        print(markers)
        inForm_spe <- format_inform_to_spe(path = path1, markers = markers, 
                                           intensity_columns_interest = cols)
        save(inForm_spe, file = "Objects/inForm_spe.Rda")
    }
}


get_colData <- function(spe_object){
    formatted_data <- data.frame(SummarizedExperiment::colData(spe_object))
    formatted_data <- cbind(formatted_data, 
                            data.frame(SpatialExperiment::spatialCoords(spe_object)))
    if (is.null(formatted_data$Cell.ID)){
        formatted_data <- formatted_data %>% tibble::rownames_to_column("Cell.ID")
    }
    
    # shouldn't delete column `sample_id`
    # formatted_data$sample_id <- NULL
    
    return(formatted_data)
}


# read marker or gene
read_markerGene <- function(format, path1, path2){
    if (format == "inForm"){
        if (is.null(path1)){
            path1 <- "Data/example_inForm/tiny_inform.txt.gz"
        }
        df <- vroom::vroom(path1, delim = "\t")
    }else if(format == "HALO"){
        df <- vroom::vroom(path)
    }else if (format == "MERSCOPE" || format == "CosMX"){
        df <- vroom::vroom(path2)
    }
    return(df)
}
