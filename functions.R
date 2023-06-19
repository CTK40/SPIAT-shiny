library(DropletUtils)
# TODO: To add selecting sample column
format_image_to_spe <- function(format = "general", intensity_matrix = NULL,
                                Sample_IDs = NULL, Cell_IDs = NULL, 
                                phenotypes = NULL, coord_x = NULL,
                                coord_y = NULL, path = NULL, markers = NULL,
                                locations = NULL,
                                intensity_columns_interest = NULL,
                                dye_columns_interest = NULL,
                                path_to_codex_cell_phenotypes = NULL){
    if (format == "general"){
        if (is.null(coord_x) || is.null(coord_y)) stop("Cell locations are missing!")
        if (is.null(phenotypes)){
            phenotypes <- rep("Undefined", length(coord_x))  
        }
        if (is.null(Cell_IDs)){
            if ("matrix" %in% class(intensity_matrix)){
                Cell_IDs <- colnames(intensity_matrix)
            }else Cell_IDs <- seq_len(length(coord_x))
        }
        if (is.null(Sample_IDs)){
            Sample_IDs <- "sample01"
        }
        if (is.null(intensity_matrix)){
            metadata_columns <- data.frame(Cell.ID = Cell_IDs,
                                           Phenotype = phenotypes,
                                           Cell.X.Position = coord_x,
                                           Cell.Y.Position = coord_y)
            spe <- SpatialExperiment::SpatialExperiment(
                assay = NULL,
                colData = metadata_columns,
                spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position"))
        }
        else{
            if ("matrix" %in% class(intensity_matrix)){
                markers <- rownames(intensity_matrix)
                intensity_columns <- t(intensity_matrix)
                Cell_IDs_update <- Cell_IDs
            }
            if ("data.frame" %in% class(intensity_matrix)){
                markers <- colnames(intensity_matrix)
                Cell_IDs_update <- Cell_IDs
                intensity_columns <- intensity_matrix
            }
            metadata_columns <- data.frame(Cell.ID = Cell_IDs,
                                           Phenotype = phenotypes,
                                           Cell.X.Position = coord_x,
                                           Cell.Y.Position = coord_y)
            #transpose the matrix so every column is a cell and every row is a marker
            assay_data_matrix <- as.matrix(intensity_columns)
            colnames(assay_data_matrix) <- NULL
            rownames(assay_data_matrix) <- NULL
            assay_data_matrix_t <- t(assay_data_matrix)
            spe <- SpatialExperiment::SpatialExperiment(
                assay = list(counts = assay_data_matrix_t),
                colData = metadata_columns,
                sample_id = Sample_IDs, 
                spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position"))
            rownames(spe) <- markers
            colnames(spe) <- Cell_IDs_update
        }
    } else if (format == "inForm") {
        spe <- format_inform_to_spe(path = path, markers = markers,
                                    locations = locations,
                                    intensity_columns_interest =
                                    intensity_columns_interest)
    }else if (format == "HALO"){
        spe <- format_halo_to_spe(path = path, markers = markers,
                                  locations = locations,
                                  dye_columns_interest = dye_columns_interest,
                                  intensity_columns_interest =
                                  intensity_columns_interest)
    } else { methods::show("Please entre a valid format!" )}
  
    return(spe)

}


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

format_image <- function(format, var_gene_select = input$var_gene_select, 
                         var_gene_ignore = input$var_gene_ignore,
                         intensity_columns = input$intensity_columns,
                         dye_columns = input$dye_columns,
                         df = markerOrGene(), df_meta = metadata(),
                         cellID_metadata = input$cellID_metadata, 
                         sample = input$sample, phenotype = input$phenotype,
                         coord_x = input$coord_x, coord_y = input$coord_y,
                         cellID_gene = input$cellID_gene,
                         fov_gene = input$fov_gene,
                         dir = input$dir, markers = NULL,
                         path1 = input$file_markers$datapath){
    if (format == "general" ){
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
      save(general_format_image, file = "Objects/general_spe.Rda")
    } else if (format == "MERSCOPE" ){
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
        
        MERSCOPE_spe <- format_image_to_spe(format = "general", 
                                            intensity_matrix = new_df,
                                            Cell_IDs = metadata_df_update$Cell_IDs,
                                            Sample_IDs = metadata_df_update$Sample_IDs, 
                                            phenotypes = metadata_df_update$phenotypes,
                                            coord_x = metadata_df_update$coord_x, 
                                            coord_y = metadata_df_update$coord_y)
        save(MERSCOPE_spe, file = "Objects/MERSCOPE_spe.Rda")
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
        
        CosMX_spe <- format_image_to_spe(format = "general", 
                                        intensity_matrix = new_df_update,
                                        Cell_IDs = metadata_df_update$new_cellID,
                                        Sample_IDs = metadata_df_update$Sample_IDs, 
                                        phenotypes = metadata_df_update$phenotypes,
                                        coord_x = metadata_df_update$coord_x, 
                                        coord_y = metadata_df_update$coord_y)
        save(CosMX_spe, file = "Objects/CosMX_spe.Rda")
    }else if (format == "Xenium"){
        Xenium_spe <- read_Xenium(samples = dir, type = "HDF5", data = "cell")
        save(Xenium_spe, file = "Objects/Xenium_spe.Rda")
    }else if (format == "Visium"){
        Visium_spe <- SpatialExperiment::read10xVisium(
            samples = dir, type = "HDF5", data = "filtered")
        spatialCoordsNames(Visium_spe) <- c("Cell.X.Position", "Cell.Y.Position")
        save(Visium_spe, file = "Objects/Visium_spe.Rda")
    }else if (format == "inForm"){
        cols <- substr(as.character(unlist(intensity_columns)),2,
                       stringr::str_length(as.character(unlist(intensity_columns)))-1)
        inForm_spe <- format_inform_to_spe(path = path1, markers = markers, 
                                           intensity_columns_interest = cols)
        save(inForm_spe, file = "Objects/inForm_spe.Rda")
    }else if (format == "HALO") {
      cols1 <- substr(as.character(unlist(intensity_columns)),2,
                     stringr::str_length(as.character(unlist(intensity_columns)))-1)
      cols2 <- substr(as.character(unlist(dye_columns)),2,
                      stringr::str_length(as.character(unlist(dye_columns)))-1)
      HALO_spe <- format_halo_to_spe(path = path1, markers = markers, 
                                         intensity_columns_interest = cols1,
                                         dye_columns_interest = cols2)
      save(HALO_spe, file = "Objects/HALO_spe.Rda")
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
            path1 <- "Data/example/tiny_inform.txt.gz"
        }
        df <- vroom::vroom(path1, delim = "\t")
    }else if(format == "HALO"){
      if (is.null(path1)){
        path1 <- "Data/example/tiny_halo.csv.gz"
      }
        df <- vroom::vroom(path1)
    }else if (format == "MERSCOPE" || format == "CosMX" || format == "general"){
        df <- vroom::vroom(path2)
    }
    return(df)
}




plot_cell_marker_levels <- function(spe_object, marker) {
  
  Cell.X.Position <- Cell.Y.Position <- NULL
  
  intensity_matrix <- SummarizedExperiment::assay(spe_object)
  markers <- rownames(intensity_matrix)
  
  #CHECK
  if (is.element(marker, markers) == FALSE) {
    stop("The marker specified is not in the data")
  }
  
  formatted_data <- bind_info(spe_object)
  
  
  #selecting the cells that have intensity for a specific marker
  column <- which(colnames(formatted_data) == marker)
  rows_non_zero <- which(formatted_data[,column] != 0)
  intensity_by_marker <- formatted_data[rows_non_zero,]
  
  if (nrow(intensity_by_marker) == 0) {
    methods::show(paste("There are no true intensity for: ", marker, sep=""))
  }
  
  #log the intensity to improve contrast
  intensity_by_marker[,marker] <- log10(intensity_by_marker[,marker])
  
  ggplot(intensity_by_marker, aes(x = Cell.X.Position, y = Cell.Y.Position,
                                  colour = eval(parse(text = marker)))) +
    geom_point(aes(colour=eval(parse(text = marker))),size = 0.1) +
    ggtitle(marker) +
    guides(alpha = "none") + scale_colour_viridis_c(direction = -1) +
    labs(colour = paste("log10","(", as.character(marker),
                        " Intensity", ")",sep="")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}



# convert spe object to a data frame with both colData and intensity matrix
bind_info <- function(spe_object){
  formatted_data <- data.frame(SummarizedExperiment::colData(spe_object))
  formatted_data <- cbind(formatted_data, 
                          data.frame(SpatialExperiment::spatialCoords(spe_object)))
  #convert rowname to column
  if (is.null(formatted_data$Cell.ID)){
    formatted_data <- formatted_data %>% tibble::rownames_to_column("Cell.ID") 
  }
  # get the intensity matrix
  intensity_matrix <- SummarizedExperiment::assay(spe_object)
  markers <- rownames(intensity_matrix)
  cell_ids <- colnames(intensity_matrix)
  rownames(intensity_matrix) <- NULL
  colnames(intensity_matrix) <- NULL
  intensity_t <- data.frame(t(intensity_matrix))
  colnames(intensity_t) <- markers
  
  # bind
  formatted_data <- cbind(formatted_data, intensity_t)
  
  # shouldn't delete column `sample_id`
  # formatted_data$sample_id <- NULL
  
  return(formatted_data)
}
library(ggplot2)
# convert spe object to a data frame with only colData
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


