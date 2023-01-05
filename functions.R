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
            View(intensity_matrix)
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
                assay = assay_data_matrix_t,
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
    } else if(format == "CODEX"){
        spe <- format_codex_to_spe(path = path, markers = markers,
                                   path_to_codex_cell_phenotypes = path_to_codex_cell_phenotypes)
    }else if(format == "cellprofiler") {
        spe <- format_cellprofiler_to_spe(path = path, markers = markers,
                                          intensity_columns_interest = intensity_columns_interest)
    } else { methods::show("Please entre a valid format!" )}
    return(spe)
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