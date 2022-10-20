#' Add Mito and Ribo percentages to LIGER
#'
#' Add Mito, Ribo, percentages to meta.data slot of LIGER Object
#'
#' @param liger_object LIGER object name.
#' @param species Species of origin for given Seurat Object.  If mouse, human, marmoset, zebrafish, rat,
#' drosophila, or rhesus macaque (name or abbreviation) are provided the function will automatically
#' generate mito_pattern and ribo_pattern values.
#' @param mito_name name to use for the new meta.data column containing percent mitochondrial counts.
#' Default is "percent_mito".
#' @param ribo_name name to use for the new meta.data column containing percent ribosomal counts.
#' Default is "percent_mito".
#' @param mito_ribo_name name to use for the new meta.data column containing percent mitochondrial+ribosomal
#' counts.  Default is "percent_mito".
#' @param mito_pattern A regex pattern to match features against for mitochondrial genes (will set automatically
#' if species is mouse or human; marmoset features list saved separately).
#' @param ribo_pattern A regex pattern to match features against for ribosomal genes (will set automatically
#' if species is mouse, human, or marmoset).
#' @param mito_features A list of mitochrondial gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param ribo_features A list of ribosomal gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param ensembl_ids logical, whether feature names in the object are gene names or
#' ensembl IDs (default is FALSE; set TRUE if feature names are ensembl IDs).
#' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
#' function will abort if columns with any one of the names provided to `mito_name` `ribo_name` or `mito_ribo_name`
#' is present in meta.data slot.
#' @param list_species_names returns list of all accepted values to use for default species names which
#' contain internal regex/feature lists (human, mouse, marmoset, zebrafish, rat, drosophila, and
#' rhesus macaque).  Default is FALSE.
#'
#' @import cli
#' @importFrom dplyr mutate select intersect
#' @importFrom magrittr "%>%"
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @return A LIGER Object
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' object <- Add_Mito_Ribo_LIGER(liger_object = object, species = "mouse")
#' }
#'

Add_Mito_Ribo_LIGER <- function(
  liger_object,
  species,
  mito_name = "percent_mito",
  ribo_name = "percent_ribo",
  mito_ribo_name = "percent_mito_ribo",
  mito_pattern = NULL,
  ribo_pattern = NULL,
  mito_features = NULL,
  ribo_features = NULL,
  ensembl_ids = FALSE,
  overwrite = FALSE,
  list_species_names = FALSE
) {
  # Accepted species names
  accepted_names <- data.frame(
    Mouse_Options = c("Mouse", "mouse", "Ms", "ms", "Mm", "mm"),
    Human_Options = c("Human", "human", "Hu", "hu", "Hs", "hs"),
    Marmoset_Options = c("Marmoset", "marmoset", "CJ", "Cj", "cj", NA),
    Zebrafish_Options = c("Zebrafish", "zebrafish", "DR", "Dr", "dr", NA),
    Rat_Options = c("Rat", "rat", "RN", "Rn", "rn", NA),
    Drosophila_Options = c("Drosophila", "drosophila", "DM", "Dm", "dm", NA),
    Macaque_Options = c("Macaque", "macaque", "Rhesus", "macaca", "mmulatta", NA)
  )

  # Return list of accepted default species name options
  if (list_species_names) {
    return(accepted_names)
    stop_quietly()
  }

  # LIGER object check
  Is_LIGER(liger_object = liger_object)

  # Overwrite check
  if (mito_name %in% colnames(x = liger_object@cell.data) || ribo_name %in% colnames(x = liger_object@cell.data) || mito_ribo_name %in% colnames(x = liger_object@cell.data)) {
    if (!overwrite) {
      cli_abort(message = c("Columns with {mito_name} and/or {ribo_name} already present in cell.data slot.",
                            "i" = "*To run function and overwrite columns set parameter `overwrite = TRUE` or change respective 'mito_name', 'ribo_name', or mito_ribo_name'*")
      )
    }
    cli_inform(message = c("Columns with {mito_name} and/or {ribo_name} already present in cell.data slot.",
                           "i" = "Overwriting those columns as overwrite = TRUE.")
    )
  }

  # Checks species
  if (is.null(x = species)) {
    cli_abort(message = c("No species name or abbreivation was provided to `species` parameter.",
                          "i" = "If not using default species please set `species = other`.")
    )
  }

  # Species Spelling Options
  mouse_options <- accepted_names$Mouse_Options
  human_options <- accepted_names$Human_Options
  marmoset_options <- accepted_names$Marmoset_Options
  zebrafish_options <- accepted_names$Zebrafish_Options
  rat_options <- accepted_names$Rat_Options
  drosophila_options <- accepted_names$Drosophila_Options
  macaque_options <- accepted_names$Macaque_Options

  # Check ensembl vs patterns
  if (ensembl_ids && species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options) && any(!is.null(x = mito_pattern), !is.null(x = ribo_pattern), !is.null(x = mito_features), !is.null(x = ribo_features))) {
    cli_warn(message = c("When using a default species and setting `ensembl_ids = TRUE` provided patterns or features are ignored.",
                         "*" = "Supplied `mito_pattern`,`ribo_pattern`, `mito_features`,`ribo_features` will be disregarded.")
    )
  }

  # Assign mito/ribo pattern to stored species
  if (species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options) && any(!is.null(x = mito_pattern), !is.null(x = ribo_pattern))) {
    cli_warn(message = c("Pattern expressions for included species (Human & Mouse) are set by default.",
                         "*" = "Supplied `mito_pattern` and `ribo_pattern` will be disregarded.",
                         "i" = "To override defaults please supply a feature list for mito and/or ribo genes.")
    )
  }

  # default patterns or features
  if (species %in% mouse_options) {
    mito_pattern <- "^mt-"
    ribo_pattern <- "^Rp[sl]"
  }
  if (species %in% human_options) {
    mito_pattern <- "^MT-"
    ribo_pattern <- "^RP[SL]"
  }
  if (species %in% c(marmoset_options, macaque_options)) {
    mito_features <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6")
    ribo_pattern <- "^RP[SL]"
  }
  if (species %in% zebrafish_options) {
    mito_pattern <- "^mt-"
    ribo_pattern <- "^rp[sl]"
  }
  if (species %in% rat_options) {
    mito_pattern <- "^Mt-"
    ribo_pattern <- "^Rp[sl]"
  }
  if (species %in% drosophila_options) {
    mito_pattern <- "^mt:"
    ribo_pattern <- "^Rp[SL]"
  }

  # Check that values are provided for mito and ribo
  if (is.null(x = mito_pattern) && is.null(x = mito_features) && is.null(x = ribo_pattern) && is.null(x = ribo_pattern)) {
    cli_abort(message = c("No features or patterns provided for mito/ribo genes.",
                          "i" = "Please provide a default species name or pattern/features."))
  }

  # Retrieve ensembl ids if TRUE
  if (ensembl_ids) {
    mito_features <- Retrieve_Ensembl_Mito(species = species)
    ribo_features <- Retrieve_Ensembl_Ribo(species = species)
  }

  # get features from patterns
  mito_features <- mito_features %||% grep(pattern = mito_pattern, x = rownames(x = liger_object@raw.data[[1]]), value = TRUE)

  ribo_features <- ribo_features %||% grep(pattern = ribo_pattern, x = rownames(x = liger_object@raw.data[[1]]), value = TRUE)

  # Check features are present in object
  length_mito_features <- length(x = intersect(x = mito_features, y = rownames(x = liger_object@raw.data[[1]])))

  length_ribo_features <- length(x = intersect(x = ribo_features, y = rownames(x = liger_object@raw.data[[1]])))

  # Check length of mito and ribo features found in object
  if (length_mito_features < 1 && length_ribo_features < 1) {
    cli_abort(message = c("No Mito or Ribo features found in object using patterns/feature list provided.",
                          "i" = "Please check pattern/feature list and/or gene names in object.")
    )
  }
  if (length_mito_features < 1) {
    cli_warn(message = c("No Mito features found in object using pattern/feature list provided.",
                         "i" = "No column will be added to meta.data.")
    )
  }
  if (length_ribo_features < 1) {
    cli_warn(message = c("No Ribo features found in object using pattern/feature list provided.",
                         "i" = "No column will be added to meta.data.")
    )
  }

  # Add mito and ribo percent
  if (length_mito_features > 0) {
    good_mito <- mito_features[mito_features %in% rownames(x = liger_object@raw.data)]
    percent_mito <- unlist(lapply(liger_object@raw.data, function(x) {
      (Matrix::colSums(x[good_mito, ])/Matrix::colSums(x))*100}))
    liger_object@cell.data[ , mito_name] <- percent_mito
  }

  if (length_ribo_features > 0){
    good_ribo <- ribo_features[ribo_features %in% rownames(x = liger_object@raw.data)]
    percent_ribo <- unlist(lapply(liger_object@raw.data, function(x) {
      (Matrix::colSums(x[good_ribo, ])/Matrix::colSums(x))*100}))
    liger_object@cell.data[ , ribo_name] <- percent_ribo
  }

  # Create combined mito ribo column if both present
  if (length_mito_features > 0 && length_ribo_features > 0) {
    object_meta <- liger_object@cell.data %>%
      rownames_to_column("barcodes")

    object_meta <- object_meta %>%
      mutate(percent_mito_ribo = percent_mito + percent_ribo)

    liger_object@cell.data[ , mito_ribo_name] <- object_meta$percent_mito_ribo
  }

  # return object
  return(liger_object)
}


#' Add Cell Complexity Value
#'
#' Add measure of cell complexity/novelty (log10PerUMI) for data QC.
#'
#' @param liger_object object name.
#' @param meta_col_name name to use for new meta data column.  Default is "log10GenesPerUMI".
#' @param overwrite Logical.  Whether to overwrite existing an meta.data column.  Default is FALSE meaning that
#' function will abort if column with name provided to `meta_col_name` is present in meta.data slot.
#'
#' @import cli
#'
#' @return A LIGER Object
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' object <- Add_Cell_Complexity_Seurat(seurat_object = object)
#' }
#'

Add_Cell_Complexity_LIGER <- function(
  liger_object,
  meta_col_name = "log10GenesPerUMI",
  overwrite = FALSE
) {
  # Check Seurat
  Is_LIGER(liger_object = liger_object)

  # Check columns for overwrite
  if (meta_col_name %in% colnames(x = liger_object@cell.data)) {
    if (!overwrite) {
      cli_abort(message = c("Column '{meta_col_name}' already present in cell.data slot.",
                            "i" = "*To run function and overwrite column, set parameter `overwrite = TRUE` or change respective 'meta_col_name'*.")
      )
    }
    cli_inform(message = c("Column '{meta_col_name}' already present in cell.data slot",
                           "i" = "Overwriting those columns as `overwrite = TRUE`.")
    )
  }

  # Add score
  liger_object@cell.data[ , meta_col_name] <- log10(liger_object@cell.data$nGene) / log10(liger_object@cell.data$nUMI)

  #return object
  return(liger_object)
}




#' Check if meta data are present
#'
#' Check if meta data columns are present in object and return vector of found columns  Return warning
#' messages for meta data columns not found.
#'
#' @param seurat_object object name.
#' @param meta_col_names vector of column names to check.
#' @param print_msg logical. Whether message should be printed if all features are found.  Default is TRUE.
#'
#' @return vector of meta data columns that are present
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' meta_variables <- Meta_Present_LIGER(seurat_object = obj_name, gene_list = DEG_list, print_msg = TRUE)
#' }
#'

Meta_Present_LIGER <- function(
  liger_object,
  meta_col_names,
  print_msg = TRUE
) {
  # Check Seurat
  Is_LIGER(liger_object = liger_object)

  # get all features
  possible_features <- colnames(x = liger_object@cell.data)

  # If any features not found
  if (any(!meta_col_names %in% possible_features)) {
    bad_meta <- meta_col_names[!meta_col_names %in% possible_features]
    found_meta <- meta_col_names[meta_col_names %in% possible_features]

    # Return message of features not found
    if (length(x = found_meta) < 1) {
      stop("No valid meta data column names found. The following @cell.data columns were omitted as they were not found",
           ": ", glue_collapse_scCustom(input_string = bad_meta, and = TRUE))
    }

    if (length(x = bad_meta) > 0) {
      warning("The following @cell.data columns were omitted as they were not found",
              ": ", glue_collapse_scCustom(input_string = bad_meta, and = TRUE))
    }

    # Return the found features omitting the not found ones.
    return(found_meta)
  }

  # Print all found message if TRUE
  if (print_msg) {
    message("All @cell.data columns present.")
  }

  # Return full input gene list.
  return(meta_col_names)
}


#' Extract top loading genes for LIGER factor
#'
#' Extract vector to the top loading genes for specified LIGER iNMF factor
#'
#' @param liger_object LIGER object name.
#' @param liger_factor LIGER factor number to pull genes from.
#' @param num_genes number of top loading genes to return as vector.
#'
#' @return A LIGER Object
#'
#' @import cli
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' top_genes_factor10 <- Top_Genes_Factor(liger_object = object, num_genes = 10)
#' }
#'

Top_Genes_Factor <- function(
  liger_object,
  liger_factor,
  num_genes = 10
) {
  # LIGER object check
  Is_LIGER(liger_object = liger_object)

  # check number of factors present
  if (!liger_factor %in% 1:dim(x = liger_object@W)[[1]]) {
    cli_abort(message = c("'liger_factor' provided: {liger_factor} not found",
                          "i" = "'liger_object' only contains {dim(x = liger_object@W)[[1]]} factors.")
    )
  }

  # Extract genes
  W <- t(liger_object@W)
  rownames(W) <- colnames(liger_object@scale.data[[1]])
  top_genes <- rownames(W)[order(W[, liger_factor], decreasing = TRUE)[1:num_genes]]
  return(top_genes)
}


#' DimPlot LIGER Version
#'
#' Standard and modified version of LIGER's plotByDatasetAndCluster
#'
#' @param liger_object Name of LIGER object.  Need to perform clustering before calling this function.
#' @param clusters Another clustering to use for coloring second plot (must have same names as
#' clusters slot) (default NULL).
#' @param shuffle Randomly shuffle points so that points from same dataset are not plotted one after
#' the other (default TRUE).
#' @param shuffle_seed Random seed for reproducibility of point shuffling (default 1).
#' @param redorder.idents logical whether to reorder the datasets from default order before plotting (default FALSE).
#' @param new.order new dataset factor order for plotting.  must set reorder.idents = TRUE.
#' @param group_by meta data varibale to group plots by
#' @param split_by meta data variable to splot plots by
#'
#' @return A data.frame with information for plotting
#'
#' @references This function is encompasses the first part of the LIGER function plotByDatasetAndCluster.
#' However, this function is modified to allow plotting other meta data variables.  In this case the function
#' just returns the data.frame needed for plotting rather than plots themselves.
#' (https://github.com/welch-lab/liger). (Licence: GPL-3).
#'
#' @noRd
#'
#' @concept liger_plotting_util
#'

Generate_Plotting_df_LIGER <- function(object,
                                       clusters = NULL,
                                       shuffle = TRUE,
                                       shuffle_seed = 1,
                                       reorder.idents = FALSE,
                                       new.order = NULL,
                                       group_by = "dataset",
                                       split_by = NULL
) {
  tsne_df <- data.frame(object@tsne.coords)
  colnames(tsne_df) <- c("tsne1", "tsne2")
  tsne_df[[group_by]] <- object@cell.data[[group_by]]
  if (!is.null(x = split_by)) {
    tsne_df[[split_by]] <- object@cell.data[[split_by]]
  }

  if (reorder.idents == TRUE){
    tsne_df[[group_by]]  <- factor(tsne_df[[group_by]], levels = new.order)
  }
  c_names <- names(object@clusters)
  if (is.null(clusters)) {
    # if clusters have not been set yet
    if (length(object@clusters) == 0) {
      clusters <- rep(1, nrow(object@tsne.coords))
      names(clusters) <- c_names <- rownames(object@tsne.coords)
    } else {
      clusters <- object@clusters
      c_names <- names(object@clusters)
    }
  }
  tsne_df[['Cluster']] <- clusters[c_names]

  if (shuffle) {
    set.seed(shuffle_seed)
    idx <- sample(1:nrow(tsne_df))
    tsne_df <- tsne_df[idx, ]
  }
  return(tsne_df)
}


#' LIGER plot by cluster.
#'
#' Modified version of LIGER's plotByDatasetAndCluster just for plotting clusters.
#'
#' @param liger_object Name of LIGER object.  Need to perform clustering before calling this function.
#' @param colors_use colors to use for plotting by cluster.  By default if number of levels plotted is
#' less than or equal to 36 it will use "polychrome" and if greater than 36 will use "varibow" with
#' shuffle = TRUE both from \code{\link{DiscretePalette_scCustomize}}.
#' @param group_by Variable to be plotted.  If `NULL` will plot clusters from `liger@clusters` slot.
#' If `combination = TRUE` will plot both clusters and meta data variable.
#' @param split_by meta data variable to split plots by (i.e. "dataset").
#' @param title plot title.
#' @param pt_size Adjust point size for plotting.
#' @param reduction_label What to label the x and y axes of resulting plots.  LIGER does not store
#' name of technique and therefore needs to be set manually.  Default is "UMAP".
#' @param num_columns Number of columns to plot by if `split_by` is not NULL.
#' @param shuffle logical. Whether to randomly shuffle the order of points. This can be useful for
#' crowded plots if points of interest are being buried. (Default is TRUE).
#' @param shuffle_seed Sets the seed if randomly shuffling the order of points.
#' @param legend.size what to set legend size to.
#' @param label logical.  Whether or not to label the clusters.  Default is TRUE.
#' @param label_size size of cluster labels.
#' @param label_repel logical.  Whether to repel cluster labels from each other if plotting by
#' cluster (if `group_by = NULL` or `group_by = "cluster`).  Default is FALSE.
#' @param label_box logical.  Whether to put a box around the label text (uses `geom_text` vs `geom_label`).
#' Default is FALSE.
#' @param label_color Color to use for cluster labels.  Default is "black".
#' @param redorder.idents logical. should the idents plotted by reordered.  Default is FALSE.
#' @param new.order What should the new ident order be if `reorder.idents = TRUE`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#'
#' @return A ggplot/patchwork object
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel geom_label_repel
#' @importFrom cowplot theme_cowplot
#' @importFrom patchwork wrap_plots
#'
#' @references This function is encompasses part of the LIGER function plotByDatasetAndCluster.
#' However, this function is modified to just return cluster plots based on `Generate_Plotting_df_LIGER`.
#' (https://github.com/welch-lab/liger). (Licence: GPL-3).
#'
#' @noRd
#'
#' @concept liger_plotting_util
#'

Plot_By_Cluster_LIGER <- function(
  liger_object,
  colors_use = NULL,
  group_by = "dataset",
  split_by = NULL,
  title = NULL,
  pt_size = NULL,
  reduction_label = "UMAP",
  num_columns = NULL,
  shuffle = TRUE,
  shuffle_seed = 1,
  legend.size = 5,
  label = TRUE,
  label_size = NA,
  label_repel = FALSE,
  label_box = FALSE,
  label_color = "black",
  reorder.idents = FALSE,
  new.order = NULL,
  raster = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123
) {
  # Create plotting data.frame
  tsne_df <- Generate_Plotting_df_LIGER(object = liger_object, group_by = group_by, split_by = split_by, reorder.idents = reorder.idents, shuffle = shuffle, shuffle_seed = shuffle_seed)

  if (!is.null(x = split_by)) {
    list_of_splits <- unique(tsne_df[[split_by]])
  }

  # Get length of meta data feature
  if (!is.null(x = split_by) && !is.null(x = num_columns)) {
    split.by_length <- length(list_of_splits)

    # Calculate number of rows for selected number of columns
    num_rows <- ceiling(split.by_length/num_columns)

    # Check column and row compatibility
    if (num_columns > split.by_length) {
      cli_abort(message = c("The number of columns specified is greater than the number of meta data variables.",
                            "*" = "'{split_by}' only contains: {split.by_length} variables.",
                            "i" = "Please adjust `num_columns` to be less than or equal t: {split.by_length}.")
      )
    }
  }

  centers <- tsne_df %>% group_by(.data[['Cluster']]) %>% summarize(
    tsne1 = median(x = .data[['tsne1']]),
    tsne2 = median(x = .data[['tsne2']])
  )

  cluster_length <- length(x = unique(x = liger_object@clusters))

  if (is.null(x = colors_use)) {
    # set default plot colors
    if (is.null(x = colors_use)) {
      colors_use <- scCustomize_Palette(num_groups = cluster_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
    }
  }

  # Create accurate axis labels
  x_axis_label <- paste0(reduction_label, "_1")
  y_axis_label <- paste0(reduction_label, "_2")

  # plot
  if (raster) {
    if (!is.null(x = split_by)) {
      p2 <- lapply(1:length(x = list_of_splits), function(x){
        p2 <- ggplot(subset(tsne_df, tsne_df[[split_by]] %in% list_of_splits[x]), aes_string(x = 'tsne1', y = 'tsne2', color = 'Cluster')) +
          theme_cowplot() +
          geom_scattermore(pointsize = pt_size) +
          guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(list_of_splits[x]) +
          scale_color_manual(values = colors_use) +
          theme(legend.position = "right",
                axis.text = element_text(size = rel(0.95)),
                plot.title = element_text(hjust = 0.5)) +
          guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
          xlab(x_axis_label) +
          ylab(y_axis_label)

        if (label_box) {
          geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes_string(label = 'Cluster', fill = 'Cluster'), size = label_size,
            show.legend = FALSE, color = label_color
          ) + scale_fill_manual(values = colors_use)
        } else if (label) {
          geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes_string(label = 'Cluster'), size = label_size, color = label_color,
            show.legend = FALSE
          )
        } else {
          p2 <- p2
        }
      })
    } else {
      p2 <- ggplot(tsne_df, aes_string(x = 'tsne1', y = 'tsne2', color = 'Cluster')) +
        theme_cowplot() +
        geom_scattermore(pointsize = pt_size) +
        guides(color = guide_legend(override.aes = list(size = legend.size))) +
        scale_color_manual(values = colors_use) +
        theme(legend.position = "right",
              axis.text = element_text(size = rel(0.95)),
              plot.title = element_text(hjust = 0.5)) +
        guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
        xlab(x_axis_label) +
        ylab(y_axis_label)

      if (label_box) {
        geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes_string(label = 'Cluster', fill = 'Cluster'), size = label_size,
          show.legend = FALSE, color = label_color
        ) + scale_fill_manual(values = colors_use)
      } else if (label) {
        geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes_string(label = 'Cluster'), size = label_size, color = label_color,
          show.legend = FALSE
        )
      } else {
        p2 <- p2
      }

    }
  } else {
    if (!is.null(x = split_by)) {
      p2 <- lapply(1:length(x = list_of_splits), function(x){
        p2 <- ggplot(subset(tsne_df, tsne_df[[split_by]] %in% list_of_splits[x]),aes_string(x = 'tsne1', y = 'tsne2', color = 'Cluster')) +
          theme_cowplot() +
          geom_point(size = pt_size) +
          guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(list_of_splits[x]) +
          scale_color_manual(values = colors_use) +
          theme(legend.position = "right",
                axis.text = element_text(size = rel(0.95)),
                plot.title = element_text(hjust = 0.5)) +
          guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
          xlab(x_axis_label) +
          ylab(y_axis_label)

        if (label_box) {
          geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes_string(label = 'Cluster', fill = 'Cluster'), size = label_size,
            show.legend = FALSE, color = label_color
          ) + scale_fill_manual(values = colors_use)
        } else if (label) {
          geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes_string(label = 'Cluster'), size = label_size, color = label_color,
            show.legend = FALSE
          )
        } else {
          p2 <- p2
        }
      })
    } else {
      p2 <- ggplot(tsne_df, aes_string(x = 'tsne1', y = 'tsne2', color = 'Cluster')) +
        theme_cowplot() +
        geom_point(size = pt_size) +
        guides(color = guide_legend(override.aes = list(size = legend.size))) +
        scale_color_manual(values = colors_use) +
        theme(legend.position = "right",
              axis.text = element_text(size = rel(0.95)),
              plot.title = element_text(hjust = 0.5)) +
        guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
        xlab(x_axis_label) +
        ylab(y_axis_label)

      if (label_box) {
        geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes_string(label = 'Cluster', fill = 'Cluster'), size = label_size,
          show.legend = FALSE, color = label_color
        ) + scale_fill_manual(values = colors_use)
      } else if (label) {
        geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes_string(label = 'Cluster'), size = label_size, color = label_color,
          show.legend = FALSE
        )
      } else {
        p2 <- p2
      }
    }
  }
  if (!is.null(x = split_by) && !is.null(x = num_columns)) {
    p2 <- wrap_plots(p2) + plot_layout(nrow = num_rows, ncol = num_columns, guides = 'collect')
    return(p2)
  }
  if (!is.null(x = split_by) && is.null(x = num_columns)) {
    p2 <- wrap_plots(p2) + plot_layout(guides = 'collect')
    return(p2)
  } else {
    return(p2)
  }
}

#' LIGER plot by meta variables.
#'
#' Modified version of LIGER's plotByDatasetAndCluster just for plotting meta variables.
#'
#' @param liger_object Name of LIGER object.  Need to perform clustering before calling this function.
#' @param colors_use colors to use for plotting by cluster.  By default if number of levels plotted is
#' less than or equal to 36 it will use "polychrome" and if greater than 36 will use "varibow" with
#' shuffle = TRUE both from \code{\link{DiscretePalette_scCustomize}}.
#' @param group_by Variable to be plotted.  If `NULL` will plot clusters from `liger@clusters` slot.
#' If `combination = TRUE` will plot both clusters and meta data variable.
#' @param split_by meta data variable to split plots by (i.e. "dataset").
#' @param title plot title.
#' @param pt_size Adjust point size for plotting.
#' @param reduction_label What to label the x and y axes of resulting plots.  LIGER does not store name
#' of technique and therefore needs to be set manually.  Default is "UMAP".
#' @param num_columns Number of columns to plot by if `split_by` is not NULL.
#' @param shuffle logical. Whether to randomly shuffle the order of points. This can be useful for
#' crowded plots if points of interest are being buried. (Default is TRUE).
#' @param shuffle_seed Sets the seed if randomly shuffling the order of points.
#' @param legend.size what to set legend size to.
#' @param redorder.idents logical. should the idents plotted by reordered.  Default is FALSE.
#' @param new.order What should the new ident order be if `reorder.idents = TRUE`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#'
#' @return A ggplot/patchwork object
#'
#' @import ggplot2
#' @importFrom cowplot theme_cowplot
#' @importFrom patchwork wrap_plots
#'
#' @references This function is encompasses part of the LIGER function plotByDatasetAndCluster.
#' However, this function is modified to just return cluster plots based on `Generate_Plotting_df_LIGER`.
#' (https://github.com/welch-lab/liger). (Licence: GPL-3).
#'
#' @noRd
#'
#' @concept liger_plotting_util
#'

Plot_By_Meta_LIGER <- function(
  liger_object,
  colors_use = NULL,
  group_by = "dataset",
  split_by = NULL,
  title = NULL,
  pt_size = NULL,
  reduction_label = "UMAP",
  num_columns = NULL,
  shuffle = TRUE,
  shuffle_seed = 1,
  legend.size = 3,
  reorder.idents = FALSE,
  new.order = NULL,
  raster = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123
) {

  tsne_df <- Generate_Plotting_df_LIGER(object = liger_object, group_by = group_by, split_by = split_by, reorder.idents = reorder.idents, shuffle = shuffle, shuffle_seed = shuffle_seed)

  if (!is.null(x = split_by)) {
    list_of_splits <- unique(tsne_df[[split_by]])
  }

  # Get length of meta data feature
  if (!is.null(x = split_by) && !is.null(x = num_columns)) {
    split.by_length <- length(list_of_splits)

    # Calculate number of rows for selected number of columns
    num_rows <- ceiling(split.by_length/num_columns)

    # Check column and row compatibility
    if (num_columns > split.by_length) {
      cli_abort(message = c("The number of columns specified is greater than the number of meta data variables.",
                            "*" = "'{split_by}' only contains: {split.by_length} variables.",
                            "i" = "Please adjust `num_columns` to be less than or equal t: {split.by_length}.")
      )
    }
  }

  meta_length <- length(x = unique(x = liger_object@cell.data[[group_by]]))

  if (is.null(x = colors_use)) {
    # set default plot colors
    if (is.null(x = colors_use)) {
      colors_use <- scCustomize_Palette(num_groups = meta_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
    }
  }

  # Create accurate axis labels
  x_axis_label <- paste0(reduction_label, "_1")
  y_axis_label <- paste0(reduction_label, "_2")

  if (raster) {
    if (!is.null(x = split_by)) {
      p1 <- lapply(1:length(x = list_of_splits), function(x){
        ggplot(subset(tsne_df, tsne_df[[split_by]] %in% list_of_splits[x]), aes_string(x = 'tsne1', y = 'tsne2', color = group_by)) +
          theme_cowplot() +
          geom_scattermore(pointsize = pt_size) +
          guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(list_of_splits[x]) +
          scale_color_manual(values = colors_use) +
          theme(legend.position = "right",
                axis.text = element_text(size = rel(0.95)),
                plot.title = element_text(hjust = 0.5)) +
          guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
          xlab(x_axis_label) +
          ylab(y_axis_label)
      })
    } else {
      p1 <- ggplot(tsne_df, aes_string(x = 'tsne1', y = 'tsne2', color = group_by)) +
        theme_cowplot() +
        geom_scattermore(pointsize = pt_size) +
        guides(color = guide_legend(override.aes = list(size = legend.size))) +
        scale_color_manual(values = colors_use) +
        theme(legend.position = "right",
              axis.text = element_text(size = rel(0.95)),
              plot.title = element_text(hjust = 0.5)) +
        guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
        xlab(x_axis_label) +
        ylab(y_axis_label)

    }
  } else {
    if (!is.null(x = split_by)) {
      p1 <- lapply(1:length(x = list_of_splits), function(x){
        ggplot(subset(tsne_df, tsne_df[[split_by]] %in% list_of_splits[x]),aes_string(x = 'tsne1', y = 'tsne2', color = group_by)) +
          theme_cowplot() +
          geom_point(size = pt_size) +
          guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(list_of_splits[x]) +
          scale_color_manual(values = colors_use) +
          theme(legend.position = "right",
                axis.text = element_text(size = rel(0.95)),
                plot.title = element_text(hjust = 0.5)) +
          guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
          xlab(x_axis_label) +
          ylab(y_axis_label)
      })
    } else {
      p1 <- ggplot(tsne_df, aes_string(x = 'tsne1', y = 'tsne2', color = group_by)) +
        theme_cowplot() +
        geom_point(size = pt_size) +
        guides(color = guide_legend(override.aes = list(size = legend.size))) +
        scale_color_manual(values = colors_use) +
        theme(legend.position = "right",
              axis.text = element_text(size = rel(0.95)),
              plot.title = element_text(hjust = 0.5)) +
        guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
        xlab(x_axis_label) +
        ylab(y_axis_label)
    }
  }
  if (!is.null(x = split_by) && !is.null(x = num_columns)) {
    p1 <- wrap_plots(p1) + plot_layout(nrow = num_rows, ncol = num_columns)
    return(p1)
  }
  if (!is.null(x = split_by) && is.null(x = num_columns)) {
    p1 <- wrap_plots(p1)
    return(p1)
  } else {
    return(p1)
  }
}


#' Perform variable gene selection over whole dataset
#'
#' Performs variable gene selection for LIGER object across the entire object instead of by
#' dataset and then taking union.
#'
#' @param liger_object LIGER object name.
#' @param num_genes Number of genes to find. Optimizes the value of `var.thresh`  to get
#' this number of genes, (Default is NULL).
#' @param var.thresh Variance threshold. Main threshold used to identify variable genes.
#' Genes with expression variance greater than threshold (relative to mean) are selected.
#' (higher threshold -> fewer selected genes).
#' @param alpha.thresh Alpha threshold. Controls upper bound for expected mean gene
#' expression (lower threshold -> higher upper bound). (default 0.99)
#' @param tol Tolerance to use for optimization if num.genes values passed in (default 0.0001).
#' @param do.plot Display log plot of gene variance vs. gene expression. Selected genes are
#' plotted in green. (Default FALSE)
#' @param pt.size Point size for plot.
#' @param chunk size of chunks in hdf5 file. (Default 1000)
#'
#' @return A LIGER Object with variable genes in correct slot.
#'
#' @import cli
#' @import rliger
#'
#' @references Matching function parameter text descriptions are taken from `rliger::selectGenes`
#' which is called by this function after creating new temporary object/dataset.
#' (https://github.com/welch-lab/liger). (Licence: GPL-3).
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' liger_obj <- Variable_Features_LIGER(liger_object = liger_obj, num_genes = 2000)
#' }
#'

Variable_Features_LIGER <- function(
  liger_object,
  num_genes = NULL,
  var.thresh = 0.3,
  alpha.thresh = 0.99,
  tol = 0.0001,
  do.plot = FALSE,
  cex.use = 0.3,
  chunk=1000
) {
  Is_LIGER(liger_object = liger_object)

  raw_data <- liger_object@raw.data

  cli_inform(message = "Creating temporary object with combined data.")

  temp_liger <- createLiger(raw.data = list("dataset" = Merge_Sparse_Data_All(raw_data)), remove.missing = FALSE)

  rm(raw_data)
  gc()

  cli_inform(message = "Normalizing and identifying variable features.")

  temp_liger <- rliger::normalize(temp_liger)
  temp_liger <- selectGenes(object = temp_liger, var.thresh = var.thresh, do.plot = do.plot, num.genes = num_genes, tol = tol, alpha.thresh = alpha.thresh, cex.use = pt.size, chunk = chunk)
  var_genes <- temp_liger@var.genes

  rm(temp_liger)
  gc()

  liger_object@var.genes <- var_genes
  return(liger_object)
}
