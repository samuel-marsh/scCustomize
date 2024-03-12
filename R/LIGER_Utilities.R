#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### DATA ACCESS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname Fetch_Meta
#' @importFrom methods slot
#' @export
#' @concept liger_object_util
#' @method Fetch_Meta liger

Fetch_Meta.liger <- function(
    object,
    ...
) {
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    object_meta <- rliger2::cellMeta(x = object, as.data.frame = TRUE)
  } else {
    object_meta <- object_meta <- slot(object = object, name = "cell.data")
  }

  # return meta
  return(object_meta)
}


#' Extract Features from LIGER Object
#'
#' Extract all unique features from LIGER object
#'
#' @param liger_object LIGER object name.
#' @param by_dataset logical, whether to return list with vector of features for each dataset in
#' LIGER object or to return single vector of unique features across all datasets in object
#' (default is FALSE; return vector of unique features)
#'
#' @return vector or list depending on `by_dataset` parameter
#'
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # return single vector of all unique features
#' all_features <- LIGER_Features(liger_object = object, by_dataset = FALSE)
#'
#' # return list of vectors containing features from each individual dataset in object
#' dataset_features <- LIGER_Features(liger_object = object, by_dataset = TRUE)
#' }
#'

LIGER_Features <- function(
    liger_object,
    by_dataset = FALSE
) {
  # check liger
  Is_LIGER(liger_object = liger_object)

  # liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    # Extract features
    features_by_dataset <- lapply(1:length(x = liger_object@datasets), function(x) {
      rownames(x = liger_object@datasets[[x]]@featureMeta)
    })

  } else {
    # Extract features
    features_by_dataset <- lapply(1:length(x = liger_object@raw.data), function(x) {
      rownames(x = liger_object@raw.data[[x]])
    })
  }

  # Return features
  if (isFALSE(x = by_dataset)) {
    features <- unique(x = unlist(x = features_by_dataset))
    return(features)
  } else {
    return(features_by_dataset)
  }
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
#' @importFrom utils packageVersion
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
  # temp liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    cli_abort(message = c("Liger functionality is currently restricted to rliger v1.0.1 or lower.",
                          "i" = "Functionality with rliger v2+ is currently in development."))
  }

  # LIGER object check
  Is_LIGER(liger_object = liger_object)

  # check number of factors present
  if (!liger_factor %in% 1:dim(x = liger_object@W)[[1]]) {
    cli_abort(message = c("{.code liger_factor} provided: {.field {liger_factor}} not found",
                          "i" = "{.code liger_object} only contains {.field {dim(x = liger_object@W)[[1]]}} factors.")
    )
  }

  # Extract genes
  W <- t(liger_object@W)
  rownames(x = W) <- colnames(x = liger_object@scale.data[[1]])
  top_genes <- rownames(x = W)[order(W[, liger_factor], decreasing = TRUE)[1:num_genes]]
  return(top_genes)
}


#' Extract dimensionality reduction coordinates from Liger object
#'
#' Extract data.frame containing dimensionality reduction coordinates from new format of
#' Liger objects
#'
#' @param liger_object LIGER object name.
#' @param reduction name of dimensionality reduction stored in cellMeta slot.  Default is
#' "UMAP")
#' @param check_only logical, return `TRUE` if valid reduction is present.
#'
#' @return dimensionality reduction coordinates in 2 column format
#'
#' @import cli
#' @importFrom methods slotNames
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' umap_coords <- LIGER_DimReduc(liger_object = object)
#' }
#'

LIGER_DimReduc <- function(
    liger_object,
    reduction = NULL,
    check_only = FALSE
) {
  # Check new liger object
  if (!"cellMeta" %in% slotNames(liger_object)) {
    cli_abort(message = "This function is only for objects created with rliger >= v2.0.0")
  }

  # check reduction in cellMeta
  if (reduction %in% names(x = liger_object@cellMeta)) {
    # check the right dims
    if (length(dim(liger_object@cellMeta[[reduction]])) != 2) {
      cli_abort(message = "The cellMeta entry {.field {reduction}} is not 2-dimensional entry.")
    } else {
      if (isTRUE(x = check_only)) {
        return(TRUE)
      }
      # get coords
      reduc_coords <- liger_object@cellMeta[[reduction]]

      # add colnames
      colnames(reduc_coords) <- paste0(reduction, "_", 1:2)
    }
  } else {
    cli_abort("The reduction {.field {reduction}} is not present in cellMeta slot.")
  }

  # return coords
  return(reduc_coords)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### QC UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Add Mito and Ribo percentages
#'
#' @param species Species of origin for given Seurat Object.  If mouse, human, marmoset, zebrafish, rat,
#' drosophila, or rhesus macaque (name or abbreviation) are provided the function will automatically
#' generate mito_pattern and ribo_pattern values.
#' @param mito_name name to use for the new meta.data column containing percent mitochondrial counts.
#' Default is "percent_mito".
#' @param ribo_name name to use for the new meta.data column containing percent ribosomal counts.
#' Default is "percent_ribo".
#' @param mito_ribo_name name to use for the new meta.data column containing percent mitochondrial+ribosomal
#' counts.  Default is "percent_mito_ribo".
#' @param mito_pattern A regex pattern to match features against for mitochondrial genes (will set automatically
#' if species is mouse or human; marmoset features list saved separately).
#' @param ribo_pattern A regex pattern to match features against for ribosomal genes (will set automatically
#' if species is mouse, human, or marmoset).
#' @param mito_features A list of mitochondrial gene names to be used instead of using regex pattern.
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
#' @importFrom rlang ":="
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom utils packageVersion
#'
#' @method Add_Mito_Ribo liger
#'
#' @export
#' @rdname Add_Mito_Ribo
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # Liger
#' liger_object <- Add_Mito_Ribo(object = liger_object, species = "human")
#' }
#'

Add_Mito_Ribo.liger <- function(
  object,
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
  list_species_names = FALSE,
  ...
) {
  # temp liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    cli_abort(message = c("Liger functionality is currently restricted to rliger v1.0.1 or lower.",
                          "i" = "Functionality with rliger v2+ is currently in development."))
  }

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
  if (isTRUE(x = list_species_names)) {
    return(accepted_names)
    stop_quietly()
  }

  # LIGER object check
  Is_LIGER(liger_object = object)

  # Check name collision
  if (any(duplicated(x = c(mito_name, ribo_name, mito_ribo_name)))) {
    cli_abort(message = "One or more of values provided to {.code mito_name}, {.code ribo_name}, {.code mito_ribo_name} are identical.")
  }

  # Overwrite check
  if (mito_name %in% colnames(x = object@cell.data) || ribo_name %in% colnames(x = object@cell.data) || mito_ribo_name %in% colnames(x = object@cell.data)) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Columns with {.val {mito_name}} and/or {.val {ribo_name}} already present in cell.data slot.",
                            "i" = "*To run function and overwrite columns set parameter {.code overwrite = TRUE} or change respective {.code mito_name}, {.code ribo_name}, and/or {.code mito_ribo_name}.*")
      )
    }
    cli_inform(message = c("Columns with {.val {mito_name}} and/or {.val {ribo_name}} already present in cell.data slot.",
                           "i" = "Overwriting those columns as {.code overwrite = TRUE}.")
    )
  }

  # Checks species
  if (is.null(x = species)) {
    cli_abort(message = c("No species name or abbreivation was provided to {.code species} parameter.",
                          "i" = "If not using default species please set {.code species = other}.")
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
  if (isTRUE(x = ensembl_ids) && species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options) && any(!is.null(x = mito_pattern), !is.null(x = ribo_pattern), !is.null(x = mito_features), !is.null(x = ribo_features))) {
    cli_warn(message = c("When using a default species and setting {.code ensembl_ids = TRUE} provided patterns or features are ignored.",
                         "*" = "Supplied {.code mito_pattern}, {.code ribo_pattern}, {.code mito_features}, {.code ribo_features} will be disregarded.")
    )
  }

  # Assign mito/ribo pattern to stored species
  if (species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options) && any(!is.null(x = mito_pattern), !is.null(x = ribo_pattern))) {
    cli_warn(message = c("Pattern expressions for included species are set by default.",
                         "*" = "Supplied {.code mito_pattern} and {.code ribo_pattern} will be disregarded.",
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
  if (is.null(x = mito_pattern) && is.null(x = mito_features) && is.null(x = ribo_pattern) && is.null(x = ribo_features)) {
    cli_abort(message = c("No features or patterns provided for mito/ribo genes.",
                          "i" = "Please provide a default species name or pattern/features."))
  }

  # Retrieve ensembl ids if TRUE
  if (isTRUE(x = ensembl_ids)) {
    mito_features <- Retrieve_Ensembl_Mito(species = species)
    ribo_features <- Retrieve_Ensembl_Ribo(species = species)
  }

  all_features <- LIGER_Features(liger_object = object)

  # get features from patterns
  mito_features <- mito_features %||% grep(pattern = mito_pattern, x = all_features, value = TRUE)

  ribo_features <- ribo_features %||% grep(pattern = ribo_pattern, x = all_features, value = TRUE)

  # Check features are present in object
  length_mito_features <- length(x = intersect(x = mito_features, y = all_features))

  length_ribo_features <- length(x = intersect(x = ribo_features, y = all_features))

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
    good_mito <- mito_features[mito_features %in% all_features]
    percent_mito <- unlist(lapply(object@raw.data, function(x) {
      (Matrix::colSums(x[good_mito, ])/Matrix::colSums(x))*100}))
    object@cell.data[ , mito_name] <- percent_mito
  }

  if (length_ribo_features > 0){
    good_ribo <- ribo_features[ribo_features %in% all_features]
    percent_ribo <- unlist(lapply(object@raw.data, function(x) {
      (Matrix::colSums(x[good_ribo, ])/Matrix::colSums(x))*100}))
    object@cell.data[ , ribo_name] <- percent_ribo
  }

  # Create combined mito ribo column if both present
  if (length_mito_features > 0 && length_ribo_features > 0) {
    object_meta <- Fetch_Meta(object = object) %>%
      rownames_to_column("barcodes")

    object_meta <- object_meta %>%
      mutate({{mito_ribo_name}} := .data[[mito_name]] + .data[[ribo_name]])

    object@cell.data[ , mito_ribo_name] <- object_meta[[mito_ribo_name]]
  }

  # return object
  return(object)
}


#' Add Cell Complexity Value
#'
#' @param meta_col_name name to use for new meta data column.  Default is "log10GenesPerUMI".
#' @param overwrite Logical.  Whether to overwrite existing an meta.data column.  Default is FALSE meaning that
#' function will abort if column with name provided to `meta_col_name` is present in meta.data slot.
#'
#' @import cli
#' @importFrom utils packageVersion
#'
#' @method Add_Cell_Complexity liger
#'
#' @export
#' @rdname Add_Cell_Complexity
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # Liger
#' liger_object <- Add_Cell_Complexity(object = liger_object)
#' }
#'

Add_Cell_Complexity.liger <- function(
  object,
  meta_col_name = "log10GenesPerUMI",
  overwrite = FALSE,
  ...
) {
  # temp liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    cli_abort(message = c("Liger functionality is currently restricted to rliger v1.0.1 or lower.",
                          "i" = "Functionality with rliger v2+ is currently in development."))
  }

  # Check Seurat
  Is_LIGER(liger_object = object)

  # Check columns for overwrite
  if (meta_col_name %in% colnames(x = object@cell.data)) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Column {.val {meta_col_name}} already present in cell.data slot.",
                            "i" = "*To run function and overwrite column, set parameter {.code overwrite = TRUE} or change respective {.code meta_col_name}*.")
      )
    }
    cli_inform(message = c("Column {.val {meta_col_name}} already present in cell.data slot",
                           "i" = "Overwriting those columns as `overwrite = TRUE`.")
    )
  }

  # Add score
  object@cell.data[ , meta_col_name] <- log10(object@cell.data$nGene) / log10(object@cell.data$nUMI)

  #return object
  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### PLOTTING UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
#' @importFrom utils packageVersion
#'
#' @references This function is encompasses the first part of the LIGER function plotByDatasetAndCluster.
#' However, this function is modified to allow plotting other meta data variables.  In this case the function
#' just returns the data.frame needed for plotting rather than plots themselves.
#' \url{https://github.com/welch-lab/liger}. (License: GPL-3).
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
  # temp liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    cli_abort(message = c("Liger functionality is currently restricted to rliger v1.0.1 or lower.",
                          "i" = "Functionality with rliger v2+ is currently in development."))
  }

  tsne_df <- data.frame(object@tsne.coords)
  colnames(x = tsne_df) <- c("tsne1", "tsne2")
  tsne_df[[group_by]] <- object@cell.data[[group_by]]
  if (!is.null(x = split_by)) {
    tsne_df[[split_by]] <- object@cell.data[[split_by]]
  }

  if (isTRUE(x = reorder.idents)) {
    tsne_df[[group_by]]  <- factor(x = tsne_df[[group_by]], levels = new.order)
  }
  c_names <- names(x = object@clusters)
  if (is.null(x = clusters)) {
    # if clusters have not been set yet
    if (length(x = object@clusters) == 0) {
      clusters <- rep(1, nrow(x = object@tsne.coords))
      names(x = clusters) <- c_names <- rownames(x = object@tsne.coords)
    } else {
      clusters <- object@clusters
      c_names <- names(x = object@clusters)
    }
  }
  tsne_df[['Cluster']] <- clusters[c_names]

  if (isTRUE(x = shuffle)) {
    set.seed(shuffle_seed)
    idx <- sample(x = 1:nrow(tsne_df))
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
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#'
#' @return A ggplot/patchwork object
#'
#' @import ggplot2
#' @importFrom cowplot theme_cowplot
#' @importFrom dplyr summarize
#' @importFrom ggrepel geom_text_repel geom_label_repel
#' @importFrom patchwork wrap_plots
#' @importFrom scattermore geom_scattermore
#' @importFrom stats median
#' @importFrom utils packageVersion
#'
#' @references This function is encompasses part of the LIGER function plotByDatasetAndCluster.
#' However, this function is modified to just return cluster plots based on `Generate_Plotting_df_LIGER`.
#' \url{https://github.com/welch-lab/liger}. (Licence: GPL-3).
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
  raster.dpi = c(512, 512),
  ggplot_default_colors = FALSE,
  color_seed = 123
) {
  # temp liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    cli_abort(message = c("Liger functionality is currently restricted to rliger v1.0.1 or lower.",
                          "i" = "Functionality with rliger v2+ is currently in development."))
  }

  # Create plotting data.frame
  tsne_df <- Generate_Plotting_df_LIGER(object = liger_object, group_by = group_by, split_by = split_by, reorder.idents = reorder.idents, shuffle = shuffle, shuffle_seed = shuffle_seed)

  if (!is.null(x = split_by)) {
    list_of_splits <- unique(x = tsne_df[[split_by]])
  }

  # Get length of meta data feature
  if (!is.null(x = split_by) && !is.null(x = num_columns)) {
    split.by_length <- length(x = list_of_splits)

    # Calculate number of rows for selected number of columns
    num_rows <- ceiling(x = split.by_length/num_columns)

    # Check column and row compatibility
    if (num_columns > split.by_length) {
      cli_abort(message = c("The number of columns specified is greater than the number of meta data variables.",
                            "*" = "{.field {split_by}} only contains: {.field {split.by_length}} variables.",
                            "i" = "Please adjust {.code num_columns} to be less than or equal to: {.field {split.by_length}}.")
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
  if (isTRUE(x = raster)) {
    if (!is.null(x = split_by)) {
      p2 <- lapply(1:length(x = list_of_splits), function(x){
        p2 <- ggplot(subset(tsne_df, tsne_df[[split_by]] %in% list_of_splits[x]), aes(x = .data[['tsne1']], y = .data[['tsne2']], color = .data[['Cluster']])) +
          theme_cowplot() +
          geom_scattermore(pointsize = pt_size, pixels = raster.dpi) +
          guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(list_of_splits[x]) +
          scale_color_manual(values = colors_use) +
          theme(legend.position = "right",
                axis.text = element_text(size = rel(0.95)),
                plot.title = element_text(hjust = 0.5)) +
          guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
          xlab(x_axis_label) +
          ylab(y_axis_label)

        if (isTRUE(x = label_box)) {
          geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes(label = .data[['Cluster']], fill = .data[['Cluster']]), size = label_size,
            show.legend = FALSE, color = label_color
          ) + scale_fill_manual(values = colors_use)
        } else if (isTRUE(x = label)) {
          geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes(label = .data[['Cluster']]), size = label_size, color = label_color,
            show.legend = FALSE
          )
        } else {
          p2 <- p2
        }
      })
    } else {
      p2 <- ggplot(tsne_df, aes(x = .data[['tsne1']], y = .data[['tsne2']], color = .data[['Cluster']])) +
        theme_cowplot() +
        geom_scattermore(pointsize = pt_size, pixels = raster.dpi) +
        guides(color = guide_legend(override.aes = list(size = legend.size))) +
        scale_color_manual(values = colors_use) +
        theme(legend.position = "right",
              axis.text = element_text(size = rel(0.95)),
              plot.title = element_text(hjust = 0.5)) +
        guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
        xlab(x_axis_label) +
        ylab(y_axis_label)

      if (isTRUE(x = label_box)) {
        geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes(label = .data[['Cluster']], fill = .data[['Cluster']]), size = label_size,
          show.legend = FALSE, color = label_color
        ) + scale_fill_manual(values = colors_use)
      } else if (isTRUE(x = label)) {
        geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes(label = .data[['Cluster']]), size = label_size, color = label_color,
          show.legend = FALSE
        )
      } else {
        p2 <- p2
      }

    }
  } else {
    if (!is.null(x = split_by)) {
      p2 <- lapply(1:length(x = list_of_splits), function(x){
        p2 <- ggplot(subset(tsne_df, tsne_df[[split_by]] %in% list_of_splits[x]),aes(x = .data[['tsne1']], y = .data[['tsne2']], color = .data[['Cluster']])) +
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

        if (isTRUE(x = label_box)) {
          geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes(label = .data[['Cluster']], fill = .data[['Cluster']]), size = label_size,
            show.legend = FALSE, color = label_color
          ) + scale_fill_manual(values = colors_use)
        } else if (isTRUE(x = label)) {
          geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes(label = .data[['Cluster']]), size = label_size, color = label_color,
            show.legend = FALSE
          )
        } else {
          p2 <- p2
        }
      })
    } else {
      p2 <- ggplot(tsne_df, aes(x = .data[['tsne1']], y = .data[['tsne2']], color = .data[['Cluster']])) +
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

      if (isTRUE(x = label_box)) {
        geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes(label = .data[['Cluster']], fill = .data[['Cluster']]), size = label_size,
          show.legend = FALSE, color = label_color
        ) + scale_fill_manual(values = colors_use)
      } else if (isTRUE(x = label)) {
        geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes(label = .data[['Cluster']]), size = label_size, color = label_color,
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
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
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
#' @importFrom rlang sym "!!"
#' @importFrom scattermore geom_scattermore
#' @importFrom utils packageVersion
#'
#' @references This function is encompasses part of the LIGER function plotByDatasetAndCluster.
#' However, this function is modified to just return cluster plots based on `Generate_Plotting_df_LIGER`.
#' \url{https://github.com/welch-lab/liger}. (Licence: GPL-3).
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
  raster.dpi = c(512, 512),
  ggplot_default_colors = FALSE,
  color_seed = 123
) {
  # temp liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    cli_abort(message = c("Liger functionality is currently restricted to rliger v1.0.1 or lower.",
                          "i" = "Functionality with rliger v2+ is currently in development."))
  }

  tsne_df <- Generate_Plotting_df_LIGER(object = liger_object, group_by = group_by, split_by = split_by, reorder.idents = reorder.idents, shuffle = shuffle, shuffle_seed = shuffle_seed)

  if (!is.null(x = split_by)) {
    list_of_splits <- unique(x = tsne_df[[split_by]])
  }

  # Get length of meta data feature
  if (!is.null(x = split_by) && !is.null(x = num_columns)) {
    split.by_length <- length(x = list_of_splits)

    # Calculate number of rows for selected number of columns
    num_rows <- ceiling(x = split.by_length/num_columns)

    # Check column and row compatibility
    if (num_columns > split.by_length) {
      cli_abort(message = c("The number of columns specified is greater than the number of meta data variables.",
                            "*" = "{.field {split_by}} only contains: {.field {split.by_length}} variables.",
                            "i" = "Please adjust {.code num_columns} to be less than or equal to: {.field {split.by_length}}.")
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

  group_by <- sym(x = group_by)

  if (isTRUE(x = raster)) {
    if (!is.null(x = split_by)) {
      p1 <- lapply(1:length(x = list_of_splits), function(x){
        ggplot(subset(tsne_df, tsne_df[[split_by]] %in% list_of_splits[x]), aes(x = .data[['tsne1']], y = .data[['tsne2']], color = !!group_by)) +
          theme_cowplot() +
          geom_scattermore(pointsize = pt_size, pixels = raster.dpi) +
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
      p1 <- ggplot(tsne_df, aes(x = .data[['tsne1']], y = .data[['tsne2']], color = !!group_by)) +
        theme_cowplot() +
        geom_scattermore(pointsize = pt_size, pixels = raster.dpi) +
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
        ggplot(subset(tsne_df, tsne_df[[split_by]] %in% list_of_splits[x]),aes(x = .data[['tsne1']], y = .data[['tsne2']], color = !!group_by)) +
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
      p1 <- ggplot(tsne_df, aes(x = .data[['tsne1']], y = .data[['tsne2']], color = !!group_by)) +
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


#' Customized version of plotFactors
#'
#' Modified and optimized version of `plotFactors` function from LIGER package.
#'
#' @param liger_object \code{liger} liger_object.  Need to perform clustering and factorization before calling this function
#' @param num_genes Number of genes to display for each factor (Default 8).
#' @param colors_use_factors colors to use for plotting factor loadings  By default datasets will be
#' plotted using "varibow" with shuffle = TRUE from both from \code{\link{DiscretePalette_scCustomize}}.
#' @param colors_use_dimreduc colors to use for plotting factor loadings on dimensionality reduction
#' coordinates (tSNE/UMAP).  Default is c('lemonchiffon', 'red'),
#' @param pt.size_factors Adjust point size for plotting in the factor plots.
#' @param pt.size_dimreduc Adjust point size for plotting in dimensionality reduction plots.
#' @param reduction Name of dimensionality reduction to use for plotting.
#' @param reduction_label `r lifecycle::badge("deprecated")` deprecated for newer style liger
#' objects.  Use `reduction` instead.
#' @param plot_legend logical, whether to plot the legend on factor loading plots, default is TRUE.
#' Helpful if number of datasets is large to avoid crowding the plot with legend.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param order logical. Whether to plot higher loading cells on top of cells with lower loading values in the
#' dimensionality reduction plots (Default = FALSE).
#' @param plot_dimreduc logical.  Whether to plot factor loadings on dimensionality reduction coordinates.  Default is TRUE.
#' @param save_plots logical.  Whether to save plots.  Default is TRUE
#' @param file_path directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix to append after sample name.
#' @param return_plots logical. Whether or not to return plots to the environment.  (Default is FALSE)
#' @param cells.highlight Names of specific cells to highlight in plot (black) (default NULL).
#' @param reorder_datasets `r lifecycle::badge("deprecated")` deprecated for newer style liger objects
#' @param ggplot_default_colors logical.  If `colors_use_factors = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "varibow" palette.
#' @param color_seed random seed for the palette shuffle if `colors_use_factors = NULL`.  Default = 123.
#'
#' @return A list of ggplot/patchwork objects and/or PDF file.
#'
#' @import cli
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom lifecycle deprecated
#' @importFrom patchwork wrap_plots
#' @importFrom scattermore geom_scattermore
#'
#' @noRd
#'
#' @concept liger_plotting
#'
#' @author Velina Kozareva (Original code for modified function), Sam Marsh (Added/modified functionality)
#' @references Based on `plotFactors` functionality from original LIGER package.
#'
#' @examples
#' \dontrun{
#' plotFactors_scCustom(liger_object = liger_obj, return_plots = FALSE, plot_dimreduc = TRUE,
#' raster = FALSE, save_plots = TRUE)
#' }
#'

plotFactors_liger2_scCustom <- function(
    liger_object,
    num_genes = 8,
    colors_use_factors = NULL,
    colors_use_dimreduc = c('lemonchiffon', 'red'),
    pt.size_factors = 1,
    pt.size_dimreduc = 1,
    reduction = "UMAP",
    reduction_label = deprecated(),
    plot_legend = TRUE,
    raster = TRUE,
    raster.dpi = c(512, 512),
    order = FALSE,
    plot_dimreduc = TRUE,
    save_plots = TRUE,
    file_path = NULL,
    file_name = NULL,
    return_plots = FALSE,
    cells.highlight = NULL,
    reorder_datasets = deprecated(),
    ggplot_default_colors = FALSE,
    color_seed = 123
) {
  # Check is slot is supplied
  if (lifecycle::is_present(reorder_datasets)) {
    lifecycle::deprecate_warn(when = "2.2.0",
                              what = "plotFactors_scCustom(reorder_datasets)",
                              details = c("i" = "The {.code reorder_datasets} parameter is deprecated for newer style Liger objects.",)
    )
  }

  # Check is slot is supplied
  if (lifecycle::is_present(reduction_label)) {
    lifecycle::deprecate_warn(when = "2.2.0",
                              what = "plotFactors_scCustom(reduction_label)",
                              details = c("v" = "The {.code reduction_label} parameter is deprecated for newer style Liger objects.",
                                          "i" = "Use {.code reduction} parameter instead")
    )
  }

  # if returning and saving
  if (isTRUE(x = save_plots)) {
    # Check file path is valid
    if (!is.null(x = file_path) && file_path != "") {
      if (!dir.exists(paths = file_path)) {
        cli_abort(message = "Provided {.code file_path}: {.val {file_path}} does not exist.")
      }
    }

    # Set file_path before path check if current dir specified as opposed to leaving set to NULL
    if (is.null(x = file_path)) {
      file_path <- ""
    }

    # Check if file name provided
    file_ext <- grep(x = file_name, pattern = ".pdf$", ignore.case = TRUE)
    if (length(x = file_ext) == 0) {
      file_name <- file_name
    } else {
      file_name <- gsub(pattern = ".pdf", replacement = "", x = file_name, ignore.case = TRUE)
    }

    if (is.null(x = file_name)) {
      cli_abort(message = c("No file name provided.",
                            "i" = "Please provide a file name using {.code file_name}.")
      )
    }
  }

  # Extract dataset number
  num_datasets <- length(x = liger_object@datasets)

  # Default Colors for Factor Plots
  if (is.null(x = colors_use_factors)) {
    if (isTRUE(x = ggplot_default_colors)) {
      colors_use_factors <- Hue_Pal(num_colors = num_datasets)
    } else {
      colors_use_factors <- DiscretePalette_scCustomize(num_colors = num_datasets, palette = "varibow", shuffle_pal = TRUE, seed = color_seed)
    }
  }

  # Check valid number of colors for tsne/UMAP
  if (length(x = colors_use_dimreduc) < 2) {
    cli_abort(message = c("Less than two values provided to {.code colors_use_dimreduc}.",
                          "i" = "Must provided either two colors to use for creating a gradient or a larger color gradient.")
    )
  }

  # Get Data and Plot Factors
  cli_inform(message = "{.field Generating plots}")
  k <- ncol(x = liger_object@H.norm)
  pb <- txtProgressBar(min = 0, max = k, style = 3)
  W <- liger_object@W
  rownames(x = W) <- rownames(x = liger_object@datasets[[1]]@scaleData)
  Hs_norm <- liger_object@H.norm
  dataset_names <- names(liger_object@datasets)
  H_raw_list <- lapply(1:num_datasets, function(x){
    H_raw <- t(liger_object@datasets[[x]]@H)
  })
  H_raw = do.call(rbind, H_raw_list)
  # Create accurate axis labels
  reduc_check <- LIGER_DimReduc(liger_object = liger_object, reduction = reduction, check_only = TRUE)
  x_axis_label <- paste0(reduction, "_1")
  y_axis_label <- paste0(reduction, "_2")
  plot_list = list()
  tsne_list = list()
  for (i in 1:k) {
    top_genes.W <- rownames(x = W)[order(W[, i], decreasing = T)[1:num_genes]]
    top_genes.W.string <- paste0(top_genes.W, collapse = ", ")
    factor_textstring <- paste0("Factor", i)
    plot_title1 <- paste(factor_textstring, "\n", top_genes.W.string, "\n")
    h_df = data.frame(x = 1:nrow(Hs_norm), h_norm = Hs_norm[, i],
                      h_raw = H_raw[, i], dataset = liger_object@cellMeta$dataset,
                      highlight = FALSE)
    if (isTRUE(x = raster)) {
      top <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_raw"]], col = .data[["dataset"]])) +
        geom_scattermore(pointsize = pt.size_factors, pixels = raster.dpi) +
        labs(x = 'Cell', y = 'Raw H Score') +
        ggtitle(plot_title1) +
        theme(legend.position = 'none') +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        top <- top + NoLegend()
      }

      bottom <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_norm"]], col = .data[["dataset"]])) +
        geom_scattermore(pointsize = pt.size_factors, pixels = raster.dpi) +
        labs(x = 'Cell', y = 'H_norm Score') +
        theme(legend.position = 'top',
              legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2))) +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        bottom <- bottom + NoLegend()
      }

    } else {
      top <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_raw"]], col = .data[["dataset"]])) +
        geom_point(size = pt.size_factors) +
        labs(x = 'Cell', y = 'Raw H Score') +
        ggtitle(plot_title1) +
        theme(legend.position = 'none') +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        top <- top + NoLegend()
      }

      bottom <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_norm"]], col = .data[["dataset"]])) +
        geom_point(size = pt.size_factors) +
        labs(x = 'Cell', y = 'H_norm Score') +
        theme(legend.position = 'top',
              legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2))) +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        bottom <- bottom + NoLegend()
      }

    }

    if (!is.null(cells.highlight)) {
      h_df[cells.highlight, 'highlight'] = TRUE
      if (isTRUE(x = raster)) {
        top <- top + geom_scattermore(data = subset(h_df, .data[["highlight"]] == TRUE),
                                      aes(.data[["x"]], .data[["h_raw"]]),
                                      col = "black",
                                      pointsize = pt.size_factors,
                                      pixels = raster.dpi)
        bottom <- bottom + geom_scattermore(data = subset(h_df, .data[["highlight"]] == TRUE),
                                            aes(.data[["x"]], .data[["h_norm"]]),
                                            col = "black",
                                            pointsize = pt.size_factors,
                                            pixels = raster.dpi)
      } else {
        top <- top + geom_point(data = subset(h_df, .data[["highlight"]] == TRUE),
                                aes(.data[["x"]], .data[["h_raw"]]),
                                col = "black",
                                size = pt.size_factors)
        bottom <- bottom + geom_point(data = subset(h_df, .data[["highlight"]] == TRUE),
                                      aes(.data[["x"]], .data[["h_norm"]]),
                                      col = "black",
                                      size = pt.size_factors)
      }
    }
    full <- wrap_plots(top, bottom, ncol = 1)
    plot_list[[i]] = full

    # plot tSNE/UMAP
    if (isTRUE(x = plot_dimreduc)) {
      tsne_df <- data.frame(Hs_norm[, i], LIGER_DimReduc(liger_object = liger_object, reduction = reduction))
      factorlab <- paste0("Factor", i)
      colnames(x = tsne_df) <- c(factorlab, x_axis_label, y_axis_label)

      if (isTRUE(x = order)) {
        tsne_df <- tsne_df[order(tsne_df[,1], decreasing = FALSE),]
      }

      if (isTRUE(x = raster)) {
        p1 <- ggplot(tsne_df, aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = .data[[factorlab]])) +
          geom_scattermore(pointsize = pt.size_dimreduc, pixels = raster.dpi) +
          ggtitle(label = paste('Factor', i)) +
          theme(legend.position = 'none') +
          xlab(x_axis_label) +
          ylab(y_axis_label) +
          if (length(x = colors_use_dimreduc) == 2) {
            scale_color_gradient(low = colors_use_dimreduc[1], high = colors_use_dimreduc[2])
          } else {
            scale_color_gradientn(colours = colors_use_dimreduc)
          }
      } else {
        p1 <- ggplot(tsne_df, aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = .data[[factorlab]])) +
          geom_point(size = pt.size_dimreduc) +
          ggtitle(label = paste('Factor', i)) +
          theme(legend.position = 'none') +
          xlab(x_axis_label) +
          ylab(y_axis_label) +
          if (length(x = colors_use_dimreduc) == 2) {
            scale_color_gradient(low = colors_use_dimreduc[1], high = colors_use_dimreduc[2])
          } else {
            scale_color_gradientn(colours = colors_use_dimreduc)
          }
      }

      tsne_list[[i]] = p1
    }
    setTxtProgressBar(pb, i)
  }

  # save plots
  if (isTRUE(x = save_plots)) {
    cli_inform(message = "{.field Saving plots to file}")
    pdf(paste(file_path, file_name, ".pdf", sep=""))
    pb <- txtProgressBar(min = 0, max = length(x = 1:k), style = 3, file = stderr())
    for (i in 1:k) {
      if (isTRUE(x = plot_dimreduc)) {
        print(plot_list[[i]])
        print(tsne_list[[i]])
        setTxtProgressBar(pb = pb, value = i)
      } else {
        print(plot_list[[i]])
        setTxtProgressBar(pb = pb, value = i)
      }
    }
    close(con = pb)
    dev.off()
  }

  # return plots
  if (isTRUE(x = return_plots)) {
    return(list(factor_plots = plot_list,
                dimreduc_plots = tsne_list))
  }
}


#' Customized version of plotFactors
#'
#' Modified and optimized version of `plotFactors` function from LIGER package.
#'
#' @param liger_object \code{liger} liger_object.  Need to perform clustering and factorization before calling this function
#' @param num_genes Number of genes to display for each factor (Default 8).
#' @param colors_use_factors colors to use for plotting factor loadings  By default datasets will be
#' plotted using "varibow" with shuffle = TRUE from both from \code{\link{DiscretePalette_scCustomize}}.
#' @param colors_use_dimreduc colors to use for plotting factor loadings on dimensionality reduction
#' coordinates (tSNE/UMAP).  Default is c('lemonchiffon', 'red'),
#' @param pt.size_factors Adjust point size for plotting in the factor plots.
#' @param pt.size_dimreduc Adjust point size for plotting in dimensionality reduction plots.
#' @param reduction_label What to label the x and y axes of resulting plots.  LIGER does not store name of
#' technique and therefore needs to be set manually.  Default is "UMAP".
#' @param plot_legend logical, whether to plot the legend on factor loading plots, default is TRUE.
#' Helpful if number of datasets is large to avoid crowding the plot with legend.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param order logical. Whether to plot higher loading cells on top of cells with lower loading values in the
#' dimensionality reduction plots (Default = FALSE).
#' @param plot_dimreduc logical.  Whether to plot factor loadings on dimensionality reduction coordinates.  Default is TRUE.
#' @param save_plots logical.  Whether to save plots.  Default is TRUE
#' @param file_path directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix to append after sample name.
#' @param return_plots logical. Whether or not to return plots to the environment.  (Default is FALSE)
#' @param cells.highlight Names of specific cells to highlight in plot (black) (default NULL).
#' @param reorder_datasets New order to plot datasets in for the factor plots if different from current
#' factor level order in cell.data slot.
#' @param ggplot_default_colors logical.  If `colors_use_factors = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "varibow" palette.
#' @param color_seed random seed for the palette shuffle if `colors_use_factors = NULL`.  Default = 123.
#'
#' @return A list of ggplot/patchwork objects and/or PDF file.
#'
#' @import cli
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom patchwork wrap_plots
#' @importFrom scattermore geom_scattermore
#'
#' @noRd
#'
#' @concept liger_plotting
#'
#' @author Velina Kozareva (Original code for modified function), Sam Marsh (Added/modified functionality)
#' @references Based on `plotFactors` functionality from original LIGER package.
#'
#' @examples
#' \dontrun{
#' plotFactors_scCustom(liger_object = liger_obj, return_plots = FALSE, plot_dimreduc = TRUE,
#' raster = FALSE, save_plots = TRUE)
#' }
#'

plotFactors_liger_scCustom <- function(
    liger_object,
    num_genes = 8,
    colors_use_factors = NULL,
    colors_use_dimreduc = c('lemonchiffon', 'red'),
    pt.size_factors = 1,
    pt.size_dimreduc = 1,
    reduction_label = "UMAP",
    plot_legend = TRUE,
    raster = TRUE,
    raster.dpi = c(512, 512),
    order = FALSE,
    plot_dimreduc = TRUE,
    save_plots = TRUE,
    file_path = NULL,
    file_name = NULL,
    return_plots = FALSE,
    cells.highlight = NULL,
    reorder_datasets = NULL,
    ggplot_default_colors = FALSE,
    color_seed = 123
) {
  # if returning and saving
  if (isTRUE(x = save_plots)) {

    # Check file path is valid
    if (!is.null(x = file_path) && file_path != "") {
      if (!dir.exists(paths = file_path)) {
        cli_abort(message = "Provided {.code file_path}: {.val {file_path}} does not exist.")
      }
    }

    # Set file_path before path check if current dir specified as opposed to leaving set to NULL
    if (is.null(x = file_path)) {
      file_path <- ""
    }

    # Check if file name provided
    file_ext <- grep(x = file_name, pattern = ".pdf$", ignore.case = TRUE)
    if (length(x = file_ext) == 0) {
      file_name <- file_name
    } else {
      file_name <- gsub(pattern = ".pdf", replacement = "", x = file_name, ignore.case = TRUE)
    }

    if (is.null(x = file_name)) {
      cli_abort(message = c("No file name provided.",
                            "i" = "Please provide a file name using {.code file_name}.")
      )
    }
  }

  if (!is.null(x = reorder_datasets)) {
    # Check new order contains same dataset names and number of datasets
    if (length(x = levels(x = liger_object@cell.data$dataset)) != length(x = reorder_datasets)) {
      cli_abort(message = c("Error reordering datasets (number mismatch).",
                            "i" = "The number of datasets provided to {.code reorder_datasets} ({.field {length(x = reorder_datasets)}}) does not match number of datasets in LIGER object ({.field {length(x = levels(x = levels(liger_object@cell.data$dataset)))}}).")
      )
    } else {
      if (!all(levels(x = liger_object@cell.data$dataset) %in% reorder_datasets)) {
        cli_abort(message = c("Error reordering datasets (name mismatch).",
                              "*" = "Dataset names provided to {.code reorder_datasets} do not match names of datasets in LIGER object.",
                              "i" = "Please check spelling.")
        )
      } else {
        liger_object@cell.data$dataset <- factor(x = liger_object@cell.data$dataset, levels = reorder_datasets)
      }
    }
  }

  # Create accurate axis labels
  x_axis_label <- paste0(reduction_label, "_1")
  y_axis_label <- paste0(reduction_label, "_2")

  # Extract dataset number
  num_datasets <- length(x = liger_object@scale.data)

  # Default Colors for Factor Plots
  if (is.null(x = colors_use_factors)) {
    if (isTRUE(x = ggplot_default_colors)) {
      colors_use_factors <- Hue_Pal(num_colors = num_datasets)
    } else {
      colors_use_factors <- DiscretePalette_scCustomize(num_colors = num_datasets, palette = "varibow", shuffle_pal = TRUE, seed = color_seed)
    }
  }

  # Check valid number of colors for tsne/UMAP
  if (length(x = colors_use_dimreduc) < 2) {
    cli_abort(message = c("Less than two values provided to {.code colors_use_dimreduc}.",
                          "i" = "Must provided either two colors to use for creating a gradient or a larger color gradient.")
    )
  }

  # Add one time dim label warning
  if (getOption(x = 'scCustomize_warn_LIGER_dim_labels_plotFactors', default = TRUE)) {
    cli_inform(message = c("",
                           "NOTE: {.field plotFactors_scCustom} uses the {.code reduction_label} parameter to set axis labels",
                           "on the dimensionality reduction plots.",
                           "By default this is set to {.val UMAP}.",
                           "Please take note of this parameter as LIGER objects do not store the name",
                           "of reduction technique used and therefore this needs to be set manually.",
                           "",
                           "-----This message will be shown once per session.-----"))
    options(scCustomize_warn_LIGER_dim_labels_plotFactors = FALSE)
  }

  # Get Data and Plot Factors
  cli_inform(message = "{.field Generating plots}")
  k <- ncol(x = liger_object@H.norm)
  pb <- txtProgressBar(min = 0, max = k, style = 3)
  W <- t(x = liger_object@W)
  rownames(x = W) <- colnames(x = liger_object@scale.data[[1]])
  Hs_norm <- liger_object@H.norm
  H_raw = do.call(rbind, liger_object@H)
  plot_list = list()
  tsne_list = list()
  for (i in 1:k) {
    top_genes.W <- rownames(x = W)[order(W[, i], decreasing = T)[1:num_genes]]
    top_genes.W.string <- paste0(top_genes.W, collapse = ", ")
    factor_textstring <- paste0("Factor", i)
    plot_title1 <- paste(factor_textstring, "\n", top_genes.W.string, "\n")
    h_df = data.frame(x = 1:nrow(Hs_norm), h_norm = Hs_norm[, i],
                      h_raw = H_raw[, i], dataset = liger_object@cell.data$dataset,
                      highlight = FALSE)
    if (isTRUE(x = raster)) {
      top <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_raw"]], col = .data[["dataset"]])) +
        geom_scattermore(pointsize = pt.size_factors, pixels = raster.dpi) +
        labs(x = 'Cell', y = 'Raw H Score') +
        ggtitle(plot_title1) +
        theme(legend.position = 'none') +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        top <- top + NoLegend()
      }

      bottom <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_norm"]], col = .data[["dataset"]])) +
        geom_scattermore(pointsize = pt.size_factors, pixels = raster.dpi) +
        labs(x = 'Cell', y = 'H_norm Score') +
        theme(legend.position = 'top',
              legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2))) +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        bottom <- bottom + NoLegend()
      }

    } else {
      top <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_raw"]], col = .data[["dataset"]])) +
        geom_point(size = pt.size_factors) +
        labs(x = 'Cell', y = 'Raw H Score') +
        ggtitle(plot_title1) +
        theme(legend.position = 'none') +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        top <- top + NoLegend()
      }

      bottom <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_norm"]], col = .data[["dataset"]])) +
        geom_point(size = pt.size_factors) +
        labs(x = 'Cell', y = 'H_norm Score') +
        theme(legend.position = 'top',
              legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2))) +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        bottom <- bottom + NoLegend()
      }

    }

    if (!is.null(cells.highlight)) {
      h_df[cells.highlight, 'highlight'] = TRUE
      if (isTRUE(x = raster)) {
        top <- top + geom_scattermore(data = subset(h_df, .data[["highlight"]] == TRUE),
                                      aes(.data[["x"]], .data[["h_raw"]]),
                                      col = "black",
                                      pointsize = pt.size_factors,
                                      pixels = raster.dpi)
        bottom <- bottom + geom_scattermore(data = subset(h_df, .data[["highlight"]] == TRUE),
                                            aes(.data[["x"]], .data[["h_norm"]]),
                                            col = "black",
                                            pointsize = pt.size_factors,
                                            pixels = raster.dpi)
      } else {
        top <- top + geom_point(data = subset(h_df, .data[["highlight"]] == TRUE),
                                aes(.data[["x"]], .data[["h_raw"]]),
                                col = "black",
                                size = pt.size_factors)
        bottom <- bottom + geom_point(data = subset(h_df, .data[["highlight"]] == TRUE),
                                      aes(.data[["x"]], .data[["h_norm"]]),
                                      col = "black",
                                      size = pt.size_factors)
      }
    }
    full <- wrap_plots(top, bottom, ncol = 1)
    plot_list[[i]] = full

    # plot tSNE/UMAP
    if (isTRUE(x = plot_dimreduc)) {
      tsne_df <- data.frame(Hs_norm[, i], liger_object@tsne.coords)
      factorlab <- paste0("Factor", i)
      colnames(x = tsne_df) <- c(factorlab, x_axis_label, y_axis_label)

      if (isTRUE(x = order)) {
        tsne_df <- tsne_df[order(tsne_df[,1], decreasing = FALSE),]
      }

      if (isTRUE(x = raster)) {
        p1 <- ggplot(tsne_df, aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = .data[[factorlab]])) +
          geom_scattermore(pointsize = pt.size_dimreduc, pixels = raster.dpi) +
          ggtitle(label = paste('Factor', i)) +
          theme(legend.position = 'none') +
          xlab(x_axis_label) +
          ylab(y_axis_label) +
          if (length(x = colors_use_dimreduc) == 2) {
            scale_color_gradient(low = colors_use_dimreduc[1], high = colors_use_dimreduc[2])
          } else {
            scale_color_gradientn(colours = colors_use_dimreduc)
          }
      } else {
        p1 <- ggplot(tsne_df, aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = .data[[factorlab]])) +
          geom_point(size = pt.size_dimreduc) +
          ggtitle(label = paste('Factor', i)) +
          theme(legend.position = 'none') +
          xlab(x_axis_label) +
          ylab(y_axis_label) +
          if (length(x = colors_use_dimreduc) == 2) {
            scale_color_gradient(low = colors_use_dimreduc[1], high = colors_use_dimreduc[2])
          } else {
            scale_color_gradientn(colours = colors_use_dimreduc)
          }
      }

      tsne_list[[i]] = p1
    }
    setTxtProgressBar(pb, i)
  }

  # save plots
  if (isTRUE(x = save_plots)) {
    cli_inform(message = "{.field Saving plots to file}")
    pdf(paste(file_path, file_name, ".pdf", sep=""))
    pb <- txtProgressBar(min = 0, max = length(x = 1:k), style = 3, file = stderr())
    for (i in 1:k) {
      if (isTRUE(x = plot_dimreduc)) {
        print(plot_list[[i]])
        print(tsne_list[[i]])
        setTxtProgressBar(pb = pb, value = i)
      } else {
        print(plot_list[[i]])
        setTxtProgressBar(pb = pb, value = i)
      }
    }
    close(con = pb)
    dev.off()
  }

  # return plots
  if (isTRUE(x = return_plots)) {
    return(list(factor_plots = plot_list,
                dimreduc_plots = tsne_list))
  }
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### ANALYSIS UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
#' @importFrom utils packageVersion
#'
#' @references Matching function parameter text descriptions are taken from `rliger::selectGenes`
#' which is called by this function after creating new temporary object/dataset.
#' \url{https://github.com/welch-lab/liger}. (License: GPL-3).
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' liger_obj <- Variable_Features_ALL_LIGER(liger_object = liger_obj, num_genes = 2000)
#' }
#'

Variable_Features_ALL_LIGER <- function(
  liger_object,
  num_genes = NULL,
  var.thresh = 0.3,
  alpha.thresh = 0.99,
  tol = 0.0001,
  do.plot = FALSE,
  pt.size = 0.3,
  chunk=1000
) {
  # temp liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    cli_abort(message = c("Liger functionality is currently restricted to rliger v1.0.1 or lower.",
                          "i" = "Functionality with rliger v2+ is currently in development."))
  }

  Is_LIGER(liger_object = liger_object)

  raw_data <- liger_object@raw.data

  cli_inform(message = "Creating temporary object with combined data.")

  temp_liger <- rliger::createLiger(raw.data = list("dataset" = Merge_Sparse_Data_All(raw_data)), remove.missing = FALSE)

  rm(raw_data)
  gc()

  cli_inform(message = "Normalizing and identifying variable features.")

  temp_liger <- rliger::normalize(object = temp_liger)
  temp_liger <- rliger::selectGenes(object = temp_liger, var.thresh = var.thresh, do.plot = do.plot, num.genes = num_genes, tol = tol, alpha.thresh = alpha.thresh, cex.use = pt.size, chunk = chunk)
  var_genes <- temp_liger@var.genes

  rm(temp_liger)
  gc()

  liger_object@var.genes <- var_genes
  return(liger_object)
}
