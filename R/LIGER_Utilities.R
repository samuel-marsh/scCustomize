#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### EXTENDED SEURAT GENERICS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Extract Features from LIGER Object
#'
#' Extract all unique features from LIGER object
#'
#' @param x LIGER object name.
#' @param by_dataset logical, whether to return list with vector of features for each dataset in
#' LIGER object or to return single vector of unique features across all datasets in object
#' (default is FALSE; return vector of unique features)
#' @param ... Arguments passed to other methods
#'
#' @method Features liger
#' @return vector or list depending on `by_dataset` parameter
#'
#' @importFrom utils packageVersion
#'
#' @export
#' @rdname Features
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # return single vector of all unique features
#' all_features <- Features(x = object, by_dataset = FALSE)
#'
#' # return list of vectors containing features from each individual dataset in object
#' dataset_features <- Features(x = object, by_dataset = TRUE)
#' }
#'

Features.liger <- function(
    x,
    by_dataset = FALSE,
    ...
) {
  # liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    # Extract features
    features_by_dataset <- lapply(1:length(x = x@datasets), function(j) {
      rownames(x = x@datasets[[j]]@featureMeta)
    })
  } else {
    # Extract features
    features_by_dataset <- lapply(1:length(x = x@raw.data), function(j) {
      rownames(x = x@raw.data[[j]])
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


#' Extract Cells from LIGER Object
#'
#' Extract all cell barcodes from LIGER object
#'
#' @param x LIGER object name.
#' @param by_dataset logical, whether to return list with vector of cell barcodes for each
#' dataset in LIGER object or to return single vector of cell barcodes across all
#' datasets in object (default is FALSE; return vector of cells).
#' @param ... Arguments passed to other methods
#'
#' @method Cells liger
#' @return vector or list depending on `by_dataset` parameter
#'
#' @importFrom utils packageVersion
#'
#' @export
#' @rdname Cells
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # return single vector of all cells
#' all_features <- Cells(x = object, by_dataset = FALSE)
#'
#' # return list of vectors containing cells from each individual dataset in object
#' dataset_features <- Cells(x = object, by_dataset = TRUE)
#' }
#'

Cells.liger <- function(
    x,
    by_dataset = FALSE,
    ...
) {
  # liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    # Extract features
    cells_by_dataset <- lapply(1:length(x = x@datasets), function(j) {
      colnames(x = x@datasets[[j]])
    })
    names(cells_by_dataset) <- names(x@datasets)
  } else {
    # Extract features
    cells_by_dataset <- lapply(1:length(x = x@raw.data), function(j) {
      colnames(x = x@raw.data[[j]])
    })
  }

  # Return features
  if (isFALSE(x = by_dataset)) {
    cells <- unlist(x = cells_by_dataset)
    return(cells)
  } else {
    return(cells_by_dataset)
  }
}


#' Extract Cells for particular identity
#'
#' Extract all cell barcodes for a specific identity
#'
#' @param object LIGER object name.
#' @param idents identities to extract cell barcodes.
#' @param ident_col name of meta data column to use when subsetting cells by identity values.
#' Default is NULL, which will use the objects default clustering as the `ident_col`.
#' @param by_dataset logical, whether to return vector with cell barcodes for all `idents` in or
#' to return list (1 entry per dataset with vector of cells) (default is FALSE; return vector).
#' @param invert logical, invert the selection of cells (default is FALSE).
#' @param ... Arguments passed to other methods
#'
#' @method WhichCells liger
#' @return vector or list depending on `by_dataset` parameter
#'
#' @concept liger_object_util
#'
#' @import cli
#' @import Seurat
#' @importFrom dplyr all_of select filter
#' @importFrom magrittr "%>%"
#'
#' @export
#' @rdname WhichCells
#'
#' @examples
#' \dontrun{
#' # Extract cells from ident =1 in current default clustering
#' ident1_cells <- WhichCells(object = liger_object, idents = 1)
#'
#' # Extract all cells from "stim" treatment from object
#' stim_cells <- WhichCells(object = liger_object, idents = "stim", ident_col = "Treatment")
#' }
#'

WhichCells.liger <- function(
    object,
    idents = NULL,
    ident_col = NULL,
    by_dataset = FALSE,
    invert = FALSE,
    ...
) {
  # Check new liger object
  if (packageVersion(pkg = 'rliger') < "2.0.0") {
    cli_abort(message = "This function is only for objects with rliger >= v2.0.0")
  }

  # Get cells data.frame
  ident_col <- ident_col %||% LIGER_Default_Cluster_Name(liger_object = object)

  cell_df <- Fetch_Meta(object = object) %>%
    select(all_of(c(ident_col, "dataset")))

  # possible idents
  if (inherits(x = cell_df[[ident_col]], what = "factor")) {
    ident_levels <- levels(x = cell_df[[ident_col]])
  } else {
    ident_levels <- unique(x = cell_df[[ident_col]])
  }

  # check idents valid
  valid_idents <- intersect(x = idents, y = ident_levels)

  if (length(x = valid_idents) == 0) {
    cli_abort(message = "None of the provided {.code idents} were found in object.")
  }
  if (length(x = valid_idents) != length(x = idents)) {
    missing_idents <- setdiff(x = idents, y = valid_idents)
    cli_warn(message = c("The following {.code idents} were not found and therefore ignored:",
                         "i" = "{.field {missing_idents}}"))
  }

  # get cells
  if (isFALSE(x = by_dataset)) {
    cells <- cell_df %>%
      filter(.data[[ident_col]] %in% valid_idents) %>%
      rownames()
    if (isTRUE(x = invert)) {
      cells <- setdiff(x = Cells(x = object, by_dataset = FALSE), y = cells)
    }
  } else {
    dataset_names <- names(x = rliger::datasets(x = object))
    cells <- lapply(dataset_names, function(x) {
      sample_cells <- cell_df %>%
        filter(.data[["dataset"]] == x & .data[[ident_col]] %in% valid_idents) %>%
        rownames()
    })
    if (isTRUE(x = invert)) {
      all_cells <- Cells(x = object, by_dataset = TRUE)
      cells <- lapply(1:length(x = cells), function(x) {
        cells_inverted <- setdiff(x = all_cells[[x]], y = cells[[x]])
      })
    }

    names(x = cells) <- dataset_names
  }

  # return cells
  return(cells)
}


#' Extract matrix of embeddings
#'
#' Extract matrix containing iNMF or dimensionality reduction embeddings.
#'
#' @param object LIGER object name.
#' @param reduction name of dimensionality reduction to pull
#' @param iNMF logical, whether to extract iNMF h.norm matrix instead of dimensionality reduction embeddings.
#' @param check_only logical, return `TRUE` if valid reduction is present.
#' @param ... Arguments passed to other methods
#'
#' @method Embeddings liger
#' @return matrix
#'
#' @concept liger_object_util
#'
#' @import cli
#' @import Seurat
#'
#' @export
#' @rdname Embeddings
#'
#' @examples
#' \dontrun{
#' # Extract embedding matrix for current dimensionality reduction
#' UMAP_coord <- Embeddings(object = liger_object)
#'
#' # Extract iNMF h.norm matrix
#' iNMF_mat <- Embeddings(object = liger_object, reduction = "iNMF")
#' }
#'

Embeddings.liger <- function(
    object,
    reduction = NULL,
    iNMF = FALSE,
    check_only = FALSE,
    ...
) {
  # Check new liger object
  if (packageVersion(pkg = 'rliger') < "2.0.0") {
    cli_abort(message = "This function is only for objects with rliger >= v2.0.0")
  }

  # return iNMF h.norm
  if (isTRUE(x = iNMF)) {
    embeddings <- object@h.norm
    return(embeddings)
  }

  # set reduction if not supplied
  reduction <- reduction %||% Default_DimReduc_LIGER(liger_object = object)

  # check reduction in cellMeta
  if (reduction %in% names(x = rliger::dimReds(x = object))) {
    if (isTRUE(x = check_only)) {
      return(TRUE)
    }
    # get coords
    embeddings <- rliger::dimReds(x = object)[[reduction]]
  } else {
    cli_abort("The reduction {.field {reduction}} is not present in dimReds slot.")
  }

  # return embeddings
  return(embeddings)
}


#' Extract or set default identities from object
#'
#' Extract default identities from object in factor form.
#'
#' @param object LIGER object name.
#' @param ... Arguments passed to other methods
#'
#' @method Idents liger
#' @return factor
#'
#' @concept liger_object_util
#'
#' @import cli
#' @import Seurat
#' @importFrom dplyr pull
#' @importFrom magrittr "%>%"
#'
#' @export
#' @rdname Idents
#'
#' @examples
#' \dontrun{
#' # Extract idents
#' object_idents <- Idents(object = liger_object)
#' }
#'

Idents.liger <- function(
    object,
    ...
) {
  # Check default cluster present
  if (is.null(x = rliger::defaultCluster(x = object))) {
    cli_abort(message = "No default cell identity/cluster present in object.")
  }

  # get current default ident name
  identity_name <- LIGER_Default_Cluster_Name(liger_object = object)

  # pull active ident column
  active_idents <- Fetch_Meta(object = object) %>%
    pull(.data[[identity_name]])

  # return active idents
  return(active_idents)
}


#' Set default identities of object
#'
#' @param object LIGER object name.
#' @param value name of column in cellMeta slot to set as new default cluster/ident
#' @param ... Arguments passed to other methods
#'
#' @method Idents<- liger
#' @return object
#'
#' @note Use of Idents<- is only for setting new default ident/cluster from column already present in cellMeta.
#' To add new column with new cluster values to cellMeta and set as default see \code{\link{Rename_Clusters}}.
#'
#' @concept liger_object_util
#'
#' @import cli
#' @import Seurat
#'
#' @export
#' @rdname Idents
#'
#' @examples
#' \dontrun{
#' # Set idents
#' Idents(object = liger_object) <- "new_annotation"
#' }
#'

"Idents<-.liger" <- function(
    object,
    ...,
    value
) {
  # Check new ident value is present in cellMeta
  new_ident_name <- Meta_Present(object = object, meta_col_names = value, print_msg = FALSE, omit_warn = FALSE)[[1]]

  if (length(x = new_ident_name) == 0) {
    cli_abort(message = c("The provided {.code value} ({.field {value}}) is not present in obect cellMeta slot.",
                          "i" = "Provide different value or use {.field Rename_Clusters} to add vector of new idents to cellMeta."))
  }

  # change defaults
  rliger::defaultCluster(x = object) <- new_ident_name

  # return object
  return(object)
}


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
    object_meta <- rliger::cellMeta(x = object, as.data.frame = TRUE)
  } else {
    object_meta <- object_meta <- slot(object = object, name = "cell.data")
  }

  # return meta
  return(object_meta)
}


#' Extract Cells by identity
#'
#' Extract all cell barcodes by identity from LIGER object
#'
#' @param liger_object LIGER object name.
#' @param group.by name of meta data column to use, default is current default clustering.
#' @param by_dataset logical, whether to return list with entries for cell barcodes for each
#' identity in `group.by`
#' or to return list of lists (1 entry per dataset and each ident within the dataset)
#' (default is FALSE; return list)
#'
#' @return list or list of lists depending on `by_dataset` parameter
#'
#' @import cli
#' @importFrom dplyr filter select all_of
#' @importFrom magrittr "%>%"
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # return single vector of all cells
#' cells_by_idents <- LIGER_Cells_by_Identities(liger_object = object, by_dataset = FALSE)
#'
#' # return list of vectors containing cells from each individual dataset in object
#' cells_by_idents_by_dataset <- LIGER_Cells_by_Identities(liger_object = object, by_dataset = TRUE)
#' }
#'

LIGER_Cells_by_Identities <- function(
    liger_object,
    group.by = NULL,
    by_dataset = FALSE
) {
  # check liger
  Is_LIGER(liger_object = liger_object)

  # Check new liger object
  if (packageVersion(pkg = 'rliger') < "2.0.0") {
    cli_abort(message = "This function is only for objects created with rliger >= v2.0.0")
  }

  # check group.by is valid
  if (!is.null(x = group.by)) {
    Meta_Present(object = liger_object, meta_col_names = group.by, print_msg = FALSE)
  }

  # set group.by if not set
  group.by <- group.by %||% LIGER_Default_Cluster_Name(liger_object = liger_object)

  # Check cluster df
  cell_df <- Fetch_Meta(object = liger_object) %>%
    select(all_of(c(group.by, "dataset")))

  if (inherits(x = cell_df[[group.by]], what = "factor")) {
    ident_levels <- levels(x = cell_df[[group.by]])
  } else {
    ident_levels <- unique(x = cell_df[[group.by]])
  }

  # Get cells for object overall
  if (isFALSE(x = by_dataset)) {
    cells_list <- lapply(ident_levels, function(x) {
      cells <- cell_df %>%
        filter(.data[[group.by]] == x) %>%
        rownames()
    })

    names(cells_list) <- ident_levels
  } else {
    # Get cells by cluster by dataset
    dataset_names <- names(x = rliger::datasets(x = liger_object))
    cells_list <- lapply(1:length(x = dataset_names), function(x) {
      sample_cells_df <- cell_df %>%
        filter(.data[["dataset"]] == dataset_names[x])

      sample_cells <- lapply(ident_levels, function(y) {
        sample_cells_df %>%
          filter(.data[[group.by]] == y) %>%
          rownames()
      })
      names(sample_cells) <- ident_levels

      return(sample_cells)
    })
    names(cells_list) <- dataset_names
  }

  return(cells_list)
}


#' @param new_idents vector of new cluster names.  Must be equal to the length of current default identity
#' of Object.  Will accept named vector (with old idents as names) or will name the new_idents vector internally.
#' @param meta_col_name `r lifecycle::badge("soft-deprecated")`. See `old_ident_name`.
#' @param old_ident_name optional, name to use for storing current object idents in object meta data slot.
#' @param new_ident_name optional, name to use for storing new object idents in object meta data slot.
#' @param overwrite logical, whether to overwrite columns in object meta data slot. if they have same
#' names as `old_ident_name` and/or `new_ident_name`.
#'
#' @method Rename_Clusters liger
#'
#' @import cli
#' @importFrom dplyr right_join
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @rdname Rename_Clusters
#' @export
#'
#' @concept marker_annotation_util
#'
#' @examples
#' \dontrun{
#' # Liger version
#' obj <- Rename_Clusters(object = obj_name, new_idents = new_idents_vec,
#' old_ident_name = "LIGER_Idents_Round01", new_ident_name = "LIGER_Idents_Round02")
#' }
#'

Rename_Clusters.liger <- function(
    object,
    new_idents,
    old_ident_name = NULL,
    new_ident_name = NULL,
    overwrite = FALSE,
    ...
) {
  # Check new liger object
  if (packageVersion(pkg = 'rliger') < "2.0.0") {
    cli_abort(message = "This function is only for objects with rliger >= v2.0.0")
  }

  # check old ident name
  if (!is.null(x = old_ident_name)) {
    if (old_ident_name %in% colnames(x = object@cellMeta)) {
      if (isFALSE(x = overwrite)) {
        cli_abort(message = c("The {.code old_ident_name}: {.field {old_ident_name}} is already a column in meta data",
                              "i" = "To overwrite current meta data column set {.code overwrite = TRUE}."))
      } else {
        cli_inform(message = "Overwriting old meta data column: {.field {old_ident_name}} as {.code overwrite = TRUE}")

      }
    } else {
      object@cellMeta[[old_ident_name]] <- rliger::defaultCluster(x = object)
    }
  }

  # check new ident name
  if (!is.null(x = new_ident_name) && new_ident_name %in% colnames(x = object@cellMeta)) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("The {.code new_ident_name}: {.field {new_ident_name}} is already a column in meta data",
                            "i" = "To overwrite current meta data column set {.code overwrite = TRUE}."))
    } else {
      cli_inform(message = "Overwriting new meta data column: {.field {new_ident_name}} as {.code overwrite = TRUE}")
    }
  }

  # Check equivalent lengths
  if (length(x = new_idents) != length(x = levels(x = rliger::defaultCluster(x = object)))) {
    cli_abort(message = c("Length of {.code new_idents} must be equal to the number of clusters in Liger Object.",
                          "i" = "{.code new_idents} length: {.field {length(x = new_idents)}} object 'defaultCluster' length: {.field {length(x = levels(x = rliger::defaultCluster(x = object)))}}.")
    )
  }

  # Name the new idents vector
  if (is.null(x = names(x = new_idents))) {
    names(x = new_idents) <- levels(x = rliger::defaultCluster(x = object))
  }

  # If named check that names are right length
  if (!is.null(x = names(x = new_idents)) && length(x = unique(x = names(x = new_idents))) != length(x = levels(x = rliger::defaultCluster(x = object)))) {
    cli_abort(message = c("The number of unique names for {.code new idents} is not equal to number of clusters.",
                          "i" = "names(new_idents) length: {.field {length(x = unique(x = names(x = new_idents)))} object 'defaultCluster' length: {length(x = levels(x = defaultCluster(x = object)))}}.")
    )
  }

  # Add new idents
  ident_df <- data.frame(rliger::defaultCluster(x = object))
  colnames(x = ident_df) <- "current_idents"
  ident_df <- ident_df %>%
    rownames_to_column("barcodes")

  new_idents_df <- data.frame("current_idents" = names(x = new_idents),
                              "new_idents" = new_idents)

  new_idents_meta <- suppressMessages(right_join(x = ident_df, y = new_idents_df)) %>%
    column_to_rownames("barcodes")

  suppressMessages(rliger::defaultCluster(x = object) <- new_idents_meta$new_idents)
  cli_inform(message = c("v" = "{.code defaultCluster} updated and stored as: {.val defaultCluster} in object cellMeta slot."))

  # Add new ident to cellMeta information if desired
  if (!is.null(x = new_ident_name)) {
    object@cellMeta[[new_ident_name]] <- rliger::defaultCluster(x = object)
    cli_inform(message = c("i" = "{.code new_idents} also stored as: {.val new_ident_name} in object cellMeta slot."))
  }

  # return object
  return(object)
}


#' Subset LIGER object
#'
#' Subset LIGER object by cluster or other meta data variable.
#'
#' @param liger_object LIGER object name.
#' @param cluster Name(s) of cluster to subset from object.
#' @param cluster_col name of `@cellMeta` column containing cluster names, default is "leiden_cluster".
#' @param ident variable within `ident_col` to use in sub-setting object.
#' @param ident_col column in `@cellMeta` that contains values provided to `ident`.
#' @param invert logical, whether to subset the inverse of the clusters or idents provided, default is FALSE.
#'
#' @return liger object
#'
#' @import cli
#' @importFrom dplyr pull filter
#' @importFrom magrittr "%>%"
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # subset clusters 3 and 5
#' sub_liger <- subset_liger(liger_object = liger_object, cluster = c(3, 5))
#'
#' # subset control samples from column "Treatment"
#' sub_liger <- subset_liger(liger_object = liger_object, ident = "control",
#' ident_col = "Treatment")
#'
#' # subset control samples from column "Treatment" in clusters 3 and 5
#' sub_liger <- subset_liger(liger_object = liger_object, ident = "control",
#' ident_col = "Treatment", cluster = c(3, 5))
#'
#' # Remove cluster 9
#' sub_liger <- subset_liger(liger_object = liger_object, cluster = 9, invert = TRUE)
#' }
#'

Subset_LIGER <- function(
    liger_object,
    cluster = NULL,
    cluster_col = "leiden_cluster",
    ident = NULL,
    ident_col = NULL,
    invert = FALSE
) {
  # check liger
  Is_LIGER(liger_object = liger_object)

  # Check new liger object
  if (packageVersion(pkg = 'rliger') < "2.0.0") {
    cli_abort(message = "This function is only for objects with rliger >= v2.0.0")
  }

  # Check value provided
  if (is.null(x = cluster) && is.null(x = ident)) {
    cli_abort(message = "No values provided to subset object")
  }

  # Check meta present
  if (!is.null(x = ident_col)) {
    ident_col <- Meta_Present(object = liger_object, meta_col_names = ident_col, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  # Check meta present
  if (!is.null(x = cluster_col)) {
    cluster_col <- Meta_Present(object = liger_object, meta_col_names = cluster_col, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  # pull meta data
  meta <- Fetch_Meta(object = liger_object)

  # check subset value ok idents
  if (!is.null(x = ident)) {
    ident_values <- meta %>%
      pull(.data[[ident_col]]) %>%
      unique()

    if (!all(ident %in% ident_values)) {
      cli_abort(message = "One or more of provided ident values ({.field {ident}}) were not found in provided ident_col ({.field {ident_col}})")
    }
  }

  # check subset value ok idents
  if (!is.null(x = cluster)) {
    cluster_values <- meta %>%
      pull(.data[[cluster_col]]) %>%
      unique()

    if (!all(cluster %in% cluster_values)) {
      cli_abort(message = "One or more of provided cluster values ({.field {cluster}}) were not found in provided cluster_col ({.field {cluster_col}})")
    }
  }

  # filter just by cluster
  if (!is.null(x = cluster) && is.null(x = ident)) {
    cells_filter <- WhichCells(object = liger_object, ident = cluster, ident_col = cluster_col)
  }

  # filter just by ident
  if (!is.null(x = ident) && is.null(x = cluster)) {
    cells_filter <- WhichCells(object = liger_object, ident = ident, ident_col = ident_col)
  }

  # Filter by ident and cluster
  if (!is.null(x = ident) && !is.null(x = cluster)) {
    cells_filter_cluster <- WhichCells(object = liger_object, ident = cluster, ident_col = cluster_col)
    cells_filter_ident <- WhichCells(object = liger_object, ident = ident, ident_col = ident_col)

    cells_filter <- intersect(x = cells_filter_cluster, y = cells_filter_ident)
  }

  # invert filtering
  if (isTRUE(x = invert)) {
    # get vector of call cells
    all_cells <- Cells(x = liger_object)

    # setdiff to get inverse
    cells_filter <- setdiff(x = all_cells, y = cells_filter)
  }

  # subset object
  sub_obj <- rliger::subsetLiger(object = liger_object, cellIdx = cells_filter)

  return(sub_obj)
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
  # LIGER object check
  Is_LIGER(liger_object = liger_object)

  # check number of factors present
  if (!liger_factor %in% 1:dim(x = liger_object@W)[[1]]) {
    cli_abort(message = c("{.code liger_factor} provided: {.field {liger_factor}} not found",
                          "i" = "{.code liger_object} only contains {.field {dim(x = liger_object@W)[[1]]}} factors.")
    )
  }

  # liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    W <- liger_object@W
    rownames(x = W) <- rownames(x = liger_object@datasets[[1]]@scaleData)
    top_genes <- rownames(x = W)[order(W[, liger_factor], decreasing = TRUE)[1:num_genes]]
    return(top_genes)
  } else {
    # Extract genes
    W <- t(x = liger_object@W)
    rownames(x = W) <- colnames(x = liger_object@scale.data[[1]])
    top_genes <- rownames(x = W)[order(W[, liger_factor], decreasing = TRUE)[1:num_genes]]
    return(top_genes)
  }
}


#' Find Factor Correlations
#'
#' Calculate correlations between gene loadings for all factors in liger object.
#'
#' @param liger_object LIGER object name.
#'
#' @return correlation matrix
#'
#' @import cli
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' factor_correlations <- Find_Factor_Cor(liger_object = object)
#' }
#'

Find_Factor_Cor <- function(
    liger_object
) {
  # check liger
  Is_LIGER(liger_object = liger_object)

  # Check new liger object
  if (packageVersion(pkg = 'rliger') < "2.0.0") {
    cli_abort(message = "This function is only for objects with rliger >= v2.0.0")
  }

  # Get loadings
  factor_loadings <- data.frame(rliger::getMatrix(x = liger_object, slot = "W"))

  # Rename is zero padding
  colnames(x = factor_loadings) <- paste0("Factor_", seq_zeros(seq_length = ncol(x = factor_loadings), num_zeros = 1))

  # Correlation
  cor_mat <- cor(x = factor_loadings)

  return(cor_mat)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### QC UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param species Species of origin for given Seurat Object.  If mouse, human, marmoset, zebrafish, rat,
#' drosophila, rhesus macaque, or chicken (name or abbreviation) are provided the function will automatically
#' generate patterns and features.
#' @param add_mito_ribo logical, whether to add percentage of counts belonging to mitochondrial/ribosomal
#' genes to object (Default is TRUE).
#' @param add_complexity logical, whether to add Cell Complexity to object (Default is TRUE).
#' @param add_top_pct logical, whether to add Top Gene Percentages to object (Default is TRUE).
#' @param add_MSigDB logical, whether to add percentages of counts belonging to genes from of mSigDB hallmark
#' gene lists: "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_APOPTOSIS", and "HALLMARK_DNA_REPAIR" to
#' object (Default is TRUE).
#' @param add_IEG logical, whether to add percentage of counts belonging to IEG genes to object (Default is TRUE).
#' @param add_hemo logical, whether to add percentage of counts belonging to homoglobin genes to object (Default is TRUE).
#' @param mito_name name to use for the new meta.data column containing percent mitochondrial counts.
#' Default is "percent_mito".
#' @param ribo_name name to use for the new meta.data column containing percent ribosomal counts.
#' Default is "percent_ribo".
#' @param mito_ribo_name name to use for the new meta.data column containing percent
#' mitochondrial+ribosomal counts.  Default is "percent_mito_ribo".
#' @param complexity_name name to use for new meta data column for `Add_Cell_Complexity`.
#' Default is "log10GenesPerUMI".
#' @param top_pct_name name to use for new meta data column for `Add_Top_Gene_Pct`.
#' Default is "percent_topXX", where XX is equal to the value provided to `num_top_genes`.
#' @param oxphos_name name to use for new meta data column for percentage of MSigDB oxidative phosphorylation
#' counts.  Default is "percent_oxphos".
#' @param apop_name name to use for new meta data column for percentage of MSigDB apoptosis counts.
#' Default is "percent_apop".
#' @param dna_repair_name name to use for new meta data column for percentage of MSigDB DNA repair
#' counts.  Default is "percent_dna_repair"..
#' @param ieg_name name to use for new meta data column for percentage of IEG counts.  Default is "percent_ieg".
#' @param hemo_name name to use for the new meta.data column containing percent hemoglobin counts.
#' Default is "percent_mito".
#' @param mito_pattern A regex pattern to match features against for mitochondrial genes (will set automatically if
#' species is mouse or human; marmoset features list saved separately).
#' @param ribo_pattern A regex pattern to match features against for ribosomal genes
#' (will set automatically if species is in default list).
#' @param hemo_pattern A regex pattern to match features against for hemoglobin genes
#' (will set automatically if species is in default list).
#' @param mito_features A list of mitochondrial gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param ribo_features A list of ribosomal gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param hemo_features A list of hemoglobin gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param ensembl_ids logical, whether feature names in the object are gene names or
#' ensembl IDs (default is FALSE; set TRUE if feature names are ensembl IDs).
#' @param num_top_genes An integer vector specifying the size(s) of the top set of high-abundance genes.
#' Used to compute the percentage of library size occupied by the most highly expressed genes in each cell.
#' @param overwrite Logical.  Whether to overwrite existing an meta.data column.  Default is FALSE meaning that
#' function will abort if column with name provided to `meta_col_name` is present in meta.data slot.
#'
#' @import cli
#'
#' @return A liger Object
#'
#' @method Add_Cell_QC_Metrics liger
#'
#' @export
#' @rdname Add_Cell_QC_Metrics
#'
#'
#' @concept qc_util
#'
#' @examples
#' \dontrun{
#' obj <- Add_Cell_QC_Metrics(object = obj, species = "Human")
#'}
#'

Add_Cell_QC_Metrics.liger <- function(
    object,
    add_mito_ribo = TRUE,
    add_complexity = TRUE,
    add_top_pct = TRUE,
    add_MSigDB = TRUE,
    add_IEG = TRUE,
    add_hemo = TRUE,
    add_cell_cycle = TRUE,
    species,
    mito_name = "percent_mito",
    ribo_name = "percent_ribo",
    mito_ribo_name = "percent_mito_ribo",
    complexity_name = "log10GenesPerUMI",
    top_pct_name = NULL,
    oxphos_name = "percent_oxphos",
    apop_name = "percent_apop",
    dna_repair_name = "percent_dna_repair",
    ieg_name = "percent_ieg",
    hemo_name = "percent_hemo",
    mito_pattern = NULL,
    ribo_pattern = NULL,
    hemo_pattern = NULL,
    mito_features = NULL,
    ribo_features = NULL,
    hemo_features = NULL,
    ensembl_ids = FALSE,
    num_top_genes = 50,
    assay = NULL,
    overwrite = FALSE,
    ...
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

  # Species Spelling Options
  mouse_options <- accepted_names$Mouse_Options
  human_options <- accepted_names$Human_Options
  marmoset_options <- accepted_names$Marmoset_Options
  zebrafish_options <- accepted_names$Zebrafish_Options
  rat_options <- accepted_names$Rat_Options
  drosophila_options <- accepted_names$Drosophila_Options
  macaque_options <- accepted_names$Macaque_Options

  # Add mito/ribo
  if (isTRUE(x = add_mito_ribo)) {
    cli_inform(message = c("*" = "Adding {.field Mito/Ribo Percentages} to meta.data."))
    liger_object <- Add_Mito_Ribo(object = liger_object, species = species, mito_name = mito_name, ribo_name = ribo_name, mito_ribo_name = mito_ribo_name, mito_pattern = mito_pattern, ribo_pattern = ribo_pattern, mito_features = mito_features, ribo_features = ribo_features, ensembl_ids = ensembl_ids, overwrite = overwrite)
  }

  # Add complexity
  if (isTRUE(x = add_complexity)) {
    cli_inform(message = c("*" = "Adding {.field Cell Complexity #1 (log10GenesPerUMI)} to meta.data."))
    liger_object <- Add_Cell_Complexity(object = liger_object, meta_col_name = complexity_name, overwrite = overwrite)
  }

  # Add top gene expression percent
  if (isTRUE(x = add_top_pct)) {
    cli_inform(message = c("*" = "Adding {.field Cell Complexity #2 (Top {num_top_genes} Percentages)} to meta.data."))
    liger_object <- Add_Top_Gene_Pct(object = liger_object, num_top_genes = num_top_genes, meta_col_name = top_pct_name, overwrite = overwrite)
  }

  # Add MSigDB
  if (isTRUE(x = add_MSigDB)) {
    if (species %in% marmoset_options) {
      cli_warn(message = c("{.val Marmoset} is not currently a part of MSigDB gene list database.",
                           "i" = "No columns will be added to object meta.data"))
    } else {
      cli_inform(message = c("*" = "Adding {.field MSigDB Oxidative Phosphorylation, Apoptosis, and DNA Repair Percentages} to meta.data."))
      liger_object <- Add_MSigDB_LIGER(liger_object = liger_object, species = species, oxphos_name = oxphos_name, apop_name = apop_name, dna_repair_name = dna_repair_name, overwrite = overwrite, ensembl_ids = ensembl_ids)
    }
  }

  # Add IEG
  if (isTRUE(x = add_IEG)) {
    if (species %in% c(marmoset_options, rat_options, zebrafish_options, macaque_options, drosophila_options)) {
      cli_warn(message = c("{.val Rat, Marmoset, Macaque, Zebrafish, and Drosophila} are not currently supported.",
                           "i" = "No column will be added to object meta.data"))
    } else {
      cli_inform(message = c("*" = "Adding {.field IEG Percentages} to meta.data."))
      liger_object <- Add_IEG_LIGER(liger_object = liger_object, species = species, ieg_name = ieg_name, overwrite = overwrite, ensembl_ids = ensembl_ids)
    }
  }

  # Add hemo
  if (isTRUE(x = add_hemo)) {
    cli_inform(message = c("*" = "Adding {.field Hemoglobin Percentages} to meta.data."))
    liger_object <- Add_Hemo(object = liger_object, species = species, hemo_name = hemo_name, hemo_pattern = hemo_pattern, hemo_features = hemo_features, overwrite = overwrite, ensembl_ids = ensembl_ids)
  }

  # return object
  return(liger_object)
}



#' Add Mito and Ribo percentages
#'
#' @param species Species of origin for given Object.  If mouse, human, marmoset, zebrafish, rat,
#' drosophila, rhesus macaque, or chicken (name or abbreviation) are provided the function will automatically
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
  # Accepted species names
  accepted_names <- data.frame(
    Mouse_Options = c("Mouse", "mouse", "Ms", "ms", "Mm", "mm"),
    Human_Options = c("Human", "human", "Hu", "hu", "Hs", "hs"),
    Marmoset_Options = c("Marmoset", "marmoset", "CJ", "Cj", "cj", NA),
    Zebrafish_Options = c("Zebrafish", "zebrafish", "DR", "Dr", "dr", NA),
    Rat_Options = c("Rat", "rat", "RN", "Rn", "rn", NA),
    Drosophila_Options = c("Drosophila", "drosophila", "DM", "Dm", "dm", NA),
    Macaque_Options = c("Macaque", "macaque", "Rhesus", "macaca", "mmulatta", NA),
    Chicken_Options = c("Chicken", "chicken", "Gallus", "gallus", "Gg", "gg")
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
  meta_names <- colnames(x = Fetch_Meta(object = object))

  if (mito_name %in% meta_names || ribo_name %in% meta_names || mito_ribo_name %in% meta_names) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Columns with {.val {mito_name}} and/or {.val {ribo_name}} already present in meta data.",
                            "i" = "*To run function and overwrite columns set parameter {.code overwrite = TRUE} or change respective {.code mito_name}, {.code ribo_name}, and/or {.code mito_ribo_name}.*")
      )
    }
    cli_inform(message = c("Columns with {.val {mito_name}} and/or {.val {ribo_name}} already present in meta data.",
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
  chicken_options <- accepted_names$Chicken_Options

  # Check ensembl vs patterns
  if (isTRUE(x = ensembl_ids) && species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options, chicken_options) && any(!is.null(x = mito_pattern), !is.null(x = ribo_pattern), !is.null(x = mito_features), !is.null(x = ribo_features))) {
    cli_warn(message = c("When using a default species and setting {.code ensembl_ids = TRUE} provided patterns or features are ignored.",
                         "*" = "Supplied {.code mito_pattern}, {.code ribo_pattern}, {.code mito_features}, {.code ribo_features} will be disregarded.")
    )
  }

  # Assign mito/ribo pattern to stored species
  if (species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options, chicken_options) && any(!is.null(x = mito_pattern), !is.null(x = ribo_pattern))) {
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
  if (species %in% c(marmoset_options, macaque_options, chicken_options)) {
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

  all_features <- Features(x = object)

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
    cli_inform(message = "Adding Percent Mitochondrial genes for {.field {species_use}} using gene symbol pattern: {.val {mito_pattern}}.")
    good_mito <- mito_features[mito_features %in% all_features]

    if (packageVersion(pkg = 'rliger') > "1.0.1") {
      object <- rliger::runGeneralQC(object = object, mito = FALSE, ribo = FALSE, hemo = FALSE, features = list(mito_name = good_mito), verbose = FALSE)
    } else {
      percent_mito <- unlist(lapply(object@raw.data, function(x) {
        (Matrix::colSums(x[good_mito, ])/Matrix::colSums(x))*100}))
      object@cell.data[ , mito_name] <- percent_mito
    }
  }

  if (length_ribo_features > 0){
    cli_inform(message = "Adding Percent Ribosomal genes for {.field {species_use}} using gene symbol pattern: {.val {ribo_pattern}}.")
    good_ribo <- ribo_features[ribo_features %in% all_features]

    if (packageVersion(pkg = 'rliger') > "1.0.1") {
      object <- rliger::runGeneralQC(object = object, mito = FALSE, ribo = FALSE, hemo = FALSE, features = list(ribo_name = good_ribo), verbose = FALSE)
    } else {
      percent_ribo <- unlist(lapply(object@raw.data, function(x) {
        (Matrix::colSums(x[good_ribo, ])/Matrix::colSums(x))*100}))
      object@cell.data[ , ribo_name] <- percent_ribo
    }
  }

  # Create combined mito ribo column if both present
  if (length_mito_features > 0 && length_ribo_features > 0) {
    cli_inform(message = "Adding Percent Mito+Ribo by adding Mito & Ribo percentages.")
    if (packageVersion(pkg = 'rliger') > "1.0.1") {
      object@cellMeta[[mito_ribo_name]] <- object@cellMeta[[mito_name]] + object@cellMeta[[ribo_name]]
    } else {
      object_meta <- Fetch_Meta(object = object) %>%
        rownames_to_column("barcodes")

      object_meta <- object_meta %>%
        mutate({{mito_ribo_name}} := .data[[mito_name]] + .data[[ribo_name]])


      object@cell.data[ , mito_ribo_name] <- object_meta[[mito_ribo_name]]
    }
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
  # Check liger
  Is_LIGER(liger_object = object)

  # Check columns for overwrite
  meta_names <- colnames(x = Fetch_Meta(object = object))

  if (meta_col_name %in% meta_names) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Column {.val {meta_col_name}} already present in meta data.",
                            "i" = "*To run function and overwrite column, set parameter {.code overwrite = TRUE} or change respective {.code meta_col_name}*.")
      )
    }
    cli_inform(message = c("Column {.val {meta_col_name}} already present in meta data slot",
                           "i" = "Overwriting those columns as `overwrite = TRUE`.")
    )
  }

  # Add score
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    object@cellMeta[[meta_col_name]] <- log10(object@cellMeta$nGene) / log10(object@cellMeta$nUMI)
  } else {
    object@cell.data[ , meta_col_name] <- log10(object@cell.data$nGene) / log10(object@cell.data$nUMI)
  }

  #return object
  return(object)
}


#' @param species Species of origin for given Object.  If mouse, human, marmoset, zebrafish, rat,
#' drosophila, rhesus macaque, or chicken (name or abbreviation) are provided the function will automatically
#' generate hemo_pattern values.
#' @param hemo_name name to use for the new meta.data column containing percent hemoglobin counts.
#' Default is "percent_hemo".
#' @param hemo_pattern A regex pattern to match features against for hemoglobin genes (will set automatically if
#' species is mouse or human; marmoset features list saved separately).
#' @param hemo_features A list of hemoglobin gene names to be used instead of using regex pattern.
#' @param ensembl_ids logical, whether feature names in the object are gene names or
#' ensembl IDs (default is FALSE; set TRUE if feature names are ensembl IDs).
#' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
#' function will abort if columns with any one of the names provided to `hemo_name` is
#' present in meta.data slot.
#' @param list_species_names returns list of all accepted values to use for default species names which
#' contain internal regex/feature lists (human, mouse, marmoset, zebrafish, rat, drosophila, and
#' rhesus macaque).  Default is FALSE.
#'
#' @import cli
#' @importFrom magrittr "%>%"
#'
#' @method Add_Hemo liger
#'
#' @export
#' @rdname Add_Hemo
#'
#' @concept qc_util
#'
#' @examples
#' \dontrun{
#' # Liger
#' liger_object <- Add_Hemo(object = liger_object, species = "human")
#'}
#'

Add_Hemo.liger <- function(
    object,
    species,
    hemo_name = "percent_hemo",
    hemo_pattern = NULL,
    hemo_features = NULL,
    ensembl_ids = FALSE,
    overwrite = FALSE,
    list_species_names = FALSE,
    ...
) {
  # Accepted species names
  accepted_names <- data.frame(
    Mouse_Options = c("Mouse", "mouse", "Ms", "ms", "Mm", "mm"),
    Human_Options = c("Human", "human", "Hu", "hu", "Hs", "hs"),
    Marmoset_Options = c("Marmoset", "marmoset", "CJ", "Cj", "cj", NA),
    Zebrafish_Options = c("Zebrafish", "zebrafish", "DR", "Dr", "dr", NA),
    Rat_Options = c("Rat", "rat", "RN", "Rn", "rn", NA),
    Drosophila_Options = c("Drosophila", "drosophila", "DM", "Dm", "dm", NA),
    Macaque_Options = c("Macaque", "macaque", "Rhesus", "macaca", "mmulatta", NA),
    Chicken_Options = c("Chicken", "chicken", "Gallus", "gallus", "Gg", "gg")
  )

  # Return list of accepted default species name options
  if (isTRUE(x = list_species_names)) {
    return(accepted_names)
    stop_quietly()
  }

  # Check liger
  Is_LIGER(liger_object = object)

  # Overwrite check
  meta_names <- colnames(x = Fetch_Meta(object = object))

  if (hemo_name %in% meta_names) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Columns with {.val {hemo_name}} already present in meta data.",
                            "i" = "*To run function and overwrite columns set parameter {.code overwrite = TRUE} or change {.code hemo_name}.*")
      )
    }
    cli_inform(message = c("Columns with {.val {hemo_name}} already present in meta data.",
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
  chicken_options <- accepted_names$Chicken_Options

  # Assign hemo pattern to stored species
  if (species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options, macaque_options, chicken_options) && any(!is.null(x = hemo_pattern))) {
    cli_warn(message = c("Pattern expressions for included species are set by default.",
                         "*" = "Supplied {.code hemo_pattern} and {.code hemo_pattern} will be disregarded.",
                         "i" = "To override defaults please supply a feature list for hemo genes.")
    )
  }

  if (species %in% mouse_options) {
    species_use <- "Mouse"
    hemo_pattern <- "^Hb[^(P)]"
  }
  if (species %in% human_options) {
    species_use <- "Human"
    hemo_pattern <- "^HB[^(P)]"
  }
  if (species %in% c(marmoset_options, macaque_options)) {
    species_use <- "Marmoset/Macaque"
    hemo_pattern <- "^^HB[^(P)]"
  }
  if (species %in% zebrafish_options) {
    species_use <- "Zebrafish"
    hemo_pattern <- "^hb[^(P)]"
  }
  if (species %in% rat_options) {
    species_use <- "Rat"
    hemo_pattern <- "^Hb[^(P)]"
  }
  if (species %in% drosophila_options) {
    species_use <- "Drosophila"
    hemo_pattern <- "^glob"
  }
  if (species %in% chicken_options) {
    species_use <- "Chicken"
    hemo_pattern <- "^HB[^(P)]"
  }

  # Check that values are provided for mito and ribo
  if (is.null(x = hemo_pattern) && is.null(x = hemo_features)) {
    cli_abort(message = c("No features or patterns provided for hemo genes.",
                          "i" = "Please provide a default species name or pattern/features."))
  }

  # get all features
  all_features <- Features(x = object, by_dataset = FALSE)

  # Retrieve ensembl ids if TRUE
  if (isTRUE(x = ensembl_ids)) {
    hemo_features <- Retrieve_Ensembl_Hemo(species = species)
  }

  # get features from patterns
  hemo_features <- hemo_features %||% grep(pattern = hemo_pattern, x = all_features, value = TRUE)

  # Check features are present in object
  length_hemo_features <- length(x = intersect(x = hemo_features, y = all_features))

  # Check length of hemo features found in object
  if (length_hemo_features < 1) {
    cli_warn(message = c("No hemoglobin features found in object using pattern/feature list provided.",
                         "i" = "No column will be added to meta.data.")
    )
  }

  # Add hemo column
  cli_inform(message = "Adding Percent Hemoglobin for {.field {species_use}} using gene symbol pattern: {.val {hemo_pattern}}.")
  if (length_hemo_features > 0) {
    good_hemo <- hemo_features[hemo_features %in% all_features]

    if (packageVersion(pkg = 'rliger') > "1.0.1") {
      object <- rliger::runGeneralQC(object = object, mito = FALSE, ribo = FALSE, hemo = FALSE, features = list(hemo_name = good_hemo), verbose = FALSE)
    } else {
      percent_hemo <- unlist(lapply(object@raw.data, function(x) {
        (Matrix::colSums(x[good_hemo, ])/Matrix::colSums(x))*100}))
      object@cell.data[ , hemo_name] <- percent_hemo
    }
  }

  # return final object
  return(object)
}


#' @param num_top_genes An integer vector specifying the size(s) of the top set of high-abundance genes.
#' Used to compute the percentage of library size occupied by the most highly expressed genes in each cell.
#' @param meta_col_name name to use for new meta data column.  Default is "percent_topXX", where XX is
#' equal to the value provided to `num_top_genes`.
#' @param overwrite Logical.  Whether to overwrite existing an meta.data column.  Default is FALSE meaning that
#' function will abort if column with name provided to `meta_col_name` is present in meta.data slot.
#' @param verbose logical, whether to print messages with status updates, default is TRUE.
#'
#' @import cli
#' @importFrom dplyr select all_of bind_rows
#' @importFrom magrittr "%>%"
#' @importFrom rlang is_installed
#'
#' @return A liger Object
#'
#' @method Add_Top_Gene_Pct liger
#'
#' @export
#' @rdname Add_Top_Gene_Pct
#'
#' @concept qc_util
#'
#' @examples
#' \dontrun{
#' liger_object <- Add_Top_Gene_Pct(object = liger_object, num_top_genes = 50)
#' }
#'

Add_Top_Gene_Pct.liger <- function(
    object,
    num_top_genes = 50,
    meta_col_name = NULL,
    overwrite = FALSE,
    verbose = TRUE,
    ...
){
  # Check for scuttle first
  scuttle_check <- is_installed(pkg = "scuttle")
  if (isFALSE(x = scuttle_check)) {
    cli_abort(message = c(
      "Please install the {.val scuttle} package to calculate/add top {num_top_genes} genes percentage.",
      "i" = "This can be accomplished with the following commands: ",
      "----------------------------------------",
      "{.field `install.packages({symbol$dquote_left}BiocManager{symbol$dquote_right})`}",
      "{.field `BiocManager::install({symbol$dquote_left}scuttle{symbol$dquote_right})`}",
      "----------------------------------------"
    ))
  }

  # Check Liger
  Is_LIGER(liger_object = object)

  # Set colnames
  scuttle_colname <- paste0("percent.top_", num_top_genes)
  if (is.null(x = meta_col_name)) {
    meta_col_name <- paste0("percent_top", num_top_genes)
  }

  # Check columns for overwrite
  meta_names <- colnames(x = Fetch_Meta(object = object))

  if (meta_col_name %in% meta_names) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Column {.val {meta_col_name}} already present in meta data.",
                            "i" = "*To run function and overwrite column, set parameter {.code overwrite = TRUE} or change respective {.code meta_col_name}*.")
      )
    }
    cli_inform(message = c("Column {.val {meta_col_name}} already present in meta data.",
                           "i" = "Overwriting those columns as {.code overwrite = TRUE}.")
    )
  }

  # Get number of datasets
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    num_datasets <- length(x = object@datasets)
  } else {
    num_datasets <- length(x = object@raw.data)
  }

  # Extract matrix
  if (isTRUE(x = verbose)) {
    cli_inform(message = "Calculating percent expressing top {num_top_genes} across all datasets.")
  }

  # apply over all datasets
  res_list <- lapply(1:num_datasets, function(x) {
    if (packageVersion(pkg = 'rliger') > "1.0.1") {
      dataset_mat <- rliger::getMatrix(x = object, slot = "rawData")[[x]]
    } else {
      dataset_mat <- object@raw.data[[x]]
    }

    # run scuttle
    dataset_res <- as.data.frame(scuttle::perCellQCMetrics(x = dataset_mat, percent.top = num_top_genes))
    # select results column
    dataset_res <- dataset_res %>%
      select(all_of(scuttle_colname))
  })

  # combine results
  if (isTRUE(x = verbose)) {
    cli_inform(message = "Combining data from all datasets.")
  }
  res <- bind_rows(res_list)

  # Add to object and return
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    object@cellMeta[[meta_col_name]] <- res
  } else {
    object@cell.data[ , meta_col_name] <- res
  }

  # return object
  return(object)
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
#' Only applicable for rliger < 2.0.0.
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
  pt.size = 1.5,
  chunk = 1000
) {
  # check liger
  Is_LIGER(liger_object = liger_object)

  # version check and split
  if (packageVersion(pkg = 'rliger') >= "2.0.0") {
    raw_data <- rliger::rawData(x = liger_object)

    cli_inform(message = "Creating temporary object with combined data.")

    temp_liger <- rliger::createLiger(rawData = list("dataset" = Merge_Sparse_Data_All(matrix_list = raw_data)), removeMissing = FALSE)

    rm(raw_data)
    gc()

    cli_inform(message = "Normalizing and identifying variable features.")

    temp_liger <- rliger::normalize(object = temp_liger)
    temp_liger <- rliger::selectGenes(object = temp_liger, thresh = var.thresh, alpha = alpha.thresh, chunk = chunk)
    if (isTRUE(x = do.plot)) {
      print(plotVarFeatures(object = temp_liger, dotSize = pt.size))
    }

    var_genes <- rliger::varFeatures(x = temp_liger)

    rm(temp_liger)
    gc()

    rliger::varFeatures(x = liger_object) <- var_genes
  } else {
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
  }
  return(liger_object)
}
