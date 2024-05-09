#' Calculate Cluster Stats
#'
#' Calculates both overall and per sample cell number and percentages per cluster based on orig.ident
#'
#' @param seurat_object Seurat object name.
#' @param group_by_var meta data column to classify samples (default = "orig.ident").
#'
#' @import cli
#' @importFrom dplyr left_join rename all_of
#' @importFrom janitor adorn_totals
#' @importFrom magrittr "%>%"
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom tidyr pivot_wider
#'
#' @return A data.frame
#'
#' @export
#'
#' @concept stats
#'
#' @examples
#' \dontrun{
#' stats <- Cluster_Stats_All_Samples(seurat_object = object, group_by_var = "orig.ident")
#' }
#'

Cluster_Stats_All_Samples <- function(
  seurat_object,
  group_by_var = "orig.ident"
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check on meta data column
  possible_meta_col <- colnames(x = seurat_object@meta.data)
  if (!group_by_var %in% possible_meta_col) {
    cli_abort(message = "{.val {group_by_var}} was not found in meta.data slot of Seurat Object.")
  }

  # Extract total percents
  total_percent <- prop.table(x = table(seurat_object@active.ident)) * 100
  total_percent <- data.frame(total_percent) %>%
    rename(Cluster = all_of("Var1")) %>%
    arrange(desc(.data[["Freq"]]))

  # Extract total cell number per cluster across all samples
  total_cells <- table(seurat_object@active.ident) %>%
    data.frame() %>%
    rename(Cluster = all_of("Var1"), Number = all_of("Freq")) %>%
    arrange(desc(.data[["Number"]]))

  # Cluster overall stats across all animals
  cluster_stats <- suppressMessages(left_join(total_cells, total_percent))

  # Extract cells per metadata column per cluster
  cells_per_cluster_2 <- table(seurat_object@active.ident, seurat_object@meta.data[, group_by_var])
  cells_per_cluster_2 <- data.frame(cells_per_cluster_2) %>%
    rename(Cluster = all_of("Var1"), group_by_var = all_of("Var2"), cell_number = all_of("Freq"))

  cells_per_cluster_2 <- cells_per_cluster_2 %>%
    pivot_wider(names_from = group_by_var, values_from = .data[["cell_number"]])

  # Merge cells per metadata column per cluster with cluster stats
  cluster_stats_2 <- suppressMessages(left_join(cluster_stats, cells_per_cluster_2))

  # Calculate and extract percents of cells per cluster per
  percent_per_cluster_2 <- prop.table(x = table(seurat_object@active.ident, seurat_object@meta.data[, group_by_var]), margin = 2) * 100
  percent_per_cluster_2 <- data.frame(percent_per_cluster_2) %>%
    rename(cluster = all_of("Var1"), group_by_var = all_of("Var2"), percent = all_of("Freq"))
  percent_per_cluster_2 <- percent_per_cluster_2 %>%
    pivot_wider(names_from = group_by_var, values_from = .data[["percent"]]) %>%
    column_to_rownames("cluster")
  colnames(x = percent_per_cluster_2) <- paste(colnames(x = percent_per_cluster_2), "%", sep = "_")

  percent_per_cluster_2 <- percent_per_cluster_2 %>%
    rownames_to_column(var = "Cluster")

  # Merge percent cells per metadata column per cluster with cluster stats and add Totals column
  cluster_stats <- suppressMessages(left_join(cluster_stats_2, percent_per_cluster_2)) %>%
    adorn_totals("row")
}


#' Calculate percent of expressing cells
#'
#' Calculates the percent of cells that express a given set of features by various grouping factors
#'
#' @param seurat_object Seurat object name.
#' @param features Feature(s) to plot.
#' @param threshold Expression threshold to use for calculation of percent expressing (default is 0).
#' @param group_by Factor to group the cells by.
#' @param split_by Factor to split the groups by.
#' @param entire_object logical (default = FALSE).  Whether to calculate percent of expressing cells
#' across the entire object as opposed to by cluster or by `group_by` variable.
#' @param assay Assay to pull feature data from.  Default is active assay.
#' @param slot `r lifecycle::badge("deprecated")` soft-deprecated.  See `layer`
#' @param layer Which layer to pull expression data from?  Default is "data".
#'
#' @return A data.frame
#'
#' @references Part of code is modified from Seurat package as used by \code{\link[Seurat]{DotPlot}}
#' to generate values to use for plotting.  Source code can be found here:
#' \url{https://github.com/satijalab/seurat/blob/4e868fcde49dc0a3df47f94f5fb54a421bfdf7bc/R/visualization.R#L3391} (License: GPL-3).
#'
#' @import cli
#'
#' @export
#'
#' @concept stats
#'
#' @examples
#' \dontrun{
#' percent_stats <- Percent_Expressing(seurat_object = object, features = "Cx3cr1", threshold = 0)
#' }
#'

Percent_Expressing <- function(
  seurat_object,
  features,
  threshold = 0,
  group_by = NULL,
  split_by = NULL,
  entire_object = FALSE,
  slot = deprecated(),
  layer = "data",
  assay = NULL
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check is slot is supplied
  if (lifecycle::is_present(slot)) {
    lifecycle::deprecate_warn(when = "2.0.0",
                              what = "Percent_Expressing(slot)",
                              with = "Percent_Expressing(layer)",
                              details = c("v" = "As of Seurat 5.0.0 the {.code slot} parameter is deprecated and replaced with {.code layer}.",
                                          "i" = "Please adjust code now to prepare for full deprecation.")
    )
    layer <- slot
  }

  # set assay (if null set to active assay)
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Check features exist in object
  features_list <- Feature_Present(data = seurat_object, features = features, print_msg = FALSE, case_check = TRUE, seurat_assay = assay)[[1]]

  # Check group_by is in object
  if (!is.null(x = group_by) && group_by == "ident") {
    group_by <- NULL
  }

  if (!is.null(x = group_by)) {
    possible_groups <- colnames(x = seurat_object@meta.data)
    if (!group_by %in% possible_groups) {
      cli_abort("Grouping variable {.val {group_by}} was not found in Seurat Object.")
    }
  }

  # Check split_by is in object
  if (!is.null(x = split_by)) {
    possible_groups <- colnames(x = seurat_object@meta.data)
    if (!split_by %in% possible_groups) {
      cli_abort("Splitting variable {.val {split_by}} was not found in Seurat Object.")
    }
  }

  # Pull Expression Info
  cells <- unlist(x = CellsByIdentities(object = seurat_object, idents = NULL))
  expression_info <- FetchData(object = seurat_object, vars = features_list, cells = cells, layer = layer)

  # Add grouping variable
  if (isTRUE(x = entire_object)) {
    expression_info$id <- "All_Cells"
  } else {
    expression_info$id <- if (is.null(x = group_by)) {
      Idents(object = seurat_object)[cells, drop = TRUE]
    } else {
      seurat_object[[group_by, drop = TRUE]][cells, drop = TRUE]
    }
  }
  if (!is.factor(x = expression_info$id)) {
    expression_info$id <- factor(x = expression_info$id)
  }
  id.levels <- levels(x = expression_info$id)
  expression_info$id <- as.vector(x = expression_info$id)

  # Split data if split.by is true
  if (!is.null(x = split_by)) {
    splits <- seurat_object[[split_by, drop = TRUE]][cells, drop = TRUE]
    expression_info$id <- paste(expression_info$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }

  # Calculate percent expressing
  percent_expressing <- lapply(
    X = unique(x = expression_info$id),
    FUN = function(ident) {
      data.use <- expression_info[expression_info$id == ident, 1:(ncol(x = expression_info) - 1), drop = FALSE]
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove_Seurat, threshold = threshold)
      return(list(pct.exp = pct.exp))
    }
  )
  names(x = percent_expressing) <- unique(x = expression_info$id)

  # Convert & return data.frame
  row_dim_names <- features_list
  col_dim_names <- names(x = percent_expressing)
  mat_dims <- list(row_dim_names, col_dim_names)
  final_df <- data.frame(matrix(unlist(percent_expressing), nrow = length(features_list), byrow = FALSE, dimnames = mat_dims), stringsAsFactors = FALSE)
  return(final_df)
}


#' Median Statistics
#'
#' Get quick values for median Genes, UMIs, %mito per cell grouped by meta.data variable.
#'
#' @param seurat_object Seurat object name.
#' @param group_by_var Column in meta.data slot to group results by (default = "orig.ident").
#' @param default_var logical.  Whether to include the default meta.data variables of: "nCount_RNA",
#' "nFeature_RNA", "percent_mito", "percent_ribo", "percent_mito_ribo", and "log10GenesPerUMI"
#' in addition to variables supplied to `median_var`.
#' @param median_var Column(s) in `@meta.data` to calculate medians for in addition to defaults.
#' Must be of `class()` integer or numeric.
#'
#' @return A data.frame.
#'
#' @importFrom dplyr group_by select summarise any_of across all_of
#' @importFrom magrittr "%>%"
#' @importFrom stats median
#'
#' @export
#'
#' @concept stats
#'
#' @examples
#' \dontrun{
#' med_stats <- Median_Stats(seurat_object - obj, group_by_var = "orig.ident")
#' }
#'

Median_Stats <- function(
  seurat_object,
  group_by_var = "orig.ident",
  default_var = TRUE,
  median_var = NULL
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  if (isTRUE(x = default_var)) {
    default_var <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo", "percent_mito_ribo", "log10GenesPerUMI")
  } else {
    default_var <- NULL
  }

  # Check group variable present
  group_by_var <- Meta_Present(object = seurat_object, meta_col_names = group_by_var, print_msg = FALSE)[[1]]

  # Check stats variables present
  all_variables <- c(default_var, median_var)

  all_variables <- Meta_Present(object = seurat_object, meta_col_names = all_variables, print_msg = FALSE)[[1]]

  # Filter meta data for columns of interest
  meta_numeric_check <- Fetch_Meta(object = seurat_object) %>%
    select(any_of(all_variables))

  all_variables <- Meta_Numeric(data = meta_numeric_check)

  # Create column names for final data frame from valid columns
  all_variable_col_names <- c(group_by_var, paste0("Median_", all_variables))

  # Calculate medians for each group_by
  meta_data <- Fetch_Meta(object = seurat_object)

  median_by_group <- meta_data %>%
    group_by(.data[[group_by_var]]) %>%
    summarise(across(all_of(all_variables), median))

  # Calculate overall medians
  median_overall <- meta_data %>%
    summarise(across(all_of(all_variables), median))

  # Create data.frame with group_by_var as column name
  meta_col_name_df <- data.frame(col_name = "Totals (All Cells)")
  colnames(x = meta_col_name_df) <- group_by_var
  # Merge with overall median data.frame
  median_overall <- cbind(meta_col_name_df, median_overall)

  # Merge by group_by_var and overall median data.frames
  median_all <- rbind(median_by_group, median_overall)

  # Rename columns and return data.frame
  colnames(x = median_all) <- all_variable_col_names

  return(median_all)
}


#' Median Absolute Deviation Statistics
#'
#' Get quick values for X x median absolute deviation for Genes, UMIs, %mito per cell grouped by meta.data variable.
#'
#' @param seurat_object Seurat object name.
#' @param group_by_var Column in meta.data slot to group results by (default = "orig.ident").
#' @param default_var logical.  Whether to include the default meta.data variables of: "nCount_RNA",
#' "nFeature_RNA", "percent_mito", "percent_ribo", "percent_mito_ribo", and "log10GenesPerUMI"
#' in addition to variables supplied to `mad_var`.
#' @param mad_var Column(s) in `@meta.data` to calculate medians for in addition to defaults.
#' Must be of `class()` integer or numeric.
#' @param mad_num integer value to multiply the MAD in returned data.frame (default is 2).
#' Often helpful when calculating a outlier range to base of of median + (X*MAD).
#'
#' @return A data.frame.
#'
#' @importFrom dplyr group_by select summarise any_of across all_of
#' @importFrom magrittr "%>%"
#' @importFrom stats mad
#'
#' @export
#'
#' @concept stats
#'
#' @examples
#' \dontrun{
#' mad_stats <- MAD_Stats(seurat_object = obj, group_by_var = "orig.ident")
#' }
#'

MAD_Stats <- function(
    seurat_object,
    group_by_var = "orig.ident",
    default_var = TRUE,
    mad_var = NULL,
    mad_num = 2
) {
  # check mad_num
  if (mad_num <= 0) {
    cli_abort(message = "The {.code mad_num} parameter must be greater than 0.")
  }

  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  if (isTRUE(x = default_var)) {
    default_var <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo", "percent_mito_ribo", "log10GenesPerUMI")
  } else {
    default_var <- NULL
  }

  # set to active ident if "ident" is provided
  if (group_by_var == "ident") {
    seurat_object[["active.ident"]] <- Idents(object = seurat_object)
    group_by_var <- "active.ident"
  }

  # Check group variable present
  group_by_var <- Meta_Present(object = seurat_object, meta_col_names = group_by_var, print_msg = FALSE)[[1]]

  # Check stats variables present
  all_variables <- c(default_var, mad_var)

  all_variables <- Meta_Present(object = seurat_object, meta_col_names = all_variables, print_msg = FALSE)[[1]]

  # Filter meta data for columns of interest
  meta_numeric_check <- Fetch_Meta(object = seurat_object) %>%
    select(any_of(all_variables))

  all_variables <- Meta_Numeric(data = meta_numeric_check)

  # Create column names for final data frame from valid columns
  all_variable_col_names <- c(group_by_var, paste0("MAD x ", mad_num, " ", all_variables))

  # Calculate medians for each group_by
  meta_data <- Fetch_Meta(object = seurat_object)

  mad_by_group <- meta_data %>%
    group_by(.data[[group_by_var]]) %>%
    summarise(across(all_of(all_variables), mad)*mad_num)

  # Calculate overall medians
  mad_overall <- meta_data %>%
    summarise(across(all_of(all_variables), mad)*mad_num)

  # Create data.frame with group_by_var as column name
  meta_col_name_df <- data.frame(col_name = "Totals (All Cells)")
  colnames(x = meta_col_name_df) <- group_by_var
  # Merge with overall median data.frame
  mad_overall <- cbind(meta_col_name_df, mad_overall)

  # Merge by group_by_var and overall median data.frames
  mad_all <- rbind(mad_by_group, mad_overall)

  # Rename columns and return data.frame
  colnames(x = mad_all) <- all_variable_col_names

  # return data.frame
  return(mad_all)
}


#' CellBender Feature Differences
#'
#' Get quick values for raw counts, CellBender counts, count differences, and percent count differences
#' per feature.
#'
#' @param seurat_object Seurat object name.
#' @param raw_assay Name of the assay containing the raw count data.
#' @param cell_bender_assay Name of the assay containing the CellBender count data.
#' @param raw_mat Name of raw count matrix in environment if not using Seurat object.
#' @param cell_bender_mat Name of CellBender count matrix in environment if not using Seurat object.
#'
#' @return A data.frame containing summed raw counts, CellBender counts, count difference, and
#' percent difference in counts.
#'
#' @import cli
#' @importFrom dplyr arrange desc left_join mutate
#' @importFrom magrittr "%>%"
#' @importFrom Matrix rowSums
#' @importFrom SeuratObject Layers JoinLayers LayerData
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @export
#'
#' @concept stats
#'
#' @examples
#' \dontrun{
#' cb_stats <- CellBender_Feature_Diff(seurat_object - obj, raw_assay = "RAW",
#' cell_bender_assay = "RNA")
#' }
#'

CellBender_Feature_Diff <- function(
  seurat_object = NULL,
  raw_assay = NULL,
  cell_bender_assay = NULL,
  raw_mat = NULL,
  cell_bender_mat = NULL
) {
  if (!is.null(x = seurat_object)) {
    # Is Seurat
    Is_Seurat(seurat_object = seurat_object)

    # Check assays present
    assays_not_found <- Assay_Present(seurat_object = seurat_object, assay_list = c(raw_assay, cell_bender_assay), print_msg = FALSE, omit_warn = TRUE)[[2]]

    if (!is.null(x = assays_not_found)) {
      stop_quietly()
    }

    # Check layers and join if split
    raw_layers_present <- Layers(object = seurat_object, assay = raw_assay, search = 'counts')

    if (length(x = raw_layers_present) > 1) {
      cli_inform(message = c("Multiple raw count layers present in {.field {raw_assay}} assay.",
                             "i" = "Plot will join layers."))
      seurat_object <- JoinLayers(object = seurat_object, assay = raw_assay)
    }

    # Check layers and join if split
    cb_layers_present <- Layers(object = seurat_object, assay = cell_bender_assay, search = 'counts')

    if (length(x = cb_layers_present) > 1) {
      cli_inform(message = c("Multiple raw count layers present in {.field {raw_assay}} assay.",
                             "i" = "Plot will join layers."))
      seurat_object <- JoinLayers(object = seurat_object, assay = cell_bender_assay)
    }

    # Pull raw counts
    raw_counts <- LayerData(object = seurat_object, assay = raw_assay, layer = "counts") %>%
      rowSums() %>%
      data.frame() %>%
      rownames_to_column("Feature_Names")

    colnames(x = raw_counts)[2] <- "Raw_Counts"

    # Pull Cell Bender Counts
    cb_counts <- LayerData(object = seurat_object, assay = cell_bender_assay, layer = "counts") %>%
      rowSums() %>%
      data.frame() %>%
      rownames_to_column("Feature_Names")

    colnames(x = cb_counts)[2] <- "CellBender_Counts"

    # Check features identical
    diff_features <- symdiff(x = raw_counts$Feature_Names, y = cb_counts$Feature_Names)

    if (length(x = diff_features > 0)) {
      cli_warn(message = c("The following features are not present in both assays:",
                           "*" = "{.field {diff_features}}",
                           "i" = "Check matrices used to create object.")
      )
    }

    # merge
    merged_counts <- suppressMessages(left_join(x = raw_counts, y = cb_counts))

    # Add diff and % diff
    merged_counts <- merged_counts %>%
      mutate(Count_Diff = .data[["Raw_Counts"]] - .data[["CellBender_Counts"]],
             Pct_Diff = 100 - ((.data[["CellBender_Counts"]] / .data[["Raw_Counts"]]) * 100)) %>%
      arrange(desc(.data[["Pct_Diff"]])) %>%
      column_to_rownames("Feature_Names")

    # return data
    return(merged_counts)
  }

  # Matrix Version (Need to update warnings and checks here)
  if (!is.null(x = raw_mat) && !is.null(x = cell_bender_mat)) {
    # Pull raw counts
    raw_counts <- raw_mat %>%
      rowSums() %>%
      data.frame() %>%
      rownames_to_column("Feature_Names")

    colnames(x = raw_counts)[2] <- "Raw_Counts"

    # Pull Cell Bender Counts
    cb_counts <- cell_bender_mat %>%
      rowSums() %>%
      data.frame() %>%
      rownames_to_column("Feature_Names")

    colnames(x = cb_counts)[2] <- "CellBender_Counts"

    # Check features identical
    diff_features <- symdiff(x = raw_counts$Feature_Names, y = cb_counts$Feature_Names)

    if (length(x = diff_features > 0)) {
      cli_warn(message = c("The following features are not present in both assays:",
                           "*" = "{.field {diff_features}}",
                           "i" = "Check matrices used to create object.")
      )
    }

    # merge
    merged_counts <- suppressMessages(left_join(x = raw_counts, y = cb_counts))

    # Add diff and % diff
    merged_counts <- merged_counts %>%
      mutate(Count_Diff = .data[["Raw_Counts"]] - .data[["CellBender_Counts"]],
             Pct_Diff = 100 - ((.data[["CellBender_Counts"]] / .data[["Raw_Counts"]]) * 100)) %>%
      arrange(desc(.data[["Pct_Diff"]])) %>%
      column_to_rownames("Feature_Names")

    # return data
    return(merged_counts)
  }
}
