#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GENERAL OBJECT UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Merge a list of Seurat Objects
#'
#' Enables easy merge of a list of Seurat Objects.  See  See \code{\link[SeuratObject]{merge}} for more information,
#'
#' @param list_seurat list composed of multiple Seurat Objects.
#' @param add.cell.ids A character vector of equal length to the number of objects in `list_seurat`.
#' Appends the corresponding values to the start of each objects' cell names.  See \code{\link[SeuratObject]{merge}}.
#' @param merge.data Merge the data slots instead of just merging the counts (which requires renormalization).
#' This is recommended if the same normalization approach was applied to all objects.
#' See \code{\link[SeuratObject]{merge}}.
#' @param project Project name for the Seurat object. See \code{\link[SeuratObject]{merge}}.
#'
#' @import cli
#' @importFrom magrittr "%>%"
#' @importFrom purrr reduce
#'
#' @return A Seurat Object
#'
#' @export
#'
#' @concept misc_util
#'
#' @examples
#' \dontrun{
#' object_list <- list(obj1, obj2, obj3, ...)
#' merged_object <- Merge_Seurat_List(list_seurat = object_list)
#' }
#'

Merge_Seurat_List <- function(
    list_seurat,
    add.cell.ids = NULL,
    merge.data = TRUE,
    project = "SeuratProject"
) {
  # Check list_seurat is list
  if (!inherits(x = list_seurat, what = "list")) {
    cli_abort(message = "{.code list_seurat} must be environmental variable of class {.val list}")
  }

  # Check list_seurat is only composed of Seurat objects
  for (i in 1:length(x = list_seurat)) {
    if (!inherits(x = list_seurat[[i]], what = "Seurat")) {
      cli_abort("One or more of entries in {.code list_seurat} are not objects of class {.val Seurat}")
    }
  }

  # Check all barcodes are unique to begin with
  duplicated_barcodes <- list_seurat %>%
    lapply(colnames) %>%
    unlist() %>%
    duplicated() %>%
    any()

  if (isTRUE(x = duplicated_barcodes) && is.null(x = add.cell.ids)) {
    cli_abort(message = c("There are overlapping cell barcodes present in the input objects",
                          "i" = "Please rename cells or provide prefixes to {.code add.cell.ids} parameter to make unique.")
    )
  }

  # Check right number of suffix/prefix ids are provided
  if (!is.null(x = add.cell.ids) && length(x = add.cell.ids) != length(x = list_seurat)) {
    cli_abort(message = "The number of prefixes in {.code add.cell.ids} must be equal to the number of objects supplied to {.code list_seurat}.")
  }

  # Rename cells if provided
  list_seurat <- lapply(1:length(x = list_seurat), function(x) {
    list_seurat[[x]] <- RenameCells(object = list_seurat[[x]], add.cell.id = add.cell.ids[x])
  })

  # Merge objects
  merged_object <- reduce(list_seurat, function(x, y) {
    merge(x = x, y = y, merge.data = merge.data, project = project)
  })
}



#' Re-filter Seurat object
#'
#' Allows for re-filtering of Seurat object based on new parameters for `min.cells` and
#' `min.features` (see \code{\link[SeuratObject]{CreateSeuratObject}} for more details)
#'
#' @param seurat_object Seurat object to filter
#' @param min.cells Include features detected in at least this many cells. Will recalculate nCount and nFeature
#' meta.data values as well.
#' @param min.features Include cells where at least this many features are detected.
#' @param override logical, override the Yes/No interactive check (see details).  Default is FALSE; don't override.
#' @param verbose logical, whether to print information on filtering parameters and number of cells/features
#' removed, Default is TRUE.
#'
#' @returns Seurat object
#'
#' @details
#' When running this function any existing reductions, graphs, and all layers except "counts" in the
#' RNA assay.  None of these aspects will be valid once cells/features are removed.
#' To ensure users understand this default behavior of function will provide interactive prompt that
#' users must select "Yes" in order to continue. To avoid this behavior users can set `override = TRUE` and
#' function will skip the interactive prompt.
#'
#'
#' @export
#'
#' @concept misc_util
#'
#' @examples
#' \dontrun{
#' # Remove features expressed in fewer than 10 cells
#' obj_fil <- ReFilter_SeuratObject(seurat_object = obj, min.cells = 10)
#'
#' # Remove cells with fewer than 1000 features
#' obj_fil <- ReFilter_SeuratObject(seurat_object = obj, min.features = 1000)
#'
#' # Filter on both parameters
#' obj_fil <- ReFilter_SeuratObject(seurat_object = obj, min.features = 1000, min.cells = 10)
#' }


ReFilter_SeuratObject <- function(
    seurat_object,
    min.cells = NULL,
    min.features = NULL,
    override = FALSE,
    verbose = TRUE
) {
  # Check Seurat
  Is_Seurat(seurat_object = pbmc)

  # get assays
  assays <- Assays(object = seurat_object)

  # check filtering parameters
  if (is.null(x = min.cells) && is.null(x = min.features)) {
    cli_abort(message = "Must provide a value to either or both {.code min.cells} or {.code min.features}")
  }

  # Set defaults if one value is NULL
  min.cells <- min.cells %||% 0
  min.features <- min.features %||% 0

  # check for Seurat 5
  assay_layer_check <- unlist(lapply(1:length(x = assays), function(x) {
    Assay5_Check(seurat_object = seurat_object, assay = x)
  }))

  # Check for multiple layers and abort if TRUE
  if (isTRUE(x = any(assay_layer_check))) {
    layers_check <- Layers(object = seurat_object, search = "counts")
    if (length(x = layers_check) > 1) {
      cli_abort(message = c("Multiple layers present {.field {head(x = layers_check, n = 2)}}.",
                            "i" = "Please run {.code JoinLayers} before running {.code ReFilter_SeuratObject}"))
    }
  }

  # Add warning about removal
  if (isFALSE(x = override)) {
    if (yesno(c("This function will remove all {.field reductions, graphs, non-RNA assays} and all layers except {.val counts} in RNA assay",  "\nDo you still want to proceed?"))) {
      return(invisible())
    }
  }

  # Diet Object
  DefaultAssay(object = seurat_object) <- "RNA"
  seurat_object <- DietSeurat(object = seurat_object, layers = "counts", assays = "RNA")

  # Get cells to keep
  nfeatures <- Matrix::colSums(x = LayerData(object = seurat_object, layer = "counts") > 0)
  cells_fil <- colnames(LayerData(object = seurat_object, layer = "counts")[, which(x = nfeatures >= min.features)])
  num_cells_rem <- length(x = Cells(x = seurat_object)) - length(x = cells_fil)

  # Get features to keep
  num.cells <- Matrix::rowSums(x = LayerData(object = seurat_object, layer = "counts", cells = cells_fil) > 0)
  features_fil <- rownames(LayerData(object = seurat_object, layer = "counts", cells = cells_fil)[which(x = num.cells >= min.cells), ])
  num_features_rem <- length(x = Features(x = seurat_object)) - length(x = features_fil)

  # print parameters of filtering
  if (isTRUE(x = verbose)) {
    cli_inform(message = c("*" = "Filtering object",
                           "i" = "Keeping cells with greater than or equal to {.field {min.features} features}.",
                           "{col_green({symbol$double_line})} Total of {.field {num_cells_rem} cells} being removed.",
                           "i" = "Keeping features expressed in greater than or equal to {.field {min.cells} cells}.",
                           "{col_green({symbol$arrow_right})} Total of {.field {num_features_rem} features} being removed."))
  }

  # subset object
  seurat_object <- subset(x = seurat_object, cells = cells_fil, features = features_fil)

  # return object
  return(seurat_object)
}





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### META DATA UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Remove meta data columns containing Seurat Defaults
#'
#' Remove any columns from new meta_data data.frame in preparation for adding back to Seurat Object
#'
#' @param meta_data data.frame containing meta data.
#' @param seurat_object object name.
#' @param barcodes_to_rownames logical, are barcodes present as column and should they be moved to
#' rownames (to be compatible with `Seurat::AddMetaData`).  Default is FALSE.
#' @param barcodes_colname name of barcodes column in meta_data.  Required if `barcodes_to_rownames = TRUE`.
#'
#' @import cli
#' @importFrom dplyr select all_of
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#'
#' @return data.frame with only new columns.
#'
#' @export
#'
#' @concept get_set_util
#'
#' @examples
#' \dontrun{
#' new_meta <- Meta_Remove_Seurat(meta_data = meta_data_df, seurat_object = object)
#' object <- AddMetaData(object = object, metadata = new_meta)
#' }
#'

Meta_Remove_Seurat <- function(
  meta_data,
  seurat_object,
  barcodes_to_rownames = FALSE,
  barcodes_colname = "barcodes"
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Get existing meta.data column names
  existing_names <- colnames(x = seurat_object@meta.data)

  # Filter meta_data to remove already existing columns
  meta_data_filtered <- meta_data %>%
    select(-all_of(x = existing_names))

  if (isTRUE(x = barcodes_to_rownames)) {
    # Check barcodes colname exists
    if (!barcodes_colname %in% colnames(x = meta_data)) {
      cli_abort(message = "{.code barcodes_colname}: {.val {barcodes_colname}} was not present in the column names of meta_data data.frame provided.")
    }
    # Move barcodes to rownames
    meta_data_filtered <- meta_data_filtered %>%
      column_to_rownames(barcodes_colname)
    # return data.frame
    return(meta_data_filtered)
  }

  # Return filtered data.frame
  return(meta_data_filtered)
}


#' Add Sample Level Meta Data
#'
#' Add meta data from ample level data.frame/tibble to cell level seurat `@meta.data` slot
#'
#' @param seurat_object object name.
#' @param meta_data data.frame/tibble containing meta data or path to file to read.  Must be formatted as
#' either data.frame or tibble.
#' @param join_by_seurat name of the column in `seurat_object@meta.data` that contains matching
#' variables to `join_by_meta` in `meta_data.`
#' @param join_by_meta name of the column in `meta_data` that contains matching
#' variables to `join_by_seurat` in `seurat_object@meta.data`.
#' @param na_ok logical, is it ok to add NA values to `seurat_object@meta.data`.  Default is FALSE.
#' Be very careful if setting TRUE because if there is error in join operation it may result in all
#' `@meta.data` values being replaced with NA.
#' @param overwrite logical, if there are shared columns between `seurat_object@meta.data` and `meta_data`
#' should the current `seurat_object@meta.data` columns be overwritten.  Default is FALSE.  This parameter
#' excludes values provided to `join_by_seurat` and `join_by_meta`.
#'
#' @import cli
#' @importFrom data.table fread
#' @importFrom dplyr select left_join all_of
#' @importFrom magrittr "%>%"
#' @importFrom stats setNames
#' @importFrom tibble column_to_rownames rownames_to_column
#'
#' @return Seurat object with new `@meta.data` columns
#'
#' @export
#'
#' @concept get_set_util
#'
#' @examples
#' \dontrun{
#' # meta_data present in environment
#' sample_level_meta <- data.frame(...)
#' obj <- Add_Sample_Meta(seurat_object = obj, meta_data = sample_level_meta,
#' join_by_seurat = "orig.ident", join_by_meta = "sample_ID")
#'
#' # from meta data file
#' obj <- Add_Sample_Meta(seurat_object = obj, meta_data = "meta_data/sample_level_meta.csv",
#' join_by_seurat = "orig.ident", join_by_meta = "sample_ID")
#' }
#'

Add_Sample_Meta <- function(
  seurat_object,
  meta_data,
  join_by_seurat,
  join_by_meta,
  na_ok = FALSE,
  overwrite = FALSE
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check variable vs. file path
  if (isTRUE(x = exists(x = deparse(expr = substitute(expr = meta_data))))) {
    meta_data <- meta_data
  } else {
    if (file.exists(meta_data)) {
      meta_data <- fread(file = meta_data, data.table = FALSE)
    } else {
      cli_abort(message = c("Could not find {.code meta_data} {.field {deparse(expr = substitute(expr = meta_data))}}.",
                            "*" = "If providing environmental variable please check {.code meta_data} name.",
                            "i" = "If providing path to file please check path is correct.")
      )
    }
  }

  # Check meta data structure
  if (!inherits(x = meta_data, what = "data.frame")) {
    cli_abort(message = c("{.code meta_data} not in correct format",
                          "*" = "{.code meta_data} must be a data.frame or tibble.",
                          "i" = "Change format and re-run function.")
    )
  }

  # Check NA in meta data
  if (anyNA(x = meta_data) && isFALSE(x = na_ok)) {
    cli_abort(message = c("{.code meta_data} contains NA values.",
                          "i" = "If you would like NA values added to Seurat meta data please set {.code na_ok = TRUE}.")
    )
  }

  # Check join variables exist
  if (!join_by_seurat %in% colnames(x = seurat_object@meta.data)) {
    cli_abort(message = "The column {.val {join_by_seurat}} was not found in object @meta.data slot."
    )
  }

  if (!join_by_meta %in% colnames(x = meta_data)) {
    cli_abort(message = "The column {.val {join_by_meta}} was not found in supplied `meta_data`."
    )
  }

  # Check if any duplicate column names
  dup_columns <- colnames(x = meta_data)[colnames(x = meta_data) %in% colnames(x = seurat_object@meta.data)]

  if (length(x = dup_columns) > 0) {
    dup_columns <- dup_columns[!dup_columns %in% c(join_by_seurat, join_by_meta)]

    if (any(dup_columns %in% colnames(x = seurat_object@meta.data)) && !overwrite) {
      cli_abort(message = c(" Duplicate {.code meta_data} contains column names in object @meta.data.",
                            "i" = "{.code meta_data} and object@meta.data both contain columns: {.field {glue_collapse_scCustom(input_string = dup_columns)}}.",
                            "*" = "To overwrite existing object @meta.data columns with those in {.code meta_data} set {.code overwrite = TRUE}.")
      )
    }
  }

  # Pull meta data
  meta_seurat <- seurat_object@meta.data %>%
    rownames_to_column("temp_barcodes")

  # remove
  if (isTRUE(x = overwrite)) {
    meta_seurat <- meta_seurat %>%
      select(-all_of(x = dup_columns))
  } else {
    meta_seurat <- meta_seurat
  }

  # join
  meta_merged <- left_join(x = meta_seurat, y = meta_data, by = setNames(join_by_meta, join_by_seurat))

  # Remove existing Seurat meta
  if (length(x = dup_columns) > 0 && isTRUE(x = overwrite)) {
    meta_merged <- meta_merged %>%
      column_to_rownames("temp_barcodes")
  } else {
    meta_merged <- Meta_Remove_Seurat(meta_data = meta_merged, seurat_object = seurat_object) %>%
      column_to_rownames("temp_barcodes")
  }

  # check NA
  if (anyNA(x = meta_merged) && !na_ok) {
    cli_abort(message = c("NAs found in new seurat meta.data.",
                          "*" = "No new meta data added to Seurat object.",
                          "i" = "Check to make sure all levels of joining factors are present in both sets of meta data.")
    )
  }

  # Add meta data
  seurat_object <- AddMetaData(object = seurat_object, metadata = meta_merged)

  # Return object
  return(seurat_object)
}


#' Extract sample level meta.data
#'
#' Returns a by identity meta.data data.frame with one row per sample.  Useful for downstream
#' quick view of sample breakdown, meta data table creation, and/or use in pseudobulk analysis
#'
#' @param object Seurat or LIGER object
#' @param sample_name meta.data column to use as sample.  Output data.frame will contain one row per
#' level or unique value in this variable.
#' @param variables_include `@meta.data` columns to keep in final data.frame.  All other columns will
#' be discarded.  Default is NULL.
#' @param variables_exclude columns to discard in final data.frame.  Many cell level columns are
#' irrelevant at the sample level (e.g., nFeature_RNA, percent_mito).
#' \itemize{
#' \item Default parameter value is `NULL` but internally will set to discard nFeature_ASSAY(s),
#' nCount_ASSAY(s), percent_mito, percent_ribo, percent_mito_ribo, and log10GenesPerUMI.
#' \item If sample level median values are desired for these type of variables the output of this
#' function can be joined with output of \code{\link[scCustomize]{Median_Stats}}.
#' \item Set parameter to `include_all = TRUE` to prevent any columns from being excluded.
#' }
#' @param include_all logical, whether or not to include all object meta data columns in output data.frame.
#' Default is FALSE.
#'
#' @return Returns a data.frame with one row per `sample_name`.
#'
#' @import cli
#' @importFrom dplyr any_of grouped_df select slice
#' @importFrom magrittr "%>%"
#'
#' @export
#'
#' @concept get_set_util
#'
#' @examples
#' library(Seurat)
#' pbmc_small[["batch"]] <- sample(c("batch1", "batch2"), size = ncol(pbmc_small), replace = TRUE)
#'
#' sample_meta <- Extract_Sample_Meta(object = pbmc_small, sample_name = "orig.ident")
#'
#' # Only return specific columns from meta data (orig.ident and batch)
#' sample_meta2 <- Extract_Sample_Meta(object = pbmc_small, sample_name = "orig.ident",
#' variables_include = "batch")
#'
#' # Return all columns from meta data
#' sample_meta3 <- Extract_Sample_Meta(object = pbmc_small, sample_name = "orig.ident",
#' include_all = TRUE)
#'

Extract_Sample_Meta <- function(
  object,
  sample_name = "orig.ident",
  variables_include = NULL,
  variables_exclude = NULL,
  include_all = FALSE
) {
  # Pull meta data
  meta_df <- Fetch_Meta(object = object)

  # Check sample name parameter is present
  if (!sample_name %in% colnames(x = meta_df)) {
    cli_abort(message = "The {.code sample_name} parameter: {.val {sample_name}} was not found in object meta.data")
  }

  if (!is.null(x = variables_include) && !is.null(x = variables_exclude) && include_all) {
    cli_abort(message = "If {.code include_all = TRUE} then {.code variables_include} and {.code variables_exclude} must be NULL.")
  }

  # Generate nCount and nFeature variable vectors for exclusion
  if (is.null(x = variables_exclude)) {
    nFeature_cols <- grep(x = colnames(x = meta_df), pattern = "^nFeature", value = TRUE)

    nCount_cols <- grep(x = colnames(x = meta_df), pattern = "^nCount", value = TRUE)

    combined_exclude <- c(nFeature_cols, nCount_cols, "percent_mito", "percent_ribo", "percent_mito_ribo", "log10GenesPerUMI")

    variables_exclude <- Meta_Present(object = object, meta_col_names = combined_exclude, omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)[[1]]
  }

  # Ensure include exclude are unique
  if (!is.null(x = variables_include) && !is.null(x = variables_exclude)) {
    variable_intersect <- intersect(x = variables_include, y = variables_exclude)
    if (length(x = variable_intersect > 0)) {
      cli_abort(message = c("Following variable is present in both {.code variables_include} and {.code variables_exclude} parameters.",
                            "i" = "{.field {variable_intersect}}")
      )
    }
  }

  # Check variables include/exclude are present
  if (!is.null(x = variables_include)) {
    include_meta_list <- Meta_Present(object = object, meta_col_names = variables_include, omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)
  } else {
    include_meta_list <- NULL
  }

  if (!is.null(x = variables_exclude)) {
    exclude_meta_list <- Meta_Present(object = object, meta_col_names = variables_exclude, omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)
  } else {
    exclude_meta_list <- NULL
  }

  all_not_found_features <- c(include_meta_list[[2]], exclude_meta_list[[2]])

  all_found_features <- c(include_meta_list[[1]], exclude_meta_list[[1]])

  # Stop if no features found
  if (length(x = all_found_features) < 1) {
    cli_abort(message = c("No features were found.",
                          "*" = "The following are not present in object meta data:",
                          "i" = "{.field {glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}}")
    )
  }

  # Return message of features not found
  if (length(x = all_not_found_features) > 0) {
    op <- options(warn = 1)
    on.exit(options(op))
    cli_warn(message = c("The following meta data variables were omitted as they were not found:",
                         "i" = "{.field {glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}}")
    )
  }

  # Create by sample data.frame
  sample_meta_df <- meta_df %>%
    grouped_df(vars = sample_name) %>%
    slice(1)

  # remove rownames
  rownames(x = sample_meta_df) <- NULL

  # Filter data.frame
  if (isTRUE(x = include_all)) {
    sample_meta_df_filtered <- sample_meta_df
  } else {
    if (length(x = include_meta_list[[1]]) > 0) {
      sample_meta_df_filtered <- sample_meta_df %>%
        select(any_of(c(include_meta_list[[1]], sample_name)))
      if (length(x = exclude_meta_list[[1]]) > 0) {
        sample_meta_df_filtered <- sample_meta_df_filtered %>%
          select(-any_of(exclude_meta_list[[1]]))
      }
    } else {
      if (length(x = exclude_meta_list[[1]]) > 0) {
        sample_meta_df_filtered <- sample_meta_df %>%
          select(-any_of(exclude_meta_list[[1]]))
      }
    }
  }

  # return data
  return(sample_meta_df_filtered)
}


#' @rdname Fetch_Meta
#' @importFrom methods slot
#' @export
#' @concept get_set_util
#' @method Fetch_Meta Seurat

Fetch_Meta.Seurat <- function(
    object,
    ...
) {
  # Pull meta data
  object_meta <- slot(object = object, name = "meta.data")

  return(object_meta)
}


#' Randomly downsample by identity
#'
#' Get a randomly downsampled set of cell barcodes with even numbers of cells for each identity class.
#' Can return either as a list (1 entry per identity class) or vector of barcodes.
#'
#' @param seurat_object Seurat object
#' @param num_cells number of cells per ident to use in down-sampling.  This value must be less than or
#' equal to the size of ident with fewest cells.  Alternatively, can set to "min" which will
#' use the maximum number of barcodes based on size of smallest group.
#' @param group.by The ident to use to group cells.  Default is NULL which use current active.ident.  .
#' @param return_list logical, whether or not to return the results as list instead of vector, default is
#' FALSE.
#' @param allow_lower logical, if number of cells in identity is lower than `num_cells` keep the
#' maximum number of cells, default is FALSE.  If FALSE will report error message if `num_cells` is
#' too high, if TRUE will subset cells with more than `num_cells` to that value and those with less
#' than `num_cells` will not be downsampled.
#' @param seed random seed to use for downsampling.  Default is 123.
#'
#' @import cli
#' @importFrom dplyr filter
#' @importFrom magrittr "%>%"
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @return either a vector or list of cell barcodes
#'
#' @export
#'
#' @concept get_set_util
#'
#' @examples
#' library(Seurat)
#'
#' # return vector of barcodes
#' random_cells <- Random_Cells_Downsample(seurat_object = pbmc_small, num_cells = 10)
#' head(random_cells)
#'
#' # return list
#' random_cells_list <- Random_Cells_Downsample(seurat_object = pbmc_small, return_list = TRUE,
#' num_cells = 10)
#' head(random_cells_list)
#'
#' # return max total number of cells (setting `num_cells = "min`)
#' random_cells_max <- Random_Cells_Downsample(seurat_object = pbmc_small, num_cells = "min")
#'

Random_Cells_Downsample <- function(
    seurat_object,
    num_cells,
    group.by = NULL,
    return_list = FALSE,
    allow_lower = FALSE,
    seed = 123
) {
  # Check seurat
  Is_Seurat(seurat_object = seurat_object)

  # set ident in case of NULL
  group.by <- group.by %||% "ident"

  # Check and set idents if not "ident"
  if (group.by != "ident") {
    group.by <- Meta_Present(object = seurat_object, meta_col_names = group.by, print_msg = FALSE, omit_warn = FALSE)[[1]]

    Idents(object = seurat_object) <- group.by
  }

  # get barcodes
  cluster_barcodes <- FetchData(object = seurat_object, vars = "ident") %>%
    rownames_to_column("barcodes")

  # get unique ident vector
  idents_all <- as.character(x = levels(x = Idents(object = seurat_object)))

  # Find minimum length ident and warn if num_cells not equal or lower
  min_cells <- CellsByIdentities(object = seurat_object)
  ident_lengths <- lengths(x = min_cells)
  min_cells <- min(ident_lengths)
  min_cell_ident <- names(x = which(x = ident_lengths ==  min(ident_lengths)))

  if (isFALSE(x = allow_lower)) {
    # set num_cells if value is "min"
    if (num_cells == "min") {
      cli_inform(message = c("The number of cells was set to {.val min}, returning {.field {min_cells}} cells per identity class (equal to size of smallest identity class(es): {.val {min_cell_ident}}).",
                             "{col_green({symbol$tick})} Total of {.field {min_cells * length(x = idents_all)}} cells across whole object."))
      num_cells <- min_cells
    }

    # check size of num_cells
    if (min_cells < num_cells) {
      cli_abort(message = c("The {.code num_cells} value ({.field {num_cells}}) must be lower than or equal to the number of cells in the smallest identity.",
                            "i" = "The identity/identities {.val {min_cell_ident}} contains {.field {min_cells}} cells)."))
    }
  }

  # set values if force
  if (isTRUE(x = allow_lower)) {
    num_cells <- unlist(x = lapply(1:length(x = ident_lengths), function(x){
      val <- ident_lengths[x]
      val <- replace(val, val < num_cells, val)
      val <- replace(val, val > num_cells, num_cells)
      val
    }))
  } else {
    num_cells <- rep(num_cells, length(x = idents_all))
  }

  # set seed and select random cells per ident
  set.seed(seed = seed)

  random_cells <- lapply(1:length(x = idents_all), function(x) {
    clus_barcodes <- cluster_barcodes %>%
      filter(.data[["ident"]] == idents_all[x]) %>%
      column_to_rownames("barcodes") %>%
      rownames() %>%
      sample(size = num_cells[x])
  })

  # return downsampled cells
  if (isTRUE(x = return_list)) {
    # return list
    return(random_cells)
  } else {
    # unlist to vector
    random_cells_unlist <- unlist(x = random_cells)

    # return vector
    return(random_cells_unlist)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### MISC OBJECT UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Store misc data in Seurat object
#'
#' Wrapper function save variety of data types to the `object@misc` slot of Seurat object.
#'
#' @param seurat_object object name.
#' @param data_to_store data to be stored in `@misc` slot.  Can be single piece of data or list.
#' If list of data see `list_as_list` parameter for control over data storage.
#' @param data_name name to give the entry in `@misc` slot.  Must be of equal length of the number
#' of data items being stored.
#' @param list_as_list logical.  If `data_to_store` is a list, this dictates whether to store in `@misc` slot
#' as list (TRUE) or whether to store each entry in the list separately (FALSE).  Default is FALSE.
#' @param overwrite Logical.  Whether to overwrite existing items with the same name.  Default is FALSE, meaning
#' that function will abort if item with `data_name` is present in misc slot.
#' @param verbose logical, whether to print messages when running function, default is TRUE.
#'
#' @return Seurat Object with new entries in the `@misc` slot.
#'
#' @import cli
#'
#' @export
#'
#' @concept get_set_util
#'
#' @examples
#' library(Seurat)
#' clu_pal <- c("red", "green", "blue")
#'
#' pbmc_small <- Store_Misc_Info_Seurat(seurat_object = pbmc_small, data_to_store = clu_pal,
#' data_name = "rd1_colors")
#'

Store_Misc_Info_Seurat <- function(
  seurat_object,
  data_to_store,
  data_name,
  list_as_list = FALSE,
  overwrite = FALSE,
  verbose = TRUE
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check if name already present
  misc_present <- names(x = seurat_object@misc)
  if (data_name %in% misc_present) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Item(s) named: {.val {data_name}} already present in @misc slot.",
                            "i" = "*To run function and overwrite items set parameter {.code overwrite = TRUE} or change {.code data_name}*")
      )
    } else {
      cli_inform(message = c("Items named {.val {data_name}} already present in @misc slot.",
                             "i" = "Overwriting those items as overwrite = TRUE.")
      )
    }
  }

  # Commenting our for now.  Not sure why you would need this check...
  # # Check length of data
  # if (length(x = data_to_store) != 1 && class(x = data_to_store) != "list") {
  #   stop("'data_to_store' must be either single vector/string or a single list of items.")
  # }

  # Check class of data
  if (inherits(x = data_to_store, what = "list")) {
    if (isTRUE(x = list_as_list)) {
      # Check length of name
      if (length(x = data_name) != 1) {
        cli_abort(message = "When storing as list the length {.code data_name} must be {.field 1 (one)}.")
      }

      # Add data
      Misc(object = seurat_object, slot = data_name) <- data_to_store
      # seurat_object@misc[[data_name]] <- data_to_store
      if (isTRUE(x = verbose)) {
        cli_inform(message = c("Seurat Object now contains the following items in @misc slot: ",
                               "i" = "{.field {paste(shQuote(names(x = seurat_object@misc)), collapse=', ')}}")
        )
      }
      return(seurat_object)
    }

    # length of list
    data_list_length <- length(x = data_to_store)

    if (length(x = data_name) != data_list_length) {
      cli_abort(message = "The lengths of {.code data_to_store} ({.field {data_list_length}}) and {.code data_name} ({.field {length(x = data_name)}}) must be equal.")
    }

    # Add data
    for (i in 1:data_list_length) {
      Misc(object = seurat_object, slot = data_name[[i]]) <- data_to_store[[i]]
      # seurat_object@misc[[data_name[i]]] <- data_to_store[[i]]
    }
    if (isTRUE(x = verbose)) {
      cli_inform(message = c("Seurat Object now contains the following items in @misc slot: ",
                             "i" = "{.field {paste(shQuote(names(x = seurat_object@misc)), collapse=', ')}}")
      )
    }

    return(seurat_object)
  } else {
    # Check length of name
    if (length(x = data_name) != 1) {
      cli_abort(message = "When storing a string/vector the length {.code data_name} must be {.field 1 (one)}.")
    }

    # Add data
    Misc(object = seurat_object, slot = data_name) <- data_to_store
    # seurat_object@misc[[data_name]] <- data_to_store
    misc_names <- shQuote(string = names(x = seurat_object@misc))
    if (isTRUE(x = verbose)) {
      cli_inform(message = c("Seurat Object now contains the following items in @misc slot: ",
                             "i" = "{.field {glue_collapse_scCustom(input_string = misc_names, and = TRUE)}}")
      )
    }

    return(seurat_object)
  }
}


#' Store color palette in Seurat object
#'
#' Wrapper function around `Store_Misc_Info_Seurat` to store color palettes.
#'
#' @param seurat_object object name.
#' @param palette vector or list of vectors containing color palettes to store.  If list of palettes
#' see `list_as_list` parameter for control over data storage.
#' @param palette_name name to give the palette(s) in `@misc` slot.  Must be of equal length to the number
#' of data items being stored.
#' @param list_as_list logical.  If `data_to_store` is a list, this dictates whether to store in `@misc` slot
#' as list (TRUE) or whether to store each entry in the list separately (FALSE).  Default is FALSE.
#' @param overwrite Logical.  Whether to overwrite existing items with the same name.  Default is FALSE, meaning
#' that function will abort if item with `data_name` is present in misc slot.
#' @param verbose logical, whether to print messages when running function, default is TRUE.
#'
#' @return Seurat Object with new entries in the `@misc` slot.
#'
#' @export
#'
#' @concept get_set_util
#'
#' @examples
#' library(Seurat)
#' clu_pal <- c("red", "green", "blue")
#'
#' pbmc_small <- Store_Misc_Info_Seurat(seurat_object = pbmc_small, data_to_store = clu_pal,
#' data_name = "rd1_colors")
#'

Store_Palette_Seurat <- function(
  seurat_object,
  palette,
  palette_name,
  list_as_list = FALSE,
  overwrite = FALSE,
  verbose = TRUE
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  seurat_object <- Store_Misc_Info_Seurat(seurat_object = seurat_object, data_to_store = palette, data_name = palette_name, list_as_list = list_as_list, overwrite = overwrite, verbose = verbose)
  return(seurat_object)
}


#' Add Alternative Feature IDs
#'
#' Add alternative feature ids data.frame to the misc slot of Seurat object.
#'
#' @param seurat_object object name.
#' @param features_tsv_file output file from Cell Ranger used for creation of Seurat object.
#' (Either provide this of `hdf5_file`)
#' @param hdf5_file output file from Cell Ranger used for creation of Seurat object.
#' (Either provide this of `features_tsv_file`)
#' @param assay name of assay(s) to add the alternative features to.  Can specify "all"
#' to add to all assays.
#' @param data_name name to use for data.frame when stored in `@misc` slot.
#' @param overwrite logical, whether to overwrite item with the same `data_name` in the
#' `@misc` slot of object (default is FALSE).
#'
#' @import cli
#' @importFrom dplyr filter
#'
#' @return Seurat Object with new entries in the `obj@misc` slot.
#'
#' @export
#'
#' @concept get_set_util
#'
#' @examples
#' \dontrun{
#' # Using features.tsv.gz file
#'    # Either file from filtered or raw outputs can be used as they are identical.
#' obj <- Add_Alt_Feature_ID(seurat_object = obj,
#' features_tsv = "sample01/outs/filtered_feature_bc_matrix/features.tsv.gz", assay = "RNA")
#'
#' #' # Using hdf5 file
#'    # Either filtered_feature_bc or raw_feature_bc can be used as the features slot is identical
#'    # Though it is faster to load filtered_feature_bc file due to droplet filtering
#' obj <- Add_Alt_Feature_ID(seurat_object = obj,
#' hdf5_file = "sample01/outs/outs/filtered_feature_bc_matrix.h5", assay = "RNA")
#'}
#'

Add_Alt_Feature_ID <- function(
    seurat_object,
    features_tsv_file = NULL,
    hdf5_file = NULL,
    assay = NULL,
    data_name = "feature_id_mapping_table",
    overwrite = FALSE
) {
  # check file
  if (is.null(x = features_tsv_file) && is.null(x = hdf5_file)) {
    cli_abort(message = "Either {.code features_tsv_file} or {.code hdf5_file} must be provided.")
  }

  if (!is.null(x = features_tsv_file) && !is.null(x = hdf5_file)) {
    cli_abort(message = "Both {.code features_tsv_file} and {.code hdf5_file} provided.  Please only supply one or the other parameter.")
  }

  # check assay
  if (is.null(x = assay)) {
    cli_abort(message = c("Must provide value to {.code assay} to add alternative featutres to assay meta.data",
                          "i" = "Value can either be name of assay or {.val all} to add to all compatible assays present."))
  }

  # set assays to use
  if (assay == "all") {
    assays_use <- Assays(object = seurat_object)
  } else {
    assays_use <- assay
  }

  # get features
  object_features <- Features(x = seurat_object, assay = assays_use[1])

  # if providing features_tsv
  if (!is.null(x = features_tsv_file)) {
    features_table <- data.table::fread(file = features_tsv_file, header = FALSE, data.table = FALSE)
    colnames(features_table) <- c("Ensembl_ID", "Symbol", "Modality")

    features_table$Symbol <- make.unique(features_table$Symbol)

    features_present <- features_table %>%
      filter(.data[["Symbol"]] %in% object_features)
  }

  if (!is.null(x = hdf5_file)) {
    h5 <- Read10X_h5(filename = hdf5_file)
    symbols <- rownames(x = h5)

    h5 <- Read10X_h5(filename = hdf5_file, use.names = F)
    ensembl <- rownames(x = h5)

    features_table <- data.frame("Ensembl_ID" = ensembl,
                                 "Symbol" = symbols)

    features_present <- features_table %>%
      filter(.data[["Symbol"]] %in% object_features)
  }

  # Add to object
  seurat_object <- Store_Misc_Info_Seurat(seurat_object = seurat_object, data_to_store = features_present, data_name = data_name, overwrite = overwrite)

  # return object
  return(seurat_object)
}
