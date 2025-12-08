#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### Operators ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Set a default value if an object is NOT null
#'
#' @param lhs An object to set if it's NOT null
#' @param rhs The value to provide if x is NOT null
#'
#' @return lhs if lhs is null, else rhs
#'
#' @author Hadley Wickham
#' @references \url{https://adv-r.hadley.nz/functions.html#missing-arguments}
#'
#' @noRd
#'

`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}


#' Set a default value if an object is null
#'
#' @param lhs An object to set if it's null
#' @param rhs The value to provide if x is null
#'
#' @return rhs if lhs is null, else lhs
#'
#' @author Hadley Wickham
#' @references \url{https://adv-r.hadley.nz/functions.html#missing-arguments}
#'
#' @noRd
#'

`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### Object/Feature Checks ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Check Seurat Object
#'
#' Checks if object is of class: Seurat and returns error message if not.
#'
#' @param seurat_object Seurat object name.
#'
#' @return error if not Seurat object.
#'
#' @noRd
#'

Is_Seurat <- function(
  seurat_object
) {
  if (isFALSE(x = inherits(what = "Seurat", x = seurat_object))) {
    cli_abort(message = "{.code seurat_object} provided is not an object of class: Seurat.")
  }
}


#' Check LIGER Object
#'
#' Checks if object is of class: liger and returns error message if not.
#'
#' @param liger_object liger object name.
#'
#' @return error is not LIGER object
#'
#' @noRd
#'

Is_LIGER <- function(
  liger_object
) {
  if (isFALSE(x = inherits(what = "liger", x = liger_object))) {
    cli_abort(message = "{.code liger_object} provided is not an object of class: liger")
  }
}


#' Check for assays present in object
#
#' Checks Seurat object for the presence of assays against list of specified assays
#'
#' @param seurat_object Seurat object name.
#' @param assay_list vector of genes to check.
#' @param print_msg logical. Whether message should be printed if all features are found.  Default is TRUE.
#' @param omit_warn logical. Whether to print message about features that are not found in current object.
#'
#' @return List of found vs not found assays.
#'
#' @noRd
#'

Assay_Present <- function(
  seurat_object,
  assay_list,
  print_msg = TRUE,
  omit_warn = TRUE
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # get all features
  possible_assays <- Assays(object = seurat_object)

  # If any features not found
  if (any(!assay_list %in% possible_assays)) {
    bad_assays <- assay_list[!assay_list %in% possible_assays]
    found_assays <- assay_list[assay_list %in% possible_assays]
    if (length(x = found_assays) == 0) {
      cli_abort(message = "No requested assays found.")
    }

    # Return message of assays not found
    if (length(x = bad_assays) > 0 && isTRUE(x = omit_warn)) {
      cli_warn(message = c("The following assays were omitted as they were not found:",
                           "i" = "{.field {glue_collapse_scCustom(input_string = bad_assays, and = TRUE)}}.")
      )
    }

    # Combine into list and return
    assay_list <- list(
      found_assays = found_assays,
      bad_assays = bad_assays
    )
    return(assay_list)
  }

  # Print all found message if TRUE
  if (isTRUE(x = print_msg)) {
    cli_inform(message = "All assays present.")
  }

  # Return full input gene list.
  # Combine into list and return
  assay_list <- list(
    found_assays = assay_list,
    bad_assays = NULL
  )
  return(assay_list)
}


#' Check whether assay is V5
#
#' Checks Seurat object to verify whether it is composed of "Assay" or "Assay5" slots.
#'
#' @param seurat_object Seurat object name.
#' @param assay name of assay to check, default is NULL.
#'
#' @return TRUE if seurat_object contains "Assay5" class.
#'
#' @noRd
#'

Assay5_Check <- function(
  seurat_object,
  assay = NULL
) {
  assay <- assay %||% DefaultAssay(object = seurat_object)

  if (inherits(x = seurat_object@assays[[assay]], what = "Assay")) {
    return(FALSE)
  }
  if (inherits(x = seurat_object@assays[[assay]], what = "Assay5")) {
    return(TRUE)
  }
}


#' Perform Feature and Meta Checks before plotting
#'
#' Wraps the `Feature_Present`, `Meta_Present`, `Reduction_Loading_Present`, and `Case_Check` into
#' single function to perform feature checks before plotting.
#'
#' @param object Seurat object
#' @param features vector of features and/or meta data variables to plot.
#' @param assay Assay to use (default all assays present).
#' @param print_missing logical, whether to print the names of missing features.  Default is TRUE.
#'
#' @return vector of features and/or meta data that were found in object.
#'
#' @noRd
#'
#' @keywords internal
#'

Feature_PreCheck <- function(
  object,
  features,
  assay = NULL,
  print_missing = TRUE
) {
  # set assay (if null set to active assay)
  assay <- assay %||% Assays(object = object)

  # Check features and meta to determine which features present
  features_list <- Feature_Present(data = object, features = features, omit_warn = FALSE, print_msg = FALSE, case_check_msg = FALSE, return_none = TRUE, seurat_assay = assay)

  meta_list <- Meta_Present(object = object, meta_col_names = features_list[[2]], omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)

  reduction_list <- Reduction_Loading_Present(seurat_object = object, reduction_names = meta_list[[2]], omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)

  all_not_found_features <- reduction_list[[2]]

  all_found_features <- c(features_list[[1]], meta_list[[1]], reduction_list[[1]])

  # Stop if no features found
  if (length(x = all_found_features) < 1) {
    cli_abort(message = c("No features were found.",
                          "*" = "The following are not present in object:",
                          "i" = "{.field {glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}}")
    )
  }

  # Return message of features not found
  if (isTRUE(x = print_missing)) {
    if (length(x = all_not_found_features) > 0) {
      op <- options(warn = 1)
      on.exit(options(op))
      cli_warn(message = c("The following features were omitted as they were not found:",
                           "i" = "{.field {glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}}")
      )
    }
  }

  # Check feature case and message if found
  Case_Check(seurat_object = object, gene_list = all_not_found_features, case_check_msg = TRUE, return_features = FALSE)

  # return all found features
  return(all_found_features)
}


#' Check for Normalized Data in Seurat object
#'
#' Runs a check to determine if all values in the "data" layer of object are whole numbers.  If so
#' functions returns error message to first normalize data.
#'
#' @param object Seurat object.
#' @param assay Name of assay to check for normalized values within.
#' @param error logical, whether the function errors if no normalized data is found (default is TRUE), or
#' if value of the check is returned.
#'
#' @returns either error message, nothing, or logical value depending on `error` parameter.
#' @noRd
#'
#' @examples
#' \dontrun{
#' Check_Normalized(object = obj)
#'}
#'

Check_Normalized <- function(
    object,
    assay = NULL,
    error = TRUE
) {
  # set assay (if null set to active assay)
  assay <- assay %||% Assays(object = object)

  # First check if data present
  if (isTRUE(x = Assay5_Check(seurat_object = object, assay = assay))) {
    if (isTRUE(x = length(grep(x = Layers(object = object), pattern = "data", value = T)) == 0)) {
      value <- FALSE
      if (isTRUE(x = error)) {
        cli_abort(message = c("Layer with normalized data not present.",
                              "i" = "Please Normalize data first."))
      } else {
        return(value)
      }
    } else {
      value <- TRUE
      return(value)
    }
  }

  # Pull data to check (limit to first 100 cells for speed and memory considerations)
  if (length(x = Cells(x = object)) > 100) {
    data <- suppressWarnings(LayerData(object = object, layer = "data", assay = "RNA")[,1:100])
  } else {
    data <- suppressWarnings(LayerData(object = object, layer = "data", assay = "RNA"))
  }

  # Check normalized by seeing if all values are whole numbers
  value <- !all(floor(x = data) == data)

  # check if erroring or returning value
  if (isTRUE(x = error)) {
    if (isFALSE(x = value)) {
      cli_abort(message = c("The data is not normalized.",
                            "*" = "The {symbol$dquote_left}data{symbol$dquote_right} layer contains all whole numbers.",
                            "i" = "Please Normalize data first."))
    } else {
      return(value)
    }
  } else {
    return(value)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### FUNCTION HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Stop function without error message
#'
#' Modifies R options within the function call only to hide the error message from `stop` while
#' keeping global options preserved outside of function.
#'
#' @return stops function without error message
#'
#' @author Stibu
#' @references \url{https://stackoverflow.com/a/42945293/15568251}
#' @details \url{https://creativecommons.org/licenses/by-sa/3.0/}
#'
#' @noRd
#'

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}


#' Custom glue collapse
#
#' Customized glue_collapse that is based on the number of items in input list
#'
#' @param input_string input string to be collapsed.
#' @param and logical.  Whether to use "and" or "or" in the collapsed string
#'
#' @return collapsed string
#'
#' @importFrom glue glue_collapse
#'
#' @noRd
#'

glue_collapse_scCustom <- function(
  input_string,
  and = TRUE
) {
  # Check length of input string
  input_length <- length(x = input_string)

  # set last seperator
  if (isTRUE(x = and)) {
    last_sep <- " and "
  } else {
    last_sep <- " or "
  }

  if (input_length <= 3) {
    glue_collapse(x = input_string, sep = ", ", last = last_sep)
  } else {
    glue_collapse(x = input_string, sep = ", ", last = paste0(",", last_sep))
  }
}


#' Ask yes/no question to proceed
#'
#' Asks the user to answer yes/no question and returns logical value depending on
#' the answer.
#'
#' @return logical
#'
#' @references function modified from function in devtools R package (License: MIT) \url{https://github.com/r-lib/devtools}.
#' @details \url{https://github.com/r-lib/devtools/blob/9f27cc3e6335e74d6f51ed331509ebda56747901/R/release.R#L147-L156}.
#'
#' @noRd
#'

yesno <- function(msg, .envir = parent.frame()) {
  cli_inform(message = msg, .envir = .envir)
  qs <- c("Yes", "No")
  rand <- c(1,2)

  utils::menu(qs[1:2]) != which(rand == 1)
}


#' Change function parameter value from NULL to NA
#'
#' Provides method to change parameter value dynamically within function to suit defaults of other functions.
#' Used in iterative plotting functions.
#'
#' @param parameter the parameter to check for NULL
#'
#' @return if NULL returns NA otherwise returns input value.
#'
#' @noRd
#'

replace_null <- function(
  parameter
) {
  # check length
  if (length(x = parameter) > 1) {
    cli_abort(message = "{.code parameter} must be single value.")
  }

  # check NULL and swap NA
  if (is.null(x = parameter)) {
    parameter <- NA
  }
  return(parameter)
}


#' Check if file name has proper extension at the end
#'
#' Checks if file name has the correct extension at the end and returns logical value
#'
#' @param file_name string containing file name to check
#' @param extension the extension to check
#'
#' @return TRUE if extension is present otherwise FALSE
#'
#' @noRd
#'

check_extension <- function(
    file_name,
    extension
) {
  # check extension
  file_ext <- grep(x = file_name, pattern = paste0(extension, "$"))

  res <- ifelse(
    test = length(x = file_ext) == 0,
    yes = FALSE,
    no = TRUE
  )
  return(res)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GENERAL HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Calculate the percentage of a vector above some threshold
#'
#' @param x Vector of values
#' @param threshold Threshold to use when calculating percentage
#'
#' @return Returns the percentage of `x` values above the given threshold
#'
#' @author Satija Lab & all co-Authors of Seurat Package
#' @references See Utilities.R in source code of Seurat \url{https://github.com/satijalab/seurat/blob/master/R/utilities.R}  (Licence: GPL-3).
#'
#' @note modified from Seurat version to return a percentage instead of proportion/decimal as part of `Percent_Expressing` function.
#'
#' @keywords internal
#'
#' @noRd
#'

PercentAbove_Seurat <- function(
  x,
  threshold
) {
  return((sum(x > threshold, na.rm = T) / length(x = x)) * 100)
}


#' Calculate percent of expressing cells from meta data category
#'
#' For internal use in calculating module score significance
#'
#' @param seurat_object Seurat object name.
#' @param features Feature(s) to plot.
#' @param threshold Expression threshold to use for calculation of percent expressing (default is 0).
#' @param group.by Factor to group the cells by.
#' @param split.by Factor to split the groups by.
#' @param entire_object logical (default = FALSE).  Whether to calculate percent of expressing cells
#' across the entire object as opposed to by cluster or by `group.by` variable.
#' @param assay Assay to pull feature data from.  Default is active assay.
#' @param layer Which layer to pull expression data from?  Default is "data".
#'
#' @return A data.frame
#'
#' @references Part of code is modified from Seurat package as used by \code{\link[Seurat]{DotPlot}}
#' to generate values to use for plotting.  Source code can be found here:
#' \url{https://github.com/satijalab/seurat/blob/4e868fcde49dc0a3df47f94f5fb54a421bfdf7bc/R/visualization.R#L3391} (License: GPL-3).
#'
#' @keywords internal
#'
#' @noRd

Percent_Expressing_Meta <- function(
  seurat_object,
  features,
  threshold = 0,
  group.by = NULL,
  split.by = NULL,
  entire_object = FALSE,
  layer = "data",
  assay = NULL
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # set assay (if null set to active assay)
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Check features exist in object
  features_list <- Feature_PreCheck(object = seurat_object, features = features, assay = assay)

  # Check group.by is in object
  if (!is.null(x = group.by) && group.by == "ident") {
    group.by <- NULL
  }

  if (!is.null(x = group.by)) {
    possible_groups <- colnames(x = seurat_object@meta.data)
    if (!group.by %in% possible_groups) {
      cli_abort("Grouping variable {.val {group.by}} was not found in Seurat Object.")
    }
  }

  # Check split.by is in object
  if (!is.null(x = split.by)) {
    possible_groups <- colnames(x = seurat_object@meta.data)
    if (!split.by %in% possible_groups) {
      cli_abort("Splitting variable {.val {split.by}} was not found in Seurat Object.")
    }
  }

  # Pull Expression Info
  cells <- unlist(x = CellsByIdentities(object = seurat_object, idents = NULL))
  expression_info <- FetchData(object = seurat_object, vars = features_list, cells = cells, layer = layer)

  # Add grouping variable
  if (isTRUE(x = entire_object)) {
    expression_info$id <- "All_Cells"
  } else {
    expression_info$id <- if (is.null(x = group.by)) {
      Idents(object = seurat_object)[cells, drop = TRUE]
    } else {
      seurat_object[[group.by, drop = TRUE]][cells, drop = TRUE]
    }
  }
  if (!is.factor(x = expression_info$id)) {
    expression_info$id <- factor(x = expression_info$id)
  }
  id.levels <- levels(x = expression_info$id)
  expression_info$id <- as.vector(x = expression_info$id)

  # Split data if split.by is true
  if (!is.null(x = split.by)) {
    splits <- seurat_object[[split.by, drop = TRUE]][cells, drop = TRUE]
    expression_info$id <- paste(expression_info$id, splits, sep = "_")
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


#' Extract delimiter information from a string.
#'
#' Parses a string (usually a cell name) and extracts fields based on a delimiter
#'
#' @param string String to parse.
#' @param field Integer(s) indicating which field(s) to extract. Can be a vector multiple numbers.
#' @param delim Delimiter to use, set to underscore by default.
#'
#' @return A new string, that parses out the requested fields, and (if multiple), rejoins them with the same delimiter
#'
#' @references See Utilities.R in source code of Seurat \url{https://github.com/satijalab/seurat/blob/master/R/utilities.R}  (Licence: GPL-3).
#'
#' @keywords internal
#'
#' @noRd
#'

ExtractField <- function(
  string,
  field = 1,
  delim = "_"
) {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}


#' Calculate the middle value between two numbers
#'
#' @param min Lower value.
#' @param max Higher value.
#'
#' @return Returns number in middle of two provided values.
#'
#' @references Code to calculate middle value from: adapted from: \url{https://stackoverflow.com/a/54147509/15568251}.
#' Renamed and wrapped into function by Samuel Marsh.
#'
#' @keywords internal
#'
#' @noRd
#'

Middle_Number <- function(
  min,
  max
) {
  min_max <- c(min, max)
  middle <- min_max[-length(min_max)] + diff(min_max) / 2
  return(middle)
}


#' Symmetrical setdiff
#'
#' tests for differences between two vectors symmetrically.
#'
#' @param x first vector to test
#' @param y second vector to test
#'
#' @return vector differences x vs. y and y vs. x
#'
#' @references Function name and code from R-bloggers post:
#' \url{https://www.r-bloggers.com/2013/06/symmetric-set-differences-in-r/}
#'
#' @keywords internal
#'
#' @noRd
#'

symdiff <- function(
  x,
  y
) {
  setdiff(x = union(x = x, y = y), intersect(x = x, y = y))
}


#' Whole number check
#'
#' Checks whether a number is whole
#'
#' @param x number to check
#'
#' @return NULL or error message if number is not whole
#'
#' @keywords internal
#'
#' @noRd
#'

check_whole_num <- function(
  x
) {
  round_x <- round(x = x)

  res <-  identical(x = x, y = round_x)

  return(res)
}


#' Remove single value columns of data.frame
#'
#' Checks all columns within data.frame and returns data.frame minus columns that have the same value in all rows.
#'
#' @param df data.frame to filter
#'
#' @references Code used in function has been slightly modified from `sceasy:::.regularise_df` function of
#' sceasy package \url{https://github.com/cellgeni/sceasy} (License: GPL-3).
#' Code modified to match scCustomize & tidyverse style, add error checks, and
#' add cli formatted messages.
#'
#' @importFrom dplyr select all_of
#' @importFrom magrittr "%>%"
#'
#' @return data.frame
#'
#' @keywords internal
#'
#' @noRd
#'

drop_single_value_cols <- function(
  df
) {
  if (!inherits(what = "data.frame", x = df)) {
    cli_abort(message = "{.code df} must of be of class data.frame.")
  }

  single_val_columns <- sapply(df, function(x) {
    length(x = unique(x = x)) == 1
  })

  col_names_single <- df %>%
    select(which(single_val_columns)) %>%
    colnames()

  if (length(x = col_names_single) > 0) {
    cli_inform(message = c("The following columns were removed as they contain identical values for all rows:",
                           "i" = "{.field {col_names_single}}"))
  }

  # filter df
  df_filtered <- df %>%
    select(-(all_of(col_names_single)))

  # return df
  return(df_filtered)
}


#' Check valid color
#'
#' Checks if input values are valid colors representations in R.
#'
#' @param colors vector of color(s) to check
#'
#' @references Code for function \url{https://stackoverflow.com/a/13290832/15568251}.
#' Renamed by Samuel Marsh.
#'
#' @importFrom grDevices col2rgb
#'
#' @return logical named vector
#'
#' @keywords internal
#'
#' @noRd
#'

Is_Color <- function(
  colors
) {
  sapply(colors, function(X) {
    tryCatch(is.matrix(x = col2rgb(col = X)),
             error = function(e) FALSE)
  })
}


#' Get current seed when needed in `RunLeiden`
#'
#' @param seed current seed
#'
#' @references from BPCells package Benjamin Parks
#' https://github.com/bnprks/BPCells/blob/f2026c3509f5e2542f7624bdaf75669d5d45d78b/r/R/utils.R#L10
#'
#' @return current seed if set
#'
#' @keywords internal
#'
#' @noRd

get_seed <- function() {
  if (exists(".Random.seed", globalenv(), mode = "integer", inherits = FALSE)) {
    return(get(".Random.seed", globalenv(), mode = "integer", inherits = FALSE))
  } else {
    return(NULL)
  }
}


#' Restore previous seed when needed in `RunLeiden`
#'
#' @param seed old seed
#
#' @references from BPCells package Benjamin Parks
#' https://github.com/bnprks/BPCells/blob/f2026c3509f5e2542f7624bdaf75669d5d45d78b/r/R/utils.R#L18
#'
#' @return old seed
#'
#' @keywords internal
#'
#' @noRd

restore_seed <- function(seed) {
  if (is.null(seed)) {
    rm(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    assign(".Random.seed", seed, envir = globalenv(), inherits = FALSE)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### METRICS HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Read Gene Expression Statistics from 10X Cell Ranger Count
#'
#' Get data.frame with all metrics from the Cell Ranger `count` analysis (present in web_summary.html)
#'
#' @param base_path path to the parent directory which contains all of the subdirectories of interest.
#' @param secondary_path path from the parent directory to count "outs/" folder which contains the
#' "metrics_summary.csv" file.
#' @param default_10X logical (default TRUE) sets the secondary path variable to the default 10X directory structure.
#' @param lib_list a list of sample names (matching directory names) to import.  If `NULL` will read
#' in all samples in parent directory.
#' @param lib_names a set of sample names to use for each sample.  If `NULL` will set names to the
#' directory name of each sample.
#'
#' @return A data frame with sample metrics produced by Cell Ranger `count` pipeline.
#'
#' @import pbapply
#' @importFrom dplyr bind_rows setdiff
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' count_metrics <- Metrics_Count_GEX(base_path = base_path, lib_list = lib_list, secondary_path = secondary_path, lib_names = lib_names)
#' }
#'

Metrics_Count_GEX <- function(
  lib_list,
  base_path,
  secondary_path,
  lib_names
) {
  cli_inform(message = "Reading {.field Gene Expression} Metrics")
  raw_data_list <- pblapply(1:length(x = lib_list), function(x) {
    if (is.null(x = secondary_path)) {
      file_path <- file.path(base_path, lib_list[x])
    } else {
      file_path <- file.path(base_path, lib_list[x], secondary_path)
    }

    raw_data <- read.csv(file = file.path(file_path, "metrics_summary.csv"), stringsAsFactors = FALSE)
    # Change format of numeric columns to due commas in data csv output.
    column_numbers <- grep(pattern = ",", x = raw_data[1, ])
    raw_data[, c(column_numbers)] <- lapply(raw_data[, c(column_numbers)], function(x) {
      as.numeric(gsub(",", "", x))
    })

    column_numbers_pct <- grep(pattern = "%", x = raw_data[1, ])
    all_columns <- 1:ncol(x = raw_data)

    column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

    raw_data[, c(column_numbers_numeric)] <- lapply(raw_data[, c(column_numbers_numeric)], function(x) {
      as.numeric(x)
    })

    return(raw_data)
  })

  # Name the list items
  if (is.null(x = lib_names)) {
    names(x = raw_data_list) <- lib_list
  } else {
    names(x = raw_data_list) <- lib_names
  }

  # Combine the list and add sample_id column
  full_data <- bind_rows(raw_data_list, .id = "sample_id")

  # Change column nams to use "_" separator instead of "." for readability
  colnames(x = full_data) <- gsub(pattern = "\\.", replacement = "_", x = colnames(x = full_data))

  rownames(x = full_data) <- full_data$sample_id

  return(full_data)
}


#' Read Gene Expression Statistics from 10X Cell Ranger (v9+) Count
#'
#' Get data.frame with all metrics from the Cell Ranger (v9+) `count` analysis (present in web_summary.html)
#'
#' @param base_path path to the parent directory which contains all of the subdirectories of interest.
#' @param secondary_path path from the parent directory to count "outs/" folder which contains the
#' "metrics_summary.csv" file.
#' @param default_10X logical (default TRUE) sets the secondary path variable to the default 10X directory structure.
#' @param lib_list a list of sample names (matching directory names) to import.  If `NULL` will read
#' in all samples in parent directory.
#' @param lib_names a set of sample names to use for each sample.  If `NULL` will set names to the
#' directory name of each sample.
#'
#' @return A data frame with sample metrics produced by Cell Ranger `count` pipeline.
#'
#' @import pbapply
#' @importFrom dplyr bind_rows setdiff filter select all_of rename
#' @importFrom tibble column_to_rownames
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' count_metrics <- Metrics_Count_GEX_v9plus(base_path = base_path, lib_list = lib_list, secondary_path = secondary_path, lib_names = lib_names)
#' }
#'

Metrics_Count_GEX_v9plus <- function(
  lib_list,
  base_path,
  secondary_path,
  lib_names
) {
  cli_inform(message = "Reading {.field Gene Expression} Metrics")
  raw_data_list <- pblapply(1:length(x = lib_list), function(x) {
    if (is.null(x = secondary_path)) {
      file_path <- file.path(base_path, lib_list[x])
    } else {
      file_path <- file.path(base_path, lib_list[x], secondary_path)
    }

    raw_data <- read.csv(file = file.path(file_path, "metrics_summary.csv"), stringsAsFactors = FALSE)

    # Change format to column based and select relevant metrics
    GEX_metrics <- raw_data %>%
      filter(.data[["Grouped.By"]] == "Physical library ID" & .data[["Library.Type"]] == "Gene Expression") %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    GEX_metrics2 <- raw_data %>%
      filter(.data[["Metric.Name"]] %in% c(c("Median UMI counts per cell",
                                             "Median genes per cell",
                                             "Median reads per cell",
                                             "Total genes detected"))) %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    raw_data_gex <- cbind(GEX_metrics, GEX_metrics2)

    # Change format of numeric columns to due commas in data csv output.
    column_numbers <- grep(pattern = ",", x = raw_data_gex[1, ])
    raw_data_gex[, c(column_numbers)] <- lapply(raw_data_gex[, c(column_numbers)], function(x) {
      as.numeric(gsub(",", "", x))
    })

    if ("Estimated.number.of.cells" %in% colnames(x = raw_data_gex)) {
      # Rename multi columns to match names from count
      names_to_replace <- c(Reads.Mapped.to.Genome = "Mapped.to.genome",
                            Reads.Mapped.Confidently.to.Genome = "Confidently.mapped.to.genome",
                            Reads.Mapped.Confidently.to.Intergenic.Regions = "Confidently.mapped.to.intergenic.regions",
                            Reads.Mapped.Confidently.to.Intronic.Regions = "Confidently.mapped.to.intronic.regions",
                            Reads.Mapped.Confidently.to.Exonic.Regions = "Confidently.mapped.to.exonic.regions",
                            Reads.Mapped.Confidently.to.Transcriptome = "Confidently.mapped.to.transcriptome",
                            Reads.Mapped.Antisense.to.Gene = "Confidently.mapped.antisense",
                            Fraction.Reads.in.Cells = "Confidently.mapped.reads.in.cells",
                            Estimated.Number.of.Cells = "Estimated.number.of.cells",
                            Mean.Reads.per.Cell = "Mean.reads.per.cell",
                            Median.Genes.per.Cell = "Median.genes.per.cell",
                            Number.of.Reads = "Number.of.reads",
                            Valid.Barcodes = "Valid.barcodes",
                            Sequencing.Saturation = "Sequencing.saturation",
                            Total.Genes.Detected = "Total.genes.detected",
                            Median.UMI.Counts.per.Cell = "Median.UMI.counts.per.cell")
    } else {
      # Rename multi columns to match names from count
      names_to_replace <- c(Reads.Mapped.to.Genome = "Mapped.to.genome",
                            Reads.Mapped.Confidently.to.Genome = "Confidently.mapped.to.genome",
                            Reads.Mapped.Confidently.to.Intergenic.Regions = "Confidently.mapped.to.intergenic.regions",
                            Reads.Mapped.Confidently.to.Intronic.Regions = "Confidently.mapped.to.intronic.regions",
                            Reads.Mapped.Confidently.to.Exonic.Regions = "Confidently.mapped.to.exonic.regions",
                            Reads.Mapped.Confidently.to.Transcriptome = "Confidently.mapped.to.transcriptome",
                            Reads.Mapped.Antisense.to.Gene = "Confidently.mapped.antisense",
                            Fraction.Reads.in.Cells = "Confidently.mapped.reads.in.cells",
                            Estimated.Number.of.Cells = "Cells",
                            Mean.Reads.per.Cell = "Mean.reads.per.cell",
                            Median.Genes.per.Cell = "Median.genes.per.cell",
                            Number.of.Reads = "Number.of.reads",
                            Valid.Barcodes = "Valid.barcodes",
                            Sequencing.Saturation = "Sequencing.saturation",
                            Total.Genes.Detected = "Total.genes.detected",
                            Median.UMI.Counts.per.Cell = "Median.UMI.counts.per.cell")
    }

    raw_data_gex <- raw_data_gex %>%
      rename(all_of(names_to_replace))

    column_numbers_pct <- grep(pattern = "%", x = raw_data_gex[1, ])
    all_columns <- 1:ncol(x = raw_data_gex)

    column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

    raw_data_gex[, c(column_numbers_numeric)] <- lapply(raw_data_gex[, c(column_numbers_numeric)], function(x) {
      as.numeric(x = x)})

    return(raw_data_gex)
  })

  # Name the list items
  if (is.null(x = lib_names)) {
    names(x = raw_data_list) <- lib_list
  } else {
    names(x = raw_data_list) <- lib_names
  }

  # Combine the list and add sample_id column
  full_data <- bind_rows(raw_data_list, .id = "sample_id")

  # Change column nams to use "_" separator instead of "." for readability
  colnames(x = full_data) <- gsub(pattern = "\\.", replacement = "_", x = colnames(x = full_data))

  rownames(x = full_data) <- full_data$sample_id

  return(full_data)
}


#' Read Gene Expression Statistics from 10X Cell Ranger multi
#'
#' Get data.frame with all gene expression metrics from the Cell Ranger `multi` analysis (present in web_summary.html)
#'
#' @param base_path path to the parent directory which contains all of the subdirectories of interest.
#' @param secondary_path path from the parent directory to count "outs/" folder which contains the
#' "metrics_summary.csv" file.
#' @param default_10X logical (default TRUE) sets the secondary path variable to the default 10X directory structure.
#' @param lib_list a list of sample names (matching directory names) to import.  If `NULL` will read
#' in all samples in parent directory.
#' @param lib_names a set of sample names to use for each sample.  If `NULL` will set names to the
#' directory name of each sample.
#'
#' @return A data frame with sample gene expression metrics produced by Cell Ranger `multi` pipeline.
#'
#' @import pbapply
#' @importFrom dplyr all_of bind_rows filter rename select setdiff
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' count_multi_metrics <- Metrics_Multi_GEX(base_path = base_path, lib_list = lib_list, secondary_path = secondary_path, lib_names = lib_names)
#' }
#'

Metrics_Multi_GEX <- function(
  lib_list,
  base_path,
  secondary_path,
  lib_names
) {
  cli_inform(message = "Reading {.field Gene Expression} Metrics")

  raw_data_list <- pblapply(1:length(x = lib_list), function(x) {
    if (is.null(x = secondary_path)) {
      file_path <- file.path(base_path, lib_list[x])
    } else {
      file_path <- file.path(base_path, lib_list[x], secondary_path, lib_list[x])
    }

    raw_data <- read.csv(file = file.path(file_path, "metrics_summary.csv"), stringsAsFactors = FALSE)

    # Change format to column based and select relevant metrics
    GEX_metrics <- raw_data %>%
      filter(.data[["Grouped.By"]] == "Physical library ID" & .data[["Library.Type"]] == "Gene Expression") %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    GEX_metrics2 <- raw_data %>%
      filter(.data[["Metric.Name"]] %in% c("Median UMI counts per cell",
                                           "Median genes per cell",
                                           "Median reads per cell",
                                           "Total genes detected")
             & .data[["Grouped.By"]]== "" & .data[["Library.Type"]] == "Gene Expression") %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    raw_data_gex <- cbind(GEX_metrics, GEX_metrics2)

    # Change format of numeric columns to due commas in data csv output.
    column_numbers <- grep(pattern = ",", x = raw_data_gex[1, ])
    raw_data_gex[, c(column_numbers)] <- lapply(raw_data_gex[, c(column_numbers)], function(x) {
      as.numeric(gsub(",", "", x))})

    # Rename multi columns to match names from count
    if ("Estimated.number.of.cells" %in% colnames(x = raw_data_gex)) {
      # Rename multi columns to match names from count
      names_to_replace <- c(Reads.Mapped.to.Genome = "Mapped.to.genome",
                            Reads.Mapped.Confidently.to.Genome = "Confidently.mapped.to.genome",
                            Reads.Mapped.Confidently.to.Intergenic.Regions = "Confidently.mapped.to.intergenic.regions",
                            Reads.Mapped.Confidently.to.Intronic.Regions = "Confidently.mapped.to.intronic.regions",
                            Reads.Mapped.Confidently.to.Exonic.Regions = "Confidently.mapped.to.exonic.regions",
                            Reads.Mapped.Confidently.to.Transcriptome = "Confidently.mapped.to.transcriptome",
                            Reads.Mapped.Antisense.to.Gene = "Confidently.mapped.antisense",
                            Fraction.Reads.in.Cells = "Confidently.mapped.reads.in.cells",
                            Estimated.Number.of.Cells = "Estimated.number.of.cells",
                            Mean.Reads.per.Cell = "Mean.reads.per.cell",
                            Median.Genes.per.Cell = "Median.genes.per.cell",
                            Number.of.Reads = "Number.of.reads",
                            Valid.Barcodes = "Valid.barcodes",
                            Sequencing.Saturation = "Sequencing.saturation",
                            Total.Genes.Detected = "Total.genes.detected",
                            Median.UMI.Counts.per.Cell = "Median.UMI.counts.per.cell")
    } else {
      # Rename multi columns to match names from count
      names_to_replace <- c(Reads.Mapped.to.Genome = "Mapped.to.genome",
                            Reads.Mapped.Confidently.to.Genome = "Confidently.mapped.to.genome",
                            Reads.Mapped.Confidently.to.Intergenic.Regions = "Confidently.mapped.to.intergenic.regions",
                            Reads.Mapped.Confidently.to.Intronic.Regions = "Confidently.mapped.to.intronic.regions",
                            Reads.Mapped.Confidently.to.Exonic.Regions = "Confidently.mapped.to.exonic.regions",
                            Reads.Mapped.Confidently.to.Transcriptome = "Confidently.mapped.to.transcriptome",
                            Reads.Mapped.Antisense.to.Gene = "Confidently.mapped.antisense",
                            Fraction.Reads.in.Cells = "Confidently.mapped.reads.in.cells",
                            Estimated.Number.of.Cells = "Cells",
                            Mean.Reads.per.Cell = "Mean.reads.per.cell",
                            Median.Genes.per.Cell = "Median.genes.per.cell",
                            Number.of.Reads = "Number.of.reads",
                            Valid.Barcodes = "Valid.barcodes",
                            Sequencing.Saturation = "Sequencing.saturation",
                            Total.Genes.Detected = "Total.genes.detected",
                            Median.UMI.Counts.per.Cell = "Median.UMI.counts.per.cell")
    }

    raw_data_gex <- raw_data_gex %>%
      rename(all_of(names_to_replace))

    column_numbers_pct <- grep(pattern = "%", x = raw_data_gex[1, ])
    all_columns <- 1:ncol(x = raw_data_gex)

    column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

    raw_data_gex[, c(column_numbers_numeric)] <- lapply(raw_data_gex[, c(column_numbers_numeric)], function(x) {
      as.numeric(x)})

    return(raw_data_gex)
  })

  # Name the list items
  if (is.null(x = lib_names)) {
    names(x = raw_data_list) <- lib_list
  } else {
    names(x = raw_data_list) <- lib_names
  }

  # Combine the list and add sample_id column
  full_data <- bind_rows(raw_data_list, .id = "sample_id")

  # Change column nams to use "_" separator instead of "." for readability
  colnames(x = full_data) <- gsub(pattern = "\\.", replacement = "_", x = colnames(x = full_data))

  rownames(x = full_data) <- full_data$sample_id

  return(full_data)
}


#' Read VDJ T Statistics from 10X Cell Ranger multi
#'
#' Get data.frame with all VDJ T metrics from the Cell Ranger `multi` analysis (present in web_summary.html)
#'
#' @param base_path path to the parent directory which contains all of the subdirectories of interest.
#' @param secondary_path path from the parent directory to count "outs/" folder which contains the
#' "metrics_summary.csv" file.
#' @param default_10X logical (default TRUE) sets the secondary path variable to the default 10X directory structure.
#' @param lib_list a list of sample names (matching directory names) to import.  If `NULL` will read
#' in all samples in parent directory.
#' @param lib_names a set of sample names to use for each sample.  If `NULL` will set names to the
#' directory name of each sample.
#'
#' @return A data frame with sample VDJ T metrics produced by Cell Ranger `multi` pipeline.
#'
#' @import pbapply
#' @importFrom dplyr all_of bind_rows filter select setdiff
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' vdj_multi_metrics <- Metrics_Multi_VDJT(base_path = base_path, lib_list = lib_list, secondary_path = secondary_path, lib_names = lib_names)
#' }
#'

Metrics_Multi_VDJT <- function(
  lib_list,
  base_path,
  secondary_path,
  lib_names
) {
  cli_inform(message = "Reading {.field VDJ T} Metrics")

  raw_data_list <- pblapply(1:length(x = lib_list), function(x) {
    if (is.null(x = secondary_path)) {
      file_path <- file.path(base_path, lib_list[x])
    } else {
      file_path <- file.path(base_path, lib_list[x], secondary_path, lib_list[x])
    }

    raw_data <- read.csv(file = file.path(file_path, "metrics_summary.csv"), stringsAsFactors = FALSE)

    VDJ_T_Metrics <- raw_data %>%
      filter(.data[["Grouped.By"]]== "Physical library ID" & .data[["Library.Type"]] == "VDJ T") %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name")

    current_metrics <- rownames(x = VDJ_T_Metrics)

    VDJ_T_Metrics <- VDJ_T_Metrics %>%
      t() %>%
      data.frame()

    remaining_metrics <- setdiff(x = c("Cells with productive TRA contig",
                                       "Cells with productive TRB contig",
                                       "Cells with productive V-J spanning (TRA, TRB) pair",
                                       "Cells with productive V-J spanning pair",
                                       "Median TRA UMIs per Cell",
                                       "Median TRB UMIs per Cell",
                                       "Number of cells with productive V-J spanning pair",
                                       "Paired clonotype diversity"),
                                 y = current_metrics)

    VDJ_T_Metrics2 <- raw_data %>%
      filter(.data[["Metric.Name"]] %in% remaining_metrics & .data[["Grouped.By"]]== "" & .data[["Library.Type"]] == "VDJ T"
      ) %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    raw_data_vdjt <- cbind(VDJ_T_Metrics, VDJ_T_Metrics2)

    column_numbers <- grep(pattern = ",", x = raw_data_vdjt[1, ])
    raw_data_vdjt[, c(column_numbers)] <- lapply(raw_data_vdjt[, c(column_numbers)], function(x) {
      as.numeric(gsub(",", "", x))})

    column_numbers_pct <- grep(pattern = "%", x = raw_data_vdjt[1, ])
    all_columns <- 1:ncol(x = raw_data_vdjt)

    column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

    raw_data_vdjt[,c(column_numbers_numeric)] <- lapply(raw_data_vdjt[,c(column_numbers_numeric)],function(x){as.numeric(x)})

    return(raw_data_vdjt)
  })

  # Name the list items
  if (is.null(x = lib_names)) {
    names(x = raw_data_list) <- lib_list
  } else {
    names(x = raw_data_list) <- lib_names
  }

  # Combine the list and add sample_id column
  full_data <- bind_rows(raw_data_list, .id = "sample_id")

  # Change column nams to use "_" separator instead of "." for readability
  colnames(x = full_data) <- gsub(pattern = "\\.", replacement = "_", x = colnames(x = full_data))

  rownames(x = full_data) <- full_data$sample_id

  return(full_data)
}


#' Read VDJ B Statistics from 10X Cell Ranger multi
#'
#' Get data.frame with all VDJ B metrics from the Cell Ranger `multi` analysis (present in web_summary.html)
#'
#' @param base_path path to the parent directory which contains all of the subdirectories of interest.
#' @param secondary_path path from the parent directory to count "outs/" folder which contains the
#' "metrics_summary.csv" file.
#' @param default_10X logical (default TRUE) sets the secondary path variable to the default 10X directory structure.
#' @param lib_list a list of sample names (matching directory names) to import.  If `NULL` will read
#' in all samples in parent directory.
#' @param lib_names a set of sample names to use for each sample.  If `NULL` will set names to the
#' directory name of each sample.
#'
#' @return A data frame with sample VDJ B metrics produced by Cell Ranger `multi` pipeline.
#'
#' @import pbapply
#' @importFrom dplyr all_of bind_rows filter select setdiff
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' vdj_multi_metrics <- Metrics_Multi_VDJB(base_path = base_path, lib_list = lib_list, secondary_path = secondary_path, lib_names = lib_names)
#' }
#'

Metrics_Multi_VDJB <- function(
    lib_list,
    base_path,
    secondary_path,
    lib_names
) {
  cli_inform(message = "Reading {.field VDJ B} Metrics")

  raw_data_list <- pblapply(1:length(x = lib_list), function(x) {
    if (is.null(x = secondary_path)) {
      file_path <- file.path(base_path, lib_list[x])
    } else {
      file_path <- file.path(base_path, lib_list[x], secondary_path, lib_list[x])
    }

    raw_data <- read.csv(file = file.path(file_path, "metrics_summary.csv"), stringsAsFactors = FALSE)

    VDJ_B_Metrics <- raw_data %>%
      filter(.data[["Grouped.By"]]== "Physical library ID" & .data[["Library.Type"]] == "VDJ B") %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name")

    current_metrics <- rownames(x = VDJ_B_Metrics)

    VDJ_B_Metrics <- VDJ_B_Metrics %>%
      t() %>%
      data.frame()

    remaining_metrics <- setdiff(x = c("Cells with productive IGH contig",
                                       "Cells with productive IGK contig",
                                       "Cells with productive IGL contig",
                                       "Cells with productive V-J spanning (IGK, IGH) pair",
                                       "Cells with productive V-J spanning (IGL, IGH) pair",
                                       "Cells with productive V-J spanning pair",
                                       "Median IGH UMIs per Cell",
                                       "Median IGK UMIs per Cell",
                                       "Median IGL UMIs per Cell",
                                       "Number of cells with productive V-J spanning pair",
                                       "Paired clonotype diversity"),
                                 y = current_metrics)

    VDJ_B_Metrics2 <- raw_data %>%
      filter(.data[["Metric.Name"]] %in% remaining_metrics & .data[["Grouped.By"]]== "" & .data[["Library.Type"]] == "VDJ B"
      ) %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    raw_data_vdjb <- cbind(VDJ_B_Metrics, VDJ_B_Metrics2)

    column_numbers <- grep(pattern = ",", x = raw_data_vdjb[1, ])
    raw_data_vdjb[, c(column_numbers)] <- lapply(raw_data_vdjb[, c(column_numbers)], function(x) {
      as.numeric(gsub(",", "", x))})

    column_numbers_pct <- grep(pattern = "%", x = raw_data_vdjb[1, ])
    all_columns <- 1:ncol(x = raw_data_vdjb)

    column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

    raw_data_vdjb[,c(column_numbers_numeric)] <- lapply(raw_data_vdjb[,c(column_numbers_numeric)],function(x){as.numeric(x)})

    return(raw_data_vdjb)
  })

  # Name the list items
  if (is.null(x = lib_names)) {
    names(x = raw_data_list) <- lib_list
  } else {
    names(x = raw_data_list) <- lib_names
  }

  # Combine the list and add sample_id column
  full_data <- bind_rows(raw_data_list, .id = "sample_id")

  # Change column nams to use "_" separator instead of "." for readability
  colnames(x = full_data) <- gsub(pattern = "\\.", replacement = "_", x = colnames(x = full_data))

  rownames(x = full_data) <- full_data$sample_id

  return(full_data)
}


#' Read Antibody Capture Statistics from 10X Cell Ranger multi
#'
#' Get data.frame with all Antibody Capture metrics from the Cell Ranger `multi` analysis (present in web_summary.html)
#'
#' @param base_path path to the parent directory which contains all of the subdirectories of interest.
#' @param secondary_path path from the parent directory to count "outs/" folder which contains the
#' "metrics_summary.csv" file.
#' @param default_10X logical (default TRUE) sets the secondary path variable to the default 10X directory structure.
#' @param lib_list a list of sample names (matching directory names) to import.  If `NULL` will read
#' in all samples in parent directory.
#' @param lib_names a set of sample names to use for each sample.  If `NULL` will set names to the
#' directory name of each sample.
#'
#' @return A data frame with sample Antibody Capture metrics produced by Cell Ranger `multi` pipeline.
#'
#' @import pbapply
#' @importFrom dplyr all_of bind_rows filter select setdiff
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' abc_multi_metrics <- Metrics_Multi_ABC(base_path = base_path, lib_list = lib_list, secondary_path = secondary_path, lib_names = lib_names)
#' }
#'

Metrics_Multi_ABC <- function(
    lib_list,
    base_path,
    secondary_path,
    lib_names
) {
  cli_inform(message = "Reading {.field Antibody Capture} Metrics")

  raw_data_list <- pblapply(1:length(x = lib_list), function(x) {
    if (is.null(x = secondary_path)) {
      file_path <- file.path(base_path, lib_list[x])
    } else {
      file_path <- file.path(base_path, lib_list[x], secondary_path, lib_list[x])
    }

    raw_data <- read.csv(file = file.path(file_path, "metrics_summary.csv"), stringsAsFactors = FALSE)

    ABC_Metrics <- raw_data %>%
      filter(.data[["Grouped.By"]]== "Physical library ID" & .data[["Library.Type"]] == "Antibody Capture") %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name")

    current_metrics <- rownames(x = ABC_Metrics)

    ABC_Metrics <- ABC_Metrics %>%
      t() %>%
      data.frame()

    remaining_metrics <- setdiff(x = c("Antibody reads in cells",
                                        "Fraction antibody reads",
                                        "Fraction antibody reads in aggregate barcodes",
                                        "Mean antibody reads usable per cell",
                                        "Median UMI counts per cell",
                                        "Number of reads from cells associated with this sample"),
                                  y = current_metrics)

    ABC_Metrics2 <- raw_data %>%
      filter(.data[["Metric.Name"]] %in% remaining_metrics & .data[["Grouped.By"]]== "" & .data[["Library.Type"]] == "Antibody Capture"
      ) %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    raw_data_abc <- cbind(ABC_Metrics, ABC_Metrics2)

    column_numbers <- grep(pattern = ",", x = raw_data_abc[1, ])
    raw_data_abc[, c(column_numbers)] <- lapply(raw_data_abc[, c(column_numbers)], function(x) {
      as.numeric(gsub(",", "", x))})

    column_numbers_pct <- grep(pattern = "%", x = raw_data_abc[1, ])
    all_columns <- 1:ncol(x = raw_data_abc)

    column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

    raw_data_abc[,c(column_numbers_numeric)] <- lapply(raw_data_abc[,c(column_numbers_numeric)],function(x){as.numeric(x)})

    return(raw_data_abc)
  })

  # Name the list items
  if (is.null(x = lib_names)) {
    names(x = raw_data_list) <- lib_list
  } else {
    names(x = raw_data_list) <- lib_names
  }

  # Combine the list and add sample_id column
  full_data <- bind_rows(raw_data_list, .id = "sample_id")

  # Change column nams to use "_" separator instead of "." for readability
  colnames(x = full_data) <- gsub(pattern = "\\.", replacement = "_", x = colnames(x = full_data))

  rownames(x = full_data) <- full_data$sample_id

  return(full_data)
}


#' Read single Summary Statistics csv from 10X Cell Ranger Count or Multi
#'
#' Get data.frame with all metrics from the Cell Ranger `count` analysis or list of data.frames if using Cell Ranger `multi`.
#'
#' @param base_path path to the metrics file
#' @param cellranger_multi logical, whether or not metrics come from Cell Ranger `count` or from Cell Ranger `multi`.  Default is FALSE.
#'
#' @return A data frame or list of data.frames with sample metrics from cell ranger.
#'
#' @importFrom dplyr all_of bind_rows filter rename select setdiff
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' count_metrics <- Metrics_Single_File(base_path = base_path)
#' }
#'

Metrics_Single_File <- function(
  base_path,
  cellranger_multi = FALSE
) {
  # Read GEX count metrics
  if (isFALSE(x = cellranger_multi)) {
    raw_data <- read.csv(file = base_path, stringsAsFactors = FALSE)
    # Change format of numeric columns to due commas in data csv output.
    column_numbers <- grep(pattern = ",", x = raw_data[1, ])
    raw_data[, c(column_numbers)] <- lapply(raw_data[, c(column_numbers)], function(x) {
      as.numeric(gsub(",", "", x))})


    column_numbers_pct <- grep(pattern = "%", x = raw_data[1, ])
    all_columns <- 1:ncol(x = raw_data)

    column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

    raw_data[, c(column_numbers_numeric)] <- lapply(raw_data[, c(column_numbers_numeric)], function(x) {
      as.numeric(x)})

    # Change column names to use "_" separator instead of "." for readability
    colnames(x = raw_data) <- gsub(pattern = "\\.", replacement = "_", x = colnames(x = raw_data))

    # return data
    return(raw_data)
  } else {
    # GEX metrics
    raw_data <- read.csv(file = base_path, stringsAsFactors = FALSE)

    # Change format to column based and select relevant metrics
    GEX_metrics <- raw_data %>%
      filter(.data[["Grouped.By"]] == "Physical library ID" & .data[["Library.Type"]] == "Gene Expression") %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    GEX_metrics2 <- raw_data %>%
      filter(.data[["Metric.Name"]] %in% c(c("Median UMI counts per cell",
                                             "Median genes per cell",
                                             "Median reads per cell",
                                             "Total genes detected"))) %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    raw_data_gex <- cbind(GEX_metrics, GEX_metrics2)

    # Change format of numeric columns to due commas in data csv output.
    column_numbers <- grep(pattern = ",", x = raw_data_gex[1, ])
    raw_data_gex[, c(column_numbers)] <- lapply(raw_data_gex[, c(column_numbers)], function(x) {
      as.numeric(gsub(",", "", x))})

    # Rename multi columns to match names from count
    names_to_replace <- c(Reads.Mapped.to.Genome = "Mapped.to.genome",
                          Reads.Mapped.Confidently.to.Genome = "Confidently.mapped.to.genome",
                          Reads.Mapped.Confidently.to.Intergenic.Regions = "Confidently.mapped.to.intergenic.regions",
                          Reads.Mapped.Confidently.to.Intronic.Regions = "Confidently.mapped.to.intronic.regions",
                          Reads.Mapped.Confidently.to.Exonic.Regions = "Confidently.mapped.to.exonic.regions",
                          Reads.Mapped.Confidently.to.Transcriptome = "Confidently.mapped.to.transcriptome",
                          Reads.Mapped.Antisense.to.Gene = "Confidently.mapped.antisense",
                          Fraction.Reads.in.Cells = "Confidently.mapped.reads.in.cells",
                          Estimated.Number.of.Cells = "Estimated.number.of.cells",
                          Mean.Reads.per.Cell = "Mean.reads.per.cell",
                          Median.Genes.per.Cell = "Median.genes.per.cell",
                          Number.of.Reads = "Number.of.reads",
                          Valid.Barcodes = "Valid.barcodes",
                          Sequencing.Saturation = "Sequencing.saturation",
                          Total.Genes.Detected = "Total.genes.detected",
                          Median.UMI.Counts.per.Cell = "Median.UMI.counts.per.cell")

    raw_data_gex <- raw_data_gex %>%
      rename(all_of(names_to_replace))

    column_numbers_pct <- grep(pattern = "%", x = raw_data_gex[1, ])
    all_columns <- 1:ncol(x = raw_data_gex)

    column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

    raw_data_gex[, c(column_numbers_numeric)] <- lapply(raw_data_gex[, c(column_numbers_numeric)], function(x){as.numeric(x)})

    # Change column nams to use "_" separator instead of "." for readability
    colnames(x = raw_data_gex) <- gsub(pattern = "\\.", replacement = "_", x = colnames(x = raw_data_gex))

    # Get VDJT metrics
    raw_data <- read.csv(file = base_path, stringsAsFactors = FALSE)

    VDJ_T_Metrics <- raw_data %>%
      filter(.data[["Grouped.By"]]== "Physical library ID" & .data[["Library.Type"]] == "VDJ T") %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    VDJ_T_Metrics2 <- raw_data %>%
      filter(.data[["Metric.Name"]] %in% c("Cells with productive TRA contig",
                                           "Cells with productive TRB contig",
                                           "Cells with productive V-J spanning (TRA, TRB) pair",
                                           "Cells with productive V-J spanning pair",
                                           "Median TRA UMIs per Cell",
                                           "Median TRB UMIs per Cell",
                                           "Number of cells with productive V-J spanning pair",
                                           "Paired clonotype diversity")
      ) %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    raw_data_vdjt <- cbind(VDJ_T_Metrics, VDJ_T_Metrics2)

    column_numbers <- grep(pattern = ",", x = raw_data_vdjt[1, ])
    raw_data_vdjt[, c(column_numbers)] <- lapply(raw_data_vdjt[, c(column_numbers)], function(x) {
      as.numeric(gsub(",", "", x))})

    column_numbers_pct <- grep(pattern = "%", x = raw_data_vdjt[1, ])
    all_columns <- 1:ncol(x = raw_data_vdjt)

    column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

    raw_data_vdjt[,c(column_numbers_numeric)] <- lapply(raw_data_vdjt[,c(column_numbers_numeric)],function(x){as.numeric(x)})

    # combine outputs into a list
    data_list <- list(
      multi_gex_metrics = raw_data_gex,
      multi_vdjt_metrics = raw_data_vdjt
    )

    # return data list
    return(data_list)
  }
}


#' Read single Summary Statistics csv from 10X Cell Ranger (v9+) Count or Multi
#'
#' Get data.frame with all metrics from the Cell Ranger (v9+) `count` analysis or list of data.frames if using Cell Ranger (v9+) `multi`.
#'
#' @param base_path path to the metrics file
#' @param cellranger_multi logical, whether or not metrics come from Cell Ranger `count` or from Cell Ranger `multi`.  Default is FALSE.
#'
#' @return A data frame or list of data.frames with sample metrics from cell ranger.
#'
#' @importFrom dplyr all_of bind_rows filter rename select setdiff
#' @importFrom magrittr "%>%"
#' @importFrom purrr compact
#' @importFrom tibble column_to_rownames
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' count_metrics <- Metrics_Single_File_v9plus(base_path = base_path)
#' }
#'

Metrics_Single_File_v9plus <- function(
  base_path,
  cellranger_multi = FALSE
){
  # read count metrics
  if (isFALSE(x = cellranger_multi)) {
    cli_inform(message = "Reading {.field Gene Expression} Metrics")
    raw_data <- read.csv(file = base_path, stringsAsFactors = FALSE)

    # Change format to column based and select relevant metrics
    GEX_metrics <- raw_data %>%
      filter(.data[["Grouped.By"]] == "Physical library ID" & .data[["Library.Type"]] == "Gene Expression") %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    GEX_metrics2 <- raw_data %>%
      filter(.data[["Metric.Name"]] %in% c(c("Median UMI counts per cell",
                                             "Median genes per cell",
                                             "Median reads per cell",
                                             "Total genes detected"))) %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    raw_data_gex <- cbind(GEX_metrics, GEX_metrics2)

    # Change format of numeric columns to due commas in data csv output.
    column_numbers <- grep(pattern = ",", x = raw_data_gex[1, ])
    raw_data_gex[, c(column_numbers)] <- lapply(raw_data_gex[, c(column_numbers)], function(x) {
      as.numeric(gsub(",", "", x))})

    # check how cells metric is named and adjust renaming
    if ("Estimated.number.of.cells" %in% colnames(x = raw_data_gex)) {
      # Rename multi columns to match names from count
      names_to_replace <- c(Reads.Mapped.to.Genome = "Mapped.to.genome",
                            Reads.Mapped.Confidently.to.Genome = "Confidently.mapped.to.genome",
                            Reads.Mapped.Confidently.to.Intergenic.Regions = "Confidently.mapped.to.intergenic.regions",
                            Reads.Mapped.Confidently.to.Intronic.Regions = "Confidently.mapped.to.intronic.regions",
                            Reads.Mapped.Confidently.to.Exonic.Regions = "Confidently.mapped.to.exonic.regions",
                            Reads.Mapped.Confidently.to.Transcriptome = "Confidently.mapped.to.transcriptome",
                            Reads.Mapped.Antisense.to.Gene = "Confidently.mapped.antisense",
                            Fraction.Reads.in.Cells = "Confidently.mapped.reads.in.cells",
                            Estimated.Number.of.Cells = "Estimated.number.of.cells",
                            Mean.Reads.per.Cell = "Mean.reads.per.cell",
                            Median.Genes.per.Cell = "Median.genes.per.cell",
                            Number.of.Reads = "Number.of.reads",
                            Valid.Barcodes = "Valid.barcodes",
                            Sequencing.Saturation = "Sequencing.saturation",
                            Total.Genes.Detected = "Total.genes.detected",
                            Median.UMI.Counts.per.Cell = "Median.UMI.counts.per.cell")
    } else {
      # Rename multi columns to match names from count
      names_to_replace <- c(Reads.Mapped.to.Genome = "Mapped.to.genome",
                            Reads.Mapped.Confidently.to.Genome = "Confidently.mapped.to.genome",
                            Reads.Mapped.Confidently.to.Intergenic.Regions = "Confidently.mapped.to.intergenic.regions",
                            Reads.Mapped.Confidently.to.Intronic.Regions = "Confidently.mapped.to.intronic.regions",
                            Reads.Mapped.Confidently.to.Exonic.Regions = "Confidently.mapped.to.exonic.regions",
                            Reads.Mapped.Confidently.to.Transcriptome = "Confidently.mapped.to.transcriptome",
                            Reads.Mapped.Antisense.to.Gene = "Confidently.mapped.antisense",
                            Fraction.Reads.in.Cells = "Confidently.mapped.reads.in.cells",
                            Estimated.Number.of.Cells = "Cells",
                            Mean.Reads.per.Cell = "Mean.reads.per.cell",
                            Median.Genes.per.Cell = "Median.genes.per.cell",
                            Number.of.Reads = "Number.of.reads",
                            Valid.Barcodes = "Valid.barcodes",
                            Sequencing.Saturation = "Sequencing.saturation",
                            Total.Genes.Detected = "Total.genes.detected",
                            Median.UMI.Counts.per.Cell = "Median.UMI.counts.per.cell")
    }

    raw_data_gex <- raw_data_gex %>%
      rename(all_of(names_to_replace))

    column_numbers_pct <- grep(pattern = "%", x = raw_data_gex[1, ])
    all_columns <- 1:ncol(x = raw_data_gex)

    column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

    raw_data_gex[, c(column_numbers_numeric)] <- lapply(raw_data_gex[, c(column_numbers_numeric)], function(x) {
      as.numeric(x = x)})

    # Change column nams to use "_" separator instead of "." for readability
    colnames(x = raw_data_gex) <- gsub(pattern = "\\.", replacement = "_", x = colnames(x = raw_data_gex))

    rownames(x = raw_data_gex) <- NULL

    return(raw_data_gex)
  } else {
    # read GEX metrics
    cli_inform(message = "Reading {.field Gene Expression} Metrics")
    raw_data <- read.csv(file = base_path, stringsAsFactors = FALSE)

    modalities <- unique(x = raw_data[["Library.Type"]])

    # Change format to column based and select relevant metrics
    GEX_metrics <- raw_data %>%
      filter(.data[["Grouped.By"]] == "Physical library ID" & .data[["Library.Type"]] == "Gene Expression") %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    GEX_metrics2 <- raw_data %>%
      filter(.data[["Metric.Name"]] %in% c(c("Median UMI counts per cell",
                                             "Median genes per cell",
                                             "Median reads per cell",
                                             "Total genes detected"))) %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    raw_data_gex <- cbind(GEX_metrics, GEX_metrics2)

    # Change format of numeric columns to due commas in data csv output.
    column_numbers <- grep(pattern = ",", x = raw_data_gex[1, ])
    raw_data_gex[, c(column_numbers)] <- lapply(raw_data_gex[, c(column_numbers)], function(x) {
      as.numeric(gsub(",", "", x))})

    if ("Estimated.number.of.cells" %in% colnames(x = raw_data_gex)) {
      # Rename multi columns to match names from count
      names_to_replace <- c(Reads.Mapped.to.Genome = "Mapped.to.genome",
                            Reads.Mapped.Confidently.to.Genome = "Confidently.mapped.to.genome",
                            Reads.Mapped.Confidently.to.Intergenic.Regions = "Confidently.mapped.to.intergenic.regions",
                            Reads.Mapped.Confidently.to.Intronic.Regions = "Confidently.mapped.to.intronic.regions",
                            Reads.Mapped.Confidently.to.Exonic.Regions = "Confidently.mapped.to.exonic.regions",
                            Reads.Mapped.Confidently.to.Transcriptome = "Confidently.mapped.to.transcriptome",
                            Reads.Mapped.Antisense.to.Gene = "Confidently.mapped.antisense",
                            Fraction.Reads.in.Cells = "Confidently.mapped.reads.in.cells",
                            Estimated.Number.of.Cells = "Estimated.number.of.cells",
                            Mean.Reads.per.Cell = "Mean.reads.per.cell",
                            Median.Genes.per.Cell = "Median.genes.per.cell",
                            Number.of.Reads = "Number.of.reads",
                            Valid.Barcodes = "Valid.barcodes",
                            Sequencing.Saturation = "Sequencing.saturation",
                            Total.Genes.Detected = "Total.genes.detected",
                            Median.UMI.Counts.per.Cell = "Median.UMI.counts.per.cell")
    } else {
      # Rename multi columns to match names from count
      names_to_replace <- c(Reads.Mapped.to.Genome = "Mapped.to.genome",
                            Reads.Mapped.Confidently.to.Genome = "Confidently.mapped.to.genome",
                            Reads.Mapped.Confidently.to.Intergenic.Regions = "Confidently.mapped.to.intergenic.regions",
                            Reads.Mapped.Confidently.to.Intronic.Regions = "Confidently.mapped.to.intronic.regions",
                            Reads.Mapped.Confidently.to.Exonic.Regions = "Confidently.mapped.to.exonic.regions",
                            Reads.Mapped.Confidently.to.Transcriptome = "Confidently.mapped.to.transcriptome",
                            Reads.Mapped.Antisense.to.Gene = "Confidently.mapped.antisense",
                            Fraction.Reads.in.Cells = "Confidently.mapped.reads.in.cells",
                            Estimated.Number.of.Cells = "Cells",
                            Mean.Reads.per.Cell = "Mean.reads.per.cell",
                            Median.Genes.per.Cell = "Median.genes.per.cell",
                            Number.of.Reads = "Number.of.reads",
                            Valid.Barcodes = "Valid.barcodes",
                            Sequencing.Saturation = "Sequencing.saturation",
                            Total.Genes.Detected = "Total.genes.detected",
                            Median.UMI.Counts.per.Cell = "Median.UMI.counts.per.cell")
    }

    raw_data_gex <- raw_data_gex %>%
      rename(all_of(names_to_replace))

    column_numbers_pct <- grep(pattern = "%", x = raw_data_gex[1, ])
    all_columns <- 1:ncol(x = raw_data_gex)

    column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

    raw_data_gex[, c(column_numbers_numeric)] <- lapply(raw_data_gex[, c(column_numbers_numeric)], function(x) {
      as.numeric(x)})

    # Change column nams to use "_" separator instead of "." for readability
    colnames(x = raw_data_gex) <- gsub(pattern = "\\.", replacement = "_", x = colnames(x = raw_data_gex))

    rownames(x = raw_data_gex) <- NULL

    if ("VDJ T" %in% modalities) {
      # Get VDJT metrics
      cli_inform(message = "Reading {.field VDJ T} Metrics")
      raw_data <- read.csv(file = base_path, stringsAsFactors = FALSE)

      VDJ_T_Metrics <- raw_data %>%
        filter(.data[["Grouped.By"]]== "Physical library ID" & .data[["Library.Type"]] == "VDJ T") %>%
        select(all_of(c("Metric.Name", "Metric.Value"))) %>%
        column_to_rownames("Metric.Name")

      current_metrics <- rownames(x = VDJ_T_Metrics)

      VDJ_T_Metrics <- VDJ_T_Metrics %>%
        t() %>%
        data.frame()

      remaining_metrics <- setdiff(x = c("Cells with productive TRA contig",
                                         "Cells with productive TRB contig",
                                         "Cells with productive V-J spanning (TRA, TRB) pair",
                                         "Cells with productive V-J spanning pair",
                                         "Median TRA UMIs per Cell",
                                         "Median TRB UMIs per Cell",
                                         "Number of cells with productive V-J spanning pair",
                                         "Paired clonotype diversity"),
                                   y = current_metrics)

      VDJ_T_Metrics2 <- raw_data %>%
        filter(.data[["Metric.Name"]] %in% remaining_metrics & .data[["Grouped.By"]]== "" & .data[["Library.Type"]] == "VDJ T"
        ) %>%
        select(all_of(c("Metric.Name", "Metric.Value"))) %>%
        column_to_rownames("Metric.Name") %>%
        t() %>%
        data.frame()

      raw_data_vdjt <- cbind(VDJ_T_Metrics, VDJ_T_Metrics2)

      column_numbers <- grep(pattern = ",", x = raw_data_vdjt[1, ])
      raw_data_vdjt[, c(column_numbers)] <- lapply(raw_data_vdjt[, c(column_numbers)], function(x) {
        as.numeric(gsub(",", "", x))})

      column_numbers_pct <- grep(pattern = "%", x = raw_data_vdjt[1, ])
      all_columns <- 1:ncol(x = raw_data_vdjt)

      column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

      raw_data_vdjt[, c(column_numbers_numeric)] <- lapply(raw_data_vdjt[, c(column_numbers_numeric)], function(x) {
        as.numeric(x)})
    } else {
      raw_data_vdjt <- NULL
    }

    # VDJ B Metrics
    if ("VDJ B" %in% modalities) {
      # Get VDJB metrics
      cli_inform(message = "Reading {.field VDJ B} Metrics")
      raw_data <- read.csv(file = base_path, stringsAsFactors = FALSE)

      VDJ_B_Metrics <- raw_data %>%
        filter(.data[["Grouped.By"]]== "Physical library ID" & .data[["Library.Type"]] == "VDJ B") %>%
        select(all_of(c("Metric.Name", "Metric.Value"))) %>%
        column_to_rownames("Metric.Name")

      current_metrics <- rownames(x = VDJ_T_Metrics)

      VDJ_B_Metrics <- VDJ_T_Metrics %>%
        t() %>%
        data.frame()

      remaining_metrics <- setdiff(x = c("Cells with productive IGH contig",
                                         "Cells with productive IGK contig",
                                         "Cells with productive IGL contig",
                                         "Cells with productive V-J spanning (IGK, IGH) pair",
                                         "Cells with productive V-J spanning (IGL, IGH) pair",
                                         "Cells with productive V-J spanning pair",
                                         "Median IGH UMIs per Cell",
                                         "Median IGK UMIs per Cell",
                                         "Median IGL UMIs per Cell",
                                         "Number of cells with productive V-J spanning pair",
                                         "Paired clonotype diversity"),
                                   y = current_metrics)

      VDJ_B_Metrics2 <- raw_data %>%
        filter(.data[["Metric.Name"]] %in% remaining_metrics & .data[["Grouped.By"]]== "" & .data[["Library.Type"]] == "VDJ B"
        ) %>%
        select(all_of(c("Metric.Name", "Metric.Value"))) %>%
        column_to_rownames("Metric.Name") %>%
        t() %>%
        data.frame()

      raw_data_vdjb <- cbind(VDJ_B_Metrics, VDJ_B_Metrics2)

      column_numbers <- grep(pattern = ",", x = raw_data_vdjb[1, ])
      raw_data_vdjb[, c(column_numbers)] <- lapply(raw_data_vdjb[, c(column_numbers)], function(x) {
        as.numeric(gsub(",", "", x))})

      column_numbers_pct <- grep(pattern = "%", x = raw_data_vdjb[1, ])
      all_columns <- 1:ncol(x = raw_data_vdjb)

      column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

      raw_data_vdjb[,c(column_numbers_numeric)] <- lapply(raw_data_vdjb[,c(column_numbers_numeric)],function(x){as.numeric(x)})
    } else {
      raw_data_vdjb <- NULL
    }

    # Antibody Capture Metrics
    if ("Antibody Capture" %in% modalities) {
      # Get ABC metrics
      cli_inform(message = "Reading {.field Antibody Capture} Metrics")
      raw_data <- read.csv(file = base_path, stringsAsFactors = FALSE)

      ABC_Metrics <- raw_data %>%
        filter(.data[["Grouped.By"]]== "Physical library ID" & .data[["Library.Type"]] == "Antibody Capture") %>%
        select(all_of(c("Metric.Name", "Metric.Value"))) %>%
        column_to_rownames("Metric.Name")

      current_metrics <- rownames(x = ABC_Metrics)

      ABC_Metrics <- ABC_Metrics %>%
        t() %>%
        data.frame()

      remaining_metrics <- setdiff(x = c("Antibody reads in cells",
                                          "Fraction antibody reads",
                                          "Fraction antibody reads in aggregate barcodes",
                                          "Mean antibody reads usable per cell",
                                          "Median UMI counts per cell",
                                          "Number of reads from cells associated with this sample"),
                                    y = current_metrics)

      ABC_Metrics2 <- raw_data %>%
        filter(.data[["Metric.Name"]] %in% remaining_metrics & .data[["Grouped.By"]]== "" & .data[["Library.Type"]] == "Antibody Capture"
        ) %>%
        select(all_of(c("Metric.Name", "Metric.Value"))) %>%
        column_to_rownames("Metric.Name") %>%
        t() %>%
        data.frame()

      raw_data_abc <- cbind(ABC_Metrics, ABC_Metrics2)

      column_numbers <- grep(pattern = ",", x = raw_data_abc[1, ])
      raw_data_abc[, c(column_numbers)] <- lapply(raw_data_abc[, c(column_numbers)], function(x) {
        as.numeric(gsub(",", "", x))})

      column_numbers_pct <- grep(pattern = "%", x = raw_data_abc[1, ])
      all_columns <- 1:ncol(x = raw_data_abc)

      column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

      raw_data_abc[,c(column_numbers_numeric)] <- lapply(raw_data_abc[,c(column_numbers_numeric)],function(x){as.numeric(x)})
    } else {
      raw_data_abc <- NULL
    }

    # combine outputs into a list
    data_list <- compact(
      list(
        multi_gex_metrics = raw_data_gex,
        multi_vdjt_metrics = raw_data_vdjt,
        multi_vdjb_metrics = raw_data_vdjb,
        multi_abc_metrics = raw_data_abc)
    )
    # return data list
    return(data_list)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GENE NAME/FILE CACHE HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' BiocFileCache Interface
#'
#' Internal function to manage/call BiocFileCache.
#'
#' @return cache
#'
#' @references \url{https://bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html#cache-to-manage-package-data}
#'
#' @noRd
#'

.get_bioc_cache <- function(
) {
  cache <- tools::R_user_dir(package = "scCustomize", which="cache")
  BiocFileCache::BiocFileCache(cache)
}


#' Download HGNC Dataset
#'
#' Internal function to download and cache the latest version of HGNC dataset for use with renaming genes.
#'
#' @param update logical, whether to manually override update parameters and download new data.
#'
#' @return path to data cache
#'
#' @references \url{https://bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html}
#'
#' @noRd
#'


download_hgnc_data <- function(
  update = NULL
) {
  # Get cache
  bfc <- .get_bioc_cache()

  # URL from https://www.genenames.org/download/statistics-and-files/
  hgnc_ftp_url <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

  # bfc <- BiocFileCache::BiocFileCache(hgnc_ftp_url)

  rid <- BiocFileCache::bfcquery(bfc, hgnc_ftp_url, "fpath")$rid
  if (!length(rid)) {               # not in cache, add but do not download
    rid <- names(BiocFileCache::bfcadd(bfc, hgnc_ftp_url, download = FALSE))
  }

  if (isTRUE(x = update)) {
    update <- update
  } else {
    update <- BiocFileCache::bfcneedsupdate(bfc, rid)  # TRUE if newly added or stale
  }

  # download & process
  if (!isFALSE(x = update)) {
    cli_inform(message = "Downloading HGNC data from: {.field {hgnc_ftp_url}}")
    BiocFileCache::bfcdownload(bfc, rid, ask = FALSE, FUN = process_hgnc_data)
  }

  rpath <- BiocFileCache::bfcrpath(bfc, rids=rid)    # path to processed result

  return(rpath)
}


#' Download MGI Dataset
#'
#' Internal function to download and cache the latest version of MGI dataset for use with renaming genes.
#'
#' @param update logical, whether to manually override update parameters and download new data.
#'
#'
#' @return path to data cache
#'
#' @references \url{https://bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html}
#'
#' @noRd
#'


download_mgi_data <- function(
  update = NULL
) {
  # Get cache
  bfc <- .get_bioc_cache()

  # URL from https://www.genenames.org/download/statistics-and-files/
  mgi_ftp_url <- "https://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt"

  # bfc <- BiocFileCache::BiocFileCache(mgi_ftp_url)

  rid <- BiocFileCache::bfcquery(bfc, mgi_ftp_url, "fpath")$rid
  if (!length(rid)) {               # not in cache, add but do not download
    rid <- names(BiocFileCache::bfcadd(bfc, mgi_ftp_url, download = FALSE))
  }

  if (isTRUE(x = update)) {
    update <- update
  } else {
    update <- BiocFileCache::bfcneedsupdate(bfc, rid)  # TRUE if newly added or stale
  }

  # download & process
  if (!isFALSE(x = update)) {
    cli_inform(message = "Downloading MGI data from: {.field {mgi_ftp_url}}")
    BiocFileCache::bfcdownload(bfc, rid, ask = FALSE, FUN = process_mgi_data)
  }

  rpath <- BiocFileCache::bfcrpath(bfc, rids=rid)    # path to processed result

  return(rpath)
}


#' Process HGNC Dataset
#'
#' Internal function process/filter and save HGNC dataset during cache process
#'
#' @param from input (cache location).
#' @param to output (cached data).
#'
#' @importFrom dplyr mutate select filter any_of contains
#' @importFrom magrittr "%>%"
#' @importFrom tidyr separate_wider_delim pivot_longer
#'
#' @return path to data cache
#'
#' @references \url{https://bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html}
#'
#' @noRd
#'

process_hgnc_data <- function(
  from,
  to
) {
  # read in data
  hgnc_full_data <- data.table::fread(file = from, data.table = FALSE, fill = TRUE)

  # filter data: Approved Genes > select relevant categories
  hgnc_filtered_data <- hgnc_full_data %>%
    filter(.data[["status"]] == "Approved") %>%
    select(any_of(c("hgnc_id", "symbol", "status", "alias_symbol", "prev_symbol", "date_symbol_changed", "entrez_id", "ensembl_gene_id")))

  # Select needed for renaming > split prev symbol column by number of additional columns needed > pivot wider without NAs > mutate
  hgnc_long_data <- hgnc_filtered_data %>%
    select(any_of(c("symbol", "prev_symbol"))) %>%
    separate_wider_delim(cols = "prev_symbol", delim = "|", names_sep = "_", names = NULL, too_few = "align_start") %>%
    pivot_longer(cols = contains("_symbol"),
                 names_to = "column",
                 values_to = "prev_symbol",
                 values_drop_na = TRUE) %>%
    mutate("prev_symbol" = ifelse(.data[["prev_symbol"]] %in% "", .data[["symbol"]], .data[["prev_symbol"]]))

  # save processed data
  saveRDS(hgnc_long_data, file = to)
  TRUE
}


#' Process MGI Dataset
#'
#' Internal function process/filter and save MGI dataset during cache process
#'
#' @param from input (cache location).
#' @param to output (cached data).
#'
#' @importFrom dplyr mutate select filter any_of contains
#' @importFrom magrittr "%>%"
#' @importFrom tidyr separate_wider_delim pivot_longer
#'
#' @return path to data cache
#'
#' @references \url{https://bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html}
#'
#' @noRd
#'

process_mgi_data <- function(
  from,
  to
) {
  # read in data
  mgi_full_data <- data.table::fread(file = from, data.table = FALSE)

  # Rename columns
  colnames(mgi_full_data) <- c("MGI Marker Accession ID", "Marker Symbol", "Status", "Marker Name", "cM Position", "Chromosome", "Type", "Secondary", "Entrez Gene ID", "Synonyms", "Feature Types", "Genome Coordinate Start", "Genome Coordinate End", "Strand", "BioTypes")

  # set accepted gene types
  accepted_biotypes <- c("protein coding gene", "lncRNA gene", "lincRNA gene", "antisense lncRNA gene")

  # filter data: Approved Genes > select relevant categories
  mgi_filtered_data <- mgi_full_data %>%
    filter(.data[["Status"]] == "O" & .data[["Type"]] == "Gene" & .data[["Chromosome"]] != "UN" & .data[["Feature Types"]] %in% accepted_biotypes) %>%
    select(any_of(c("MGI Marker Accession ID", "Marker Symbol", "Status", "Synonyms", "Entrez Gene ID", "Type", "Feature Types")))


  mgi_long_data <- mgi_filtered_data %>%
    select(any_of(c("Marker Symbol", "Synonyms"))) %>%
    separate_wider_delim(cols = "Synonyms", delim = "|", names_sep = "_", names = NULL, too_few = "align_start") %>%
    pivot_longer(cols = contains("Synonyms"),
                 names_to = "column",
                 values_to = "Synonyms",
                 values_drop_na = TRUE) %>%
    mutate("Synonyms" = ifelse(.data[["Synonyms"]] %in% "", .data[["Marker Symbol"]], .data[["Synonyms"]]))

  colnames(mgi_long_data) <- c("symbol", "column", "prev_symbol")
  # save processed data
  saveRDS(mgi_long_data, file = to)
  TRUE
}
