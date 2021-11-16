#' Set a default value if an object is NOT null
#'
#' @param lhs An object to set if it's NOT null
#' @param rhs The value to provide if x is NOT null
#'
#' @return lhs if lhs is null, else rhs
#'
#' @author Hadley Wickham
#' @references https://adv-r.hadley.nz/functions.html#missing-arguments
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
#' @references https://adv-r.hadley.nz/functions.html#missing-arguments
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

#' Stop function without error message
#'
#' Modifies R options within the function call only to hide the error message from `stop` while
#' keeping global options preserved outside of function.
#'
#' @return stops function without error message
#'
#' @author Stibu
#' @references https://stackoverflow.com/a/42945293/15568251
#' @details https://creativecommons.org/licenses/by-sa/3.0/
#'
#' @noRd
#'

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}


#' Check Seurat Object
#'
#' Checks if object is of class: Seurat and returns error message if not.
#'
#' @param seurat_object Seurat object name.
#'
#' @return stops function without error message
#'
#' @noRd
#'

Is_Seurat <- function(
  seurat_object
) {
  if (class(x = seurat_object)[[1]] != "Seurat") {
    stop("'seurat_object' provided is not an object of class: Seurat.")
  }
}


#' Check LIGER Object
#'
#' Checks if object is of class: liger and returns error message if not.
#'
#' @param liger_object liger object name.
#'
#' @return stops function without error message
#'
#' @noRd
#'

Is_LIGER <- function(
  liger_object
) {
  if (class(x = liger_object)[[1]] != "liger") {
    stop("'liger_object' provided is not an object of class: liger")
  }
}


#' Automatically calculate a point size for ggplot2-based scatter plots
#
#' It happens to look good
#'
#' @param data a single value length vector corresponding to the number of cells.
#' @param raster If TRUE, point size is set to 1
#'
#' @return The "optimal" point size for visualizing these data
#'
#' @noRd
#'
#' @references This function and documentation text are modified versions of the `AutoPointSize` function
#' and documentation from Seurat (https://github.com/satijalab/seurat/blob/master/R/visualization.R) (Licence: GPL-3).
#' This version has been modified to take single value length input instead of data.frame input.
#'

AutoPointSize_scCustom <- function(data, raster = NULL) {
  # for single value
  if (is.null(x = nrow(x = data)) && length(x = data) == 1 && is.numeric(x = data)) {
    return(ifelse(
      test = isTRUE(x = raster),
      yes = 1,
      no = min(1583 / data, 1)
    ))
  } else {
    # for data frame/object based values (from Seurat, see documentation)
    return(ifelse(
      test = isTRUE(x = raster),
      yes = 1,
      no = min(1583 / nrow(x = data), 1)
    ))
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
      stop("No requested assays found.")
    }

    # Return message of assays not found
    if (length(x = bad_assays) > 0 && omit_warn) {
      warning("The following assays were omitted as they were not found",
              ": ", glue_collapse_scCustom(input_string = bad_assays, and = TRUE))
    }

    # Combine into list and return
    assay_list <- list(
      found_assays = found_assays,
      bad_assays = bad_assays
    )
    return(assay_list)
  }

  # Print all found message if TRUE
  if (print_msg) {
    message("All assays present.")
  }

  # Return full input gene list.
  # Combine into list and return
  assay_list <- list(
    found_assays = assay_list,
    bad_assays = NULL
  )
  return(assay_list)
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
  if (and) {
    last_sep <- " and "
  } else {
    last_sep <- " or "
  }

  if (input_length < 3) {
    glue_collapse(x = input_string, sep = ", ", last = last_sep)
  } else {
    glue_collapse(x = input_string, sep = ", ", last = paste0(",", last_sep))
  }
}


#' Calculate the percentage of a vector above some threshold
#'
#' @param x Vector of values
#' @param threshold Threshold to use when calculating percentage
#'
#' @return Returns the percentage of `x` values above the given threshold
#'
#' @author Satija Lab & all co-Authors of Seurat Package
#' @references See Utilities.R in source code of Seurat https://github.com/satijalab/seurat/blob/master/R/utilities.R  (Licence: GPL-3).
#'
#' @note modified from Seurat version to return a percentage instead of proportion/decimal as part of `Percent_Expressing` function.  To be replaced following next Seurat version update.
#'
#' @keywords internal
#'
#' @noRd
#'

PercentAbove_Seurat <- function(x, threshold) {
  return((length(x = x[x > threshold]) / length(x = x))*100)
}
