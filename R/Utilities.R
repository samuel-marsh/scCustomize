#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Check if genes/features are present
#'
#' Check if genes are present in object and return vector of found genes.  Return warning messages for
#' genes not found.
#'
#' @param data Name of input data.  Currently only data of classes: Seurat, liger, data.frame,
#' dgCMatrix, dgTMatrix, tibble are accepted.  Gene_IDs must be present in rownames of the data.
#' @param gene_list vector of genes to check.
#' @param case_check logical. Whether or not to check if features are found if the case is changed from the
#' input list (Sentence case to Upper and vice versa).  Default is TRUE.
#' @param case_check_msg logical. Whether to print message to console if alternate case features are found
#' in addition to inclusion in returned list.  Default is TRUE.
#' @param print_msg logical. Whether message should be printed if all features are found.  Default is TRUE.
#' @param omit_warn logical. Whether to print message about features that are not found in current object.
#'  Default is TRUE.
#' @param return_none logical. Whether list of found vs. bad features should still be returned if no
#' features are found.  Default is FALSE.
#' @param seurat_assay Name of assay to pull feature names from if `data` is Seurat Object.
#' Defaults to  `DefaultAssay(OBJ)` if NULL.
#'
#' @importFrom purrr reduce
#' @importFrom stringr str_to_upper str_to_sentence
#'
#' @return A list of length 3 containing 1) found features, 2) not found features, 3) features found if
#' case was modified.
#'
#' @export
#'
#' @concept helper_util
#'
#' @examples
#' \dontrun{
#' features <- Gene_Present(data = obj_name, gene_list = DEG_list, print_msg = TRUE, case_check = TRUE)
#' found_features <- features[[1]]
#' }
#'

Gene_Present <- function(
  data,
  gene_list,
  case_check = TRUE,
  case_check_msg = TRUE,
  print_msg = TRUE,
  omit_warn = TRUE,
  return_none = FALSE,
  seurat_assay = NULL
) {
  # Check object type
  # Seurat
  accepted_types <- c("data.frame", "dgCMatrix", "dgTMatrix", "tibble")
  if (inherits(x = data, what = "Seurat")) {
    # set assay (if null set to active assay)
    assay <- seurat_assay %||% DefaultAssay(object = data)

    possible_features <- rownames(x = GetAssayData(object = data, assay = assay))
  } else if ((class(x = data)[[1]] == "liger")) {
    # get complete gene list
    length_liger <- length(data@raw.data)

    list_genes <- lapply(1:length_liger, function(x){
      rownames(x = data@raw.data[[x]])
    })

    possible_features <- reduce(list_genes, function(x, y) {
      union(x = x, y = y)})
  } else if ((class(x = data) %in% accepted_types)) {
    possible_features <- rownames(x = data)
  } else {
    all_accepted <- c(accepted_types, "Seurat", "liger")
    cli_abort(message = c("Input data is currently accepted only in the following formats:",
                          "i" = "{.field {glue_collapse_scCustom(input_string = all_accepted, and = FALSE)}}.")
    )
  }

  # If any features not found
  if (any(!gene_list %in% possible_features)) {
    bad_features <- gene_list[!gene_list %in% possible_features]
    found_features <- gene_list[gene_list %in% possible_features]
    if (length(x = found_features) == 0) {
      if (return_none) {
        # Combine into list and return
        feature_list <- list(
          found_features = NULL,
          bad_features = bad_features,
          wrong_case_found_features = NULL
        )
        return(feature_list)
      } else {
        cli_abort(message ="No requested features found.")
      }
    }

    # Return message of features not found
    if (length(x = bad_features) > 0 && omit_warn) {
      cli_warn(message = c("The following features were omitted as they were not found:",
                            "i" = "{.field {glue_collapse_scCustom(input_string = bad_features, and = TRUE)}}")
      )
    }

    # Check if features found if case is changed.
    if (case_check) {
      upper_bad_features <- str_to_upper(string = bad_features)
      upper_found_features <- upper_bad_features[upper_bad_features %in% possible_features]

      sentence_bad_features <- str_to_sentence(string = bad_features)
      sentence_found_features <- sentence_bad_features[sentence_bad_features %in% possible_features]

      # Combine case check
      wrong_case_found_features <- c(upper_found_features, sentence_found_features)

      # Additional messages if found.
      if (length(x = wrong_case_found_features) > 0) {
        if (case_check_msg) {
          cli_abort(message = c("NOTE: However, the following features were found: {.field {glue_collapse_scCustom(input_string = wrong_case_found_features, and = TRUE)}}",
                                "i" = "Please check intended case of features provided.")
          )
        }
        # Combine into list and return
        feature_list <- list(
          found_features = found_features,
          bad_features = bad_features,
          wrong_case_found_features = wrong_case_found_features
        )
        return(feature_list)
      }
    }
    # Combine into list and return
    feature_list <- list(
      found_features = found_features,
      bad_features = bad_features,
      wrong_case_found_features = "NA (check not performed.  Set 'case_check = TRUE' to perform check."
    )
    return(feature_list)
  }

  # Print all found message if TRUE
  if (print_msg) {
    cli_inform(message = "All features present.")
  }

  # Return full input gene list.
  # Combine into list and return
  feature_list <- list(
    found_features = gene_list,
    bad_features = NULL,
    wrong_case_found_features = NULL
  )
  return(feature_list)
}


#' Check for alternate case features
#
#' Checks Seurat object for the presence of features with the same spelling but alternate case.
#'
#' @param seurat_object Seurat object name.
#' @param gene_list vector of genes to check.
#' @param case_check_msg logical. Whether to print message to console if alternate case features are
#' found in addition to inclusion in returned list.  Default is TRUE.
#' @param return_features logical. Whether to return vector of alternate case features.  Default is TRUE.
#'
#' @return If features found returns vector of found alternate case features and prints message depending on
#' parameters specified.
#' @export
#'
#' @concept helper_util
#'
#' @examples
#' \dontrun{
#' alt_features <- Case_Check(seurat_object = obj_name, gene_list = DEG_list)
#' }
#'

Case_Check <- function(
  seurat_object,
  gene_list,
  case_check_msg = TRUE,
  return_features = TRUE
) {
  # get all features
  possible_features <- rownames(x = GetAssayData(object = seurat_object))

  upper_bad_features <- str_to_upper(string = gene_list)
  upper_found_features <- upper_bad_features[upper_bad_features %in% possible_features]

  sentence_bad_features <- str_to_sentence(string = gene_list)
  sentence_found_features <- sentence_bad_features[sentence_bad_features %in% possible_features]

  # Combine case check
  wrong_case_found_features <- c(upper_found_features, sentence_found_features)

  # Additional messages if found.
  if (length(x = wrong_case_found_features) > 0) {
    if (case_check_msg) {
      cli_inform(message = c("NOTE: However, the following features were found: {.field {glue_collapse_scCustom(input_string = wrong_case_found_features, and = TRUE)}}",
                            "i" = "Please check intended case of features provided.")
      )
    }
    if (return_features) {
      return(wrong_case_found_features)
    }
  }
}


#' Check if meta data are present
#'
#' Check if meta data columns are present in object and return vector of found columns
#' Return warning messages for meta data columns not found.
#'
#' @param seurat_object object name.
#' @param meta_col_names vector of column names to check.
#' @param print_msg logical. Whether message should be printed if all features are found.  Default is TRUE.
#' @param omit_warn logical. Whether to print message about features that are not found in current object. Default is TRUE.
#' @param abort logical. Whether or not to stop function and print stop message if no input `meta_col_names` are found.  Default is TRUE.
#'
#' @return vector of meta data columns that are present
#'
#' @import cli
#'
#' @export
#'
#' @concept helper_util
#'
#' @examples
#' \dontrun{
#' meta_variables <- Meta_Present(seurat_object = obj_name, gene_list = DEG_list, print_msg = TRUE)
#' }
#'

Meta_Present <- function(
  seurat_object,
  meta_col_names,
  print_msg = TRUE,
  omit_warn = TRUE,
  abort = TRUE
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # get all features
  possible_features <- colnames(x = seurat_object@meta.data)

  # If any features not found
  if (any(!meta_col_names %in% possible_features)) {
    bad_meta <- meta_col_names[!meta_col_names %in% possible_features]
    found_meta <- meta_col_names[meta_col_names %in% possible_features]

    if (abort) {
      if (length(x = found_meta) < 1) {
        cli_abort(message = c("No meta data columns found.",
                              "i" = "The following @meta.data columns were not found: {glue_collapse_scCustom(input_string = bad_meta, and = TRUE)}")
        )
      }
    }

    # Return message of features not found
    if (length(x = bad_meta) > 0 && omit_warn) {
      cli_warn(message = c("The following @meta.data columns were omitted as they were not found:",
                            "i" = "{.field {glue_collapse_scCustom(input_string = bad_meta, and = TRUE)}}")
      )
    }

    # Return the found features omitting the not found ones.
    meta_list <- list(
      found_meta = found_meta,
      bad_meta = bad_meta
    )

    return(meta_list)
  }

  # Print all found message if TRUE
  if (print_msg) {
    cli_inform(message = "All @meta.data columns present.")
  }

  # Return full input gene list.
  meta_list <- list(
    found_meta = meta_col_names,
    bad_meta = NULL
  )
  return(meta_list)
}

#' Check if meta data columns are numeric
#'
#' Check if any present meta data columns are numeric and returns vector of valid numeric columns.
#' Issues warning message if any columns not in numeric form.
#'
#' @param data a data.frame contain meta.data.
#'
#' @return vector of meta data columns that are numeric.
#'
#' @importFrom dplyr filter pull
#' @importFrom magrittr "%>%"
#' @importFrom tibble rownames_to_column
#'
#' @export
#'
#' @concept helper_util
#'
#' @examples
#' \dontrun{
#' numeric_meta_columns <- Meta_Numeric(data = meta_data)
#' }
#'

Meta_Numeric <- function(
  data
) {
  # Logical data.frame of is.numeric results
  all_numeric <- data.frame(sapply(colnames(x = data), function(x) {
    is.numeric(x = data[[x]])
  }))

  colnames(all_numeric) <- "Is_Numeric"

  # Pull results into vectors
  invalid_variables <- all_numeric %>%
    rownames_to_column("variables") %>%
    filter(.data[["Is_Numeric"]] == FALSE) %>%
    pull(.data[["variables"]])

  valid_variables <- all_numeric %>%
    rownames_to_column("variables") %>%
    filter(.data[["Is_Numeric"]] == TRUE) %>%
    pull(.data[["variables"]])

  # Warn if columns are not numeric
  if (length(x = invalid_variables) > 0) {
    cli_warn(message = c("Some of the meta.data columns provided are not in numeric form and will be excluded from results:",
                         "i" = "{.field {paste(shQuote(invalid_variables, type = 'cmd'), collapse=', ')}}"))
  }

  # Return valid column names
  return(valid_variables)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### DATA ACCESS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Get meta data from object
#'
#' Quick function to properly pull meta.data from objects.
#'
#' @param object list of matrices to merge.
#'
#' @importFrom methods slot
#'
#' @return A data.frame
#'
#' @export
#'
#' @concept helper_util
#'
#' @examples
#' \dontrun{
#' meta_data <- Fetch_Meta(object = pbmc)
#' }
#'

Fetch_Meta <- function(
  object
) {
  # Check Seurat
  Is_Seurat(seurat_object = object)

  # Pull meta data
  object_meta <- slot(object = object, name = "meta.data")

  return(object_meta)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### MATRIX HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Merge a list of Sparse Matrices
#'
#' Enables easy merge of a list of sparse matrices
#'
#' @param matrix_list list of matrices to merge.
#' @param add_cell_ids a vector of sample ids to add as prefix to cell barcode during merge.
#' @param prefix logical.  Whether `add_cell_ids` should be added as prefix to current cell barcodes/names
#' or as suffix to current cell barcodes/names.  Default is TRUE, add as prefix.
#' @param cell_id_delimiter The delimiter to use when adding cell id prefix/suffix.  Default is "_".
#'
#' @references Original function is part of LIGER package as non-exported function
#' \url{https://github.com/welch-lab/liger/blob/master/R/utilities.R} (License: GPL-3).
#' Function was modified for use in scCustomize (add progress bar, prefix vs. suffix, and delimiter options).
#'
#' @import Matrix
#' @importFrom dplyr intersect
#' @importFrom magrittr "%>%"
#'
#' @return A sparse Matrix
#'
#' @export
#'
#' @concept helper_util
#'
#' @examples
#' \dontrun{
#' data_list <- Read10X_GEO(...)
#' merged <- Merge_Sparse_Data_All(matrix_list = data_list, add_cell_ids = names(data_list),
#' prefix = TRUE, cell_id_delimiter = "_")
#' }
#'

Merge_Sparse_Data_All <- function(
  matrix_list,
  add_cell_ids = NULL,
  prefix = TRUE,
  cell_id_delimiter = "_"
) {
  # Check all barcodes are unique to begin with
  duplicated_barcodes <- matrix_list %>%
    lapply(colnames) %>%
    unlist() %>%
    duplicated() %>%
    any()

  if (duplicated_barcodes && is.null(x = add_cell_ids)) {
    stop("There are overlapping cell barcodes present in the input matrices.  Please provide prefixes/suffixes to 'add_cell_ids' parameter to make unique.")
  }

  # Check right number of suffix/prefix ids are provided
  if (!is.null(x = add_cell_ids) && length(x = add_cell_ids) != length(x = matrix_list)) {
    stop("The number of prefixes in `add_cell_ids` must be equal to the number of matrices supplied to `matrix_list`.")
  }

  if (!is.null(x = add_cell_ids)) {
    # check barcodes will be unique after adding prefixes/suffixes
    all_names <- lapply(1:length(x = matrix_list), function(i){
      cell_names <- colnames(matrix_list[[i]])
    })

    new_names <- lapply(X = 1:length(x = matrix_list), function(x){
      colnames(matrix_list[[x]]) <- paste0(add_cell_ids[x], cell_id_delimiter, colnames(matrix_list[[x]]))
    })

    are_duplicates <- unlist(new_names) %>%
      duplicated() %>%
      any()

    if (are_duplicates) {
      stop("Supplied 'add_cell_ids' will result in overlapping barcodes names.  If provided cell prefixes/suffixes are not unique please change and re-run.")
    }
  }

  # Use summary to convert the sparse matrices into three-column indexes where i are the
  # row numbers, j are the column numbers, and x are the nonzero entries
  col_offset <- 0
  allGenes <- unique(unlist(lapply(matrix_list, rownames)))
  allCells <- c()
  message("Preparing & merging matrices.")
  pb <- txtProgressBar(min = 0, max = length(x = matrix_list), style = 3, file = stderr())
  for (i in 1:length(matrix_list)) {
    curr <- matrix_list[[i]]
    curr_s <- summary(curr)

    # Now, alter the indexes so that the two 3-column matrices can be properly merged.
    # First, make the current and full column numbers non-overlapping.
    curr_s[, 2] <- curr_s[, 2] + col_offset

    # Update full cell names
    if (!is.null(x = add_cell_ids)) {
      if (prefix) {
        cellnames <- paste0(add_cell_ids [i], cell_id_delimiter, colnames(curr))
      } else {
        cellnames <- paste0(colnames(curr), cell_id_delimiter, add_cell_ids [i])
      }
    } else {
      cellnames <- colnames(curr)
    }
    allCells <- c(allCells, cellnames)

    # Next, change the row (gene) indexes so that they index on the union of the gene sets,
    # so that proper merging can occur.
    idx <- match(rownames(curr), allGenes)
    newgenescurr <- idx[curr_s[, 1]]
    curr_s[, 1] <- newgenescurr

    # Now bind the altered 3-column matrices together, and convert into a single sparse matrix.
    if (!exists("full_mat")) {
      full_mat <- curr_s
    } else {
      full_mat <- rbind(full_mat, curr_s)
    }
    col_offset <- length(x = allCells)
    setTxtProgressBar(pb = pb, value = i)
  }
  close(con = pb)
  message("Creating final sparse matrix.")
  M <- sparseMatrix(
    i = full_mat[, 1],
    j = full_mat[, 2],
    x = full_mat[, 3],
    dims = c(
      length(x = allGenes),
      length(x = allCells)
    ),
    dimnames = list(
      allGenes,
      allCells
    )
  )
  return(M)
}


#' Check Matrix Validity
#'
#' Native implementation of SeuratObjects CheckMatrix but with modified warning messages.
#'
#' @param object A matrix
#' @param checks Type of checks to perform, choose one or more from:
#' \itemize{
#'  \item \dQuote{\code{infinite}}: Emit a warning if any value is infinite
#'  \item \dQuote{\code{logical}}: Emit a warning if any value is a logical
#'  \item \dQuote{\code{integer}}: Emit a warning if any value is \emph{not}
#'   an integer
#'  \item \dQuote{\code{na}}: Emit a warning if any value is an \code{NA}
#'   or \code{NaN}
#' }
#'
#' @return Emits warnings for each test and invisibly returns \code{NULL}
#'
#' @import cli
#' @importFrom methods slot
#'
#' @references Re-implementing `CheckMatrix` only for sparse matrices with modified warning messages.  Original function from SeuratObject \url{https://github.com/mojaveazure/seurat-object/blob/9c0eda946e162d8595696e5280a6ecda6284db39/R/utils.R#L625-L650} (License: MIT).
#'
#' @export
#'
#' @concept helper_util
#'
#' @examples
#' \dontrun{
#' mat <- Read10X(...)
#' CheckMatrix_scCustom(object = mat)
#' }
#'

CheckMatrix_scCustom <- function(
  object,
  checks = c('infinite', 'logical', 'integer', 'na')
) {
  checks <- match.arg(arg = checks, several.ok = TRUE)
  x <- slot(object = object, name = 'x')
  for (i in checks) {
    switch(
      EXPR = i,
      'infinite' = if (any(is.infinite(x = x))) {
        cli_warn(message = "Input matrix contains infinite values")
      },
      'logical' = if (any(is.logical(x = x))) {
        cli_warn(message = "Input matrix contains logical values")
      },
      'integer' = if (!all(round(x = x) == x, na.rm = TRUE)) {
        cli_warn(message = c("Input matrix contains non-integer values.",
                             "*" = "Data may represent normalized or scaled values.",
                             "i" = "Take into account when performing analysis.")
        )
      },
      'na' = if (anyNA(x = x)) {
        cli_warn(message = "Input matrix contains NA/NaN values")
      },
    )
  }
  return(invisible(x = NULL))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### BARCODE UTILS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Replace barcode suffixes
#'
#' Replace barcode suffixes in matrix, data.frame, or list of matrices/data.frames
#'
#' @param data Either matrix/data.frame or list of matrices/data.frames with the cell barcodes in the column names.
#' @param current_suffix a single value or vector of values representing current barcode suffix.
#' If suffix is the same for all matrices/data.frames in list only single value is required.
#' @param new_suffix a single value or vector of values representing new barcode suffix to be added.
#' If desired suffix is the same for all matrices/data.frames in list only single value is required.
#' If no suffix is desired set `new_suffix = ""`.`
#'
#' @export
#'
#' @return matrix or data.frame with new column names.
#'
#' @concept helper_util
#'
#' @examples
#' \dontrun{
#' dge_matrix <- Replace_Suffix(data = dge_matrix, current_suffix = "-1", new_suffix = "-2")
#' }
#'

Replace_Suffix <- function(
  data,
  current_suffix,
  new_suffix
) {
  # is data a list
  if (inherits(x = data, what = "list")) {
    # Make list of current names
    current_cell_names <- lapply(X = 1:length(x = data), function(x) {
      cell_names <- colnames(x = data[[x]])
    })
    # Get current suffix regexps
    if (length(x = current_suffix) == 1) {
      current_suffix_regexp <- rep(paste0("\\", current_suffix, "$"), length(x = data))
    }
    if (length(x = current_suffix) == length(x = data)) {
      current_suffix_regexp <- sapply(X = 1:length(x = data), FUN = function(i){
        current_suffix_regexp <- paste0("\\", current_suffix[[i]], "$")
      })
    }
    if (length(x = current_suffix) != length(x = data) && length(x = current_suffix) != 1) {
      stop("'current_suffix' must either be single value or a vector of values equal in length to number of datasets in 'data'.")
    }

    # Check new suffix
    if (length(x = new_suffix) == 1) {
      new_suffix <- rep(new_suffix, length(x = data))
    }
    if (length(x = new_suffix) == length(x = data)) {
      new_suffix <- new_suffix
    }
    if (length(x = new_suffix) != length(x = data) && length(x = new_suffix) != 1) {
      stop("'new_suffix' must either be single value or a vector of values equal in length to number of datasets in 'data'.")
    }

    # Is current suffix found in all cell names
    check_suffixes <- sapply(1:length(data), FUN = function(j){
      all(grepl(pattern = current_suffix_regexp[[j]], x = current_cell_names[[j]]))
    })

    if (all(check_suffixes) != TRUE) {
      stop("One or more 'current_suffixes' do not match cell names in data.  Check inputs.")
    }

    # Create list of string of new names
    new_cell_names_list <- lapply(1:length(x = data), function(k){
      new_cell_names <- gsub(pattern = current_suffix_regexp[[k]], replacement = new_suffix[[k]], x = current_cell_names[[k]])
    })


    # replace names and return data
    data_mod <- lapply(1:length(x = data), function(m){
      data_single <- data[[m]]
      colnames(x = data_single) <- new_cell_names_list[[m]]
      return(data_single)
    })
    # Add names back to output
    names(data_mod) <- names(data)
    return(data_mod)

  } else {
    # for data.frames and individual matrices
    current_cell_names <- colnames(x = data)
    current_suffix_regexp <- paste0("\\", current_suffix, "$")
    # Is current suffix found in cell names
    if (all(grepl(pattern = current_suffix_regexp, x = current_cell_names)) == FALSE) {
      stop("Supplied 'current_suffix': ", current_suffix, " was not found in the cell names of data provided.")
    }
    # Create string of new names
    new_cell_names <- gsub(pattern = current_suffix_regexp, replacement = new_suffix, x = current_cell_names)

    # replace names and return data
    colnames(x = data) <- new_cell_names
    return(data)
  }
}


#' Change barcode suffix delimiter
#'
#' Change barcode suffix delimiter from list of data.frames/matrices or single data.frame/matrix
#'
#' @param data Either matrix/data.frame or list of matrices/data.frames with the cell barcodes in the column names.
#' @param current_delim a single value of current delimiter.
#' @param new_delim a single value of new delimiter desired.
#'
#' @importFrom stringi stri_replace_last_fixed
#'
#' @export
#'
#' @return matrix or data.frame with new column names.
#'
#' @concept helper_util
#'
#' @examples
#' \dontrun{
#' dge_matrix <- Change_Delim_Suffix(data = dge_matrix, current_delim = ".", new_delim = "-")
#' }
#'

Change_Delim_Suffix <- function(
  data,
  current_delim,
  new_delim
) {
  # is data a list
  if (inherits(x = data, what = "list")) {
    # Make list of current names
    current_cell_names <- lapply(X = 1:length(x = data), function(x) {
      cell_names <- colnames(x = data[[x]])
    })

    # Is current suffix delim found in all cell names
    check_suffix_delim <- sapply(1:length(data), FUN = function(j){
      all(grepl(pattern = current_delim, x = current_cell_names[[j]], fixed = TRUE))
    })

    if (all(check_suffix_delim) != TRUE) {
      stop("One or more 'current_delim' does not match cell names in data.  Check inputs.")
    }

    # Create list of string of new names
    new_cell_names_list <- lapply(1:length(x = data), function(k){
      new_cell_names <- stri_replace_last_fixed(pattern = current_delim, str = current_cell_names[[k]], replacement = new_delim)
    })

    # replace names and return data
    data_mod <- lapply(1:length(x = data), function(m){
      data_single <- data[[m]]
      colnames(x = data_single) <- new_cell_names_list[[m]]
      return(data_single)
    })
    # Add names back to output
    names(data_mod) <- names(data)
    return(data_mod)
  } else {
    # for data.frames and individual matrices
    current_cell_names <- colnames(x = data)

    # Is current suffix found in cell names
    if (all(grepl(pattern = current_delim, x = current_cell_names, fixed = TRUE)) != TRUE) {
      stop("Supplied 'current_delim': ", current_delim, " was not found in the cell names of data provided.")
    }

    # Create string of new names
    new_cell_names <-  stri_replace_last_fixed(pattern = current_delim, str = current_cell_names, replacement = new_delim)

    # replace names and return data
    colnames(x = data) <- new_cell_names
    return(data)
  }
}


#' Change barcode prefix delimiter
#'
#' Change barcode prefix delimiter from list of data.frames/matrices or single data.frame/matrix
#'
#' @param data Either matrix/data.frame or list of matrices/data.frames with the cell barcodes in the column names.
#' @param current_delim a single value of current delimiter.
#' @param new_delim a single value of new delimiter desired.
#'
#' @importFrom stringi stri_replace_first_fixed
#'
#' @export
#'
#' @return matrix or data.frame with new column names.
#'
#' @concept helper_util
#'
#' @examples
#' \dontrun{
#' dge_matrix <- Change_Delim_Prefix(data = dge_matrix, current_delim = ".", new_delim = "-")
#' }
#'

Change_Delim_Prefix <- function(
  data,
  current_delim,
  new_delim
) {
  # is data a list
  if (inherits(x = data, what = "list")) {
    # Make list of current names
    current_cell_names <- lapply(X = 1:length(x = data), function(x) {
      cell_names <- colnames(x = data[[x]])
    })

    # Is current prefix delim found in all cell names
    check_prefix_delim <- sapply(1:length(data), FUN = function(j){
      all(grepl(pattern = current_delim, x = current_cell_names[[j]], fixed = TRUE))
    })

    if (all(check_prefix_delim) != TRUE) {
      stop("One or more 'current_delim' does not match cell names in data.  Check inputs.")
    }

    # Create list of string of new names
    new_cell_names_list <- lapply(1:length(x = data), function(k){
      new_cell_names <- stri_replace_first_fixed(pattern = current_delim, str = current_cell_names[[k]], replacement = new_delim)
    })

    # replace names and return data
    data_mod <- lapply(1:length(x = data), function(m){
      data_single <- data[[m]]
      colnames(x = data_single) <- new_cell_names_list[[m]]
      return(data_single)
    })
    # Add names back to output
    names(data_mod) <- names(data)
    return(data_mod)
  } else {
    # for data.frames and individual matrices
    current_cell_names <- colnames(x = data)

    # Is current prefix found in cell names
    if (all(grepl(pattern = current_delim, x = current_cell_names, fixed = TRUE)) != TRUE) {
      stop("Supplied 'current_delim': ", current_delim, " was not found in the cell names of data provided.")
    }

    # Create string of new names
    new_cell_names <-  stri_replace_first_fixed(pattern = current_delim, str = current_cell_names, replacement = new_delim)

    # replace names and return data
    colnames(x = data) <- new_cell_names
    return(data)
  }
}


#' Change all delimiters in cell name
#'
#' Change all instances of delimiter in cell names from list of data.frames/matrices or single data.frame/matrix
#'
#' @param data Either matrix/data.frame or list of matrices/data.frames with the cell barcodes in the column names.
#' @param current_delim a single value of current delimiter.
#' @param new_delim a single value of new delimiter desired.
#'
#' @export
#'
#' @return matrix or data.frame with new column names.
#'
#' @concept helper_util
#'
#' @examples
#' \dontrun{
#' dge_matrix <- Change_Delim_All(data = dge_matrix, current_delim = ".", new_delim = "-")
#' }
#'

Change_Delim_All <- function(
  data,
  current_delim,
  new_delim
) {
  # is data a list
  if (inherits(x = data, what = "list")) {
    # Make list of current names
    current_cell_names <- lapply(X = 1:length(x = data), function(x) {
      cell_names <- colnames(x = data[[x]])
    })

    # Is current prefix delim found in all cell names
    check_prefix_delim <- sapply(1:length(data), FUN = function(j){
      all(grepl(pattern = current_delim, x = current_cell_names[[j]], fixed = TRUE))
    })

    if (all(check_prefix_delim) != TRUE) {
      stop("One or more 'current_delim' does not match cell names in data.  Check inputs.")
    }

    # Create list of string of new names
    new_cell_names_list <- lapply(1:length(x = data), function(k){
      new_cell_names <- gsub(pattern = current_delim, x = current_cell_names[[k]], replacement = new_delim, fixed = TRUE)
    })

    # replace names and return data
    data_mod <- lapply(1:length(x = data), function(m){
      data_single <- data[[m]]
      colnames(x = data_single) <- new_cell_names_list[[m]]
      return(data_single)
    })
    # Add names back to output
    names(data_mod) <- names(data)
    return(data_mod)
  } else {
    # for data.frames and individual matrices
    current_cell_names <- colnames(x = data)

    # Is current prefix found in cell names
    if (all(grepl(pattern = current_delim, x = current_cell_names, fixed = TRUE)) != TRUE) {
      stop("Supplied 'current_delim': ", current_delim, " was not found in the cell names of data provided.")
    }

    # Create string of new names
    new_cell_names <-  gsub(pattern = current_delim, x = current_cell_names, replacement = new_delim, fixed = TRUE)

    # replace names and return data
    colnames(x = data) <- new_cell_names
    return(data)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### CLUSTER MARKERS/ANNOTATION ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Add percentage difference to DE results
#'
#' Adds new column labeled "pct_diff" to the data.frame output of \code{\link[Seurat]{FindMarkers}},  \code{\link[Seurat]{FindAllMarkers}}, or other DE test data.frames.
#'
#' @param marker_dataframe data.frame containing the results of \code{\link[Seurat]{FindMarkers}},  \code{\link[Seurat]{FindAllMarkers}}, or other DE test data.frame.
#' @param pct.1_name the name of data.frame column corresponding to percent expressed in group 1.
#' Default is Seurat default "pct.1".
#' @param pct.2_name the name of data.frame column corresponding to percent expressed in group 2.
#' Default is Seurat default "pct.2".
#' @param overwrite logical.  If the `marker_dataframe` already contains column named "pct_diff" whether to
#'  overwrite or return error message.  Default is FALSE.
#'
#' @import cli
#' @importFrom dplyr mutate
#' @importFrom magrittr "%>%"
#'
#' @return Returns input `marker_dataframe` with additional "pct_diff" column.
#'
#' @export
#'
#' @concept marker_annotation_util
#'
#' @examples
#' \dontrun{
#' marker_df <- FindAllMarkers(object = obj_name)
#' marker_df <- Add_Pct_Diff(marker_dataframe = marker_df)
#' # or piped with function
#' marker_df <- FindAllMarkers(object = obj_name) %>%
#'   Add_Pct_Diff()
#' }
#'

Add_Pct_Diff <- function(
  marker_dataframe,
  pct.1_name = "pct.1",
  pct.2_name = "pct.2",
  overwrite = FALSE
) {
  # Check if percent difference exists already
  if ("pct_diff" %in% colnames(marker_dataframe)) {
    df_name <- deparse(expr = substitute(expr = marker_dataframe))
    if (!overwrite) {
      cli_abort(message = c("'pct_diff' column already present in `marker_dataframe: '{name}'.",
                            "i" = "To overwrite previous results set `overwrite = TRUE`.")
      )
    } else {
      cli_inform(message = c("'pct_diff' column already present in `marker_dataframe: '{name}'.",
                            "i" = "Overwriting column as overwrite = TRUE.")
      )
    }
  }
  # Add percentage difference
  pct_diff_df <- marker_dataframe %>%
    mutate(pct_diff = .data[[pct.1_name]] - .data[[pct.2_name]])
  return(pct_diff_df)
}


#' Extract Top N Marker Genes
#'
#' Extract vector gene list (or named gene vector) from data.frame results of \code{\link[Seurat]{FindAllMarkers}}
#'  or similar analysis.
#'
#' @param marker_dataframe data.frame output from \code{\link[Seurat]{FindAllMarkers}} or similar analysis.
#' @param num_genes number of genes per group (e.g., cluster) to include in output list.
#' @param group_by column name of `marker_dataframe` to group data by.  Default is "cluster" based on
#'  \code{\link[Seurat]{FindAllMarkers}}.
#' @param rank_by column name of `marker_dataframe` to rank data by when selecting `num_genes` per `group_by`.
#' Default is "avg_log2FC" based on \code{\link[Seurat]{FindAllMarkers}}.
#' @param gene_column column name of `marker_dataframe` that contains the gene IDs.  Default is "gene"
#' based on \code{\link[Seurat]{FindAllMarkers}}.
#' @param gene_rownames_to_column logical. Whether gene IDs are stored in rownames and should be moved to
#' column.  Default is FALSE.
#' @param data_frame Logical, whether or not to return filtered data.frame of the original `markers_dataframe` or
#' to return a vector of gene IDs.  Default is FALSE.
#' @param named_vector Logical, whether or not to name the vector of gene names that is returned by the function.
#' If `TRUE` will name the vector using the column provided to `group_by`.  Default is TRUE.
#' @param make_unique Logical, whether an unnamed vector should return only unique values.  Default is FALSE.
#' Not applicable when `data_frame = TRUE` or `named_vector = TRUE`.
#'
#' @import cli
#' @importFrom dplyr group_by slice_max
#' @importFrom magrittr "%>%"
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @return filtered data.frame, vector, or named vector containing gene IDs.
#'
#' @export
#'
#' @concept marker_annotation_util
#'
#' @examples
#' \dontrun{
#' top10_genes <- Extract_Top_Markers(marker_dataframe = markers_results, num_genes = 10,
#' group_by = "cluster", rank_by = "avg_log2FC")
#' }
#'

Extract_Top_Markers <- function(
  marker_dataframe,
  num_genes = 10,
  group_by = "cluster",
  rank_by = "avg_log2FC",
  gene_column = "gene",
  gene_rownames_to_column = FALSE,
  data_frame = FALSE,
  named_vector = TRUE,
  make_unique = FALSE
) {
  # Check ranking factor in marker data.frame
  if (!rank_by %in% colnames(x = marker_dataframe)) {
    cli_abort(message = "`rank_by`: '{rank_by}' not found in column names of `marker_dataframe`.")
  }

  # Check grouping factor in marker data.frame
  if (!is.null(x = group_by)) {
    if (!group_by %in% colnames(x = marker_dataframe)) {
      cli_abort(message = "`group_by`: '{group_by}' not found in column names of `marker_dataframe`.")
    }
  }

  # Check gene column is present
  if (!gene_column %in% colnames(x = marker_dataframe) && !gene_rownames_to_column) {
    cli_abort(message = c("`gene_column`: '{gene_column}' not found in column names of `marker_dataframe.",
                          "i" = "Set `gene_rownames_to_column` to move genes from rownames to column.")
    )
  }


  # create filtered data.frame
  if (is.null(x = group_by)) {
    filtered_markers <- marker_dataframe %>%
      rownames_to_column("rownames") %>%
      slice_max(n = num_genes, order_by = .data[[rank_by]]) %>%
      column_to_rownames("rownames")
  } else {
    filtered_markers <- marker_dataframe %>%
      rownames_to_column("rownames") %>%
      group_by(.data[[group_by]]) %>%
      slice_max(n = num_genes, order_by = .data[[rank_by]]) %>%
      column_to_rownames("rownames")
  }

  if (gene_rownames_to_column) {
    filtered_markers <- filtered_markers %>%
      rownames_to_column(gene_column)
  }

  # return data.frame
  if (data_frame) {
    return(filtered_markers)
  }

  # pull gene list
  gene_list <- filtered_markers[[gene_column]]

  # should gene list be named
  # check naming
  if (named_vector && is.null(x = group_by)) {
    cli_warn(message = c("Cannot return named vector if `group_by` is NULL.",
                         "i" = "Returning unnamed vector.")
    )
  }

  if (named_vector && !is.null(x = group_by)) {
    if (make_unique) {
      cli_abort(message = "Cannot return unique list if 'named_vector = TRUE'.")
    }
    names(gene_list) <- filtered_markers[[group_by]]
    return(gene_list)
  }

  # make unique
  if (make_unique) {
    gene_list <- unique(x = gene_list)
  }

  return(gene_list)
}


#' Create cluster annotation csv file
#'
#' create annotation file
#'
#' @param file_path path to directory to save file.  Default is current working directory.
#' @param file_name name to use for annotation file.  Function automatically adds file type ".csv" suffix.
#' Default is "cluster_annotation".
#'
#' @import cli
#' @importFrom utils write.csv
#'
#' @export
#'
#' @return No value returned.  Creates .csv file.
#'
#' @concept marker_annotation_util
#'
#' @examples
#' \dontrun{
#' Create_Cluster_Annotation_File(file_path = "cluster_annotation_folder_name")
#' }
#'

Create_Cluster_Annotation_File <- function(
  file_path = NULL,
  file_name = "cluster_annotation"
) {
  # Set file path is parameter is NULL
  if (is.null(x = file_path)) {
    dir_path <- getwd()
  } else {
    dir_path <- file_path
  }
  # Check directory path is exists
  if (!dir.exists(paths = dir_path)) {
    cli_abort(message = c("Target directory '{dir_path}' does not exist.",
                          "i" = "Please create directory or fix `file_path` and re-run function.")
    )
  }
  # Confirm no files with same name in the same directory path.
  full_path <- file.path(dir_path, paste0(file_name, ".csv"))
  if (file.exists(full_path)) {
    cli_abort(message = c("File with name {file_name} already exists in directory directory.",
                          "i" = "Please supply a different file_name.")
    )
  }
  # Save `Cluster_Annotation_Tibble`
  write.csv(Cluster_Annotation_Tibble(), full_path, row.names = F)
  cli_inform("Cluster annotation file created in: {dir_path}.")
}


#' Cluster Annotation Tibble
#'
#' Basic cluster annotation tibble for use in `Create_Cluster_Annotation_File`.  Contains columns:
#' "cluster" with values 0-32, and blank columns: "cell_type", "sub_type", "notes".
#'
#' @importFrom tibble tribble
#'
#' @keywords internal
#'
#' @noRd
#'

Cluster_Annotation_Tibble <- function(
) {
  annotation_tibble <- tribble(
                         ~cluster, ~cell_type, ~sub_type, ~notes,
                               0L,         "",        "",     "",
                               1L,         "",        "",     "",
                               2L,         "",        "",     "",
                               3L,         "",        "",     "",
                               4L,         "",        "",     "",
                               5L,         "",        "",     "",
                               6L,         "",        "",     "",
                               7L,         "",        "",     "",
                               8L,         "",        "",     "",
                               9L,         "",        "",     "",
                              10L,         "",        "",     "",
                              11L,         "",        "",     "",
                              12L,         "",        "",     "",
                              13L,         "",        "",     "",
                              14L,         "",        "",     "",
                              15L,         "",        "",     "",
                              16L,         "",        "",     "",
                              17L,         "",        "",     "",
                              18L,         "",        "",     "",
                              19L,         "",        "",     "",
                              20L,         "",        "",     "",
                              21L,         "",        "",     "",
                              22L,         "",        "",     "",
                              23L,         "",        "",     "",
                              24L,         "",        "",     "",
                              25L,         "",        "",     "",
                              26L,         "",        "",     "",
                              27L,         "",        "",     "",
                              28L,         "",        "",     "",
                              29L,         "",        "",     "",
                              30L,         "",        "",     "",
                              31L,         "",        "",     "",
                              32L,         "",        "",     ""
                         )
  return(annotation_tibble)
}


#' Pull cluster information from annotation csv file.
#'
#' shortcut filter and pull function compatible with annotation files created by `Create_Cluster_Annotation_File`
#' by default but also any other csv file.
#'
#' @param annotation name of the data.frame/tibble or path to CSV file containing cluster annotation.
#' @param cluster_name_col name of column containing cluster names/numbers (default is "cluster").
#' @param cell_type_col name of column contain the cell type annotation (default is "cell_type").
#'
#' @return a list of named vectors for every cell type in the `cell_type_col` column of the annotation table and
#' vectors new cluster names (for use with `Rename_Clusters` function or manual identity renaming).
#'
#' @import cli
#' @importFrom dplyr filter pull
#' @importFrom magrittr "%>%"
#' @importFrom utils read.csv
#'
#' @export
#'
#' @concept marker_annotation_util
#'
#' @examples
#' \dontrun{
#' # If pulling from a data.frame/tibble
#' cluster_annotation <- Pull_Cluster_Annotation(annotation = annotation_df,
#' cluster_name_col = "cluster", cell_type_col = "cell_type")
#'
#' # If pulling from csv file
#' cluster_annotation <- Pull_Cluster_Annotation(annotation = "file_path/file_name.csv",
#' cluster_name_col = "cluster", cell_type_col = "cell_type")
#' }
#'

Pull_Cluster_Annotation <- function(
  annotation = NULL,
  cluster_name_col = "cluster",
  cell_type_col = "cell_type"
) {
  # Check that annotation is in environment or a file that exists.
  if (!exists(x = deparse(expr = substitute(expr = annotation))) && !file.exists(annotation)) {
    cli_abort(message = "No file or environmental variable: {annotation} exists.")
  }
  # Read or specify annotation table
  if (exists(x = deparse(expr = substitute(expr = annotation)))) {
    annotation_table <- annotation
  } else {
    annotation_table <- read.csv(file = annotation, stringsAsFactors = FALSE)
  }

  # Check that cluster and cell type columns are present
  if (!cluster_name_col %in% colnames(x = annotation_table)) {
    cli_abort(message = "`cluster_name_col`: '{cluster_name_col}' not found in annotation data.frame.")
  }

  if (!cell_type_col %in% colnames(x = annotation_table)) {
    cli_abort(message = "`cell_type_col`: '{cell_type_col}' not found in annotation data.frame.")
  }

  # Create list elements per cluster
  cell_type_list <- unique(annotation_table[[cell_type_col]])
  cluster_annotation_list <- lapply(c(1:length(cell_type_list)), function(x){
    cluster <- annotation_table %>%
      filter(.data[[cell_type_col]] == cell_type_list[x]) %>%
      pull(cluster_name_col)
  })
  names(cluster_annotation_list) <- cell_type_list

  # Create list elements for renaming idents
  new_cluster_ids <- annotation_table %>%
    pull(cell_type_col)
  secondary_ids <- annotation_table[, 3]

  new_cluster_ids_list <- list(new_cluster_ids)
  secondary_ids_list <- list(secondary_ids)
  # Name the new cluster ids list
  names(new_cluster_ids_list) <- "new_cluster_idents"
  names(secondary_ids_list) <- colnames(annotation_table)[[3]]

  # Combine and return both lists as single list
  final_cluster_annotation_list <- c(cluster_annotation_list, new_cluster_ids_list, secondary_ids_list)
  return(final_cluster_annotation_list)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### PROJECT ORGANIZATION ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Setup project directory structure
#'
#' Create reproducible project directory organization when initiating a new analysis.
#'
#' @param custom_dir_file file to file containing desired directory structure.  Default is NULL and
#' will provide generic built-in directory structure.
#' @param cluster_annotation_path path to place cluster annotation file using \code{\link{Create_Cluster_Annotation_File}}.
#' @param cluster_annotation_file_name name to use for annotation file if created (optional).
#'
#' @import cli
#' @importFrom data.table fread
#' @importFrom dplyr pull
#' @importFrom magrittr "%>%"
#'
#' @export
#'
#' @return no return value.  Creates system folders.
#'
#' @concept organization_util
#'
#' @examples
#' \dontrun{
#' # If using built-in directory structure.
#' Setup_scRNAseq_Project()
#' }
#'

Setup_scRNAseq_Project <- function(
  custom_dir_file = NULL,
  cluster_annotation_path = NULL,
  cluster_annotation_file_name = "cluster_annotation.csv"
) {

  if (is.null(x = custom_dir_file)) {
    # File paths setup
    output_dirs <- list(
      scripts_path = "01_scripts/",
      data_path = "02_raw_data/",
      metadata_path = "03_meta_data/",
      object_path = "04_data_objects/",
      plots_path = "05_plots/",
      plots_qc_path = "05_plots/01_QC_plots/",
      plots_Rd01_path = "05_plots/02_Round01_plots/",
      plots_Rd02_path = "05_plots/03_Round02_plots/",
      annotation_path = "06_cluster_annotation/",
      outputs_path = "07_csv_outputs/",
      final_plots_path = "08_final_plots_for_figures/"
    )
  } else {
    # Check custom paths file exists
    if (!file.exists(custom_dir_file)) {
      cli_abort(message = "`custom_dir_file` not found.  Please check file path and name provided.")
    }

    # Read file and create directory list
    output_dirs <- fread(input = custom_dir_file, data.table = FALSE) %>%
      pull() %>%
      as.list()
  }

  # Check for directories and create new ones
  lapply(output_dirs, function(dir_path){
    if (!dir.exists(dir_path)){
      dir.create(dir_path)
    } else {
      cli_warn(message = "The directory {dir_path} aleady exists.  No new directory created.")
    }
  })

  # Create annotation file if desired.
  if (!is.null(x = cluster_annotation_path)) {
    # Add Cluster annotation file
    Create_Cluster_Annotation_File(file_path = cluster_annotation_path, file_name = cluster_annotation_file_name)
  }

  # Print completion message
  cli_inform(message = "scRNA-seq R project setup complete.")
}


#' Copy folder to GCP bucket from R Console
#'
#' Run command from R console without moving to terminal to copy folder to GCP bucket
#'
#' @param folder_file_path folder to be copied to GCP bucket.
#' @param gcp_bucket_path GCP bucket path to copy to files.
#'
#' @import cli
#'
#' @export
#'
#' @return No return value.  Performs system copy to GCP bucket.
#'
#' @concept organization_util
#'
#' @examples
#' \dontrun{
#' Copy_To_GCP(folder_file_path = "plots/", gcp_bucket_path = "gs://bucket_name_and_folder_path")
#' }
#'

Copy_To_GCP <- function(
  folder_file_path,
  gcp_bucket_path
) {
  # Check directory path is exists
  if (!dir.exists(paths = folder_file_path)) {
    cli_abort(message = c("Target directory '{folder_file_path}' does not exist.",
                          "i" = "Please create directory or fix `file_path` and re-run function.")
    )
  }

  # Copy files
  system(paste0("gsutil -m cp -r ", folder_file_path, " ", gcp_bucket_path))
}


#' Copy folder from GCP bucket from R Console
#'
#' Run command from R console without moving to terminal to copy folder from GCP bucket to local storage
#'
#' @param folder_file_path folder to be copied to GCP bucket.
#' @param gcp_bucket_path GCP bucket path to copy to files.
#'
#' @import cli
#'
#' @export
#'
#' @return No return value.  Performs system copy from GCP bucket.
#'
#' @concept organization_util
#'
#' @examples
#' \dontrun{
#' Copy_From_GCP(folder_file_path = "plots/", gcp_bucket_path = "gs://bucket_name_and_folder_path")
#' }
#'

Copy_From_GCP <- function(
  folder_file_path,
  gcp_bucket_path
) {
  # Check directory path is exists
  if (!dir.exists(paths = folder_file_path)) {
    cli_abort(message = c("Target directory '{folder_file_path}' does not exist.",
                          "i" = "Please create directory or fix `file_path` and re-run function.")
    )
  }

  # Copy files
  system(paste0("gsutil -m cp -r ", gcp_bucket_path, " ", folder_file_path))
}
