#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### OBJECT HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Check if genes/features are present
#'
#' Check if genes are present in object and return vector of found genes.  Return warning messages for
#' genes not found.
#'
#' @param data Name of input data.  Currently only data of classes: Seurat, liger, data.frame,
#' dgCMatrix, dgTMatrix, tibble are accepted.  Gene_IDs must be present in rownames of the data.
#' @param features vector of features to check.
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
#' Default is NULL which will check against features from all assays present.
#'
#' @import cli
#' @importFrom purrr reduce
#' @importFrom SeuratObject Features
#' @importFrom stringr str_to_upper str_to_sentence
#'
#' @return A list of length 3 containing 1) found features, 2) not found features, 3) features found if
#' case was modified.
#'
#' @export
#'
#' @concept check_util
#'
#' @examples
#' \dontrun{
#' features <- Feature_Present(data = obj_name, features = DEG_list, print_msg = TRUE,
#' case_check = TRUE)
#' found_features <- features[[1]]
#' }
#'

Feature_Present <- function(
    data,
    features,
    case_check = TRUE,
    case_check_msg = TRUE,
    print_msg = TRUE,
    omit_warn = TRUE,
    return_none = FALSE,
    seurat_assay = NULL
) {
  # Check object type
  # Seurat
  accepted_types <- c("data.frame", "dgCMatrix", "dgTMatrix", "tibble", "ligerDataset")
  if (inherits(x = data, what = "Seurat")) {
    # set assay (if null set to active assay)
    assays_present <- seurat_assay %||% Assays(object = data)

    possible_features <- lapply(assays_present, function(j) {
      Features(x = data, assay = j)
    })

    possible_features <- unlist(possible_features)
  } else if ((class(x = data)[[1]] == "liger")) {
      possible_features <- Features(x = data, by_dataset = FALSE)
  } else if ((class(x = data) %in% accepted_types)) {
    possible_features <- rownames(x = data)
  } else {
    all_accepted <- c(accepted_types, "Seurat", "liger")
    cli_abort(message = c("Input data is currently accepted only in the following formats:",
                          "i" = "{.field {glue_collapse_scCustom(input_string = all_accepted, and = FALSE)}}.")
    )
  }

  # If any features not found
  if (any(!features %in% possible_features)) {
    bad_features <- features[!features %in% possible_features]
    found_features <- features[features %in% possible_features]
    if (length(x = found_features) == 0) {
      if (isTRUE(x = return_none)) {
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
    if (length(x = bad_features) > 0 && isTRUE(x = omit_warn)) {
      cli_warn(message = c("The following features were omitted as they were not found:",
                           "i" = "{.field {glue_collapse_scCustom(input_string = bad_features, and = TRUE)}}")
      )
    }

    # Check if features found if case is changed.
    if (isTRUE(x = case_check)) {
      upper_bad_features <- str_to_upper(string = bad_features)
      upper_found_features <- upper_bad_features[upper_bad_features %in% possible_features]

      sentence_bad_features <- str_to_sentence(string = bad_features)
      sentence_found_features <- sentence_bad_features[sentence_bad_features %in% possible_features]

      # Combine case check
      wrong_case_found_features <- c(upper_found_features, sentence_found_features)

      # Additional messages if found.
      if (length(x = wrong_case_found_features) > 0) {
        if (isTRUE(x = case_check_msg)) {
          cli_warn(message = c("NOTE: However, the following features were found: {.field {glue_collapse_scCustom(input_string = wrong_case_found_features, and = TRUE)}}",
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
  if (isTRUE(x = print_msg)) {
    cli_inform(message = "All features present.")
  }

  # Return full input gene list.
  # Combine into list and return
  feature_list <- list(
    found_features = features,
    bad_features = NULL,
    wrong_case_found_features = NULL
  )
  return(feature_list)
}


#' Check for alternate case features
#'
#' Checks Seurat object for the presence of features with the same spelling but alternate case.
#'
#' @param seurat_object Seurat object name.
#' @param gene_list vector of genes to check.
#' @param case_check_msg logical. Whether to print message to console if alternate case features are
#' found in addition to inclusion in returned list.  Default is TRUE.
#' @param return_features logical. Whether to return vector of alternate case features.  Default is TRUE.
#' @param assay Name of assay to pull feature names from. If NULL will use the result of `DefaultAssay(seurat_object)`.
#'
#' @import cli
#' @importFrom SeuratObject Features
#' @importFrom stringr str_to_sentence str_to_upper
#'
#' @return If features found returns vector of found alternate case features and prints message depending on
#' parameters specified.
#'
#' @export
#'
#' @concept check_util
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
  return_features = TRUE,
  assay = NULL
) {
  # set assay (if null set to active assay)
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # get all features
  possible_features <- Features(x = seurat_object, assay = assay)

  upper_bad_features <- str_to_upper(string = gene_list)
  upper_found_features <- upper_bad_features[upper_bad_features %in% possible_features]

  sentence_bad_features <- str_to_sentence(string = gene_list)
  sentence_found_features <- sentence_bad_features[sentence_bad_features %in% possible_features]

  # Combine case check
  wrong_case_found_features <- c(upper_found_features, sentence_found_features)

  # Additional messages if found.
  if (length(x = wrong_case_found_features) > 0) {
    if (isTRUE(x = case_check_msg)) {
      cli_inform(message = c("{col_cyan('*NOTE*')}: However, the following features were found: {.field {glue_collapse_scCustom(input_string = wrong_case_found_features, and = TRUE)}}",
                             "i" = "Please check intended case of features provided.")
      )
    }
    if (isTRUE(x = return_features)) {
      return(wrong_case_found_features)
    }
  }
}


#' Check if meta data are present
#'
#' Check if meta data columns are present in object and return vector of found columns
#' Return warning messages for meta data columns not found.
#'
#' @param object Seurat or Liger object name.
#' @param meta_col_names vector of column names to check.
#' @param print_msg logical. Whether message should be printed if all features are found.  Default is TRUE.
#' @param omit_warn logical. Whether to print message about features that are not found in current object. Default is TRUE.
#' @param return_none logical. Whether list of found vs. bad features should still be returned if no
#' `meta_col_names` are found.  Default is FALSE.
#'
#' @return vector of meta data columns that are present
#'
#' @import cli
#'
#' @export
#'
#' @concept check_util
#'
#' @examples
#' \dontrun{
#' meta_variables <- Meta_Present(object = obj_name, meta_col_names = "percent_mito", print_msg = TRUE)
#' }
#'

Meta_Present <- function(
  object,
  meta_col_names,
  print_msg = TRUE,
  omit_warn = TRUE,
  return_none = FALSE
) {
  # Set possible variables based on object type
  possible_features <- colnames(x = Fetch_Meta(object = object))

  # If any features not found
  if (any(!meta_col_names %in% possible_features)) {
    bad_meta <- meta_col_names[!meta_col_names %in% possible_features]
    found_meta <- meta_col_names[meta_col_names %in% possible_features]

    if (isFALSE(return_none)) {
      if (length(x = found_meta) < 1) {
        cli_abort(message = c("No meta data columns found.",
                              "i" = "The following meta data columns were not found: {.field {glue_collapse_scCustom(input_string = bad_meta, and = TRUE)}}")
        )
      }
    }

    # Return message of features not found
    if (length(x = bad_meta) > 0 && isTRUE(x = omit_warn)) {
      cli_warn(message = c("The following meta data columns were omitted as they were not found:",
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
  if (isTRUE(x = print_msg)) {
    cli_inform(message = "All meta data columns present.")
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
#' @concept check_util
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

  colnames(x = all_numeric) <- "Is_Numeric"

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


#' Check if reduction loadings are present
#'
#' Check if reduction loadings are present in object and return vector of found loading names.  Return
#' warning messages for genes not found.
#'
#' @param seurat_object object name.
#' @param reduction_names vector of genes to check.
#' @param print_msg logical. Whether message should be printed if all features are found.  Default is TRUE.
#' @param omit_warn logical. Whether to print message about features that are not found in current object.
#'  Default is TRUE.
#' @param return_none logical. Whether list of found vs. bad features should still be returned if no
#' features are found.  Default is FALSE.
#'
#' @importFrom purrr reduce
#' @importFrom stringr str_to_upper str_to_sentence
#'
#' @return A list of length 3 containing 1) found features, 2) not found features.
#'
#' @export
#'
#' @concept check_util
#'
#' @examples
#' \dontrun{
#' reductions <- Reduction_Loading_Present(seurat_object = obj_name, reduction_name = "PC_1")
#' found_reductions <- reductions[[1]]
#' }
#'

Reduction_Loading_Present <- function(
    seurat_object,
    reduction_names,
    print_msg = TRUE,
    omit_warn = TRUE,
    return_none = FALSE
) {
  # If no reductions are present
  if (length(x = seurat_object@reductions) == 0) {
    if (isTRUE(x = return_none)) {
      # Combine into list and return
      feature_list <- list(
        found_features = NULL,
        bad_features = NULL
      )
      return(feature_list)
    } else {
      cli_abort(message ="No reductions present in object.")
    }
  }

  # Get all reduction names
  possible_reduction_names <- unlist(x = lapply(1:length(x = seurat_object@reductions), function(z) {
    names <- names(x = seurat_object@reductions[[z]])
  })
  )

  # If any features not found
  if (any(!reduction_names %in% possible_reduction_names)) {
    bad_reductions <- reduction_names[!reduction_names %in% possible_reduction_names]
    found_reductions <- reduction_names[reduction_names %in% possible_reduction_names]
    if (length(x = found_reductions) == 0) {
      if (isTRUE(x = return_none)) {
        # Combine into list and return
        reduction_list <- list(
          found_reductions = NULL,
          bad_reductions = bad_reductions
        )
        return(reduction_list)
      } else {
        cli_abort(message ="No requested features found.")
      }
    }

    # Return message of features not found
    if (length(x = bad_reductions) > 0 && isTRUE(x = omit_warn)) {
      cli_warn(message = c("The following features were omitted as they were not found:",
                           "i" = "{.field {glue_collapse_scCustom(input_string = bad_features, and = TRUE)}}")
      )
    }

    # Combine into list and return
    reduction_list <- list(
      found_reductions = found_reductions,
      bad_reductions = bad_reductions
    )
    return(reduction_list)
  }

  # Print all found message if TRUE
  if (isTRUE(x = print_msg)) {
    cli_inform(message = "All features present.")
  }

  # Return full input gene list.
  # Combine into list and return
  reduction_list <- list(
    found_reductions = reduction_names,
    bad_reductions = NULL
  )
  return(reduction_list)
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
#' @references Original function is part of LIGER package
#' \url{https://github.com/welch-lab/liger/blob/master/R/mergeObject.R} (License: GPL-3).
#' Function was modified for use in scCustomize (add progress bar, prefix vs. suffix, and delimiter options).
#'
#' @import cli
#' @import Matrix
#' @importFrom dplyr intersect
#' @importFrom magrittr "%>%"
#'
#' @return A sparse Matrix
#'
#' @export
#'
#' @concept read_merge_util
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

  if (isTRUE(x = duplicated_barcodes) && is.null(x = add_cell_ids)) {
    cli_abort(message = c("There are overlapping cell barcodes present in the input matrices.",
                          "i" = "Please provide prefixes/suffixes to {.code add_cell_ids} parameter to make unique.")
    )
  }

  # Check right number of suffix/prefix ids are provided
  if (!is.null(x = add_cell_ids) && length(x = add_cell_ids) != length(x = matrix_list)) {
    cli_abort(message = "The number of prefixes in {.code add_cell_ids} must be equal to the number of matrices supplied to {.code matrix_list}.")
  }

  if (!is.null(x = add_cell_ids)) {
    # check barcodes will be unique after adding prefixes/suffixes
    all_names <- lapply(1:length(x = matrix_list), function(i){
      cell_names <- colnames(x = matrix_list[[i]])
    })

    new_names <- lapply(X = 1:length(x = matrix_list), function(x){
      colnames(x = matrix_list[[x]]) <- paste0(add_cell_ids[x], cell_id_delimiter, colnames(x = matrix_list[[x]]))
    })

    are_duplicates <- unlist(x = new_names) %>%
      duplicated() %>%
      any()

    if (isTRUE(x = are_duplicates)) {
      cli_abort(message = c("Supplied {.code add_cell_ids} will result in overlapping barcodes names if provided cell prefixes/suffixes are not unique.",
                            "i" = "Please change and re-run.")
      )
    }
  }

  # Use summary to convert the sparse matrices into three-column indexes where i are the
  # row numbers, j are the column numbers, and x are the nonzero entries
  col_offset <- 0
  allGenes <- unique(x = unlist(x = lapply(matrix_list, rownames)))
  allCells <- c()
  cli_inform(message = "{.field Preparing & merging matrices.}")
  pb <- txtProgressBar(min = 0, max = length(x = matrix_list), style = 3, file = stderr())
  for (i in 1:length(x = matrix_list)) {
    curr <- matrix_list[[i]]
    curr_s <- summary(curr)

    # Now, alter the indexes so that the two 3-column matrices can be properly merged.
    # First, make the current and full column numbers non-overlapping.
    curr_s[, 2] <- curr_s[, 2] + col_offset

    # Update full cell names
    if (!is.null(x = add_cell_ids)) {
      if (isTRUE(x = prefix)) {
        cellnames <- paste0(add_cell_ids [i], cell_id_delimiter, colnames(x = curr))
      } else {
        cellnames <- paste0(colnames(x = curr), cell_id_delimiter, add_cell_ids [i])
      }
    } else {
      cellnames <- colnames(x = curr)
    }
    allCells <- c(allCells, cellnames)

    # Next, change the row (gene) indexes so that they index on the union of the gene sets,
    # so that proper merging can occur.
    idx <- match(x = rownames(x = curr), allGenes)
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
  cli_inform(message = "{.field Creating final sparse matrix.}")
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


#' Extract multi-modal data into list by modality
#'
#' Reorganize multi-modal data after import with `Read10X()` or scCustomize read functions.
#' Organizes sub-lists by data modality instead of by sample.
#'
#' @param matrix_list list of matrices to split by modality
#'
#' @return list of lists, with one sublist per data modality.  Sub-list contain 1 matrix entry per sample
#'
#' @import cli
#'
#' @export
#'
#' @concept read_merge_util
#'
#' @examples
#' \dontrun{
#' multi_mat <- Read10X(...)
#' new_multi_mat <- Extract_Modality(matrix_list = multi_mat)
#' }
#'

Extract_Modality <- function(
    matrix_list
) {
  modality_names <- names(x = matrix_list[[1]])

  unlist_mat <- unlist(x = matrix_list)

  index_list <- lapply(1:length(x = modality_names), function(x) {
    modality_index <- grep(x = names(x = unlist_mat), pattern = modality_names[x])
  })

  split_list <- lapply(1:length(x = modality_names), function(i) {
    modality_list <- unlist_mat[index_list[[i]]]
    sample_name <- gsub(pattern = paste0("_.", modality_names[i]), x = names(x = modality_list), replacement = "")
    names(x = modality_list) <- sample_name
    return(modality_list)
  })

  names(x = split_list) <- modality_names
  return(split_list)
}


#' Merge a list of Sparse Matrices contain multi-modal data.
#'
#' Enables easy merge of a list of sparse matrices for multi-modal data.
#'
#' @param matrix_list list of matrices to merge.
#' @param add_cell_ids a vector of sample ids to add as prefix to cell barcode during merge.
#' @param prefix logical.  Whether `add_cell_ids` should be added as prefix to current cell barcodes/names
#' or as suffix to current cell barcodes/names.  Default is TRUE, add as prefix.
#' @param cell_id_delimiter The delimiter to use when adding cell id prefix/suffix.  Default is "_".
#'
#' @import cli
#'
#' @return A list containing one sparse matrix for each modality
#'
#' @export
#'
#' @concept read_merge_util
#'
#' @examples
#' \dontrun{
#' data_list <- Read10X_GEO(...)
#' merged_list <- Merge_Sparse_Multimodal_All(matrix_list = data_list, add_cell_ids = names(data_list),
#' prefix = TRUE, cell_id_delimiter = "_")
#' }
#'

Merge_Sparse_Multimodal_All <- function(
    matrix_list,
    add_cell_ids = NULL,
    prefix = TRUE,
    cell_id_delimiter = "_"
) {
  # Check matrix_list is list of lists
  if (!inherits(x = matrix_list[[1]], what = "list")) {
    cli_abort(message = "{.code matrix_list} is not multimodal, please use {.field Merge_Sparse_Data_All}.")
  }

  # Extract matrices
  mat_list <- Extract_Modality(matrix_list = matrix_list)

  # Merge and return
  modality_names <- names(x = mat_list)

  merged_list <- lapply(1:length(x = modality_names), function(x) {
    cli_inform(message = "Merging {.val {modality_names[x]}} matrices.")
    merged <- Merge_Sparse_Data_All(matrix_list = mat_list[[x]], add_cell_ids = add_cell_ids, prefix = prefix, cell_id_delimiter = cell_id_delimiter)
  })

  names(x = merged_list) <- modality_names

  return(merged_list)
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
#' @references Re-implementing `CheckMatrix` only for sparse matrices with modified warning messages.  Original function from SeuratObject \url{https://github.com/satijalab/seurat-object/blob/9c0eda946e162d8595696e5280a6ecda6284db39/R/utils.R#L625-L650} (License: MIT).
#'
#' @export
#'
#' @concept check_util
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
#' @concept barcode_util
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
      cli_abort(message = "{.code current_suffix} must either be single value or a vector of values equal in length to number of datasets in {.code data}.")
    }

    # Check new suffix
    if (length(x = new_suffix) == 1) {
      new_suffix <- rep(new_suffix, length(x = data))
    }
    if (length(x = new_suffix) == length(x = data)) {
      new_suffix <- new_suffix
    }
    if (length(x = new_suffix) != length(x = data) && length(x = new_suffix) != 1) {
      cli_abort(message = "{.code new_suffix} must either be single value or a vector of values equal in length to number of datasets in {.code data}.")
    }

    # Is current suffix found in all cell names
    check_suffixes <- sapply(1:length(x = data), FUN = function(j){
      all(grepl(pattern = current_suffix_regexp[[j]], x = current_cell_names[[j]]))
    })

    if (all(check_suffixes) != TRUE) {
      cli_abort(message = c("One or more {.code current_suffixes} do not match cell names in data.",
                            "i" = "Check inputs.")
                )
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
    names(x = data_mod) <- names(x = data)
    return(data_mod)

  } else {
    # for data.frames and individual matrices
    current_cell_names <- colnames(x = data)
    current_suffix_regexp <- paste0("\\", current_suffix, "$")
    # Is current suffix found in cell names
    if (all(grepl(pattern = current_suffix_regexp, x = current_cell_names)) == FALSE) {
      cli_abort(message = "Supplied {.code current_suffix}: {.field {current_suffix}} was not found in the cell names of data provided.")
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
#' @concept barcode_util
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
    check_suffix_delim <- sapply(1:length(x = data), FUN = function(j){
      all(grepl(pattern = current_delim, x = current_cell_names[[j]], fixed = TRUE))
    })

    if (all(check_suffix_delim) != TRUE) {
      cli_abort(message = c("One or more {.code current_delim} do not match cell names in data.",
                            "i" = "Check inputs.")
      )
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
    names(x = data_mod) <- names(x = data)
    return(data_mod)
  } else {
    # for data.frames and individual matrices
    current_cell_names <- colnames(x = data)

    # Is current suffix found in cell names
    if (all(grepl(pattern = current_delim, x = current_cell_names, fixed = TRUE)) != TRUE) {
      cli_abort(message = "Supplied {.code current_delim}: {.field {current_delim}} was not found in the cell names of data provided.")
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
#' @concept barcode_util
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
    check_prefix_delim <- sapply(1:length(x = data), FUN = function(j){
      all(grepl(pattern = current_delim, x = current_cell_names[[j]], fixed = TRUE))
    })

    if (all(check_prefix_delim) != TRUE) {
      cli_abort(message = c("One or more {.code current_delim} do not match cell names in data.",
                            "i" = "Check inputs.")
      )
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
    names(x = data_mod) <- names(x = data)
    return(data_mod)
  } else {
    # for data.frames and individual matrices
    current_cell_names <- colnames(x = data)

    # Is current prefix found in cell names
    if (all(grepl(pattern = current_delim, x = current_cell_names, fixed = TRUE)) != TRUE) {
      cli_abort(message = "Supplied {.code current_delim}: {.field {current_delim}} was not found in the cell names of data provided.")
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
#' @concept barcode_util
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
    check_prefix_delim <- sapply(1:length(x = data), FUN = function(j){
      all(grepl(pattern = current_delim, x = current_cell_names[[j]], fixed = TRUE))
    })

    if (all(check_prefix_delim) != TRUE) {
      cli_abort(message = c("One or more {.code current_delim} do not match cell names in data.",
                            "i" = "Check inputs.")
      )
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
    names(x = data_mod) <- names(x = data)
    return(data_mod)
  } else {
    # for data.frames and individual matrices
    current_cell_names <- colnames(x = data)

    # Is current prefix found in cell names
    if (all(grepl(pattern = current_delim, x = current_cell_names, fixed = TRUE)) != TRUE) {
      cli_abort(message = "Supplied {.code current_delim}: {.field {current_delim}} was not found in the cell names of data provided.")
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
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("{.val pct_diff} column already present in {.code marker_dataframe}: {.val {df_name}}.",
                            "i" = "To overwrite previous results set `overwrite = TRUE`.")
      )
    } else {
      cli_inform(message = c("{.val pct_diff} column already present in {.code marker_dataframe}: {.val {df_name}}.",
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
    cli_abort(message = "{.code rank_by}: {.val {rank_by}} not found in column names of {.code marker_dataframe}.")
  }

  # Check grouping factor in marker data.frame
  if (!is.null(x = group_by)) {
    if (!group_by %in% colnames(x = marker_dataframe)) {
      cli_abort(message = "{.code group_by}: {.val {group_by}} not found in column names of {.code marker_dataframe}.")
    }
  }

  # Check gene column is present
  if (!gene_column %in% colnames(x = marker_dataframe) && isFALSE(x = gene_rownames_to_column)) {
    cli_abort(message = c("{.code gene_column}: '{gene_column}' not found in column names of {.code marker_dataframe}.",
                          "i" = "Set {.code gene_rownames_to_column} to move genes from rownames to column.")
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

  if (isTRUE(x = gene_rownames_to_column)) {
    filtered_markers <- filtered_markers %>%
      rownames_to_column(gene_column)
  }

  # return data.frame
  if (isTRUE(x = data_frame)) {
    return(filtered_markers)
  }

  # pull gene list
  gene_list <- filtered_markers[[gene_column]]

  # should gene list be named
  # check naming
  if (isTRUE(x = named_vector) && is.null(x = group_by)) {
    cli_warn(message = c("Cannot return named vector if {.code group_by} is NULL.",
                         "i" = "Returning unnamed vector.")
    )
  }

  if (isTRUE(x = named_vector) && !is.null(x = group_by)) {
    if (isTRUE(x = make_unique)) {
      cli_abort(message = "Cannot return unique list if {.code named_vector = TRUE}.")
    }
    names(x = gene_list) <- filtered_markers[[group_by]]
    return(gene_list)
  }

  # make unique
  if (isTRUE(x = make_unique)) {
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
    if (file_path == "") {
      dir_path <- getwd()
    } else {
      dir_path <- file_path
    }
  }
  # Check directory path is exists
  if (!dir.exists(paths = dir_path)) {
    cli_abort(message = c("Target directory {.val {dir_path}} does not exist.",
                          "i" = "Please create directory or fix {.code file_path} and re-run function.")
    )
  }

  # Check extension
  file_ext <- grep(x = file_name, pattern = ".csv$")

  if (length(x = file_ext) == 0) {
    file_name <- paste0(file_name, ".csv")
  }


  # Confirm no files with same name in the same directory path.
  full_path <- file.path(dir_path, file_name)
  if (file.exists(full_path)) {
    cli_abort(message = c("File with name {.val {file_name}} already exists in directory directory.",
                          "i" = "Please supply a different {.code file_name}.")
    )
  }
  # Save `Cluster_Annotation_Tibble`
  write.csv(Cluster_Annotation_Tibble(), full_path, row.names = F)
  cli_inform("Cluster annotation file created in: {.val {dir_path}}.")
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
    cli_abort(message = "No file or environmental variable: {.field {annotation}} exists.")
  }
  # Read or specify annotation table
  if (exists(x = deparse(expr = substitute(expr = annotation)))) {
    annotation_table <- annotation
  } else {
    annotation_table <- read.csv(file = annotation, stringsAsFactors = FALSE)
  }

  # Check that cluster and cell type columns are present
  if (!cluster_name_col %in% colnames(x = annotation_table)) {
    cli_abort(message = "{.code cluster_name_col}: {.val {cluster_name_col}} not found in annotation data.frame.")
  }

  if (!cell_type_col %in% colnames(x = annotation_table)) {
    cli_abort(message = "{.code cell_type_col}: {.val {cell_type_col}} not found in annotation data.frame.")
  }

  # Create list elements per cluster
  cell_type_list <- unique(x = annotation_table[[cell_type_col]])
  cluster_annotation_list <- lapply(c(1:length(x = cell_type_list)), function(x){
    cluster <- annotation_table %>%
      filter(.data[[cell_type_col]] == cell_type_list[x]) %>%
      pull(cluster_name_col)
  })
  names(x = cluster_annotation_list) <- cell_type_list

  # Create list elements for renaming idents
  new_cluster_ids <- annotation_table %>%
    pull(cell_type_col)
  secondary_ids <- annotation_table[, 3]

  new_cluster_ids_list <- list(new_cluster_ids)
  secondary_ids_list <- list(secondary_ids)
  # Name the new cluster ids list
  names(x = new_cluster_ids_list) <- "new_cluster_idents"
  names(x = secondary_ids_list) <- colnames(x = annotation_table)[[3]]

  # Combine and return both lists as single list
  final_cluster_annotation_list <- c(cluster_annotation_list, new_cluster_ids_list, secondary_ids_list)
  return(final_cluster_annotation_list)
}


#' @param new_idents vector of new cluster names.  Must be equal to the length of current default identity
#' of Object.  Will accept named vector (with old idents as names) or will name the new_idents vector internally.
#' @param meta_col_name `r lifecycle::badge("soft-deprecated")`. See `old_ident_name`.
#' @param old_ident_name optional, name to use for storing current object idents in object meta data slot.
#' @param new_ident_name optional, name to use for storing new object idents in object meta data slot.
#' @param overwrite logical, whether to overwrite columns in object meta data slot. if they have same
#' names as `old_ident_name` and/or `new_ident_name`.
#'
#' @method Rename_Clusters Seurat
#'
#' @import cli
#' @importFrom lifecycle deprecated
#'
#' @rdname Rename_Clusters
#' @export
#'
#' @concept marker_annotation_util
#'
#' @examples
#' \dontrun{
#' obj <- Rename_Clusters(seurat_object = obj_name, new_idents = new_idents_vec,
#' old_ident_name = "Seurat_Idents_Round01", new_ident_name = "Round01_Res0.6_Idents")
#' }
#'

Rename_Clusters.Seurat <- function(
  object,
  new_idents,
  old_ident_name = NULL,
  new_ident_name = NULL,
  meta_col_name = deprecated(),
  overwrite = FALSE,
  ...
) {
  # Deprecation warning
  if (lifecycle::is_present(meta_col_name)) {
    lifecycle::deprecate_stop(when = "2.2.0",
                              what = "Rename_Clusters(meta_col_name)",
                              with = "Rename_Clusters(old_ident_name)",
                              details = c("i" = "To store old idents please provide name to `old_ident_name`",
                                          "i" = "To store new idents please provide name to `new_ident_name`")
    )
  }

  # Check Seurat
  Is_Seurat(seurat_object = object)

  # check old ident name
  if (!is.null(x = old_ident_name)) {
    if (old_ident_name %in% colnames(x = object@meta.data)) {
      if (isFALSE(x = overwrite)) {
        cli_abort(message = c("The {.code old_ident_name}: {.field {old_ident_name}} is already a column in meta.data",
                              "i" = "To overwrite current meta.data column set {.code overwrite = TRUE}."))
      } else {
        cli_inform(message = "Overwriting old meta.data column: {.field {old_ident_name}} as {.code overwrite = TRUE}")

      }
    } else {
      object[[old_ident_name]] <- Idents(object = object)
    }
  }

  # check new ident name
  if (!is.null(x = new_ident_name) && new_ident_name %in% colnames(x = object@meta.data)) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("The {.code new_ident_name}: {.field {new_ident_name}} is already a column in meta.data",
                            "i" = "To overwrite current meta.data column set {.code overwrite = TRUE}."))
    } else {
      cli_inform(message = "Overwriting new meta.data column: {.field {new_ident_name}} as {.code overwrite = TRUE}")
    }
  }

  # Check equivalent lengths
  if (length(x = new_idents) != length(x = levels(x = object))) {
    cli_abort(message = c("Length of {.code new_idents} must be equal to the number of active.idents in Seurat Object.",
                          "i" = "{.code new_idents} length: {.field {length(x = new_idents)}} Object@active.idents length: {.field {length(x = levels(x = object))}}.")
    )
  }

  # Name the new idents vector
  if (is.null(x = names(x = new_idents))) {
    names(x = new_idents) <- levels(x = object)
  }

  # If named check that names are right length
  if (!is.null(x = names(x = new_idents)) && length(x = unique(x = names(x = new_idents))) != length(x = levels(x = object))) {
    cli_abort(message = c("The number of unique names for {.code new idents} is not equal to number of active.idents.",
                          "i" = "names(new_idents) length: {.field {length(x = unique(x = names(x = new_idents)))} Object@active.idents length: {length(x = levels(x = object))}}.")
    )
  }

  # Add new idents
  object <- RenameIdents(object = object, new_idents)

  # Add new ident to meta.data information if desired
  if (!is.null(x = new_ident_name)) {
    object[[new_ident_name]] <- Idents(object = object)
  }

  # return object
  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GENERAL HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Split vector into list
#'
#' Splits vector into chunks of x sizes
#'
#' @param x vector to split
#' @param chunk_size size of chunks for vector to be split into, default is NULL.  Only valid if
#' `num_chunk` is NULL.
#' @param num_chunk number of chunks to split the vector into, default is NULL.  Only valid if
#' `chunk_size` is NULL.
#' @param verbose logical, print details of vector and split, default is FALSE.
#'
#' @return list with vector of X length
#'
#' @import cli
#'
#' @export
#'
#' @references Base code from stackoverflow post:
#' \url{https://stackoverflow.com/a/3321659/15568251}
#'
#' @concept misc_util
#'
#' @examples
#' vector <- c("gene1", "gene2", "gene3", "gene4", "gene5", "gene6")
#'
#' vector_list <- Split_Vector(x = vector, chunk_size = 3)
#'

Split_Vector <- function(
    x,
    chunk_size = NULL,
    num_chunk = NULL,
    verbose = FALSE
) {
  if (!is.null(x = chunk_size) && !is.null(x = num_chunk)) {
    cli_abort(message = "Cannot specify both {.code chunk_size} and {.code num_chunk}, use one or the other.")
  }

  # set chunk size
  chunk_size <- chunk_size %||% (length(x = x) / num_chunk)

  # Split vector
  vector_list <- split(x, ceiling(x = seq_along(x)/chunk_size))

  # Report info
  if (isTRUE(x = verbose)) {
    cli_inform(message = c("Original vector length: ({.field {length(x = x)}}).",
                           "Split into {.field {length(x = vector_list)}} vectors of {.field {chunk_size}} items." ))
  }

  # return list
  return(vector_list)
}


#' Create sequence with zeros
#'
#' Create sequences of numbers like `seq()` or `seq_len()` but with zeros prefixed to
#' keep numerical order
#'
#' @param seq_length a seqeunce or numbers of numbers to create sequence.
#' Users can provide sequence (1:XX) or number of values to add in sequence (will
#' be used as second number in `seq_len`; 1:XX).
#' @param num_zeros number of zeros to prefix sequence, default is  (e.g, 01, 02, 03, ...)
#'
#' @return vector of numbers in sequence
#'
#' @import cli
#' @importFrom stringr str_pad
#'
#' @export
#'
#' @references Base code from stackoverflow post:
#' \url{https://stackoverflow.com/a/38825614}
#'
#' @concept misc_util
#'
#' @examples
#' # Using sequence
#' new_seq <- seq_zeros(seq_length = 1:15, num_zeros = 1)
#' new_seq
#'
#' # Using number
#' new_seq <- seq_zeros(seq_length = 15, num_zeros = 1)
#' new_seq
#'
#' # Sequence with 2 zeros
#' new_seq <- seq_zeros(seq_length = 1:15, num_zeros = 2)
#' new_seq
#'

seq_zeros <- function(
    seq_length,
    num_zeros = NULL
) {
  # Check seq_length
  if (is.null(x = num_zeros)) {
    if (length(x = seq_length) == 1) {
      if (nchar(x = seq_length) == 1) {
        num_zeros <- 1
      } else {
        num_zeros <- nchar(x = seq_length) - 1
      }
    } else {
      if (nchar(x = tail(x = seq_length, n = 1)) == 1) {
        num_zeros <- 1
      } else {
        num_zeros <- nchar(x = tail(x = seq_length, n = 1)) - 1
      }
    }
  }

  # set padding
  padding <- 1 + num_zeros

  # make sequence if single number
  if (length(x = seq_length) == 1) {
    seq_length <- seq_len(seq_length)
  }

  # make sequence
  new_seq <- str_pad(string = seq_length, pad = 0, width = padding, side = "left")

  return(new_seq)
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
      cli_abort(message = c("{.code custom_dir_file} {.val {custom_dir_file}} not found.",
                            "i" = "Please check file path and name provided."))
    }

    # Read file and create directory list
    output_dirs <- fread(input = custom_dir_file, data.table = FALSE) %>%
      pull() %>%
      as.list()
  }

  # Check for directories and create new ones
  lapply(output_dirs, function(dir_path){
    if (!dir.exists(dir_path)){
      dir.create(path = dir_path)
    } else {
      cli_warn(message = "The directory {.val {dir_path}} aleady exists.  No new directory created.")
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
    cli_abort(message = c("Target directory {.val {folder_file_path}} does not exist.",
                          "i" = "Please create directory or fix {.code file_path} and re-run function.")
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
    cli_abort(message = c("Target directory {.val {folder_file_path}} does not exist.",
                          "i" = "Please create directory or fix `file_path` and re-run function.")
    )
  }

  # Copy files
  system(paste0("gsutil -m cp -r ", gcp_bucket_path, " ", folder_file_path))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GENE NAME HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Update HGNC Gene Symbols
#'
#' Update human gene symbols using data from HGNC.  This function will store cached data in package directory using (BiocFileCache).  Use of this function requires internet connection on first use (or if setting `update_symbol_data = TRUE`).  Subsequent use does not require connection and will pull from cached data.
#'
#' @param input_data Data source containing gene names.  Accepted formats are:
#' \itemize{
#'  \item \code{charcter vector}
#'  \item \code{Seurat Objects}
#'  \item \code{data.frame}: genes as rownames
#'  \item \code{dgCMatrix/dgTMatrix}: genes as rownames
#'  \item \code{tibble}: genes in first column
#' }
#' @param update_symbol_data logical, whether to update cached HGNC data, default is NULL.
#' If `NULL` BiocFileCache will check and prompt for update if cache is stale.
#' If `FALSE` the BiocFileCache stale check will be skipped and current cache will be used.
#' If `TRUE` the BiocFileCache stale check will be skipped and HGNC data will be downloaded.
#' @param case_check_as_warn logical, whether case checking of features should cause abort or
#' only warn, default is FALSE (abort).  Set to TRUE if atypical names (i.e. old LOC naming) are
#' present in input_data.
#' @param verbose logical, whether to print results detailing numbers of symbols, found, updated,
#' and not found; default is TRUE.
#'
#' @return data.frame containing columns: input_features, Approved_Symbol (already approved; output unchanged), Not_Found_Symbol (symbol not in HGNC; output unchanged), Updated_Symbol (new symbol from HGNC; output updated).
#'
#' @import cli
#' @importFrom dplyr mutate filter select across left_join join_by
#' @importFrom magrittr "%>%"
#' @importFrom stats complete.cases
#' @importFrom stringr str_to_upper str_replace_na str_c str_replace
#' @importFrom tidyr drop_na everything
#'
#' @export
#'
#' @concept misc_util
#'
#' @examples
#' \dontrun{
#' new_names <- Updated_HGNC_Symbols(input_data = Seurat_Object)
#' }
#'

Updated_HGNC_Symbols <- function(
    input_data,
    update_symbol_data = NULL,
    case_check_as_warn = FALSE,
    verbose = TRUE
) {
  # Check BiocFileCache installed
  BiocFileCache_check <- is_installed(pkg = "BiocFileCache")
  if (isFALSE(x = BiocFileCache_check)) {
    cli_abort(message = c(
      "Please install the {.val BiocFileCache} package to use {.code Updated_HGNC_Symbols}",
      "i" = "This can be accomplished with the following commands: ",
      "----------------------------------------",
      "{.field `install.packages({symbol$dquote_left}BiocManager{symbol$dquote_right})`}",
      "{.field `BiocManager::install({symbol$dquote_left}BiocFileCache{symbol$dquote_right})`}",
      "----------------------------------------"
    ))
  }

  # Check input data type
  accepted_types <- c("data.frame", "dgCMatrix", "dgTMatrix")

  if (inherits(x = input_data, what = "Seurat")) {
    input_symbols <- Features(input_data)
  }
  if ((class(x = input_data) %in% accepted_types)) {
    input_symbols <- rownames(x = input_data)
  }
  if (inherits(x = input_data, what = "tibble")) {
    input_symbols <-   input_data[, 1]
  }
  if (inherits(x = input_data, what = "character")) {
    input_symbols <- input_data
  }

  # Check for duplicates
  num_duplicated <- length(x = unique(x = input_symbols[duplicated(x = input_symbols)]))

  if (num_duplicated > 0) {
    cli_abort(message = c("Input data contains duplicate gene symbols.",
                          "i" = "Check input data and/or make unique."))
  }

  # Check input symbols have correct case
  case_check <- str_to_upper(input_symbols)
  case_check <- gsub(pattern = "(.*C[0-9XY]+)ORF(.+)", replacement = "\\1orf\\2", x = case_check)
  # Currently two genes that are case anomalies so correcting them here
  case_check <- gsub(pattern = "HSA-MIR-", replacement = "hsa-mir-", x = case_check)

  if (isFALSE(x = identical(x = input_symbols, y = case_check)) && isFALSE(x = case_check_as_warn)) {
    cli_abort("Uncovered potential errors in case/capitalization of input symbols.  Please check case is correct.")
  }
  if (isFALSE(x = identical(x = input_symbols, y = case_check)) && isTRUE(x = case_check_as_warn)) {
    cli_warn(c("Uncovered potential errors in case/capitalization of input symbols.  This may cause errors in updating gene symbols.",
               "i" = "Please check case is correct and re-run if errors are found."))
  }

  # Download and process HGNC dataset if not already cached
  hgnc_data_path <- download_hgnc_data(update = update_symbol_data)

  hgnc_long_data <- readRDS(hgnc_data_path)

  input_features_df <- data.frame("input_features" = input_symbols)

  symbols_not_approved <- input_symbols[!input_symbols %in% hgnc_long_data$symbol]
  symbols_approved <- input_symbols[input_symbols %in% hgnc_long_data$symbol]

  input_features_df_approved <- input_features_df %>%
    mutate("Approved_Symbol" = ifelse(.data[["input_features"]] %in% symbols_approved, .data[["input_features"]], NA)) %>%
    drop_na()


  input_features_updated_df <- hgnc_long_data %>%
    filter(.data[["prev_symbol"]] %in% symbols_not_approved) %>%
    mutate("Updated_Symbol" = symbol) %>%
    select(any_of(c("prev_symbol", "Updated_Symbol"))) %>%
    rename("input_features" = any_of("prev_symbol")) %>%
    drop_na()

  symbols_not_found <- data.frame("input_features" = symbols_not_approved[!symbols_not_approved %in% input_features_updated_df$input_features]) %>%
    mutate("Not_Found_Symbol" = .data[["input_features"]])

  merged_df <- left_join(input_features_df, y = input_features_df_approved, by = join_by("input_features")) %>%
    left_join(symbols_not_found, by = join_by("input_features")) %>%
    left_join(input_features_updated_df, by = join_by("input_features")) %>%
    mutate(across(everything(), ~str_replace_na(string = .x, replacement = ""))) %>%
    mutate(Output_Features = str_c(.data[["Approved_Symbol"]], .data[["Not_Found_Symbol"]], .data[["Updated_Symbol"]])) %>%
    mutate(across(everything(), ~str_replace(string = .x, pattern = "^$", replacement = NA_character_))) %>%
    filter(!(.data[["input_features"]] == "QARS" & .data[["Updated_Symbol"]] == "EPRS1"))

  # Report the results
  if (isTRUE(x = verbose)) {
    num_features <- length(input_symbols)

    num_updated <- sum(complete.cases(merged_df$Updated_Symbol))
    num_not_found <- sum(complete.cases(merged_df$Not_Found_Symbol))
    num_approved <- sum(complete.cases(merged_df$Approved_Symbol))

    cli_inform(message = c("Input features contained {.field {format(x = num_features, big.mark = ',')}} gene symbols",
                           "{col_green({symbol$tick})} {.field {format(x = num_approved, big.mark = ',')}} were already approved symbols.",
                           "{col_blue({symbol$arrow_right})} {.field {format(x = num_updated, big.mark = ',')}} were updated to approved symbol.",
                           "{col_red({symbol$cross})} {.field {format(x = num_not_found, big.mark = ',')}} were not found in HGNC dataset and remain unchanged."))
  }

  # Return results
  return(merged_df)
}


#' Update MGI Gene Symbols
#'
#' Update mouse gene symbols using data from MGI  This function will store cached data in package directory using (BiocFileCache).  Use of this function requires internet connection on first use (or if setting `update_symbol_data = TRUE`).  Subsequent use does not require connection and will pull from cached data.
#'
#' @param input_data Data source containing gene names.  Accepted formats are:
#' \itemize{
#'  \item \code{charcter vector}
#'  \item \code{Seurat Objects}
#'  \item \code{data.frame}: genes as rownames
#'  \item \code{dgCMatrix/dgTMatrix}: genes as rownames
#'  \item \code{tibble}: genes in first column
#' }
#' @param update_symbol_data logical, whether to update cached MGI data, default is NULL.
#' If `NULL` BiocFileCache will check and prompt for update if cache is stale.
#' If `FALSE` the BiocFileCache stale check will be skipped and current cache will be used.
#' If `TRUE` the BiocFileCache stale check will be skipped and MGI data will be downloaded.
#' @param verbose logical, whether to print results detailing numbers of symbols, found, updated,
#' and not found; default is TRUE.
#'
#' @return data.frame containing columns: input_features, Approved_Symbol (already approved; output unchanged), Not_Found_Symbol (symbol not in MGI; output unchanged), Updated_Symbol (new symbol from MGI; output updated).
#'
#' @import cli
#' @importFrom dplyr mutate filter select across left_join join_by
#' @importFrom magrittr "%>%"
#' @importFrom stats complete.cases
#' @importFrom stringr str_to_upper str_replace_na str_c str_replace
#' @importFrom tidyr drop_na everything
#'
#' @export
#'
#' @concept misc_util
#'
#' @examples
#' \dontrun{
#' new_names <- Updated_MGI_Symbols(input_data = Seurat_Object)
#' }
#'

Updated_MGI_Symbols <- function(
    input_data,
    update_symbol_data = NULL,
    verbose = TRUE
) {
  # Check BiocFileCache installed
  BiocFileCache_check <- is_installed(pkg = "BiocFileCache")
  if (isFALSE(x = BiocFileCache_check)) {
    cli_abort(message = c(
      "Please install the {.val BiocFileCache} package to use {.code Updated_HGNC_Symbols}",
      "i" = "This can be accomplished with the following commands: ",
      "----------------------------------------",
      "{.field `install.packages({symbol$dquote_left}BiocManager{symbol$dquote_right})`}",
      "{.field `BiocManager::install({symbol$dquote_left}BiocFileCache{symbol$dquote_right})`}",
      "----------------------------------------"
    ))
  }

  # Check input data type
  accepted_types <- c("data.frame", "dgCMatrix", "dgTMatrix")

  if (inherits(x = input_data, what = "Seurat")) {
    input_symbols <- Features(input_data)
  }
  if ((class(x = input_data) %in% accepted_types)) {
    input_symbols <- rownames(x = input_data)
  }
  if (inherits(x = input_data, what = "tibble")) {
    input_symbols <-   input_data[, 1]
  }
  if (inherits(x = input_data, what = "character")) {
    input_symbols <- input_data
  }

  # Check for duplicates
  num_duplicated <- length(x = unique(x = input_symbols[duplicated(x = input_symbols)]))

  if (num_duplicated > 0) {
    cli_abort(message = c("Input data contains duplicate gene symbols.",
                          "i" = "Check input data and/or make unique."))
  }

  # Download and process HGNC dataset if not already cached
  mgi_data_path <- download_mgi_data(update = update_symbol_data)

  mgi_long_data <- readRDS(mgi_data_path)

  input_features_df <- data.frame("input_features" = input_symbols)

  symbols_not_approved <- input_symbols[!input_symbols %in% mgi_long_data$symbol]
  symbols_approved <- input_symbols[input_symbols %in% mgi_long_data$symbol]

  input_features_df_approved <- input_features_df %>%
    mutate("Approved_Symbol" = ifelse(.data[["input_features"]] %in% symbols_approved, .data[["input_features"]], NA)) %>%
    drop_na()


  input_features_updated_df <- mgi_long_data %>%
    filter(.data[["prev_symbol"]] %in% symbols_not_approved) %>%
    mutate("Updated_Symbol" = symbol) %>%
    select(any_of(c("prev_symbol", "Updated_Symbol"))) %>%
    rename("input_features" = any_of("prev_symbol")) %>%
    drop_na()

  symbols_not_found <- data.frame("input_features" = symbols_not_approved[!symbols_not_approved %in% input_features_updated_df$input_features]) %>%
    mutate("Not_Found_Symbol" = .data[["input_features"]])

  merged_df <- left_join(input_features_df, y = input_features_df_approved, by = join_by("input_features")) %>%
    left_join(symbols_not_found, by = join_by("input_features")) %>%
    left_join(input_features_updated_df, by = join_by("input_features")) %>%
    mutate(across(everything(), ~str_replace_na(string = .x, replacement = ""))) %>%
    mutate(Output_Features = str_c(.data[["Approved_Symbol"]], .data[["Not_Found_Symbol"]], .data[["Updated_Symbol"]])) %>%
    mutate(across(everything(), ~str_replace(string = .x, pattern = "^$", replacement = NA_character_))) %>%
    filter(!(.data[["input_features"]] == "QARS" & .data[["Updated_Symbol"]] == "EPRS1"))

  # Report the results
  if (isTRUE(x = verbose)) {
    num_features <- length(input_symbols)

    num_updated <- sum(complete.cases(merged_df$Updated_Symbol))
    num_not_found <- sum(complete.cases(merged_df$Not_Found_Symbol))
    num_approved <- sum(complete.cases(merged_df$Approved_Symbol))

    cli_inform(message = c("Input features contained {.field {format(x = num_features, big.mark = ',')}} gene symbols",
                           "{col_green({symbol$tick})} {.field {format(x = num_approved, big.mark = ',')}} were already approved symbols.",
                           "{col_blue({symbol$arrow_right})} {.field {format(x = num_updated, big.mark = ',')}} were updated to approved symbol.",
                           "{col_red({symbol$cross})} {.field {format(x = num_not_found, big.mark = ',')}} were not found in MGI dataset and remain unchanged."))
  }

  # Return results
  return(merged_df)
}
