#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### WRITE/CREATE ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create H5 from 10X Outputs
#'
#' Creates HDF5 formatted output analogous to the outputs created by Cell Ranger and can be read into
#' Seurat, LIGER, or SCE class object.  Requires DropletUtils package from Bioconductor.
#'
#' @param raw_data_file_path file path to raw data file(s).
#' @param source_type type of source data (Default is "10X").  Alternatively can provide "Matrix" or "data.frame".
#' @param save_file_path file path to directory to save file.
#' @param save_name name prefix for output H5 file.
#'
#' @importFrom Matrix readMM
# #' @importFrom DropletUtils write10xCounts
#' @importFrom Seurat Read10X
#'
#' @return A HDF5 format file that will be recognized as 10X Cell Ranger formatted file by Seurat or LIGER.
#'
#' @export
#'
#' @concept read_&_write
#'
#' @examples
#' \dontrun{
#' Create_10X_H5(raw_data_file_path = "file_path", save_file_path = "file_path2", save_name = "NAME")
#' }
#'

Create_10X_H5 <- function(
  raw_data_file_path,
  source_type = "10X",
  save_file_path,
  save_name
) {
  DropletUtils_check <- is_installed(pkg = "DropletUtils")
  if (isFALSE(DropletUtils_check)) {
    cli_abort(message = c(
      "Please install the {.val DropletUtils} package to use {.code Create_10X_H5}",
      "i" = "This can be accomplished with the following commands: ",
      "----------------------------------------",
      "{.field `install.packages({symbol$dquote_left}BiocManager{symbol$dquote_right})`}",
      "{.field `BiocManager::install({symbol$dquote_left}DropletUtils{symbol$dquote_right})`}",
      "----------------------------------------"
    ))
  }
  valid_source_types <- list(
    source_10X = c("10X", "10x"),
    source_matrix = c("Matrix", "matrix"),
    source_dataframe = c("Dataframe", "dataframe", "DataFrame", "data.frame", "Data.Frame", "Data.frame", "data.Frame")
    )
  if (source_type %in% valid_source_types[["source_10X"]]) {
    count_matrix <- Read10X(data.dir = raw_data_file_path)
  }
  if (source_type %in% valid_source_types[["source_matrix"]]) {
    count_matrix <- readMM(file = raw_data_file_path)
  }
  if (source_type %in% valid_source_types[["source_dataframe"]]) {
    count_matrix <- read.delim(file = raw_data_file_path,
                               header = TRUE,
                               stringsAsFactors = FALSE)
  }

  cli_inform(message = "{.field Import complete. Start write to H5}")
  # check extension
  if (isTRUE(x = check_extension(file_name = save_name, extension = ".h5"))) {
    save_name <- gsub(pattern = ".h5", replacement = "", x = save_name, fixed = TRUE)
  }

  temp_file <- tempfile(pattern = paste(save_name, "_", sep = ""),
                        tmpdir = save_file_path,
                        fileext=".h5")
  DropletUtils::write10xCounts(path = temp_file,
                 x = count_matrix,
                 barcodes = colnames(x = count_matrix),
                 gene.symbol = rownames(x = count_matrix),
                 gene.type = "Gene Expression",
                 type = "HDF5",
                 version = "3")
}


#' Create Seurat Object with Cell Bender and Raw data
#'
#' Enables easy creation of Seurat object which contains both cell bender data and raw count data as
#' separate assays within the object.
#'
#' @param raw_cell_bender_matrix matrix file containing the cell bender correct counts.
#' @param raw_counts_matrix matrix file contain the uncorrected Cell Ranger (or other) counts.
#' @param raw_assay_name a key value to use specifying the name of assay.  Default is "RAW".
#' @param min_cells `r lifecycle::badge("deprecated")` soft-deprecated. See `min.cells`.
#' @param min_features `r lifecycle::badge("deprecated")` soft-deprecated. See `min.features`.
#' @param min.cells value to supply to min.cells parameter of \code{\link[SeuratObject]{CreateSeuratObject}}.
#' Default is 5.
#' @param min.features value to supply to min.features parameter of \code{\link[SeuratObject]{CreateSeuratObject}}.
#'  Default is 200.
#' @param ... Extra parameters passed to \code{\link[SeuratObject]{CreateSeuratObject}}.
#'
#' @import Seurat
#' @importFrom dplyr intersect
#'
#' @return A Seurat Object contain both the Cell Bender corrected counts ("RNA" assay) and uncorrected
#' counts ("RAW" assay; or other name specified to `raw_assay_name`).
#'
#' @export
#'
#' @concept read_&_write
#'
#' @examples
#' \dontrun{
#' seurat_obj <- Create_CellBender_Merged_Seurat(raw_cell_bender_matrix = cb_matrix,
#' raw_counts_matrix = cr_matrix)
#' }
#'

Create_CellBender_Merged_Seurat <- function(
  raw_cell_bender_matrix = NULL,
  raw_counts_matrix = NULL,
  raw_assay_name = "RAW",
  min_cells = deprecated(),
  min_features = deprecated(),
  min.cells = 5,
  min.features = 200,
  ...
) {
  if (is_present(min_cells)) {
    deprecate_warn(when = "3.3.0",
                              what = "Create_CellBender_Merged_Seurat(min_cells)",
                              details = c("i" = "The {.code min_cells} parameter is soft-deprecated.  Please update code to use `min.cells` instead.")
    )
    min.cells <- min_cells
  }

  if (is_present(min_features)) {
    deprecate_warn(when = "3.3.0",
                              what = "Create_CellBender_Merged_Seurat(min_features)",
                              details = c("i" = "The {.code min_features} parameter is soft-deprecated.  Please update code to use `min.features` instead.")
    )
    min.features <- min_features
  }



  # Filter Cell Bender matrix for Cell Ranger cells
  cell_intersect <- intersect(x = colnames(x = raw_counts_matrix), y = colnames(x = raw_cell_bender_matrix))

  cli_inform(message = "{.field Filtering Cell Bender matrix for cells present in raw counts matrix.}")

  raw_cell_bender_matrix <- raw_cell_bender_matrix[, cell_intersect]

  # Create Seurat Object
  cli_inform(message = "{.field Creating Seurat Object from Cell Bender matrix.}")
  cell_bender_seurat <- CreateSeuratObject(counts = raw_cell_bender_matrix, min.cells = min.cells, min.features = min.features, ...)

  # Pull cell and gene names
  cell_names_seurat <- colnames(x = cell_bender_seurat)
  gene_names_seurat <- rownames(x = cell_bender_seurat)

  # Create raw counts assay object
  cli_inform(message = "{.field Creating raw counts Seurat Assay Object.}")
  counts <- CreateAssayObject(counts = raw_counts_matrix, min.cells = 0, min.features = 0)

  # Filter raw counts by created Seurat parameters
  cli_inform(message = "{.field Filtering raw counts Assay Object to match Seurat Object.}")
  counts <- subset(x = counts, cells = Cells(x = cell_bender_seurat), features = rownames(x = cell_bender_seurat))

  # Add counts assay to Seurat Object
  cli_inform(message = "{.field Adding assay to Seurat Object.}")
  cell_bender_seurat[[raw_assay_name]] <- counts

  return(cell_bender_seurat)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### READ 10X DATA ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Load in NCBI GEO data from 10X
#'
#' Enables easy loading of sparse data matrices provided by 10X genomics. That have file prefixes
#' added to them by NCBI GEO or other repos.
#'
#' @param data_dir Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv
#' files provided by 10X.
#' @param sample_list A vector of file prefixes/names if specific samples are desired.  Default is `NULL` and
#'  will load all samples in given directory.
#' @param sample_names a set of sample names to use for each sample entry in returned list.  If `NULL` will
#'  set names to the file name of each sample.
#' @param gene.column Specify which column of genes.tsv or features.tsv to use for gene names; default is 2.
#' @param cell.column Specify which column of barcodes.tsv to use for cell names; default is 1.
#' @param unique.features Make feature names unique (default TRUE).
#' @param strip.suffix Remove trailing "-1" if present in all cell barcodes.
#' @param parallel logical (default FALSE).  Whether to use multiple cores when reading in data.
#' Only possible on Linux based systems.
#' @param num_cores if `parallel = TRUE` indicates the number of cores to use for multicore processing.
#' @param merge logical (default FALSE) whether or not to merge samples into a single matrix or return
#' list of matrices.  If TRUE each sample entry in list will have cell barcode prefix added.  The prefix
#' will be taken from `sample_names`.
#'
#' @return If features.csv indicates the data has multiple data types, a list
#'   containing a sparse matrix of the data from each type will be returned.
#'   Otherwise a sparse matrix containing the expression data will be returned.
#'
#' @references Code used in function has been slightly modified from `Seurat::Read10X` function of
#' Seurat package \url{https://github.com/satijalab/seurat} (License: GPL-3).  Function was modified to
#' support file prefixes and altered loop by Samuel Marsh for scCustomize (also previously posted as
#' potential PR to Seurat GitHub).
#'
#' @import parallel
#' @import pbapply
#' @importFrom Matrix readMM
#' @importFrom utils read.delim txtProgressBar setTxtProgressBar read.table
#'
#' @export
#'
#' @concept read_&_write
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/data/directory'
#' expression_matrices <- Read10X_GEO(data_dir = data_dir)
#' # To create object from single file
#' seurat_object = CreateSeuratObject(counts = expression_matrices[[1]])
#' }
#'

Read10X_GEO <- function(
  data_dir = NULL,
  sample_list = NULL,
  sample_names = NULL,
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE,
  parallel = FALSE,
  num_cores = NULL,
  merge = FALSE
) {
  if (!dir.exists(paths = data_dir)) {
    cli_abort(message = "Directory provided does not exist")
  }
  if (length(x = data_dir) > 1) {
    cli_abort(message = "{.code Read10X_GEO} only supports reading from single data directory at a time.")
  }

  # Confirm num_cores specified
  if (isTRUE(x = parallel) && is.null(x = num_cores)) {
    cli_abort("If {.code parallel = TRUE} then {.code num_cores} must be specified.")
  }

  if (is.null(x = sample_list)) {
    # pull all file names from directory
    file.list <- list.files(path = data_dir, pattern = "barcodes.tsv", full.names = FALSE)
    # Remove "barcodes.tsv.gz" file suffix
    sample_list <- gsub(pattern = "barcodes.tsv.gz", x = file.list, replacement = "")
  }

  if (!is.null(x = sample_names)) {
    if (length(x = sample_names) != length(x = sample_list)) {
      cli_abort(message = "Length of {.code sample_names} provided {.field {length(x = sample_names)}} does not equal length of {.code sample_list} {.field {length(x = sample_list)}}.")
    }
  }


  cli_inform(message = "{.field Reading 10X files from directory}")
  pboptions(char = "=")
  if (isTRUE(x = parallel)) {
    cli_inform(message = c("NOTE: Progress bars not currently supported for parallel processing.",
                           "NOTE: Parallel processing will not report informative error messages.", "
                           If function fails set {.code parallel = FALSE} and re-run for informative error reporting.\n"))
    raw_data_list <- mclapply(mc.cores = num_cores, 1:length(sample_list), function(i) {
      barcode.loc <- file.path(data_dir, paste0(sample_list[i], 'barcodes.tsv.gz'))
      gene.loc <- file.path(data_dir, paste0(sample_list[i], 'genes.tsv.gz'))
      features.loc <- file.path(data_dir, paste0(sample_list[i], 'features.tsv.gz'))
      matrix.loc <- file.path(data_dir, paste0(sample_list[i], 'matrix.mtx.gz'))
      # Flag to indicate if this data is from CellRanger >= 3.0
      pre_ver_3 <- file.exists(gene.loc)
      if (!file.exists(barcode.loc)) {
        cli_abort(message = "Barcode file missing. Expecting {val {basename(path = barcode.loc)}}")
      }
      if (isFALSE(x = pre_ver_3) && !file.exists(features.loc) ) {
        cli_abort(message = "Gene name or features file missing. Expecting {val {basename(path = features.loc)}}")
      }
      if (!file.exists(matrix.loc)) {
        cli_abort(message = "Expression matrix file missing. Expecting {val {basename(path = matrix.loc)}}")
      }
      data <- readMM(file = matrix.loc)
      cell.barcodes <- read.table(file = barcode.loc, header = FALSE, sep = '\t', row.names = NULL)
      if (ncol(x = cell.barcodes) > 1) {
        cell.names <- cell.barcodes[, cell.column]
      } else {
        cell.names <- readLines(con = barcode.loc)
      }
      if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
        cell.names <- as.vector(x = as.character(x = sapply(
          X = cell.names,
          FUN = ExtractField,
          field = 1,
          delim = "-"
        )))
      }
      if (is.null(x = names(x = data_dir))) {
        if (i < 2) {
          colnames(x = data) <- cell.names
        } else {
          colnames(x = data) <- paste0(i, "_", cell.names)
        }
      } else {
        colnames(x = data) <- paste0(names(x = data_dir)[i], "_", cell.names)
      }
      feature.names <- read.delim(
        file = ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc),
        header = FALSE,
        stringsAsFactors = FALSE
      )
      if (any(is.na(x = feature.names[, gene.column]))) {
        warning(
          'Some features names are NA. Replacing NA names with ID from the opposite column requested',
          call. = FALSE,
          immediate. = TRUE
        )
        na.features <- which(x = is.na(x = feature.names[, gene.column]))
        replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
        feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
      }
      if (isTRUE(x = unique.features)) {
        fcols = ncol(x = feature.names)
        if (fcols < gene.column) {
          cli_abort(message = c("{.code gene.column} was set to {.val {gene.column}}, but feature.tsv.gz (or genes.tsv) only has {.field {cols}} columns.",
                                "i" = "Try setting the {.code gene.column} argument to a value <= to {.field {cols}}."))
        }
        rownames(x = data) <- make.unique(names = feature.names[, gene.column])
      }
      # In cell ranger 3.0, a third column specifying the type of data was added
      # and we will return each type of data as a separate matrix
      if (ncol(x = feature.names) > 2) {
        data_types <- factor(x = feature.names$V3)
        lvls <- levels(x = data_types)
        # if (length(x = lvls) > 1 && length(x = full.data) == 0) {
        #   message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
        # }
        expr_name <- "Gene Expression"
        if (expr_name %in% lvls) { # Return Gene Expression first
          lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
        }
        data <- lapply(
          X = lvls,
          FUN = function(l) {
            return(data[data_types == l, , drop = FALSE])
          }
        )
        names(x = data) <- lvls
      } else{
        data <- list(data)
      }

      if (length(x = data) == 1) {
        return(data[[1]])
      } else {
        return(data)
      }
    })
  } else {
    raw_data_list <- pblapply(1:length(sample_list), function(i) {
      barcode.loc <- file.path(data_dir, paste0(sample_list[i], 'barcodes.tsv.gz'))
      gene.loc <- file.path(data_dir, paste0(sample_list[i], 'genes.tsv.gz'))
      features.loc <- file.path(data_dir, paste0(sample_list[i], 'features.tsv.gz'))
      matrix.loc <- file.path(data_dir, paste0(sample_list[i], 'matrix.mtx.gz'))
      # Flag to indicate if this data is from CellRanger >= 3.0
      pre_ver_3 <- file.exists(gene.loc)
      if (!file.exists(barcode.loc)) {
        cli_abort(message = "Barcode file missing. Expecting {.val {basename(path = barcode.loc)}}")
      }
      if (isFALSE(x = pre_ver_3) && !file.exists(features.loc) ) {
        cli_abort(message = "Gene name or features file missing. Expecting {.val {basename(path = features.loc)}}")
      }
      if (!file.exists(matrix.loc)) {
        cli_abort(message = "Expression matrix file missing. Expecting {.val {basename(path = matrix.loc)}}")
      }
      data <- readMM(file = matrix.loc)
      cell.barcodes <- read.table(file = barcode.loc, header = FALSE, sep = '\t', row.names = NULL)
      if (ncol(x = cell.barcodes) > 1) {
        cell.names <- cell.barcodes[, cell.column]
      } else {
        cell.names <- readLines(con = barcode.loc)
      }
      if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
        cell.names <- as.vector(x = as.character(x = sapply(
          X = cell.names,
          FUN = ExtractField,
          field = 1,
          delim = "-"
        )))
      }
      if (is.null(x = names(x = data_dir))) {
        if (i < 2) {
          colnames(x = data) <- cell.names
        } else {
          colnames(x = data) <- paste0(i, "_", cell.names)
        }
      } else {
        colnames(x = data) <- paste0(names(x = data_dir)[i], "_", cell.names)
      }
      feature.names <- read.delim(
        file = ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc),
        header = FALSE,
        stringsAsFactors = FALSE
      )
      if (any(is.na(x = feature.names[, gene.column]))) {
        warning(
          'Some features names are NA. Replacing NA names with ID from the opposite column requested',
          call. = FALSE,
          immediate. = TRUE
        )
        na.features <- which(x = is.na(x = feature.names[, gene.column]))
        replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
        feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
      }
      if (isTRUE(x = unique.features)) {
        fcols = ncol(x = feature.names)
        if (fcols < gene.column) {
          cli_abort(message = c("{.code gene.column} was set to {.val {gene.column}}, but feature.tsv.gz (or genes.tsv) only has {.field {cols}} columns.",
                                "i" = "Try setting the {.code gene.column} argument to a value <= to {.field {cols}}."))
        }
        rownames(x = data) <- make.unique(names = feature.names[, gene.column])
      }
      # In cell ranger 3.0, a third column specifying the type of data was added
      # and we will return each type of data as a separate matrix
      if (ncol(x = feature.names) > 2) {
        data_types <- factor(x = feature.names$V3)
        lvls <- levels(x = data_types)
        # if (length(x = lvls) > 1 && length(x = full.data) == 0) {
        #   message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
        # }
        expr_name <- "Gene Expression"
        if (expr_name %in% lvls) { # Return Gene Expression first
          lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
        }
        data <- lapply(
          X = lvls,
          FUN = function(l) {
            return(data[data_types == l, , drop = FALSE])
          }
        )
        names(x = data) <- lvls
      } else{
        data <- list(data)
      }

      if (length(x = data) == 1) {
        return(data[[1]])
      } else {
        return(data)
      }
    })
  }
  # # Combine all the data from different directories into one big matrix, note this
  # # assumes that all data directories essentially have the same features files
  # if (single.matrix) {
  #   message("Creating combined matrix output")
  #   list_of_data <- list()
  #   pb <- txtProgressBar(min = 0, max = length(x = full.data[[1]]), style = 3, file = stderr())
  #   for (j in 1:length(x = full.data[[1]])) {
  #     list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, FUN = `[[`, j))
  #     # Fix for Issue #913
  #     list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "dgCMatrix")
  #     setTxtProgressBar(pb = pb, value = j)
  #   }
  #   close(con = pb)
  #   names(x = list_of_data) <- names(x = full.data[[1]])
  #   # If multiple features, will return a list, otherwise
  #   # a matrix.
  #   if (length(x = list_of_data) == 1) {
  #     return(list_of_data[[1]])
  #   } else {
  #     return(list_of_data)
  #   }
  # }

  # Name the list
  if (!is.null(x = sample_names)) {
    names(x = raw_data_list) <- sample_names
  } else {
    names(x = raw_data_list) <- sample_list
  }

  # Merge data
  if (merge) {
    raw_data_merged <- Merge_Sparse_Data_All(matrix_list = raw_data_list, add_cell_ids = names(x = raw_data_list))
    return(raw_data_merged)
  }

  # return list
  return(raw_data_list)
}


#' Load in NCBI GEO data from 10X in HDF5 file format
#'
#' Enables easy loading of HDF5 data matrices provided by 10X genomics. That have file prefixes added to
#' them by NCBI GEO or other repos or programs (i.e. Cell Bender)
#'
#' @param data_dir Directory containing the .h5 files provided by 10X.
#' @param sample_list A vector of file prefixes/names if specific samples are desired.  Default is `NULL` and
#' will load all samples in given directory.
#' @param sample_names a set of sample names to use for each sample entry in returned list.  If `NULL`
#' will set names to the file name of each sample.
#' @param shared_suffix a suffix and file extension shared by all samples.
#' @param parallel logical (default FALSE).  Whether to use multiple cores when reading in data.
#' Only possible on Linux based systems.
#' @param num_cores if `parallel = TRUE` indicates the number of cores to use for multicore processing.
#' @param merge logical (default FALSE) whether or not to merge samples into a single matrix or return
#' list of matrices.  If TRUE each sample entry in list will have cell barcode prefix added.  The prefix
#' will be taken from `sample_names`.
#' @param ... Additional arguments passed to \code{\link[Seurat]{Read10X_h5}}
#'
#' @return If the data has multiple data types, a list
#'   containing a sparse matrix of the data from each type will be returned.
#'   Otherwise a sparse matrix containing the expression data will be returned.
#'
#' @import parallel
#' @import pbapply
#' @importFrom Matrix readMM
#' @importFrom utils read.delim txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @concept read_&_write
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/data/directory'
#' expression_matrices <- Read10X_h5_GEO(data_dir = data_dir)
#' # To create object from single file
#' seurat_object = CreateSeuratObject(counts = expression_matrices[[1]])
#' }
#'

Read10X_h5_GEO <- function(
  data_dir = NULL,
  sample_list = NULL,
  sample_names = NULL,
  shared_suffix = NULL,
  parallel = FALSE,
  num_cores = NULL,
  merge = FALSE,
  ...
) {
  if (!dir.exists(paths = data_dir)) {
    cli_abort(message = "Directory provided does not exist")
  }
  if (length(x = data_dir) > 1) {
    cli_abort(message = "{.code Read10X_h5_GEO} only supports reading from single data directory at a time.")
  }

  # Confirm num_cores specified
  if (isTRUE(x = parallel) && is.null(x = num_cores)) {
    cli_abort("If {.code parallel = TRUE} then {.code num_cores} must be specified.")
  }

  file.list <- list.files(path = data_dir, pattern = ".h5", full.names = FALSE)

  # Remove file suffix if provided
  if (!is.null(x = shared_suffix)) {
    shared_suffix <- gsub(pattern = ".h5", replacement = "", x = shared_suffix)
  }

  if (is.null(x = sample_list)) {
    if (is.null(x = shared_suffix)) {
      sample_list <- gsub(pattern = ".h5", x = file.list, replacement = "")
    } else {
      sample_list <- gsub(pattern = paste0(shared_suffix, ".h5"), x = file.list, replacement = "")
    }
  }

  # Check sample_names length is ok
  if (!is.null(x = sample_names) && length(x = sample_names) != length(x = sample_list)) {
    cli_abort(message = "Length of {.code sample_names} {.field {length(x = sample_names)}} must be equal to number of samples {.field {length(x = sample_list)}}.")
  }

  cli_inform(message = "{.field Reading 10X H5 files from directory}")
  pboptions(char = "=")
  if (isTRUE(x = parallel)) {
    cli_inform(message = c("NOTE: Progress bars not currently supported for parallel processing.",
                           "NOTE: Parallel processing will not report informative error messages.", "
                           If function fails set {.code parallel = FALSE} and re-run for informative error reporting.\n"))
    raw_data_list <- mclapply(mc.cores = num_cores, 1:length(sample_list), function(i) {
      h5_loc <- file.path(data_dir, paste0(sample_list[i], shared_suffix, ".h5"))
      data <- Read10X_h5(filename = h5_loc, ...)
    })
  } else {
    raw_data_list <- pblapply(1:length(x = sample_list), function(i) {
      h5_loc <- file.path(data_dir, paste0(sample_list[i], shared_suffix, ".h5"))
      data <- Read10X_h5(filename = h5_loc, ...)
    })
  }

  # Name the matrices
  if (is.null(x = sample_names)) {
    names(x = raw_data_list) <- sample_list
  } else {
    names(x = raw_data_list) <- sample_names
  }

  # Merge data
  if (isTRUE(x = merge)) {
    raw_data_merged <- Merge_Sparse_Data_All(matrix_list = raw_data_list, add_cell_ids = names(x = raw_data_list))
    return(raw_data_merged)
  }

  # return object
  return(raw_data_list)
}


#' Load 10X count matrices from multiple directories
#'
#' Enables easy loading of sparse data matrices provided by 10X genomics that are present in multiple
#' subdirectories.  Can function with either default output directory structure of Cell Ranger or
#' custom directory structure.
#'
#' @param base_path path to the parent directory which contains all of the subdirectories of interest.
#' @param secondary_path path from the parent directory to count matrix files for each sample.
#' @param default_10X_path logical (default TRUE) sets the secondary path variable to the default 10X
#' directory structure.
#' @param cellranger_multi logical, whether samples were processed with Cell Ranger `multi`, default is FALSE.
#' @param sample_list a vector of sample directory names if only specific samples are desired.  If `NULL` will
#' read in subdirectories in parent directory.
#' @param sample_names a set of sample names to use for each sample entry in returned list.  If `NULL` will
#' set names to the subdirectory name of each sample.
#' @param parallel logical (default FALSE) whether or not to use multi core processing to read in matrices.
#' @param num_cores how many cores to use for parallel processing.
#' @param merge logical (default FALSE) whether or not to merge samples into a single matrix or return
#' list of matrices.  If TRUE each sample entry in list will have cell barcode prefix added.  The prefix
#' will be taken from `sample_names`.
#' @param ... Extra parameters passed to \code{\link[Seurat]{Read10X}}.
#'
#' @return a list of sparse matrices (merge = FALSE) or a single sparse matrix (merge = TRUE).
#'
#' @import parallel
#' @import pbapply
#' @importFrom Seurat Read10X
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @concept read_&_write
#'
#' @examples
#' \dontrun{
#' base_path <- 'path/to/data/directory'
#' expression_matrices <- Read10X_Multi_Directory(base_path = base_path)
#' }
#'

Read10X_Multi_Directory <- function(
  base_path,
  secondary_path = NULL,
  default_10X_path = TRUE,
  cellranger_multi = FALSE,
  sample_list = NULL,
  sample_names = NULL,
  parallel = FALSE,
  num_cores = NULL,
  merge = FALSE,
  ...
) {
  # Confirm num_cores specified
  if (parallel && is.null(x = num_cores)) {
    cli_abort("If {.code parallel = TRUE} then {.code num_cores} must be specified.")
  }
  # Confirm directory exists
  if (dir.exists(paths = base_path) == FALSE) {
    cli_abort(message = "Directory: {.val {base_path}} specified by {.code base_path} does not exist.")
  }
  # Detect libraries if sample_list is NULL
  if (is.null(x = sample_list)) {
    sample_list <- Pull_Directory_List(base_path = base_path)
  }
  # Add file path for 10X default directories
  if (isTRUE(x = default_10X_path) && !is.null(x = secondary_path)) {
    cli_abort(message = "If {.code default_10X_path = TRUE} then {.code secondary_path} must be NULL.")
  }

  if (isFALSE(x = default_10X_path) && !is.null(x = secondary_path) && isTRUE(x = cellranger_multi)) {
    cli_abort(message = "If {.code cellranger_multi = TRUE} then {.code default_10X_path} must be TRUE")
  }

  if (isTRUE(x = default_10X_path)) {
    if (isTRUE(x = cellranger_multi)) {
      secondary_path <- "/outs/per_sample_outs/"
      multi_extra_path <- "count/sample_filtered_feature_bc_matrix"
    } else {
      secondary_path <- "outs/filtered_feature_bc_matrix/"
    }
  }
  if (is.null(x = secondary_path)) {
    secondary_path <- ""
  }
  # Check if full directory path exists
  for (i in 1:length(x = sample_list)) {
    full_directory_path <- file.path(base_path, sample_list[i], secondary_path)
    if (dir.exists(paths = full_directory_path) == FALSE) {
      cli_abort(message = "Full Directory does not exist {.val {full_directory_path}} was not found.")
    }
  }
  # read data
  cli_inform(message = "{.field Reading gene expression files.}")
  if (isTRUE(x = parallel)) {
    cli_inform(message = c("NOTE: Progress bars not currently supported for parallel processing.",
                           "NOTE: Parallel processing will not report informative error messages.", "
                           If function fails set {.code parallel = FALSE} and re-run for informative error reporting.\n"))
    # *** Here is where the swap of mclapply or pbmclapply is occuring ***
    raw_data_list <- mclapply(mc.cores = num_cores, 1:length(x = sample_list), function(x) {
      if (isTRUE(x = cellranger_multi)) {
        file_path <- file.path(base_path, sample_list[x], secondary_path, sample_list[x], multi_extra_path)
      } else {
        file_path <- file.path(base_path, sample_list[x], secondary_path)
      }
      raw_data <- Read10X(data.dir = file_path, ...)
      return(raw_data)
    })
  } else {
    raw_data_list <- pblapply(1:length(x = sample_list), function(x) {
      if (is.null(x = secondary_path)) {
        file_path <- file.path(base_path, sample_list[x])
      } else {
        if (isTRUE(x = cellranger_multi)) {
          file_path <- file.path(base_path, sample_list[x], secondary_path, sample_list[x], multi_extra_path)
        } else {
          file_path <- file.path(base_path, sample_list[x], secondary_path)
        }
      }
      raw_data <- Read10X(data.dir = file_path, ...)
    })
  }
  # Name the list items
  if (is.null(x = sample_names)) {
    names(x = raw_data_list) <- sample_list
  } else {
    names(x = raw_data_list) <- sample_names
  }
  # Merge data
  if (isTRUE(x = merge)) {
    raw_data_merged <- Merge_Sparse_Data_All(matrix_list = raw_data_list, add_cell_ids = names(x = raw_data_list))
    return(raw_data_merged)
  }
  return(raw_data_list)
}


#' Load 10X h5 count matrices from multiple directories
#'
#' Enables easy loading of sparse data matrices provided by 10X genomics that are present in multiple
#' subdirectories.  Can function with either default output directory structure of Cell Ranger or
#' custom directory structure.
#'
#' @param base_path path to the parent directory which contains all of the subdirectories of interest.
#' @param secondary_path path from the parent directory to count matrix files for each sample.
#' @param default_10X_path logical (default TRUE) sets the secondary path variable to the default 10X
#' directory structure.
#' @param cellranger_multi logical, whether samples were processed with Cell Ranger `multi`, default is FALSE.
#' @param h5_filename name of h5 file (including .h5 suffix).  If all h5 files have same name (i.e. Cell Ranger output)
#' then use full file name.  By default function uses Cell Ranger name: "filtered_feature_bc_matrix.h5".
#' If h5 files have sample specific prefixes (i.e. from Cell Bender) then use only the shared part of file
#' name (e.g., "_filtered_out.h5").
#' @param sample_list a vector of sample directory names if only specific samples are desired.  If `NULL` will
#' read in subdirectories in parent directory.
#' @param sample_names a set of sample names to use for each sample entry in returned list.  If `NULL` will
#' set names to the subdirectory name of each sample.
#' @param replace_suffix logical (default FALSE).  Whether or not to replace the barcode suffixes of matrices
#' using \code{\link{Replace_Suffix}}.
#' @param new_suffix_list a vector of new suffixes to replace existing suffixes if `replace_suffix = TRUE`.
#' See \code{\link{Replace_Suffix}} for more information.  To remove all suffixes set `new_suffix_list = ""`.
#' @param parallel logical (default FALSE) whether or not to use multi core processing to read in matrices.
#' @param num_cores how many cores to use for parallel processing.
#' @param merge logical (default FALSE) whether or not to merge samples into a single matrix or return
#' list of matrices.  If TRUE each sample entry in list will have cell barcode prefix added.  The prefix
#' will be taken from `sample_names`.
#' @param ... Extra parameters passed to \code{\link[Seurat]{Read10X_h5}}.
#'
#' @return a list of sparse matrices (merge = FALSE) or a single sparse matrix (merge = TRUE).
#'
#' @import parallel
#' @import pbapply
#' @importFrom Seurat Read10X_h5
#' @importFrom stringr str_extract
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @concept read_&_write
#'
#' @examples
#' \dontrun{
#' base_path <- 'path/to/data/directory'
#' expression_matrices <- Read10X_h5_Multi_Directory(base_path = base_path)
#' }
#'

Read10X_h5_Multi_Directory <- function(
  base_path,
  secondary_path = NULL,
  default_10X_path = TRUE,
  cellranger_multi = FALSE,
  h5_filename = "filtered_feature_bc_matrix.h5",
  sample_list = NULL,
  sample_names = NULL,
  replace_suffix = FALSE,
  new_suffix_list = NULL,
  parallel = FALSE,
  num_cores = NULL,
  merge = FALSE,
  ...
) {
  # Confirm num_cores specified
  if (isTRUE(x = parallel) && is.null(x = num_cores)) {
    cli_abort("If {.code parallel = TRUE} then {.code num_cores} must be specified.")
  }
  # Confirm directory exists
  if (dir.exists(paths = base_path) == FALSE) {
    cli_abort(message = "Directory: {.val {base_path}} specified by {.code base_path} does not exist.")
  }
  # Detect libraries if sample_list is NULL
  if (is.null(x = sample_list)) {
    sample_list <- Pull_Directory_List(base_path = base_path)
  }

  # Add file path for 10X default directories
  if (isTRUE(x = default_10X_path) && !is.null(x = secondary_path)) {
    cli_abort(message = "If {.code default_10X_path = TRUE} then {.code secondary_path} must be NULL.")
  }

  if (isFALSE(x = default_10X_path) && !is.null(x = secondary_path) && isTRUE(x = cellranger_multi)) {
    cli_abort(message = "If {.code cellranger_multi = TRUE} then {.code default_10X_path} must be TRUE")
  }

  if (isTRUE(x = default_10X_path)) {
    if (isTRUE(x = cellranger_multi)) {
      secondary_path <- "/outs/per_sample_outs/"
      multi_extra_path <- "count/"
    } else {
      secondary_path <- "outs/"
    }
  }

  if (is.null(x = secondary_path)) {
    secondary_path <- ""
  }
  # Check if full directory path exists
  for (i in 1:length(x = sample_list)) {
    full_directory_path <- file.path(base_path, sample_list[i], secondary_path)
    if (dir.exists(paths = full_directory_path) == FALSE) {
      cli_abort(message = "Full Directory does not exist {.val {full_directory_path}} was not found.")
    }
  }

  # read data
  cli_inform(message = "{.field Reading gene expression files.}")
  if (isTRUE(x = parallel)) {
    cli_inform(message = c("NOTE: Progress bars not currently supported for parallel processing.",
                           "NOTE: Parallel processing will not report informative error messages.", "
                           If function fails set {.code parallel = FALSE} and re-run for informative error reporting.\n"))
    # *** Here is where the swap of mclapply or pbmclapply is occuring ***
    raw_data_list <- mclapply(mc.cores = num_cores, 1:length(x = sample_list), function(x) {
      if (isTRUE(x = cellranger_multi)) {
        file_path <- file.path(base_path, sample_list[x], secondary_path, sample_list[x], multi_extra_path, h5_filename)
      } else {
        file_path <- file.path(base_path, sample_list[x], secondary_path, h5_filename)
      }
      raw_data <- Read10X_h5(filename = file_path, ...)
      return(raw_data)
    })
  } else {
    raw_data_list <- pblapply(1:length(x = sample_list), function(x) {
      if (isTRUE(x = cellranger_multi)) {
        file_path <- file.path(base_path, sample_list[x], secondary_path, sample_list[x], multi_extra_path, h5_filename)
      } else {
        file_path <- file.path(base_path, sample_list[x], secondary_path, h5_filename)
      }
      raw_data <- Read10X_h5(filename = file_path, ...)
    })
  }
  # Name the list items
  if (is.null(x = sample_names)) {
    names(x = raw_data_list) <- sample_list
  } else {
    names(x = raw_data_list) <- sample_names
  }

  # Replace Suffixes
  if (isTRUE(x = replace_suffix)) {
    if (is.null(x = new_suffix_list)) {
      cli_abort(message = "No values provided to {.code new_suffix_list} but {.code replace_suffix = TRUE}.")
    }

    current_suffix_list <- sapply(1:length(x = raw_data_list), function(x) {
      unique(str_extract(string = colnames(x = raw_data_list[[x]]), pattern = "-.$"))
    })

    if (length(x = new_suffix_list) != 1 & length(x = new_suffix_list) != length(x = current_suffix_list)) {
      cli_abort(message = c("`new_suffix_list` must be either single value or list of values equal to the number of samples.",
                            "i" = "Number of samples is: {.field {length(current_suffix_list)}} and number of new_suffixes provided is: {.field {length(x = new_suffix_list)}}."))
    }

    raw_data_list <- Replace_Suffix(data = raw_data_list, current_suffix = current_suffix_list, new_suffix = new_suffix_list)
  }

  # Merge data
  if (isTRUE(x = merge)) {
    raw_data_merged <- Merge_Sparse_Data_All(matrix_list = raw_data_list, add_cell_ids = names(x = raw_data_list))
    return(raw_data_merged)
  }
  return(raw_data_list)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### READ DELIM DATA ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Load in NCBI GEO data formatted as single file per sample
#'
#' Can read delimited file types (i.e. csv, tsv, txt)
#'
#' @param data_dir Directory containing the files.
#' @param file_suffix The file suffix of the individual files.  Must be the same across all files being
#' imported.  This is used to detect files to import and their GEO IDs.
#' @param move_genes_rownames logical.  Whether gene IDs are present in first column or in row names of
#' delimited file.  If TRUE will move the first column to row names before creating final matrix.
#' Default is TRUE.
#' @param sample_list a vector of samples within directory to read in (can be either with or
#' without `file_suffix` see `full_names`).  If NULL will read in all subdirectories.
#' @param full_names logical (default FALSE).  Whether or not the `sample_list` vector includes the file suffix.
#' If `FALSE` the function will add suffix based on `file_suffix` parameter.
#' @param sample_names a set of sample names to use for each sample entry in returned list.
#' If `NULL` will set names to the directory name of each sample.
#' @param barcode_suffix_period Is the barcode suffix a period and should it be changed to "-".  Default (FALSE;
#' barcodes will be left identical to their format in input files.).  If TRUE "." in barcode suffix will
#' be changed to "-".
#' @param parallel logical (default FALSE).  Whether to use multiple cores when reading in data.
#' Only possible on Linux based systems.
#' @param num_cores if `parallel = TRUE` indicates the number of cores to use for multicore processing.
#' @param merge logical (default FALSE) whether or not to merge samples into a single matrix or return
#' list of matrices.  If TRUE each sample entry in list will have cell barcode prefix added.  The prefix
#' will be taken from `sample_names`.
#'
#' @return List of gene x cell matrices in list format named by sample name.
#'
#' @import Matrix
#' @import parallel
#' @import pbapply
#' @importFrom data.table fread
#' @importFrom magrittr "%>%"
#' @importFrom methods as
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom utils read.delim
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @concept read_&_write
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/data/directory'
#' expression_matrices <- Read_GEO_Delim(data_dir = data_dir)
#' }
#'

Read_GEO_Delim <- function(
  data_dir,
  file_suffix,
  move_genes_rownames = TRUE,
  sample_list = NULL,
  full_names = FALSE,
  sample_names = NULL,
  barcode_suffix_period = FALSE,
  parallel = FALSE,
  num_cores = NULL,
  merge = FALSE
) {
  # Create list of all files in directory
  possible_file_list <- list.files(path = data_dir, pattern = file_suffix, full.names = FALSE)

  # Check files found
  if (is.null(x = possible_file_list)) {
    cli_abort(message = "No files found.  Check that {.code data_dir} and {.code file_suffix} are correct.")
  }

  # Set all files to be used if sample_list is NULL
  if (is.null(sample_list)) {
    file_list <- possible_file_list
  }

  # Confirm num_cores specified
  if (isTRUE(x = parallel) && is.null(x = num_cores)) {
    cli_abort("If {.code parallel = TRUE} then {.code num_cores} must be specified.")
  }

  # Read in subset of files
  if (!is.null(x = sample_list)) {
    # Add suffix
    if (isTRUE(x = full_names)) {
      file_list <- sample_list
    } else {
      file_list <- paste0(sample_list, file_suffix)
    }
    file_list <- file_list

    if (any(!file_list %in% possible_file_list)) {
      bad_file_list <- file_list[!file_list %in% possible_file_list]
      file_list <- file_list[file_list %in% possible_file_list]
      if (length(x = file_list) == 0) {
        cli_abort(message = c("No requested files found.",
                              "i" = "Check that {.code data_dir} and {.code file_suffix} are correct and {.code full_names} parameter is accurate."))
      }
      cli_warn(message = c("The following files were not imported as they were not found in specified directory:",
                           "i" = "{.field {glue_collapse_scCustom(input_string = bad_file_list, and = TRUE)}}"))
    }
  }

  # Get sample names
  if (is.null(x = sample_names)) {
    sample_names <- gsub(pattern = file_suffix, x = file_list, replacement = "")
  } else {
    sample_names <- sample_names
  }

  # Read in files
  cli_inform(message = "{.field Reading gene expression files from directory}")
  pboptions(char = "=")
  if (isTRUE(x = parallel)) {
    cli_inform(message = c("NOTE: Progress bars not currently supported for parallel processing.",
                           "NOTE: Parallel processing will not report informative error messages.", "
                           If function fails set {.code parallel = FALSE} and re-run for informative error reporting.\n"))
    raw_data_list <- mclapply(mc.cores = num_cores, 1:length(x = file_list), function(i) {
      dge_loc <- file.path(data_dir, file_list[i])
      data <- fread(file = dge_loc, data.table = F)
      if (isTRUE(x = move_genes_rownames)) {
        first_col_name <- colnames(x = data[1])
        data <- data %>%
          column_to_rownames(first_col_name)
      }
      if (isTRUE(x = barcode_suffix_period)) {
        colnames(x = data) <- gsub("\\.", "-", colnames(x = data))
      }
      data_sparse <- as(data, "Matrix")
      return(data_sparse)
    })
  } else {
    raw_data_list <- pblapply(1:length(x = file_list), function(i) {
      dge_loc <- file.path(data_dir, file_list[i])
      data <- fread(file = dge_loc, data.table = F)
      if (isTRUE(x = move_genes_rownames)) {
        first_col_name <- colnames(x = data[1])
        data <- data %>%
          column_to_rownames(first_col_name)
      }
      # Check all columns numeric
      col_data_numeric <- sapply(data, is.numeric)
      if (!all(col_data_numeric)) {
        cli_abort(message = c("One or more columns in the file: {.val {dge_loc}} contains non-numeric data.",
                              "i" = "Please check original file and/or that parameter {.code move_genes_rownames} is set appropriately."))
      }
      if (isTRUE(x = barcode_suffix_period)) {
        colnames(x = data) <- gsub("\\.", "-", colnames(x = data))
      }
      data_sparse <- as(data, "Matrix")
      return(data_sparse)
    })
  }

  # Name the items in list
  names(x = raw_data_list) <- sample_names

  # Check matrices
  for (i in 1:length(x = raw_data_list)) {
    CheckMatrix_scCustom(object = raw_data_list[[i]])
  }

  # Merge data
  if (isTRUE(x = merge)) {
    raw_data_merged <- Merge_Sparse_Data_All(matrix_list = raw_data_list, add_cell_ids = names(x = raw_data_list))
    return(raw_data_merged)
  }

  # return list
  return(raw_data_list)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### READ CellBender DATA ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Load CellBender h5 matrices (corrected)
#'
#' Extract sparse matrix with corrected counts from CellBender h5 output file.
#'
#' @param file_name Path to h5 file.
#' @param use.names Label row names with feature names rather than ID numbers (default TRUE).
#' @param unique.features Make feature names unique (default TRUE).
#' @param h5_group_name Name of the group within H5 file that contains count data.  This is only
#' required if H5 file contains multiple subgroups and non-default names.  Default is `NULL`.
#' @param feature_slot_name Name of the slot contain feature names/ids.  Must be one of:
#' "features"(Cell Ranger v3+) or "genes" (Cell Ranger v1/v2 or STARsolo).  Default is "features".
#'
#' @return sparse matrix
#'
#' @references Code used in function has been modified from `Seurat::Read10X_h5` function of
#' Seurat package \url{https://github.com/satijalab/seurat} (License: GPL-3).
#'
#' @import Matrix
#'
#' @export
#'
#' @concept read_&_write
#'
#' @examples
#' \dontrun{
#' mat <- Read_CellBender_h5_Mat(file_name = "/SampleA_out_filtered.h5")
#' }
#'

Read_CellBender_h5_Mat <- function(
    file_name,
    use.names = TRUE,
    unique.features = TRUE,
    h5_group_name = NULL,
    feature_slot_name = "features"
) {
  # Check hdf5r installed
  hdf5r_check <- is_installed(pkg = "hdf5r")
  if (isFALSE(x = hdf5r_check)) {
    cli_abort(message = c(
      "Please install the {.val hdf5r} package to use {.code Read_CellBender_h5_Mat} and read HDF5 files.",
      "i" = "This can be accomplished with the following commands: ",
      "----------------------------------------",
      "{.field `install.packages({symbol$dquote_left}hdf5r{symbol$dquote_right})`}",
      "----------------------------------------"
    ))
  }

  # Check file
  if (!file.exists(file_name)) {
    cli_abort(message = "File: {.val {file_name}} not found.")
  }

  # Check feature_slot_name is acceptable
  if (!feature_slot_name %in% c("features", "genes")) {
    cli_abort(message = c("{.code feature_slot_name} must be one of {.val features} or {.val genes}.",
                               "i" = "If unsure, check contents of H5 file {.code rhdf5::h5ls('{file_name}')}."))
  }

  # Read file
  infile <- hdf5r::H5File$new(filename = file_name, mode = "r")

  # Get list of H5 contents
  h5_dataset_list <- hdf5r::list.datasets(infile)

  # Check feature_slot_name is correct
  if (!length(x = grep(pattern = feature_slot_name, x = h5_dataset_list, value = TRUE)) > 0) {
    cli_abort(message = c("{.code feature_slot_name}: {.val {feature_slot_name}} not found in H5 file.",
                               "i" = "Check contents of H5 file {.code rhdf5::h5ls('{file_name}')} to confirm correct {.code feature_slot_name}."))
  }

  # Assign feature slot name
  if (feature_slot_name == "features") {
    if (isTRUE(x = use.names)) {
      feature_slot <- 'features/name'
    }
    else {
      feature_slot <- 'features/id'
    }
  }

  if (feature_slot_name == "genes") {
    if (isTRUE(x = use.names)) {
      feature_slot <- 'gene_names'
    }
    else {
      feature_slot <- 'genes'
    }
  }

  # add name check
  group_names <- names(x = infile)

  if (!is.null(x = h5_group_name) && !h5_group_name %in% group_names) {
    cli_abort(message = c("{.code h5_group_name} {.val {h5_group_name}} not found.",
                               "i" = "Check H5 file group names {.code rhdf5::h5ls('{file_name}')}."))
  }

  # Read in data
  if ("matrix" %in% group_names) {
    counts <- infile[["matrix/data"]]
    indices <- infile[["matrix/indices"]]
    indptr <- infile[["matrix/indptr"]]
    shp <- infile[["matrix/shape"]]
    features <- infile[[paste0("matrix/", feature_slot)]][]
    barcodes <- infile[["matrix/barcodes"]]
  } else {
    if (length(x = group_names) == 1) {
      counts <- infile[[paste0(group_names, '/data')]]
      indices <- infile[[paste0(group_names, '/indices')]]
      indptr <- infile[[paste0(group_names, '/indptr')]]
      shp <- infile[[paste0(group_names, '/shape')]]
      features <- infile[[paste0(group_names, '/', feature_slot)]][]
      barcodes <- infile[[paste0(group_names, '/barcodes')]]
    } else {
      # check subgroups
      if (is.null(x = h5_group_name)) {
        cli_abort(message = c("H5 file contains multiple sub-groups.",
                                   "i" = "Please provide {.code h5_group_name} specifying which subgroup contains count data."))
      } else {
        counts <- infile[[paste0(h5_group_name, '/data')]]
        indices <- infile[[paste0(h5_group_name, '/indices')]]
        indptr <- infile[[paste0(h5_group_name, '/indptr')]]
        shp <- infile[[paste0(h5_group_name, '/shape')]]
        features <- infile[[paste0(h5_group_name, '/', feature_slot)]][]
        barcodes <- infile[[paste0(h5_group_name, '/barcodes')]]
      }
    }
  }

  # Create sparse matrix
  sparse.mat <- sparseMatrix(
    i = indices[] + 1,
    p = indptr[],
    x = as.numeric(x = counts[]),
    dims = shp[],
    repr = "T"
  )

  if (isTRUE(x = unique.features)) {
    features <- make.unique(names = features)
  }

  rownames(x = sparse.mat) <- features
  colnames(x = sparse.mat) <- barcodes[]
  sparse.mat <- as.sparse(x = sparse.mat)

  infile$close_all()

  return(sparse.mat)
}


#' Load CellBender h5 matrices (corrected) from multiple directories
#'
#' Extract sparse matrix with corrected counts from CellBender h5 output file across multiple sample
#' subdirectories.
#'
#' @param base_path path to the parent directory which contains all of the subdirectories of interest.
#' @param secondary_path path from the parent directory to count matrix files for each sample.
#' @param filtered_h5 logical (default TRUE).  Will set the shared file name suffix `custom_name` is NULL.
#' @param custom_name if file name was customized in CellBender then this parameter should contain the portion
#' of file name that is shared across all samples.  Must included the ".h5" extension as well.
#' @param sample_list a vector of sample directory names if only specific samples are desired.  If `NULL` will
#' read in subdirectories in parent directory.
#' @param sample_names a set of sample names to use for each sample entry in returned list.  If `NULL` will
#' set names to the subdirectory name of each sample.  NOTE: unless `sample_list` is specified this will
#' rename files in the order they are read which will be alphabetical.
#' @param no_file_prefix logical, whether or not the file has prefix identical to folder name. Default is TRUE.
#' @param h5_group_name Name of the group within H5 file that contains count data.  This is only
#' required if H5 file contains multiple subgroups and non-default names.  Default is `NULL`.
#' @param feature_slot_name Name of the slot contain feature names/ids.  Must be one of:
#' "features"(Cell Ranger v3+) or "genes" (Cell Ranger v1/v2 or STARsolo).  Default is "features".
#' @param replace_suffix logical (default FALSE).  Whether or not to replace the barcode suffixes of matrices
#' using \code{\link{Replace_Suffix}}.
#' @param new_suffix_list a vector of new suffixes to replace existing suffixes if `replace_suffix = TRUE`.
#' See \code{\link{Replace_Suffix}} for more information.  To remove all suffixes set `new_suffix_list = ""`.
#' @param parallel logical (default FALSE) whether or not to use multi core processing to read in matrices.
#' @param num_cores how many cores to use for parallel processing.
#' @param merge logical (default FALSE) whether or not to merge samples into a single matrix or return
#' list of matrices.  If TRUE each sample entry in list will have cell barcode prefix added.  The prefix
#' will be taken from `sample_names`.
#' @param ... Extra parameters passed to \code{\link[scCustomize]{Read_CellBender_h5_Mat}}.
#'
#' @return list of sparse matrices
#'
#' @import parallel
#' @import pbapply
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @concept read_&_write
#'
#' @examples
#' \dontrun{
#' base_path <- 'path/to/data/directory'
#' mat_list <- Read_CellBender_h5_Multi_Directory(base_path = base_path)
#' }
#'

Read_CellBender_h5_Multi_Directory <- function(
  base_path,
  secondary_path = NULL,
  filtered_h5 = TRUE,
  custom_name = NULL,
  sample_list = NULL,
  sample_names = NULL,
  no_file_prefix = FALSE,
  h5_group_name = NULL,
  feature_slot_name = "features",
  replace_suffix = FALSE,
  new_suffix_list = NULL,
  parallel = FALSE,
  num_cores = NULL,
  merge = FALSE,
  ...
) {
  # Confirm num_cores specified
  if (isTRUE(x = parallel) && is.null(x = num_cores)) {
    cli_abort("If {.code parallel = TRUE} then {.code num_cores} must be specified.")
  }
  # Confirm directory exists
  if (dir.exists(paths = base_path) == FALSE) {
    cli_abort(message = "Directory: {.val {base_path}} specified by {.code base_path} does not exist.")
  }
  # Detect libraries if sample_list is NULL
  if (is.null(x = sample_list)) {
    sample_list <- Pull_Directory_List(base_path = base_path)
  }

  # Add file suffix
  if (!is.null(x = custom_name)) {
    file_suffix <- custom_name

    # check suffix
    file_ext <- grep(x = file_suffix, pattern = ".h5$")
    if (length(x = file_ext) == 0) {
      cli_abort(message = "'custom_name' must end with file extension '.h5'.")
    }
  } else if (isTRUE(x = filtered_h5)) {
    file_suffix <- "_out_filtered.h5"
  } else {
    file_suffix <- "_out.h5"
  }

  # Edit secondary path if NULL
  if (is.null(x = secondary_path)) {
    secondary_path <- ""
  }

  # Check if full directory path exists
  for (i in 1:length(x = sample_list)) {
    full_directory_path <- file.path(base_path, sample_list[i], secondary_path)
    if (dir.exists(paths = full_directory_path) == FALSE) {
      cli_abort(message = "Full Directory does not exist {.val {full_directory_path}} was not found.")
    }
  }

  # read data
  cli_inform(message = "{.field Reading gene expression files.}")
  if (isTRUE(x = parallel)) {
    cli_inform(message = c("NOTE: Progress bars not currently supported for parallel processing.",
                           "NOTE: Parallel processing will not report informative error messages.", "
                           If function fails set {.code parallel = FALSE} and re-run for informative error reporting.\n"))
    # *** Here is where the swap of mclapply or pbmclapply is occuring ***
    raw_data_list <- mclapply(mc.cores = num_cores, 1:length(x = sample_list), function(x) {
      # Create file path
      if (isFALSE(x = no_file_prefix)) {
        file_path <- file.path(base_path, sample_list[x], secondary_path, paste0(sample_list[x], file_suffix))
      } else {
        file_path <- file.path(base_path, sample_list[x], secondary_path, file_suffix)
      }

      # read and return data
      raw_data <- Read_CellBender_h5_Mat(file_name = file_path, h5_group_name = h5_group_name, feature_slot_name = feature_slot_name, ...)
      return(raw_data)
    })
  } else {
    raw_data_list <- pblapply(1:length(x = sample_list), function(x) {
      # Create file path
      if (isFALSE(x = no_file_prefix)) {
        file_path <- file.path(base_path, sample_list[x], secondary_path, paste0(sample_list[x], file_suffix))
      } else {
        file_path <- file.path(base_path, sample_list[x], secondary_path, file_suffix)
      }

      # read and return data
      raw_data <- Read_CellBender_h5_Mat(file_name = file_path, h5_group_name = h5_group_name, feature_slot_name = feature_slot_name, ...)
      return(raw_data)
    })
  }
  # Name the list items
  if (is.null(x = sample_names)) {
    names(x = raw_data_list) <- sample_list
  } else {
    names(x = raw_data_list) <- sample_names
  }

  # Replace Suffixes
  if (isTRUE(x = replace_suffix)) {
    if (is.null(x = new_suffix_list)) {
      cli_abort(message = "No values provided to {.code new_suffix_list} but {.code replace_suffix = TRUE}.")
    }

    current_suffix_list <- sapply(1:length(raw_data_list), function(x) {
      unique(str_extract(string = colnames(x = raw_data_list[[x]]), pattern = "-.$"))
    })

    if (length(x = new_suffix_list) != 1 & length(x = new_suffix_list) != length(x = current_suffix_list)) {
      cli_abort(message = c("`new_suffix_list` must be either single value or list of values equal to the number of samples.",
                            "i" = "Number of samples is: {.field {length(current_suffix_list)}} and number of new_suffixes provided is: {.field {length(x = new_suffix_list)}}."))
    }

    raw_data_list <- Replace_Suffix(data = raw_data_list, current_suffix = current_suffix_list, new_suffix = new_suffix_list)
  }

  # Merge data
  if (isTRUE(x = merge)) {
    raw_data_merged <- Merge_Sparse_Data_All(matrix_list = raw_data_list, add_cell_ids = names(x = raw_data_list))
    return(raw_data_merged)
  }
  return(raw_data_list)
}


#' Load CellBender h5 matrices (corrected) from multiple files
#'
#' Extract sparse matrix with corrected counts from CellBender h5 output file across multiple samples
#' within the same directory.
#'
#' @param data_dir Directory containing the .h5 files output by CellBender.
#' @param filtered_h5 logical (default TRUE).  Will set the shared file name suffix if `custom_name` is NULL.
#' @param custom_name if file name was customized in CellBender then this parameter should contain the portion
#' of file name that is shared across all samples.  Must included the ".h5" extension as well.
#' @param sample_list a vector of sample names if only specific samples are desired.  If `NULL` will
#' read in all files within `data_dir` directory.
#' @param sample_names a set of sample names to use for each sample entry in returned list.  If `NULL` will
#' set names to the subdirectory name of each sample.
#' @param h5_group_name Name of the group within H5 file that contains count data.  This is only
#' required if H5 file contains multiple subgroups and non-default names.  Default is `NULL`.
#' @param feature_slot_name Name of the slot contain feature names/ids.  Must be one of:
#' "features"(Cell Ranger v3+) or "genes" (Cell Ranger v1/v2 or STARsolo).  Default is "features".
#' @param parallel logical (default FALSE) whether or not to use multi core processing to read in matrices
#' @param num_cores how many cores to use for parallel processing.
#' @param merge logical (default FALSE) whether or not to merge samples into a single matrix or return
#' list of matrices.  If TRUE each sample entry in list will have cell barcode prefix added.  The prefix
#' will be taken from `sample_names`.
#' @param ... Extra parameters passed to \code{\link[scCustomize]{Read_CellBender_h5_Mat}}.
#'
#' @return list of sparse matrices
#'
#' @import parallel
#' @import pbapply
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @concept read_&_write
#'
#' @examples
#' \dontrun{
#' base_path <- 'path/to/data/directory'
#' mat_list <- Read_CellBender_h5_Multi_File(data_dir = base_path)
#' }
#'

Read_CellBender_h5_Multi_File <- function(
  data_dir = NULL,
  filtered_h5 = TRUE,
  custom_name = NULL,
  sample_list = NULL,
  sample_names = NULL,
  h5_group_name = NULL,
  feature_slot_name = "features",
  parallel = FALSE,
  num_cores = NULL,
  merge = FALSE,
  ...
) {
  if (!dir.exists(paths = data_dir)) {
    cli_abort(message = "Directory provided does not exist")
  }
  if (length(x = data_dir) > 1) {
    cli_abort(message = "{.code Read_CellBender_h5_Multi_File} only supports reading from single data directory at a time.")
  }

  # Confirm num_cores specified
  if (isTRUE(x = parallel) && is.null(x = num_cores)) {
    cli_abort("If {.code parallel = TRUE} then {.code num_cores} must be specified.")
  }

  # Add file suffix
  if (!is.null(x = custom_name)) {
    file_suffix <- custom_name

    # check suffix
    file_ext <- check_extension(file_name = file_suffix, extension = ".h5")
    if (isFALSE(x = file_ext)) {
      cli_abort(message = "'custom_name' must end with file extension '.h5'.")
    }
  } else if (isTRUE(x = filtered_h5)) {
    file_suffix <- "_out_filtered.h5"
  } else {
    file_suffix <- "_out.h5"
  }

  file.list <- list.files(path = data_dir, pattern = file_suffix, full.names = FALSE)
  # Remove "barcodes.tsv.gz" file suffix
  if (is.null(x = sample_list)) {
    sample_list <- gsub(pattern = file_suffix, x = file.list, replacement = "")
  }

  # Check sample_names length is OK
  if (!is.null(x = sample_names) && length(x = sample_names) != length(x = sample_list)) {
    cli_abort(message = "Length of {.code sample_names} {.field {length(x = sample_names)}} must be equal to number of samples {.field {length(x = sample_list)}}.")
  }

  cli_inform(message = "{.field Reading Cell Bender H5 files from directory}")
  pboptions(char = "=")
  if (isTRUE(x = parallel)) {
    cli_inform(message = c("NOTE: Progress bars not currently supported for parallel processing.",
                           "NOTE: Parallel processing will not report informative error messages.", "
                           If function fails set {.code parallel = FALSE} and re-run for informative error reporting.\n"))
    raw_data_list <- mclapply(mc.cores = num_cores, 1:length(sample_list), function(i) {
      h5_loc <- file.path(data_dir, paste0(sample_list[i], file_suffix))
      data <- Read_CellBender_h5_Mat(file_name = h5_loc, h5_group_name = h5_group_name, feature_slot_name = feature_slot_name, ...)
    })
  } else {
    raw_data_list <- pblapply(1:length(x = sample_list), function(i) {
      h5_loc <- file.path(data_dir, paste0(sample_list[i], file_suffix))
      data <- Read_CellBender_h5_Mat(file_name = h5_loc, h5_group_name = h5_group_name, feature_slot_name = feature_slot_name, ...)
    })
  }

  # Name the matrices
  if (is.null(x = sample_names)) {
    names(x = raw_data_list) <- sample_list
  } else {
    names(x = raw_data_list) <- sample_names
  }

  # Merge data
  if (isTRUE(x = merge)) {
    raw_data_merged <- Merge_Sparse_Data_All(matrix_list = raw_data_list, add_cell_ids = names(x = raw_data_list))
    return(raw_data_merged)
  }

  # return object
  return(raw_data_list)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### READ OTHER DATA ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Read Overall Statistics from 10X Cell Ranger Count
#'
#' Get data.frame with all metrics from the Cell Ranger count analysis (present in web_summary.html)
#'
#' @param base_path path to the parent directory which contains all of the sub-directories of interest or
#' alternatively can provide single csv file to read and format identically to reading multiple files.
#' @param secondary_path path from the parent directory to count "outs/" folder which contains the
#' "metrics_summary.csv" file.
#' @param default_10X logical (default TRUE) sets the secondary path variable to the default 10X directory structure.
#' @param cellranger_multi logical, whether or not metrics come from Cell Ranger `count` or from Cell Ranger `multi`.  Default is FALSE.
#' @param lib_list a list of sample names (matching directory names) to import.  If `NULL` will read
#' in all samples in parent directory.
#' @param lib_names a set of sample names to use for each sample.  If `NULL` will set names to the
#' directory name of each sample.
#'
#' @return A data frame or list of data.frames with sample metrics from cell ranger.
#'
#' @import pbapply
#' @importFrom dplyr bind_rows
#' @importFrom magrittr "%>%"
#' @importFrom purrr compact
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv
#'
#' @export
#'
#' @concept read_&_write
#'
#' @examples
#' \dontrun{
#' metrics <- Read_Metrics_10X(base_path = "/path/to/directories", default_10X = TRUE)
#' }
#'

Read_Metrics_10X <- function(
  base_path,
  secondary_path = NULL,
  default_10X = TRUE,
  cellranger_multi = FALSE,
  lib_list = NULL,
  lib_names = NULL
) {
  # Check if single file
  file_ending <- grep(pattern = ".csv$", x = base_path, value = TRUE)
  if (length(x = file_ending) == 1) {
    temp_csv <- read.csv(file = base_path)
    if (ncol(x = temp_csv) > nrow(x = temp_csv)) {
      metrics_data <- Metrics_Single_File(base_path = base_path, cellranger_multi = cellranger_multi)
    } else {
      metrics_data <- Metrics_Single_File_v9plus(base_path = base_path, cellranger_multi = cellranger_multi)
    }
    return(metrics_data)
  }

  # Confirm directory exists
  if (dir.exists(paths = base_path) == FALSE) {
    cli_abort(message = "Directory: {.val {base_path}} specified by {.code base_path} does not exist.")
  }
  # Detect libraries if lib_list is NULL
  if (is.null(x = lib_list)) {
    lib_list <- list.dirs(path = base_path, full.names = F, recursive = F)
  }

  # Add file path for 10X default directories
  if (isTRUE(x = default_10X) && !is.null(x = secondary_path)) {
    cli_abort(message = "If {.code default_10X_path = TRUE} then {.code secondary_path} must be NULL.")
  }
  if (isTRUE(x = default_10X)) {
    if (isTRUE(x = cellranger_multi)) {
      secondary_path <- "outs/per_sample_outs/"
    } else {
      secondary_path <- "outs/"
    }
  }
  if (is.null(x = secondary_path)) {
    secondary_path <- ""
  }
  # Check if full directory path exists
  for (i in 1:length(x = lib_list)) {
    full_directory_path <- file.path(base_path, lib_list[i], secondary_path)
    if (dir.exists(paths = full_directory_path) == FALSE) {
      cli_abort(message = "Full Directory does not exist {.val {full_directory_path}} was not found.")
    }
  }

  if (isTRUE(x = cellranger_multi)) {
    if (is.null(x = secondary_path)) {
      s1_file_path <- file.path(base_path, lib_list[1])
    } else {
      s1_file_path <- file.path(base_path, lib_list[1], secondary_path, lib_list[1])
    }

    modalities <- read.csv(file = file.path(s1_file_path, "metrics_summary.csv"), stringsAsFactors = F)$Library.Type %>%
      unique()

    if ("Gene Expression" %in% modalities) {
      multi_gex_metrics <- Metrics_Multi_GEX(lib_list = lib_list, base_path = base_path, secondary_path = secondary_path, lib_names = lib_names)
    }

    if ("VDJ T" %in% modalities) {
      multi_vdjt_metrics <- Metrics_Multi_VDJT(lib_list = lib_list, base_path = base_path, secondary_path = secondary_path, lib_names = lib_names)
    } else {
      multi_vdjt_metrics <- NULL
    }

    if ("VDJ B" %in% modalities) {
      multi_vdjb_metrics <- Metrics_Multi_VDJB(lib_list = lib_list, base_path = base_path, secondary_path = secondary_path, lib_names = lib_names)
    } else {
      multi_vdjb_metrics <- NULL
    }

    if ("Antibody Capture" %in% modalities) {
      multi_abc_metrics <- Metrics_Multi_ABC(lib_list = lib_list, base_path = base_path, secondary_path = secondary_path, lib_names = lib_names)
    } else {
      multi_abc_metrics <- NULL
    }

    # Return data
    data_list <- purrr::compact(
      list(
        multi_gex_metrics = multi_gex_metrics,
        multi_vdjt_metrics = multi_vdjt_metrics,
        multi_vdjb_metrics = multi_vdjb_metrics,
        multi_abc_metrics = multi_abc_metrics)
    )

    return(data_list)
  } else {
    temp_csv <- read.csv(file = file.path(base_path, lib_list[1], secondary_path, "metrics_summary.csv"))
    if (ncol(x = temp_csv) > nrow(x = temp_csv)) {
      count_gex_metrics <- Metrics_Count_GEX(lib_list = lib_list, base_path = base_path, secondary_path = secondary_path, lib_names = lib_names)
    } else {
      count_gex_metrics <- Metrics_Count_GEX_v9plus(lib_list = lib_list, base_path = base_path, secondary_path = secondary_path, lib_names = lib_names)
    }
    return(count_gex_metrics)
  }
}


#' Read Overall Statistics from CellBender
#'
#' Get data.frame with all metrics from the CellBender `remove-background` analysis.
#'
#' @param base_path path to the parent directory which contains all of the sub-directories of interest or
#' path to single metrics csv file.
#' @param lib_list a list of sample names (matching directory names) to import.  If `NULL` will read
#' in all samples in parent directory.
#' @param lib_names a set of sample names to use for each sample.  If `NULL` will set names to the
#' directory name of each sample.
#'
#' @return A data frame with sample metrics from CellBender.
#'
#' @import pbapply
#' @importFrom dplyr bind_rows
#' @importFrom magrittr "%>%"
#' @importFrom utils read.csv
#'
#' @export
#'
#' @concept read_&_write
#'
#' @examples
#' \dontrun{
#' CB_metrics <- Read_Metrics_CellBender(base_path = "/path/to/directories")
#' }
#'

Read_Metrics_CellBender <- function(
    base_path,
    lib_list = NULL,
    lib_names = NULL
) {
  # single file vs. multi-sample
  if (length(x = grep(pattern = "\\.csv", x = base_path, value = TRUE)) > 0) {
    if (file.exists(base_path) == FALSE) {
      cli_abort(message = "Metrics file: {.val {base_path}} specified by {.code base_path} does not exist.")
    } else {
      # read in metrics file
      raw_data_single <- read.csv(file = base_path, stringsAsFactors = FALSE, header = FALSE)

      # Move statistic names to rownames and transpose
      raw_data_single <- raw_data_single %>%
        column_to_rownames("V1") %>%
        t() %>%
        data.frame()

      # Add sample name
      file_name <- basename(path = base_path)
      sample_name <- gsub(pattern = "_out_metrics.csv", replacement = "", x = file_name)

      rownames(raw_data_single) <- sample_name

      raw_data_single <- cbind(sample_id = sample_name, raw_data_single)

      # return the new raw_data data.frame
      return(raw_data_single)
    }
  }

  # Confirm directory exists
  if (dir.exists(paths = base_path) == FALSE) {
    cli_abort(message = "Directory: {.val {base_path}} specified by {.code base_path} does not exist.")
  }
  # Detect libraries if lib_list is NULL
  if (is.null(x = lib_list)) {
    lib_list <- list.dirs(path = base_path, full.names = F, recursive = F)
  }

  # Check if full directory path exists
  for (i in 1:length(x = lib_list)) {
    full_directory_path <- file.path(base_path, lib_list[i])
    if (dir.exists(paths = full_directory_path) == FALSE) {
      cli_abort(message = "Full Directory does not exist {.val {full_directory_path}} was not found.")
    }
  }

  cli_inform(message = "Reading {.field CellBender} Metrics for {.field {length(lib_list)} samples}.")
  raw_data_list <- pblapply(1:length(x = lib_list), function(x) {
    # get directory path
    file_path <- file.path(base_path, lib_list[x])

    # full path with file name
    full_path <- file.path(file_path, paste0(lib_list[x], "_out_metrics.csv"))

    # read in metrics file
    raw_data <- read.csv(file = full_path, stringsAsFactors = FALSE, header = FALSE)

    # Move statistic names to rownames and transpose
    raw_data <- raw_data %>%
      column_to_rownames("V1") %>%
      t() %>%
      data.frame()

    # return the new raw_data data.frame
    raw_data
  })

  # Name the list items
  if (is.null(x = lib_names)) {
    names(x = raw_data_list) <- lib_list
  } else {
    names(x = raw_data_list) <- lib_names
  }

  # Combine the list and add sample_id column
  full_data <- bind_rows(raw_data_list, .id = "sample_id")

  # replace nonsense with sample_id as well
  rownames(x = full_data) <- full_data$sample_id

  return(full_data)
}


#' Read and add results from cNMF
#'
#' Reads the usage and spectra files from cNMF results and adds them as dimensionality
#' reduction to seurat object.
#'
#' @param seurat_object Seurat object name to add cNMF reduction
#' @param usage_file path and name of cNMF usage file
#' @param spectra_file path and name of cNMF spectra file
#' @param reduction_name name to use for reduction to be added, default is "cnmf".
#' @param reduction_key key to use for reduction to be added, default is "cNMF_".
#' @param normalize logical, whether to normalize the cNMF usage data, default is TRUE
#' @param assay assay to add reduction.  Default is NULL and will use current
#' active assay.
#' @param overwrite logical, whether to overwrite a reduction with the name `reduction_name` already
#' present in reduction slot of given Seurat object.
#'
#' @return Seurat object with new dimensionality reduction "cnmf"
#'
#' @importFrom data.table fread
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#'
#' @references For more information about cNMF and usage see \url{https://github.com/dylkot/cNMF}
#'
#' @export
#'
#' @concept read_&_write
#'
#' @examples
#' \dontrun{
#' object <- Read_cNMF(seurat_object = object,
#' usage_file = "example_cNMF/example_cNMF.usages.k_27.dt_0_01.consensus.txt",
#' spectra_file = "example_cNMF/example_cNMF.gene_spectra_score.k_27.dt_0_01.txt")
#' }
#'

Read_Add_cNMF <- function(
    seurat_object,
    usage_file,
    spectra_file,
    reduction_name = "cnmf",
    reduction_key = "cNMF_",
    normalize = TRUE,
    assay = NULL,
    overwrite = FALSE
) {
  # check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # check reduction present
  if (reduction_name %in% Reductions(object = seurat_object) && isFALSE(x = overwrite)) {
    cli_abort(message = c("A reduction with name {.field {reduction_name}} is already present in Seurat object.",
                          "i" = "To run function and overwrite existing reduction set parameter {.code overwrite = TRUE}",
                          "i" = "To keep existing reduction and add new reduction change respective {.code reduction_name} and {.code reduction_key} parameters."))
  }

  # set assay (if null set to active assay)
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Read and transform usage data
  usage_data <- fread(usage_file, data.table = FALSE, header = TRUE) %>%
    column_to_rownames("V1")

  colnames(x = usage_data) <- paste0("Factor_", colnames(x = usage_data))

  # normalize usage data
  if (isTRUE(x = normalize)) {
    usage_data <- as.data.frame(t(apply(usage_data, 1, function(x) x / sum(x))))
  }

  # transform to matrix for dimreduc creation
  usage_data <- as.matrix(x = usage_data)

  # read spectra data
  spectra_data <- fread(spectra_file, data.table = FALSE, header = TRUE) %>%
    column_to_rownames("V1") %>%
    t()

  colnames(spectra_data) <- paste0("Factor_", colnames(spectra_data))

  # create and add dimreduc object
  cnmf_dimreduc <- CreateDimReducObject(embeddings = usage_data, loadings = spectra_data, assay = assay, key = reduction_key)

  seurat_object[[reduction_name]] <- cnmf_dimreduc

  # return object
  return(seurat_object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### READ Utilities ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Pull Directory List
#'
#' Enables easy listing of all sub-directories for use as input library lists in Read10X multi functions.
#'
#' @param base_path path to the parent directory which contains all of the subdirectories of interest.
#'
#' @return A vector of sub-directories within `base_path`.
#'
#' @export
#'
#' @concept read_&_write
#'
#' @examples
#' \dontrun{
#' data_dir <- 'path/to/data/directory'
#' library_list <- Pull_Directory_List(base_path = data_dir)
#' }
#'

Pull_Directory_List <- function(
  base_path
) {
  # Confirm directory exists
  if (dir.exists(paths = base_path) == FALSE) {
    cli_abort(message = "Directory: {.val base_path} specified by {.code base_path} does not exist.")
  }

  # Pull sub-directory list
  dir_list <- list.dirs(path = base_path, full.names = F, recursive = F)
  return(dir_list)
}
