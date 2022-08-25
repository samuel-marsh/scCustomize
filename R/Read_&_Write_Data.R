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
  DropletUtils_check <- PackageCheck("DropletUtils", error = FALSE)
  if (!DropletUtils_check[1]) {
    stop(
      "Please install the DropletUtils package to use Create_10X_H5",
      "\nThis can be accomplished with the following commands: ",
      "\n----------------------------------------",
      "\ninstall.packages('BiocManager')",
      "\nBiocManager::install('DropletUtils')",
      "\n----------------------------------------",
      call. = FALSE
    )
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
  message("Import complete. Start write to H5")
  temp_file <- tempfile(pattern = paste(save_name, "_", sep = ""),
                        tmpdir = save_file_path,
                        fileext=".h5")
  write10xCounts(path = temp_file,
                 x = count_matrix,
                 barcodes = colnames(count_matrix),
                 gene.symbol = rownames(count_matrix),
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
#' @param min_cells value to supply to min.cells parameter of \code{\link[SeuratObject]{CreateSeuratObject}}.
#' Default is 5.
#' @param min_features value to supply to min.features parameter of \code{\link[SeuratObject]{CreateSeuratObject}}.
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
#' seurat_obj <- Create_CellBender_Merged_Seurat(raw_cell_bender_matrix = cb_matrix, raw_counts_matrix = cr_matrix)
#' }
#'

Create_CellBender_Merged_Seurat <- function(
  raw_cell_bender_matrix = NULL,
  raw_counts_matrix = NULL,
  raw_assay_name = "RAW",
  min_cells = 5,
  min_features = 200,
  ...
) {
  # Filter Cell Bender matrix for Cell Ranger cells
  cell_intersect <- intersect(x = colnames(x = raw_counts_matrix), y = colnames(raw_cell_bender_matrix))

  message("Filtering Cell Bender matrix for cells present in raw counts matrix.")

  raw_cell_bender_matrix <- raw_cell_bender_matrix[, cell_intersect]

  # Create Seurat Object
  message("Creating Seurat Object from Cell Bender matrix.")
  cell_bender_seurat <- CreateSeuratObject(counts = raw_cell_bender_matrix, min.cells = min_cells, min.features = min_features, ...)

  # Pull cell and gene names
  cell_names_seurat <- colnames(x = cell_bender_seurat)
  gene_names_seurat <- rownames(x = cell_bender_seurat)

  # Filter raw counts by created Seurat parameters
  message("Filtering raw counts matrix to match Seurat Object.")
  raw_counts_matrix <- raw_counts_matrix[gene_names_seurat, cell_names_seurat]

  # Create raw counts assay object
  message("Creating raw counts Seurat Assay Object.")
  counts <- CreateAssayObject(counts = raw_counts_matrix, min.cells = 0, min.features = 0)

  # Add counts assay to Seurat Object
  message("Adding assay to Seurat Object.")
  cell_bender_seurat[[raw_assay_name]] <- counts

  return(cell_bender_seurat)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### READ 10X DATA ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Load in NCBI GEO data from 10X
#'
#' Enables easy loading of sparse data matrices provided by 10X genomics. That have file prefixes added to them by NCBI GEO or other repos.
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
#'
#' @return If features.csv indicates the data has multiple data types, a list
#'   containing a sparse matrix of the data from each type will be returned.
#'   Otherwise a sparse matrix containing the expression data will be returned.
#'
#' @references Code used in function has been slightly modified from `Seurat::Read10X` function of
#' Seurat package (https://github.com/satijalab/seurat) (Licence: GPL-3).  Function was modified to
#' support file prefixes and altered loop by Samuel Marsh for scCustomize (also previously posted as
#' potential PR to Seurat GitHub).
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
  num_cores = NULL#,
  #single.matrix = FALSE
) {
  if (!dir.exists(paths = data_dir)) {
    stop("Directory provided does not exist")
  }
  if (length(x = data_dir) > 1) {
    stop("Read10X_GEO only supports reading from single data directory at a time.")
  }

  # Confirm num_cores specified
  if (parallel && is.null(x = num_cores)) {
    stop("If 'parallel = TRUE' then 'num_cores' must be specified.")
  }

  if (is.null(x = sample_list)) {
    # pull all file names from directory
    file.list <- list.files(path = data_dir, pattern = "barcodes.tsv", full.names = FALSE)
    # Remove "barcodes.tsv.gz" file suffix
    sample_list <- gsub(pattern = "barcodes.tsv.gz", x = file.list, replacement = "")
  }

  if (!is.null(x = sample_names)) {
    if (length(x = sample_names) != length(x = sample_list)) {
      stop("The length of 'sample_names' provided (", length(x = sample_names), ") does not equal length of 'sample_list' (", length(x = sample_list), ").")
    }
  }


  message("Reading 10X files from directory")
  pboptions(char = "=")
  if (parallel) {
    message("NOTE: Progress bars not currently supported for parallel processing.\n",
            "NOTE: Parallel processing will not report informative error messages.  If function fails set 'parallel = FALSE' and re-run for informative error reporting.\n")
    raw_data_list <- mclapply(mc.cores = num_cores, 1:length(sample_list), function(i) {
      barcode.loc <- file.path(data_dir, paste0(sample_list[i], 'barcodes.tsv.gz'))
      gene.loc <- file.path(data_dir, paste0(sample_list[i], 'genes.tsv.gz'))
      features.loc <- file.path(data_dir, paste0(sample_list[i], 'features.tsv.gz'))
      matrix.loc <- file.path(data_dir, paste0(sample_list[i], 'matrix.mtx.gz'))
      # Flag to indicate if this data is from CellRanger >= 3.0
      pre_ver_3 <- file.exists(gene.loc)
      if (!file.exists(barcode.loc)) {
        stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
      }
      if (!pre_ver_3 && !file.exists(features.loc) ) {
        stop("Gene name or features file missing. Expecting ", basename(path = features.loc))
      }
      if (!file.exists(matrix.loc)) {
        stop("Expression matrix file missing. Expecting ", basename(path = matrix.loc))
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
      if (unique.features) {
        fcols = ncol(x = feature.names)
        if (fcols < gene.column) {
          stop(paste0("gene.column was set to ", gene.column,
                      " but feature.tsv.gz (or genes.tsv) only has ", fcols, " columns.",
                      " Try setting the gene.column argument to a value <= to ", fcols, "."))
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
        stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
      }
      if (!pre_ver_3 && !file.exists(features.loc) ) {
        stop("Gene name or features file missing. Expecting ", basename(path = features.loc))
      }
      if (!file.exists(matrix.loc)) {
        stop("Expression matrix file missing. Expecting ", basename(path = matrix.loc))
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
      if (unique.features) {
        fcols = ncol(x = feature.names)
        if (fcols < gene.column) {
          stop(paste0("gene.column was set to ", gene.column,
                      " but feature.tsv.gz (or genes.tsv) only has ", fcols, " columns.",
                      " Try setting the gene.column argument to a value <= to ", fcols, "."))
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
    names(raw_data_list) <- sample_names
  } else {
    names(raw_data_list) <- sample_list
  }

  # return list
  return(raw_data_list)
}


#' Load in NCBI GEO data from 10X in HDF5 file format
#'
#' Enables easy loading of HDF5 data matrices provided by 10X genomics. That have file prefixes added to
#' them by NCBI GEO or other repos or programs (i.e. Cell Bender)
#'
#' @param data_dir Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv
#' files provided by 10X.
#' @param sample_list A vector of file prefixes/names if specific samples are desired.  Default is `NULL` and
#' will load all samples in given directory.
#' @param sample_names a set of sample names to use for each sample entry in returned list.  If `NULL`
#' will set names to the file name of each sample.
#' @param parallel logical (default FALSE).  Whether to use multiple cores when reading in data.
#' Only possible on Linux based systems.
#' @param num_cores if `parallel = TRUE` indicates the number of cores to use for multicore processing.
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
  ...
) {
  if (!dir.exists(paths = data_dir)) {
    stop("Directory provided does not exist")
  }
  if (length(x = data_dir) > 1) {
    stop("Read10X_GEO only supports reading from single data directory at a time.")
  }

  # Confirm num_cores specified
  if (parallel && is.null(x = num_cores)) {
    stop("If 'parallel = TRUE' then 'num_cores' must be specified.")
  }

  file.list <- list.files(path = data_dir, pattern = ".h5", full.names = FALSE)
  # Remove "barcodes.tsv.gz" file suffix
  if (is.null(x = sample_list)) {
    if (is.null(x = shared_suffix)) {
      sample_list <- gsub(pattern = ".h5", x = file.list, replacement = "")
    } else {
      sample_list <- gsub(pattern = paste0(shared_suffix, ".h5"), x = file.list, replacement = "")
    }
  }

  # Check sample_names length is ok
  if (!is.null(x = sample_names) && length(x = sample_names) != length(x = sample_list)) {
    stop("Length of `sample_names` must be equal to number of samples.")
  }

  message("Reading 10X H5 files from directory")
  pboptions(char = "=")
  if (parallel) {
    message("NOTE: Progress bars not currently supported for parallel processing.\n",
            "NOTE: Parallel processing will not report informative error messages.  If function fails set 'parallel = FALSE' and re-run for informative error reporting.\n")
    raw_data_list <- mclapply(mc.cores = num_cores, 1:length(sample_list), function(i) {
      h5_loc <- file.path(data_dir, paste0(sample_list[i], shared_suffix, ".h5"))
      data <- Read10X_h5(filename = h5_loc)
    })
  } else {
    raw_data_list <- pblapply(1:length(x = sample_list), function(i) {
      h5_loc <- file.path(data_dir, paste0(sample_list[i], shared_suffix, ".h5"))
      data <- Read10X_h5(filename = h5_loc)
    })
  }

  # Name the matrices
  if (is.null(x = sample_names)) {
    names(raw_data_list) <- sample_list
  } else {
    names(raw_data_list) <- sample_names
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
#' @param sample_list a vector of sample directory names if only specific samples are desired.  If `NULL` will
#' read in subdirectories in parent directory.
#' @param sample_names a set of sample names to use for each sample entry in returned list.  If `NULL` will
#' set names to the subdirectory name of each sample.
#' @param parallel logical (default FALSE) whether or not to use multi core processing to read in matrices.
#' @param num_cores how many cores to use for parallel processing.
#' @param merge logical (default FALSE) whether or not to merge samples into a single matrix or return list of
#' matrices.
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
  sample_list = NULL,
  sample_names = NULL,
  parallel = FALSE,
  num_cores = NULL,
  merge = FALSE,
  ...
) {
  # Confirm num_cores specified
  if (parallel && is.null(x = num_cores)) {
    stop("If 'parallel = TRUE' then 'num_cores' must be specified.")
  }
  # Confirm directory exists
  if (dir.exists(paths = base_path) == FALSE) {
    stop(paste0("Directory: ", base_path, "specified by 'base_path' does not exist."))
  }
  # Detect libraries if sample_list is NULL
  if (is.null(x = sample_list)) {
    sample_list <- Pull_Directory_List(base_path = base_path)
  }
  # Add file path for 10X default directories
  if (default_10X_path && !is.null(x = secondary_path)) {
    stop("If 'default_10X_path = TRUE' then 'secondary_path' must be NULL.")
  }
  if (default_10X_path) {
    secondary_path <- "outs/filtered_feature_bc_matrix/"
  }
  if (is.null(x = secondary_path)) {
    secondary_path <- ""
  }
  # Check if full directory path exists
  for (i in 1:length(x = sample_list)) {
    full_directory_path <- file.path(base_path, sample_list[i], secondary_path)
    if (dir.exists(paths = full_directory_path) == FALSE) {
      stop(paste0("Full Directory does not exist: ", full_directory_path, " was not found."))
    }
  }
  # read data
  message("Reading gene expression files.")
  if (parallel) {
    message("NOTE: Progress bars not currently supported for parallel processing.\n",
            "NOTE: Parallel processing will not report informative error messages.  If function fails set 'parallel = FALSE' and re-run for informative error reporting.\n")
    # *** Here is where the swap of mclapply or pbmclapply is occuring ***
    raw_data_list <- mclapply(mc.cores = num_cores, 1:length(x = sample_list), function(x) {
      file_path <- file.path(base_path, sample_list[x], secondary_path)
      raw_data <- Read10X(data.dir = file_path, ...)
      return(raw_data)
    })
  } else {
    raw_data_list <- pblapply(1:length(x = sample_list), function(x) {
      if (is.null(x = secondary_path)) {
        file_path <- file.path(base_path, sample_list[x])
      } else {
        file_path <- file.path(base_path, sample_list[x], secondary_path)
      }
      raw_data <- Read10X(data.dir = file_path, ...)
    })
  }
  # Name the list items
  if (is.null(x = sample_names)) {
    names(raw_data_list) <- sample_list
  } else {
    names(raw_data_list) <- sample_names
  }
  # Merge data
  if (merge) {
    raw_data_merged <- Merge_Sparse_Data_All(matrix_list = raw_data_list, add_cell_ids = names(raw_data_list))
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
#' @param h5_filename name of h5 file (including .h5 suffix).  If all h5 files have same name (i.e. Cell Ranger output)
#' then use full file name.  By default function uses Cell Ranger name: "filtered_feature_bc_matrix.h5".
#' If h5 files have sample specific prefixes (i.e. from Cell Bender) then use only the shared part of file
#' name (e.g., "_filtered_out.h5").
#' @param cell_bender logical (default FALSE).  Is the h5 file from cell bender output, needed to set correct file names.
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
#' list of matrices.  Will use the `sample_names` parameter to add prefix to cell barcodes.
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
  h5_filename = "filtered_feature_bc_matrix.h5",
  cell_bender = FALSE,
  sample_list = NULL,
  sample_names = NULL,
  replace_suffix = FALSE,
  new_suffix_list = NULL,
  parallel = FALSE,
  num_cores = NULL,
  merge = FALSE,
  ...
) {
  # Check cell bender or default 10X
  if (cell_bender && default_10X_path) {
    stop("Both `cell_bender` and `default_10X_path` cannot be simultaneously set to TRUE.")
  }

  # Confirm num_cores specified
  if (parallel && is.null(x = num_cores)) {
    stop("If 'parallel = TRUE' then 'num_cores' must be specified.")
  }
  # Confirm directory exists
  if (dir.exists(paths = base_path) == FALSE) {
    stop(paste0("Directory: ", base_path, "specified by 'base_path' does not exist."))
  }
  # Detect libraries if sample_list is NULL
  if (is.null(x = sample_list)) {
    sample_list <- Pull_Directory_List(base_path = base_path)
  }

  # Add file path for 10X default directories
  if (default_10X_path && !is.null(x = secondary_path)) {
    stop("If 'default_10X_path = TRUE' then 'secondary_path' must be NULL.")
  }
  if (default_10X_path) {
    secondary_path <- "outs/"
  }
  if (is.null(x = secondary_path)) {
    secondary_path <- ""
  }
  # Check if full directory path exists
  for (i in 1:length(x = sample_list)) {
    full_directory_path <- file.path(base_path, sample_list[i], secondary_path)
    if (dir.exists(paths = full_directory_path) == FALSE) {
      stop(paste0("Full Directory does not exist: ", full_directory_path, " was not found."))
    }
  }

  # read data
  message("Reading gene expression files.")
  if (parallel) {
    message("NOTE: Progress bars not currently supported for parallel processing.\n",
            "NOTE: Parallel processing will not report informative error messages.  If function fails set 'parallel = FALSE' and re-run for informative error reporting.\n")
    # *** Here is where the swap of mclapply or pbmclapply is occuring ***
    raw_data_list <- mclapply(mc.cores = num_cores, 1:length(x = sample_list), function(x) {
      if (cell_bender) {
        file_path <- file.path(base_path, sample_list[x], secondary_path, paste0(sample_list[x], h5_filename))
      } else {
        file_path <- file.path(base_path, sample_list[x], secondary_path, h5_filename)
      }
      raw_data <- Read10X_h5(filename = file_path, ...)
      return(raw_data)
    })
  } else {
    raw_data_list <- pblapply(1:length(x = sample_list), function(x) {
      if (is.null(x = secondary_path)) {
        if (cell_bender) {
          file_path <- file.path(base_path, sample_list[x], paste0(sample_list[x], h5_filename))
        } else {
          file_path <- file.path(base_path, sample_list[x], h5_filename)
        }
      } else {
        if (cell_bender) {
          file_path <- file.path(base_path, sample_list[x], secondary_path, paste0(sample_list[x], h5_filename))
        } else {
          file_path <- file.path(base_path, sample_list[x], secondary_path, h5_filename)
        }
      }
      raw_data <- Read10X_h5(filename = file_path, ...)
    })
  }
  # Name the list items
  if (is.null(x = sample_names)) {
    names(raw_data_list) <- sample_list
  } else {
    names(raw_data_list) <- sample_names
  }

  # Replace Suffixes
  if (replace_suffix) {
    if (is.null(x = new_suffix_list)) {
      stop("No values provided to `new_suffix_list` but `replace_suffix = TRUE`.")
    }

    current_suffix_list <- sapply(1:length(raw_data_list), function(x) {
      unique(str_extract(string = colnames(x = raw_data_list[[x]]), pattern = "-.$"))
    })

    if (length(x = new_suffix_list) != 1 & length(x = new_suffix_list) != length(x = current_suffix_list)) {
      stop("`new_suffix_list` must be either single value or list of values equal to the number of samples.  Number of samples is: ", length(current_suffix_list), " and number of new_suffixes provided is:", length(x = new_suffix_list), ".")
    }

    raw_data_list <- Replace_Suffix(data = raw_data_list, current_suffix = current_suffix_list, new_suffix = new_suffix_list)
  }

  # Merge data
  if (merge) {
    raw_data_merged <- Merge_Sparse_Data_All(matrix_list = raw_data_list, add_cell_ids = names(raw_data_list))
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
#'
#' @return List of gene x cell matrices in list format named by sample name.
#'
#' @import Matrix
#' @import parallel
#' @import pbapply
#' @importFrom data.table fread
#' @importFrom magrittr "%>%"
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
  num_cores = NULL
) {
  # Create list of all files in directory
  possible_file_list <- list.files(path = data_dir, pattern = file_suffix, full.names = FALSE)

  # Check files found
  if (is.null(x = possible_file_list)) {
    stop("No files found.  Check that `data_dir` and `file_suffix` are correct.")
  }

  # Set all files to be used if sample_list is NULL
  if (is.null(sample_list)) {
    file_list <- possible_file_list
  }

  # Confirm num_cores specified
  if (parallel && is.null(x = num_cores)) {
    stop("If 'parallel = TRUE' then 'num_cores' must be specified.")
  }

  # Read in subset of files
  if (!is.null(x = sample_list)) {
    # Add suffix
    if (full_names) {
      file_list <- sample_list
    } else {
      file_list <- paste0(sample_list, file_suffix)
    }
    file_list <- file_list

    if (any(!file_list %in% possible_file_list)) {
      bad_file_list <- file_list[!file_list %in% possible_file_list]
      file_list <- file_list[file_list %in% possible_file_list]
      if (length(x = file_list) == 0) {
        stop("No requested files found. Check that 'data_dir' and file_suffix' are correct \n
             and `full_names` parameter is accurate.")
      }
      warning("The following files were not imported as they were not found in specified directory",
              ": ", glue_collapse_scCustom(input_string = bad_file_list, and = TRUE))
    }
  }

  # Get sample names
  if (is.null(x = sample_names)) {
    sample_names <- gsub(pattern = file_suffix, x = file_list, replacement = "")
  } else {
    sample_names <- sample_names
  }

  # Read in files
  message("Reading gene expression files from directory")
  pboptions(char = "=")
  if (parallel) {
    message("NOTE: Progress bars not currently supported for parallel processing.\n",
            "NOTE: Parallel processing will not report informative error messages.  If function fails set 'parallel = FALSE' and re-run for informative error reporting.\n")
    raw_data_list <- mclapply(mc.cores = num_cores, 1:length(x = file_list), function(i) {
      dge_loc <- file.path(data_dir, file_list[i])
      data <- fread(file = dge_loc, data.table = F)
      if (move_genes_rownames) {
        first_col_name <- colnames(data[1])
        data <- data %>%
          column_to_rownames(first_col_name)
      }
      if (barcode_suffix_period) {
        colnames(data) <- gsub("\\.", "-", colnames(data))
      }
      data_sparse <- as(data, "Matrix")
      return(data_sparse)
    })
  } else {
    raw_data_list <- pblapply(1:length(x = file_list), function(i) {
      dge_loc <- file.path(data_dir, file_list[i])
      data <- fread(file = dge_loc, data.table = F)
      if (move_genes_rownames) {
        first_col_name <- colnames(data[1])
        data <- data %>%
          column_to_rownames(first_col_name)
      }
      # Check all columns numeric
      col_data_numeric <- sapply(data, is.numeric)
      if (!all(col_data_numeric)) {
        stop("One or more columns in the file: ", '"', dge_loc, '"', " contains non-numeric data.  Please check original file and/or that parameter `move_genes_rownames` is set appropriately.")
      }
      if (barcode_suffix_period) {
        colnames(data) <- gsub("\\.", "-", colnames(data))
      }
      data_sparse <- as(data, "Matrix")
      return(data_sparse)
    })
  }

  # Name the items in list
  names(raw_data_list) <- sample_names
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
#'
#' @return sparse matrix
#'
#' @references Code used in function has been modified from `Seurat::Read10X_h5` function of
#' Seurat package (https://github.com/satijalab/seurat) (Licence: GPL-3).
#'
#' @import cli
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
  unique.features = TRUE
) {
  # Check hdf5r installed
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    cli_abort(message = c("Please install hdf5r to read HDF5 files",
                          "i" = "`install.packages('hdf5r')`")
    )
  }
  # Check file
  if (!file.exists(file_name)) {
    cli_abort(message = "File: {file_name} not found.")
  }

  if (use.names) {
    feature_slot <- 'features/name'
  } else {
    feature_slot <- 'features/id'
  }

  # Read file
  infile <- hdf5r::H5File$new(filename = file_name, mode = "r")

  counts <- infile[["matrix/data"]]
  indices <- infile[["matrix/indices"]]
  indptr <- infile[["matrix/indptr"]]
  shp <- infile[["matrix/shape"]]
  features <- infile[[paste0("matrix/", feature_slot)]][]
  barcodes <- infile[["matrix/barcodes"]]


  sparse.mat <- sparseMatrix(
    i = indices[] + 1,
    p = indptr[],
    x = as.numeric(x = counts[]),
    dims = shp[],
    repr = "T"
  )

  if (unique.features) {
    features <- make.unique(names = features)
  }

  rownames(x = sparse.mat) <- features
  colnames(x = sparse.mat) <- barcodes[]
  sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")

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
#' @param filtered_h5 BLANK
#' @param custom_name BLANK
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
#' list of matrices.  Will use the `sample_names` parameter to add prefix to cell barcodes.
#' @param ... Extra parameters passed to \code{\link[scCustomize]{Read_CellBender_h5_Mat}}.
#'
#' @return list of sparse matrices
#'
#' @import cli
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
#' mat <- Read_CellBender_h5_Multi_Directory(base_path = base_path)
#' }
#'

Read_CellBender_h5_Multi_Directory <- function(
  base_path,
  secondary_path = NULL,
  filtered_h5 = TRUE,
  custom_name = NULL,
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
  if (parallel && is.null(x = num_cores)) {
    stop("If 'parallel = TRUE' then 'num_cores' must be specified.")
  }
  # Confirm directory exists
  if (dir.exists(paths = base_path) == FALSE) {
    stop(paste0("Directory: ", base_path, "specified by 'base_path' does not exist."))
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
  } else if (filtered_h5) {
    file_suffix <- "_out_filtered.h5"
  } else {
    file_suffix <- "_out.h5"
  }

  # Check if full directory path exists
  for (i in 1:length(x = sample_list)) {
    full_directory_path <- file.path(base_path, sample_list[i], secondary_path)
    if (dir.exists(paths = full_directory_path) == FALSE) {
      stop(paste0("Full Directory does not exist: ", full_directory_path, " was not found."))
    }
  }

  # read data
  message("Reading gene expression files.")
  if (parallel) {
    message("NOTE: Progress bars not currently supported for parallel processing.\n",
            "NOTE: Parallel processing will not report informative error messages.  If function fails set 'parallel = FALSE' and re-run for informative error reporting.\n")
    # *** Here is where the swap of mclapply or pbmclapply is occuring ***
    raw_data_list <- mclapply(mc.cores = num_cores, 1:length(x = sample_list), function(x) {
      # Create file path
      file_path <- file.path(base_path, sample_list[x], secondary_path, paste0(sample_list[x], file_suffix))

      # read and return data
      raw_data <- Read_CellBender_h5_Mat(file_name = file_path, ...)
      return(raw_data)
    })
  } else {
    raw_data_list <- pblapply(1:length(x = sample_list), function(x) {
      # Create file path
      file_path <- file.path(base_path, sample_list[x], secondary_path, paste0(sample_list[x], file_suffix))

      # read and return data
      raw_data <- Read_CellBender_h5_Mat(file_name = file_path, ...)
      return(raw_data)
    })
  }
  # Name the list items
  if (is.null(x = sample_names)) {
    names(raw_data_list) <- sample_list
  } else {
    names(raw_data_list) <- sample_names
  }

  # Replace Suffixes
  if (replace_suffix) {
    if (is.null(x = new_suffix_list)) {
      stop("No values provided to `new_suffix_list` but `replace_suffix = TRUE`.")
    }

    current_suffix_list <- sapply(1:length(raw_data_list), function(x) {
      unique(str_extract(string = colnames(x = raw_data_list[[x]]), pattern = "-.$"))
    })

    if (length(x = new_suffix_list) != 1 & length(x = new_suffix_list) != length(x = current_suffix_list)) {
      stop("`new_suffix_list` must be either single value or list of values equal to the number of samples.  Number of samples is: ", length(current_suffix_list), " and number of new_suffixes provided is:", length(x = new_suffix_list), ".")
    }

    raw_data_list <- Replace_Suffix(data = raw_data_list, current_suffix = current_suffix_list, new_suffix = new_suffix_list)
  }

  # Merge data
  if (merge) {
    raw_data_merged <- Merge_Sparse_Data_All(matrix_list = raw_data_list, add_cell_ids = names(raw_data_list))
    return(raw_data_merged)
  }
  return(raw_data_list)
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### READ OTHER DATA ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Read Overall Statistics from 10X Cell Ranger Count
#'
#' Get data.frame with all metrics from the Cell Ranger count analysis (present in web_summary.html)
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
#' @return A data frame with sample metrics from cell ranger.
#'
#' @import pbapply
#' @importFrom dplyr bind_rows
#' @importFrom utils txtProgressBar setTxtProgressBar
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
  lib_list = NULL,
  lib_names = NULL
) {
  # Confirm directory exists
  if (dir.exists(paths = base_path) == FALSE) {
    stop(paste0("Directory: ", base_path, "specified by 'base_path' does not exist."))
  }
  # Detect libraries if lib_list is NULL
  if (is.null(x = lib_list)) {
    lib_list <- list.dirs(path = base_path, full.names = F, recursive = F)
  }

  # Add file path for 10X default directories
  if (default_10X && !is.null(x = secondary_path)) {
    stop("If 'default_10X = TRUE' then 'secondary_path' must be NULL.")
  }
  if (default_10X) {
    secondary_path <- "outs/"
  }
  if (is.null(x = secondary_path)) {
    secondary_path <- ""
  }
  # Check if full directory path exists
  for (i in 1:length(x = lib_list)) {
    full_directory_path <- file.path(base_path, lib_list[i], secondary_path)
    if (dir.exists(paths = full_directory_path) == FALSE) {
      stop(paste0("Full Directory does not exist: ", full_directory_path, " was not found."))
    }
  }

  # Read in raw data
  raw_data_list <- pblapply(1:length(x = lib_list), function(x) {
    if (is.null(x = secondary_path)) {
      file_path <- file.path(base_path, lib_list[x])
    } else {
      file_path <- file.path(base_path, lib_list[x], secondary_path)
    }

    raw_data <- read.csv(file = paste0(file_path, "metrics_summary.csv"), stringsAsFactors = F)
    # Change format of numeric columns to due commas in data csv output.
    column_numbers <- grep(pattern = ",", x = raw_data[1, ])
    raw_data[,c(column_numbers)] <- lapply(raw_data[,c(column_numbers)],function(x){as.numeric(gsub(",", "", x))})
    return(raw_data)
  })

  # Name the list items
  if (is.null(x = lib_names)) {
    names(raw_data_list) <- lib_list
  } else {
    names(raw_data_list) <- lib_names
  }

  # Combine the list and add sample_id column
  full_data <- bind_rows(raw_data_list, .id = "sample_id")

  # Change column nams to use "_" separator instead of "." for readability
  colnames(full_data) <- gsub(pattern = "\\.", replacement = "_", x = colnames(x = full_data))

  return(full_data)
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
    stop(paste0("Directory: ", base_path, "specified by 'base_path' does not exist."))
  }

  # Pull sub-directory list
  dir_list <- list.dirs(path = base_path, full.names = F, recursive = F)
  return(dir_list)
}
