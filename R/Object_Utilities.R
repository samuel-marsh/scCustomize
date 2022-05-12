#' Add Mito and Ribo percentages
#'
#' Add Mito, Ribo, & Mito+Ribo percentages to meta.data slot of Seurat Object
#'
#' @param seurat_object object name.
#' @param species Species of origin for given Seurat Object.  If mouse, human, marmoset, zebrafish, rat,
#' drosophila, or rhesus macaque (name or abbreviation) are provided the function will automatically
#' generate mito_pattern and ribo_pattern values.
#' @param mito_name name to use for the new meta.data column containing percent mitochondrial counts.
#' Default is "percent_mito".
#' @param ribo_name name to use for the new meta.data column containing percent ribosomal counts.
#' Default is "percent_mito".
#' @param mito_ribo_name name to use for the new meta.data column containing percent
#' mitochondrial+ribosomal counts.  Default is "percent_mito".
#' @param mito_pattern A regex pattern to match features against for mitochondrial genes (will set automatically if
#' species is mouse or human; marmoset features list saved separately).
#' @param ribo_pattern A regex pattern to match features against for ribosomal genes
#' (will set automatically if species is mouse, human, or marmoset).
#' @param mito_features A list of mitochrondial gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param ribo_features A list of ribosomal gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param assay Assay to use (default is the current object default assay).
#' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
#' function will abort if columns with any one of the names provided to `mito_name` `ribo_name` or
#' `mito_ribo_name` is present in meta.data slot.
#' @param list_species_names returns list of all accepted values to use for default species names which
#' contain internal regex/feature lists (human, mouse, marmoset, zebrafish, rat, drosophila, and
#' rhesus macaque).  Default is FALSE.
#'
#' @importFrom dplyr mutate select intersect
#' @importFrom magrittr "%>%"
#' @importFrom Seurat PercentageFeatureSet AddMetaData
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @return A Seurat Object
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' \dontrun{
#' object <- Add_Mito_Ribo_Seurat(seurat_object = object, species = "mouse")
#' }
#'

Add_Mito_Ribo_Seurat <- function(
  seurat_object,
  species,
  mito_name = "percent_mito",
  ribo_name = "percent_ribo",
  mito_ribo_name = "percent_mito_ribo",
  mito_pattern = NULL,
  ribo_pattern = NULL,
  mito_features = NULL,
  ribo_features = NULL,
  assay = NULL,
  overwrite = FALSE,
  list_species_names = FALSE
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

  # Return list of accepted default species name options
  if (list_species_names) {
    return(accepted_names)
    stop_quietly()
  }

  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Overwrite check
  if (mito_name %in% colnames(x = seurat_object@meta.data) || ribo_name %in% colnames(x = seurat_object@meta.data) || mito_ribo_name %in% colnames(x = seurat_object@meta.data)) {
    if (!overwrite) {
      cli_abort(message = c("Columns with {mito_name} and/or {ribo_name} already present in meta.data slot.",
                            "i" = "*To run function and overwrite columns set parameter `overwrite = TRUE` or change respective 'mito_name', 'ribo_name', or mito_ribo_name'*")
      )
    }
    cli_inform(message = c("Columns with {mito_name} and/or {ribo_name} already present in meta.data slot.",
                           "i" = "Overwriting those columns as overwrite = TRUE.")
    )
  }

  # Checks species
  if (is.null(x = species)) {
    cli_abort(message = c("No species name or abbreivation was provided to `species` parameter.",
                          "i" = "If not using default species please set `species = other`.")
    )
  }

  # Set default assay
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Species Spelling Options
  mouse_options <- accepted_names$Mouse_Options
  human_options <- accepted_names$Human_Options
  marmoset_options <- accepted_names$Marmoset_Options
  zebrafish_options <- accepted_names$Zebrafish_Options
  rat_options <- accepted_names$Rat_Options
  drosophila_options <- accepted_names$Drosophila_Options
  macaque_options <- accepted_names$Macaque_Options

  # Assign mito/ribo pattern to stored species
  if (species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options) && any(!is.null(x = mito_pattern), !is.null(x = ribo_pattern))) {
    cli_warn(message = c("Pattern expressions for included species (Human & Mouse) are set by default.",
                         "*" = "Supplied `mito_pattern` and `ribo_pattern` will be disregarded.",
                         "i" = "To override defaults please supply a feature list for mito and/or ribo genes.")
    )
  }

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
  if (is.null(x = mito_pattern) && is.null(x = mito_features) && is.null(x = ribo_pattern) && is.null(x = ribo_pattern)) {
    cli_abort(message = c("No features or patterns provided for mito/ribo genes.",
                          "i" = "Please provide a default species name or pattern/features."))
  }

  mito_features <- mito_features %||% grep(pattern = mito_pattern, x = rownames(x = seurat_object[[assay]]), value = TRUE)

  ribo_features <- ribo_features %||% grep(pattern = ribo_pattern, x = rownames(x = seurat_object[[assay]]), value = TRUE)

  # Check features are present in object
  length_mito_features <- length(x = intersect(x = mito_features, y = rownames(x = seurat_object[[assay]])))

  length_ribo_features <- length(x = intersect(x = ribo_features, y = rownames(x = seurat_object[[assay]])))

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

  # Add mito and ribo columns
  if (length_mito_features > 0) {
    seurat_object[["percent_mito"]] <- PercentageFeatureSet(object = seurat_object, features = mito_features, assay = assay)
  }
  if (length_ribo_features > 0) {
    seurat_object[["percent_ribo"]] <- PercentageFeatureSet(object = seurat_object, features = ribo_features, assay = assay)
  }

  # Create combined mito ribo column if both present
  if (length_mito_features > 0 && length_ribo_features > 0) {
    object_meta <- seurat_object@meta.data %>%
      rownames_to_column("barcodes")

    object_meta <- object_meta %>%
      mutate(percent_mito_ribo = percent_mito + percent_ribo)

    object_meta <- object_meta %>%
      select(barcodes, percent_mito_ribo) %>%
      column_to_rownames("barcodes")

    seurat_object <- AddMetaData(object = seurat_object, metadata = object_meta)
  }

  # return final object
  return(seurat_object)
}


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
#' @importFrom dplyr select
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#' @importFrom tidyselect all_of
#'
#' @return data.frame with only new columns.
#'
#' @export
#'
#' @concept object_util
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

  if (barcodes_to_rownames) {
    # Check barcodes colname exists
    if (!barcodes_colname %in% colnames(x = meta_data)) {
      stop("barcodes_colname: ", '"', barcodes_colname, '"', " was not present in the column names of meta_data data.frame provided.")
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


#' Calculate and add differences post-cell bender analysis
#'
#' Calculate the difference in features and UMIs per cell when both cell bender and raw assays are present.
#'
#' @param seurat_object object name.
#' @param raw_assay_name name of the assay containing the raw data.
#' @param cell_bender_assay_name name of the assay containing the Cell Bender'ed data.
#'
#' @importFrom magrittr "%>%"
#'
#' @return Seurat object with 2 new columns in the meta.data slot.
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' \dontrun{
#' object <- Add_Cell_Bender_Diff(seurat_object = obj, raw_assay_name = "RAW", cell_bender_assay_name = "RNA")
#' }
#'

Add_Cell_Bender_Diff <- function(
  seurat_object,
  raw_assay_name,
  cell_bender_assay_name
) {
  # Is Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check assays present
  assays_not_found <- Assay_Present(seurat_object = seurat_object, assay_list = c(raw_assay_name, cell_bender_assay_name), print_msg = FALSE, omit_warn = TRUE)[[2]]

  if (!is.null(x = assays_not_found)) {
    stop_quietly()
  }

  # pull meta
  meta_seurat <- seurat_object@meta.data

  # Set variable names
  raw_nFeature <- paste0("nFeature_", raw_assay_name)
  raw_nCount <- paste0("nCount_", raw_assay_name)

  cb_nFeature <- paste0("nFeature_", cell_bender_assay_name)
  cb_nCount <- paste0("nCount_", cell_bender_assay_name)

  # Mutate meta data
  meta_modified <- meta_seurat %>%
    mutate(nFeature_Diff = .data[[raw_nFeature]] - .data[[cb_nFeature]],
           nCount_Diff = .data[[raw_nCount]] - .data[[cb_nCount]])

  # Remove already in object meta data
  meta_modified <- Meta_Remove_Seurat(meta_data = meta_modified, seurat_object = seurat_object)

  # Add back to Seurat Object
  seurat_object <- AddMetaData(object = seurat_object, metadata = meta_modified)

  return(seurat_object)
}


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
#'
#' @return Seurat Object with new entries in the `@misc` slot.
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' \dontrun{
#' obj <- Store_Misc_Info_Seurat(seurat_object = obj_name, data_to_store = data, data_name = "rd1_colors")
#' }
#'

Store_Misc_Info_Seurat <- function(
  seurat_object,
  data_to_store,
  data_name,
  list_as_list = FALSE
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check length of data
  if (length(x = data_to_store) != 1 && class(x = data_to_store) != "list") {
    stop("'data_to_store' must be either single vector/string or a single list of items.")
  }

  # Check class of data
  if (class(x = data_to_store) == "list") {
    if (list_as_list) {
      # Check length of name
      if (length(x = data_name) != 1) {
        stop("The lengths of 'data_to_store' (", length(x = data_to_store), ") and 'data_name' (", length(x = data_name), ") must be equal.")
      }

      # Add data
      seurat_object@misc[[data_name]] <- data_to_store
      message("Seurat Object now contains the following items in @misc slot: ",
              "\n", paste(shQuote(names(x = seurat_object@misc)), collapse=", "))
      return(seurat_object)
    }

    # length of list
    data_list_length <- length(x = data_to_store)

    if (length(x = data_name) != data_list_length) {
      stop("The lengths of 'data_to_store' (", data_list_length, ") and 'data_name' (", length(x = data_name), ") must be equal.")
    }

    # Add data
    for (i in 1:data_list_length) {
      seurat_object@misc[[data_name[i]]] <- data_to_store[[i]]
    }
    message("Seurat Object now contains the following items in @misc slot: ",
            "\n", paste(shQuote(names(x = seurat_object@misc)), collapse=", "))
    return(seurat_object)
  } else {
    # Check length of name
    if (length(x = data_name) != 1) {
      stop("The lengths of 'data_to_store' (", length(x = data_to_store), ") and 'data_name' (", length(x = data_name), ") must be equal.")
    }

    # Add data
    seurat_object@misc[[data_name]] <- data_to_store
    message("Seurat Object now contains the following items in @misc slot: ",
            "\n", paste(shQuote(names(x = seurat_object@misc)), collapse=", "))
    return(seurat_object)
  }
}


#' Store color palette in Seurat object
#'
#' Wrapper function around `Store_Misc_Info_Seurat` to store color palettes.
#'
#' @param seurat_object object name.
#' @param palette(s) vector or list of vectors containing color palettes to store.  If list of palettes
#' see `list_as_list` parameter for control over data storage.
#' @param palette_name name to give the palette(s) in `@misc` slot.  Must be of equal length to the number
#' of data items being stored.
#' @param list_as_list logical.  If `data_to_store` is a list, this dictates whether to store in `@misc` slot
#' as list (TRUE) or whether to store each entry in the list separately (FALSE).  Default is FALSE.
#'
#' @return Seurat Object with new entries in the `@misc` slot.
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' \dontrun{
#' obj <- Store_Palette_Seurat(seurat_object = obj_name, data_to_store = colors_vector, data_name = "rd1_colors")
#' }
#'

Store_Palette_Seurat <- function(
  seurat_object,
  palette,
  palette_name,
  list_as_list = FALSE
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  seurat_object <- Store_Misc_Info_Seurat(seurat_object = seurat_object, data_to_store = palette, data_name = palette_name, list_as_list = list_as_list)
  return(seurat_object)
}


#' Rename Cluster Seurat
#'
#' Wrapper function to rename active identities in Seurat Object with new idents.
#'
#' @param seurat_object object name.
#' @param new_idents vector of new cluster names.  Must be equal to the length of current active.ident
#' in Seurat Object.  Will accept named vector (with old idents as names) or will name the new_idents vector internally.
#' @param meta_col_name (Optional).  Whether or not to create new named column in `Object@meta.data`
#' to store the old identities.
#' @param ... Extra parameters passed to \code{\link[SeuratObject]{RenameIdents}}.
#'
#' @return Seurat Object with new identities placed in active.ident slot.
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' \dontrun{
#' obj <- Rename_Clusters(seurat_object = obj_name, new_idents = new_idents_vec, meta_col_name = "Round01_Res0.6_Idents")
#' }
#'

Rename_Clusters <- function(
  seurat_object,
  new_idents,
  meta_col_name = NULL,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check equivalent lengths
  if (length(x = new_idents) != length(x = levels(x = seurat_object))) {
    stop("Length of `new_idents` must be equal to the number of active.idents in Seurat Object.
         `new_idents` length: ", length(x = new_idents), ", Object@active.idents length: ", length(levels(x = seurat_object)), ".")
  }

  # Name the new idents vector
  if (is.null(x = names(x = new_idents))) {
    names(new_idents) <- levels(seurat_object)
  }
  # If named check that names are right length
  if (!is.null(x = names(x = new_idents)) && length(x = unique(x = names(x = new_idents))) != length(x = levels(x = seurat_object))) {
    stop("The number of unique names for `new idents is not equal to number of active.idents.
         names(new_idents) length: ", length(x = unique(x = names(x = new_idents))), " Object@active.idents length: ", length(levels(x = seurat_object)), ".")
  }

  # Rename meta column for old ident information if desired
  if (!is.null(x = meta_col_name)) {
    seurat_object[[meta_col_name]] <- Idents(seurat_object)
  }

  # Add new idents & return object
  seurat_object <- RenameIdents(object = seurat_object, new_idents)
  return(seurat_object)
}


#' Merge a list of Seurat Objects
#'
#' Enables easy merge of a list of Seurat Objects.  See  See \code{\link[Seurat]{merge}} for more information,
#'
#' @param list_seurat list composed of multiple Seurat Objects.
#' @param add.cell.ids A character vector of length(x = c(x, y)). Appends the corresponding values
#' to the start of each objects' cell names.  See \code{\link[Seurat]{merge}}.
#' @param merge.data Merge the data slots instead of just merging the counts (which requires renormalization).
#' This is recommended if the same normalization approach was applied to all objects.
#' See \code{\link[Seurat]{merge}}.
#' @param project Project name for the Seurat object. See \code{\link[Seurat]{merge}}.
#'
#' @importFrom purrr reduce
#'
#' @return A Seurat Object
#'
#' @export
#'
#' @concept object_util
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
  merged_object <- reduce(list_seurat, function(x, y) {
    merge(x = x, y = y, add.cell.ids = add.cell.ids, merge.data = merge.data, project = project)
  })
}


#' Create a Seurat object containing the data from a liger object
#'
#' Merges raw.data and scale.data of object, and creates Seurat object with these values along with
#' tsne.coords, iNMF factorization, and cluster assignments. Supports Seurat V2 and V3.
#'
#' Stores original dataset identity by default in new object metadata if dataset names are passed
#' in nms. iNMF factorization is stored in dim.reduction object with key "iNMF".
#'
#' @param object \code{liger} object.
#' @param nms By default, labels cell names with dataset of origin (this is to account for cells in
#' different datasets which may have same name). Other names can be passed here as vector, must have
#' same length as the number of datasets. (default names(H)).
#' @param renormalize Whether to log-normalize raw data using Seurat defaults (default TRUE).
#' @param use.liger.genes Whether to carry over variable genes (default TRUE).
#' @param by.dataset Include dataset of origin in cluster identity in Seurat object (default FALSE).
#' @param keep.meta logical. Whether to transfer additional metadata (nGene/nUMI/dataset already transferred)
#' to new Seurat Object.  Default is TRUE.
#' @param reduction_label Name of dimensionality reduction technique used.  Enables accurate transfer
#' or name to Seurat object instead of defaulting to "tSNE".
#' @param seurat_assay Name to set for assay in Seurat Object.  Default is "RNA".
#'
#' @return Seurat object with raw.data, scale.data, dr$reduction_label, dr$inmf, and ident slots set.
#'
#' @references Original function is part of LIGER package (https://github.com/welch-lab/liger) (Licence: GPL-3).
#' Function was slightly modified for use in scCustomize with keep.meta parameter.  Also posted as
#' PR to liger GitHub.
#'
#' @import Matrix
#' @importFrom dplyr pull select
#' @importFrom methods new
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' \dontrun{
#' seurat_object <- Liger_to_Seurat(liger_object = LIGER_OBJ, reduction_label = "UMAP")
#' }

Liger_to_Seurat <- function(
  liger_object,
  nms = names(liger_object@H),
  renormalize = TRUE,
  use.liger.genes = TRUE,
  by.dataset = FALSE,
  keep_meta = TRUE,
  reduction_label = NULL,
  seurat_assay = "RNA"
) {
  if (is.null(x = reduction_label)) {
    stop("`reduction_label` parameter was not set.  LIGER objects do not store name of dimensionality
         reduction technique used.  In order to retain proper labels in Seurat object please set
         `reduction_label` to", '"', "tSNE", '"', ", ", '"', "UMAP", '"', " etc.")
  }

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  # get Seurat version
  maj_version <- packageVersion('Seurat')$major
  if (class(liger_object@raw.data[[1]])[1] != 'dgCMatrix') {
    # mat <- as(x, 'CsparseMatrix')
    liger_object@raw.data <- lapply(liger_object@raw.data, function(x) {
      as(x, 'CsparseMatrix')
    })
  }

  key_name <- paste0(reduction_label, "_")

  raw.data <- Merge_Sparse_Data_All(liger_object@raw.data, nms)
  scale.data <- do.call(rbind, liger_object@scale.data)
  rownames(scale.data) <- colnames(raw.data)
  if (maj_version < 3) {
    var.genes <- liger_object@var.genes
    inmf.obj <- new(
      Class = "dim.reduction", gene.loadings = t(liger_object@W),
      cell.embeddings = liger_object@H.norm, key = "iNMF_"
    )
    rownames(inmf.obj@gene.loadings) <- var.genes

    tsne.obj <- new(
      Class = "dim.reduction", cell.embeddings = liger_object@tsne.coords,
      key = key_name
    )
  } else {
    var.genes <- liger_object@var.genes
    if (any(grepl('_', var.genes))) {
      print("Warning: Seurat v3 genes cannot have underscores, replacing with dashes ('-')")
      var.genes <- gsub("_", replacement = "-", var.genes)
    }
    inmf.loadings <- t(x = liger_object@W)
    inmf.embeddings <- liger_object@H.norm
    ncol_Hnorm <- ncol(x = liger_object@H.norm)
    colnames(inmf.embeddings) <- paste0("iNMF_", 1:ncol_Hnorm)

    tsne.embeddings <- liger_object@tsne.coords
    colnames(tsne.embeddings) <- paste0(key_name, 1:2)
    rownames(x = inmf.loadings) <- var.genes
    rownames(x = inmf.embeddings) <-
      rownames(x = tsne.embeddings) <-
      rownames(x = scale.data)
    inmf.obj <- CreateDimReducObject(
      embeddings = inmf.embeddings,
      loadings = inmf.loadings,
      key = "iNMF_",
      global = TRUE,
      assay = seurat_assay
    )
    tsne.obj <- CreateDimReducObject(
      embeddings = tsne.embeddings,
      key = key_name,
      global = TRUE,
      assay = seurat_assay
    )
  }
  new.seurat <- CreateSeuratObject(raw.data)
  if (renormalize) {
    new.seurat <- NormalizeData(new.seurat)
  }
  if (by.dataset) {
    ident.use <- as.character(unlist(lapply(1:length(liger_object@raw.data), function(i) {
      dataset.name <- names(liger_object@raw.data)[i]
      paste0(dataset.name, as.character(liger_object@clusters[colnames(liger_object@raw.data[[i]])]))
    })))
  } else {
    if (maj_version < 3) {
      ident.use <- as.character(liger_object@clusters)
    } else {
      ident.use <- liger_object@clusters
    }
  }

  if (maj_version < 3) {
    if (use.liger.genes) {
      new.seurat@var.genes <- var.genes
    }
    new.seurat@scale.data <- t(scale.data)
    new.seurat@dr[[reduction_label]] <- tsne.obj
    new.seurat@dr$inmf <- inmf.obj
    new.seurat <- SetIdent(new.seurat, ident.use = ident.use)

  } else {
    if (use.liger.genes) {
      VariableFeatures(new.seurat) <- var.genes
    }
    SetAssayData(new.seurat, slot = "scale.data",  t(scale.data), assay = "RNA")
    new.seurat[[reduction_label]] <- tsne.obj
    new.seurat[['inmf']] <- inmf.obj
    Idents(new.seurat) <- ident.use
  }
  if (keep_meta){
    # extract meta data from liger object
    liger_meta <- liger_object@cell.data
    # remove meta data values already transferred
    liger_meta <- liger_meta %>%
      select(-nUMI, -nGene, -dataset)
    # extract meta data names
    meta_names <- colnames(liger_meta)
    # add meta data to new seurat object
    for (meta_var in meta_names){
      meta_transfer <- liger_meta %>%
        pull(meta_var)
      names(meta_transfer) <- colnames(x = new.seurat)
      new.seurat <- AddMetaData(object = new.seurat,
                                metadata = meta_transfer,
                                col.name = meta_var)
    }
  }

  return(new.seurat)
}

