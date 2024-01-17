#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### CONVERT TO LIGER ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Create liger object from one Seurat Object
#'
#' This function is part of generic `as.LIGER` for conversion of Seurat Objects.  Updated for Seurat5 compatibility and enhanced functionality.
#'
#' @param group.by Variable in meta data which contains variable to split data by, (default is "orig.ident").
#' @param assay Assay containing raw data to use, (default is "RNA").
#' @param remove_missing logical, whether to remove missing genes with no counts when converting to
#' LIGER object (default is FALSE).
#' @param renormalize logical, whether to perform normalization after LIGER object creation (default is TRUE).
#' @param use_seurat_var_genes logical, whether to transfer variable features from Seurat object to
#' new LIGER object (default is FALSE).
#' @param use_seurat_dimreduc logical, whether to transfer dimensionality reduction coordinates from
#' Seurat to new LIGER object (default is FALSE).
#' @param reduction Name of Seurat reduction to transfer if `use_seurat_dimreduc = TRUE`.
#' @param keep_meta logical, whether to transfer columns in Seurat meta.data slot to LIGER cell.data
#' slot (default is TRUE).
#' @param verbose logical, whether to print status messages during object conversion (default is TRUE).
#'
#'
#' @references modified and enhanced version of `rliger::seuratToLiger`.
#'
#' @method as.LIGER Seurat
#' @return liger object.
#'
#' @concept object_conversion
#'
#' @import cli
#' @import Seurat
#' @importFrom dplyr left_join join_by select any_of
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @export
#' @rdname as.LIGER
#'
#' @examples
#' \dontrun{
#' liger_object <- as.LIGER(x = seurat_object)
#' }
#'

as.LIGER.Seurat <- function(
    x,
    group.by = "orig.ident",
    assay = "RNA",
    remove_missing = FALSE,
    renormalize = TRUE,
    use_seurat_var_genes = FALSE,
    use_seurat_dimreduc = FALSE,
    reduction = NULL,
    keep_meta = TRUE,
    verbose = TRUE,
    ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = x)

  # Run update to ensure functionality
  if (isTRUE(x = verbose)) {
    cli_inform(message = "Checking Seurat object validity")
  }

  x <- suppressMessages(UpdateSeuratObject(object = x))

  # Check Assay5 for multiple layers
  if (isTRUE(x = Assay5_Check(seurat_object = x, assay = assay))) {
    layers_check <- Layers(object = x, search = "counts")
    if (length(x = layers_check) > 1) {
      cli_abort(message = c("Multiple layers containing raw counts present {.field {head(x = layers_check, n = 2)}}.",
                            "i" = "Please run {.code JoinLayers} before converting to LIGER object."))
    }
  }

  # Check meta data
  group.by <- Meta_Present(seurat_object = x, meta_col_names = group.by, omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)[[1]]

  # stop if none found
  if (length(x = group.by) == 0) {
    cli_abort(message = c("{.code group.by} was not found.",
                          "i" = "No column found in object meta.data named: {.val {group.by}}.")
    )
  }

  # Set ident to grouping variable
  Idents(x) <- group.by

  # Check & Set Assay
  if (!assay %in% Assays(object = x)) {
    cli_abort(message = "Provided assay {.field {assay}} not found in Seurat object.")
  }

  if (assay != DefaultAssay(object = x)) {
    cli_inform("Changing object DefaultAssay from ({.field {DefaultAssay(object = x)}}) to provided assay ({.field {assay}}).")
    DefaultAssay(x) <- assay
  }

  # Check & Pull other relevant data
  if (isTRUE(x = use_seurat_dimreduc)) {
    # Extract default reduction
    reduction <- reduction %||% DefaultDimReduc(object = x)

    if (!reduction %in% Reductions(object = x)) {
      cli_abort(message = "Provided reduction: {.field {reduction}} was not found in Seurat Object.")
    }

    reduc_coords <- Embeddings(object = x, reduction = reduction)
  }

  if (isTRUE(x = use_seurat_var_genes)) {
    var_genes <- VariableFeatures(object = x)

    if (!length(x = var_genes) > 0) {
      cli_abort(message ="{.code use_seurat_var_genes = TRUE}, but no variable features found in Seurat object.")
    }
  }

  # Get raw data & cells
  raw_data_full <- LayerData(object = x, layer = layers_check)

  cells_per_dataset <- CellsByIdentities(object = x)

  # Split data by dataset
  idents <- names(x = cells_per_dataset)

  raw_data_list <- lapply(idents, function(x){
    raw_data_full[, cells_per_dataset[[x]]]
  })

  names(raw_data_list) <- idents

  # Create LIGER Object
  if (isTRUE(x = verbose)) {
    cli_inform(message = "Creating LIGER object.")
  }

  liger_object <- rliger::createLiger(raw.data = raw_data_list, remove.missing = remove_missing)

  if (isTRUE(x = renormalize)) {
    if (isTRUE(x = verbose)) {
      cli_inform(message = "Normalizing data.")
    }
    liger_object <- rliger::normalize(object = liger_object, remove.missing = remove_missing)
  }

  # Add var genes
  if (isTRUE(x = use_seurat_var_genes)) {
    liger_object@var.genes <- var_genes
  }

  # Add dim reduc
  if (isTRUE(x = use_seurat_dimreduc)) {
    liger_object@tsne.coords <- reduc_coords

    # Add new attribute to enable more accurate scCustomize plotting
    attributes(liger_object)$reduction_key <- reduction
  }

  # transfer meta
  if (isTRUE(x = keep_meta)) {
    # extract meta data from liger object
    seurat_meta <- Fetch_Meta(object = x)
    # remove meta data values already transferred
    seurat_meta <- seurat_meta %>%
      select(-any_of(c("nFeature_RNA", "nCount_RNA", group.by))) %>%
      rownames_to_column("barcodes")

    # pull current liger meta
    liger_meta <- Fetch_Meta(object = liger_object) %>%
      rownames_to_column("barcodes")

    # join meta
    new_liger_meta <- suppressMessages(left_join(x = liger_meta, y = seurat_meta, by = join_by("barcodes"))) %>%
      column_to_rownames("barcodes")

    # Add to LIGER object
    liger_object@cell.data <- new_liger_meta
  }

  # return object
  return(liger_object)
}


#' Create liger object from one Seurat Object
#'
#' This function is part of generic `as.LIGER` for conversion of Seurat Objects.  Updated for Seurat5 compatibility and enhanced functionality.
#'
#' @param group.by Variable in meta data which contains variable to split data by, (default is "orig.ident").
#' @param dataset_names optional, vector of names to use for naming datasets.
#' @param assay Assay containing raw data to use, (default is "RNA").
#' @param remove_missing logical, whether to remove missing genes with no counts when converting to
#' LIGER object (default is FALSE).
#' @param renormalize logical, whether to perform normalization after LIGER object creation (default is TRUE).
#' @param use_seurat_var_genes logical, whether to transfer variable features from Seurat object to
#' new LIGER object (default is FALSE).
#' @param var_genes_method how variable genes should be selected from Seurat objects if `use_seurat_var_genes = TRUE`.  Can be either "intersect" or "union", (default is "intersect").
#' @param keep_meta logical, whether to transfer columns in Seurat meta.data slot to LIGER cell.data
#' slot (default is TRUE).
#' @param verbose logical, whether to print status messages during object conversion (default is TRUE).
#'
#'
#' @references modified and enhanced version of `rliger::seuratToLiger`.
#'
#' @method as.LIGER list
#' @return liger object.
#'
#' @concept object_conversion
#'
#' @import cli
#' @import Seurat
#' @importFrom dplyr left_join join_by select any_of bind_rows union intersect
#' @importFrom stringr str_to_lower
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @export
#' @rdname as.LIGER
#'
#' @examples
#' \dontrun{
#' liger_object <- as.LIGER(x = seurat_object_list)
#' }
#'

as.LIGER.list <- function(
    x,
    group.by = "orig.ident",
    dataset_names = NULL,
    assay = "RNA",
    remove_missing = FALSE,
    renormalize = TRUE,
    use_seurat_var_genes = FALSE,
    var_genes_method = "intersect",
    keep_meta = TRUE,
    verbose = TRUE,
    ...
) {
  # Check Seurat
  seurat_check <- unlist(lapply(x, function(x) {
    inherits(x = x, what = "Seurat")
  }))

  if (any(seurat_check) == "FALSE") {
    cli_abort(message = "One or more of items in list are not Seurat Objects.")
  }

  # Run update to ensure functionality
  if (isTRUE(x = verbose)) {
    cli_inform(message = "Checking Seurat object validity")
  }

  x <- lapply(x, function(y) {
    suppressMessages(UpdateSeuratObject(object = y))
  })

  # Check Assay5 for multiple layers
  for (i in x) {
    if (isTRUE(x = Assay5_Check(seurat_object = i, assay = assay))) {
      layers_check <- Layers(object = i, search = "counts")
      if (length(x = layers_check) > 1) {
        cli_abort(message = c("Multiple layers containing raw counts present {.field {head(x = layers_check, n = 2)}}.",
                              "i" = "Please run {.code JoinLayers} before converting to LIGER object."))
      }
    }
  }

  # Check meta data
  if (is.null(x = dataset_names)) {
    for (j in x) {
      group.by <- Meta_Present(seurat_object = j, meta_col_names = group.by, omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)[[1]]

      # stop if none found
      if (length(x = group.by) == 0) {
        cli_abort(message = c("{.code group.by} was not found in all objects in list.",
                              "i" = "All objects must contain column in meta.data named: {.val {group.by}}.")
        )
      }
    }

  } else {
    if (length(x = dataset_names) != length(x = x)) {
      cli_abort(message = "The number of {.code dataset_names} provided ({.field {length(x = dataset_names)}}) does not match number of Seurat objects in list ({.field {length(x = x)}}).")
    }
  }

  # Check & Set Assay
  for (k in x) {
    if (!assay %in% Assays(object = k)) {
      cli_abort(message = "Provided assay {.field {assay}} not found in all Seurat objects in list.")
    }
  }

  for (l in x) {
    if (assay != DefaultAssay(object = l)) {
      cli_inform("Changing object DefaultAssay from ({.field {DefaultAssay(object = x)}}) to provided assay ({.field {assay}}).")
      DefaultAssay(l) <- assay
    }
  }

  if (isTRUE(x = use_seurat_var_genes)) {
    var_genes <- lapply(x, function(z) {
      VariableFeatures(object = z)
    })

    for (m in var_genes) {
      if (!length(x = m) > 0) {
        cli_abort(message ="{.code use_seurat_var_genes = TRUE}, but not all objects in list have variable features.")
      }
    }

    var_genes_method <- str_to_lower(string = var_genes_method)
    if (!var_genes_method %in% c("intersect", "union")) {
      cli_abort(message = "{.code var_genes_method} must be either {.field intersect} or {.field union}.")
    }

    if (var_genes_method == "union") {
      var_genes <- reduce(var_genes, function(a, b) {
        union(x = a, y = b)})
    }
    if (var_genes_method == "intersect") {
      var_genes <- reduce(var_genes, function(c, d) {
        intersect(x = c, y = d)
      })
    }
  }

  # Get raw data & cells
  raw_data_list <- lapply(x, function(e){
    LayerData(object = e, layer = "counts")
  })

  if (is.null(x = dataset_names)) {
    group_names <- unique(x = sapply(1:length(x = x), function(f) {
      obj_meta <- Fetch_Meta(object = x[[f]]) %>%
        select(any_of(group.by)) %>%
        unique()
      if (length(x = obj_meta) > 1) {
        cli_abort(message = c("Some objects in list have multiple values within the {.field {group.by}} column.",
                              "i" = "This column must only contain one value per object"))
      }
    }))

    if (length(x = group_names) != length(x = x)) {
      cli_abort(message = c("Some objects in list have the same values within the {.field {group.by}} column.",
                            "i" = "All objects must have unique value in this column."))
    }

    names(x = raw_data_list) <- group_names
  } else {
    names(x = raw_data_list) <- dataset_names
  }


  # Create LIGER Object
  if (isTRUE(x = verbose)) {
    cli_inform(message = "Creating LIGER object.")
  }

  liger_object <- rliger::createLiger(raw.data = raw_data_list, remove.missing = remove_missing)

  if (isTRUE(x = renormalize)) {
    if (isTRUE(x = verbose)) {
      cli_inform(message = "Normalizing data.")
    }
    liger_object <- rliger::normalize(object = liger_object, remove.missing = remove_missing)
  }

  # Add var genes
  if (isTRUE(x = use_seurat_var_genes)) {
    liger_object@var.genes <- var_genes
  }

  # transfer meta
  if (isTRUE(x = keep_meta)) {
    # extract meta data from seurat object
    seurat_meta <- lapply(x, function(g) {
      obj_meta <- Fetch_Meta(object = g) %>%
        select(-any_of(c("nFeature_RNA", "nCount_RNA", group.by)))
    })

    seurat_meta <- bind_rows(seurat_meta) %>%
      rownames_to_column("barcodes")

    # pull current liger meta
    liger_meta <- Fetch_Meta(object = liger_object) %>%
      rownames_to_column("barcodes")

    # join meta
    new_liger_meta <- suppressMessages(left_join(x = liger_meta, y = seurat_meta, by = join_by("barcodes"))) %>%
      column_to_rownames("barcodes")

    # Add to LIGER object
    liger_object@cell.data <- new_liger_meta
  }

  # return object
  return(liger_object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### CONVERT TO SEURAT ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Create a Seurat object containing the data from a liger object
#'
#' Merges raw.data and scale.data of object, and creates Seurat object with these values along with
#' tsne.coords, iNMF factorization, and cluster assignments. Supports Seurat V2 and V3.
#'
#' Stores original dataset identity by default in new object metadata if dataset names are passed
#' in nms. iNMF factorization is stored in dim.reduction object with key "iNMF".
#'
#' @param liger_object \code{liger} object.
#' @param nms By default, labels cell names with dataset of origin (this is to account for cells in
#' different datasets which may have same name). Other names can be passed here as vector, must have
#' same length as the number of datasets. (default names(H)).
#' @param renormalize Whether to log-normalize raw data using Seurat defaults (default TRUE).
#' @param use.liger.genes Whether to carry over variable genes (default TRUE).
#' @param by.dataset Include dataset of origin in cluster identity in Seurat object (default FALSE).
#' @param keep_meta logical. Whether to transfer additional metadata (nGene/nUMI/dataset already transferred)
#' to new Seurat Object.  Default is TRUE.
#' @param reduction_label Name of dimensionality reduction technique used.  Enables accurate transfer
#' or name to Seurat object instead of defaulting to "tSNE".
#' @param seurat_assay Name to set for assay in Seurat Object.  Default is "RNA".
#' @param add_barcode_names logical, whether to add dataset names to the cell barcodes when
#' creating Seurat object, default is FALSE.
#' @param barcode_prefix logical, if `add_barcode_names = TRUE` should the names be added as
#' prefix to current cell barcodes/names or a suffix (default is TRUE; prefix).
#' @param barcode_cell_id_delimiter The delimiter to use when adding dataset id to barcode
#' prefix/suffix.  Default is "_".
#'
#' @return Seurat object with raw.data, scale.data, reduction_label, iNMF, and ident slots set.
#'
#' @references Original function is part of LIGER package \url{https://github.com/welch-lab/liger} (Licence: GPL-3).
#' Function was slightly modified for use in scCustomize with keep.meta parameter.  Also posted as
#' PR to liger GitHub.
#'
#' @import cli
#' @import Matrix
#' @importFrom dplyr any_of pull select
#' @importFrom methods as new
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
    reduction_label = "UMAP",
    seurat_assay = "RNA",
    add_barcode_names = FALSE,
    barcode_prefix = TRUE,
    barcode_cell_id_delimiter = "_"
) {
  if (is.null(x = reduction_label)) {
    cli_abort(message = c("{.code reduction_label} parameter was not set.",
                          "*" = "LIGER objects do not store name of dimensionality reduction technique used.",
                          "i" = "In order to retain proper labels in Seurat object please set {.code reduction_label} to {.val tSNE}, {.val UMAP}, {.val etc}."))
  }

  # Adjust name for dimreduc key
  key_name <- paste0(reduction_label, "_")

  # adjust raw data slot if needed
  if (!inherits(x = liger_object@raw.data[[1]], what = 'dgCMatrix')) {
    liger_object@raw.data <- lapply(liger_object@raw.data, as, Class = "CsparseMatrix")
  }

  # merge raw data
  if (isTRUE(x = add_barcode_names)) {
    raw.data <- Merge_Sparse_Data_All(matrix_list = liger_object@raw.data, add_cell_ids = nms, prefix = barcode_prefix, cell_id_delimiter = barcode_cell_id_delimiter)
  } else {
    raw.data <- Merge_Sparse_Data_All(matrix_list = liger_object@raw.data)
  }

  # create object
  new.seurat <- CreateSeuratObject(counts = raw.data, assay = seurat_assay)

  # normalize data
  if (isTRUE(x = renormalize)) {
    new.seurat <- Seurat::NormalizeData(new.seurat)
  } else {
    if (length(x = liger_object@norm.data) > 0) {
      if (isTRUE(x = add_barcode_names)) {
        norm.data <- Merge_Sparse_Data_All(matrix_list = liger_object@norm.data, add_cell_ids = nms, prefix = barcode_prefix, cell_id_delimiter = barcode_cell_id_delimiter)
      } else {
        norm.data <- Merge_Sparse_Data_All(matrix_list = liger_object@norm.data)
      }

      new.seurat <- SetAssayData(object = new.seurat, layer = "data", slot = "data", new.data = norm.data)
    }
  }

  if (length(x = liger_object@var.genes) > 0 && isTRUE(x = use.liger.genes)) {
    VariableFeatures(object = new.seurat) <- liger_object@var.genes
  }
  if (length(x = liger_object@scale.data) > 0) {
    scale.data <- t(x = Reduce(rbind, liger_object@scale.data))
    colnames(x = scale.data) <- colnames(x = raw.data)
    new.seurat <- SetAssayData(object = new.seurat, layer = "scale.data", slot = "scale.data", new.data = scale.data)
  }


  if (all(dim(x = liger_object@W) > 0) && all(dim(x = liger_object@H.norm) > 0)) {
    inmf.loadings <- t(x = liger_object@W)
    rinmf.loadings <- inmf.loadings

    dimnames(x = inmf.loadings) <- list(liger_object@var.genes,
                                        paste0("iNMF_", seq_len(ncol(inmf.loadings))))
    dimnames(x = rinmf.loadings) <- list(liger_object@var.genes,
                                         paste0("rawiNMF_", seq_len(ncol(rinmf.loadings))))

    inmf.embeddings <- liger_object@H.norm
    rinmf.embeddings <- do.call(what = 'rbind', args = liger_object@H)

    dimnames(x = inmf.embeddings) <- list(unlist(x = lapply(liger_object@scale.data, rownames), use.names = FALSE),
                                          paste0("iNMF_", seq_len(ncol(inmf.loadings))))
    dimnames(x = rinmf.embeddings) <- list(unlist(x = lapply(liger_object@scale.data, rownames), use.names = FALSE),
                                           paste0("rawiNMF_", seq_len(ncol(x = inmf.loadings))))


    inmf.obj <- CreateDimReducObject(
      embeddings = inmf.embeddings,
      loadings = inmf.embeddings,
      assay = seurat_assay,
      global = TRUE,
      key = "iNMF_"
    )
    new.seurat[["iNMF"]] <- inmf.obj

    rinmf.obj <- CreateDimReducObject(
      embeddings = rinmf.embeddings,
      loadings = rinmf.loadings,
      key = "rawiNMF_",
      global = TRUE,
      assay = seurat_assay
    )
  }


  if (all(dim(x = liger_object@tsne.coords) > 0)) {
    dimreduc.embeddings <- liger_object@tsne.coords
    dimnames(x = dimreduc.embeddings) <- list(rownames(liger_object@H.norm),
                                              paste0(key_name, 1:2))

    dimreduc.obj <- CreateDimReducObject(
      embeddings = dimreduc.embeddings,
      assay = seurat_assay,
      global = TRUE,
      key = key_name
    )
    new.seurat[[reduction_label]] <- dimreduc.obj
  }

  new.seurat$orig.ident <- liger_object@cell.data$dataset

  idents <- liger_object@clusters

  if (length(x = idents) == 0 || isTRUE(x = by.dataset)) idents <- liger_object@cell.data$dataset
  Idents(object = new.seurat) <- idents

  # transfer meta
  if (isTRUE(x = keep_meta)) {
    # extract meta data from liger object
    liger_meta <- Fetch_Meta(object = liger_object)
    # remove meta data values already transferred
    liger_meta <- liger_meta %>%
      select(-any_of(c("nUMI", "nGene", "dataset")))
    # extract meta data names
    meta_names <- colnames(x = liger_meta)
    # add meta data to new seurat object
    for (meta_var in meta_names){
      meta_transfer <- liger_meta %>%
        pull(meta_var)
      names(x = meta_transfer) <- colnames(x = new.seurat)
      new.seurat <- AddMetaData(object = new.seurat,
                                metadata = meta_transfer,
                                col.name = meta_var)
    }
  }

  # return object
  return(new.seurat)
}
