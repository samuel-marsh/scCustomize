#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### CONVERT TO LIGER ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Create liger object from one Seurat Object
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
#'
#' @concept object_conversion
#'
#' @import cli
#' @import Seurat
#' @importFrom dplyr left_join join_by select any_of
#' @importFrom magrittr "%>%"
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
    cli_inform(message = c("*" = "Checking Seurat object validity"))
  }

  x <- suppressMessages(UpdateSeuratObject(object = x))

  # Check & Set Assay
  if (!assay %in% Assays(object = x)) {
    cli_abort(message = "Provided assay {.field {assay}} not found in Seurat object.")
  }

  if (assay != DefaultAssay(object = x)) {
    cli_inform(c("*" = "Changing object DefaultAssay from ({.field {DefaultAssay(object = x)}}) to provided assay ({.field {assay}})."))
    DefaultAssay(x) <- assay
  }

  # Check Assay5 for multiple layers
  count_layers <- Layers(object = x, search = "counts")

  if (isTRUE(x = Assay5_Check(seurat_object = x, assay = assay))) {
    if (length(x = count_layers) > 1 && group.by != "orig.ident") {
      cli_abort(message = c("Multiple layers containing raw counts present ({.field {count_layers[1]}}, {.field {count_layers[2]}}, {.field ...}) and a non-default value was provided to {.code group.by}.",
                            "i" = "To group LIGER object by assay layers do not change value of {.code group.by}",
                            "i" = "If {.code group.by} variable is different than the split layers please run {.code JoinLayers} first."))
    }
  }

  # Check meta data
  group.by <- Meta_Present(object = x, meta_col_names = group.by, omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)[[1]]

  # stop if none found
  if (length(x = group.by) == 0) {
    cli_abort(message = c("{.code group.by} was not found.",
                          "i" = "No column found in object meta.data named: {.val {group.by}}.")
    )
  }

  # Set ident to grouping variable
  if (length(x = count_layers) == 1) {
    Idents(object = x) <- group.by
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
  if (length(x = count_layers) == 1) {
    raw_data_full <- LayerData(object = x, layer = count_layers)

    cells_per_dataset <- CellsByIdentities(object = x)

    # Split data by dataset
    idents <- names(x = cells_per_dataset)

    raw_data_list <- lapply(idents, function(x){
      raw_data_full[, cells_per_dataset[[x]]]
    })

    names(raw_data_list) <- idents
  }

  # If multiple layers
  if (length(x = count_layers) > 1) {
    raw_data_list <- lapply(count_layers, function (i){
      counts <- LayerData(object = x, layer = i)
    })

    new_names <- gsub(pattern = "counts.", replacement = "", x = count_layers)

    names(raw_data_list) <- new_names
  }

  # Create LIGER Object
  if (isTRUE(x = verbose)) {
    cli_inform(message = c("*" = "Creating LIGER object."))
  }

  liger_object <- rliger::createLiger(raw.data = raw_data_list, remove.missing = remove_missing)

  if (isTRUE(x = renormalize)) {
    if (isTRUE(x = verbose)) {
      cli_inform(message = c("*" = "Normalizing data."))
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
#' @method as.LIGER list
#'
#' @concept object_conversion
#'
#' @import cli
#' @import Seurat
#' @importFrom dplyr left_join join_by select any_of bind_rows union intersect
#' @importFrom magrittr "%>%"
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
    cli_inform(message = c("*" = "Checking Seurat object validity"))
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
      group.by <- Meta_Present(object = j, meta_col_names = group.by, omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)[[1]]

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
      cli_inform(c("*" = "Changing object DefaultAssay from ({.field {DefaultAssay(object = x)}}) to provided assay ({.field {assay}})."))
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
    counts_layer <- Layers(object = e, search = "counts")
    LayerData(object = e, layer = counts_layer)
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
    cli_inform(message = c("*" = "Creating LIGER object."))
  }

  liger_object <- rliger::createLiger(raw.data = raw_data_list, remove.missing = remove_missing)

  if (isTRUE(x = renormalize)) {
    if (isTRUE(x = verbose)) {
      cli_inform(message = c("*" = "Normalizing data."))
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


#' Convert objects to \code{Seurat} objects
#'
#' Merges raw.data and scale.data of object, and creates Seurat object with these values along with slots
#' containing dimensionality reduction coordinates, iNMF factorization, and cluster assignments.
#' Supports Seurat V3/4 and V4.
#'
#' Stores original dataset identity by default in new object metadata if dataset names are passed
#' in nms. iNMF factorization is stored in dim.reduction object with key "iNMF".
#'
#' @param x \code{liger} object.
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
#' @param assay_type what type of Seurat assay to create in new object (Assay vs Assay5).
#' Default is NULL which will default to the current user settings.
#' See \code{\link{Convert_Assay}} parameter `convert_to` for acceptable values.
#' @param add_barcode_names logical, whether to add dataset names to the cell barcodes when
#' creating Seurat object, default is FALSE.
#' @param barcode_prefix logical, if `add_barcode_names = TRUE` should the names be added as
#' prefix to current cell barcodes/names or a suffix (default is TRUE; prefix).
#' @param barcode_cell_id_delimiter The delimiter to use when adding dataset id to barcode
#' prefix/suffix.  Default is "_".
#' @param ... unused.
#'
#' @return Seurat object with raw.data, scale.data, reduction_label, iNMF, and ident slots set.
#'
#' @references Original function is part of LIGER package \url{https://github.com/welch-lab/liger} (Licence: GPL-3).
#' Function was modified for use in scCustomize with additional parameters/functionality.
#'
#' @method as.Seurat liger
#' @return Seurat object.
#'
#' @concept object_conversion
#'
#' @import cli
#' @import Matrix
#' @import Seurat
#' @importFrom dplyr any_of pull select
#' @importFrom magrittr "%>%"
#' @importFrom methods as new
#' @importFrom utils packageVersion
#'
#' @export
#' @rdname as.Seurat
#'
#' @examples
#' \dontrun{
#' seurat_object <- as.Seurat(x = liger_object)
#' }
#'

as.Seurat.liger <- function(
    x,
    nms = names(liger_object@H),
    renormalize = TRUE,
    use.liger.genes = TRUE,
    by.dataset = FALSE,
    keep_meta = TRUE,
    reduction_label = "UMAP",
    seurat_assay = "RNA",
    assay_type = NULL,
    add_barcode_names = FALSE,
    barcode_prefix = TRUE,
    barcode_cell_id_delimiter = "_",
    ...
) {
  if (is.null(x = reduction_label)) {
    cli_abort(message = c("{.code reduction_label} parameter was not set.",
                          "*" = "LIGER objects do not store name of dimensionality reduction technique used.",
                          "i" = "In order to retain proper labels in Seurat object please set {.code reduction_label} to {.val tSNE}, {.val UMAP}, {.val etc}."))
  }

  # Adjust name for dimreduc key
  key_name <- paste0(reduction_label, "_")

  # adjust raw data slot if needed
  if (!inherits(x = x@raw.data[[1]], what = 'dgCMatrix')) {
    x@raw.data <- lapply(x@raw.data, as, Class = "CsparseMatrix")
  }

  # check assay_type is ok
  if (!is.null(x = assay_type)) {
    # Check accepted
    accepted_V3 <- c("Assay", "assay", "V3", "v3")
    accepted_V5 <- c("Assay5", "assay5", "V5", "v5")

    if (!convert_to %in% c(accepted_V5, accepted_V3)) {
      cli_abort(message = c("Value provided to {.code convert_to} ({.field {convert_to}}) was not accepted value.",
                            "i" = "Accepted values to convert to V3/4 are: {.field {accepted_V3}}",
                            "i" = "Accepted values to convert to V5 are: {.field {accepted_V5}}"))
    }

    # set assay value
    if (convert_to %in% accepted_V5) {
      if (packageVersion(pkg = 'Seurat') < 5) {
        cli_abort(message = "Seurat version must be v5.0.0 or greater to create To create Assay5.")
      }

      convert_to <- "v5"
    }

    if (convert_to %in% accepted_V3) {
      convert_to <- "v3"
    }
  }

  # merge raw data
  if (isTRUE(x = add_barcode_names)) {
    raw.data <- Merge_Sparse_Data_All(matrix_list = x@raw.data, add_cell_ids = nms, prefix = barcode_prefix, cell_id_delimiter = barcode_cell_id_delimiter)
  } else {
    raw.data <- Merge_Sparse_Data_All(matrix_list = x@raw.data)
  }

  # create object
  new.seurat <- CreateSeuratObject(counts = raw.data, assay = seurat_assay)

  # normalize data
  if (isTRUE(x = renormalize)) {
    new.seurat <- Seurat::NormalizeData(new.seurat)
  } else {
    if (length(x = x@norm.data) > 0) {
      if (isTRUE(x = add_barcode_names)) {
        norm.data <- Merge_Sparse_Data_All(matrix_list = x@norm.data, add_cell_ids = nms, prefix = barcode_prefix, cell_id_delimiter = barcode_cell_id_delimiter)
      } else {
        norm.data <- Merge_Sparse_Data_All(matrix_list = x@norm.data)
      }

      new.seurat <- SetAssayData(object = new.seurat, layer = "data", slot = "data", new.data = norm.data)
    }
  }

  if (length(x = x@var.genes) > 0 && isTRUE(x = use.liger.genes)) {
    VariableFeatures(object = new.seurat) <- x@var.genes
  }
  if (length(x = x@scale.data) > 0) {
    scale.data <- t(x = Reduce(rbind, x@scale.data))
    colnames(x = scale.data) <- colnames(x = raw.data)
    new.seurat <- SetAssayData(object = new.seurat, layer = "scale.data", slot = "scale.data", new.data = scale.data)
  }


  if (all(dim(x = x@W) > 0) && all(dim(x = x@H.norm) > 0)) {
    inmf.loadings <- t(x = x@W)
    rinmf.loadings <- inmf.loadings

    dimnames(x = inmf.loadings) <- list(x@var.genes,
                                        paste0("iNMF_", seq_len(ncol(inmf.loadings))))
    dimnames(x = rinmf.loadings) <- list(x@var.genes,
                                         paste0("rawiNMF_", seq_len(ncol(rinmf.loadings))))

    inmf.embeddings <- x@H.norm
    rinmf.embeddings <- do.call(what = 'rbind', args = x@H)

    dimnames(x = inmf.embeddings) <- list(unlist(x = lapply(x@scale.data, rownames), use.names = FALSE),
                                          paste0("iNMF_", seq_len(ncol(inmf.loadings))))
    dimnames(x = rinmf.embeddings) <- list(unlist(x = lapply(x@scale.data, rownames), use.names = FALSE),
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


  if (all(dim(x = x@tsne.coords) > 0)) {
    dimreduc.embeddings <- x@tsne.coords
    dimnames(x = dimreduc.embeddings) <- list(rownames(x@H.norm),
                                              paste0(key_name, 1:2))

    dimreduc.obj <- CreateDimReducObject(
      embeddings = dimreduc.embeddings,
      assay = seurat_assay,
      global = TRUE,
      key = key_name
    )
    new.seurat[[reduction_label]] <- dimreduc.obj
  }

  new.seurat$orig.ident <- x@cell.data$dataset

  idents <- x@clusters

  if (length(x = idents) == 0 || isTRUE(x = by.dataset)) idents <- x@cell.data$dataset
  Idents(object = new.seurat) <- idents

  # transfer meta
  if (isTRUE(x = keep_meta)) {
    # extract meta data from liger object
    liger_meta <- Fetch_Meta(object = x)
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


  if (!is.null(x = assay_type)) {
    options_list <- options()
    if (options_list$Seurat.object.assay.version != convert_to) {
      new.seurat <- Convert_Assay(seurat_object = new.seurat, convert_to = convert_to)
    }
  }

  # return object
  return(new.seurat)
}


#' Create a Seurat object containing the data from a liger object `r lifecycle::badge("soft-deprecated")`
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
#' @param assay_type what type of Seurat assay to create in new object (Assay vs Assay5).
#' Default is NULL which will default to the current user settings.
#' See \code{\link{Convert_Assay}} parameter `convert_to` for acceptable values.
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
#' @importFrom magrittr "%>%"
#' @importFrom methods as new
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @concept object_conversion
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
    assay_type = NULL,
    add_barcode_names = FALSE,
    barcode_prefix = TRUE,
    barcode_cell_id_delimiter = "_"
) {
  lifecycle::deprecate_soft(when = "2.1.0",
                            what = "Liger_to_Seurat()",
                            with = "as.Seurat()",
                            details = c("i" = "Please adjust code now to prepare for full deprecation.")
  )



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

  # check assay_type is ok
  if (!is.null(x = assay_type)) {
    # Check accepted
    accepted_V3 <- c("Assay", "assay", "V3", "v3")
    accepted_V5 <- c("Assay5", "assay5", "V5", "v5")

    if (!convert_to %in% c(accepted_V5, accepted_V3)) {
      cli_abort(message = c("Value provided to {.code convert_to} ({.field {convert_to}}) was not accepted value.",
                            "i" = "Accepted values to convert to V3/4 are: {.field {accepted_V3}}",
                            "i" = "Accepted values to convert to V5 are: {.field {accepted_V5}}"))
    }

    # set assay value
    if (convert_to %in% accepted_V5) {
      if (packageVersion(pkg = 'Seurat') < 5) {
        cli_abort(message = "Seurat version must be v5.0.0 or greater to create To create Assay5.")
      }

      convert_to <- "v5"
    }

    if (convert_to %in% accepted_V3) {
      convert_to <- "v3"
    }
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


  if (!is.null(x = assay_type)) {
    options_list <- options()
    if (options_list$Seurat.object.assay.version != convert_to) {
      new.seurat <- Convert_Assay(seurat_object = new.seurat, convert_to = convert_to)
    }
  }

  # return object
  return(new.seurat)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### CONVERT TO ANNDATA ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Create & Save Anndata Object
#'
#' @param file_path directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name file name.
#' @param assay Assay containing data to use, (default is "RNA").
#' @param main_layer the layer of data to become default layer in anndata object (default is "data").
#' @param other_layers other data layers to transfer to anndata object (default is "counts").
#' @param transer_dimreduc logical, whether to transfer dimensionality reduction coordinates from
#' Seurat to anndata object (default is TRUE).
#' @param verbose logical, whether to print status messages during object conversion (default is TRUE).
#'
#'
#' @references Seurat version modified and enhanced version of `sceasy::seurat2anndata` (sceasy package: \url{https://github.com/cellgeni/sceasy}; License: GPL-3.  Function has additional checks and supports Seurat V3 and V5 object structure.
#'
#' @method as.anndata Seurat
#'
#' @concept object_conversion
#'
#' @import cli
#' @import Seurat
#' @importFrom stringr str_to_lower
#'
#' @export
#' @rdname as.anndata
#'
#' @examples
#' \dontrun{
#' as.anndata(x = seurat_object, file_path = "/folder_name", file_name = "anndata_converted.h5ad")
#' }
#'

as.anndata.Seurat <- function(
    x,
    file_path,
    file_name,
    assay = "RNA",
    main_layer = "data",
    other_layers = "counts",
    transer_dimreduc = TRUE,
    verbose = TRUE,
    ...
) {
  # Check reticulate installed
  reticulate_check <- is_installed(pkg = "reticulate")
  if (isFALSE(x = reticulate_check)) {
    cli_abort(message = c(
      "Please install the {.val reticulate} package to use {.code as.anndata}.",
      "i" = "This can be accomplished with the following commands: ",
      "----------------------------------------",
      "{.field `install.packages({symbol$dquote_left}reticulate{symbol$dquote_right})`}",
      "----------------------------------------"
    ))
  }

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path)) {
    if (!dir.exists(paths = file_path)) {
      cli_abort(message = "Provided {.code file_path}: {symbol$dquote_left}{.field {file_path}}{symbol$dquote_right} does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name)) {
    cli_abort(message = "No file name provided.  Please provide a file name using {.code file_name}.")
  }

  file_ext <- grep(x = file_name, pattern = ".h5ad$")

  if (length(x = file_ext) == 0) {
    file_name <- paste0(file_name, ".h5ad")
  }

  if (!is.null(x = file_path)) {
    norm_path <- normalizePath(path = file_path)
    full_path_name <- file.path(norm_path, file_name)
  } else {
    full_path_name <- file.path(file_name)
  }

  # Run update to ensure functionality
  if (isTRUE(x = verbose)) {
    cli_inform(message = c("*" = "Checking Seurat object validity & Extracting Data"))
  }

  # Check Seurat
  Is_Seurat(seurat_object = x)

  # Run update to ensure functionality
  x <- suppressMessages(UpdateSeuratObject(object = x))

  # Check Assay5 for multiple layers
  if (isTRUE(x = Assay5_Check(seurat_object = x, assay = assay))) {
    layers_check <- Layers(object = x, search = main_layer)
    if (length(x = layers_check) > 1) {
      cli_abort(message = c("Multiple data layers present {.field {head(x = layers_check, n = 2)}}.",
                            "i" = "Please run {.code JoinLayers} before converting to anndata object."))
    }
  }

  main_approved_slots <- Layers(object = x, search = c("counts", "data"))

  if (!main_layer %in% main_approved_slots) {
    cli_abort(message = "{.code main_layer} must be one of {.field {main_approved_slots}}")
  }

  if (main_layer %in% other_layers) {
    cli_abort(message = "{.code main_layer} and {.code other_layers} cannot overlap.")
  }

  if (isFALSE(x = all(other_layers %in% Layers(object = x)))) {
    cli_abort(message = "One or more of {.field {other_layers}} were not found in Seurat object.")
  }

  # Extract Data
  main_layer_data <- LayerData(object = x, assay = assay, layer = main_layer)

  meta_data <- Fetch_Meta(object = x)

  meta_data <- drop_single_value_cols(df = meta_data)

  if (isTRUE(x = Assay5_Check(seurat_object = x, assay = assay))) {
    seurat_var_info <- drop_single_value_cols(df = x[[assay]]@meta.data)
  } else {
    seurat_var_info <- drop_single_value_cols(df = x[[assay]]@meta.features)
  }

  if (isTRUE(x = transer_dimreduc)) {
    dim_reducs_present <- Reductions(object = x)
    if (length(x = dim_reducs_present) > 0) {
      dim_reducs_list <- lapply(dim_reducs_present, function(z) {
        as.matrix(x = Embeddings(object = x, reduction = z))
      })
      names(x = dim_reducs_list) <- paste0("X_", str_to_lower(string = dim_reducs_present))
    } else {
      dim_reducs_present <- NULL
    }
  } else {
    dim_reducs_present <- NULL
  }

  if (length(x = other_layers) > 0) {
    other_layers_list <- lapply(other_layers, function(i) {
      Matrix::t(LayerData(object = x, layer = i, assay = assay))
    })
    names(x = other_layers_list) <- other_layers
  } else {
    other_layers_list <- list()
  }

  # convert
  if (isTRUE(x = verbose)) {
    cli_inform(message = c("*" = "Creating anndata object."))
  }
  anndata <- reticulate::import("anndata", convert = FALSE)

  adata <- anndata$AnnData(
    X = Matrix::t(main_layer_data),
    obs = meta_data,
    var = seurat_var_info,
    obsm = dim_reducs_list,
    layers = other_layers_list
  )

  if (isTRUE(x = verbose)) {
    cli_inform(message = c("*" = "Writing anndata file: {.val {full_path_name}}"))
  }
  adata$write(full_path_name, compression = "gzip")

  adata
}


#' Create & Save Anndata Object
#'
#' @param file_path directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name file name.
#' @param transfer_norm.data logical, whether to transfer the norm.data in addition to
#' raw.data, default is FALSE.
#' @param reduction_label What to label the visualization dimensionality reduction.
#' LIGER does not store name of technique and therefore needs to be set manually.
#' @param add_barcode_names logical, whether to add dataset names to the cell barcodes when
#' merging object data, default is FALSE.
#' @param barcode_prefix logical, if `add_barcode_names = TRUE` should the names be added as
#' prefix to current cell barcodes/names or a suffix (default is TRUE; prefix).
#' @param barcode_cell_id_delimiter The delimiter to use when adding dataset id to barcode
#' prefix/suffix.  Default is "_".
#' @param verbose logical, whether to print status messages during object conversion (default is TRUE).
#'
#' @references LIGER version inspired by `sceasy::seurat2anndata` modified and updated to apply to LIGER objects (sceasy package: \url{https://github.com/cellgeni/sceasy}; License: GPL-3.
#'
#' @method as.anndata liger
#'
#' @concept object_conversion
#'
#' @import cli
#' @importFrom dplyr mutate
#' @importFrom magrittr "%>%"
#' @importFrom stringr str_to_lower
#' @importFrom tibble column_to_rownames
#'
#' @export
#' @rdname as.anndata
#'
#' @examples
#' \dontrun{
#' as.anndata(x = liger_object, file_path = "/folder_name", file_name = "anndata_converted.h5ad")
#' }
#'

as.anndata.liger <- function(
    x,
    file_path,
    file_name,
    transfer_norm.data = FALSE,
    reduction_label = NULL,
    add_barcode_names = FALSE,
    barcode_prefix = TRUE,
    barcode_cell_id_delimiter = "_",
    verbose = TRUE,
    ...
) {
  # Check hdf5r installed
  reticulate_check <- is_installed(pkg = "reticulate")
  if (isFALSE(x = reticulate_check)) {
    cli_abort(message = c(
      "Please install the {.val reticulate} package to use {.code as.anndata}.",
      "i" = "This can be accomplished with the following commands: ",
      "----------------------------------------",
      "{.field `install.packages({symbol$dquote_left}reticulate{symbol$dquote_right})`}",
      "----------------------------------------"
    ))
  }

  # Check all barcodes are unique to begin with
  duplicated_barcodes <- x@raw.data %>%
    lapply(colnames) %>%
    unlist() %>%
    duplicated() %>%
    any()

  if (isTRUE(x = duplicated_barcodes) && is.null(x = add_barcode_names)) {
    cli_abort(message = c("There are overlapping cell barcodes present in the input data",
                          "i" = "Please set {.code add_barcode_names = TRUE} to make all cell barcodes unique.")
    )
  }

  if (is.null(x = reduction_label)) {
    cli_abort(message = c("{.code reduction_label} parameter was not set.",
                          "*" = "LIGER objects do not store name of dimensionality reduction technique used.",
                          "i" = "In order to retain proper labels in Seurat object please set {.code reduction_label} to {.val tSNE}, {.val UMAP}, {.val etc}."))
  }

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path)) {
    if (!dir.exists(paths = file_path)) {
      cli_abort(message = "Provided {.code file_path}: {symbol$dquote_left}{.field {file_path}}{symbol$dquote_right} does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name)) {
    cli_abort(message = "No file name provided.  Please provide a file name using {.code file_name}.")
  }

  file_ext <- grep(x = file_name, pattern = ".h5ad$")

  if (length(x = file_ext) == 0) {
    file_name <- paste0(file_name, ".h5ad")
  }

  if (!is.null(x = file_path)) {
    norm_path <- normalizePath(path = file_path)
    full_path_name <- file.path(norm_path, file_name)
  } else {
    full_path_name <- file.path(file_name)
  }

  if (isTRUE(x = verbose)) {
    cli_inform(message = c("*" = "Creating main layer from {.field raw.data}"))
  }
  if (isTRUE(x = add_barcode_names)) {
    nms <-  names(x = x@H)
    main_layer_data <- Merge_Sparse_Data_All(matrix_list = x@raw.data, add_cell_ids = nms, prefix = barcode_prefix, cell_id_delimiter = barcode_cell_id_delimiter)
  } else {
    main_layer_data <- Merge_Sparse_Data_All(matrix_list = x@raw.data)
  }

  # merge norm data
  if (isTRUE(x = transfer_norm.data)) {
    cli_inform(message = c("*" = "Creating other layer from {.field norm.data}"))
    if (isTRUE(x = add_barcode_names)) {
      nms <-  names(x = x@H)
      norm_data <- Merge_Sparse_Data_All(matrix_list = x@norm.data, add_cell_ids = nms, prefix = barcode_prefix, cell_id_delimiter = barcode_cell_id_delimiter)
    } else {
      norm_data <- Merge_Sparse_Data_All(matrix_list = x@norm.data)
    }

    other_layers <- list("norm.data" = Matrix::t(norm_data)
    )
  } else {
    other_layers <- list()
  }

  # pull var genes
  liger_var_genes <- x@var.genes
  total_features <- data.frame("all_genes" = LIGER_Features(liger_object = x))

  liger_var_df <- total_features %>%
    mutate("variable_genes" = ifelse(.data[["all_genes"]] %in% liger_var_genes, .data[["all_genes"]], NA)) %>%
    column_to_rownames("all_genes")

  # Prep reductions
  if (all(dim(x = x@W) > 0) && all(dim(x = x@H.norm) > 0)) {
    inmf.loadings <- Matrix::t(x = x@W)
    rinmf.loadings <- inmf.loadings

    dimnames(x = inmf.loadings) <- list(x@var.genes,
                                        paste0("iNMF_", seq_len(ncol(inmf.loadings))))
    dimnames(x = rinmf.loadings) <- list(x@var.genes,
                                         paste0("rawiNMF_", seq_len(ncol(rinmf.loadings))))

    inmf.embeddings <- x@H.norm
    rinmf.embeddings <- do.call(what = 'rbind', args = x@H)

    dimnames(x = inmf.embeddings) <- list(unlist(x = lapply(x@scale.data, rownames), use.names = FALSE),
                                          paste0("iNMF_", seq_len(ncol(inmf.loadings))))
    dimnames(x = rinmf.embeddings) <- list(unlist(x = lapply(x@scale.data, rownames), use.names = FALSE),
                                           paste0("rawiNMF_", seq_len(ncol(x = inmf.loadings))))

    inmf.obj <- CreateDimReducObject(
      embeddings = inmf.embeddings,
      loadings = inmf.embeddings,
      global = TRUE,
      key = "iNMF_"
    )

    inmf_anndata <- as.matrix(x = Embeddings(object = inmf.obj))

    rinmf.obj <- CreateDimReducObject(
      embeddings = rinmf.embeddings,
      loadings = rinmf.loadings,
      key = "rawiNMF_",
      global = TRUE
    )

    rinmf_anndata <- as.matrix(x = Embeddings(object = rinmf.obj))
  }

  # prep visualization reduction
  dimreduc.embeddings <- x@tsne.coords
  dimnames(x = dimreduc.embeddings) <- list(rownames(x@H.norm),
                                            paste0(reduction_label, 1:2))

  # create reducs list
  reducs <- list(inmf_anndata,
                 rinmf_anndata,
                 dimreduc.embeddings)

  names(x = reducs) <- paste0("X_", str_to_lower(c("inmf", "rinmf", reduction_label)))

  # get meta and drop single value columns
  liger_meta <- Fetch_Meta(object = x)

  liger_meta <- drop_single_value_cols(df = liger_meta)

  # Create anndata
  if (isTRUE(x = verbose)) {
    cli_inform(message = c("*" = "Creating anndata object."))
  }
  anndata <- reticulate::import("anndata", convert = FALSE)

  adata <- anndata$AnnData(
    X = Matrix::t(main_layer_data),
    obs = liger_meta,
    var = liger_var_genes,
    obsm = reducs,
    layers = other_layers
  )

  # write file
  if (isTRUE(x = verbose)) {
    cli_inform(message = c("*" = "Writing anndata file: {.val {full_path_name}}"))
  }
  adata$write(full_path_name, compression = "gzip")

  adata
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### CONVERT SEURAT ASSAYS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Convert between Seurat Assay types
#'
#' Will convert assays within a Seurat object between "Assay" and "Assay5" types.
#'
#' @param seurat_object Seurat object name.
#' @param assay name(s) of assays to convert.  Default is NULL and will check with users
#' which assays they want to convert.
#' @param convert_to value of what assay type to convert current assays to.
#' #' \itemize{
#'       \item Accepted values for V3/4 are: "Assay", "assay", "V3", or "v3".
#'       \item Accepted values for V5 are: "Assay5", "assay5", "V5", or "v5".
#'       }
#'
#' @concept object_conversion
#'
#' @import cli
#' @importFrom dplyr mutate
#' @importFrom magrittr "%>%"
#' @importFrom stringr str_to_lower
#' @importFrom tibble column_to_rownames
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert to V3/4 assay
#' obj <- Convert_Assay(seurat_object = obj, convert_to = "V3")
#'
#' # Convert to 5 assay
#' obj <- Convert_Assay(seurat_object = obj, convert_to = "V5")
#' }
#'

Convert_Assay <- function(
    seurat_object,
    assay = NULL,
    convert_to
) {
  # Check accepted
  accepted_V3 <- c("Assay", "assay", "V3", "v3", "3")
  accepted_V5 <- c("Assay5", "assay5", "V5", "v5", "5")

  # convert to character in case numeric provided
  convert_to <- as.character(x = convert_to)

  if (!convert_to %in% c(accepted_V5, accepted_V3)) {
    cli_abort(message = c("Value provided to {.code convert_to} ({.field {convert_to}}) was not accepted value.",
                          "i" = "Accepted values to convert to V3/4 are: {.field {accepted_V3}}",
                          "i" = "Accepted values to convert to V5 are: {.field {accepted_V5}}"))
  }

  # set assay value
  if (convert_to %in% accepted_V5) {
    if (packageVersion(pkg = 'Seurat') < 5) {
      cli_abort(message = "Seurat version must be v5.0.0 or greater to create Assay5.")
    }

    convert_to <- "Assay5"
    convert_from <- "Assay"
  }
  if (convert_to %in% accepted_V3) {
    convert_to <- "Assay"
    convert_from <- "Assay5"
  }

  if (is.null(x = assay)) {
    num_assays <- length(x = Assays(object = seurat_object))
    if (num_assays > 1) {
      if (yesno("Multiple assays ({.field {Assays(object = seurat_object)}}) are present.  Should all assays be converted to assay type: {.field {convert_to}}?")) {
        cli_inform(message = c("!" = "To only convert specific assays please specify assay names using {.code assay} parameter."))
        return(invisible())
      }
    }
  }

  # Check assays are present if provided
  if (!is.null(x = assay)) {
    assays_not_found <- Assay_Present(seurat_object = seurat_object, assay_list = assay, print_msg = FALSE, omit_warn = TRUE)[[2]]

    if (!is.null(x = assays_not_found)) {
      stop_quietly()
    }
  }

  # set assays to convert
  assays_convert <- assay %||% Assays(object = seurat_object)

  # Check against current assay class
  current_assay_classes <- sapply(assays_convert, function(x) {
    class(x = seurat_object[[x]])
  })

  if (convert_to %in% current_assay_classes) {
    cli_abort(message = c("Attempting to assays to {.field {convert_to}}, however one or more of current assays is already of that type",
                          "i" = "Check assay type and/or whether {.code {convert_to}} value is correct."))
  }

  if ("SCTAssay" %in% current_assay_classes) {
    cli_abort(message = "Cannot convert assay of class {.field SCTAssay}.")
  }

  # convert assays
  for (i in assays_convert) {
    cli_inform(message = "Converting assay {.val {i}} from {.field {convert_from}} to {.field {convert_to}}.")
    suppressWarnings(seurat_object[[i]] <- as(seurat_object[[i]], convert_to))
  }

  # return object
  return(seurat_object)
}


#' Split Seurat object into layers
#'
#' Split Assay5 of Seurat object into layers by variable in meta.data
#'
#' @param seurat_object Seurat object name.
#' @param assay name(s) of assays to convert.  Defaults to current active assay.
#' @param split.by Variable in meta.data to use for splitting layers.
#'
#' @concept object_conversion
#'
#' @import cli
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Split object by "treatment"
#' obj <- Split_Layers(object = obj, assay = "RNA", split.by = "treatment")
#' }
#'

Split_Layers <- function(
    seurat_object,
    assay = "RNA",
    split.by
) {
  # Make sure single assay
  if (length(x = assay) > 1) {
    cli_abort(message = c("Multiple assays specified ({.field {assay}}).",
                          "i" = "{.code Split_Layers} only supports splitting one assay at a time."))
  }

  # Check assay present
  assay_present <- Assay_Present(seurat_object = seurat_object, assay_list = assay, print_msg = FALSE, omit_warn = TRUE)[[1]]

  # check split is valid and length
  split.by <- Meta_Present(object = seurat_object, meta_col_names = split.by, print_msg = FALSE, omit_warn = FALSE)[[1]]

  length_split <- length(x = unique(x = seurat_object@meta.data[[split.by]]))

  # Check for already split layers
  check_split <- Layers(object = seurat_object, search = "counts", assay = assay_present)

  if (length(x = check_split) > 1) {
    cli_warn(message = "Layers in the assay: {.field {assay_present}} already appear split.  Skipping assay.")
  } else {
    cli_inform(message = c("*" = "Splitting layers within assay: {.field {assay_present}} into {.field {length_split} parts} by {.val {split.by}}"))
    # Check v3 vs. v5 and convert if needed
    if (isFALSE(x = Assay5_Check(seurat_object = seurat_object, assay = assay_present))) {
      cli_inform(message = c("i" = "{.field {assay_present}} is not Assay5, converting to Assay5 before splitting."))

      seurat_object <- suppressMessages(Convert_Assay(seurat_object = seurat_object, assay = assay_present, convert_to = "V5"))
    }

    # split layers
    seurat_object[[assay_present]] <- split(seurat_object[[assay_present]], f = seurat_object@meta.data[[split.by]])
  }

  # return object
  return(seurat_object)
}
