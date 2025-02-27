#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GENERICS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert objects to LIGER objects
#'
#' Convert objects (Seurat & lists of Seurat Objects) to anndata objects
#'
#' @param x An object to convert to class `liger`
#' @param ... Arguments passed to other methods
#'
#' @return a liger object generated from `x`
#'
#' @rdname as.LIGER
#' @export as.LIGER
#'

as.LIGER <- function(x, ...) {
  UseMethod(generic = "as.LIGER", object = x)
}


#' Convert objects to anndata objects
#'
#' Convert objects (Seurat & LIGER) to anndata objects
#'
#' @param x Seurat or LIGER object
#' @param ... Arguments passed to other methods
#'
#' @return an anndata object generated from `x`, saved at path provided.
#'
#' @rdname as.anndata
#' @export as.anndata
#'

as.anndata <- function(x, ...) {
  UseMethod(generic = "as.anndata", object = x)
}


#' Add Mito and Ribo percentages
#'
#' Add Mito, Ribo, & Mito+Ribo percentages to meta.data slot of Seurat Object or
#' cell.data slot of Liger object
#'
#' @param object Seurat or LIGER object
#' @param ... Arguments passed to other methods
#'
#' @return An object of the same class as `object` with columns added to object meta data.
#'
#' @rdname Add_Mito_Ribo
#' @export Add_Mito_Ribo
#'

Add_Mito_Ribo <- function(object, ...) {
  UseMethod(generic = "Add_Mito_Ribo", object = object)
}


#' Add Hemoglobin percentages
#'
#' Add hemoglobin percentages to meta.data slot of Seurat Object or
#' cell.data/cellMeta slot of Liger object
#'
#' @param object Seurat or LIGER object
#' @param ... Arguments passed to other methods
#'
#' @return An object of the same class as `object` with columns added to object meta data.
#'
#' @rdname Add_Hemo
#' @export Add_Hemo
#'

Add_Hemo <- function(object, ...) {
  UseMethod(generic = "Add_Hemo", object = object)
}


#' Add Cell Complexity
#'
#' Add measure of cell complexity/novelty (log10GenesPerUMI) for data QC.
#'
#' @param object Seurat or LIGER object
#' @param ... Arguments passed to other methods
#'
#' @return An object of the same class as `object` with columns added to object meta data.
#'
#' @rdname Add_Cell_Complexity
#' @export Add_Cell_Complexity
#'

Add_Cell_Complexity <- function(object, ...) {
  UseMethod(generic = "Add_Cell_Complexity", object = object)
}


#' Add Percent of High Abundance Genes
#'
#' Add the percentage of counts occupied by the top XX most highly expressed genes in each cell.
#'
#' @param object Seurat or LIGER object.
#' @param ... Arguments passed to other methods
#'
#' @rdname Add_Top_Gene_Pct
#' @export Add_Top_Gene_Pct
#'

Add_Top_Gene_Pct <- function(object, ...) {
  UseMethod(generic = "Add_Top_Gene_Pct", object = object)
}


#' Add MALAT1 QC Threshold
#'
#' Adds TRUE/FALSE values to each cell based on calculation of MALAT1 threshold.
#' This function incorporates a threshold calculation and procedure as described in
#' Clarke & Bader (2024). bioRxiv \url{doi.org/10.1101/2024.07.14.603469}.  Please cite this preprint
#' whenever using this function.
#'
#' @param object Seurat or LIGER object
#' @param ... Arguments passed to other methods
#'
#' @rdname Add_MALAT1_Threshold
#' @export Add_MALAT1_Threshold
#'

Add_MALAT1_Threshold <- function(object, ...) {
  UseMethod(generic = "Add_MALAT1_Threshold", object = object)
}


#' Add Multiple Cell Quality Control Values with Single Function
#'
#' Add Mito/Ribo %, Cell Complexity (log10GenesPerUMI), Top Gene Percent with single
#' function call to Seurat or liger objects.
#'
#' @param object Seurat or LIGER object
#' @param ... Arguments passed to other methods
#'
#' @rdname Add_Cell_QC_Metrics
#' @export Add_Cell_QC_Metrics
#'

Add_Cell_QC_Metrics <- function(object, ...) {
  UseMethod(generic = "Add_Cell_QC_Metrics", object = object)
}


#' Get meta data from object
#'
#' Quick function to properly pull meta.data from objects.
#'
#' @param object Object of class Seurat or liger.
#' @param ... Arguments passed to other methods
#'
#' @return A data.frame containing cell-level meta data
#'
#' @export
#'
#' @concept get_set_util
#'
#' @rdname Fetch_Meta
#'
#' @examples
#' library(Seurat)
#' meta_data <- Fetch_Meta(object = pbmc_small)
#' head(meta_data, 5)
#'

Fetch_Meta <- function(object, ...) {
  UseMethod(generic = "Fetch_Meta", object = object)
}


#' Rename Clusters
#'
#' Wrapper function to rename active cluster identity in Seurat or Liger Object with new idents.
#'
#' @param object Object of class Seurat or liger.
#' @param ... Arguments passed to other methods
#'
#' @return An object of the same class as `object` with updated default identities.
#'
#' @export
#'
#' @concept get_set_util
#'
#' @rdname Rename_Clusters
#'

Rename_Clusters <- function(object, ...) {
  UseMethod(generic = "Rename_Clusters", object = object)
}
