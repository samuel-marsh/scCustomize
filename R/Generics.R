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
  UseMethod(generic = 'as.LIGER', object = x)
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
  UseMethod(generic = 'as.anndata', object = x)
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
  UseMethod(generic = 'Add_Mito_Ribo', object = object)
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
  UseMethod(generic = 'Add_Cell_Complexity', object = object)
}
