#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GENERICS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert objects to LIGER objects
#'
#' @param x Seurat object
#' @param ... Arguments passed to other methods
#'
#' @rdname as.LIGER
#' @export as.LIGER
#'

as.LIGER <- function(x, ...) {
  UseMethod(generic = 'as.LIGER', object = x)
}


#' Convert objects to anndata objects
#'
#' This function is part of generic `as.anndata` for conversion of Seurat Objects to anndata objects.
#'
#' @param x Seurat or LIGER object
#' @param ... Arguments passed to other methods
#'
#' @return saved anndata object to at path provided.
#'
#' @rdname as.anndata
#' @export as.anndata
#'

as.anndata <- function(x, ...) {
  UseMethod(generic = 'as.anndata', object = x)
}


#' Add Mitochondrial & Ribosomal Percentages
#'
#' @param object Seurat or LIGER object
#' @param ... Arguments passed to other methods
#'
#' @rdname Add_Mito_Ribo
#' @export Add_Mito_Ribo
#'

Add_Mito_Ribo <- function(object, ...) {
  UseMethod(generic = 'Add_Mito_Ribo', object = object)
}


#' Add Cell Complexity
#'
#' @param object Seurat or LIGER object
#' @param ... Arguments passed to other methods
#'
#' @rdname Add_Cell_Complexity
#' @export Add_Cell_Complexity
#'

Add_Cell_Complexity <- function(object, ...) {
  UseMethod(generic = 'Add_Cell_Complexity', object = object)
}
