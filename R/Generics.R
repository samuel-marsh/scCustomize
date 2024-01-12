#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GENERICS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert objects to LIGER objects
#'
#' @param x An object(s) to convert to class \code{LIGER}
#' @param ... Arguments passed to other methods
#'
#' @rdname as.LIGER
#' @export as.LIGER
#'

as.LIGER <- function(x, ...) {
  UseMethod(generic = 'as.LIGER', object = x)
}

#' Convert objects to Seurat objects
#'
#' @param x An object(s) to convert to class \code{Seurat}
#' @param ... Arguments passed to other methods
#'
#' @rdname as.Seurat
#' @importMethodsFrom SeuratObject as.Seurat
#' @export
#'

as.Seurat <- function(x, ...) {
  UseMethod(generic = 'as.Seurat', object = x)
}
