#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### TEMP DEPRECATED ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Deprecated functions
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' Use [Add_Mito_Ribo()] instead of `Add_Mito_Ribo_Seurat()`.
#'
#' @export
#' @keywords internal
#' @name deprecated

Add_Mito_Ribo_Seurat <- function(...) {
  lifecycle::deprecate_stop(when = "2.1.0", what = "Add_Mito_Ribo_Seurat()", with = "Add_Mito_Ribo()")
}
