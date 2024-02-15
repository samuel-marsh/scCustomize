#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### DEPRECATED FUNCTIONS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Deprecated functions `r lifecycle::badge("deprecated")`
#'
#'
#' @description
#' Use [FeatureScatter_scCustom()] instead of `Split_FeatureScatter()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated

Split_FeatureScatter <- function(...) {
  lifecycle::deprecate_stop(when = "2.0.0", what = "Split_FeatureScatter()", with = "FeatureScatter_scCustom()")
}


#' @description
#' Use [Add_Mito_Ribo()] instead of `Add_Mito_Ribo_Seurat()`.
#'
#' @export
#' @keywords internal
#' @name deprecated

Add_Mito_Ribo_Seurat <- function(...) {
  lifecycle::deprecate_stop(when = "2.1.0", what = "Add_Mito_Ribo_Seurat()", with = "Add_Mito_Ribo()")
}


#' @description
#' Use [Add_Mito_Ribo()] instead of `Add_Mito_Ribo_LIGER()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated

Add_Mito_Ribo_LIGER <- function(...) {
  lifecycle::deprecate_stop(when = "2.1.0", what = "Add_Mito_Ribo_LIGER()", with = "Add_Mito_Ribo()")
}


#' @description
#' Use [Add_Cell_Complexity()] instead of `Add_Cell_Complexity_Seurat()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated

Add_Cell_Complexity_Seurat <- function(...) {
  lifecycle::deprecate_stop(when = "2.1.0", what = "Add_Cell_Complexity_Seurat()", with = "Add_Cell_Complexity()")
}


#' @description
#' Use [Add_Cell_Complexity()] instead of `Add_Cell_Complexity_LIGER()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated

Add_Cell_Complexity_LIGER <- function(...) {
  lifecycle::deprecate_stop(when = "2.1.0", what = "Add_Cell_Complexity_LIGER()", with = "Add_Cell_Complexity()")
}


#' @description
#' Use [Meta_Present()] instead of `Meta_Present_LIGER()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated

Meta_Present_LIGER <- function(...) {
  lifecycle::deprecate_stop(when = "2.1.0", what = "Meta_Present_LIGER()", with = "Meta_Present()")
}
