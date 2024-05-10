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
  lifecycle::deprecate_stop(when = "2.0.0", what = "Split_FeatureScatter()", with = "FeatureScatter_scCustom()", details = "Deprecation error when calling function will be removed in v2.3.0+")
}


#' @description
#' Use [Add_Mito_Ribo()] instead of `Add_Mito_Ribo_Seurat()`.
#'
#' @export
#' @keywords internal
#' @name deprecated

Add_Mito_Ribo_Seurat <- function(...) {
  lifecycle::deprecate_stop(when = "2.1.0", what = "Add_Mito_Ribo_Seurat()", with = "Add_Mito_Ribo()", details = "Deprecation error when calling function will be removed in v2.3.0+")
}


#' @description
#' Use [Add_Mito_Ribo()] instead of `Add_Mito_Ribo_LIGER()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated

Add_Mito_Ribo_LIGER <- function(...) {
  lifecycle::deprecate_stop(when = "2.1.0", what = "Add_Mito_Ribo_LIGER()", with = "Add_Mito_Ribo()", details = "Deprecation error when calling function will be removed in v2.3.0+")
}


#' @description
#' Use [Add_Cell_Complexity()] instead of `Add_Cell_Complexity_Seurat()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated

Add_Cell_Complexity_Seurat <- function(...) {
  lifecycle::deprecate_stop(when = "2.1.0", what = "Add_Cell_Complexity_Seurat()", with = "Add_Cell_Complexity()", details = "Deprecation error when calling function will be removed in v2.3.0+")
}


#' @description
#' Use [Add_Cell_Complexity()] instead of `Add_Cell_Complexity_LIGER()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated

Add_Cell_Complexity_LIGER <- function(...) {
  lifecycle::deprecate_stop(when = "2.1.0", what = "Add_Cell_Complexity_LIGER()", with = "Add_Cell_Complexity()", details = "Deprecation error when calling function will be removed in v2.3.0+")
}


#' @description
#' Use [Meta_Present()] instead of `Meta_Present_LIGER()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated

Meta_Present_LIGER <- function(...) {
  lifecycle::deprecate_stop(when = "2.1.0", what = "Meta_Present_LIGER()", with = "Meta_Present()", details = "Deprecation error when calling function will be removed in v2.3.0+")
}


#' @description
#' Use [Add_Top_Gene_Pct()] instead of `Add_Top_Gene_Pct_Seurat()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated

Add_Top_Gene_Pct_Seurat <- function(...) {
  lifecycle::deprecate_stop(when = "2.2.0", what = "Add_Top_Gene_Pct_Seurat()", with = "Add_Top_Gene_Pct()", details = "Deprecation error when calling function will be removed in v2.4.0+")
}

#' @description
#' Use [Feature_Present()] instead of `Gene_Present()`.
#'
#' @export
#' @keywords internal
#' @rdname deprecated

Gene_Present <- function(...) {
  lifecycle::deprecate_stop(when = "2.2.0", what = "Gene_Present()", with = "Feature_Present()", details = "Deprecation error when calling function will be removed in v2.3.0+")
}
