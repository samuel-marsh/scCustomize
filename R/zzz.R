#' @section Package options:
#'
#' scCustomize uses the following [options()] to configure behavior:
#'
#' \describe{
#'   \item{\code{scCustomize_warn_raster_iterative}}{Show message about setting `raster` parameter
#'   in \code{\link{Iterate_FeaturePlot_scCustom}} if `raster = FALSE` and `single_pdf = TRUE`
#'   due to large file sizes.}
#'   \item{\code{scCustomize_warn_raster_LIGER}}{Show warning about rasterization of points in
#'   \code{\link{DimPlot_LIGER}} due to new functionality compared to LIGER.}
#'   \item{\code{scCustomize_warn_na_cutoff}}{Show message about properly setting `na_cutoff` parameter in \code{\link{FeaturePlot_scCustom}}.}
#'   #'   \item{\code{scCustomize_warn_zero_na_cutoff}}{Show message about properly setting `na_cutoff` parameter in \code{\link{FeaturePlot_scCustom}} if `na_cutoff` is set to exactly zero.}
#'   \item{\code{scCustomize_warn_vln_raster_iterative}}{Show message about \code{\link{Iterate_VlnPlot_scCustom}}
#'    when `pt.size > 0` due to current lack of raster support in \code{\link[Seurat]{VlnPlot}}}
#'    \item{\code{scCustomize_warn_LIGER_dim_labels}}{Show message about \code{\link{DimPlot_LIGER}}
#'    parameter `reduction_label` as LIGER objects do not store dimensionality reduction name and
#'    and therefore needs to be set manually.}
#'    \item{\code{scCustomize_warn_DimPlot_split_type}}{Show message about \code{\link{DimPlot_scCustom}}
#'    parameter `split.by` and `split_seurat` to alert user to difference in returned plots between
#'    scCustomize and Seurat.}
#'    \item{\code{scCustomize_warn_FeatureScatter_split_type}}{Show message about \code{\link{FeatureScatter_scCustom}}
#'    parameter `split.by` and `split_seurat` to alert user to difference in returned plots between
#'    scCustomize and Seurat.}
#'    \item{\code{scCustomize_warn_LIGER_dim_labels_plotFactors}}{Show message about \code{\link{plotFactors_scCustom}}
#'    parameter `reduction_label` as LIGER objects do not store dimensionality reduction name and
#'    and therefore needs to be set manually.}
#' }
#'
#' @import cli
#' @importFrom lifecycle deprecated deprecate_soft deprecate_stop
#' deprecate_warn is_present
#' @importFrom rlang check_installed is_installed sym "!!"
#'
#' @keywords internal
#' @docType package
#' @rdname scCustomize-package
#' @name scCustomize-package
#'
"_PACKAGE"


scCustomize_default_options <- list(
  scCustomize_warn_raster_iterative = TRUE,
  scCustomize_warn_raster_LIGER = TRUE,
  scCustomize_warn_na_cutoff = TRUE,
  scCustomize_warn_zero_na_cutoff = TRUE,
  scCustomize_warn_vln_raster_iterative = TRUE,
  scCustomize_warn_LIGER_dim_labels = TRUE,
  scCustomize_warn_LIGER_dim_labels_plotFactors = TRUE,
  scCustomize_warn_DimPlot_split_type = TRUE,
  scCustomize_warn_FeatureScatter_split_type = TRUE
)



.onAttach <- function(libname, pkgname) {
  packageStartupMessage(cli::format_message(c("scCustomize v{packageVersion(pkg = 'scCustomize')}",
                                              "If you find the scCustomize useful please cite.",
                                              "See 'samuel-marsh.github.io/scCustomize/articles/FAQ.html' for citation info.")))
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(x = scCustomize_default_options) %in% names(x = op))
  if (any(toset)) options(scCustomize_default_options[toset])
  invisible(x = NULL)
}

