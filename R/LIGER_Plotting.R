#' DimPlot LIGER Version
#'
#' Standard and modified version of LIGER's plotByDatasetAndCluster
#'
#' @param liger_object \code{liger} liger_object.  Need to perform clustering before calling this function
#' @param group_by `r lifecycle::badge("deprecated")` soft-deprecated. See `group.by`.
#' @param group.by Variable to be plotted.  If `NULL` will plot clusters from `liger@clusters` slot.
#' If `combination = TRUE` will plot both clusters and meta data variable.
#' If `combination = TRUE` will plot both clusters and meta data variable.
#' @param split_by `r lifecycle::badge("deprecated")` soft-deprecated. See `split.by`.
#' @param split.by Variable to split plots by.
#' @param colors_use_cluster colors to use for plotting by clusters.  By default if number of levels plotted is
#' less than or equal to 36 will use "polychrome" and if greater than 36 will use "varibow" with shuffle = TRUE
#' both from \code{\link{DiscretePalette_scCustomize}}.
#' @param colors_use_meta colors to use for plotting by meta data (cell.data) variable.  By default if number
#' of levels plotted is less than or equal to 36 it will use "polychrome" and if greater than 36 will use
#' "varibow" with shuffle = TRUE both from DiscretePalette_scCustomize.
#' @param pt_size Adjust point size for plotting.
#' @param shuffle logical. Whether to randomly shuffle the order of points. This can be useful for crowded plots
#' if points of interest are being buried. (Default is TRUE).
#' @param shuffle_seed Sets the seed if randomly shuffling the order of points.
#' @param reduction_label What to label the x and y axes of resulting plots.  LIGER does not store name of
#' technique and therefore needs to be set manually.  Default is "UMAP". (only valid for
#' rliger < 2.0.0).
#' @param reduction specify reduction to use when plotting.  Default is current object
#' default reduction (only valid for rliger v2.0.0 or greater).
#' @param aspect_ratio Control the aspect ratio (y:x axes ratio length).  Must be numeric value;
#' Default is NULL.
#' @param label logical.  Whether or not to label the clusters.  ONLY applies to plotting by cluster.  Default is TRUE.
#' @param label_size size of cluster labels.
#' @param label_repel logical.  Whether to repel cluster labels from each other if plotting by
#' cluster (if `group.by = NULL` or `group.by = "cluster`).  Default is FALSE.
#' @param label_box logical.  Whether to put a box around the label text (uses `geom_text` vs `geom_label`).
#' Default is FALSE.
#' @param label_color Color to use for cluster labels.  Default is "black".
#' @param combination logical, whether to return patchwork displaying both plots side by side.  (Default is FALSE).
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param num_columns Number of columns in plot layout.  Only valid if `split.by != NULL`.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#'
#' @return A ggplot/patchwork object
#'
#' @import ggplot2

#' @importFrom patchwork wrap_plots
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @concept liger_plotting
#'
#' @examples
#' \dontrun{
#' DimPlot_LIGER(liger_object = obj_name, reduction_label = "UMAP")
#' }
#'

DimPlot_LIGER <- function(
  liger_object,
  group_by = deprecated(),
  group.by = NULL,
  split_by = deprecated(),
  split.by = NULL,
  colors_use_cluster = NULL,
  colors_use_meta = NULL,
  pt_size = NULL,
  shuffle = TRUE,
  shuffle_seed = 1,
  reduction_label = "UMAP",
  reduction = NULL,
  aspect_ratio = NULL,
  label = TRUE,
  label_size = NA,
  label_repel = FALSE,
  label_box = FALSE,
  label_color = "black",
  combination = FALSE,
  raster = NULL,
  raster.dpi = c(512, 512),
  num_columns = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123
) {
  # Check LIGER
  Is_LIGER(liger_object = liger_object)

  # check deprecation
  if (is_present(group_by)) {
    deprecate_warn(when = "3.1.0",
                              what = "DimPlot_LIGER(group_by)",
                              details = c("i" = "The {.code group_by} parameter is soft-deprecated.  Please update code to use `group.by` instead.")
    )
    group.by <- group_by
  }

  # Set group.by defaults
  if (isFALSE(x = combination) && is.null(x = group.by)) {
    group.by <- "cluster"
  }

  if (isTRUE(x = combination) && is.null(x = group.by)) {
    group.by <- "dataset"
  }

  # Group by cluster options
  cluster_options <- c("cluster", "Cluster", "clusters", "Clusters")
  if (group.by %in% cluster_options) {
    group.by <- "cluster"
  }

  # Check group.by parameter
  if (!group.by == "cluster")
    group_by_var <- Meta_Present(object = liger_object, meta_col_names = group.by, print_msg = FALSE, omit_warn = FALSE)[[1]]

  if (!is.null(x = split.by)) {
    group_by_var <- Meta_Present(object = liger_object, meta_col_names = split.by, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  if (packageVersion(pkg = 'rliger') < "2.0.0") {
    # Add one time dim label warning
    if (getOption(x = 'scCustomize_warn_LIGER_dim_labels', default = TRUE)) {
      cli_inform(message = c("",
                             "NOTE: {.field DimPlot_LIGER} uses the {.code reduction_label} parameter to set axis labels ",
                             "on the plot.",
                             "By default this is set to {.val UMAP}.",
                             "Please take note of this parameter as LIGER objects do not store the name",
                             "of reduction technique used and therefore this needs to be set manually.",
                             "",
                             "-----This message will be shown once per session.-----"))
      options(scCustomize_warn_LIGER_dim_labels = FALSE)
    }
  }

  # cells in object
  cells_total <- Cells(x = liger_object)

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = cells_total) > 2e5)

  if (isTRUE(x = raster) && (length(x = cells_total) > 2e5) && getOption(x = 'scCustomize_warn_raster_LIGER', default = TRUE)) {
    cli_inform(message = c("",
                           "Rasterizing points since number of points exceeds 200,000.",
                           "To disable this behavior set {.code raster = FALSE}",
                           "",
                           "-----This message will be shown once per session.-----"))
    options(scCustomize_warn_raster_LIGER = FALSE)
  }

  # Add point size
  if (is.null(x = pt_size)) {
    # modified version of the AutoPointSize() function from Seurat
    pt_size <- AutoPointSize_scCustom(data = cells_total, raster = raster)
  }

  # liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    plots <- LIGER2_DimPlot(liger_object = liger_object,
                           group.by = group.by,
                           split.by = split.by,
                           colors_use_cluster = colors_use_cluster,
                           colors_use_meta = colors_use_meta,
                           pt_size = pt_size,
                           shuffle = shuffle,
                           shuffle_seed = shuffle_seed,
                           reduction = reduction,
                           aspect_ratio = aspect_ratio,
                           label = label,
                           label_size = label_size,
                           label_repel = label_repel,
                           label_box = label_box,
                           label_color = label_color,
                           combination = combination,
                           raster = raster,
                           raster.dpi = raster.dpi,
                           num_columns = num_columns,
                           ggplot_default_colors = ggplot_default_colors,
                           color_seed = color_seed)
  } else {
    plots <- LIGER_DimPlot(liger_object = liger_object,
                           group.by = group.by,
                           split.by = split.by,
                           colors_use_cluster = colors_use_cluster,
                           colors_use_meta = colors_use_meta,
                           pt_size = pt_size,
                           shuffle = shuffle,
                           shuffle_seed = shuffle_seed,
                           reduction_label = reduction_label,
                           aspect_ratio = aspect_ratio,
                           label = label,
                           label_size = label_size,
                           label_repel = label_repel,
                           label_box = label_box,
                           label_color = label_color,
                           combination = combination,
                           raster = raster,
                           raster.dpi = raster.dpi,
                           num_columns = num_columns,
                           ggplot_default_colors = ggplot_default_colors,
                           color_seed = color_seed)
  }
  # return plots
  return(plots)
}


#' Customized version of plotFactors
#'
#' Modified and optimized version of `plotFactors` function from LIGER package.
#'
#' @param liger_object \code{liger} liger_object.  Need to perform clustering and factorization before calling this function
#' @param num_genes Number of genes to display for each factor (Default 8).
#' @param colors_use_factors colors to use for plotting factor loadings  By default datasets will be
#' plotted using "varibow" with shuffle = TRUE from both from \code{\link{DiscretePalette_scCustomize}}.
#' @param colors_use_dimreduc colors to use for plotting factor loadings on dimensionality reduction
#' coordinates (tSNE/UMAP).  Default is c('lemonchiffon', 'red'),
#' @param pt.size_factors Adjust point size for plotting in the factor plots.
#' @param pt.size_dimreduc Adjust point size for plotting in dimensionality reduction plots.
#' @param reduction Name of dimensionality reduction to use for plotting.  Default is "UMAP".
#' Only for newer style liger objects.
#' @param reduction_label What to label the x and y axes of resulting plots.  LIGER does not store name of
#' technique and therefore needs to be set manually.  Default is "UMAP".
#' Only for older style liger objects.
#' @param plot_legend logical, whether to plot the legend on factor loading plots, default is TRUE.
#' Helpful if number of datasets is large to avoid crowding the plot with legend.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param order logical. Whether to plot higher loading cells on top of cells with lower loading values in the
#' dimensionality reduction plots (Default = FALSE).
#' @param plot_dimreduc logical.  Whether to plot factor loadings on dimensionality reduction coordinates.  Default is TRUE.
#' @param save_plots logical.  Whether to save plots.  Default is TRUE
#' @param file_path directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix to append after sample name.
#' @param return_plots logical. Whether or not to return plots to the environment.  (Default is FALSE)
#' @param cells.highlight Names of specific cells to highlight in plot (black) (default NULL).
#' @param reorder_datasets New order to plot datasets in for the factor plots if different from current
#' factor level order in cell.data slot. Only for older style liger objects.
#' @param ggplot_default_colors logical.  If `colors_use_factors = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "varibow" palette.
#' @param color_seed random seed for the palette shuffle if `colors_use_factors = NULL`.  Default = 123.
#'
#' @return A list of ggplot/patchwork objects and/or PDF file.
#'
#' @import cli
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom patchwork wrap_plots
#' @importFrom scattermore geom_scattermore
#'
#' @export
#'
#' @concept liger_plotting
#'
#' @author Velina Kozareva (Original code for modified function), Sam Marsh (Added/modified functionality)
#' @references Based on `plotFactors` functionality from original LIGER package.
#'
#' @examples
#' \dontrun{
#' plotFactors_scCustom(liger_object = liger_obj, return_plots = FALSE, plot_dimreduc = TRUE,
#' raster = FALSE, save_plots = TRUE)
#' }
#'

plotFactors_scCustom <- function(
  liger_object,
  num_genes = 8,
  colors_use_factors = NULL,
  colors_use_dimreduc = c('lemonchiffon', 'red'),
  pt.size_factors = 1,
  pt.size_dimreduc = 1,
  reduction = "UMAP",
  reduction_label = "UMAP",
  plot_legend = TRUE,
  raster = TRUE,
  raster.dpi = c(512, 512),
  order = FALSE,
  plot_dimreduc = TRUE,
  save_plots = TRUE,
  file_path = NULL,
  file_name = NULL,
  return_plots = FALSE,
  cells.highlight = NULL,
  reorder_datasets = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123
) {
  # Check LIGER
  Is_LIGER(liger_object = liger_object)

  # rliger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    plotFactors_liger2_scCustom(liger_object = liger_object,
                                num_genes = num_genes,
                                colors_use_factors = colors_use_factors,
                                colors_use_dimreduc = colors_use_dimreduc,
                                pt.size_factors = pt.size_factors,
                                pt.size_dimreduc = pt.size_dimreduc,
                                reduction = reduction,
                                plot_legend = plot_legend,
                                raster = raster,
                                raster.dpi = raster.dpi,
                                order = order,
                                plot_dimreduc = plot_dimreduc,
                                save_plots = save_plots,
                                file_path = file_path,
                                file_name = file_name,
                                return_plots = return_plots,
                                cells.highlight = cells.highlight,
                                ggplot_default_colors = ggplot_default_colors,
                                color_seed = color_seed
    )
  } else {
    plotFactors_liger_scCustom(liger_object = liger_object,
                                num_genes = num_genes,
                                colors_use_factors = colors_use_factors,
                                colors_use_dimreduc = colors_use_dimreduc,
                                pt.size_factors = pt.size_factors,
                                pt.size_dimreduc = pt.size_dimreduc,
                                reduction_label = reduction_label,
                                plot_legend = plot_legend,
                                raster = raster,
                                raster.dpi = raster.dpi,
                                order = order,
                                plot_dimreduc = plot_dimreduc,
                                save_plots = save_plots,
                                file_path = file_path,
                                file_name = file_name,
                                return_plots = return_plots,
                                cells.highlight = cells.highlight,
                                reorder_datasets = reorder_datasets,
                                ggplot_default_colors = ggplot_default_colors,
                                color_seed = color_seed
    )
  }
}
