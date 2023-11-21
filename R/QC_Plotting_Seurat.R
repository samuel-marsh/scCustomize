#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### OBJECT QC VLN ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' QC Plots Genes
#'
#' Custom VlnPlot for initial QC checks including lines for thresholding
#'
#' @param seurat_object Seurat object name.
#' @param plot_title Plot Title.
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident);
#' default is the current active.ident of the object.
#' @param y_axis_label Label for y axis.
#' @param x_axis_label Label for x axis.
#' @param low_cutoff Plot line a potential low threshold for filtering.
#' @param high_cutoff Plot line a potential high threshold for filtering.
#' @param pt.size Point size for plotting.
#' @param plot_median logical, whether to plot median for each ident on the plot (Default is FALSE).
#' @param median_size Shape size for the median is plotted.
#' @param plot_boxplot logical, whether to plot boxplot inside of violin (Default is FALSE).
#' @param colors_use vector of colors to use for plot.
#' @param x_lab_rotate Rotate x-axis labels 45 degrees (Default is TRUE).
#' @param y_axis_log logical. Whether to change y axis to log10 scale (Default is FALSE).
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 100,000 total points plotted (# Cells x # of features).
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ... Extra parameters passed to \code{\link[Seurat]{VlnPlot}}.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom Seurat VlnPlot
#'
#' @export
#'
#' @concept object_qc_plotting
#'
#' @examples
#' library(Seurat)
#' QC_Plots_Genes(seurat_object = pbmc_small, plot_title = "Genes per Cell", low_cutoff = 40,
#' high_cutoff = 85)
#'

QC_Plots_Genes <- function(
  seurat_object,
  plot_title = "Genes Per Cell/Nucleus",
  group.by = NULL,
  x_axis_label = NULL,
  y_axis_label = "Features",
  low_cutoff = NULL,
  high_cutoff = NULL,
  pt.size = NULL,
  plot_median = FALSE,
  plot_boxplot = FALSE,
  median_size = 15,
  colors_use = NULL,
  x_lab_rotate = TRUE,
  y_axis_log = FALSE,
  raster = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  plot <- VlnPlot_scCustom(seurat_object = seurat_object, features = "nFeature_RNA", group.by = group.by, colors_use = colors_use, pt.size = pt.size, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, plot_median = plot_median, plot_boxplot = plot_boxplot, median_size = median_size, ...) +
    geom_hline(yintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "red") +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    ggtitle(plot_title) +
    theme(plot.subtitle = element_text(hjust = 0.5), legend.position = "none")

  # Rotate x axis label
  if (isFALSE(x = x_lab_rotate)) {
    plot <- plot + UnRotate_X()
  }

  # return log10 y axis
  if (isTRUE(x = y_axis_log)) {
    plot <- plot + scale_y_log10()
  }

  return(plot)
}


#' QC Plots UMIs
#'
#' #' Custom VlnPlot for initial QC checks including lines for thresholding
#'
#' @param seurat_object Seurat object name.
#' @param plot_title Plot Title.
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident);
#' default is the current active.ident of the object.
#' @param y_axis_label Label for y axis.
#' @param x_axis_label Label for x axis.
#' @param low_cutoff Plot line a potential low threshold for filtering.
#' @param high_cutoff Plot line a potential high threshold for filtering.
#' @param pt.size Point size for plotting.
#' @param plot_median logical, whether to plot median for each ident on the plot (Default is FALSE).
#' @param median_size Shape size for the median is plotted.
#' @param plot_boxplot logical, whether to plot boxplot inside of violin (Default is FALSE).
#' @param colors_use vector of colors to use for plot.
#' @param x_lab_rotate Rotate x-axis labels 45 degrees (Default is TRUE).
#' @param y_axis_log logical. Whether to change y axis to log10 scale (Default is FALSE).
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 100,000 total points plotted (# Cells x # of features).
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ... Extra parameters passed to \code{\link[Seurat]{VlnPlot}}.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom Seurat VlnPlot
#'
#' @export
#'
#' @concept object_qc_plotting
#'
#' @examples
#' library(Seurat)
#' QC_Plots_UMIs(seurat_object = pbmc_small, plot_title = "UMIs per Cell", low_cutoff = 75,
#' high_cutoff = 600)
#'

QC_Plots_UMIs <- function(
  seurat_object,
  plot_title = "UMIs per Cell/Nucleus",
  group.by = NULL,
  x_axis_label = NULL,
  y_axis_label = "UMIs",
  low_cutoff = NULL,
  high_cutoff = NULL,
  pt.size = NULL,
  plot_median = FALSE,
  median_size = 15,
  plot_boxplot = FALSE,
  colors_use = NULL,
  x_lab_rotate = TRUE,
  y_axis_log = FALSE,
  raster = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  plot <- VlnPlot_scCustom(seurat_object = seurat_object, features = "nCount_RNA", group.by = group.by, colors_use = colors_use, pt.size = pt.size, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, plot_median = plot_median, plot_boxplot = plot_boxplot, median_size = median_size, ...) +
    geom_hline(yintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "red") +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    ggtitle(plot_title) +
    theme(plot.subtitle = element_text(hjust = 0.5), legend.position = "none")

  # Rotate x axis label
  if (isFALSE(x = x_lab_rotate)) {
    plot <- plot + UnRotate_X()
  }

  # return log10 y axis
  if (isTRUE(x = y_axis_log)) {
    plot <- plot + scale_y_log10()
  }

  return(plot)
}


#' QC Plots Mito
#'
#' #' Custom VlnPlot for initial QC checks including lines for thresholding
#'
#' @param seurat_object Seurat object name.
#' @param mito_name The column name containing percent mitochondrial counts information.  Default value is
#' "percent_mito" which is default value created when using `Add_Mito_Ribo_Seurat()`.
#' @param plot_title Plot Title.
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident);
#' default is the current active.ident of the object.
#' @param y_axis_label Label for y axis.
#' @param x_axis_label Label for x axis.
#' @param low_cutoff Plot line a potential low threshold for filtering.
#' @param high_cutoff Plot line a potential high threshold for filtering.
#' @param pt.size Point size for plotting.
#' @param plot_median logical, whether to plot median for each ident on the plot (Default is FALSE).
#' @param median_size Shape size for the median is plotted.
#' @param plot_boxplot logical, whether to plot boxplot inside of violin (Default is FALSE).
#' @param colors_use vector of colors to use for plot.
#' @param x_lab_rotate Rotate x-axis labels 45 degrees (Default is TRUE).
#' @param y_axis_log logical. Whether to change y axis to log10 scale (Default is FALSE).
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 100,000 total points plotted (# Cells x # of features).
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ... Extra parameters passed to \code{\link[Seurat]{VlnPlot}}.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom Seurat VlnPlot
#'
#' @export
#'
#' @concept object_qc_plotting
#'
#' @examples
#' \dontrun{
#' QC_Plots_Mito(seurat_object = object, plot_title = "Percent Mito per Cell", high_cutoff = 10)
#' }
#'

QC_Plots_Mito <- function(
  seurat_object,
  mito_name = "percent_mito",
  plot_title = "Mito Gene % per Cell/Nucleus",
  group.by = NULL,
  x_axis_label = NULL,
  y_axis_label = "% Mitochondrial Gene Counts",
  low_cutoff = NULL,
  high_cutoff = NULL,
  pt.size = NULL,
  plot_median = FALSE,
  median_size = 15,
  plot_boxplot = FALSE,
  colors_use = NULL,
  x_lab_rotate = TRUE,
  y_axis_log = FALSE,
  raster = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  plot <- VlnPlot_scCustom(seurat_object = seurat_object, features = mito_name, group.by = group.by, colors_use = colors_use, pt.size = pt.size, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, plot_median = plot_median, plot_boxplot = plot_boxplot, median_size = median_size, ...) +
    geom_hline(yintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "red") +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    ggtitle(plot_title) +
    theme(plot.subtitle = element_text(hjust = 0.5), legend.position = "none")

  # Rotate x axis label
  if (isFALSE(x = x_lab_rotate)) {
    plot <- plot + UnRotate_X()
  }

  # return log10 y axis
  if (isTRUE(x = y_axis_log)) {
    plot <- plot + scale_y_log10()
  }

  return(plot)
}


#' QC Plots Feature
#'
#' Custom VlnPlot for initial QC checks including lines for thresholding
#'
#' @param seurat_object Seurat object name.
#' @param feature Feature from Meta Data to plot.
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident);
#' default is the current active.ident of the object.
#' @param y_axis_label Label for y axis.
#' @param x_axis_label Label for x axis.
#' @param plot_title Plot Title.
#' @param low_cutoff Plot line a potential low threshold for filtering.
#' @param high_cutoff Plot line a potential high threshold for filtering.
#' @param pt.size Point size for plotting.
#' @param plot_median logical, whether to plot median for each ident on the plot (Default is FALSE).
#' @param median_size Shape size for the median is plotted.
#' @param plot_boxplot logical, whether to plot boxplot inside of violin (Default is FALSE).
#' @param colors_use vector of colors to use for plot.
#' @param x_lab_rotate Rotate x-axis labels 45 degrees (Default is TRUE).
#' @param y_axis_log logical. Whether to change y axis to log10 scale (Default is FALSE).
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 100,000 total points plotted (# Cells x # of features).
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ... Extra parameters passed to \code{\link[Seurat]{VlnPlot}}.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom Seurat VlnPlot
#'
#' @export
#'
#' @concept object_qc_plotting
#'
#' @examples
#' \dontrun{
#' QC_Plots_Feature(seurat_object = object, feature = "FEATURE_NAME",
#' y_axis_label = "FEATURE per Cell", plot_title = "FEATURE per Cell", high_cutoff = 10,
#' low_cutoff = 2)
#' }
#'

QC_Plots_Feature <- function(
  seurat_object,
  feature,
  group.by = NULL,
  x_axis_label = NULL,
  y_axis_label = NULL,
  plot_title = NULL,
  low_cutoff = NULL,
  high_cutoff = NULL,
  pt.size = NULL,
  plot_median = FALSE,
  median_size = 15,
  plot_boxplot = FALSE,
  colors_use = NULL,
  x_lab_rotate = TRUE,
  y_axis_log = FALSE,
  raster = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  if (is.null(x = plot_title)) {
    plot_title <- paste0(feature, " per Cell/Nucleus")
  }
  plot <- VlnPlot_scCustom(seurat_object = seurat_object, features = feature, group.by = group.by, colors_use = colors_use, pt.size = pt.size, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, plot_median = plot_median, plot_boxplot = plot_boxplot, median_size = median_size, ...) +
    geom_hline(yintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "red") +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    ggtitle(plot_title) +
    theme(plot.subtitle = element_text(hjust = 0.5), legend.position = "none")

  # Rotate x axis label
  if (isFALSE(x = x_lab_rotate)) {
    plot <- plot + UnRotate_X()
  }

  # return log10 y axis
  if (isTRUE(x = y_axis_log)) {
    plot <- plot + scale_y_log10()
  }

  return(plot)
}


#' QC Plots Cell "Complexity"
#'
#' Custom VlnPlot for initial QC checks including lines for thresholding
#'
#' @param seurat_object Seurat object name.
#' @param feature Feature from Meta Data to plot.
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident);
#' default is the current active.ident of the object.
#' @param y_axis_label Label for y axis.
#' @param x_axis_label Label for x axis.
#' @param plot_title Plot Title.
#' @param low_cutoff Plot line a potential low threshold for filtering.
#' @param high_cutoff Plot line a potential high threshold for filtering.
#' @param pt.size Point size for plotting
#' @param plot_median logical, whether to plot median for each ident on the plot (Default is FALSE).
#' @param median_size Shape size for the median is plotted.
#' @param plot_boxplot logical, whether to plot boxplot inside of violin (Default is FALSE).
#' @param colors_use vector of colors to use for plot.
#' @param x_lab_rotate Rotate x-axis labels 45 degrees (Default is TRUE).
#' @param y_axis_log logical. Whether to change y axis to log10 scale (Default is FALSE).
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 100,000 total points plotted (# Cells x # of features).
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ... Extra parameters passed to \code{\link[Seurat]{VlnPlot}}.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom Seurat VlnPlot
#'
#' @export
#'
#' @concept object_qc_plotting
#'
#' @examples
#' library(Seurat)
#' pbmc_small <- Add_Cell_Complexity_Seurat(pbmc_small)
#'
#' QC_Plots_Complexity(seurat_object = pbmc_small)
#'

QC_Plots_Complexity <- function(
  seurat_object,
  feature = "log10GenesPerUMI",
  group.by = NULL,
  x_axis_label = NULL,
  y_axis_label = "log10(Genes) / log10(UMIs)",
  plot_title = "Cell Complexity",
  low_cutoff = NULL,
  high_cutoff = NULL,
  pt.size = NULL,
  plot_median = FALSE,
  plot_boxplot = FALSE,
  median_size = 15,
  colors_use = NULL,
  x_lab_rotate = TRUE,
  y_axis_log = FALSE,
  raster = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123,
  ...
) {
  plot <- QC_Plots_Feature(seurat_object = seurat_object, feature = feature, group.by = group.by, x_axis_label = x_axis_label, y_axis_label = y_axis_label, plot_title = plot_title, low_cutoff = low_cutoff, high_cutoff = high_cutoff, pt.size = pt.size, colors_use = colors_use, x_lab_rotate = x_lab_rotate, y_axis_log = y_axis_log, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, plot_median = plot_median, median_size = median_size, plot_boxplot = plot_boxplot, ...)

  return(plot)
}


#' QC Plots Genes, UMIs, & % Mito
#'
#' Custom VlnPlot for initial QC checks including lines for thresholding
#'
#' @param seurat_object Seurat object name.
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident);
#' default is the current active.ident of the object.
#' @param feature_cutoffs Numeric vector of length 1 or 2 to plot lines for  potential low/high threshold for filtering.
#' @param UMI_cutoffs Numeric vector of length 1 or 2 to plot lines for  potential low/high threshold for filtering.
#' @param mito_cutoffs Numeric vector of length 1 or 2 to plot lines for  potential low/high threshold for filtering.
#' @param mito_name The column name containing percent mitochondrial counts information.  Default value is
#' "percent_mito" which is default value created when using `Add_Mito_Ribo_Seurat()`.
#' @param pt.size Point size for plotting
#' @param plot_median logical, whether to plot median for each ident on the plot (Default is FALSE).
#' @param median_size Shape size for the median is plotted.
#' @param plot_boxplot logical, whether to plot boxplot inside of violin (Default is FALSE).
#' @param colors_use vector of colors to use for plot.
#' @param x_lab_rotate Rotate x-axis labels 45 degrees (Default is TRUE).
#' @param y_axis_log logical. Whether to change y axis to log10 scale (Default is FALSE).
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 100,000 total points plotted (# Cells x # of features).
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ... Extra parameters passed to \code{\link[Seurat]{VlnPlot}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom Seurat VlnPlot
#' @importFrom patchwork wrap_plots
#'
#' @export
#'
#' @concept object_qc_plotting
#'
#' @examples
#' \dontrun{
#' QC_Plots_Combined_Vln(seurat_object = object)
#' }
#'

QC_Plots_Combined_Vln <- function(
  seurat_object,
  group.by = NULL,
  feature_cutoffs = NULL,
  UMI_cutoffs = NULL,
  mito_cutoffs = NULL,
  mito_name = "percent_mito",
  pt.size = NULL,
  plot_median = FALSE,
  median_size = 15,
  plot_boxplot = FALSE,
  colors_use = NULL,
  x_lab_rotate = TRUE,
  y_axis_log = FALSE,
  raster = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Setup cutoff values
  if (length(x = feature_cutoffs) > 2 || length(x = UMI_cutoffs) > 2 || length(x = mito_cutoffs) > 2) {
    cli_abort(message = "Length of each cutoff vector cannot be greater than {.field 2 (two)}.")
  }

  if (length(x = feature_cutoffs) == 1) {
    feature_cutoffs <- c(NULL, feature_cutoffs)
  }

  if (length(x = UMI_cutoffs) == 1) {
    UMI_cutoffs <- c(NULL, UMI_cutoffs)
  }

  if (length(x = mito_cutoffs) == 1) {
    mito_cutoffs <- c(NULL, mito_cutoffs)
  }

  # Create Individual Plots
  feature_plot <- QC_Plots_Genes(seurat_object = seurat_object, group.by = group.by, low_cutoff = feature_cutoffs[1], high_cutoff = feature_cutoffs[2], pt.size = pt.size, colors_use = colors_use, x_lab_rotate = x_lab_rotate, y_axis_log = y_axis_log, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, plot_median = plot_median, median_size = median_size, plot_boxplot = plot_boxplot, ...)

  UMI_plot <- QC_Plots_UMIs(seurat_object = seurat_object, group.by = group.by, low_cutoff = UMI_cutoffs[1], high_cutoff = UMI_cutoffs[2], pt.size = pt.size, colors_use = colors_use, x_lab_rotate = x_lab_rotate, y_axis_log = y_axis_log, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, plot_median = plot_median, median_size = median_size, plot_boxplot = plot_boxplot, ...)

  mito_plot <- QC_Plots_Mito(seurat_object = seurat_object, group.by = group.by, mito_name = mito_name, low_cutoff = mito_cutoffs[1], high_cutoff = mito_cutoffs[2], pt.size = pt.size, colors_use = colors_use, x_lab_rotate = x_lab_rotate, y_axis_log = y_axis_log, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, plot_median = plot_median, median_size = median_size, plot_boxplot = plot_boxplot, ...)

  # wrap plots
  plots <- wrap_plots(feature_plot, UMI_plot, mito_plot, ncol = 3)

  return(plots)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### OBJECT QC HISTOGRAM ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' QC Histogram Plots
#'
#' Custom histogram for initial QC checks including lines for thresholding
#'
#' @param seurat_object Seurat object name.
#' @param features Feature from meta.data, assay features, or feature name shortcut to plot.
#' @param low_cutoff Plot line a potential low threshold for filtering.
#' @param high_cutoff Plot line a potential high threshold for filtering.
#' @param split.by Feature to split plots by (i.e. "orig.ident").
#' @param bins number of bins to plot default is 250.
#' @param colors_use color to fill histogram bars, default is "dodgerblue".
#' @param num_columns Number of columns in plot layout.
#' @param plot_title optional, vector to use for plot title.  Default is the name of the
#' variable being plotted.
#' @param assay assay to pull features from, default is active assay.
#' @param print_defaults return list of accepted default shortcuts to provide to `features` instead
#' of full name.
#'
#' @return A patchwork object
#'
#' @import cli
#' @import ggplot2
#' @importFrom cowplot theme_cowplot
#' @importFrom dplyr filter
#' @importFrom magrittr "%>%"
#' @importFrom patchwork wrap_plots plot_annotation
#'
#' @export
#'
#' @concept object_qc_plotting
#'
#' @examples
#' \dontrun{
#' QC_Histogram(seurat_object = object, features = "nFeature_RNA")
#' }
#'

QC_Histogram <- function(
    seurat_object,
    features,
    low_cutoff = NULL,
    high_cutoff = NULL,
    split.by = NULL,
    bins = 250,
    colors_use = "dodgerblue",
    num_columns = NULL,
    plot_title = NULL,
    assay = NULL,
    print_defaults = FALSE
){
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # default features
  found_defaults <- Return_QC_Defaults(seurat_object = seurat_object, features = features, print_defaults = print_defaults)

  # set assay
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Check split valid
  if (!is.null(x = split.by)) {
    split.by <- Meta_Present(seurat_object = seurat_object, meta_col_names = split.by, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  # Check feature length if split.by provided
  if (!is.null(x = split.by)) {
    if (length(x = features) != 1) {
      cli_abort(message = "Only 1 feature can be plotted when {.code split.by = TRUE}.")
    }
  }

  # Check against object
  found_features <- Gene_Present(data = seurat_object, gene_list = found_defaults[[2]], omit_warn = FALSE, print_msg = FALSE, case_check_msg = FALSE, return_none = TRUE, seurat_assay = assay)

  found_meta <- Meta_Present(seurat_object = seurat_object, meta_col_names = found_features[[2]], omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)

  # Combine lists
  all_not_found_features <- found_meta[[2]]

  all_found_features <- c(found_defaults[[1]], found_features[[1]], found_meta[[1]])

  # Warn not found
  if (length(x = all_not_found_features > 0)) {
    cli_warn(message = c("The following features were omitted as they not found in default values or in Seurat object:",
                         "i" = "{.field {glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}}"))
  }

  # Check and set titles
  if (!is.null(x = plot_title) && length(x = plot_title) != length(x = features)) {
    cli_abort(message = "The number of {.code plot_title} (.field {length(x = plot_title)}}) does not equal number of features ({.field {length(x = all_found_features)}})")
  } else {
    plot_titles <- plot_title
  }

  if (is.null(x = plot_title) && is.null(x = split.by)) {
    plot_titles <- all_found_features
  }

  # Plot
  if (is.null(x = split.by)) {
    plot_list <- lapply(1:length(x = all_found_features), function(x) {
      plot <- ggplot(data = seurat_object@meta.data, aes(x = .data[[all_found_features[x]]])) +
        geom_histogram(color = "black", fill = colors_use, bins = bins) +
        theme_cowplot() +
        geom_vline(xintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "red") +
        ggtitle(plot_titles[x])
    })

    # wrap and return plots
    plots <- wrap_plots(plot_list, ncol = num_columns)

    return(plots)

  } else {
    # Pull required data
    data_to_plot <- FetchData(object = seurat_object, vars = c(all_found_features, split.by))

    # Extract split.by list of values
    if (inherits(x = seurat_object@meta.data[, split.by], what = "factor")) {
      meta_sample_list <- as.character(x = levels(x = seurat_object@meta.data[, split.by]))
    } else {
      meta_sample_list <- as.character(x = unique(x = seurat_object@meta.data[, split.by]))
    }

    if (length(x = colors_use) != length(x = meta_sample_list)) {
      if (length(x = colors_use == 1)) {
        if (colors_use == "dodgerblue") {
          colors_use <- scCustomize_Palette(num_groups = length(x = meta_sample_list))
        }
      } else {
        cli_abort(message = c("The number of colors must match the number of variables in {.code split.by}.",
                              "i" = "The length of {.code colors_use} is {.field {length(x = colors_use)}} but the number of variables in {.code spliut.by} is {.field {length(x = split.by)}}"))
      }
    }

    # Plot
    plot_list <- lapply(1:length(x = meta_sample_list), function(x) {
      sub_data <- data_to_plot %>%
        filter(.data[[split.by]] == meta_sample_list[x])

      plot <- ggplot(data = sub_data, aes(x = .data[[all_found_features]])) +
        geom_histogram(color = "black", fill = colors_use[x], bins = bins) +
        theme_cowplot() +
        geom_vline(xintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "red") +
        ggtitle(meta_sample_list[x])
    })

    # wrap and return plots
    plots <- wrap_plots(plot_list, ncol = num_columns) + plot_annotation(title = all_found_features, theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5))))

    return(plots)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### OBJECT QC SCATTER ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' QC Plots Genes vs UMIs
#'
#' Custom FeatureScatter for initial QC checks including lines for thresholding
#'
#' @param seurat_object Seurat object name.
#' @param y_axis_label Label for y axis.
#' @param x_axis_label Label for x axis.
#' @param low_cutoff_gene Plot line a potential low threshold for filtering genes per cell.
#' @param high_cutoff_gene Plot line a potential high threshold for filtering genes per cell.
#' @param low_cutoff_UMI Plot line a potential low threshold for filtering UMIs per cell.
#' @param high_cutoff_UMI Plot line a potential high threshold for filtering UMIs per cell.
#' @param colors_use vector of colors to use for plotting by identity.
#' @param meta_gradient_name Name of continuous meta data variable to color points in plot by.
#' (MUST be continuous variable i.e. "percent_mito").
#' @param meta_gradient_color The gradient color palette to use for plotting of meta variable (default is
#' viridis "Plasma" palette with dark colors high).
#' @param meta_gradient_low_cutoff Value to use as threshold for plotting.  `meta_gradient_name` values
#' below this value will be plotted using `meta_gradient_na_color`.
#' @param meta_gradient_na_color Color to use for plotting values when a `meta_gradient_low_cutoff` is
#' set (default is "lightgray").
#' @param cells Cells to include on the scatter plot (default is all cells).
#' @param combination logical (default FALSE).  Whether or not to return a plot layout with both the
#' plot colored by identity and the meta data gradient plot.
#' @param ident_legend logical, whether to plot the legend containing identities (left plot) when
#' `combination = TRUE`.  Default is TRUE.
#' @param pt.size Passes size of points to both \code{\link[Seurat]{FeatureScatter}} and `geom_point`.
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident).
#' Default is `@active.ident`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 100,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed Random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param shuffle_seed Sets the seed if randomly shuffling the order of points (Default is 1).
#' @param ... Extra parameters passed to \code{\link[Seurat]{FeatureScatter}}.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom cowplot theme_cowplot
#' @importFrom dplyr filter arrange
#' @importFrom magrittr "%>%"
#' @importFrom scattermore geom_scattermore
#' @importFrom Seurat FeatureScatter
#' @importFrom stats cor
#'
#' @export
#'
#' @concept object_qc_plotting
#'
#' @examples
#' library(Seurat)
#' QC_Plot_UMIvsGene(seurat_object = pbmc_small, x_axis_label = "UMIs per Cell/Nucleus",
#' y_axis_label = "Genes per Cell/Nucleus")
#'

QC_Plot_UMIvsGene <- function(
  seurat_object,
  x_axis_label = "UMIs per Cell/Nucleus",
  y_axis_label = "Genes per Cell/Nucleus",
  low_cutoff_gene = -Inf,
  high_cutoff_gene = Inf,
  low_cutoff_UMI = -Inf,
  high_cutoff_UMI = Inf,
  colors_use = NULL,
  meta_gradient_name = NULL,
  meta_gradient_color = viridis_plasma_dark_high,
  meta_gradient_na_color = "lightgray",
  meta_gradient_low_cutoff = NULL,
  cells = NULL,
  combination = FALSE,
  ident_legend = TRUE,
  pt.size = 1,
  group.by = NULL,
  raster = NULL,
  raster.dpi = c(512, 512),
  ggplot_default_colors = FALSE,
  color_seed = 123,
  shuffle_seed = 1,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Default raster check
  if (isTRUE(x = combination)) {
    raster <- raster %||% (length(x = Cells(x = seurat_object)) > 1e5)
  } else {
    raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)
  }

  # select color palette if not specified
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }
  if (is.null(x = colors_use)) {
    if (isTRUE(x = ggplot_default_colors)) {
      colors_use <- Hue_Pal(group_by_length)
    } else {
      if (group_by_length <= 2) {
        colors_use <- NavyAndOrange()
      }
      if (group_by_length > 2 && group_by_length <= 4) {
        colors_use <- JCO_Four()
      }
      if (group_by_length > 4 && group_by_length <= 36) {
        colors_use <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
      }
      if (group_by_length > 36) {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "varibow", shuffle_pal = TRUE, seed = color_seed)
      }
    }
  }

  if (isFALSE(x = ident_legend) && isFALSE(x = combination)) {
    cli_warn(message = "{.code ident_legend} parameter ignored as {.code combination = FALSE}")
  }


  # Pull meta data
  featurescatter_data <- Fetch_Meta(object = seurat_object) %>%
    rownames_to_column("barcodes")
  # Check valid meta variable
  if (!is.null(x = meta_gradient_name)) {
    meta_names <- colnames(x = featurescatter_data)
    if (meta_gradient_name %in% meta_names == FALSE) {
      cli_abort(message = "The meta data variable {.val {meta_gradient_name}} could not be found in object@metadata.")
    }
  }

  # filter meta data by cells plotting if valid
  if (!is.null(x = cells)) {
    featurescatter_data <- featurescatter_data %>%
      filter(.data[["barcodes"]] %in% cells)
  }

  # Add reported value for plot subtitle
  if (is.null(x = meta_gradient_low_cutoff)) {
    meta_cutoff_reported <- 0
  } else {
    meta_cutoff_reported <- meta_gradient_low_cutoff
  }

  # Sort meta data by chosen variable
  if (!is.null(x = meta_gradient_name)) {
    featurescatter_data_sort <- featurescatter_data %>%
      arrange(.data[[meta_gradient_name]])
    if (is.null(x = meta_gradient_low_cutoff)) {
      meta_gradient_low_cutoff <- min(featurescatter_data_sort[[meta_gradient_name]])
    }
  } else {
    featurescatter_data_sort <- featurescatter_data
  }

  # Calculate Correlation for all data
  plot_cor_full <- round(x = cor(x = featurescatter_data_sort[, "nCount_RNA"], y = featurescatter_data_sort[, "nFeature_RNA"]), digits = 2)

  if (is.null(x = meta_gradient_name)) {
    featurescatter_data_sort_filter <- featurescatter_data_sort %>%
      filter(.data[["nCount_RNA"]] > low_cutoff_UMI & .data[["nCount_RNA"]] < high_cutoff_UMI & .data[["nFeature_RNA"]] > low_cutoff_gene & .data[["nFeature_RNA"]] < high_cutoff_gene)
  } else {
    featurescatter_data_sort_filter <- featurescatter_data_sort %>%
      filter(.data[["nCount_RNA"]] > low_cutoff_UMI & .data[["nCount_RNA"]] < high_cutoff_UMI & .data[["nFeature_RNA"]] > low_cutoff_gene & .data[["nFeature_RNA"]] < high_cutoff_gene & .data[[meta_gradient_name]] < meta_gradient_low_cutoff)
  }

  # Calculate correlation based on cutoffs
  plot_cor_filtered <- round(x = cor(x = featurescatter_data_sort_filter[, "nCount_RNA"], y = featurescatter_data_sort_filter[, "nFeature_RNA"]), digits = 2)

  # Plot with meta gradient
  if (!is.null(x = meta_gradient_name) && isFALSE(x = combination)) {
    if (isTRUE(x = raster)) {
      p1 <- ggplot(data = featurescatter_data_sort, mapping = aes(x = .data[["nCount_RNA"]], y = .data[["nFeature_RNA"]])) +
        geom_scattermore(mapping = aes(color = .data[[meta_gradient_name]]), pointsize = pt.size) +
        scale_color_gradientn(colors = meta_gradient_color, limits = c(meta_gradient_low_cutoff, NA), na.value = meta_gradient_na_color) +
        theme_cowplot() +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_hline(yintercept = c(if(is.finite(x = low_cutoff_gene)) {low_cutoff_gene}, if(is.finite(x = high_cutoff_gene)) {high_cutoff_gene}), linetype = "dashed", color = "red") +
        geom_vline(xintercept = c(if(is.finite(x = low_cutoff_UMI)) {low_cutoff_UMI}, if(is.finite(x = high_cutoff_UMI)) {high_cutoff_UMI}), linetype = "dashed", color = "blue") +
        xlab(x_axis_label) +
        ylab(y_axis_label) +
        ggtitle("Genes vs. UMIs per Cell/Nucleus", subtitle = c(paste0("Correlation of full dataset is: ", plot_cor_full, ".", "\nCorrelation of filtered dataset would be: ", plot_cor_filtered, ".  ", "\nThe low cutoff for plotting ", meta_gradient_name, " is: ", meta_cutoff_reported)))
      return(p1)
    }
    p1 <- ggplot(data = featurescatter_data_sort, mapping = aes(x = .data[["nCount_RNA"]], y = .data[["nFeature_RNA"]])) +
      geom_point(mapping = aes(color = .data[[meta_gradient_name]]), size = pt.size) +
      scale_color_gradientn(colors = meta_gradient_color, limits = c(meta_gradient_low_cutoff, NA), na.value = meta_gradient_na_color) +
      theme_cowplot() +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_hline(yintercept = c(if(is.finite(x = low_cutoff_gene)) {low_cutoff_gene}, if(is.finite(x = high_cutoff_gene)) {high_cutoff_gene}), linetype = "dashed", color = "red") +
      geom_vline(xintercept = c(if(is.finite(x = low_cutoff_UMI)) {low_cutoff_UMI}, if(is.finite(x = high_cutoff_UMI)) {high_cutoff_UMI}), linetype = "dashed", color = "blue") +
      xlab(x_axis_label) +
      ylab(y_axis_label) +
      ggtitle("Genes vs. UMIs per Cell/Nucleus", subtitle = c(paste0("Correlation of full dataset is: ", plot_cor_full, ".", "\nCorrelation of filtered dataset would be: ", plot_cor_filtered, ".  ", "\nThe low cutoff for plotting ", meta_gradient_name, " is: ", meta_cutoff_reported)))
    return(p1)
  }
  # Plot by identity
  if (is.null(x = meta_gradient_name) && isFALSE(x = combination)) {
    p1 <- FeatureScatter(object = seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cells = cells, pt.size = pt.size, shuffle = TRUE,  raster = raster, raster.dpi = raster.dpi, cols = colors_use, group.by = group.by, seed = shuffle_seed, ...) +
      geom_hline(yintercept = c(if(is.finite(x = low_cutoff_gene)) {low_cutoff_gene}, if(is.finite(x = high_cutoff_gene)) {high_cutoff_gene}), linetype = "dashed", color = "red") +
      geom_vline(xintercept = c(if(is.finite(x = low_cutoff_UMI)) {low_cutoff_UMI}, if(is.finite(x = high_cutoff_UMI)) {high_cutoff_UMI}), linetype = "dashed", color = "blue") +
      xlab(x_axis_label) +
      ylab(y_axis_label) +
      ggtitle("Genes vs. UMIs per Cell/Nucleus", subtitle = c(paste0("Correlation of full dataset is: ", plot_cor_full, ".", "\nCorrelation of filtered dataset would be: ", plot_cor_filtered, ".")))
    return(p1)
  }

  if (isTRUE(x = combination)) {
    # Plot by identity
    p1 <- FeatureScatter(object = seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cells = cells, pt.size = pt.size, shuffle = TRUE, raster = raster, raster.dpi = raster.dpi, cols = colors_use, group.by = group.by, seed = shuffle_seed, ...) +
      geom_hline(yintercept = c(if(is.finite(x = low_cutoff_gene)) {low_cutoff_gene}, if(is.finite(x = high_cutoff_gene)) {high_cutoff_gene}), linetype = "dashed", color = "red") +
      geom_vline(xintercept = c(if(is.finite(x = low_cutoff_UMI)) {low_cutoff_UMI}, if(is.finite(x = high_cutoff_UMI)) {high_cutoff_UMI}), linetype = "dashed", color = "blue") +
      xlab(x_axis_label) +
      ylab(y_axis_label) + ggtitle("")

    if (isFALSE(x = ident_legend)) {
      p1 <- p1 + NoLegend()
    }

    # Plot with meta gradient
    if (isTRUE(x = raster)) {
      p2 <- ggplot(data = featurescatter_data_sort, mapping = aes(x = .data[["nCount_RNA"]], y = .data[["nFeature_RNA"]])) +
        geom_scattermore(mapping = aes(color = .data[[meta_gradient_name]]), pointsize = pt.size) +
        scale_color_gradientn(colors = meta_gradient_color, limits = c(meta_gradient_low_cutoff, NA), na.value = meta_gradient_na_color) +
        theme_cowplot() +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_hline(yintercept = c(if(is.finite(x = low_cutoff_gene)) {low_cutoff_gene}, if(is.finite(x = high_cutoff_gene)) {high_cutoff_gene}), linetype = "dashed", color = "red") +
        geom_vline(xintercept = c(if(is.finite(x = low_cutoff_UMI)) {low_cutoff_UMI}, if(is.finite(x = high_cutoff_UMI)) {high_cutoff_UMI}), linetype = "dashed", color = "blue") +
        xlab(x_axis_label) +
        ylab(y_axis_label)
    } else {
      p2 <- ggplot(data = featurescatter_data_sort, mapping = aes(x = .data[["nCount_RNA"]], y = .data[["nFeature_RNA"]])) +
        geom_point(mapping = aes(color = .data[[meta_gradient_name]]), size = pt.size) +
        scale_color_gradientn(colors = meta_gradient_color, limits = c(meta_gradient_low_cutoff, NA), na.value = meta_gradient_na_color) +
        theme_cowplot() +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_hline(yintercept = c(if(is.finite(x = low_cutoff_gene)) {low_cutoff_gene}, if(is.finite(x = high_cutoff_gene)) {high_cutoff_gene}), linetype = "dashed", color = "red") +
        geom_vline(xintercept = c(if(is.finite(x = low_cutoff_UMI)) {low_cutoff_UMI}, if(is.finite(x = high_cutoff_UMI)) {high_cutoff_UMI}), linetype = "dashed", color = "blue") +
        xlab(x_axis_label) +
        ylab(y_axis_label)
    }

    # plot together in layout
    plots <- wrap_plots(p1, p2) + plot_annotation(title = "Genes vs. UMIs per Cell/Nucleus", subtitle = c(paste0("Correlation of full dataset is: ", plot_cor_full, ".", "\nCorrelation of filtered dataset would be: ", plot_cor_filtered, ".  ", "\nThe low cutoff for plotting ", meta_gradient_name, " is: ", meta_cutoff_reported))) & theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    return(plots)
  }
}


#' QC Plots Genes vs Misc
#'
#' Custom FeatureScatter for initial QC checks including lines for thresholding
#'
#' @param seurat_object Seurat object name.
#' @param feature1 First feature to plot.
#' @param y_axis_label Label for y axis.
#' @param x_axis_label Label for x axis.
#' @param low_cutoff_gene Plot line a potential low threshold for filtering genes per cell.
#' @param high_cutoff_gene Plot line a potential high threshold for filtering genes per cell.
#' @param low_cutoff_feature Plot line a potential low threshold for filtering feature1 per cell.
#' @param high_cutoff_feature Plot line a potential high threshold for filtering feature1 per cell.
#' @param colors_use vector of colors to use for plotting by identity.
#' @param pt.size Adjust point size for plotting.
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident).
#'   Default is `@active.ident`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if greater
#' than 100,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using default
#' ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param shuffle_seed Sets the seed if randomly shuffling the order of points (Default is 1).
#' @param ... Extra parameters passed to \code{\link[Seurat]{FeatureScatter}}.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom Seurat FeatureScatter
#'
#' @export
#'
#' @concept object_qc_plotting
#'
#' @examples
#' \dontrun{
#' QC_Plot_GenevsFeature(seurat_object = obj, y_axis_label = "Feature per Cell")
#' }
#'

QC_Plot_GenevsFeature <- function(
  seurat_object,
  feature1,
  x_axis_label = NULL,
  y_axis_label = "Genes per Cell/Nucleus",
  low_cutoff_gene = NULL,
  high_cutoff_gene = NULL,
  low_cutoff_feature = NULL,
  high_cutoff_feature = NULL,
  colors_use = NULL,
  pt.size = 1,
  group.by = NULL,
  raster = NULL,
  raster.dpi = c(512, 512),
  ggplot_default_colors = FALSE,
  color_seed = 123,
  shuffle_seed = 1,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Label x axis
  if (is.null(x = x_axis_label)) {
    x_axis_label <- paste0(feature1, " per Cell/Nucleus")
  }

  # select color palette if not specified
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }
  if (is.null(x = colors_use)) {
    if (isTRUE(x = ggplot_default_colors)) {
      colors_use <- Hue_Pal(group_by_length)
    } else {
      if (group_by_length <= 2) {
        colors_use <- NavyAndOrange()
      }
      if (group_by_length > 2 && group_by_length <= 4) {
        colors_use <- JCO_Four()
      }
      if (group_by_length > 4 && group_by_length <= 36) {
        colors_use <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
      }
      if (group_by_length > 36) {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "varibow", shuffle_pal = TRUE, seed = color_seed)
      }
    }
  }

  # Plot
  FeatureScatter(object = seurat_object, feature1 = feature1, feature2 = "nFeature_RNA", pt.size = pt.size, shuffle = TRUE, raster = raster, raster.dpi = raster.dpi, cols = colors_use, group.by = group.by, seed = shuffle_seed, ...) +
    geom_hline(yintercept = c(low_cutoff_gene, high_cutoff_gene), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(low_cutoff_feature, high_cutoff_feature), linetype = "dashed", color = "blue") +
    xlab(x_axis_label) +
    ylab(y_axis_label)
}


#' QC Plots UMI vs Misc
#'
#' Custom FeatureScatter for initial QC checks including lines for thresholding
#'
#' @param seurat_object Seurat object name.
#' @param feature1 First feature to plot.
#' @param y_axis_label Label for y axis.
#' @param x_axis_label Label for x axis.
#' @param low_cutoff_UMI Plot line a potential low threshold for filtering UMI per cell.
#' @param high_cutoff_UMI Plot line a potential high threshold for filtering UMI per cell.
#' @param low_cutoff_feature Plot line a potential low threshold for filtering feature1 per cell.
#' @param high_cutoff_feature Plot line a potential high threshold for filtering feature1 per cell.
#' @param colors_use vector of colors to use for plotting by identity.
#' @param pt.size Adjust point size for plotting.
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident).
#' Default is `@active.ident`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if greater
#' than 100,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param shuffle_seed Sets the seed if randomly shuffling the order of points (Default is 1).
#' @param ... Extra parameters passed to \code{\link[Seurat]{FeatureScatter}}.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom Seurat FeatureScatter
#'
#' @export
#'
#' @concept object_qc_plotting
#'
#' @examples
#' \dontrun{
#' QC_Plot_UMIvsFeature(seurat_object = obj, y_axis_label = "Feature per Cell")
#' }
#'

QC_Plot_UMIvsFeature <- function(
  seurat_object,
  feature1,
  x_axis_label = NULL,
  y_axis_label = "UMIs per Cell/Nucleus",
  low_cutoff_UMI = NULL,
  high_cutoff_UMI = NULL,
  low_cutoff_feature = NULL,
  high_cutoff_feature = NULL,
  colors_use = NULL,
  pt.size = 1,
  group.by = NULL,
  raster = NULL,
  raster.dpi = c(512, 512),
  ggplot_default_colors = FALSE,
  color_seed = 123,
  shuffle_seed = 1,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Label x axis
  if (is.null(x = x_axis_label)) {
    x_axis_label <- paste0(feature1, " per Cell/Nucleus")
  }

  # select color palette if not specified
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }
  if (is.null(x = colors_use)) {
    if (isTRUE(x = ggplot_default_colors)) {
      colors_use <- Hue_Pal(group_by_length)
    } else {
      if (group_by_length <= 2) {
        colors_use <- NavyAndOrange()
      }
      if (group_by_length > 2 && group_by_length <= 4) {
        colors_use <- JCO_Four()
      }
      if (group_by_length > 4 && group_by_length <= 36) {
        colors_use <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
      }
      if (group_by_length > 36) {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "varibow", shuffle_pal = TRUE, seed = color_seed)
      }
    }
  }

  # Plot
  FeatureScatter(object = seurat_object, feature1 = feature1, feature2 = "nCount_RNA", pt.size = pt.size, shuffle = TRUE, raster = raster, raster.dpi = raster.dpi, cols = colors_use, group.by = group.by, seed = shuffle_seed, ...) +
    geom_hline(yintercept = c(low_cutoff_UMI, high_cutoff_UMI), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(low_cutoff_feature, high_cutoff_feature), linetype = "dashed", color = "blue") +
    xlab(x_axis_label) +
    ylab(y_axis_label)
}
