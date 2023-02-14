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
#' @param pt.size Point size for plotting
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

  # Add pt.size check
  pt.size <- pt.size %||% AutoPointSize_scCustom(data = seurat_object)

  plot <- VlnPlot_scCustom(seurat_object = seurat_object, features = "nFeature_RNA", group.by = group.by, colors_use = colors_use, pt.size = pt.size, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, ...) +
    geom_hline(yintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "red") +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    ggtitle(plot_title) +
    theme(plot.subtitle = element_text(hjust = 0.5), legend.position = "none")

  # Rotate x axis label
  if (!x_lab_rotate) {
    plot <- plot + UnRotate_X()
  }

  # return log10 y axis
  if (y_axis_log) {
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
#' @param pt.size Point size for plotting
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

  # Add pt.size check
  pt.size <- pt.size %||% AutoPointSize_scCustom(data = seurat_object)

  plot <- VlnPlot_scCustom(seurat_object = seurat_object, features = "nCount_RNA", group.by = group.by, colors_use = colors_use, pt.size = pt.size, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, ...) +
    geom_hline(yintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "red") +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    ggtitle(plot_title) +
    theme(plot.subtitle = element_text(hjust = 0.5), legend.position = "none")

  # Rotate x axis label
  if (!x_lab_rotate) {
    plot <- plot + UnRotate_X()
  }

  # return log10 y axis
  if (y_axis_log) {
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
#' @param pt.size Point size for plotting
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

  # Add pt.size check
  pt.size <- pt.size %||% AutoPointSize_scCustom(data = seurat_object)

  plot <- VlnPlot_scCustom(seurat_object = seurat_object, features = mito_name, group.by = group.by, colors_use = colors_use, pt.size = pt.size, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, ...) +
    geom_hline(yintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "red") +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    ggtitle(plot_title) +
    theme(plot.subtitle = element_text(hjust = 0.5), legend.position = "none")

  # Rotate x axis label
  if (!x_lab_rotate) {
    plot <- plot + UnRotate_X()
  }

  # return log10 y axis
  if (y_axis_log) {
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
#' @param pt.size Point size for plotting
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

  # Add pt.size check
  pt.size <- pt.size %||% AutoPointSize_scCustom(data = seurat_object)

  if (is.null(x = plot_title)) {
    plot_title <- paste0(feature, " per Cell/Nucleus")
  }
  plot <- VlnPlot_scCustom(seurat_object = seurat_object, features = mito_name, group.by = group.by, colors_use = colors_use, pt.size = pt.size, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, ...) +
    geom_hline(yintercept = c(low_cutoff, high_cutoff), linetype = "dashed", color = "red") +
    xlab(x_axis_label) +
    ylab(y_axis_label) +
    ggtitle(plot_title) +
    theme(plot.subtitle = element_text(hjust = 0.5), legend.position = "none")

  # Rotate x axis label
  if (!x_lab_rotate) {
    plot <- plot + UnRotate_X()
  }

  # return log10 y axis
  if (y_axis_log) {
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
  colors_use = NULL,
  x_lab_rotate = TRUE,
  y_axis_log = FALSE,
  raster = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123,
  ...
) {
  QC_Plots_Feature(seurat_object = seurat_object, feature = feature, group.by = group.by, x_axis_label = x_axis_label, y_axis_label = y_axis_label, plot_title = plot_title, low_cutoff = low_cutoff, high_cutoff = high_cutoff, pt.size = pt.size, colors_use = colors_use, x_lab_rotate = x_lab_rotate, y_axis_log = y_axis_log, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, ...)
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

  # Add pt.size check
  pt.size <- pt.size %||% AutoPointSize_scCustom(data = seurat_object)

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
  feature_plot <- QC_Plots_Genes(seurat_object = seurat_object, group.by = group.by, low_cutoff = feature_cutoffs[1], high_cutoff = feature_cutoffs[2], pt.size = pt.size, colors_use = colors_use, x_lab_rotate = x_lab_rotate, y_axis_log = y_axis_log, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, ...)

  UMI_plot <- QC_Plots_UMIs(seurat_object = seurat_object, group.by = group.by, low_cutoff = UMI_cutoffs[1], high_cutoff = UMI_cutoffs[2], pt.size = pt.size, colors_use = colors_use, x_lab_rotate = x_lab_rotate, y_axis_log = y_axis_log, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, ...)

  mito_plot <- QC_Plots_Mito(seurat_object = seurat_object, group.by = group.by, mito_name = mito_name, low_cutoff = mito_cutoffs[1], high_cutoff = mito_cutoffs[2], pt.size = pt.size, colors_use = colors_use, x_lab_rotate = x_lab_rotate, y_axis_log = y_axis_log, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, ...)

  # wrap plots
  plots <- wrap_plots(feature_plot, UMI_plot, mito_plot, ncol = 3)

  return(plots)
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
  if (combination) {
    raster <- raster %||% (length(x = colnames(x = seurat_object)) > 1e5)
  } else {
    raster <- raster %||% (length(x = colnames(x = seurat_object)) > 2e5)
  }

  # select color palette if not specified
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }
  if (is.null(x = colors_use)) {
    if (ggplot_default_colors) {
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

  # Pull meta data
  featurescatter_data <- Fetch_Meta(object = seurat_object) %>%
    rownames_to_column("barcodes")
  # Check valid meta variable
  if (!is.null(x = meta_gradient_name)) {
    meta_names <- colnames(featurescatter_data)
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
  if (!is.null(x = meta_gradient_name) && combination == FALSE) {
    if (raster) {
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
  if (is.null(x = meta_gradient_name) && combination == FALSE) {
    p1 <- FeatureScatter(object = seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cells = cells, pt.size = pt.size, shuffle = TRUE,  raster = raster, raster.dpi = raster.dpi, cols = colors_use, group.by = group.by, seed = shuffle_seed, ...) +
      geom_hline(yintercept = c(if(is.finite(x = low_cutoff_gene)) {low_cutoff_gene}, if(is.finite(x = high_cutoff_gene)) {high_cutoff_gene}), linetype = "dashed", color = "red") +
      geom_vline(xintercept = c(if(is.finite(x = low_cutoff_UMI)) {low_cutoff_UMI}, if(is.finite(x = high_cutoff_UMI)) {high_cutoff_UMI}), linetype = "dashed", color = "blue") +
      xlab(x_axis_label) +
      ylab(y_axis_label) +
      ggtitle("Genes vs. UMIs per Cell/Nucleus", subtitle = c(paste0("Correlation of full dataset is: ", plot_cor_full, ".", "\nCorrelation of filtered dataset would be: ", plot_cor_filtered, ".")))
    return(p1)
  }

  if (combination) {
    # Plot by identity
    p1 <- FeatureScatter(object = seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cells = cells, pt.size = pt.size, shuffle = TRUE, raster = raster, raster.dpi = raster.dpi, cols = colors_use, group.by = group.by, seed = shuffle_seed, ...) +
      geom_hline(yintercept = c(if(is.finite(x = low_cutoff_gene)) {low_cutoff_gene}, if(is.finite(x = high_cutoff_gene)) {high_cutoff_gene}), linetype = "dashed", color = "red") +
      geom_vline(xintercept = c(if(is.finite(x = low_cutoff_UMI)) {low_cutoff_UMI}, if(is.finite(x = high_cutoff_UMI)) {high_cutoff_UMI}), linetype = "dashed", color = "blue") +
      xlab(x_axis_label) +
      ylab(y_axis_label) + ggtitle("")

    # Plot with meta gradient
    if (raster) {
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
    if (ggplot_default_colors) {
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
    if (ggplot_default_colors) {
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
