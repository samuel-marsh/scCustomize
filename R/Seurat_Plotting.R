#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GENE EXPRESSION PLOTTING (2D) ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Customize FeaturePlot
#'
#' Create Custom FeaturePlots and preserve scale (no binning)
#'
#' @param seurat_object Seurat object name.
#' @param features Feature(s) to plot.
#' @param colors_use list of colors or color palette to use.
#' @param na_color color to use for points below lower limit.
#' @param order whether to move positive cells to the top (default = TRUE).
#' @param pt.size Adjust point size for plotting.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default).
#' @param na_cutoff Value to use as minimum expression cutoff.  This will be lowest value plotted use
#' palette provided to `colors_use`.  Leave as default value to plot only positive non-zero values using
#' color scale and zero/negative values as NA.  To plot all values using color palette set to `NA`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param split.by Variable in `@meta.data` to split the plot by.
#' @param split_collect logical, whether to collect the legends/guides when plotting with `split.by`.
#' Default is TRUE if one value is provided to `features` otherwise is set to FALSE.
#' @param aspect_ratio Control the aspect ratio (y:x axes ratio length).  Must be numeric value;
#' Default is NULL.
#' @param figure_plot logical.  Whether to remove the axes and plot with legend on left of plot denoting
#' axes labels.  (Default is FALSE).  Requires `split_seurat = TRUE`.
#' @param num_columns Number of columns in plot layout.
#' @param slot `r lifecycle::badge("deprecated")` soft-deprecated.  See `layer`
#' @param layer Which layer to pull expression data from?  Default is "data".
#' @param alpha_exp new alpha level to apply to expressing cell color palette (`colors_use`).  Must be
#' value between 0-1.
#' @param alpha_na_exp new alpha level to apply to non-expressing cell color palette (`na_color`).  Must be
#' value between 0-1.
#' @param label logical, whether to label the clusters.  Default is FALSE.
#' @param label_feature_yaxis logical, whether to place feature labels on secondary y-axis as opposed to
#' above legend key.  Default is FALSE.  When setting `label_feature_yaxis = TRUE` the number of columns
#' in plot output will automatically be set to the number of levels in `split.by'`
#' @param combine Combine plots into a single \code{\link[patchwork]{patchwork}ed} ggplot object.
#' If FALSE, return a list of ggplot objects.
#' @param ... Extra parameters passed to \code{\link[Seurat]{FeaturePlot}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @import patchwork
#' @importFrom methods hasArg
#' @importFrom scales alpha
#' @importFrom Seurat FeaturePlot
#' @importFrom SeuratObject DefaultDimReduc
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' library(Seurat)
#' FeaturePlot_scCustom(seurat_object = pbmc_small, features = "CD3E",
#' colors_use = viridis_plasma_dark_high, na_color = "lightgray")
#'

FeaturePlot_scCustom <- function(
  seurat_object,
  features,
  colors_use = viridis_plasma_dark_high,
  na_color = "lightgray",
  order = TRUE,
  pt.size = NULL,
  reduction = NULL,
  na_cutoff = 0.000000001,
  raster = NULL,
  raster.dpi = c(512, 512),
  split.by = NULL,
  split_collect = NULL,
  aspect_ratio = NULL,
  figure_plot = FALSE,
  num_columns = NULL,
  slot = deprecated(),
  layer = "data",
  alpha_exp = NULL,
  alpha_na_exp = NULL,
  label = FALSE,
  label_feature_yaxis = FALSE,
  combine = TRUE,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check is slot is supplied
  if (lifecycle::is_present(slot)) {
    lifecycle::deprecate_warn(when = "2.0.0",
                              what = "slot",
                              with = "layer",
                              details = c("v" = "As of Seurat 5.0.0 the {.code slot} parameter is deprecated and replaced with {.code layer}.",
                                          "i" = "Please adjust code now to prepare for full deprecation.")
    )
    layer <- slot
  }

  # Check meta
  if (!is.null(x = split.by)) {
    split.by <- Meta_Present(seurat_object = seurat_object, meta_col_names = split.by, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  # Set or check split_collect values
  if (is.null(x = split_collect)) {
    if (length(x = features) == 1) {
      split_collect <- TRUE
    } else {
      split_collect <- FALSE
    }
  }

  if (!is.null(x = split_collect)) {
    if (length(x = features) > 1 && isTRUE(x = split_collect)) {
      cli_abort(message = "{.code split_collect} cannot be set to {.field TRUE} if the number of features is greater than 1.")
    }
  }

  # Check features and meta to determine which features present
  all_found_features <- Feature_PreCheck(object = seurat_object, features = features)

  # Get length of meta data feature
  if (is.null(x = split.by) && isTRUE(x = label_feature_yaxis)) {
    cli_abort(message = "Setting {.code label_feature_yaxis = TRUE} is only supported when also setting {.code split.by}.")
  }

  if (!is.null(x = split.by)) {
    split.by_length <- length(x = unique(x = seurat_object@meta.data[[split.by]]))

    if (!is.null(x = num_columns) && isTRUE(x = label_feature_yaxis)) {

      cli_warn(message = c("Setting number of columns is not permitted if {.code label_feature_yaxis = TRUE}",
                           "i" = "Number of columns be automatically set to number of levels in `split.by` ({.field {split.by_length}}).")
      )
      num_columns <- split.by_length
    }

    if (is.null(x = num_columns)) {
      num_columns <- split.by_length
    }

    # Calculate number of rows for selected number of columns
    num_rows <- ceiling(split.by_length/num_columns)

    # Check column and row compatibility
    if (num_columns > split.by_length) {
      cli_abort(message = c("The number of columns specified is greater than the number of meta data variables.",
                        "*" = "{.val {split.by}} only contains {.field {split.by_length}} variables.",
                        "i" = "Please adjust {.code num_columns} to be less than or equal to {.field {split.by_length}}.")
      )
    }
  }

  if (any(all_found_features %in% colnames(x = seurat_object@meta.data))) {
    cli_warn(message = c("Some of the plotted features are from meta.data slot.",
                         "*" = "Please check that {.code na_cutoff} param is being set appropriately for those features.")
    )
  }

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)

  # Set uniform poist size is pt.size = NULL (based on plot with most cells)
  if (is.null(x = pt.size)) {
    if (!is.null(x = split.by)) {
      # cells per meta data
      cells_by_meta <- data.frame(table(seurat_object@meta.data[, split.by]))
      # Identity with greatest number of cells
      max_cells <- max(cells_by_meta$Freq)
      # modified version of the autopointsize function from Seurat
      pt.size <- AutoPointSize_scCustom(data = max_cells, raster = raster)
    } else {
      # Total cells
      cells_total <- nrow(x = seurat_object@meta.data)
      # modified version of the autopointsize function from Seurat
      pt.size <- AutoPointSize_scCustom(data = cells_total, raster = raster)
    }
  }

  # set na_cutoff if NULL is provided to ensure proper plotting
  if (is.null(x = na_cutoff)) {
    na_cutoff <- NA
  }

  # Extract default reduction
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)

  # Get Seurat version
  seurat_version <- packageVersion("Seurat")

  # Add alpha to color scales
  if (!is.null(x = alpha_exp) && seurat_version < 5) {
    colors_use <- alpha(colors_use, alpha_exp)
  }

  if (!is.null(x = alpha_na_exp) && seurat_version < 5) {
    na_color <- alpha(na_color, alpha_na_exp)
  }

  if (!is.null(x = alpha_na_exp) && seurat_version >= 5) {
    cli_warn(message = "{.code alpha_na_exp} is not currently supported for Seurat v5+")
  }

  # Set alpha if NULL
  if (is.null(x = alpha_exp) && seurat_version >= 5) {
      alpha_exp <- 1
  }

  # plot no split & combined
  if (is.null(x = split.by) && isTRUE(x = combine)) {
    # Keep until Seurat version required is > 5
    if (seurat_version >= 5) {
      plot <- suppressMessages(FeaturePlot(object = seurat_object, features = all_found_features, order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, ncol = num_columns, combine = combine, raster.dpi = raster.dpi, label = label, alpha = alpha_exp, ...) & scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, NA), na.value = na_color))
    } else {
      plot <- suppressMessages(FeaturePlot(object = seurat_object, features = all_found_features, order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, ncol = num_columns, combine = combine, raster.dpi = raster.dpi, label = label, ...) & scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, NA), na.value = na_color))
    }
  }

  # plot no split & combined
  if (is.null(x = split.by) && isFALSE(x = combine)) {
    if (seurat_version >= 5) {
      plot_list <- suppressMessages(FeaturePlot(object = seurat_object, features = all_found_features, order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, ncol = num_columns, combine = combine, raster.dpi = raster.dpi, label = label, alpha = alpha_exp, ...))

      plot <- lapply(1:length(x = plot_list), function(i) {
        p[[i]] <- suppressMessages(p[[i]] + scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, NA), na.value = na_color))
      })
    } else {
      plot_list <- suppressMessages(FeaturePlot(object = seurat_object, features = all_found_features, order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, ncol = num_columns, combine = combine, raster.dpi = raster.dpi, label = label, ...))

      plot <- lapply(1:length(x = plot_list), function(i) {
        p[[i]] <- suppressMessages(p[[i]] + scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, NA), na.value = na_color))
      })
    }
  }


  # plotting split with single feature (allows column number setting)
  if (!is.null(x = split.by) && length(x = all_found_features) == 1) {
    # Until Seurat is fixed pull feature data to set separately
    feature_data <- FetchData(
      object = seurat_object,
      vars = all_found_features,
      layer = layer)
    # Pull min and max values
    max_exp_value <- max(feature_data)
    min_exp_value <- min(feature_data)

    if (seurat_version >= 5) {
      plot <- suppressMessages(FeaturePlot(object = seurat_object, features = all_found_features, order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, raster.dpi = raster.dpi, label = label, alpha = alpha_exp, ...) & scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, max_exp_value), na.value = na_color, name = all_found_features)) & RestoreLegend() & theme(axis.title.y.right = element_blank())
    } else {
      plot <- suppressMessages(FeaturePlot(object = seurat_object, features = all_found_features, order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, raster.dpi = raster.dpi, label = label, ...) & scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, max_exp_value), na.value = na_color, name = all_found_features)) & RestoreLegend() & theme(axis.title.y.right = element_blank())
    }

    if (isTRUE(x = label_feature_yaxis)) {
      plot <- plot + plot_layout(nrow = num_rows, ncol = num_columns)
      plot <- plot & theme(legend.title=element_blank())
      plot <- suppressMessages(plot + scale_y_continuous(sec.axis = dup_axis(name = all_found_features))) + No_Right()
    } else {
      if (isTRUE(x = split_collect)) {
        if (hasArg("keep.scale")) {
          cli_abort(message = "The parameter {.code keep.scale} cannot be set different from default if {.code split_collect - TRUE}.")
        }
        plot <- plot + plot_layout(nrow = num_rows, ncol = num_columns, guides = "collect")
      } else {
        plot <- plot + plot_layout(nrow = num_rows, ncol = num_columns)
      }
    }
  }

  # plotting split multiple features
  if (!is.null(x = split.by) && length(x = all_found_features) > 1) {

    plot_list <- lapply(1:length(x = all_found_features), function(i){
      feature_data <- FetchData(
        object = seurat_object,
        vars = all_found_features[i],
        layer = layer)
      # Pull min and max values
      max_exp_value <- max(feature_data)
      min_exp_value <- min(feature_data)


      if (seurat_version >= 5) {
        single_plot <- suppressMessages(FeaturePlot(object = seurat_object, features = all_found_features[i], order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, raster.dpi = raster.dpi, label = label, alpha = alpha_exp, ...) & scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, max_exp_value), na.value = na_color, name = all_found_features[i])) & RestoreLegend() & theme(axis.title.y.right = element_blank())
      } else {
        single_plot <- suppressMessages(FeaturePlot(object = seurat_object, features = all_found_features[i], order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, raster.dpi = raster.dpi, label = label, ...) & scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, max_exp_value), na.value = na_color, name = features[i])) & RestoreLegend() & theme(axis.title.y.right = element_blank())
      }

      if (isTRUE(x = label_feature_yaxis)) {
        single_plot <- single_plot + plot_layout(nrow = num_rows, ncol = num_columns)
        single_plot <- single_plot & theme(legend.title=element_blank())
        single_plot <- suppressMessages(single_plot + scale_y_continuous(sec.axis = dup_axis(name = all_found_features[i]))) + No_Right()
      } else {
        single_plot <- single_plot + plot_layout(nrow = num_rows, ncol = num_columns)
      }
    })
    plot <- wrap_plots(plot_list) + plot_layout(ncol = 1)
  }

  # Add one time na_cutoff warning
  if (getOption(x = 'scCustomize_warn_na_cutoff', default = TRUE) && !is.na(x = na_cutoff) && na_cutoff == 0.000000001) {
    cli_inform(message = c("",
                           "NOTE: {.field FeaturePlot_scCustom} uses a specified {.code na_cutoff} when plotting to",
                           "color cells with no expression as background color separate from color scale.",
                           "Please ensure {.code na_cutoff} value is appropriate for feature being plotted.",
                           "Default setting is appropriate for use when plotting from 'RNA' assay.\n",
                           "When {.code na_cutoff} not appropriate (e.g., module scores) set to NULL to \n",
                           "plot all cells in gradient color palette.",
                           "",
                           "-----This message will be shown once per session.-----"))
    options(scCustomize_warn_na_cutoff = FALSE)
  }

  if (getOption(x = 'scCustomize_warn_zero_na_cutoff', default = TRUE) && !is.na(x = na_cutoff) && na_cutoff == 0) {
    cli_inform(message = c("",
                           "NOTE: Specified {.code na_cutoff} is set to {.field zero (0)}. This means that only cells/nuclei",
                           "with expression less than zero will be plotted with {.code na_color}: {.val {na_color}}.",
                           "To plot cells with expression values of zero using {.code na_color} leave",
                           "default {.code na_cutoff} value. If you want to plot full spectrum without",
                           "{.code na_cutoff} (e.g., for module scores) then set {.code na_cutoff = NULL`}.",
                           "",
                           "-----This message will be shown once per session.-----"))
    options(scCustomize_warn_na_cutoff = FALSE)
  }

  # Aspect ratio changes
  if (!is.null(x = aspect_ratio)) {
    if (!is.numeric(x = aspect_ratio)) {
      cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
    }
    plot <- plot & theme(aspect.ratio = aspect_ratio)
  }

  # Figure plot
  if (isTRUE(x = figure_plot)) {
    if (length(x = all_found_features) == 1) {
      plot <- Figure_Plot(plot = plot)
    } else {
      plot_list <- lapply(1:length(x = all_found_features), function(j) {
        fig_plot <- Figure_Plot(plot = plot[[j]])
      })

      plot <- wrap_plots(plot_list, ncol = num_columns)
    }
  }

  return(plot)
}


#' Customize FeaturePlot of two assays
#'
#' Create Custom FeaturePlots and preserve scale (no binning) from same features in two assays
#' simultaneously.  Intended for plotting same modality present in two assays.
#'
#' @param seurat_object Seurat object name.
#' @param features Feature(s) to plot.
#' @param assay1 name of assay one.  Default is "RAW" as featured in \code{\link{Create_CellBender_Merged_Seurat}}
#' @param assay2 name of assay two  Default is "RNA" as featured in \code{\link{Create_CellBender_Merged_Seurat}}
#' @param colors_use list of colors or color palette to use.
#' @param na_color color to use for points below lower limit.
#' @param order whether to move positive cells to the top (default = TRUE).
#' @param pt.size Adjust point size for plotting.
#' @param aspect_ratio Control the aspect ratio (y:x axes ratio length).  Must be numeric value;
#' Default is NULL.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default).
#' @param na_cutoff Value to use as minimum expression cutoff.  To set no cutoff set to `NA`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param slot `r lifecycle::badge("deprecated")` soft-deprecated.  See `layer`
#' @param layer Which layer to pull expression data from?  Default is "data".
#' @param num_columns Number of columns in plot layout.  If number of features > 1 then `num_columns`
#' dictates the number of columns in overall layout (`num_columns = 1` means stacked layout & `num_columns = 2`
#' means adjacent layout).
#' @param alpha_exp new alpha level to apply to expressing cell color palette (`colors_use`).  Must be
#' value between 0-1.
#' @param alpha_na_exp new alpha level to apply to non-expressing cell color palette (`na_color`).  Must be
#' value between 0-1.
#' @param ... Extra parameters passed to \code{\link[Seurat]{FeaturePlot}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @import patchwork
#' @importFrom Seurat FeaturePlot
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' \dontrun{
#' FeaturePlot_DualAssay(seurat_object = object, features = "Cx3cr1", assay1 = "RAW", assay2 = "RNA",
#' colors_use = viridis_plasma_dark_high, na_color = "lightgray")
#' }
#'

FeaturePlot_DualAssay <- function(
  seurat_object,
  features,
  assay1 = "RAW",
  assay2 = "RNA",
  colors_use = viridis_plasma_dark_high,
  na_color = "lightgray",
  order = TRUE,
  pt.size = NULL,
  aspect_ratio = NULL,
  reduction = NULL,
  na_cutoff = 0.000000001,
  raster = NULL,
  raster.dpi = c(512, 512),
  slot = deprecated(),
  layer = "data",
  num_columns = NULL,
  alpha_exp = NULL,
  alpha_na_exp = NULL,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check is slot is supplied
  if (lifecycle::is_present(slot)) {
    lifecycle::deprecate_warn(when = "2.0.0",
                              what = "slot",
                              with = "layer",
                              details = c("v" = "As of Seurat 5.0.0 the {.code slot} parameter is deprecated and replaced with {.code layer}.",
                                          "i" = "Please adjust code now to prepare for full deprecation.")
    )
    layer <- slot
  }

  # Check assays present
  assays_not_found <- Assay_Present(seurat_object = seurat_object, assay_list = c(assay1, assay2), print_msg = FALSE, omit_warn = TRUE)[[2]]

  if (!is.null(x = assays_not_found)) {
    stop_quietly()
  }

  # Find commands
  commands <- Command(object = seurat_object)

  # Raw normalize check
  if (!paste0("NormalizeData.", assay1) %in% commands) {
    cli_abort(message = c("Assay 1: {.val {assay1}} has not been normalized.",
                          "i" = "Please run {.code NormalizeData} on this assay before proceeding to visualization."))
  }

  # Cell Bender normalize check
  if (!paste0("NormalizeData.", assay2) %in% commands) {
    cli_abort(message = c("Assay 2: {.val {assay2}} has not been normalized.",
                      "i" = "Please run {.code NormalizeData} on this assay before proceeding to visualization."))
  }

  # Set columns if single feature
  num_features <- length(x = features)
  if (is.null(x = num_columns) && num_features == 1) {
    num_columns <- 2
  }

  # Set num columns if more than 1 feature
  if (is.null(x = num_columns) && num_features > 1) {
    num_columns <- 1
  }

  if (num_columns > 2 && num_features > 1) {
    cli_warn("When plotting more than one feature {.code num_columns} refers to patchwork columns and must either be 1 (vertical) or 2 (horizontal).")
  }

  # Change assay and plot raw
  DefaultAssay(object = seurat_object) <- assay1

  plot_raw <- FeaturePlot_scCustom(seurat_object = seurat_object, features = features, layer = layer, colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, order = order, pt.size = pt.size, reduction = reduction, raster = raster, alpha_exp = alpha_exp, alpha_na_exp = alpha_na_exp, raster.dpi = raster.dpi, ...) & labs(color = assay1)

  # Change to cell bender and plot
  DefaultAssay(object = seurat_object) <- assay2

  plot_cell_bender <- FeaturePlot_scCustom(seurat_object = seurat_object, features = features, layer = layer, colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, order = order, pt.size = pt.size, reduction = reduction, raster = raster, alpha_exp = alpha_exp, alpha_na_exp = alpha_na_exp, raster.dpi = raster.dpi, ...) & labs(color = assay2)

  # Assemble plots & return plots
  plots <- wrap_plots(plot_raw, plot_cell_bender, ncol = num_columns)

  # Aspect ratio changes
  if (!is.null(x = aspect_ratio)) {
    if (!is.numeric(x = aspect_ratio)) {
      cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
    }
    plots <- plots & theme(aspect.ratio = aspect_ratio)
  }

  return(plots)
}


#' Split FeatureScatter
#'
#' `r lifecycle::badge("deprecated")`
#' Create FeatureScatter using split.by
#'
#' @param seurat_object Seurat object name.
#' @param feature1 First feature to plot.
#' @param feature2 Second feature to plot.
#' @param split.by Feature to split plots by (i.e. "orig.ident").
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident).
#' Use 'ident' to group.by active.ident class.
#' @param colors_use color for the points on plot.
#' @param pt.size Adjust point size for plotting.
#' @param aspect_ratio Control the aspect ratio (y:x axes ratio length).  Must be numeric value;
#' Default is NULL.
#' @param title_size size for plot title labels.
#' @param num_columns number of columns in final layout plot.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 100,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ... Extra parameters passed to \code{\link[Seurat]{FeatureScatter}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @import patchwork
#' @importFrom magrittr "%>%"
#' @importFrom Seurat FeatureScatter
#' @importFrom stats cor
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' \dontrun{
#' # Function now DEPRECATED.
#' library(Seurat)
#' pbmc_small$sample_id <- sample(c("sample1", "sample2"), size = ncol(pbmc_small), replace = TRUE)
#'
#' # OLD Code
#' Split_FeatureScatter(seurat_object = pbmc_small, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
#' split.by = "sample_id")
#'
#' # NEW Code
#' FeatureScatter_scCustom(seurat_object = pbmc_small, feature1 = "nCount_RNA",
#' feature2 = "nFeature_RNA", split.by = "sample_id")
#'}
#'

Split_FeatureScatter <- function(
  seurat_object,
  feature1 = NULL,
  feature2 = NULL,
  split.by = NULL,
  group.by = NULL,
  colors_use = NULL,
  pt.size = NULL,
  aspect_ratio = NULL,
  title_size = 15,
  num_columns = NULL,
  raster = NULL,
  raster.dpi = c(512, 512),
  ggplot_default_colors = FALSE,
  color_seed = 123,
  ...
) {
  lifecycle::deprecate_stop(when = "2.0.0",
                            what = "Split_FeatureScatter()",
                            with = "FeatureScatter_scCustom()",
                            details = c("i" = "The functionality is now contained within `FeatureScatter_scCustom`")
  )

  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # split.by present
  if (is.null(x = split.by)) {
    cli_abort(message = "No value supplied to {.code split.by}.")
  }

  # Check split.by is valid
  if (split.by %in% colnames(seurat_object@meta.data) == FALSE) {
    cli_abort(message = c("The meta data variable: {.val {split.by}} could not be found in object@meta.data.",
                          "i" = "Please check the spelling and column names of meta.data slot.")
    )
  }

  # Set column and row lengths
  split.by_length <- length(x = unique(x = seurat_object@meta.data[[split.by]]))

   if (is.null(x = num_columns)) {
    num_columns <- split.by_length
  }
  # Calculate number of rows for selected number of columns
  num_rows <- ceiling(x = split.by_length/num_columns)

  # Check column and row compatibility
  if (num_columns > split.by_length) {
    cli_abort(message = c("The number of columns specified is greater than the number of meta data variables.",
                          "*" = "{.val {split.by}} only contains {.field {split.by_length}} variables.",
                          "i" = "Please adjust {.code num_columns} to be less than or equal to {.field {split.by_length}}.")
    )
  }

  # Check features are present
  possible_features <- c(rownames(x = seurat_object), colnames(x = seurat_object@meta.data))
  check_features <- setdiff(x = c(feature1, feature2), y = possible_features)
  if (length(x = check_features) > 0) {
    cli_abort(message = "The following feature(s) were not present in Seurat object: '{.field {check_features}}'")
  }

  # Extract min/maxes of features
  data_to_plot <- FetchData(object = seurat_object, vars = c(feature1, feature2))
  cor_data_features <- c("nCount_RNA", "nFeature_RNA")
  if (feature1 %in% cor_data_features && feature2 %in% cor_data_features) {
    min_feature1 <- min(data_to_plot[, feature1])-1
    max_feature1 <- max(data_to_plot[, feature1])+1
    min_feature2 <- min(data_to_plot[, feature2])-1
    max_feature2 <- max(data_to_plot[, feature2])+1
  } else {
    min_feature1 <- min(data_to_plot[, feature1])-0.05
    max_feature1 <- max(data_to_plot[, feature1])+0.05
    min_feature2 <- min(data_to_plot[, feature2])-0.05
    max_feature2 <- max(data_to_plot[, feature2])+0.05
  }

  # Extract split.by list of values
  if (inherits(x = seurat_object@meta.data[, split.by], what = "factor")) {
    meta_sample_list <- as.character(x = levels(x = seurat_object@meta.data[, split.by]))
  } else {
    meta_sample_list <- as.character(x = unique(x = seurat_object@meta.data[, split.by]))
  }

  # Extract cell names per meta data list of values
  cell_names <- lapply(meta_sample_list, function(x) {
    row.names(x = seurat_object@meta.data)[which(x = seurat_object@meta.data[, split.by] == x)]})

  # raster check
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)

  # Set uniform point size is pt.size = NULL (based on plot with most cells)
  if (is.null(x = pt.size)) {
    # cells per meta data
    cells_by_meta <- data.frame(table(seurat_object@meta.data[, split.by]))
    # Identity with greatest number of cells
    max_cells <- max(cells_by_meta$Freq)
    # modified version of the autopointsize function from Seurat
    pt.size <- AutoPointSize_scCustom(data = max_cells, raster = raster)
  }

  # Add correlations if applicable
  cor_data_features <- c("nCount_RNA", "nFeature_RNA")
  if (feature1 %in% cor_data_features && feature2 %in% cor_data_features) {
    plot_cor <- TRUE
    cor_data <- FetchData(object = seurat_object, vars = c("nCount_RNA", "nFeature_RNA", split.by))

    cor_values <- lapply(1:length(x = meta_sample_list), function(i) {
      cor_data_filtered <- cor_data %>%
        filter(.data[[split.by]] == meta_sample_list[[i]])
      round(x = cor(x = cor_data_filtered[, "nCount_RNA"], y = cor_data_filtered[, "nFeature_RNA"]), digits = 2)
    })
  } else {
    plot_cor <- FALSE
  }

  # Set colors
  group.by <- group.by %||% 'ident'

  if (group.by == "ident") {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }

  if (is.null(x = colors_use)) {
    # set default plot colors
    if (is.null(x = colors_use)) {
      colors_use <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
    }
  }

  # Plots
  plots <- lapply(1:length(x = meta_sample_list), function(j) {
    plot <- FeatureScatter(seurat_object, feature1 = feature1, feature2 = feature2, cells = cell_names[[j]], group.by = group.by, cols = colors_use, pt.size = pt.size, raster = raster, raster.dpi = raster.dpi, ...) +
      theme(plot.title = element_text(hjust = 0.5, size = title_size),
            legend.position = "right") +
      xlim(min_feature1, max_feature1) +
      ylim(min_feature2, max_feature2)
    if (isTRUE(x = plot_cor)) {
      plot + ggtitle(paste(meta_sample_list[[j]]), subtitle = paste0("Correlation: ", cor_values[j]))
    } else {
      plot + ggtitle(paste(meta_sample_list[[j]]))
    }
  })

  # Wrap Plots into single output
  plot_comb <- wrap_plots(plots, ncol = num_columns, nrow = num_rows) + plot_layout(guides = 'collect')

  # Aspect ratio changes
  if (!is.null(x = aspect_ratio)) {
    if (!is.numeric(x = aspect_ratio)) {
      cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
    }
    plot_comb <- plot_comb & theme(aspect.ratio = aspect_ratio)
  }

  return(plot_comb)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GENE EXPRESSION PLOTTING (NON-2D) ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' VlnPlot with modified default settings
#'
#' Creates DimPlot with some of the settings modified from their Seurat defaults (colors_use, shuffle, label).
#'
#' @param seurat_object Seurat object name.
#' @param features Feature(s) to plot.
#' @param colors_use color palette to use for plotting.  By default if number of levels plotted is less than
#' or equal to 36 it will use "polychrome" and if greater than 36 will use "varibow" with shuffle = TRUE
#' both from `DiscretePalette_scCustomize`.
#' @param pt.size Adjust point size for plotting.
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident);
#' default is the current active.ident of the object.
#' @param split.by Feature to split plots by (i.e. "orig.ident").
#' @param plot_median logical, whether to plot median for each ident on the plot (Default is FALSE).
#' @param plot_boxplot logical, whether to plot boxplot inside of violin (Default is FALSE).
#' @param median_size Shape size for the median is plotted.
#' @param idents Which classes to include in the plot (default is all).
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 100,000 total points plotted (# Cells x # of features).
#' @param add.noise logical, determine if adding a small noise for plotting (Default is TRUE).
#' @param num_columns Number of columns in plot layout.  Only valid if `split.by != NULL`.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ... Extra parameters passed to \code{\link[Seurat]{VlnPlot}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import patchwork
#' @import ggrastr
#' @importFrom Seurat VlnPlot
#'
#' @export
#'
#' @references Many of the param names and descriptions are from Seurat to facilitate ease of use as
#' this is simply a wrapper to alter some of the default parameters \url{https://github.com/satijalab/seurat/blob/master/R/visualization.R} (License: GPL-3).
#'
#' @concept seurat_plotting
#'
#' @examples
#' library(Seurat)
#' VlnPlot_scCustom(seurat_object = pbmc_small, features = "CD3E")
#'

VlnPlot_scCustom <- function(
  seurat_object,
  features,
  colors_use = NULL,
  pt.size = NULL,
  group.by = NULL,
  split.by = NULL,
  plot_median = FALSE,
  plot_boxplot = FALSE,
  median_size = 15,
  idents = NULL,
  num_columns = NULL,
  raster = NULL,
  add.noise = TRUE,
  ggplot_default_colors = FALSE,
  color_seed = 123,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check split valid
  if (!is.null(x = split.by)) {
    split.by <- Meta_Present(seurat_object = seurat_object, meta_col_names = split.by, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  # Add check for group.by before getting to colors
  if (!is.null(x = group.by) && group.by != "ident") {
    Meta_Present(seurat_object = seurat_object, meta_col_names = group.by, print_msg = FALSE)
  }

  # Check features and meta to determine which features present
  all_found_features <- Feature_PreCheck(object = seurat_object, features = features)

  # Check boxplot vs median
  if (isTRUE(x = plot_median) && isTRUE(x = plot_boxplot)) {
    cli_abort(message = c("Incompatible settings.",
                          "{.code plot_median} and {.code plot_boxplot} cannot both be set to TRUE."))
  }

  # set size if NULL
  if (isTRUE(x = plot_boxplot)) {
    if (!is.null(x = pt.size)) {
      cli_warn(message = c("Provided value for {.code pt.size} ({.field {pt.size}}) will be ignored.",
                                "When setting {.field plot_boxplot = TRUE}, {.code pt.size} is automatically set to 0."))
    }
    pt.size <- 0
  } else {
    pt.size <- pt.size %||% AutoPointSize_scCustom(data = seurat_object)
  }

  # check median vs boxplot
  if (isTRUE(x = plot_median) && isTRUE(x = plot_boxplot)) {
    cli_abort(message = "{.code plot_median} and {.code plot_boxplot} cannot both be set {.field TRUE}")
  }

  # Add raster check for scCustomize
  # num_cells <- unlist(x = CellsByIdentities(object = seurat_object, idents = idents))
  num_cells <- length(x = Cells(x = seurat_object))

  if (is.null(x = raster)) {
    if (pt.size == 0) {
      raster <- FALSE
    } else {
      if (length(x = num_cells) * length(x = all_found_features) > 100000 && pt.size != 0) {
        raster <- TRUE
        cli_inform(message = c("NOTE: Rasterizing points since total number of points across all plots exceeds 100,000.",
                               "i" = "To plot in vector form set {.code raster=FALSE}")
        )
      } else {
        raster <- FALSE
      }
    }
  }

  # Set default color palette based on number of levels being plotted
  if (is.null(x = split.by)) {
    if (is.null(x = group.by)) {
      group_by_length <- length(x = unique(x = seurat_object@active.ident))
    } else {
      group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
    }
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[split.by]]))
  }


  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use) && isTRUE(x = ggplot_default_colors)) {
    cli_abort(message = "Cannot provide both custom palette to {.code colors_use} and specify {.code ggplot_default_colors = TRUE}.")
  }

  # set default plot colors
  if (is.null(x = colors_use)) {
    colors_use <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
  }

  # Plot
  plot <- VlnPlot(object = seurat_object, features = all_found_features, cols = colors_use, pt.size = pt.size, idents = idents, group.by = group.by, split.by = split.by, ncol = num_columns, raster = raster, add.noise = add.noise, ...)

  # Add add median plot
  if (isTRUE(x = plot_median) && is.null(x = split.by)) {
    plot <- plot + stat_summary(fun = median, geom='point', size = median_size, colour = "white", shape = 95)
  }

  if (isTRUE(x = plot_median) && !is.null(x = split.by)) {
    cli_abort(message = "Cannot add median ({.field plot_median = TRUE}) when {.code split.by} is not NULL.")
  }

  if (isTRUE(x = plot_boxplot) && is.null(x = split.by)) {
    plot <- plot + geom_boxplot(fill='#A4A4A4', color="black", width = 0.1)
  }
  if (isTRUE(x = plot_boxplot) && !is.null(x = split.by)) {
    cli_abort(message = "Cannot add median ({.field plot_boxplot = TRUE}) when {.code split.by} is not NULL.")
  }

  return(plot)
}


#' Stacked Violin Plot
#'
#' Code for creating stacked violin plot gene expression.
#'
#' @param seurat_object Seurat object name.
#' @param features Features to plot.
#' @param group.by Group (color) cells in different ways (for example, orig.ident).
#' @param split.by A variable to split the violin plots by,
#' @param idents Which classes to include in the plot (default is all).
#' @param x_lab_rotate logical or numeric.  If logical whether to rotate x-axis labels 45 degrees
#'  (Default is FALSE).  If numeric must be either 45 or 90.  Setting 45 is equivalent to setting TRUE.
#' @param plot_legend logical.  Adds plot legend containing `idents` to the returned plot.
#' @param colors_use specify color palette to used in \code{\link[Seurat]{VlnPlot}}.  By default if
#' number of levels plotted is less than or equal to 36 it will use "polychrome" and if greater than 36
#' will use "varibow" with shuffle = TRUE both from `DiscretePalette_scCustomize`.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param plot_spacing Numerical value specifying the vertical spacing between each plot in the stack.
#' Default is 0.15 ("cm").  Spacing dependent on unit provided to `spacing_unit`.
#' @param spacing_unit Unit to use in specifying vertical spacing between plots.  Default is "cm".
#' @param vln_linewidth Adjust the linewidth of violin outline.  Must be numeric.
#' @param pt.size Adjust point size for plotting.  Default for `Stacked_VlnPlot` is 0 to avoid issues with
#' rendering so many points in vector form.  Alternatively, see `raster` parameter.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 100,000 total points plotted (# Cells x # of features).
#' @param add.noise logical, determine if adding a small noise for plotting (Default is TRUE).
#' @param ... Extra parameters passed to \code{\link[Seurat]{VlnPlot}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @import patchwork
#' @importFrom purrr map map_dbl map2
#' @importFrom Seurat VlnPlot
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @author Ming Tang (Original Code), Sam Marsh (Wrap single function, added/modified functionality)
#' @references \url{https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/}
#' @seealso \url{https://twitter.com/tangming2005}
#'
#' @examples
#' library(Seurat)
#' Stacked_VlnPlot(seurat_object = pbmc_small, features = c("CD3E", "CD8", "GZMB", "MS4A1"),
#' x_lab_rotate = TRUE)
#'

Stacked_VlnPlot <- function(
  seurat_object,
  features,
  group.by = NULL,
  split.by = NULL,
  idents = NULL,
  x_lab_rotate = FALSE,
  plot_legend = FALSE,
  colors_use = NULL,
  color_seed = 123,
  ggplot_default_colors = FALSE,
  plot_spacing = 0.15,
  spacing_unit = "cm",
  vln_linewidth = NULL,
  pt.size = NULL,
  raster = NULL,
  add.noise = TRUE,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Add check for group.by before getting to colors
  if (!is.null(x = group.by) && group.by != "ident") {
    Meta_Present(seurat_object = seurat_object, meta_col_names = group.by, print_msg = FALSE)
  }

  # Check features and meta to determine which features present
  all_found_features <- Feature_PreCheck(object = seurat_object, features = features)

  # set pt.size (default is no points)
  if (is.null(x = pt.size)) {
    pt.size <- 0
    if (isTRUE(x = raster)) {
      cli_inform(message = "Default pt.size is 0, setting {.code raster = FALSE}.")
    }
    raster <- FALSE
  }

  # Set rasterization
  num_cells <- length(x = Cells(x = seurat_object))

  if (num_cells * length(x = all_found_features) > 100000 && is.null(x = raster) && pt.size != 0) {
    raster <- TRUE
    cli_inform(message = c("NOTE: Rasterizing points since total number of points across all plots exceeds 100,000.",
                           "i" = "To plot in vector form set {.code raster=FALSE}")
    )
  }

  # Set default color palette based on number of levels being plotted
  if (is.null(x = split.by)) {
    if (is.null(x = group.by)) {
      group_by_length <- length(x = unique(x = seurat_object@active.ident))
    } else {
      group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
    }
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[split.by]]))
  }

  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use) && isTRUE(x = ggplot_default_colors)) {
    cli_abort(message = "Cannot provide both custom palette to {.code colors_use} and specify {.code ggplot_default_colors = TRUE}.")
  }
  if (is.null(x = colors_use)) {
    # set default plot colors
    if (is.null(x = colors_use)) {
      colors_use <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
    }
  }

  # Setup plot spacing/margin parameter
  plot_margin <- margin(t = plot_spacing, r = 0, b = plot_spacing, l = 0, unit = spacing_unit)

  # Create plots
  plot_list <- map(all_found_features, function(x) Modify_VlnPlot(seurat_object = seurat_object, features = x, cols = colors_use, group.by = group.by, split.by = split.by, idents = idents, plot_margin = plot_margin, pt.size = pt.size, raster = raster, add.noise = add.noise, ...))

  # Add back x-axis title to bottom plot. patchwork is going to support this?
  # Add ability to rotate the X axis labels to the function call
  if (isTRUE(x = x_lab_rotate) || x_lab_rotate == 45) {
    plot_list[[length(x = plot_list)]] <- plot_list[[length(x = plot_list)]] +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.ticks.x = element_line())
  }

  if (x_lab_rotate == 90) {
    plot_list[[length(x = plot_list)]] <- plot_list[[length(x = plot_list)]] +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.ticks.x = element_line())
  }

  if (!is.logical(x = x_lab_rotate) && !x_lab_rotate %in% c(45, 90)) {
    cli_abort(message = "{.code x_lab_rotate} must be either a logical or a numeric value of 45 or 90.")
  }

  plot_list[[length(x = plot_list)]] <- plot_list[[length(x = plot_list)]] +
    theme(axis.text.x = element_text(), axis.ticks.x = element_line())

  # change the y-axis tick to only max value
  ymaxs <- map_dbl(plot_list, Extract_Max)
  plot_list <- suppressMessages(map2(plot_list, ymaxs, function(x,y) x +
                                       scale_y_continuous(breaks = c(y)) +
                                       expand_limits(y = y)))

  # Create final plot patchwork to return
  plot_return <- wrap_plots(plotlist = plot_list, ncol = 1)

  # add single legend
  if (plot_legend) {
    plot_return <- plot_return & theme(legend.position = "right")
    plot_return <- plot_return + plot_layout(guides = 'collect')
  }

  if (!is.null(x = vln_linewidth)) {
    if (!is.numeric(x = vln_linewidth)) {
      cli_abort(message = "{.code vln_linewidth} parameter must be numeric.")
    }
    for (j in 1:length(x = plot_list)) {
      plot_return[[j]]$layers[[1]]$aes_params$linewidth <- vln_linewidth
    }
  }

  # return plot
  return(plot_return)
}


#' Customized DotPlot
#'
#' Code for creating customized DotPlot
#'
#' @param seurat_object Seurat object name.
#' @param features Features to plot.
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident);
#' default is the current active.ident of the object.
#' @param colors_use specify color palette to used.  Default is viridis_plasma_dark_high.
#' @param remove_axis_titles logical. Whether to remove the x and y axis titles.  Default = TRUE.
#' @param x_lab_rotate Rotate x-axis labels 45 degrees (Default is FALSE).
#' @param y_lab_rotate Rotate x-axis labels 45 degrees (Default is FALSE).
#' @param facet_label_rotate Rotate facet labels on grouped `DotPlots` by 45 degrees (Default is FALSE).
#' @param flip_axes whether or not to flip and X and Y axes (Default is FALSE).
#' @param ... Extra parameters passed to \code{\link[Seurat]{DotPlot}}.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom Seurat DotPlot
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' \donttest{
#' library(Seurat)
#' DotPlot_scCustom(seurat_object = pbmc_small, features = c("CD3E", "CD8", "GZMB", "MS4A1"))
#'}
#'

DotPlot_scCustom <- function(
  seurat_object,
  features,
  group.by = NULL,
  colors_use = viridis_plasma_dark_high,
  remove_axis_titles = TRUE,
  x_lab_rotate = FALSE,
  y_lab_rotate = FALSE,
  facet_label_rotate = FALSE,
  flip_axes = FALSE,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Add check for group.by before getting to colors
  if (!is.null(x = group.by) && group.by != "ident") {
    Meta_Present(seurat_object = seurat_object, meta_col_names = group.by, print_msg = FALSE)
  }

  # Check features and meta to determine which features present
  all_found_features <- Feature_PreCheck(object = seurat_object, features = features)

  # Plot
  plot <- suppressMessages(DotPlot(object = seurat_object, features = all_found_features, ...) +
                             scale_color_gradientn(colors = colors_use)
  )
  # Modify plot
  if (isTRUE(x = remove_axis_titles)) {
    plot <- plot +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank()
      )
  }

  if (isTRUE(x = flip_axes)) {
    plot <- plot & coord_flip()
  }
  if (isTRUE(x = x_lab_rotate)) {
    plot <- plot +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  if (isTRUE(x = y_lab_rotate)) {
    plot <- plot +
      theme(axis.text.y = element_text(angle = 45, hjust = 1))
  }

  if (is.list(x = features)) {
    plot + theme(strip.text.x = element_text(angle = 45))
  }

  return(plot)
}


#' Clustered DotPlot
#'
#' Clustered DotPlots using ComplexHeatmap
#'
#' @param seurat_object Seurat object name.
#' @param features Features to plot.
#' @param split.by Variable in `@meta.data` to split the identities plotted by.
#' @param colors_use_exp Color palette to use for plotting expression scale.  Default is `viridis::plasma(n = 20, direction = -1)`.
#' @param exp_color_min Minimum scaled average expression threshold (everything smaller will be set to this).
#' Default is -2.
#' @param exp_color_middle What scaled expression value to use for the middle of the provided `colors_use_exp`.
#' By default will be set to value in middle of `exp_color_min` and `exp_color_max`.
#' @param exp_color_max Minimum scaled average expression threshold (everything smaller will be set to this).
#' Default is 2.
#' @param exp_value_type Whether to plot average normalized expression or
#' scaled average normalized expression.  Only valid when `split.by` is provided.
#' @param print_exp_quantiles Whether to print the quantiles of expression data in addition to plots.
#' Default is FALSE.  NOTE: These values will be altered by choices of `exp_color_min` and `exp_color_min`
#' if there are values below or above those cutoffs, respectively.
#' @param colors_use_idents specify color palette to used for identity labels.  By default if
#' number of levels plotted is less than or equal to 36 it will use "polychrome" and if greater than 36
#' will use "varibow" with shuffle = TRUE both from `DiscretePalette_scCustomize`.
#' @param x_lab_rotate How to rotate column labels.  By default set to `TRUE` which rotates labels 45 degrees.
#' If set `FALSE` rotation is set to 0 degrees.  Users can also supply custom angle for text rotation.
#' @param flip logical, whether to flip the axes of final plot.  Default is FALSE; rows = features and
#' columns = idents.
#' @param k Value to use for k-means clustering on features  Sets (km) parameter in `ComplexHeatmap::Heatmap()`.
#' From `ComplexHeatmap::Heatmap()`: Apply k-means clustering on rows. If the value is larger than 1, the
#' heatmap will be split by rows according to the k-means clustering. For each row slice, hierarchical
#' clustering is still applied with parameters above.
#' @param feature_km_repeats Number of k-means runs to get a consensus k-means clustering for features.
#' Note if `feature_km_repeats` is set to value greater than one, the final number of groups might be
#' smaller than row_km, but this might mean the original row_km is not a good choice.  Default is 1000.
#' @param ident_km_repeats Number of k-means runs to get a consensus k-means clustering. Similar to
#' `feature_km_repeats`.  Default is 1000.
#' @param row_label_size Size of the feature labels.  Provided to `row_names_gp` in Heatmap call.
#' @param row_label_fontface Fontface to use for row labels.  Provided to `row_names_gp` in Heatmap call.
#' @param grid_color color to use for heatmap grid.  Default is NULL which "removes" grid by using NA color.
#' @param cluster_feature logical, whether to cluster and reorder feature axis.  Default is TRUE.
#' @param cluster_ident logical, whether to cluster and reorder identity axis.  Default is TRUE.
#' @param column_label_size Size of the feature labels.  Provided to `column_names_gp` in Heatmap call.
#' @param legend_label_size Size of the legend text labels.  Provided to `labels_gp` in Heatmap legend call.
#' @param legend_title_size Sise of the legend title text labels.  Provided to `title_gp` in Heatmap legend call.
#' @param raster Logical, whether to render in raster format (faster plotting, smaller files).  Default is FALSE.
#' @param plot_km_elbow Logical, whether or not to return the Sum Squared Error Elbow Plot for k-means clustering.
#' Estimating elbow of this plot is one way to determine "optimal" value for `k`.
#' Based on: \url{https://stackoverflow.com/a/15376462/15568251}.
#' @param elbow_kmax The maximum value of k to use for `plot_km_elbow`.  Suggest setting larger value so the
#' true shape of plot can be observed.  Value must be 1 less than number of features provided.  If NULL parameter
#' will be set dependent on length of feature list up to `elbow_kmax = 20`.
#' @param assay Name of assay to use, defaults to the active assay.
#' @param group.by Group (color) cells in different ways (for example, orig.ident).
#' @param idents Which classes to include in the plot (default is all).
#' @param show_parent_dend_line Logical, Sets parameter of same name in `ComplexHeatmap::Heatmap()`.
#' From `ComplexHeatmap::Heatmap()`: When heatmap is split, whether to add a dashed line to mark parent
#' dendrogram and children dendrograms.  Default is TRUE.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param seed Sets seed for reproducible plotting (ComplexHeatmap plot).
#'
#' @return A ComplexHeatmap or if plot_km_elbow = TRUE a list containing ggplot2 object and ComplexHeatmap.
#'
#' @import cli
#' @import ggplot2
#' @importFrom circlize colorRamp2
#' @importFrom dplyr any_of filter select
#' @importFrom grid grid.circle grid.rect gpar
#' @importFrom magrittr "%>%"
#' @importFrom rlang is_installed
#' @importFrom Seurat DotPlot
#' @importFrom stats quantile
#' @importFrom tidyr pivot_wider
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @author Ming Tang (Original Code), Sam Marsh (Wrap single function, added/modified functionality)
#' @references \url{https://divingintogeneticsandgenomics.rbind.io/post/clustered-dotplot-for-single-cell-rnaseq/}
#' @seealso \url{https://twitter.com/tangming2005}
#'
#' @examples
#' \donttest{
#' library(Seurat)
#' Clustered_DotPlot(seurat_object = pbmc_small, features = c("CD3E", "CD8", "GZMB", "MS4A1"))
#'}
#'

Clustered_DotPlot <- function(
  seurat_object,
  features,
  split.by = NULL,
  colors_use_exp = viridis_plasma_dark_high,
  exp_color_min = -2,
  exp_color_middle = NULL,
  exp_color_max = 2,
  exp_value_type = "scaled",
  print_exp_quantiles = FALSE,
  colors_use_idents = NULL,
  x_lab_rotate = TRUE,
  flip = FALSE,
  k = 1,
  feature_km_repeats = 1000,
  ident_km_repeats = 1000,
  row_label_size = 8,
  row_label_fontface = "plain",
  grid_color = NULL,
  cluster_feature = TRUE,
  cluster_ident = TRUE,
  column_label_size = 8,
  legend_label_size = 10,
  legend_title_size = 10,
  raster = FALSE,
  plot_km_elbow = TRUE,
  elbow_kmax = NULL,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  show_parent_dend_line = TRUE,
  ggplot_default_colors = FALSE,
  plot_width = NULL,
  plot_height = NULL,
  color_seed = 123,
  seed = 123
) {
  # check split
  if (is.null(x = split.by)) {
    Clustered_DotPlot_Single_Group(seurat_object = seurat_object,
                                   features = features,
                                   colors_use_exp = colors_use_exp,
                                   exp_color_min = exp_color_min,
                                   exp_color_middle = exp_color_middle,
                                   exp_color_max = exp_color_max,
                                   print_exp_quantiles = print_exp_quantiles,
                                   colors_use_idents = colors_use_idents,
                                   x_lab_rotate = x_lab_rotate,
                                   flip = flip,
                                   k = k,
                                   feature_km_repeats = feature_km_repeats,
                                   ident_km_repeats = ident_km_repeats,
                                   row_label_size = row_label_size,
                                   row_label_fontface = row_label_fontface,
                                   grid_color = grid_color,
                                   cluster_feature = cluster_feature,
                                   cluster_ident = cluster_ident,
                                   column_label_size = column_label_size,
                                   legend_label_size = legend_label_size,
                                   legend_title_size = legend_title_size,
                                   raster = raster,
                                   plot_km_elbow = plot_km_elbow,
                                   elbow_kmax = elbow_kmax,
                                   assay = assay,
                                   group.by = group.by,
                                   idents = idents,
                                   show_parent_dend_line = show_parent_dend_line,
                                   ggplot_default_colors = ggplot_default_colors,
                                   plot_width = plot_width,
                                   plot_height = plot_height,
                                   color_seed = color_seed,
                                   seed = seed)
  } else {
    Clustered_DotPlot_Multi_Group(seurat_object = seurat_object,
                                  features = features,
                                  split.by = split.by,
                                  colors_use_exp = colors_use_exp,
                                  exp_color_min = exp_color_min,
                                  exp_color_middle = exp_color_middle,
                                  exp_color_max = exp_color_max,
                                  exp_value_type = exp_value_type,
                                  print_exp_quantiles = print_exp_quantiles,
                                  x_lab_rotate = x_lab_rotate,
                                  flip = flip,
                                  k = k,
                                  feature_km_repeats = feature_km_repeats,
                                  ident_km_repeats = ident_km_repeats,
                                  row_label_size = row_label_size,
                                  row_label_fontface = row_label_fontface,
                                  grid_color = grid_color,
                                  cluster_feature = cluster_feature,
                                  cluster_ident = cluster_ident,
                                  column_label_size = column_label_size,
                                  legend_label_size = legend_label_size,
                                  legend_title_size = legend_title_size,
                                  raster = raster,
                                  plot_km_elbow = plot_km_elbow,
                                  elbow_kmax = elbow_kmax,
                                  assay = assay,
                                  group.by = group.by,
                                  idents = idents,
                                  show_parent_dend_line = show_parent_dend_line,
                                  plot_width = plot_width,
                                  plot_height = plot_height,
                                  seed = seed)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### 2D NON-GENE EXPRESSION PLOTTING ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Cluster Highlight Plot
#'
#' Create Plot with cluster of interest highlighted
#'
#' @param seurat_object Seurat object name.
#' @param cluster_name Name(s) (or number(s)) identity of cluster to be highlighted.
#' @param highlight_color Color(s) to highlight cells.  The default is NULL and plot will use
#' `scCustomize_Palette()`.
#' @param background_color non-highlighted cell colors.
#' @param pt.size point size for both highlighted cluster and background.
#' @param aspect_ratio Control the aspect ratio (y:x axes ratio length).  Must be numeric value;
#' Default is NULL.
#' @param figure_plot logical.  Whether to remove the axes and plot with legend on left of plot denoting
#' axes labels.  (Default is FALSE).  Requires `split_seurat = TRUE`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param label Whether to label the highlighted cluster(s).  Default is FALSE.
#' @param split.by Feature to split plots by (i.e. "orig.ident").
#' @param split_seurat logical.  Whether or not to display split plots like Seurat (shared y axis) or as
#' individual plots in layout.  Default is FALSE.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param ... Extra parameters passed to \code{\link[Seurat]{DimPlot}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @import patchwork
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' Cluster_Highlight_Plot(seurat_object = pbmc_small, cluster_name = "1", highlight_color = "gold",
#' background_color = "lightgray",  pt.size = 2)
#'

Cluster_Highlight_Plot <- function(
  seurat_object,
  cluster_name,
  highlight_color = NULL,
  background_color = "lightgray",
  pt.size = NULL,
  aspect_ratio = NULL,
  figure_plot = FALSE,
  raster = NULL,
  raster.dpi = c(512, 512),
  label = FALSE,
  split.by = NULL,
  split_seurat = FALSE,
  ggplot_default_colors = FALSE,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)

  # Perform Idents check and report errors when when length(cluster_name) > 1
  if (length(x = cluster_name) > 1) {
    idents_list <- levels(x = Idents(object = seurat_object))

    good_idents <- cluster_name[cluster_name %in% idents_list]
    bad_idents <- cluster_name[!cluster_name %in% idents_list]

    if (length(x = bad_idents) > 0) {
      cli_warn("The following `cluster_name{?s}` were omitted as they were not found the active.ident slot: {.field {bad_idents}}")
    }
  }

  # pull cells to highlight in plot
  cells_to_highlight <- CellsByIdentities(seurat_object, idents = cluster_name)

  # set point size
  if (is.null(x = pt.size)) {
    pt.size <- AutoPointSize_scCustom(data = sum(lengths(x = cells_to_highlight)), raster = raster)
  }

  # Set colors
  # Adjust colors if needed when length(cluster_name) > 1
  if (length(x = highlight_color) == 1 && length(x = cluster_name) > 1) {
    highlight_color <- rep(x = highlight_color, length(x = cluster_name))
    cli_inform(message = c("NOTE: Only one color provided to but {.field {length(x = cluster_name)}} clusters were provided.",
                           "i" = "Using the same color ({.val {highlight_color[1]}}) for all clusters."))
  }

  # If NULL set using scCustomize_Palette
  if (is.null(x = highlight_color)) {
    highlight_color <- scCustomize_Palette(num_groups = length(x = cells_to_highlight), ggplot_default_colors = ggplot_default_colors)
  }

  # plot
  plot <- DimPlot_scCustom(seurat_object = seurat_object,
          cells.highlight = cells_to_highlight,
          cols.highlight = highlight_color,
          colors_use = background_color,
          sizes.highlight = pt.size,
          pt.size = pt.size,
          order = TRUE,
          raster = raster,
          raster.dpi = raster.dpi,
          split.by = split.by,
          split_seurat = split_seurat,
          label = label,
          ...)

  # Edit plot legend
  plot <- suppressMessages(plot & scale_color_manual(breaks = names(x = cells_to_highlight), values = c(highlight_color, background_color), na.value = background_color))

  # Aspect ratio changes
  if (!is.null(x = aspect_ratio)) {
    if (!is.numeric(x = aspect_ratio)) {
      cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
    }
    plot <- plot & theme(aspect.ratio = aspect_ratio)
  }

  # Figure plot
  if (isTRUE(x = figure_plot)) {
    plot <- Figure_Plot(plot = plot)
  }

  return(plot)
}


#' Meta Highlight Plot
#'
#' Create Plot with meta data variable of interest highlighted
#'
#' @param seurat_object Seurat object name.
#' @param meta_data_column Name of the column in `seurat_object@meta.data` slot to pull value from for highlighting.
#' @param meta_data_highlight Name of variable(s) within `meta_data_name` to highlight in the plot.
#' @param highlight_color Color to highlight cells (default "navy").
#' @param background_color non-highlighted cell colors.
#' @param pt.size point size for both highlighted cluster and background.
#' @param aspect_ratio Control the aspect ratio (y:x axes ratio length).  Must be numeric value;
#' Default is NULL.
#' @param figure_plot logical.  Whether to remove the axes and plot with legend on left of plot denoting
#' axes labels.  (Default is FALSE).  Requires `split_seurat = TRUE`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param label Whether to label the highlighted meta data variable(s).  Default is FALSE.
#' @param split.by Variable in `@meta.data` to split the plot by.
#' @param split_seurat logical.  Whether or not to display split plots like Seurat (shared y axis) or as
#' individual plots in layout.  Default is FALSE.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param ... Extra parameters passed to\code{\link[Seurat]{DimPlot}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @import patchwork
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' library(Seurat)
#' pbmc_small$sample_id <- sample(c("sample1", "sample2"), size = ncol(pbmc_small), replace = TRUE)
#'
#' Meta_Highlight_Plot(seurat_object = pbmc_small, meta_data_column = "sample_id",
#' meta_data_highlight = "sample1", highlight_color = "gold", background_color = "lightgray",
#' pt.size = 2)
#'

Meta_Highlight_Plot <- function(
  seurat_object,
  meta_data_column,
  meta_data_highlight,
  highlight_color = NULL,
  background_color = "lightgray",
  pt.size = NULL,
  aspect_ratio = NULL,
  figure_plot = FALSE,
  raster = NULL,
  raster.dpi = c(512, 512),
  label = FALSE,
  split.by = NULL,
  split_seurat = FALSE,
  ggplot_default_colors = FALSE,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check meta data
  good_meta_data_column <- Meta_Present(seurat_object = seurat_object, meta_col_names = meta_data_column, omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)[[1]]

  # stop if none found
  if (length(x = good_meta_data_column) == 0) {
    cli_abort(message = c("{.code meta_data_column} was not found.",
              "i" = "No column found in object meta.data named: {.val {meta_data_column}}.")
    )
  }

  # Check that meta data is factor or character
  accepted_meta_types <- c("factor", "character", "logical")

  if (!class(x = seurat_object@meta.data[[good_meta_data_column]]) %in% accepted_meta_types) {
    cli_abort(message = c("The {.code good_meta_data_column}: {.field {good_meta_data_column}} is of class: {.val {class(x = seurat_object@meta.data[[good_meta_data_column]])}}.",
                          "i" = "Meta data variables must be of classes: factor, character, or logical to be used with {.code Meta_Highlight_Plot()}.")
              )
  }

  # Check meta_data_highlight
  meta_var_list <- as.character(x = unique(x = seurat_object@meta.data[, good_meta_data_column]))

  # Check good and bad highlight values
  bad_meta_highlight <- meta_var_list[!meta_var_list %in% meta_data_highlight]
  found_meta_highlight <- meta_var_list[meta_var_list %in% meta_data_highlight]

  # Abort if no meta_data_highlight found
  if (length(x = found_meta_highlight) == 0) {
    cli_abort(message = c("No 'meta_data_highlight' value(s) were not found.",
                          "i" = "The following {.code meta_data_highlight} variables were not found in {.field {good_meta_data_column}} and were omitted: {.field {bad_meta_highlight}}")
    )
  }

  # warn if some meta_data_highlight not found
  if (length(x = found_meta_highlight) != length(x = meta_data_highlight)) {
    cli_warn(message = c("Some {.code meta_data_highlight} value(s) were not found.",
                          "i" = "The following {.code meta_data_highlight} variables were not found in {.field {good_meta_data_column}} and were omitted: {.field {bad_meta_highlight}}")
    )
  }

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)

  # Change default ident and pull cells to highlight in plot
  Idents(object = seurat_object) <- good_meta_data_column

  cells_to_highlight <- CellsByIdentities(seurat_object, idents = found_meta_highlight)

  # set point size
  if (is.null(x = pt.size)) {
    pt.size <- AutoPointSize_scCustom(data = sum(lengths(x = cells_to_highlight)), raster = raster)
  }

  # Set colors
  # Adjust colors if needed when length(meta_data_highlight) > 1
  if (length(x = highlight_color) == 1 && length(x = found_meta_highlight) > 1) {
    highlight_color <- rep(x = highlight_color, length(x = found_meta_highlight))
    cli_inform(message = c("NOTE: Only one color provided to but {length(x = found_meta_highlight) `meta_data_highlight` variables were provided.}",
                           "i" = "Using the same color ({highlight_color[1]}) for all variables"))
  }

  # If NULL set using scCustomize_Palette
  if (is.null(x = highlight_color)) {
    highlight_color <- scCustomize_Palette(num_groups = length(x = cells_to_highlight), ggplot_default_colors = ggplot_default_colors)
  }

  # plot
  plot <- DimPlot_scCustom(seurat_object = seurat_object,
          cells.highlight = cells_to_highlight,
          cols.highlight = highlight_color,
          colors_use = background_color,
          sizes.highlight = pt.size,
          pt.size = pt.size,
          order = TRUE,
          raster = raster,
          raster.dpi = raster.dpi,
          split.by = split.by,
          split_seurat = split_seurat,
          label = label,
          ...)

  # Update legend and return plot
  plot <- suppressMessages(plot & scale_color_manual(breaks = names(x = cells_to_highlight), values = c(highlight_color, background_color), na.value = background_color))

  # Aspect ratio changes
  if (!is.null(x = aspect_ratio)) {
    if (!is.numeric(x = aspect_ratio)) {
      cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
    }
    plot <- plot & theme(aspect.ratio = aspect_ratio)
  }

  # Figure plot
  if (isTRUE(x = figure_plot)) {
    plot <- Figure_Plot(plot = plot)
  }

  return(plot)
}


#' Meta Highlight Plot
#'
#' Create Plot with meta data variable of interest highlighted
#'
#' @param seurat_object Seurat object name.
#' @param cells_highlight Cell names to highlight in named list.
#' @param highlight_color Color to highlight cells.
#' @param background_color non-highlighted cell colors (default is "lightgray")..
#' @param pt.size point size for both highlighted cluster and background.
#' @param aspect_ratio Control the aspect ratio (y:x axes ratio length).  Must be numeric value;
#' Default is NULL.
#' @param figure_plot logical.  Whether to remove the axes and plot with legend on left of plot denoting
#' axes labels.  (Default is FALSE).  Requires `split_seurat = TRUE`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param label Whether to label the highlighted meta data variable(s).  Default is FALSE.
#' @param split.by Variable in `@meta.data` to split the plot by.
#' @param split_seurat logical.  Whether or not to display split plots like Seurat (shared y axis) or as
#' individual plots in layout.  Default is FALSE.
#' @param ggplot_default_colors logical.  If `highlight_color = NULL`, Whether or not to return plot
#' using default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param ... Extra parameters passed to\code{\link[Seurat]{DimPlot}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @import patchwork
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' library(Seurat)
#'
#' # Creating example non-overlapping vectors of cells
#' MS4A1 <- WhichCells(object = pbmc_small, expression = MS4A1 > 4)
#' GZMB <- WhichCells(object = pbmc_small, expression = GZMB > 4)
#'
#' # Format as named list
#' cells <- list("MS4A1" = MS4A1,
#'               "GZMB" = GZMB)
#'
#' Cell_Highlight_Plot(seurat_object = pbmc_small, cells_highlight = cells)
#'

Cell_Highlight_Plot <- function(
  seurat_object,
  cells_highlight,
  highlight_color = NULL,
  background_color = "lightgray",
  pt.size = NULL,
  aspect_ratio = NULL,
  figure_plot = FALSE,
  raster = NULL,
  raster.dpi = c(512, 512),
  label = FALSE,
  split.by = NULL,
  split_seurat = FALSE,
  ggplot_default_colors = FALSE,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  if (!inherits(x = cells_highlight, what = "list")) {
    cli_abort(message = "{.code cells_highlight} must be of class: {.val list()}.")
  }

  if (is.null(x = names(x = cells_highlight))) {
    cli_abort(message = "Entries in {.code cells_highlight} list must be named.")
  }

  # Check duplicates
  if (any(duplicated(x = unlist(x = cells_highlight)))) {
    cli_abort(message = c("The list of {.code cells_highlight} contains duplicate cell names.",
                          "i" = "Ensure all cell names are unique before plotting."
                          )
              )
  }

  # Check all cells are present in object
  if (!all(unlist(x = cells_highlight) %in% Cells(x = seurat_object))) {
    cli_abort(message = c("Some of cells in {.code cells_highlight} are not present in object.",
                          "i" = "Ensure all cells are present in object before plotting."
                          )
              )
  }

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)

  # set point size
  if (is.null(x = pt.size)) {
    pt.size <- AutoPointSize_scCustom(data = sum(lengths(x = cells_highlight)), raster = raster)
  }

  # Check right number of colors provided
  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = highlight_color) && isTRUE(x = ggplot_default_colors)) {
    cli_abort(message = "Cannot provide both {.code highlight_color} and specify {.code ggplot_default_colors = TRUE}.")
  }

  if (!is.null(x = highlight_color)) {
    if (length(x = highlight_color) != length(x = cells_highlight)) {
      cli_abort(message = c("Incorrect number of highlight colors provided. Number of colors and groups must be equal.",
                            "i" = "{.code cells_highlight} contains: {.field {length(x = cells_highlight)}} groups but {.code highlight_color} contains: {.field {length(x = highlight_color)}} colors."
                            )
                )
    }
  } else {
    highlight_color <- scCustomize_Palette(num_groups = length(x = cells_highlight), ggplot_default_colors = ggplot_default_colors)
  }

  # plot
  plot <- DimPlot_scCustom(seurat_object = seurat_object,
                           cells.highlight = cells_highlight,
                           cols.highlight = highlight_color,
                           colors_use = background_color,
                           sizes.highlight = pt.size,
                           pt.size = pt.size,
                           order = TRUE,
                           raster = raster,
                           raster.dpi = raster.dpi,
                           split.by = split.by,
                           split_seurat = split_seurat,
                           label = label,
                           ...)

  # Edit plot legend
  plot <- suppressMessages(plot & scale_color_manual(breaks = names(x = cells_highlight), values = c(highlight_color, background_color), na.value = background_color))

  # Aspect ratio changes
  if (!is.null(x = aspect_ratio)) {
    if (!is.numeric(x = aspect_ratio)) {
      cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
    }
    plot <- plot & theme(aspect.ratio = aspect_ratio)
  }

  # Figure plot
  if (isTRUE(x = figure_plot)) {
    plot <- Figure_Plot(plot = plot)
  }

  return(plot)
}





#' DimPlot with modified default settings
#'
#' Creates DimPlot with some of the settings modified from their Seurat defaults (colors_use, shuffle, label).
#'
#' @param seurat_object Seurat object name.
#' @param colors_use color palette to use for plotting.  By default if number of levels plotted is less than
#' or equal to 36 it will use "polychrome" and if greater than 36 will use "varibow" with shuffle = TRUE
#' both from `DiscretePalette_scCustomize`.
#' @param pt.size Adjust point size for plotting.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default).
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident);
#' default is the current active.ident of the object.
#' @param split.by Feature to split plots by (i.e. "orig.ident").
#' @param split_seurat logical.  Whether or not to display split plots like Seurat (shared y axis) or as
#' individual plots in layout.  Default is FALSE.
#' @param figure_plot logical.  Whether to remove the axes and plot with legend on left of plot denoting
#' axes labels.  (Default is FALSE).  Requires `split_seurat = TRUE`.
#' @param aspect_ratio Control the aspect ratio (y:x axes ratio length).  Must be numeric value;
#' Default is NULL.
#' @param shuffle logical. Whether to randomly shuffle the order of points. This can be useful for crowded
#' plots if points of interest are being buried. (Default is TRUE).
#' @param seed Sets the seed if randomly shuffling the order of points.
#' @param label Whether to label the clusters.  By default if `group.by = NULL` label = TRUE, and
#' otherwise it is FALSE.
#' @param label.size Sets size of labels.
#' @param label.color Sets the color of the label text.
#' @param label.box Whether to put a box around the label text (geom_text vs geom_label).
#' @param dims Which dimensions to plot.  Defaults to c(1,2) if not specified.
#' @param repel Repel labels.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param num_columns Number of columns in plot layout.  Only valid if `split.by != NULL`.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ... Extra parameters passed to \code{\link[Seurat]{DimPlot}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @import patchwork
#' @importFrom Seurat DimPlot
#' @importFrom SeuratObject DefaultDimReduc
#'
#' @export
#'
#' @references Many of the param names and descriptions are from Seurat to facilitate ease of use as
#' this is simply a wrapper to alter some of the default parameters \url{https://github.com/satijalab/seurat/blob/master/R/visualization.R} (License: GPL-3).
#' `figure_plot` parameter/code modified from code by Tim Stuart via twitter: \url{https://twitter.com/timoast/status/1526237116035891200?s=20&t=foJOF81aPSjr1t7pk1cUPg}.
#'
#' @concept seurat_plotting
#'
#' @examples
#' library(Seurat)
#' DimPlot_scCustom(seurat_object = pbmc_small)
#'

DimPlot_scCustom <- function(
  seurat_object,
  colors_use = NULL,
  pt.size = NULL,
  reduction = NULL,
  group.by = NULL,
  split.by = NULL,
  split_seurat = FALSE,
  figure_plot = FALSE,
  aspect_ratio = NULL,
  shuffle = TRUE,
  seed = 1,
  label = NULL,
  label.size = 4,
  label.color = 'black',
  label.box = FALSE,
  dims = c(1, 2),
  repel = FALSE,
  raster = NULL,
  raster.dpi = c(512, 512),
  num_columns = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Change label if label.box
  if (isTRUE(x = label.box) && is.null(x = label)) {
    label <- TRUE
  }

  if (!is.null(x = split.by)) {
    split.by <- Meta_Present(seurat_object = seurat_object, meta_col_names = split.by, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  # Add check for group.by before getting to colors
  if (length(x = group.by) > 1) {
    Meta_Present(seurat_object = seurat_object, meta_col_names = group.by, print_msg = FALSE)
  } else {
    if (!is.null(x = group.by) && group.by != "ident") {
      Meta_Present(seurat_object = seurat_object, meta_col_names = group.by, print_msg = FALSE)
    }
  }

  # Add one time split_seurat warning
  if (!is.null(x = split.by) && isFALSE(x = split_seurat) && getOption(x = 'scCustomize_warn_DimPlot_split_type', default = TRUE)) {
    cli_inform(c("",
                 "NOTE: {.field DimPlot_scCustom} returns split plots as layout of all plots each",
                 "with their own axes as opposed to Seurat which returns with shared x or y axis.",
                 "To return to Seurat behvaior set {.code split_seurat = TRUE}.",
                 "",
                 "-----This message will be shown once per session.-----"))
    options(scCustomize_warn_DimPlot_split_type = FALSE)
  }

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)

  label <- label %||% (is.null(x = group.by))

  # if split.by is null set split_seurat to TRUE
  if (is.null(x = split.by)) {
    split_seurat <- TRUE
  }

  # # figure_plot check
  # if (isTRUE(x = figure_plot) && isTRUE(x = split_seurat)) {
  #   cli_abort(message = "{.code figure_plot} can only be TRUE if {.code split_seurat = FALSE}.")
  # }

  # Set default color palette based on number of levels being plotted
  if (length(x = group.by) > 1) {
    all_length <- lapply(group.by, function(x) {
      num_var <- length(x = unique(x = seurat_object@meta.data[[x]]))
    })
    group_by_length <- max(unlist(x = all_length))
  } else {
    if (is.null(x = group.by)) {
      group_by_length <- length(x = unique(x = seurat_object@active.ident))
    } else {
      group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
    }
  }

  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use) && isTRUE(x = ggplot_default_colors)) {
    cli_abort(message = "Cannot provide both custom palette to {.code colors_use} and specify {.code ggplot_default_colors = TRUE}.")
  }

  # set default plot colors
  if (is.null(x = colors_use)) {
    colors_use <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
  }

  # Set uniform point size is pt.size = NULL (based on plot with most cells)
  if (is.null(x = pt.size) && !is.null(split.by)) {
    # cells per meta data
    cells_by_split <- data.frame(table(seurat_object@meta.data[, split.by]))
    # Identity with greatest number of cells
    max_cells <- max(cells_by_split$Freq)
    # modified version of the autopointsize function from Seurat
    pt.size <- AutoPointSize_scCustom(data = max_cells, raster = raster)
  }

  # set size otherwise
  pt.size <- pt.size %||% AutoPointSize_scCustom(data = seurat_object)

  # Plot
  if (is.null(x = split.by)) {
    plot <- DimPlot(object = seurat_object, cols = colors_use, pt.size = pt.size, reduction = reduction, group.by = group.by, split.by = split.by, shuffle = shuffle, seed = seed, label = label, label.size = label.size, label.color = label.color, repel = repel, raster = raster, raster.dpi = raster.dpi, ncol = num_columns, dims = dims, label.box = label.box, ...)
    if (isTRUE(x = figure_plot)) {

      # pull axis labels
      x_lab_reduc <- plot$labels$x
      y_lab_reduc <- plot$labels$y

      plot <- plot & NoAxes()

      axis_plot <- ggplot(data.frame(x= 100, y = 100), aes(x = .data[["x"]], y = .data[["y"]])) +
        geom_point() +
        xlim(c(0, 10)) + ylim(c(0, 10)) +
        theme_classic() +
        ylab(y_lab_reduc) + xlab(x_lab_reduc) +
        theme(plot.background = element_rect(fill = "transparent", colour = NA),
              panel.background = element_rect(fill = "transparent"),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_line(
                arrow = arrow(angle = 15, length = unit(.5, "cm"), type = "closed")
              )
        )

      figure_layout <- c(
        area(t = 1, l = 2, b = 11, r = 11),
        area(t = 10, l = 1, b = 12, r = 2))

      plot_figure <- plot + axis_plot +
        plot_layout(design = figure_layout)

      # Aspect ratio changes
      if (!is.null(x = aspect_ratio)) {
        if (!is.numeric(x = aspect_ratio)) {
          cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
        }
        plot_figure <- plot_figure & theme(aspect.ratio = aspect_ratio)
      }

      return(plot_figure)
    } else {
      # Aspect ratio changes
      if (!is.null(x = aspect_ratio)) {
        if (!is.numeric(x = aspect_ratio)) {
          cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
        }
        plot <- plot & theme(aspect.ratio = aspect_ratio)
      }

      return(plot)
    }

  } else {
    if (isTRUE(x = split_seurat)) {
      # Plot Seurat Splitting
      plot <- DimPlot(object = seurat_object, cols = colors_use, pt.size = pt.size, reduction = reduction, group.by = group.by, split.by = split.by, shuffle = shuffle, seed = seed, label = label, label.size = label.size, label.color = label.color, repel = repel, raster = raster, raster.dpi = raster.dpi, ncol = num_columns, dims = dims, label.box = label.box, ...)
      if (isTRUE(x = figure_plot)) {

        # pull axis labels
        x_lab_reduc <- plot$labels$x
        y_lab_reduc <- plot$labels$y

        plot <- plot & NoAxes()

        axis_plot <- ggplot(data.frame(x= 100, y = 100), aes(x = .data[["x"]], y = .data[["y"]])) +
          geom_point() +
          xlim(c(0, 10)) + ylim(c(0, 10)) +
          theme_classic() +
          ylab(y_lab_reduc) + xlab(x_lab_reduc) +
          theme(plot.background = element_rect(fill = "transparent", colour = NA),
                panel.background = element_rect(fill = "transparent"),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_line(
                  arrow = arrow(angle = 15, length = unit(.5, "cm"), type = "closed")
                )
          )

        figure_layout <- c(
          area(t = 1, l = 2, b = 11, r = 11),
          area(t = 10, l = 1, b = 12, r = 2))

        plot_figure <- plot + axis_plot +
          plot_layout(design = figure_layout)

        # Aspect ratio changes
        if (!is.null(x = aspect_ratio)) {
          if (!is.numeric(x = aspect_ratio)) {
            cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
          }
          plot_figure <- plot_figure & theme(aspect.ratio = aspect_ratio)
        }

        return(plot_figure)
      } else {
        # Aspect ratio changes
        if (!is.null(x = aspect_ratio)) {
          if (!is.numeric(x = aspect_ratio)) {
            cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
          }
          plot <- plot & theme(aspect.ratio = aspect_ratio)
        }

        return(plot)
      }
    } else {
      if (is.null(x = group.by)) {
        group_by_vars <- as.character(unique(x = seurat_object@active.ident))
      } else {
        group_by_vars <- as.character(unique(x = seurat_object@meta.data[[group.by]]))
      }
      # Extract reduction coordinates
      reduction <- reduction %||% DefaultDimReduc(object = seurat_object)
      all_cells <- Cells(x = seurat_object)
      reduc_coordinates <- Embeddings(object = seurat_object[[reduction]])[all_cells, dims]
      reduc_coordinates <- as.data.frame(x = reduc_coordinates)
      x_axis <- c(min(reduc_coordinates[, 1]),
                  max(reduc_coordinates[, 1]))
      y_axis <- c(min(reduc_coordinates[, 2]),
                  max(reduc_coordinates[, 2]))

      # Extract cell names per meta data list of values
      # Extract split.by list of values
      if (inherits(x = seurat_object@meta.data[, split.by], what = "factor")) {
        split_by_list <- as.character(x = levels(x = seurat_object@meta.data[, split.by]))
      } else {
        split_by_list <- as.character(x = unique(x = seurat_object@meta.data[, split.by]))
      }

      cell_names <- lapply(split_by_list, function(x) {
        row.names(x = seurat_object@meta.data)[which(seurat_object@meta.data[, split.by] == x)]})

      # Unify colors across plots
      if (is.null(x = group.by)) {
        levels_overall <- levels(x = Idents(object = seurat_object))
      } else {
        levels_overall <- levels(x = seurat_object@meta.data[[group.by]])
      }

      colors_overall <- colors_use

      names(x = colors_overall) <- levels_overall

      # plot
      plots <- lapply(1:length(x = split_by_list), function(x) {
        plot <- DimPlot(object = seurat_object, cells = cell_names[[x]], group.by = group.by, cols = colors_use, reduction = reduction, pt.size = pt.size, raster = raster, raster.dpi = raster.dpi, shuffle = shuffle, seed = seed, label = label, label.size = label.size, label.color = label.color, repel = repel, dims = dims, label.box = label.box, ...) +
          ggtitle(paste(split_by_list[[x]])) +
          theme(plot.title = element_text(hjust = 0.5),
                legend.position = "right") +
          xlim(x_axis) +
          ylim(y_axis)

        # Normalize the colors across all plots
        plot <- suppressMessages(plot + scale_color_manual(values = colors_overall, drop = FALSE))

        if (!is.null(x = group.by)) {
          plot <- plot + labs(color=group.by)
        } else {
          plot <- plot
        }
      })

      # Wrap Plots into single output
      plots <- wrap_plots(plots, ncol = num_columns) + plot_layout(guides = 'collect')

      # Aspect ratio changes
      if (!is.null(x = aspect_ratio)) {
        if (!is.numeric(x = aspect_ratio)) {
          cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
        }
        plots <- plots & theme(aspect.ratio = aspect_ratio)
      }

      return(plots)
    }
  }
}


#' DimPlot by Meta Data Column
#'
#' Creates DimPlot layout containing all samples within Seurat Object from orig.ident column
#'
#' @param seurat_object Seurat object name.
#' @param meta_data_column Meta data column to split plots by.
#' @param colors_use single color to use for all plots or a vector of colors equal to the number of plots.
#' @param pt.size Adjust point size for plotting.
#' @param aspect_ratio Control the aspect ratio (y:x axes ratio length).  Must be numeric value;
#' Default is NULL.
#' @param title_size size for plot title labels.
#' @param num_columns number of columns in final layout plot.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default).
#' @param dims Which dimensions to plot.  Defaults to c(1,2) if not specified.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param ... Extra parameters passed to \code{\link[Seurat]{DimPlot}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import patchwork
#' @importFrom Seurat DimPlot
#' @importFrom SeuratObject DefaultDimReduc
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' library(Seurat)
#'
#' pbmc_small$sample_id <- sample(c("sample1", "sample2"), size = ncol(pbmc_small), replace = TRUE)
#'
#' DimPlot_All_Samples(seurat_object = pbmc_small, meta_data_column = "sample_id", color = "black",
#' num_columns = 2)
#'

DimPlot_All_Samples <- function(
  seurat_object,
  meta_data_column = "orig.ident",
  colors_use = "black",
  pt.size = NULL,
  aspect_ratio = NULL,
  title_size = 15,
  num_columns = NULL,
  reduction = NULL,
  dims = c(1, 2),
  raster = NULL,
  raster.dpi = c(512, 512),
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check meta data column exists
  if (meta_data_column %in% colnames(seurat_object@meta.data) == FALSE) {
    cli_abort(message = c("Meta data variable not found.",
                          "i" = "The {.code meta_data_column}: {.val {meta_data_column}}, could not be found in object@meta.data.")
    )
  }

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)

  # Extract reduction coordinates
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)
  all_cells <- Cells(x = seurat_object)
  reduc_coordinates <- Embeddings(object = seurat_object[[reduction]])[all_cells, dims]
  reduc_coordinates <- as.data.frame(x = reduc_coordinates)
  x_axis <- c(min(reduc_coordinates[, 1]),
              max(reduc_coordinates[, 1]))
  y_axis <- c(min(reduc_coordinates[, 2]),
              max(reduc_coordinates[, 2]))

  # Extract meta_data_column list of values
  if (inherits(x = seurat_object@meta.data[, meta_data_column], what = "factor")) {
    meta_sample_list <- levels(x = seurat_object@meta.data[, meta_data_column])
  } else {
    meta_sample_list <- as.character(x = unique(x = seurat_object@meta.data[, meta_data_column]))
  }

  # Extract cell names per meta data list of values
  cell_names <- lapply(meta_sample_list, function(x) {
    row.names(x = seurat_object@meta.data)[which(seurat_object@meta.data[, meta_data_column] == x)]})

  # Set uniform point size is pt.size = NULL (based on plot with most cells)
  if (is.null(x = pt.size)) {
    # cells per meta data
    cells_by_meta <- data.frame(table(seurat_object@meta.data[, meta_data_column]))
    # Identity with greatest number of cells
    max_cells <- max(cells_by_meta$Freq)
    # modified version of the autopointsize function from Seurat
    pt.size <- AutoPointSize_scCustom(data = max_cells, raster = raster)
  }

  # Check colors use
  num_plots <- length(x = meta_sample_list)

  if (length(x = colors_use) == 1) {
    colors_use <- rep(x = colors_use, num_plots)
  } else if (length(x = colors_use) != num_plots) {
    cli_abort(message = c("Not enough colors provided.",
                          "i" = "Length of {.code colors_use} ({.field {length(x = colors_use)}}) does not equal number of plots ({.field {num_plots}}).")
    )
  }

  # Plots
  plots <- lapply(1:length(x = meta_sample_list), function(x) {
    plot <- DimPlot(seurat_object, cells = cell_names[[x]], group.by = meta_data_column, cols = colors_use[[x]], reduction = reduction, pt.size = pt.size, raster = raster, raster.dpi = raster.dpi, ...) +
      ggtitle(paste(meta_sample_list[[x]])) +
      theme(plot.title = element_text(hjust = 0.5, size = title_size),
            legend.position = "none") +
      xlim(x_axis) +
      ylim(y_axis)
  })

  # Wrap Plots into single output
  plot_comb <- wrap_plots(plots, ncol = num_columns)

  # Aspect ratio changes
  if (!is.null(x = aspect_ratio)) {
    if (!is.numeric(x = aspect_ratio)) {
      cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
    }
    plot_comb <- plot_comb & theme(aspect.ratio = aspect_ratio)
  }

  return(plot_comb)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### OTHER PLOTTING ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Custom Labeled Variable Features Plot
#'
#' Creates variable features plot with N number of features already labeled by default.
#'
#' @param seurat_object Seurat object name.
#' @param num_features Number of top variable features to highlight by color/label.
#' @param custom_features A vector of custom feature names to label on plot instead of labeling top
#' variable genes.
#' @param label logical. Whether to label the top features.  Default is TRUE.
#' @param pt.size Adjust point size for plotting.
#' @param colors_use colors to use for plotting.  Default is "black" and "red".
#' @param repel logical (default TRUE).  Whether or not to repel the feature labels on plot.
#' @param y_axis_log logical. Whether to change y axis to log10 scale (Default is FALSE).
#' @param assay Assay to pull variable features from.
#' @param selection.method If more then one method use to calculate variable features specify which
#' method to use for plotting.  See `selection.method` parameter in \code{\link[Seurat]{VariableFeaturePlot}}
#' for list of options.
#' @param ... Extra parameters passed to \code{\link[Seurat]{VariableFeaturePlot}}.
#'
#' @return A ggplot object
#'
#' @import patchwork
#' @importFrom Seurat VariableFeaturePlot
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' library(Seurat)
#' VariableFeaturePlot_scCustom(seurat_object = pbmc_small, num_features = 10)
#'

VariableFeaturePlot_scCustom <- function(
  seurat_object,
  num_features = 10,
  custom_features = NULL,
  label = TRUE,
  pt.size = 1,
  colors_use = c("black", "red"),
  repel = TRUE,
  y_axis_log = FALSE,
  assay = NULL,
  selection.method = NULL,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # set assay (if null set to active assay)
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Extract num of desired features
  top_features <- head(x = VariableFeatures(object = seurat_object, assay = assay, selection.method = selection.method), num_features)

  # Plot
  plot <- VariableFeaturePlot(object = seurat_object, pt.size = pt.size, assay = assay, selection.method = selection.method, cols = colors_use, ...)

  # Label points
  if (isFALSE(x = label) && !is.null(x = custom_features)) {
    cli_warn(message = "The provided values provided to {.field custom_features} were not labeled as {.code label = FALSE} was also set.")
  }

  if (isTRUE(x = label)) {
    if (is.null(x = custom_features)) {
      plot <- LabelPoints(plot = plot, points = top_features, repel = repel)
    } else {
      # check all custom features are present
      all_found_features <- Feature_PreCheck(object = seurat_object, features = custom_features)
      plot <- LabelPoints(plot = plot, points = all_found_features, repel = repel)
    }
  }

  # return log10 y axis
  if (isTRUE(x = y_axis_log)) {
    plot <- plot + scale_y_log10()
    return(plot)
  }

  # Return plot
  return(plot)
}


#' Modified version of FeatureScatter
#'
#' Create customized FeatureScatter plots with scCustomize defaults.
#'
#' @param seurat_object Seurat object name.
#' @param feature1 First feature to plot.
#' @param feature2 Second feature to plot.
#' @param colors_use color for the points on plot.
#' @param pt.size Adjust point size for plotting.
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident).
#' Default is active ident.
#' @param split.by Feature to split plots by (i.e. "orig.ident").
#' @param split_seurat logical.  Whether or not to display split plots like Seurat (shared y axis) or as
#' individual plots in layout.  Default is FALSE.
#' @param shuffle logical, whether to randomly shuffle the order of points. This can be useful for crowded plots if points of interest are being buried. Default is TRUE.
#' @param aspect_ratio Control the aspect ratio (y:x axes ratio length).  Must be numeric value;
#' Default is NULL.
#' @param title_size size for plot title labels. Does NOT apply if `split_seurat = TRUE`.
#' @param plot.cor Display correlation in plot subtitle (or title if `split_seurat = TRUE`).
#' @param num_columns number of columns in final layout plot.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ... Extra parameters passed to \code{\link[Seurat]{FeatureScatter}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @import patchwork
#' @importFrom magrittr "%>%"
#' @importFrom Seurat FeatureScatter
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' \donttest{
#' library(Seurat)
#'
#' pbmc_small$sample_id <- sample(c("sample1", "sample2"), size = ncol(pbmc_small), replace = TRUE)
#'
#' FeatureScatter_scCustom(seurat_object = pbmc_small, feature1 = "nCount_RNA",
#' feature2 = "nFeature_RNA", split.by = "sample_id")
#' }
#'

FeatureScatter_scCustom <- function(
    seurat_object,
    feature1 = NULL,
    feature2 = NULL,
    colors_use = NULL,
    pt.size = NULL,
    group.by = NULL,
    split.by = NULL,
    split_seurat = FALSE,
    shuffle = TRUE,
    aspect_ratio = NULL,
    title_size = 15,
    plot.cor = TRUE,
    num_columns = NULL,
    raster = NULL,
    raster.dpi = c(512, 512),
    ggplot_default_colors = FALSE,
    color_seed = 123,
    ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  if (!is.null(x = split.by)) {
    split.by <- Meta_Present(seurat_object = seurat_object, meta_col_names = split.by, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  # Add check for group.by before getting to colors
  if (length(x = group.by) > 1) {
    Meta_Present(seurat_object = seurat_object, meta_col_names = group.by, print_msg = FALSE)
  } else {
    if (!is.null(x = group.by) && group.by != "ident") {
      Meta_Present(seurat_object = seurat_object, meta_col_names = group.by, print_msg = FALSE)
    }
  }

  # Add one time split_seurat warning
  if (!is.null(x = split.by) && isFALSE(x = split_seurat) && getOption(x = 'scCustomize_warn_FeatureScatter_split_type', default = TRUE)) {
    cli_inform(c("",
                 "NOTE: {.field FeatureScatter_scCustom} returns split plots as layout of all plots each",
                 "with their own axes as opposed to Seurat which returns with shared x or y axis.",
                 "To return to Seurat behvaior set {.code split_seurat = TRUE}.",
                 "",
                 "-----This message will be shown once per session.-----"))
    options(scCustomize_warn_FeatureScatter_split_type = FALSE)
  }

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)

  # Set default color palette based on number of levels being plotted
  if (length(x = group.by) > 1) {
    all_length <- lapply(group.by, function(x) {
      num_var <- length(x = unique(x = seurat_object@meta.data[[x]]))
    })
    group_by_length <- max(unlist(x = all_length))
  } else {
    if (is.null(x = group.by)) {
      group_by_length <- length(x = unique(x = seurat_object@active.ident))
    } else {
      group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
    }
  }

  # set default plot colors
  if (is.null(x = colors_use)) {
    colors_use <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
  }

  # Set uniform point size is pt.size = NULL (based on plot with most cells)
  if (is.null(x = pt.size) && !is.null(split.by)) {
    # cells per meta data
    cells_by_split <- data.frame(table(seurat_object@meta.data[, split.by]))
    # Identity with greatest number of cells
    max_cells <- max(cells_by_split$Freq)
    # modified version of the autopointsize function from Seurat
    pt.size <- AutoPointSize_scCustom(data = max_cells, raster = raster)
  }

  # set size otherwise
  pt.size <- pt.size %||% AutoPointSize_scCustom(data = seurat_object)

  # Plot
  if (is.null(x = split.by)) {
    plot <- FeatureScatter(object = seurat_object, feature1 = feature1, feature2 = feature2, cols = colors_use, pt.size = pt.size, group.by = group.by, split.by = split.by, shuffle = shuffle, plot.cor = plot.cor, raster = raster, raster.dpi = raster.dpi, ncol = num_columns, ...)

    # Change title
    plot <- plot +
      theme(plot.title = element_text(hjust = 0.5, size = title_size), legend.position = "right") +
      ggtitle(paste0(feature1, " vs. ", feature2), subtitle = paste0("Correlation: ", plot$labels$title))


    # Aspect ratio changes
    if (!is.null(x = aspect_ratio)) {
      if (!is.numeric(x = aspect_ratio)) {
        cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
      }
      plot <- plot & theme(aspect.ratio = aspect_ratio)
    }

    # return plot
    return(plot)
  } else {
    # Plot with Seurat splitting
    if (isTRUE(x = split_seurat)) {
      plot <- FeatureScatter(object = seurat_object, feature1 = feature1, feature2 = feature2, cols = colors_use, pt.size = pt.size, group.by = group.by, split.by = split.by, shuffle = shuffle, plot.cor = plot.cor, raster = raster, raster.dpi = raster.dpi, ncol = num_columns, ...)

      # Aspect ratio changes
      if (!is.null(x = aspect_ratio)) {
        if (!is.numeric(x = aspect_ratio)) {
          cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
        }
        plot <- plot & theme(aspect.ratio = aspect_ratio)
      }

      # return plot
      return(plot)
    } else {
      plot <- scCustomze_Split_FeatureScatter(seurat_object = seurat_object, feature1 = feature1, feature2 = feature2, split.by = split.by, group.by = group.by, colors_use = colors_use, pt.size = pt.size, aspect_ratio = aspect_ratio, title_size = title_size, num_columns = num_columns, raster = raster, raster.dpi = raster.dpi, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, ...)

      return(plot)
    }
  }
}
