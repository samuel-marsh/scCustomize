#' Customize FeaturePlot
#'
#' Create Custom FeaturePlots and preserve scale (no binning)
#'
#' @param seurat_object Seurat object name.
#' @param features Feature(s) to plot.
#' @param colors_use list of colors or color palette to use.
#' @param na_value color to use for points below lower limit.
#' @param order whether to move positive cells to the top (default = TRUE).
#' @param pt.size Adjust point size for plotting.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default).
#' @param na_cutoff Value to use as minimum expression cutoff.  This will be lowest value plotted use
#' palette provided to `colors_use`.  Leave as default value to plot only positive non-zero values using
#' color scale and zero/negative values as NA.  To plot all values using color palette set to `NA`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param split.by Variable in `@meta.data` to split the plot by.
#' @param num_columns Number of columns in plot layout.
#' @param slot Which slot to pull expression data from?  Default is "data".
#' @param ... Extra parameters passed to \code{\link[Seurat]{FeaturePlot}}.
#'
#' @return A ggplot object
#'
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
#' FeaturePlot_scCustom(seurat_object = object, features = "Cx3cr1", colors_use = viridis_plasma_dark_high,
#' na_color = "lightgray")
#' }
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
  split.by = NULL,
  num_columns = NULL,
  slot = "data",
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check meta
  if (!is.null(x = split.by)) {
    split.by <- Meta_Present(seurat_object = seurat_object, meta_col_names = split.by, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  # Check num_columns and split.by with feature length
  # if (!is.null(x = num_columns) && !is.null(x = split.by) && length(x = features) > 1) {
  #   stop("'num_columns' parameter cannot currently be used with split.by if number of features plotted is greater than 1.")
  # }

  # Get length of meta data feature
  if (!is.null(x = split.by)) {
    split.by_length <- length(unique(seurat_object@meta.data[[split.by]]))
    if (is.null(x = num_columns)) {
      num_columns <- split.by_length
    }
    # Calculate number of rows for selected number of columns
    num_rows <- ceiling(split.by_length/num_columns)

    # Check column and row compatibility
    if (num_columns > split.by_length) {
      stop("The number of columns specified is greater than the number of meta data variables.  ", paste0('"', split_by, '"', " only contains ", split.by_length, " variables.  "), "Please adjust `num_columns` to be less than or equal to", ": ", paste(split.by_length), ".")
    }
  }

  if (any(features %in% colnames(x = seurat_object@meta.data))) {
    warning("Some of the plotted features are from meta.data slot.  Please check that `na_cutoff`` param is being set appropriately for those features.")
  }

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = colnames(x = seurat_object)) > 2e5)

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

  # plot no split
  if (is.null(x = split.by)) {
    plot <- suppressMessages(FeaturePlot(object = seurat_object, features = features, order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, ncol = num_columns, ...) & scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, NA), na.value = na_color))
  }

  # plotting split with single feature (allows column number setting)
  if (!is.null(x = split.by) && length(x = features) == 1) {
    # Until Seurat is fixed pull feature data to set separately
    feature_data <- FetchData(
      object = seurat_object,
      vars = features,
      slot = slot)
    # Pull min and max values
    max_exp_value <- max(feature_data)
    min_exp_value <- min(feature_data)

    plot <- suppressMessages(FeaturePlot(object = seurat_object, features = features, order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, ...) & scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, max_exp_value), na.value = na_color, name = features)) & RestoreLegend() & theme(axis.title.y.right = element_blank())

    plot <- plot + plot_layout(nrow = num_rows, ncol = num_columns)
  }

  # plotting split multiple features
  if (!is.null(x = split.by) && length(x = features) > 1) {

    plot_list <- lapply(1:length(x = features), function(i){
      feature_data <- FetchData(
        object = seurat_object,
        vars = features[i],
        slot = slot)
      # Pull min and max values
      max_exp_value <- max(feature_data)
      min_exp_value <- min(feature_data)

      single_plot <- suppressMessages(FeaturePlot(object = seurat_object, features = features[i], order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, ...) & scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, max_exp_value), na.value = na_color, name = features[i])) & RestoreLegend() & theme(axis.title.y.right = element_blank())

      single_plot <- single_plot + plot_layout(nrow = num_rows, ncol = num_columns)
    })
    plot <- wrap_plots(plot_list) + plot_layout(ncol = 1)
  }

  # Add one time na_cutoff warning
  if (getOption(x = 'scCustomize_warn_na_cutoff', default = TRUE) && !is.na(x = na_cutoff) && na_cutoff == 0.000000001) {
    message(
      "NOTE: FeaturePlot_scCustom uses a specified `na_cutoff` when plotting to \n",
      "color cells with no expression differently.\n",
      "Please ensure `na_cutoff` value is appropriate for feature being plotted.\n",
      "Default setting is appropriate for use when plotting from 'RNA' assay.\n",
      "When `na_cutoff` not appropriate (e.g., module scores) set to NULL to \n",
      "plot all cells in gradiant color palette.
       \nThis message will be shown once per session.\n"
    )
    options(scCustomize_warn_na_cutoff = FALSE)
  }

  if (getOption(x = 'scCustomize_warn_zero_na_cutoff', default = TRUE) && !is.na(x = na_cutoff) && na_cutoff == 0) {
    message(
      "NOTE: Specified `na_cutoff` is set to zero. This means that only cells/nuclei\n",
      "with expression less than zero will be plotted with `na_color`.\n",
      "To plot cells with expression values of zero using `na_color` leave \n",
      "default `na_cutoff` value. If you want to plot full spectrum without \n",
      "`na_cutoff` (e.g., for module scores) then set `na_cutoff = NULL`.
       \nThis message will be shown once per session.\n"
    )
    options(scCustomize_warn_na_cutoff = FALSE)
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
#' @param na_value color to use for points below lower limit.
#' @param order whether to move positive cells to the top (default = TRUE).
#' @param pt.size Adjust point size for plotting.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default).
#' @param na_cutoff Value to use as minimum expression cutoff.  To set no cutoff set to `NA`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param slot Which slot to pull expression data from?  Default is "data".
#' @param num_columns Number of columns in plot layout.  If number of features > 1 then `num_columns`
#' dictates the number of columns in overall layout (`num_columns = 1` means stacked layout & `num_columns = 2`
#' means adjacent layout).
#' @param ... Extra parameters passed to \code{\link[Seurat]{FeaturePlot}}.
#'
#' @return A ggplot object
#'
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
  reduction = NULL,
  na_cutoff = 0.000000001,
  raster = NULL,
  slot = "data",
  num_columns = NULL,
  ...
) {
  # Check assays present
  assays_not_found <- Assay_Present(seurat_object = seurat_object, assay_list = c(assay1, assay2), print_msg = FALSE, omit_warn = TRUE)[[2]]

  if (!is.null(x = assays_not_found)) {
    stop_quietly()
  }

  # Find commands
  commands <- Command(object = seurat_object)

  # Raw normalize check
  if (!paste0("NormalizeData.", assay1) %in% commands) {
    stop("Assay 1: ", assay1, " has not been normalized.  Please run `NormalizeData` on this assay before proceeding to visualization.")
  }

  # Cell Bender normalize check
  if (!paste0("NormalizeData.", assay2) %in% commands) {
    stop("Assay 2: ", assay2, " has not been normalized.  Please run `NormalizeData` on this assay before proceeding to visualization.")
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
    warning("When plotting more than one feature `num_columns` refers to patchwork columns and must either be 1 (vertical) or 2 (horizontal).")
  }

  # Change assay and plot raw
  DefaultAssay(seurat_object) <- assay1

  plot_raw <- FeaturePlot_scCustom(seurat_object = seurat_object, features = features, slot = slot, colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, order = order, pt.size = pt.size, reduction = reduction, raster = raster, ...) & labs(color = assay1)

  # Change to cell bender and plot
  DefaultAssay(seurat_object) <- assay2

  plot_cell_bender <- FeaturePlot_scCustom(seurat_object = seurat_object, features = features, slot = slot, colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, order = order, pt.size = pt.size, reduction = reduction, raster = raster, ...) & labs(color = assay2)

  # Assemble plots & return plots
  plots <- wrap_plots(plot_raw, plot_cell_bender, ncol = num_columns)

  return(plots)
}


#' Split FeatureScatter
#'
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
#' @param title_size size for plot title labels.
#' @param num_columns number of columns in final layout plot.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 100,000 cells.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ... Extra parameters passed to \code{\link[Seurat]{FeatureScatter}}.
#'
#' @return A ggplot object
#'
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
#' \dontrun{
#' Split_FeatureScatter(seurat_object = object, features1= "nCount_RNA", feature2 = "nFeature_RNA",
#' split.by = "orig.ident", colors_use = "navy")
#' }
#'

Split_FeatureScatter <- function(
  seurat_object,
  feature1 = NULL,
  feature2 = NULL,
  split.by = NULL,
  group.by = NULL,
  colors_use = NULL,
  pt.size = NULL,
  title_size = 15,
  num_columns = NULL,
  raster = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # split.by present
  if (is.null(x = split.by)) {
    stop("No value supplied to `split.by`.")
  }

  # Set columna and row lengths
  split.by_length <- length(unique(seurat_object@meta.data[[split.by]]))

   if (is.null(x = num_columns)) {
    num_columns <- split.by_length
  }
  # Calculate number of rows for selected number of columns
  num_rows <- ceiling(split.by_length/num_columns)

  # Check column and row compatibility
  if (num_columns > split.by_length) {
    stop("The number of columns specified is greater than the number of meta data variables.  ", paste0('"', split.by, '"', " only contains ", split.by_length, " variables.  "), "Please adjust `num_columns` to be less than or equal to", ": ", paste(split.by_length), ".")
  }

  if (split.by %in% colnames(seurat_object@meta.data) == FALSE) {
    stop("The meta data variable: ", '"', split.by, '"', " could not be found in object@meta.data.")
  }

  # Check features are present
  possible_features <- c(rownames(seurat_object), colnames(seurat_object@meta.data))
  check_features <- setdiff(x = c(feature1, feature2), y = possible_features)
  if (length(x = check_features) > 0) {
    stop("Feature: ", '"', check_features, '"', " is not present in Seurat object.")
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
  meta_sample_list <- as.character(unique(seurat_object@meta.data[, split.by]))

  # Extract cell names per meta data list of values
  cell_names <- lapply(meta_sample_list, function(x) {
    row.names(seurat_object@meta.data)[which(seurat_object@meta.data[, split.by] == x)]})

  # raster check
  raster <- raster %||% (length(x = colnames(x = seurat_object)) > 2e5)

  # Set uniform poist size is pt.size = NULL (based on plot with most cells)
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

    cor_values <- lapply(1:length(x = meta_sample_list_test), function(i) {
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
    if (ggplot_default_colors) {
      colors_use <- Hue_Pal(num_colors = group_by_length)
    } else {
      if (group_by_length <= 36) {
        colors_use <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
      } else {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "varibow", shuffle_pal = TRUE, seed = color_seed)
      }
    }
  }

  # Plots
  plots <- lapply(1:length(x = meta_sample_list), function(j) {
    plot <- FeatureScatter(seurat_object, feature1 = feature1, feature2 = feature2, cells = cell_names[[j]], group.by = group.by, cols = colors_use, pt.size = pt.size, raster = raster, ...) +
      theme(plot.title = element_text(hjust = 0.5, size = title_size),
            legend.position = "right") +
      xlim(min_feature1, max_feature1) +
      ylim(min_feature2, max_feature2)
    if (plot_cor) {
      plot + ggtitle(paste(meta_sample_list[[j]]), subtitle = paste0("Correlation: ", cor_values[i]))
    } else {
      plot + ggtitle(paste(meta_sample_list[[j]]))
    }
  })

  # Wrap Plots into single output
  wrap_plots(plots, ncol = num_columns, nrow = num_rows) + plot_layout(guides = 'collect')
}


#' Cluster Highlight Plot
#'
#' Create Plot with cluster of interest highlighted
#'
#' @param seurat_object Seurat object name.
#' @param cluster_name Name (or number) identity of cluster to be highlighted.
#' @param highlight_color Color to highlight cells (default "navy").
#' @param background_color non-highlighted cell colors.
#' @param pt.size point size for both highlighted cluster and background.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param ... Extra parameters passed to\code{\link[Seurat]{DimPlot}}.
#'
#' @return A ggplot object
#'
#' @importFrom Seurat DimPlot
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' \dontrun{
#' Cluster_Highlight_Plot(seurat_object = object, cluster_name = "Microglia", highlight_color = "gold",
#' background_color = "lightgray",  pt.size = 2)
#' }
#'

Cluster_Highlight_Plot <- function(
  seurat_object,
  cluster_name,
  highlight_color = "navy",
  background_color = "lightgray",
  pt.size = NULL,
  raster = NULL,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = colnames(x = seurat_object)) > 2e5)

  # pull cells to highlight in plot
  cells_to_highlight <- CellsByIdentities(seurat_object, idents = cluster_name)

  # set point size
  if (is.null(x = pt.size)) {
    pt.size <- AutoPointSize_scCustom(data = sum(lengths(cells_to_highlight)), raster = raster)
  }

  # plot
  plot <- DimPlot(object = seurat_object,
          cells.highlight = cells_to_highlight,
          cols.highlight = highlight_color,
          cols = background_color,
          sizes.highlight = pt.size,
          pt.size = pt.size,
          order = TRUE,
          raster = raster,
          ...)

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
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param ... Extra parameters passed to\code{\link[Seurat]{DimPlot}}.
#'
#' @return A ggplot object
#'
#' @importFrom Seurat DimPlot
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' \dontrun{
#' Meta_Highlight_Plot(seurat_object = object, meta_data_column = "orig.ident", meta_data_highlight = "sample_01",
#' highlight_color = "gold", background_color = "lightgray",  pt.size = 2)
#' }
#'

Meta_Highlight_Plot <- function(
  seurat_object,
  meta_data_column,
  meta_data_highlight,
  highlight_color = "navy",
  background_color = "lightgray",
  pt.size = NULL,
  raster = NULL,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check meta data
  meta_data_column <- Meta_Present(seurat_object = seurat_object, meta_col_names = meta_data_column, omit_warn = FALSE, print_msg = FALSE)[[1]]

  # stop if none found
  if (length(x = meta_data_column) == 0) {
    stop("The 'meta_data_column': ", meta_data_column, " was not found in object meta.data slot.")
  }

  # Check that meta data is factor or character
  accepted_meta_types <- c("factor", "character", "logical")

  if (!class(x = seurat_object@meta.data[[meta_data_column]]) %in% accepted_meta_types) {
    stop("The 'meta_data_column': ", meta_data_column, " is of class: ", '"', class(x = seurat_object@meta.data[[meta_data_column]]), '"', " only meta data variables of classes: factor, character, or logical can be used with Meta_Highlight_Plot().")
  }

  # Check meta_data_highlight
  meta_var_list <- as.character(unique(seurat_object@meta.data[, meta_data_column]))

  if (!meta_data_highlight %in% meta_var_list) {
    stop("The 'meta_data_highlight': ", meta_data_highlight, " was not found in the meta.data column: ", meta_data_column, ".")
  }

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = colnames(x = seurat_object)) > 2e5)

  # Change default ident and pull cells to highlight in plot
  Idents(seurat_object) <- meta_data_column

  cells_to_highlight <- CellsByIdentities(seurat_object, idents = meta_data_highlight)

  # set point size
  if (is.null(x = pt.size)) {
    pt.size <- AutoPointSize_scCustom(data = sum(lengths(cells_to_highlight)), raster = raster)
  }

  # plot
  plot <- DimPlot(object = seurat_object,
          cells.highlight = cells_to_highlight,
          cols.highlight = highlight_color,
          cols = background_color,
          sizes.highlight = pt.size,
          pt.size = pt.size,
          order = TRUE,
          raster = raster,
          ...)

  return(plot)
}


#' Stacked Violin Plot
#'
#' Code for creating stacked violin plot gene expression.
#'
#' @param seurat_object Seurat object name.
#' @param features Features to plot.
#' @param x_lab_rotate Rotate x-axis labels 45 degrees (Default is FALSE).
#' @param colors_use specify color palette to used in \code{\link[Seurat]{VlnPlot}}.  By default if
#' number of levels plotted is less than or equal to 36 it will use "polychrome" and if greater than 36
#' will use "varibow" with shuffle = TRUE both from `DiscretePalette_scCustomize`.
#' @param group.by Group (color) cells in different ways (for example, orig.ident).
#' @param split.by A variable to split the violin plots by,
#' @param idents Which classes to include in the plot (default is all).
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param ... Extra parameters passed to \code{\link[Seurat]{VlnPlot}}.
#'
#' @return A ggplot object
#'
#' @import patchwork
#' @importFrom purrr map map_dbl map2
#' @importFrom Seurat VlnPlot
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @author Ming Tang (Original Code), Sam Marsh (Wrap single function, added/modified functionality)
#' @references https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
#' @seealso https://twitter.com/tangming2005
#'
#' @examples
#' \dontrun{
#' Stacked_VlnPlot(seurat_object = object, features = gene_list, x_lab_rotate = TRUE)
#' }
#'

Stacked_VlnPlot <- function(
  seurat_object,
  features,
  group.by = NULL,
  split.by = NULL,
  idents = NULL,
  x_lab_rotate = FALSE,
  colors_use = NULL,
  color_seed = 123,
  ggplot_default_colors = FALSE,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check features and meta to determine which features present
  features_list <- Gene_Present(data = seurat_object, gene_list = features, omit_warn = FALSE, print_msg = FALSE, case_check_msg = FALSE, return_none = TRUE)

  meta_list <- Meta_Present(seurat_object = seurat_object, meta_col_names = features_list[[2]], omit_warn = FALSE, print_msg = FALSE, abort = FALSE)

  all_not_found_features <- meta_list[[2]]

  all_found_features <- c(features_list[[1]], meta_list[[1]])

  # Return message of features not found
  if (length(x = all_not_found_features) > 0) {
    op <- options(warn = 1)
    on.exit(options(op))
    warning("The following features were omitted as they were not found",
            ": ", glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE))
  }

  # Check feature case correct
  Case_Check(seurat_object = seurat_object, gene_list = all_not_found_features, case_check_msg = TRUE, return_features = FALSE)

  # Set default color palette based on number of levels being plotted
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }

  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use) && ggplot_default_colors) {
    stop("Cannot provide both custom palette to `colors_use` and specify `ggplot_default_colors = TRUE`.")
  }
  if (is.null(x = colors_use)) {
    if (ggplot_default_colors) {
      colors_use <- Hue_Pal(num_colors = group_by_length)
    } else {
      if (group_by_length <= 36) {
        colors_use <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
      } else {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "varibow", shuffle_pal = TRUE, seed = color_seed)
      }
    }
  }

  # Create plots
  plot_list <- map(all_found_features, function(x) Modify_VlnPlot(seurat_object = seurat_object, features = x, cols = colors_use, group.by = group.by, split.by = split.by, idents = idents, ...))

  # Add back x-axis title to bottom plot. patchwork is going to support this?
  # Add ability to rotate the X axis labels to the function call
  if (x_lab_rotate) {
    plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.ticks.x = element_line())
  }
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(), axis.ticks.x = element_line())

  # change the y-axis tick to only max value
  ymaxs <- map_dbl(plot_list, Extract_Max)
  plot_list <- suppressMessages(map2(plot_list, ymaxs, function(x,y) x +
                                       scale_y_continuous(breaks = c(y)) +
                                       expand_limits(y = y)))

  wrap_plots(plotlist = plot_list, ncol = 1)
}


#' Customized DotPlot
#'
#' Code for creating customized DotPlot
#'
#' @param seurat_object Seurat object name.
#' @param features Features to plot.
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
#' \dontrun{
#' DotPlot_scCustom(seurat_object = object, features = gene_list)
#' }
#'

DotPlot_scCustom <- function(
  seurat_object,
  features,
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

  # Check Genes
  features_list <- Gene_Present(data = seurat_object, gene_list = features, omit_warn = TRUE, print_msg = FALSE, case_check_msg = TRUE)[[1]]

  # Plot
  plot <- suppressMessages(DotPlot(object = seurat_object, features = features_list, ...) +
                             scale_color_gradientn(colors = colors_use)
  )
  # Modify plot
  if (remove_axis_titles) {
    plot <- plot +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank()
    )
  }

  if (flip_axes) {
    plot <- plot & coord_flip()
  }
  if (x_lab_rotate) {
    plot <- plot +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  if (y_lab_rotate) {
    plot <- plot +
      theme(axis.text.y = element_text(angle = 45, hjust = 1))
  }

  if (is.list(x = features)) {
    plot + theme(strip.text.x = element_text(angle = 45))
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
#' pass 'ident' to group by identity class.
#' @param split.by Feature to split plots by (i.e. "orig.ident").
#' @param split_seurat logical.  Whether or not to display split plots like Seurat (shared y axis) or as
#' individual plots in layout.  Default is FALSE.
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
#' @param num_columns Number of columns in plot layout.  Only valid if `split.by != NULL`.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ... Extra parameters passed to \code{\link[Seurat]{DimPlot}}.
#'
#' @return A ggplot object
#'
#' @import patchwork
#' @importFrom Seurat DimPlot
#'
#' @export
#'
#' @references Many of the param names and descriptions are from Seurat to facilitate ease of use as
#' this is simply a wrapper to alter some of the default parameters (https://github.com/satijalab/seurat/blob/master/R/visualization.R) (Licence: GPL-3).
#'
#' @concept seurat_plotting
#'
#' @examples
#' \dontrun{
#' DimPlot_scCustom(seurat_object = object)
#' }
#'

DimPlot_scCustom <- function(
  seurat_object,
  colors_use = NULL,
  pt.size = NULL,
  reduction = NULL,
  group.by = NULL,
  split.by = NULL,
  split_seurat = FALSE,
  shuffle = TRUE,
  seed = 1,
  label = NULL,
  label.size = 4,
  label.color = 'black',
  label.box = FALSE,
  dims = c(1, 2),
  repel = FALSE,
  raster = NULL,
  num_columns = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  if (!is.null(x = split.by)) {
    split.by <- Meta_Present(seurat_object = seurat_object, meta_col_names = split.by, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  # Add one time split_seurat warning
  if (!is.null(x = split.by) && !split_seurat && getOption(x = 'scCustomize_warn_DimPlot_split_type', default = TRUE)) {
    message(
      "NOTE: DimPlot_scCustom returns split plots as layout of all plots each \n",
      "with their own axes as opposed to Seurat which returns with shared x or y axis.\n",
      "To return to Seurat behvaior set `split_seurat = TRUE`.
       \nThis message will be shown once per session.\n"
    )
    options(scCustomize_warn_DimPlot_split_type = FALSE)
  }

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = colnames(x = seurat_object)) > 2e5)

  label <- label %||% (is.null(x = group.by))

  # Set default color palette based on number of levels being plotted
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }

  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use) && ggplot_default_colors) {
    stop("Cannot provide both custom palette to `colors_use` and specify `ggplot_default_colors = TRUE`.")
  }

  # set default plot colors
  if (is.null(x = colors_use)) {
    if (ggplot_default_colors) {
      colors_use <- Hue_Pal(num_colors = group_by_length)
    } else {
      if (group_by_length <= 36) {
        colors_use <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
      } else {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "varibow", shuffle_pal = TRUE, seed = color_seed)
      }
    }
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
    DimPlot(object = seurat_object, cols = colors_use, pt.size = pt.size, reduction = reduction, group.by = group.by, split.by = split.by, shuffle = shuffle, seed = seed, label = label, label.size = label.size, label.color = label.color, repel = repel, raster = raster, ncol = num_columns, dims = dims, label.box = label.box, ...)
  } else {
    if (split_seurat) {
      # Plot Seurat Splitting
      DimPlot(object = seurat_object, cols = colors_use, pt.size = pt.size, reduction = reduction, group.by = group.by, split.by = split.by, shuffle = shuffle, seed = seed, label = label, label.size = label.size, label.color = label.color, repel = repel, raster = raster, ncol = num_columns, dims = dims, label.box = label.box, ...)
    } else {
      if (is.null(x = group.by)) {
        group_by_vars <- as.character(unique(x = seurat_object@active.ident))
      } else {
        group_by_vars <- as.character(unique(x = seurat_object@meta.data[[group.by]]))
      }
      # Extract reduction coordinates
      reduction <- reduction %||% DefaultDimReduc(object = seurat_object)
      all_cells <- colnames(x = seurat_object)
      reduc_coordinates <- Embeddings(object = seurat_object[[reduction]])[all_cells, dims]
      reduc_coordinates <- as.data.frame(x = reduc_coordinates)
      x_axis <- c(min(reduc_coordinates[, 1]),
                  max(reduc_coordinates[, 1]))
      y_axis <- c(min(reduc_coordinates[, 2]),
                  max(reduc_coordinates[, 2]))

      # Extract cell names per meta data list of values
      split_by_list <- as.character(unique(x = seurat_object@meta.data[[split.by]]))

      cell_names <- lapply(split_by_list, function(x) {
        row.names(seurat_object@meta.data)[which(seurat_object@meta.data[, split.by] == x)]})

      # plot
      plots <- lapply(1:length(x = split_by_list), function(x) {
        plot <- DimPlot(object = seurat_object, cells = cell_names[[x]], group.by = group.by, cols = colors_use, reduction = reduction, pt.size = pt.size, raster = raster, shuffle = shuffle, seed = seed, label = label, label.size = label.size, label.color = label.color, repel = repel, dims = dims, label.box = label.box, ...) +
          ggtitle(paste(split_by_list[[x]])) +
          theme(plot.title = element_text(hjust = 0.5),
                legend.position = "right") +
          xlim(x_axis) +
          ylim(y_axis)
        if (!is.null(x = group.by)) {
          plot <- plot + labs(color=group.by)
        } else {
          plot <- plot
        }
      })

      # Wrap Plots into single output
      wrap_plots(plots, ncol = num_columns) + plot_layout(guides = 'collect')
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
#' @param title_size size for plot title labels.
#' @param num_columns number of columns in final layout plot.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default).
#' @param dims Which dimensions to plot.  Defaults to c(1,2) if not specified.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param ... Extra parameters passed to \code{\link[Seurat]{DimPlot}}.
#'
#' @return A ggplot object
#'
#' @import patchwork
#' @importFrom Seurat DimPlot
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' \dontrun{
#' DimPlot_All_Samples(seurat_object = object, meta_data_column = "orig.ident", color = "black",
#' num_columns = 2, reduction = "tsne")
#' }
#'

DimPlot_All_Samples <- function(
  seurat_object,
  meta_data_column = "orig.ident",
  colors_use = "black",
  pt.size = NULL,
  title_size = 15,
  num_columns = NULL,
  reduction = NULL,
  dims = c(1, 2),
  raster = NULL,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check meta data column exists
  if (meta_data_column %in% colnames(seurat_object@meta.data) == FALSE) {
    stop("The meta data variable: ", '"', meta_data_column, '"', " could not be found in object@meta.data.")
  }

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = colnames(x = seurat_object)) > 2e5)

  # Extract reduction coordinates
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)
  all_cells <- colnames(x = seurat_object)
  reduc_coordinates <- Embeddings(object = seurat_object[[reduction]])[all_cells, dims]
  reduc_coordinates <- as.data.frame(x = reduc_coordinates)
  x_axis <- c(min(reduc_coordinates[, 1]),
              max(reduc_coordinates[, 1]))
  y_axis <- c(min(reduc_coordinates[, 2]),
              max(reduc_coordinates[, 2]))

  # Extract meta_data_column list of values
  meta_sample_list <- as.character(unique(seurat_object@meta.data[, meta_data_column]))

  # Extract cell names per meta data list of values
  cell_names <- lapply(meta_sample_list, function(x) {
    row.names(seurat_object@meta.data)[which(seurat_object@meta.data[, meta_data_column] == x)]})

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
    stop("Length of `colors_use` (", length(x = colors_use), ") does not equal number of plots (", num_plots, ").")
  }

  # Plots
  plots <- lapply(1:length(x = meta_sample_list), function(x) {
    plot <- DimPlot(seurat_object, cells = cell_names[[x]], group.by = meta_data_column, cols = colors_use[[x]], reduction = reduction, pt.size = pt.size, raster = raster, ...) +
      ggtitle(paste(meta_sample_list[[x]])) +
      theme(plot.title = element_text(hjust = 0.5, size = title_size),
            legend.position = "none") +
      xlim(x_axis) +
      ylim(y_axis)
  })

  # Wrap Plots into single output
  wrap_plots(plots, ncol = num_columns)
}


#' Custom Labeled Variable Features Plot
#'
#' Creates variable features plot with N number of features already labeled by default.
#'
#' @param seurat_object Seurat object name.
#' @param num_features Number of top variable features to highlight by color/label.
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
#' \dontrun{
#' VariableFeaturePlot_scCustom(seurat_object = object, num_features = 10)
#' }
#'

VariableFeaturePlot_scCustom <- function(
  seurat_object,
  num_features = 10,
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
  if (label) {
    plot <- LabelPoints(plot = plot, points = top_features, repel = repel)
  }

  # return log10 y axis
  if (y_axis_log) {
    plot <- plot + scale_y_log10()
    return(plot)
  }

  # Return plot
  return(plot)
}
