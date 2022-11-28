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
#' @param num_columns Number of columns in plot layout.
#' @param slot Which slot to pull expression data from?  Default is "data".
#' @param alpha_exp new alpha level to apply to expressing cell color palette (`colors_use`).  Must be
#' value between 0-1.
#' @param alpha_na_exp new alpha level to apply to non-expressing cell color palette (`na_color`).  Must be
#' value between 0-1.
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
#' @importFrom scales alpha
#' @importFrom Seurat FeaturePlot
#' @importFrom SeuratObject DefaultDimReduc
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' \dontrun{
#' FeaturePlot_scCustom(seurat_object = object, features = "Cx3cr1",
#' colors_use = viridis_plasma_dark_high, na_color = "lightgray")
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
  raster.dpi = c(512, 512),
  split.by = NULL,
  num_columns = NULL,
  slot = "data",
  alpha_exp = NULL,
  alpha_na_exp = NULL,
  label_feature_yaxis = FALSE,
  combine = TRUE,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check meta
  if (!is.null(x = split.by)) {
    split.by <- Meta_Present(seurat_object = seurat_object, meta_col_names = split.by, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  # Check features and meta to determine which features present to plot to avoid error in split plots
  if (!is.null(x = split.by)) {
    features_list <- Gene_Present(data = seurat_object, gene_list = features, omit_warn = FALSE, print_msg = FALSE, case_check_msg = FALSE, return_none = TRUE)

    meta_list <- Meta_Present(seurat_object = seurat_object, meta_col_names = features_list[[2]], omit_warn = FALSE, print_msg = FALSE, abort = FALSE)

    all_not_found_features <- meta_list[[2]]

    all_found_features <- c(features_list[[1]], meta_list[[1]])

    # Stop if no features found
    if (length(x = all_found_features) < 1) {
      cli_abort(message = c("No features were found.",
                            "*" = "The following are not present in object:",
                            "i" = "{glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}")
      )
    }

    # Return message of features not found
    if (length(x = all_not_found_features) > 0) {
      op <- options(warn = 1)
      on.exit(options(op))
      cli_warn(message = c("The following features were omitted as they were not found:",
                           "i" = "{glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}")
      )
    }

    # set to features to match remainder code
    features <- all_found_features
  }

  # Get length of meta data feature
  if (is.null(x = split.by) && label_feature_yaxis) {
    cli_abort(message = "Setting `label_feature_yaxis = TRUE` is only supported when also setting `split.by`.")
  }

  if (!is.null(x = split.by)) {
    split.by_length <- length(unique(seurat_object@meta.data[[split.by]]))

    if (!is.null(x = num_columns) && label_feature_yaxis) {

      cli_warn(message = c("Setting number of columns is not permitted if `label_feature_yaxis = TRUE`",
                           "i" = "Number of columns be automatically set to number of levels in `split.by` ({split.by_length}).")
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
                        "*" = "{split.by} only contains {split.by_length} variables.",
                        "i" = "Please adjust `num_columns` to be less than or equal to {split.by_length}.")
      )
    }
  }

  if (any(features %in% colnames(x = seurat_object@meta.data))) {
    cli_warn(message = c("Some of the plotted features are from meta.data slot.",
                         "*" = "Please check that `na_cutoff` param is being set appropriately for those features.")
    )
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

  # Add alpha to color scales
  if (!is.null(x = alpha_exp)) {
    colors_use <- alpha(colors_use, alpha_exp)
  }

  if (!is.null(x = alpha_na_exp)) {
    na_color <- alpha(na_color, alpha_exp)
  }

  # plot no split & combined
  if (is.null(x = split.by) && combine) {
    plot <- suppressMessages(FeaturePlot(object = seurat_object, features = features, order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, ncol = num_columns, combine = combine, raster.dpi = raster.dpi, ...) & scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, NA), na.value = na_color))
  }

  # plot no split & combined
  if (is.null(x = split.by) && !combine) {
    plot_list <- suppressMessages(FeaturePlot(object = seurat_object, features = features, order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, ncol = num_columns, combine = combine, raster.dpi = raster.dpi, ...))

    plot <- lapply(1:length(x = plot_list), function(i) {
      p[[i]] <- suppressMessages(p[[i]] + scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, NA), na.value = na_color))
    })
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

    plot <- suppressMessages(FeaturePlot(object = seurat_object, features = features, order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, raster.dpi = raster.dpi, ...) & scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, max_exp_value), na.value = na_color, name = features)) & RestoreLegend() & theme(axis.title.y.right = element_blank())

    if (label_feature_yaxis) {
      plot <- plot + plot_layout(nrow = num_rows, ncol = num_columns)
      plot <- plot & theme(legend.title=element_blank())
      plot <- suppressMessages(plot + scale_y_continuous(sec.axis = dup_axis(name = features))) + No_Right()
    } else {
      plot <- plot + plot_layout(nrow = num_rows, ncol = num_columns)
    }
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

      single_plot <- suppressMessages(FeaturePlot(object = seurat_object, features = features[i], order = order, pt.size = pt.size, reduction = reduction, raster = raster, split.by = split.by, raster.dpi = raster.dpi, ...) & scale_color_gradientn(colors = colors_use, limits = c(na_cutoff, max_exp_value), na.value = na_color, name = features[i])) & RestoreLegend() & theme(axis.title.y.right = element_blank())

      if (label_feature_yaxis) {
        single_plot <- single_plot + plot_layout(nrow = num_rows, ncol = num_columns)
        single_plot <- single_plot & theme(legend.title=element_blank())
        single_plot <- suppressMessages(single_plot + scale_y_continuous(sec.axis = dup_axis(name = features[i]))) + No_Right()
      } else {
        single_plot <- single_plot + plot_layout(nrow = num_rows, ncol = num_columns)
      }
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
#' @param na_color color to use for points below lower limit.
#' @param order whether to move positive cells to the top (default = TRUE).
#' @param pt.size Adjust point size for plotting.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default).
#' @param na_cutoff Value to use as minimum expression cutoff.  To set no cutoff set to `NA`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param slot Which slot to pull expression data from?  Default is "data".
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
  reduction = NULL,
  na_cutoff = 0.000000001,
  raster = NULL,
  raster.dpi = c(512, 512),
  slot = "data",
  num_columns = NULL,
  alpha_exp = NULL,
  alpha_na_exp = NULL,
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
    cli_abort(message = c("Assay 1: '{assay1}' has not been normalized.",
                          "i" = "Please run `NormalizeData` on this assay before proceeding to visualization."))
  }

  # Cell Bender normalize check
  if (!paste0("NormalizeData.", assay2) %in% commands) {
    cli_abort(message = c("Assay 2: '{assay2}' has not been normalized.",
                      "i" = "Please run `NormalizeData` on this assay before proceeding to visualization."))
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
    cli_warn("When plotting more than one feature `num_columns` refers to patchwork columns and must either be 1 (vertical) or 2 (horizontal).")
  }

  # Change assay and plot raw
  DefaultAssay(seurat_object) <- assay1

  plot_raw <- FeaturePlot_scCustom(seurat_object = seurat_object, features = features, slot = slot, colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, order = order, pt.size = pt.size, reduction = reduction, raster = raster, alpha_exp = alpha_exp, alpha_na_exp = alpha_na_exp, raster.dpi = raster.dpi, ...) & labs(color = assay1)

  # Change to cell bender and plot
  DefaultAssay(seurat_object) <- assay2

  plot_cell_bender <- FeaturePlot_scCustom(seurat_object = seurat_object, features = features, slot = slot, colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, order = order, pt.size = pt.size, reduction = reduction, raster = raster, alpha_exp = alpha_exp, alpha_na_exp = alpha_na_exp, raster.dpi = raster.dpi, ...) & labs(color = assay2)

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
  raster.dpi = c(512, 512),
  ggplot_default_colors = FALSE,
  color_seed = 123,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # split.by present
  if (is.null(x = split.by)) {
    cli_abort(message = "No value supplied to `split.by`.")
  }

  # Check split.by is valid
  if (split.by %in% colnames(seurat_object@meta.data) == FALSE) {
    cli_abort(message = c("The meta data variable: '{split.by}' could not be found in object@meta.data.",
                          "i" = "Please check the spelling and column names of meta.data slot.")
    )
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
    cli_abort(message = c("The number of columns specified is greater than the number of meta data variables.",
                          "*" = "'{split.by}' only contains {split.by_length} variables.",
                          "i" = "Please adjust `num_columns` to be less than or equal to {split.by_length}.")
    )
  }

  # Check features are present
  possible_features <- c(rownames(seurat_object), colnames(seurat_object@meta.data))
  check_features <- setdiff(x = c(feature1, feature2), y = possible_features)
  if (length(x = check_features) > 0) {
    cli_abort(message = "The following feature(s) were not present in Seurat object: '{check_features}'")
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
  if (class(x = seurat_object@meta.data[, split.by]) == "factor") {
    meta_sample_list <- as.character(x = levels(seurat_object@meta.data[, split.by]))
  } else {
    meta_sample_list <- as.character(unique(seurat_object@meta.data[, split.by]))
  }

  # Extract cell names per meta data list of values
  cell_names <- lapply(meta_sample_list, function(x) {
    row.names(seurat_object@meta.data)[which(seurat_object@meta.data[, split.by] == x)]})

  # raster check
  raster <- raster %||% (length(x = colnames(x = seurat_object)) > 2e5)

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
#' @param cluster_name Name(s) (or number(s)) identity of cluster to be highlighted.
#' @param highlight_color Color(s) to highlight cells (default "navy").
#' @param background_color non-highlighted cell colors.
#' @param pt.size point size for both highlighted cluster and background.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param label Whether to label the highlighted cluster(s).  Default is FALSE.
#' @param split.by Feature to split plots by (i.e. "orig.ident").
#' @param split_seurat logical.  Whether or not to display split plots like Seurat (shared y axis) or as
#' individual plots in layout.  Default is FALSE.
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
  raster.dpi = c(512, 512),
  label = FALSE,
  split.by = NULL,
  split_seurat = FALSE,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = colnames(x = seurat_object)) > 2e5)

  # Perform Idents check and report errors when when length(cluster_name) > 1
  if (length(x = cluster_name) > 1) {
    idents_list <- levels(x = Idents(object = seurat_object))

    good_idents <- cluster_name[cluster_name %in% idents_list]
    bad_idents <- cluster_name[!cluster_name %in% idents_list]

    if (length(x = bad_idents) > 0) {
      cli_warn("The following 'cluster_name(s)' were not found the active.ident slot: {bad_idents}")
    }
  }

  # pull cells to highlight in plot
  cells_to_highlight <- CellsByIdentities(seurat_object, idents = cluster_name)

  # set point size
  if (is.null(x = pt.size)) {
    pt.size <- AutoPointSize_scCustom(data = sum(lengths(cells_to_highlight)), raster = raster)
  }

  # Adjust colors if needed when length(cluster_name) > 1
  if (length(x = highlight_color) == 1 && length(x = cluster_name) > 1) {
    highlight_color <- rep(x = highlight_color, length(x = cluster_name))
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
  plot <- suppressMessages(plot & scale_color_manual(breaks = names(cells_to_highlight), values = c(highlight_color, background_color), na.value = background_color))

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
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param label Whether to label the highlighted meta data variable(s).  Default is FALSE.
#' @param split.by Variable in `@meta.data` to split the plot by.
#' @param split_seurat logical.  Whether or not to display split plots like Seurat (shared y axis) or as
#' individual plots in layout.  Default is FALSE.
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
#' \dontrun{
#' Meta_Highlight_Plot(seurat_object = object, meta_data_column = "orig.ident",
#' meta_data_highlight = "sample_01", highlight_color = "gold", background_color = "lightgray",
#' pt.size = 2)
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
  raster.dpi = c(512, 512),
  label = FALSE,
  split.by = NULL,
  split_seurat = FALSE,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check meta data
  good_meta_data_column <- Meta_Present(seurat_object = seurat_object, meta_col_names = meta_data_column, omit_warn = FALSE, print_msg = FALSE, abort = FALSE)[[1]]

  # stop if none found
  if (length(x = good_meta_data_column) == 0) {
    cli_abort(message = c("No 'meta_data_column' was not found.",
              "i" = "No column found in object meta.data named: {meta_data_column}.")
    )
  }

  # Check that meta data is factor or character
  accepted_meta_types <- c("factor", "character", "logical")

  if (!class(x = seurat_object@meta.data[[good_meta_data_column]]) %in% accepted_meta_types) {
    stop("The 'good_meta_data_column': ", good_meta_data_column, " is of class: ", '"', class(x = seurat_object@meta.data[[good_meta_data_column]]), '"', " only meta data variables of classes: factor, character, or logical can be used with Meta_Highlight_Plot().")
  }

  # Check meta_data_highlight
  meta_var_list <- as.character(unique(seurat_object@meta.data[, good_meta_data_column]))

  # Check good and bad highlight values
  bad_meta_highlight <- meta_var_list[!meta_var_list %in% meta_data_highlight]
  found_meta_highlight <- meta_var_list[meta_var_list %in% meta_data_highlight]

  # Abort if no meta_data_highlight found
  if (length(x = found_meta_highlight) == 0) {
    cli_abort(message = c("No 'meta_data_highlight' value(s) were not found.",
                          "i" = "The following 'meta_data_highlight' variables were not found in {good_meta_data_column}: {bad_meta_highlight}")
    )
  }

  # warn if some meta_data_highlight not found
  if (length(x = found_meta_highlight) != length(x = meta_data_highlight)) {
    cli_warn(message = c("Some 'meta_data_highlight' value(s) were not found.",
                          "i" = "The following 'meta_data_highlight' variables were not found in {good_meta_data_column} and were omitted: {bad_meta_highlight}")
    )
  }

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = colnames(x = seurat_object)) > 2e5)

  # Change default ident and pull cells to highlight in plot
  Idents(seurat_object) <- good_meta_data_column

  cells_to_highlight <- CellsByIdentities(seurat_object, idents = found_meta_highlight)

  # set point size
  if (is.null(x = pt.size)) {
    pt.size <- AutoPointSize_scCustom(data = sum(lengths(cells_to_highlight)), raster = raster)
  }

  # Adjust colors if needed when length(meta_data_highlight) > 1
  if (length(x = highlight_color) == 1 && length(x = found_meta_highlight) > 1) {
    highlight_color <- rep(x = highlight_color, length(x = found_meta_highlight))
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
  plot <- suppressMessages(plot & scale_color_manual(breaks = names(cells_to_highlight), values = c(highlight_color, background_color), na.value = background_color))

  return(plot)
}


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
#' @importFrom Seurat VlnPlot
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
#' VlnPlot_scCustom(seurat_object = object, features = "Cx3cr1")
#' }
#'

VlnPlot_scCustom <- function(
  seurat_object,
  features,
  colors_use = NULL,
  pt.size = NULL,
  group.by = NULL,
  split.by = NULL,
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

  # set size if NULL
  pt.size <- pt.size %||% AutoPointSize_scCustom(data = seurat_object)

  # Add raster check for scCustomize
  num_cells <- unlist(CellsByIdentities(object = seurat_object, idents = idents))

  if (is.null(x = raster)) {
    if (pt.size == 0) {
      raster <- FALSE
    } else {
      if (length(x = num_cells) * length(x = features) > 100000 && pt.size != 0) {
        raster <- TRUE
        cli_inform(message = c("NOTE: Rasterizing points since total number of points across all plots exceeds 100,000.",
                               "i" = "To plot in vector form set `raster=FALSE`")
        )
      } else {
        raster <- FALSE
      }
    }
  }

  # Set default color palette based on number of levels being plotted
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }

  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use) && ggplot_default_colors) {
    cli_abort(message = "Cannot provide both custom palette to `colors_use` and specify `ggplot_default_colors = TRUE`.")
  }

  # set default plot colors
  if (is.null(x = colors_use)) {
    colors_use <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
  }

  # Plot
  plot <- VlnPlot(object = seurat_object, features = features, cols = colors_use, pt.size = pt.size, idents = idents, group.by = group.by, split.by = split.by, ncol = num_columns, raster = raster, add.noise = add.noise, ...)

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
#' @param x_lab_rotate Rotate x-axis labels 45 degrees (Default is FALSE).
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
#' @param pt.size Adjust point size for plotting.  Default for `StackedVlnPlot` is 0 to avoid issues with
#' rendering so many points in vector form.  Alteratively, see `raster` parameter.
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
  plot_legend = FALSE,
  colors_use = NULL,
  color_seed = 123,
  ggplot_default_colors = FALSE,
  plot_spacing = 0.15,
  spacing_unit = "cm",
  pt.size = NULL,
  raster = NULL,
  add.noise = TRUE,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check features and meta to determine which features present
  features_list <- Gene_Present(data = seurat_object, gene_list = features, omit_warn = FALSE, print_msg = FALSE, case_check_msg = FALSE, return_none = TRUE)

  meta_list <- Meta_Present(seurat_object = seurat_object, meta_col_names = features_list[[2]], omit_warn = FALSE, print_msg = FALSE, abort = FALSE)

  all_not_found_features <- meta_list[[2]]

  all_found_features <- c(features_list[[1]], meta_list[[1]])

  # Stop if no features found
  if (length(x = all_found_features) < 1) {
    cli_abort(message = c("No features were found.",
                          "*" = "The following are not present in object:",
                          "i" = "{glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}")
    )
  }

  # Return message of features not found
  if (length(x = all_not_found_features) > 0) {
    op <- options(warn = 1)
    on.exit(options(op))
    cli_warn(message = c("The following features were omitted as they were not found:",
                         "i" = "{glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}")
    )
  }

  # Check feature case correct
  Case_Check(seurat_object = seurat_object, gene_list = all_not_found_features, case_check_msg = TRUE, return_features = FALSE)

  # set pt.size (default is no points)
  if (is.null(x = pt.size)) {
    pt.size <- 0
  }

  # Set rasterization
  num_cells <- unlist(CellsByIdentities(object = seurat_object, idents = idents))

  if (length(x = num_cells) * length(x = all_found_features) > 100000 && is.null(x = raster) && pt.size != 0) {
    raster <- TRUE
    cli_inform(message = c("NOTE: Rasterizing points since total number of points across all plots exceeds 100,000.",
                           "i" = "To plot in vector form set `raster=FALSE`")
    )
  }

  # Set default color palette based on number of levels being plotted
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }

  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use) && ggplot_default_colors) {
    cli_abort(message = "Cannot provide both custom palette to `colors_use` and specify `ggplot_default_colors = TRUE`.")
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

  # Create final plot patchwork to return
  plot_return <- wrap_plots(plotlist = plot_list, ncol = 1)

  # add single legend
  if (plot_legend) {
    plot_return <- plot_return & theme(legend.position = "right")
    plot_return <- plot_return + plot_layout(guides = 'collect')
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

  # Check features and meta to determine which features present
  features_list <- Gene_Present(data = seurat_object, gene_list = features, omit_warn = FALSE, print_msg = FALSE, case_check_msg = FALSE, return_none = TRUE)

  meta_list <- Meta_Present(seurat_object = seurat_object, meta_col_names = features_list[[2]], omit_warn = FALSE, print_msg = FALSE, abort = FALSE)

  all_not_found_features <- meta_list[[2]]

  all_found_features <- c(features_list[[1]], meta_list[[1]])

  # Stop if no features found
  if (length(x = all_found_features) < 1) {
    cli_abort(message = c("No features were found.",
                               "*" = "The following are not present in object:",
                               "i" = "{glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}")
    )
  }

  # Return message of features not found
  if (length(x = all_not_found_features) > 0) {
    op <- options(warn = 1)
    on.exit(options(op))
    cli_warn(message = c("The following features were omitted as they were not found:",
                              "i" = "{glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}")
    )
  }

  # Check feature case correct
  Case_Check(seurat_object = seurat_object, gene_list = all_not_found_features, case_check_msg = TRUE, return_features = FALSE)

  # Plot
  plot <- suppressMessages(DotPlot(object = seurat_object, features = all_found_features, ...) +
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


#' Clustered DotPlot
#'
#' Clustered DotPlots using ComplexHeatmap
#'
#' @param seurat_object Seurat object name.
#' @param features Features to plot.
#' @param colors_use_exp Color palette to use for plotting expression scale.  Default is `viridis::plasma(n = 20, direction = -1)`.
#' @param exp_color_min Minimum scaled average expression threshold (everything smaller will be set to this).
#' Default is -2.
#' @param exp_color_middle What scaled expression value to use for the middle of the provided `colors_use_exp`.
#' By default will be set to value in middle of `exp_color_min` and `exp_color_max`.
#' @param exp_color_max Minimum scaled average expression threshold (everything smaller will be set to this).
#' Default is 2.
#' @param print_exp_quantiles Whether to print the quantiles of expression data in addition to plots.
#' Default is FALSE.  NOTE: These values will be altered by choices of `exp_color_min` and `exp_color_min`
#' if there are values below or above those cutoffs, respectively.
#' @param colors_use_idents specify color palette to used for identity labels.  By default if
#' number of levels plotted is less than or equal to 36 it will use "polychrome" and if greater than 36
#' will use "varibow" with shuffle = TRUE both from `DiscretePalette_scCustomize`.
#' @param x_lab_rotate How to rotate column labels.  By default set to `TRUE` which rotates labels 45 degrees.
#' If set `FALSE` rotation is set to 0 degrees.  Users can also supply custom angle for text rotation.
#' @param k Value to use for k-means clustering on rows.  Sets (km) parameter in `ComplexHeatmap::Heatmap()`.
#' From `ComplexHeatmap::Heatmap()`: Apply k-means clustering on rows. If the value is larger than 1, the
#' heatmap will be split by rows according to the k-means clustering. For each row slice, hierarchical
#' clustering is still applied with parameters above.
#' @param row_km_repeats Number of k-means runs to get a consensus k-means clustering. Note if row_km_repeats
#' is set to more than one, the final number of groups might be smaller than row_km, but this might
#' mean the original row_km is not a good choice.  Default is 1000.
#' @param column_km_repeats Number of k-means runs to get a consensus k-means clustering. Similar as row_km_repeats.
#' Default is 100.
#' @param row_label_size Size of the feature labels.  Provided to `row_names_gp` in Heatmap call.
#' @param raster Logical, whether to render in raster format (faster plotting, smaller files).  Default is FALSE.
#' @param plot_km_elbow Logical, whether or not to return the Sum Squared Error Elbow Plot for k-means clustering.
#' Estimating elbow of this plot is one way to determine "optimal" value for `k`.
#' Based on: https://stackoverflow.com/a/15376462/15568251.
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
#' @importFrom dplyr filter select
#' @importFrom grid grid.circle grid.rect gpar
#' @importFrom magrittr "%>%"
#' @importFrom Seurat DotPlot PackageCheck
#' @importFrom tidyr pivot_wider
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @author Ming Tang (Original Code), Sam Marsh (Wrap single function, added/modified functionality)
#' @references https://divingintogeneticsandgenomics.rbind.io/post/clustered-dotplot-for-single-cell-rnaseq/
#' @seealso https://twitter.com/tangming2005
#'
#' @examples
#' \dontrun{
#' Clustered_DotPlot(seurat_object = object, features = gene_list)
#' }
#'

Clustered_DotPlot <- function(
  seurat_object,
  features,
  colors_use_exp = viridis_plasma_dark_high,
  exp_color_min = -2,
  exp_color_middle = NULL,
  exp_color_max = 2,
  print_exp_quantiles = FALSE,
  colors_use_idents = NULL,
  x_lab_rotate = TRUE,
  k = 1,
  row_km_repeats = 1000,
  column_km_repeats = 1000,
  row_label_size = 8,
  raster = FALSE,
  plot_km_elbow = TRUE,
  elbow_kmax = NULL,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  show_parent_dend_line = TRUE,
  ggplot_default_colors = FALSE,
  color_seed = 123,
  seed = 123
) {
  # Check for packages
  ComplexHeatmap_check <- PackageCheck("ComplexHeatmap", error = FALSE)
  if (!ComplexHeatmap_check[1]) {
    stop(
      "Please install the ComplexHeatmap package to use Clustered_DotPlot",
      "\nThis can be accomplished with the following commands: ",
      "\n----------------------------------------",
      "\ninstall.packages('BiocManager')",
      "\nBiocManager::install('ComplexHeatmap')",
      "\n----------------------------------------",
      call. = FALSE
    )
  }

  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check unique features
  features_unique <- unique(x = features)

  if (length(x = features_unique) != length(x = features)) {
    cli_warn("Feature list contains duplicates, making unique.")
  }

  # Check exp min/max set correctly
  if (!exp_color_min < exp_color_max) {
    cli_abort(message = c("Expression color min/max values are not compatible.",
                          "i" = "The value for 'exp_color_min': {exp_color_min} must be less than the value for 'exp_color_max': {exp_color_max}.")
    )
  }

  # Get DotPlot data
  seurat_plot <- DotPlot(object = seurat_object, features = features_unique, assay = assay, group.by = group.by, scale = TRUE, idents = idents, col.min = NULL, col.max = NULL)

  data <- seurat_plot$data

  # Get expression data
  exp_mat <- data %>%
    select(-pct.exp, -avg.exp) %>%
    pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
    as.data.frame()

  row.names(x = exp_mat) <- exp_mat$features.plot

  # Check NAs if idents
  if (!is.null(x = idents)) {
    # Find NA features and print warning
    excluded_features <- exp_mat[rowSums(is.na(x = exp_mat)) > 0,] %>%
      rownames()
    cli_warn(message = c("Some scaled data missing.",
                         "*" = "The following features were removed as there is no scaled expression present in subset (`idents`) of object provided:",
                         "i" = "{glue_collapse_scCustom(input_string = excluded_features, and = TRUE)}.")
    )

    # Extract good features
    good_features <- rownames(exp_mat)

    # Remove rows with NAs
    exp_mat <- exp_mat %>%
      filter(features.plot %in% good_features)
  }

  exp_mat <- exp_mat[,-1] %>%
    as.matrix()

  # Get percent expressed data
  percent_mat <- data %>%
    select(-avg.exp, -avg.exp.scaled) %>%
    pivot_wider(names_from = id, values_from = pct.exp) %>%
    as.data.frame()

  row.names(x = percent_mat) <- percent_mat$features.plot

  # Subset dataframe for NAs if idents so that exp_mat and percent_mat match
  if (!is.null(x = idents)) {
    percent_mat <- percent_mat %>%
      filter(features.plot %in% good_features)
  }

  percent_mat <- percent_mat[,-1] %>%
    as.matrix()

  # print quantiles
  if (print_exp_quantiles) {
    cli_inform(message = "Quantiles of gene expression data are:")
    print(quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99)))
  }

  # set assay (if null set to active assay)
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Set default color palette based on number of levels being plotted
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }

  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use_idents) && ggplot_default_colors) {
    cli_abort(message = "Cannot provide both custom palette to `colors_use` and specify `ggplot_default_colors = TRUE`.")
  }
  if (is.null(x = colors_use_idents)) {
    # set default plot colors
    colors_use_idents <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
  }

  # Reduce color length list due to naming requirement
  colors_use_idents <- colors_use_idents[1:group_by_length]

  # Modify if class = "colors"
  if (class(x = colors_use_idents) == "colors") {
    colors_use_idents <- as.vector(colors_use_idents)
  }

  # Pull Annotation and change colors to ComplexHeatmap compatible format
  Identity <- colnames(exp_mat)

  identity_colors <- colors_use_idents
  names(identity_colors) <- Identity
  identity_colors_list <- list(Identity = identity_colors)

  # Create identity annotation
  column_ha <- ComplexHeatmap::HeatmapAnnotation(Identity = Identity,
                                                 col =  identity_colors_list,
                                                 na_col = "grey",
                                                 name = "Identity",
                                                 show_legend = FALSE
  )

  # Set middle of color scale if not specified
  if (is.null(x = exp_color_middle)) {
    exp_color_middle <- Middle_Number(min = exp_color_min, max = exp_color_max)
  }

  palette_length <- length(colors_use_exp)
  palette_middle <- Middle_Number(min = 0, max = palette_length)

  # Create palette
  col_fun = colorRamp2(c(exp_color_min, exp_color_middle, exp_color_max), colors_use_exp[c(1,palette_middle, palette_length)])

  # Calculate and plot Elbow
  if (plot_km_elbow) {
    # if elbow_kmax not NULL check it is usable
    if (!is.null(x = elbow_kmax) && elbow_kmax > (nrow(x = exp_mat) - 1)) {
      elbow_kmax <- nrow(x = exp_mat) - 1
      cli_warn(message = c("The value provided for 'elbow_kmax' is too large.",
                           "i" = "Changing to (length(x = features)-1): {elbow_kmax}")
      )
    }

    # if elbow_kmax is NULL set value based on input feature list
    if (is.null(x = elbow_kmax)) {
      # set to (length(x = features)-1) if less than 21 features OR to 20 if greater than 21 features
      if (nrow(x = exp_mat) > 21) {
        elbow_kmax <- 20
      } else {
        elbow_kmax <- nrow(x = exp_mat) - 1
      }
    }

    km_elbow_plot <- kMeans_Elbow(data = exp_mat, k_max = elbow_kmax)
  }

  # prep heatmap
  if (raster) {
    layer_fun = function(j, i, x, y, w, h, fill) {
      grid.rect(x = x, y = y, width = w, height = h,
                      gp = gpar(col = NA, fill = NA))
      grid.circle(x=x,y=y,r= sqrt(ComplexHeatmap::pindex(percent_mat, i, j)/100)  * unit(2, "mm"),
                        gp = gpar(fill = col_fun(ComplexHeatmap::pindex(exp_mat, i, j)), col = NA))
    }
  } else {
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.rect(x = x, y = y, width = w, height = h,
                      gp = gpar(col = NA, fill = NA))
      grid.circle(x=x,y=y,r= sqrt(percent_mat[i, j]/100) * unit(2, "mm"),
                        gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))
    }
  }

  # Create legend for point size
  lgd_list = list(
    ComplexHeatmap::Legend(at = Identity, title = "Identity", legend_gp = gpar(fill = identity_colors_list[[1]])),
    ComplexHeatmap::Legend(labels = c(0.25,0.5,0.75,1), title = "Percent Expressing",
                           graphics = list(
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(2, "mm"),
                                                                    gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(2, "mm"),
                                                                    gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(2, "mm"),
                                                                    gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                                                    gp = gpar(fill = "black")))
    )
  )

  # Set x label roration
  if (is.numeric(x = x_lab_rotate)) {
    x_lab_rotate <- x_lab_rotate
  } else if (isTRUE(x = x_lab_rotate)) {
    x_lab_rotate <- 45
  } else {
    x_lab_rotate <- 0
  }

  # Create Plot
  set.seed(seed = seed)
  if (raster) {
    cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                                heatmap_legend_param=list(title="Expression"),
                                                col=col_fun,
                                                rect_gp = gpar(type = "none"),
                                                layer_fun = layer_fun,
                                                row_names_gp = gpar(fontsize = row_label_size),
                                                row_km = k,
                                                row_km_repeats = row_km_repeats,
                                                border = "black",
                                                top_annotation = column_ha,
                                                column_km_repeats = column_km_repeats,
                                                show_parent_dend_line = show_parent_dend_line,
                                                column_names_rot = x_lab_rotate)
  } else {
    cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                                heatmap_legend_param=list(title="Expression"),
                                                col=col_fun,
                                                rect_gp = gpar(type = "none"),
                                                cell_fun = cell_fun,
                                                row_names_gp = gpar(fontsize = row_label_size),
                                                row_km = k,
                                                row_km_repeats = row_km_repeats,
                                                border = "black",
                                                top_annotation = column_ha,
                                                column_km_repeats = column_km_repeats,
                                                show_parent_dend_line = show_parent_dend_line,
                                                column_names_rot = x_lab_rotate)
  }

  # Add pt.size legend & return plots
  if (plot_km_elbow) {
    return(list(km_elbow_plot, ComplexHeatmap::draw(cluster_dot_plot, annotation_legend_list = lgd_list)))
  }
  return(ComplexHeatmap::draw(cluster_dot_plot, annotation_legend_list = lgd_list))
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
#' this is simply a wrapper to alter some of the default parameters (https://github.com/satijalab/seurat/blob/master/R/visualization.R) (Licence: GPL-3).
#' `figure_plot` parameter/code modified from code by Tim Stuart via twitter: (https://twitter.com/timoast/status/1526237116035891200?s=20&t=foJOF81aPSjr1t7pk1cUPg).
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
  figure_plot = FALSE,
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

  # if split.by is null set split_seurat to TRUE
  if (is.null(x = split.by)) {
    split_seurat <- TRUE
  }

  # figure_plot check
  if (figure_plot && !split_seurat) {
    cli_abort(message = "'figure_plot' can only be TRUE if split_seurat is FALSE.")
  }

  # Set default color palette based on number of levels being plotted
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }

  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use) && ggplot_default_colors) {
    cli_abort(message = "Cannot provide both custom palette to `colors_use` and specify `ggplot_default_colors = TRUE`.")
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
    if (figure_plot) {

      # pull axis labels
      x_lab_reduc <- plot$labels$x
      y_lab_reduc <- plot$labels$y

      plot <- plot & NoAxes()

      axis_plot <- ggplot(data.frame(x= 100, y = 100), aes(x = x, y = y)) +
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
      return(plot_figure)
    } else {
      return(plot)
    }

  } else {
    if (split_seurat) {
      # Plot Seurat Splitting
      plot <- DimPlot(object = seurat_object, cols = colors_use, pt.size = pt.size, reduction = reduction, group.by = group.by, split.by = split.by, shuffle = shuffle, seed = seed, label = label, label.size = label.size, label.color = label.color, repel = repel, raster = raster, raster.dpi = raster.dpi, ncol = num_columns, dims = dims, label.box = label.box, ...)
      if (figure_plot) {

        # pull axis labels
        x_lab_reduc <- plot$labels$x
        y_lab_reduc <- plot$labels$y

        plot <- plot & NoAxes()

        axis_plot <- ggplot(data.frame(x= 100, y = 100), aes(x = x, y = y)) +
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
        return(plot_figure)
      } else {
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
      all_cells <- colnames(x = seurat_object)
      reduc_coordinates <- Embeddings(object = seurat_object[[reduction]])[all_cells, dims]
      reduc_coordinates <- as.data.frame(x = reduc_coordinates)
      x_axis <- c(min(reduc_coordinates[, 1]),
                  max(reduc_coordinates[, 1]))
      y_axis <- c(min(reduc_coordinates[, 2]),
                  max(reduc_coordinates[, 2]))

      # Extract cell names per meta data list of values
      # Extract split.by list of values
      if (class(x = seurat_object@meta.data[, split.by]) == "factor") {
        split_by_list <- as.character(x = levels(seurat_object@meta.data[, split.by]))
      } else {
        split_by_list <- as.character(unique(seurat_object@meta.data[, split.by]))
      }

      cell_names <- lapply(split_by_list, function(x) {
        row.names(seurat_object@meta.data)[which(seurat_object@meta.data[, split.by] == x)]})

      # Unify colors across plots
      if (is.null(x = group.by)) {
        levels_overall <- levels(x = Idents(object = seurat_object))
      } else {
        levels_overall <- levels(x = seurat_object@meta.data[[group.by]])
      }

      colors_overall <- colors_use

      names(colors_overall) <- levels_overall

      # plot
      plots <- lapply(1:length(x = split_by_list), function(x) {
        plot <- DimPlot(object = seurat_object, cells = cell_names[[x]], group.by = group.by, cols = colors_use, reduction = reduction, pt.size = pt.size, raster = raster, raster.dpi = raster.dpi, shuffle = shuffle, seed = seed, label = label, label.size = label.size, label.color = label.color, repel = repel, dims = dims, label.box = label.box, ...) +
          ggtitle(paste(split_by_list[[x]])) +
          theme(plot.title = element_text(hjust = 0.5),
                legend.position = "right") +
          xlim(x_axis) +
          ylim(y_axis)

        # Normalize the colors across all plots
        plot <- suppressMessages(plot + scale_color_manual(values = colors_overall))

        if (!is.null(x = group.by)) {
          plot <- plot + labs(color=group.by)
        } else {
          plot <- plot
        }
      })

      # Wrap Plots into single output
      plots <- wrap_plots(plots, ncol = num_columns) + plot_layout(guides = 'collect')
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
  raster.dpi = c(512, 512),
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check meta data column exists
  if (meta_data_column %in% colnames(seurat_object@meta.data) == FALSE) {
    cli_abort(message = c("Meta data variable not found.",
                          "i" = "The `meta_data_column`: {meta_data_column}, could not be found in object@meta.data.")
    )
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
  if (class(x = seurat_object@meta.data[, meta_data_column]) == "factor") {
    meta_sample_list <- levels(x = seurat_object@meta.data[, meta_data_column])
  } else {
    meta_sample_list <- as.character(unique(seurat_object@meta.data[, meta_data_column]))
  }

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
    cli_abort(message = c("Not enough colors provided.",
                          "i" = "Length of `colors_use` ({length(x = colors_use)}) does not equal number of plots ({num_plots}).")
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
