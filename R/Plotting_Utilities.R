#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### MODDED PLOTS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' PC Plots
#'
#' Plot PC Heatmaps and Dim Loadings for exploratory analysis.  Plots a single Heatmap and Gene Loading Plot.
#'   Used for PC_Loading_Plots function.
#'
#' @param seurat_object Seurat Object.
#' @param dim_number A single dim to plot (integer).
#'
#' @return A plot of PC heatmap and gene loadings for single
#'
#' @importFrom Seurat PCHeatmap VizDimLoadings
#' @import patchwork
#' @import ggplot2
#'
#' @seealso \code{\link[Seurat]{PCHeatmap}} and \code{\link[Seurat]{VizDimLoadings}}
#'
#' @export
#'
#' @concept seurat_plotting
#'
#' @examples
#' library(Seurat)
#' PC_Plotting(seurat_object = pbmc_small, dim_number = 1)
#'

PC_Plotting <- function(
  seurat_object,
  dim_number
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  p1 <- PCHeatmap(seurat_object,
                  cells = 500,
                  dims = dim_number,
                  balanced = TRUE,
                  fast = FALSE) +
    ggtitle(label = paste("PC", dim_number, sep="")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 22))
  p2 <- VizDimLoadings(seurat_object,
                       dims = dim_number,
                       reduction = "pca")
  p1 / p2
}


#' Modified Violin Plot
#'
#' Custom Violin Plot for stacked violin function.
#'
#' @param seurat_object Seurat object name.
#' @param features Features to plot.
#' @param plot_margin Modify white space in between plots.
#' @param ... takes arguments from \code{\link[Seurat]{VlnPlot}}.
#'
#' @return A modified violin plot
#'
#' @import ggplot2
#' @importFrom Seurat VlnPlot
#'
#' @noRd
#'
#' @author Ming Tang (Original Code), Sam Marsh (Modified function for use in scCustomtize)
#' @references \url{https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/}.  Solution for re-enabling plot spacing modification by Abdenour ABBAS (comment on original blog post; \url{http://disq.us/p/2b54qh2}).
#' @seealso \url{https://twitter.com/tangming2005}
#'

Modify_VlnPlot <- function(
  seurat_object,
  features,
  pt.size = NULL,
  cols = NULL,
  plot_margin = NULL,
  raster = NULL,
  add.noise = TRUE,
  ...
) {
  #remove the x-axis text and tick
  #plot_margin to adjust the white space between each plot.
  VlnPlot(seurat_object, features = features, pt.size = pt.size, cols = cols, raster = raster, add.noise = add.noise, ...)  +
    xlab("") +
    ylab(features) +
    ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0),
          axis.text.y = element_text(size = rel(1)),
          plot.margin = plot_margin,
          plot.title= element_blank(),
          axis.title.x = element_blank())
}


#' Sum Squared Error Elbow Plot
#'
#' Sum Squared Error Elbow Plot as method of estimating optimal k value for k-means clustering.
#'
#' @param data Expression data.
#' @param k_max Maximum number of k values to test.
#' @param plot_title Title of the plot.
#' @param cutoff_value Value to use for adding dashed line visualizing choice of k.
#'
#' @return A ggplot2 object.
#'
#' @import ggplot2
#' @importFrom magrittr "%>%"
#' @importFrom stats kmeans var
#' @importFrom tibble rownames_to_column
#'
#' @noRd
#'
#' @references Code to calculate wss values from: \url{https://stackoverflow.com/a/15376462/15568251}
#'

kMeans_Elbow <- function(
  data,
  k_max = 15,
  plot_title = "Sum of Squared Error (SSE) Plot",
  cutoff_value = NULL
) {
  # Calculate the within squares
  # code from @Ben https://stackoverflow.com/a/15376462/15568251
  wss <- (nrow(x = data)-1)*sum(apply(data,2,var))
  for (i in 2:k_max) wss[i] <- sum(kmeans(data,
                                          centers=i)$withinss)

  # Reformat for ggplot2 plotting
  plot_data <- data.frame(wss) %>%
    rownames_to_column("k")

  plot_data$k <- as.numeric(x = plot_data$k)

  # Plot data
  plot <- ggplot(data = plot_data, mapping = aes(y = wss, x = .data[["k"]])) +
    geom_point() +
    geom_path() +
    scale_x_continuous(n.breaks = k_max) +
    theme_ggprism_mod() +
    xlab("k (Number of Clusters)") +
    ylab("Within groups sum of squares") +
    ggtitle(plot_title) +
    geom_vline(xintercept = cutoff_value, linetype = "dashed", color = "red")

  return(plot)
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
#' @importFrom dplyr filter
#' @importFrom magrittr "%>%"
#' @importFrom Seurat FeatureScatter
#' @importFrom stats cor
#'
#' @noRd
#'

scCustomze_Split_FeatureScatter <- function(
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


#' Figure Plots
#'
#' Removes the axes from 2D DR plots and makes them into plot label.
#' Used for `figure_plot` parameter in plotting functions.
#'
#' @param plot 2D DR plot
#'
#' @return A modified plot
#'
#' @import ggplot2
#' @import patchwork
#'
#' @references parameter/code modified from code by Tim Stuart via twitter: \url{https://twitter.com/timoast/status/1526237116035891200?s=20&t=foJOF81aPSjr1t7pk1cUPg}.
#'
#' @noRd
#'

Figure_Plot <- function(
    plot
){
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

  return(plot_figure)
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
#' @noRd
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

Clustered_DotPlot_Single_Group <- function(
    seurat_object,
    features,
    colors_use_exp = viridis_plasma_dark_high,
    exp_color_min = -2,
    exp_color_middle = NULL,
    exp_color_max = 2,
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
    color_seed = 123,
    seed = 123
) {
  # Check for packages
  ComplexHeatmap_check <- is_installed(pkg = "ComplexHeatmap")
  if (isFALSE(x = ComplexHeatmap_check)) {
    cli_abort(message = c(
      "Please install the {.val ComplexHeatmap} package to use {.code Clustered_DotPlot}",
      "i" = "This can be accomplished with the following commands: ",
      "----------------------------------------",
      "{.field `install.packages({symbol$dquote_left}BiocManager{symbol$dquote_right})`}",
      "{.field `BiocManager::install({symbol$dquote_left}ComplexHeatmap{symbol$dquote_right})`}",
      "----------------------------------------"
    ))
  }

  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # set assay (if null set to active assay)
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Check acceptable fontface
  if (!row_label_fontface %in% c("plain", "bold", "italic", "oblique", "bold.italic")) {
    cli_abort(message = c("{.code row_label_face} {.val {row_label_face}} not recognized.",
                          "i" = "Must be one of {.val plain}, {.val bold}, {.val italic}, {.val olique}, or {.val bold.italic}."))
  }

  # Check unique features
  features_unique <- unique(x = features)

  if (length(x = features_unique) != length(x = features)) {
    cli_warn("Feature list contains duplicates, making unique.")
  }

  # Check features and meta to determine which features present
  all_found_features <- Feature_PreCheck(object = seurat_object, features = features_unique, assay = assay)

  # Check exp min/max set correctly
  if (!exp_color_min < exp_color_max) {
    cli_abort(message = c("Expression color min/max values are not compatible.",
                          "i" = "The value for {.code exp_color_min}: {.field {exp_color_min}} must be less than the value for {.code exp_color_max}: {.field {exp_color_max}}.")
    )
  }

  # Get DotPlot data
  seurat_plot <- DotPlot(object = seurat_object, features = all_found_features, assay = assay, group.by = group.by, scale = TRUE, idents = idents, col.min = NULL, col.max = NULL)

  data <- seurat_plot$data

  # Get expression data
  exp_mat <- data %>%
    select(-any_of(c("pct.exp", "avg.exp"))) %>%
    pivot_wider(names_from = any_of("id"), values_from = any_of("avg.exp.scaled")) %>%
    as.data.frame()

  row.names(x = exp_mat) <- exp_mat$features.plot

  # Check NAs if idents
  if (!is.null(x = idents)) {
    # Find NA features and print warning
    excluded_features <- exp_mat[rowSums(is.na(x = exp_mat)) > 0,] %>%
      rownames()
    cli_warn(message = c("Some scaled data missing.",
                         "*" = "The following features were removed as there is no scaled expression present in subset (`idents`) of object provided:",
                         "i" = "{.field {glue_collapse_scCustom(input_string = excluded_features, and = TRUE)}}.")
    )

    # Extract good features
    good_features <- rownames(x = exp_mat)

    # Remove rows with NAs
    exp_mat <- exp_mat %>%
      filter(.data[["features.plot"]] %in% good_features)
  }

  exp_mat <- exp_mat[,-1] %>%
    as.matrix()

  # Get percent expressed data
  percent_mat <- data %>%
    select(-any_of(c("avg.exp", "avg.exp.scaled"))) %>%
    pivot_wider(names_from = any_of("id"), values_from = any_of("pct.exp")) %>%
    as.data.frame()

  row.names(x = percent_mat) <- percent_mat$features.plot

  # Subset dataframe for NAs if idents so that exp_mat and percent_mat match
  if (!is.null(x = idents)) {
    percent_mat <- percent_mat %>%
      filter(.data[["features.plot"]] %in% good_features)
  }

  percent_mat <- percent_mat[,-1] %>%
    as.matrix()

  # print quantiles
  if (isTRUE(x = print_exp_quantiles)) {
    cli_inform(message = "Quantiles of gene expression data are:")
    print(quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99)))
  }

  # Set default color palette based on number of levels being plotted
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }

  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use_idents) && isTRUE(x = ggplot_default_colors)) {
    cli_abort(message = "Cannot provide both custom palette to {.code colors_use} and specify {.code ggplot_default_colors = TRUE}.")
  }
  if (is.null(x = colors_use_idents)) {
    # set default plot colors
    colors_use_idents <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
  }

  # Reduce color length list due to naming requirement
  colors_use_idents <- colors_use_idents[1:group_by_length]

  # Modify if class = "colors"
  if (inherits(x = colors_use_idents, what = "colors")) {
    colors_use_idents <- as.vector(x = colors_use_idents)
  }

  # Pull Annotation and change colors to ComplexHeatmap compatible format
  Identity <- colnames(x = exp_mat)

  identity_colors <- colors_use_idents
  names(x = identity_colors) <- Identity
  identity_colors_list <- list(Identity = identity_colors)

  # check grid color
  if (is.null(x = grid_color)) {
    grid_color <- NA
  } else {
    if (length(x = grid_color) > 1) {
      cli_abort(message = "{.code grid_color} can only be a single value.")
    }
    if (isTRUE(x = Is_Color(colors = grid_color))) {
      grid_color <- grid_color
    } else {
      cli_abort(message = "Value provided to {.code grid_color} ({.field {grid_color}}) is not valid value for color in R.")
    }
  }

  # Create identity annotation
  if (isTRUE(x = flip)) {
    column_ha <- ComplexHeatmap::rowAnnotation(Identity = Identity,
                                               col =  identity_colors_list,
                                               na_col = "grey",
                                               name = "Identity",
                                               show_legend = FALSE
    )
  } else {
    column_ha <- ComplexHeatmap::HeatmapAnnotation(Identity = Identity,
                                                   col =  identity_colors_list,
                                                   na_col = "grey",
                                                   name = "Identity",
                                                   show_legend = FALSE
    )
  }

  # Set middle of color scale if not specified
  if (is.null(x = exp_color_middle)) {
    exp_color_middle <- Middle_Number(min = exp_color_min, max = exp_color_max)
  }

  palette_length <- length(x = colors_use_exp)
  palette_middle <- Middle_Number(min = 0, max = palette_length)

  # Create palette
  col_fun = colorRamp2(c(exp_color_min, exp_color_middle, exp_color_max), colors_use_exp[c(1,palette_middle, palette_length)])

  # Calculate and plot Elbow
  if (isTRUE(x = plot_km_elbow)) {
    # if elbow_kmax not NULL check it is usable
    if (!is.null(x = elbow_kmax) && elbow_kmax > (nrow(x = exp_mat) - 1)) {
      elbow_kmax <- nrow(x = exp_mat) - 1
      cli_warn(message = c("The value provided for {.code elbow_kmax} is too large.",
                           "i" = "Changing to (length(x = features)-1): {.field {elbow_kmax}}.")
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
  if (isTRUE(x = flip)) {
    if (isTRUE(x = raster)) {
      layer_fun_flip = function(i, j, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(col = grid_color, fill = NA))
        grid.circle(x=x,y=y,r= sqrt(ComplexHeatmap::pindex(percent_mat, i, j)/100)  * unit(2, "mm"),
                    gp = gpar(fill = col_fun(ComplexHeatmap::pindex(exp_mat, i, j)), col = NA))
      }
    } else {
      cell_fun_flip = function(i, j, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(col = grid_color, fill = NA))
        grid.circle(x=x,y=y,r= sqrt(percent_mat[i, j]/100) * unit(2, "mm"),
                    gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))
      }
    }
  } else {
    if (isTRUE(x = raster)) {
      layer_fun = function(j, i, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(col = grid_color, fill = NA))
        grid.circle(x=x,y=y,r= sqrt(ComplexHeatmap::pindex(percent_mat, i, j)/100)  * unit(2, "mm"),
                    gp = gpar(fill = col_fun(ComplexHeatmap::pindex(exp_mat, i, j)), col = NA))
      }
    } else {
      cell_fun = function(j, i, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(col = grid_color, fill = NA))
        grid.circle(x=x,y=y,r= sqrt(percent_mat[i, j]/100) * unit(2, "mm"),
                    gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))
      }
    }
  }

  # Create legend for point size
  lgd_list = list(
    ComplexHeatmap::Legend(at = Identity, title = "Identity", legend_gp = gpar(fill = identity_colors_list[[1]]), labels_gp = gpar(fontsize = legend_label_size), title_gp = gpar(fontsize = legend_title_size, fontface = "bold")),
    ComplexHeatmap::Legend(labels = c(10,25,50,75,100), title = "Percent Expressing",
                           graphics = list(
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.1) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.50) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                                              gp = gpar(fill = "black"))),
                           labels_gp = gpar(fontsize = legend_label_size),
                           title_gp = gpar(fontsize = legend_title_size, fontface = "bold")
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
  if (isTRUE(x = raster)) {
    if (isTRUE(x = flip)) {
      cluster_dot_plot <- ComplexHeatmap::Heatmap(t(exp_mat),
                                                  heatmap_legend_param=list(title="Expression", labels_gp = gpar(fontsize = legend_label_size), title_gp = gpar(fontsize = legend_title_size, fontface = "bold")),
                                                  col=col_fun,
                                                  rect_gp = gpar(type = "none"),
                                                  layer_fun = layer_fun,
                                                  row_names_gp = gpar(fontsize = row_label_size, fontface = row_label_fontface),
                                                  column_names_gp = gpar(fontsize = column_label_size),
                                                  column_km = k,
                                                  row_km_repeats = ident_km_repeats,
                                                  border = "black",
                                                  left_annotation = column_ha,
                                                  column_km_repeats = feature_km_repeats,
                                                  show_parent_dend_line = show_parent_dend_line,
                                                  column_names_rot = x_lab_rotate,
                                                  cluster_rows = cluster_ident,
                                                  cluster_columns = cluster_feature)
    } else {
      cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                                  heatmap_legend_param=list(title="Expression", labels_gp = gpar(fontsize = legend_label_size), title_gp = gpar(fontsize = legend_title_size, fontface = "bold")),
                                                  col=col_fun,
                                                  rect_gp = gpar(type = "none"),
                                                  layer_fun = layer_fun,
                                                  row_names_gp = gpar(fontsize = row_label_size, fontface = row_label_fontface),
                                                  column_names_gp = gpar(fontsize = column_label_size),
                                                  row_km = k,
                                                  row_km_repeats = feature_km_repeats,
                                                  border = "black",
                                                  top_annotation = column_ha,
                                                  column_km_repeats = ident_km_repeats,
                                                  show_parent_dend_line = show_parent_dend_line,
                                                  column_names_rot = x_lab_rotate,
                                                  cluster_rows = cluster_feature,
                                                  cluster_columns = cluster_ident)
    }
  } else {
    if (isTRUE(x = flip)) {
      cluster_dot_plot <- ComplexHeatmap::Heatmap(t(exp_mat),
                                                  heatmap_legend_param=list(title="Expression", labels_gp = gpar(fontsize = legend_label_size), title_gp = gpar(fontsize = legend_title_size, fontface = "bold")),
                                                  col=col_fun,
                                                  rect_gp = gpar(type = "none"),
                                                  cell_fun = cell_fun_flip,
                                                  row_names_gp = gpar(fontsize = row_label_size, fontface = row_label_fontface),
                                                  column_names_gp = gpar(fontsize = column_label_size),
                                                  column_km = k,
                                                  row_km_repeats = ident_km_repeats,
                                                  border = "black",
                                                  left_annotation = column_ha,
                                                  column_km_repeats = feature_km_repeats,
                                                  show_parent_dend_line = show_parent_dend_line,
                                                  column_names_rot = x_lab_rotate,
                                                  cluster_rows = cluster_ident,
                                                  cluster_columns = cluster_feature)
    } else {
      cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                                  heatmap_legend_param=list(title="Expression", labels_gp = gpar(fontsize = legend_label_size), title_gp = gpar(fontsize = legend_title_size, fontface = "bold")),
                                                  col=col_fun,
                                                  rect_gp = gpar(type = "none"),
                                                  cell_fun = cell_fun,
                                                  row_names_gp = gpar(fontsize = row_label_size, fontface = row_label_fontface),
                                                  column_names_gp = gpar(fontsize = column_label_size),
                                                  row_km = k,
                                                  row_km_repeats = feature_km_repeats,
                                                  border = "black",
                                                  top_annotation = column_ha,
                                                  column_km_repeats = ident_km_repeats,
                                                  show_parent_dend_line = show_parent_dend_line,
                                                  column_names_rot = x_lab_rotate,
                                                  cluster_rows = cluster_feature,
                                                  cluster_columns = cluster_ident)
    }
  }

  # Add pt.size legend & return plots
  if (isTRUE(x = plot_km_elbow)) {
    return(list(km_elbow_plot, ComplexHeatmap::draw(cluster_dot_plot, annotation_legend_list = lgd_list)))
  }
  return(ComplexHeatmap::draw(cluster_dot_plot, annotation_legend_list = lgd_list))
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
#' @param seed Sets seed for reproducible plotting (ComplexHeatmap plot).
#'
#' @return A ComplexHeatmap or if plot_km_elbow = TRUE a list containing ggplot2 object and ComplexHeatmap.
#'
#' @import cli
#' @import ggplot2
#' @importFrom circlize colorRamp2
#' @importFrom dplyr any_of filter select pull
#' @importFrom grid grid.circle grid.rect gpar
#' @importFrom magrittr "%>%"
#' @importFrom rlang is_installed
#' @importFrom Seurat DotPlot
#' @importFrom stats quantile
#' @importFrom stringr str_to_lower
#' @importFrom tidyr pivot_wider
#'
#' @noRd
#'
#' @concept seurat_plotting
#'
#' @author Ming Tang (Original Code), Sam Marsh (Wrap single function, added/modified functionality)
#' @references \url{https://divingintogeneticsandgenomics.com/post/how-to-make-a-multi-group-dotplot-for-single-cell-rnaseq-data/}
#' @seealso \url{https://twitter.com/tangming2005}
#'
#' @examples
#' \donttest{
#' library(Seurat)
#' Clustered_DotPlot(seurat_object = pbmc_small, features = c("CD3E", "CD8", "GZMB", "MS4A1"))
#'}
#'

Clustered_DotPlot_Multi_Group <- function(
    seurat_object,
    features,
    split.by,
    colors_use_exp = viridis_plasma_dark_high,
    exp_color_min = -2,
    exp_color_middle = NULL,
    exp_color_max = 2,
    exp_value_type = "scaled",
    print_exp_quantiles = FALSE,
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
    seed = 123
) {
  # Check for packages
  ComplexHeatmap_check <- is_installed(pkg = "ComplexHeatmap")
  if (isFALSE(x = ComplexHeatmap_check)) {
    cli_abort(message = c(
      "Please install the {.val ComplexHeatmap} package to use {.code Clustered_DotPlot}",
      "i" = "This can be accomplished with the following commands: ",
      "----------------------------------------",
      "{.field `install.packages({symbol$dquote_left}BiocManager{symbol$dquote_right})`}",
      "{.field `BiocManager::install({symbol$dquote_left}ComplexHeatmap{symbol$dquote_right})`}",
      "----------------------------------------"
    ))
  }

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

  # set assay (if null set to active assay)
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Check expression value type
  accepted_exp_types <- c("scaled", "average")

  exp_value_type <- str_to_lower(string = exp_value_type)

  if (!exp_value_type %in% accepted_exp_types) {
    cli_abort(message = "{.code exp_value_type}, must be one of {.field {accepted_exp_types}}")
  }

  # Ignore exp_min and exp_max colors
  if (exp_value_type == "average") {
    if (exp_color_min != -2 || exp_color_max != 2 || !is.null(x = exp_color_middle)) {
      ignored_params <- c("exp_color_min", "exp_color_max", "exp_color_middle")
      cli_warn(message = c("One or more of the following parameters were set to a non-default value but are ignored when {.code exp_value_type = 'avergae'}",
                           "i" = "{.field {glue_collapse_scCustom(input_string = ignored_params, and = TRUE)}}."))
    }
  }

  # Check acceptable fontface
  if (!row_label_fontface %in% c("plain", "bold", "italic", "oblique", "bold.italic")) {
    cli_abort(message = c("{.code row_label_face} {.val {row_label_face}} not recognized.",
                          "i" = "Must be one of {.val plain}, {.val bold}, {.val italic}, {.val olique}, or {.val bold.italic}."))
  }

  # Check unique features
  features_unique <- unique(x = features)

  if (length(x = features_unique) != length(x = features)) {
    cli_warn("Feature list contains duplicates, making unique.")
  }

  # Check features and meta to determine which features present
  all_found_features <- Feature_PreCheck(object = seurat_object, features = features_unique, assay = assay)

  # Check exp min/max set correctly
  if (!exp_color_min < exp_color_max) {
    cli_abort(message = c("Expression color min/max values are not compatible.",
                          "i" = "The value for {.code exp_color_min}: {.field {exp_color_min}} must be less than the value for {.code exp_color_max}: {.field {exp_color_max}}.")
    )
  }

  # set group.by value
  group.by <- group.by %||% "ident"

  # Get data
  exp_mat_df <- suppressMessages(data.frame(AverageExpression(object = seurat_object, features = all_found_features, group.by = c(group.by, split.by), assays = assay, layer = "data")[[assay]]))

  # Data is returned in non-log space after averaging, return to log space for plotting
  exp_mat <- data.frame(lapply(exp_mat_df, function(x){
    log1p(x)
  }))

  exp_mat <- as.matrix(exp_mat)
  rownames(exp_mat) <- rownames(exp_mat_df)

  # scale data
  if (exp_value_type == "scaled") {
    exp_mat <- FastRowScale(mat = exp_mat)
    rownames(exp_mat) <- rownames(exp_mat_df)
  }

  # check underscore present in split.by and replace if so
  split_by_names <- Fetch_Meta(object = seurat_object) %>%
    select(any_of(split.by)) %>%
    pull()

  under_score <- grep(pattern = "_", x = split_by_names, value = TRUE)

  if (length(x = under_score) > 0) {
    split_by_names <- gsub(pattern = "_", replacement = ".", x = split_by_names)
    seurat_object[[split.by]] <- split_by_names
  }

  percent_mat <- Percent_Expressing(seurat_object = seurat_object, features = all_found_features, split_by = split.by, group_by = group.by, assay = assay)

  # reorder columns to match
  idx <- match(colnames(x = exp_mat), colnames(x = percent_mat))
  idx

  percent_mat <- percent_mat[, idx]
  percent_mat <- as.matrix(percent_mat)

  # print quantiles
  if (isTRUE(x = print_exp_quantiles)) {
    cli_inform(message = "Quantiles of gene expression data are:")
    print(quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99)))
  }

  # check grid color
  if (is.null(x = grid_color)) {
    grid_color <- NA
  } else {
    if (length(x = grid_color) > 1) {
      cli_abort(message = "{.code grid_color} can only be a single value.")
    }
    if (isTRUE(x = Is_Color(colors = grid_color))) {
      grid_color <- grid_color
    } else {
      cli_abort(message = "Value provided to {.code grid_color} ({.field {grid_color}}) is not valid value for color in R.")
    }
  }

  # Set middle of color scale if not specified
  if (exp_value_type == "scaled") {
    if (is.null(x = exp_color_middle)) {
      exp_color_middle <- Middle_Number(min = exp_color_min, max = exp_color_max)
    }

    palette_length <- length(x = colors_use_exp)
    palette_middle <- Middle_Number(min = 0, max = palette_length)

    # Create palette
    col_fun <-  colorRamp2(c(exp_color_min, exp_color_middle, exp_color_max), colors_use_exp[c(1,palette_middle, palette_length)])
  }

  if (exp_value_type == "average") {
    if (is.null(x = exp_color_middle)) {
      avg_color_max <- max(apply(exp_mat, 2, function(x) max(x, na.rm = TRUE)))
      avg_color_min <- 0
      avg_color_middle <- Middle_Number(min = 0, max = avg_color_max)

      palette_length <- length(x = colors_use_exp)
      palette_middle <- Middle_Number(min = 0, max = palette_length)

      # Create palette
      col_fun <- colorRamp2(c(avg_color_min, avg_color_middle, avg_color_max), colors_use_exp[c(1,palette_middle, palette_length)])

    }
  }

  # Calculate and plot Elbow
  if (isTRUE(x = plot_km_elbow)) {
    # if elbow_kmax not NULL check it is usable
    if (!is.null(x = elbow_kmax) && elbow_kmax > (nrow(x = exp_mat) - 1)) {
      elbow_kmax <- nrow(x = exp_mat) - 1
      cli_warn(message = c("The value provided for {.code elbow_kmax} is too large.",
                           "i" = "Changing to (length(x = features)-1): {.field {elbow_kmax}}.")
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
  if (isTRUE(x = flip)) {
    if (isTRUE(x = raster)) {
      layer_fun_flip = function(i, j, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(col = grid_color, fill = NA))
        grid.circle(x=x,y=y,r= sqrt(ComplexHeatmap::pindex(percent_mat, i, j)/100)  * unit(2, "mm"),
                    gp = gpar(fill = col_fun(ComplexHeatmap::pindex(exp_mat, i, j)), col = NA))
      }
    } else {
      cell_fun_flip = function(i, j, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(col = grid_color, fill = NA))
        grid.circle(x=x,y=y,r= sqrt(percent_mat[i, j]/100) * unit(2, "mm"),
                    gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))
      }
    }
  } else {
    if (isTRUE(x = raster)) {
      layer_fun = function(j, i, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(col = grid_color, fill = NA))
        grid.circle(x=x,y=y,r= sqrt(ComplexHeatmap::pindex(percent_mat, i, j)/100)  * unit(2, "mm"),
                    gp = gpar(fill = col_fun(ComplexHeatmap::pindex(exp_mat, i, j)), col = NA))
      }
    } else {
      cell_fun = function(j, i, x, y, w, h, fill) {
        grid.rect(x = x, y = y, width = w, height = h,
                  gp = gpar(col = grid_color, fill = NA))
        grid.circle(x=x,y=y,r= sqrt(percent_mat[i, j]/100) * unit(2, "mm"),
                    gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))
      }
    }
  }

  # Create legend for point size
  lgd_list = list(
    ComplexHeatmap::Legend(labels = c(10,25,50,75,100), title = "Percent Expressing",
                           graphics = list(
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.1) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.50) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(2, "mm"),
                                                              gp = gpar(fill = "black")),
                             function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                                              gp = gpar(fill = "black"))),
                           labels_gp = gpar(fontsize = legend_label_size),
                           title_gp = gpar(fontsize = legend_title_size, fontface = "bold")
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
  if (isTRUE(x = raster)) {
    if (isTRUE(x = flip)) {
      cluster_dot_plot <- ComplexHeatmap::Heatmap(t(exp_mat),
                                                  heatmap_legend_param=list(title="Expression", labels_gp = gpar(fontsize = legend_label_size), title_gp = gpar(fontsize = legend_title_size, fontface = "bold")),
                                                  col=col_fun,
                                                  rect_gp = gpar(type = "none"),
                                                  layer_fun = layer_fun,
                                                  row_names_gp = gpar(fontsize = row_label_size, fontface = row_label_fontface),
                                                  column_names_gp = gpar(fontsize = column_label_size),
                                                  column_km = k,
                                                  row_km_repeats = ident_km_repeats,
                                                  border = "black",
                                                  column_km_repeats = feature_km_repeats,
                                                  show_parent_dend_line = show_parent_dend_line,
                                                  column_names_rot = x_lab_rotate,
                                                  cluster_rows = cluster_ident,
                                                  cluster_columns = cluster_feature,
                                                  heatmap_width = heatmap_width,
                                                  heatmap_height = heatmap_height)
    } else {
      cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                                  heatmap_legend_param=list(title="Expression", labels_gp = gpar(fontsize = legend_label_size), title_gp = gpar(fontsize = legend_title_size, fontface = "bold")),
                                                  col=col_fun,
                                                  rect_gp = gpar(type = "none"),
                                                  layer_fun = layer_fun,
                                                  row_names_gp = gpar(fontsize = row_label_size, fontface = row_label_fontface),
                                                  column_names_gp = gpar(fontsize = column_label_size),
                                                  row_km = k,
                                                  row_km_repeats = feature_km_repeats,
                                                  border = "black",
                                                  column_km_repeats = ident_km_repeats,
                                                  show_parent_dend_line = show_parent_dend_line,
                                                  column_names_rot = x_lab_rotate,
                                                  cluster_rows = cluster_feature,
                                                  cluster_columns = cluster_ident,
                                                  heatmap_width = heatmap_width,
                                                  heatmap_height = heatmap_height)
    }
  } else {
    if (isTRUE(x = flip)) {
      cluster_dot_plot <- ComplexHeatmap::Heatmap(t(exp_mat),
                                                  heatmap_legend_param=list(title="Expression", labels_gp = gpar(fontsize = legend_label_size), title_gp = gpar(fontsize = legend_title_size, fontface = "bold")),
                                                  col=col_fun,
                                                  rect_gp = gpar(type = "none"),
                                                  cell_fun = cell_fun_flip,
                                                  row_names_gp = gpar(fontsize = row_label_size, fontface = row_label_fontface),
                                                  column_names_gp = gpar(fontsize = column_label_size),
                                                  column_km = k,
                                                  row_km_repeats = ident_km_repeats,
                                                  border = "black",
                                                  column_km_repeats = feature_km_repeats,
                                                  show_parent_dend_line = show_parent_dend_line,
                                                  column_names_rot = x_lab_rotate,
                                                  cluster_rows = cluster_ident,
                                                  cluster_columns = cluster_feature,
                                                  heatmap_width = heatmap_width,
                                                  heatmap_height = heatmap_height)
    } else {
      cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                                  heatmap_legend_param=list(title="Expression", labels_gp = gpar(fontsize = legend_label_size), title_gp = gpar(fontsize = legend_title_size, fontface = "bold")),
                                                  col=col_fun,
                                                  rect_gp = gpar(type = "none"),
                                                  cell_fun = cell_fun,
                                                  row_names_gp = gpar(fontsize = row_label_size, fontface = row_label_fontface),
                                                  column_names_gp = gpar(fontsize = column_label_size),
                                                  row_km = k,
                                                  row_km_repeats = feature_km_repeats,
                                                  border = "black",
                                                  column_km_repeats = ident_km_repeats,
                                                  show_parent_dend_line = show_parent_dend_line,
                                                  column_names_rot = x_lab_rotate,
                                                  cluster_rows = cluster_feature,
                                                  cluster_columns = cluster_ident,
                                                  heatmap_width = heatmap_width,
                                                  heatmap_height = heatmap_height)
    }
  }

  # Add pt.size legend & return plots
  if (isTRUE(x = plot_km_elbow)) {
    return(list(km_elbow_plot, ComplexHeatmap::draw(cluster_dot_plot, annotation_legend_list = lgd_list, merge_legend = TRUE)))
  }
  return(ComplexHeatmap::draw(cluster_dot_plot, annotation_legend_list = lgd_list, merge_legend = TRUE))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### TEST/HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Automatically calculate a point size for ggplot2-based scatter plots
#
#' It happens to look good
#'
#' @param data a single value length vector corresponding to the number of cells.
#' @param raster If TRUE, point size is set to 1
#'
#' @return The "optimal" point size for visualizing these data
#'
#' @noRd
#'
#' @references This function and documentation text are modified versions of the `AutoPointSize` function
#' and documentation from Seurat \url{https://github.com/satijalab/seurat/blob/master/R/visualization.R} (License: GPL-3).
#' This version has been modified to take single value length input instead of data.frame input.
#'

AutoPointSize_scCustom <- function(data, raster = NULL) {
  # for single value
  if (is.null(x = nrow(x = data)) && length(x = data) == 1 && is.numeric(x = data)) {
    return(ifelse(
      test = isTRUE(x = raster),
      yes = 1,
      no = min(1583 / data, 1)
    ))
  }
  if (inherits(what = "Seurat", x = data)) {

    return(ifelse(
      test = isTRUE(x = raster),
      yes = 1,
      no = min(1583 / length(x = Cells(x = data)), 1)
    ))
  } else {
    # for data frame/object based values (from Seurat, see documentation)
    return(ifelse(
      test = isTRUE(x = raster),
      yes = 1,
      no = min(1583 / nrow(x = data), 1)
    ))
  }
}


#' Extract max value for stacked violin plot
#'
#' extract max expression value
#'
#' @param p ggplot plot build to extract values from.
#'
#' @return max expression value
#'
#' @import ggplot2
#'
#' @noRd
#'
#' @author Ming Tang (Original Code), Sam Marsh (Modified function for use in scCustomtize)
#' @references \url{https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/}
#' @seealso \url{https://twitter.com/tangming2005}
#'

# extract the max value of the y axis
Extract_Max <- function(
  p
) {
  ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


#' Check integer
#'
#' Check if resulting number is integer
#'
#' @param x number to test if integer.
#'
#' @return logical value
#'
#' @noRd
#'
#' @author Iterator (StackOverflow)
#' @source  \url{https://stackoverflow.com/a/7798235}
#' @details https://creativecommons.org/licenses/by-sa/3.0/
#'

Test_Integer <- function(
  x
  ) {
  test <- all.equal(x, as.integer(x), check.attributes = FALSE)
  if (isTRUE(x = test)) {
    return(TRUE)
  } else {
      return(FALSE)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GGPLOT2/THEMES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Unrotate x axis on VlnPlot
#'
#' Shortcut for thematic modification to unrotate the x axis (e.g., for Seurat VlnPlot is rotated by default).
#'
#' @param ... extra arguments passed to `ggplot2::theme()`.
#'
#' @importFrom ggplot2 theme
#'
#' @export
#'
#' @return Returns a list-like object of class _theme_.
#'
#' @concept themes
#'
#' @examples
#' library(Seurat)
#' p <- VlnPlot(object = pbmc_small, features = "CD3E")
#' p + UnRotate_X()
#'

UnRotate_X <- function(...) {
  unrotate_x_theme <- theme(
    axis.text.x =
      element_text(angle = 0,
                   hjust = 0.5),
    validate = TRUE,
    ...
  )
  return(unrotate_x_theme)
}


#' Blank Theme
#'
#' Shortcut for thematic modification to remove all axis labels and grid lines
#'
#' @param ... extra arguments passed to `ggplot2::theme()`.
#'
#' @importFrom ggplot2 theme
#'
#' @export
#'
#' @return Returns a list-like object of class _theme_.
#'
#' @concept themes
#'
#' @examples
#' # Generate a plot and customize theme
#' library(ggplot2)
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + Blank_Theme()
#'

Blank_Theme <- function(...) {
  blank_theme <- theme(
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    validate = TRUE,
    ...
  )
  return(blank_theme)
}


#' Move Legend Position
#'
#' Shortcut for thematic modification to move legend position.
#'
#' @param position valid position to move legend.  Default is "right".
#' @param ... extra arguments passed to `ggplot2::theme()`.
#'
#' @importFrom ggplot2 theme
#'
#' @export
#'
#' @return Returns a list-like object of class _theme_.
#'
#' @concept themes
#'
#' @examples
#' # Generate a plot and customize theme
#' library(ggplot2)
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + Move_Legend("left")
#'

Move_Legend <- function(
  position = "right",
  ...
) {
  move_legend_theme <- theme(
    legend.position = position,
    validate = TRUE,
    ...
  )
  return(move_legend_theme)
}


#' Modified ggprism theme
#'
#' Modified ggprism theme which restores the legend title.
#'
#' @param palette `string`. Palette name, use
#' `names(ggprism_data$themes)` to show all valid palette names.
#' @param base_size `numeric`. Base font size, given in `"pt"`.
#' @param base_family `string`. Base font family, default is `"sans"`.
#' @param base_fontface `string`. Base font face, default is `"bold"`.
#' @param base_line_size `numeric`. Base linewidth for line elements
#' @param base_rect_size `numeric`. Base linewidth for rect elements
#' @param axis_text_angle `integer`. Angle of axis text in degrees.
#' One of: `0, 45, 90, 270`.
#' @param border `logical`. Should a border be drawn around the plot?
#' Clipping will occur unless e.g. `coord_cartesian(clip = "off")` is used.
#'
#' @references theme is a modified version of `theme_prism` from ggprism package \url{https://github.com/csdaw/ggprism}
#' (License: GPL-3).  Param text is from `ggprism:theme_prism()` documentation \code{\link[ggprism]{theme_prism}}.
#' Theme adaptation based on ggprism vignette
#' \url{https://csdaw.github.io/ggprism/articles/themes.html#make-your-own-ggprism-theme-1}.
#'
#' @import ggplot2
#' @importFrom ggprism theme_prism
#'
#' @export
#'
#' @return Returns a list-like object of class _theme_.
#'
#' @concept themes
#'
#' @examples
#' # Generate a plot and customize theme
#' library(ggplot2)
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + theme_ggprism_mod()
#'

theme_ggprism_mod <- function(
  palette = "black_and_white",
  base_size = 14,
  base_family = "sans",
  base_fontface = "bold",
  base_line_size = base_size / 20,
  base_rect_size = base_size / 20,
  axis_text_angle = 0,
  border = FALSE
) {
  theme_prism(palette = palette,
              base_size = base_size,
              base_family = base_family,
              base_fontface = base_fontface,
              base_line_size = base_line_size,
              base_rect_size = base_rect_size,
              axis_text_angle = axis_text_angle,
              border = border) %+replace%
    theme(legend.title = element_text(hjust = 0),
          axis.text = element_text(size = rel(0.95), face = "plain")
    )
}


#' Remove Right Y Axis
#'
#' Shortcut for removing right y axis from ggplot2 object
#'
#' @importFrom ggplot2 theme
#'
#' @references Shortcut slightly modified from Seurat \url{https://github.com/satijalab/seurat/blob/c4638730d0639d770ad12c35f50d19108e0491db/R/visualization.R#L1039-L1048}
#'
#' @keywords internal
#'
#' @return Returns a list-like object of class _theme_.
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' # Generate a plot without axes, labels, or grid lines
#' library(ggplot2)
#' p <- FeaturePlot(object = obj, features = "Cx3cr1")
#' p + No_Right()
#' }

No_Right <- function() {
  no.right <- theme(
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y.right = element_text(
      face = "bold",
      size = 14,
      margin = margin(r = 7),
      angle = 270
    )
  )
  return(no.right)
}


