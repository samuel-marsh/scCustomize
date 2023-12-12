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
      no = min(1583 / Cells(x = data), 1)
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


