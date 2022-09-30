#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### MODDED PLOTS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' PC Plots
#'
#' Plot PC Heatmaps and Dim Loadings for exploratory analysis.  Plots a single Heatmap and Gene Loading Plot.
#'   Used for PC_Loading_Plots function.
#'
#' @param seurat_object Seurat Object.
#' @param dims_number A single dim to plot (integer).
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
#' \dontrun{
#' PC_Plotting(seurat_object = seurat, dims_plot = 25, "plots/")
#' }
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
#' @references https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/.  Solution for re-enabling plot spacing modification by Abdenour ABBAS (comment on original blog post; http://disq.us/p/2b54qh2).
#' @seealso https://twitter.com/tangming2005
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
#' @importFrom tibble rownames_to_column
#'
#' @noRd
#'
#' @references Code to calculate wss values from: https://stackoverflow.com/a/15376462/15568251
#'

kMeans_Elbow <- function(
  data,
  k_max = 15,
  plot_title = "Sum of Squared Error (SSE) Plot",
  cutoff_value = NULL
) {
  # Calculate the within squares
  # code from @Ben https://stackoverflow.com/a/15376462/15568251
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:k_max) wss[i] <- sum(kmeans(data,
                                          centers=i)$withinss)

  # Reformat for ggplot2 plotting
  plot_data <- data.frame(wss) %>%
    rownames_to_column("k")

  plot_data$k <- as.numeric(plot_data$k)

  # Plot data
  plot <- ggplot(data = plot_data, mapping = aes(y = wss, x = k)) +
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### TEST/HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
#' @references https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
#' @seealso https://twitter.com/tangming2005
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
#' @source  https://stackoverflow.com/a/7798235
#' @details https://creativecommons.org/licenses/by-sa/3.0/
#'

Test_Integer <- function(
  x
  ) {
  test <- all.equal(x, as.integer(x), check.attributes = FALSE)
  if (test == TRUE) {
    return(TRUE)
  } else {
      return(FALSE)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GGPLOT2/THEMES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Unrotate x axis on VlnPLot
#'
#' Shortcut for thematic modification to unrotate the x axis (e.g., for Seurat VlnPlot is rotated by default).
#'
#' @importFrom ggplot2 theme
#'
#' @export
#'
#' @concept themes
#'
#' @examples
#' \dontrun{
#' # Generate a plot and unrotate the x-axis label
#' library(ggplot2)
#' p <- VlnPlot(object = obj, features = "Cx3cr1")
#' p + UnRotate_X
#' }
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
#' @importFrom ggplot2 theme
#'
#' @export
#'
#' @concept themes
#'
#' @examples
#' \dontrun{
#' # Generate a plot without axes, labels, or grid lines
#' library(ggplot2)
#' p <- FeaturePlot(object = obj, features = "Cx3cr1")
#' p + Blank_Theme()
#' }

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
#'
#' @importFrom ggplot2 theme
#'
#' @export
#'
#' @concept themes
#'
#' @examples
#' \dontrun{
#' # Generate a plot and move the legend.
#' library(ggplot2)
#' p <- VlnPlot(object = obj, features = "Cx3cr1")
#' p + Move_Legend("left")
#' }
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
#' @references theme is a modified version of `theme_prism` from ggprism package (https://github.com/csdaw/ggprism)
#' (Licence: GPL-3).  Param text is from `ggprism:theme_prism()` documentation \code{\link[ggprism]{theme_prism}}.
#' Theme adaptation based on ggprism vignette
#' (https://csdaw.github.io/ggprism/articles/themes.html#make-your-own-ggprism-theme-1).
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
#' \dontrun{
#' library(ggplot2)
#' p <- FeaturePlot(object = obj, features = "Cx3cr1")
#' p + theme_ggprism_mod()
#' }

theme_ggprism_mod <- function(base_size = 14,
                              base_family = "sans",
                              base_fontface = "bold",
                              base_line_size = base_size / 20,
                              base_rect_size = base_size / 20,
                              axis_text_angle = 0,
                              border = FALSE) {
  theme_prism( base_size = base_size,
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
