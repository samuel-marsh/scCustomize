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
  mod_theme <- theme_prism(palette = palette,
                           base_size = base_size,
                           base_family = base_family,
                           base_fontface = base_fontface,
                           base_line_size = base_line_size,
                           base_rect_size = base_rect_size,
                           axis_text_angle = axis_text_angle,
                           border = border)

  mod_theme %+replace% theme(legend.title = element_text(hjust = 0),
                             axis.text = element_text(size = rel(0.95), face = "plain")
                             )

  mod_theme[c("legend.text.align", "legend.title.align")] <- NULL

  return(mod_theme)
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
