#' Viridis Shortcuts
#'
#' Quick shortcuts to access viridis palettes
#'
#' @return A color palette for plotting
#'
#' @import viridis
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'
#' @examples
#' \dontrun{
#' FeaturePlot_scCustom(object = seurat_object, features = "Cx3cr1", colors_use = viridis_plasma_dark_high, na_color = "lightgray")
#' }
#'

viridis_plasma_dark_high <- viridis(n = 30, option = "C", direction = -1)

#' Viridis Shortcuts
#'
#' @import viridis
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'

viridis_plasma_light_high <- viridis(n = 30, option = "C", direction = 1)

#' Viridis Shortcuts
#'
#' @import viridis
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'

viridis_inferno_dark_high <- viridis(n = 30, option = "B", direction = -1)

#' Viridis Shortcuts
#'
#' @import viridis
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'

viridis_inferno_light_high <- viridis(n = 30, option = "B", direction = 1)

#' Viridis Shortcuts
#'
#' @import viridis
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'

viridis_magma_dark_high <- viridis(n = 30, option = "A", direction = -1)

#' Viridis Shortcuts
#'
#' @import viridis
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'

viridis_magma_light_high <- viridis(n = 30, option = "A", direction = 1)

#' Viridis Shortcuts
#'
#' @import viridis
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'

viridis_dark_high <- viridis(n = 30, option = "D", direction = -1)

#' Viridis Shortcuts
#'
#' @import viridis
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'

viridis_light_high <- viridis(n = 30, option = "D", direction = 1)


#' Discrete color palettes
#'
#' Helper function to return a number of discrete color palettes.
#'
#' @param num_colors Number of colors to be generated.
#' @param palette Options are
#' "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", "ditto_seq", "varibow".
#' @param shuffle_pal randomly shuffle the outputted palette.  Most useful for `varibow` palette as
#' that is normally an ordered palette.
#' @param seed random seed for the palette shuffle.  Default = 123.
#'
#' @import cli
#' @importFrom colorway varibow
#' @importFrom dittoSeq dittoColors
#' @importFrom paletteer paletteer_d
#'
#' @return A vector of colors
#'
#' @references
#' This function uses the paletteer package (https://github.com/EmilHvitfeldt/paletteer) to
#' provide simplified access to color palettes from many different R package sources while
#' minimizing scCustomize current and future dependencies.
#'
#' The following packages & palettes are called by this function (see individual packages for
#' palette references/citations):
#' \enumerate{
#'   \item pals (via paletteer) \url{https://cran.r-project.org/web/packages/pals/index.html}
#'     \itemize{
#'       \item alphabet, alphabet2, glasbey, polychrome, and stepped.
#'       }
#'   \item dittoSeq \url{https://bioconductor.org/packages/release/bioc/html/dittoSeq.html}
#'     \itemize{
#'       \item dittoColors.
#'       }
#'   \item colorway \url{https://github.com/hypercompetent/colorway}
#'     \itemize{
#'       \item varibow
#'       }
#' }
#'
#' Function name and implementation modified from Seurat (Licence: GPL-3).
#' \url{https://github.com/satijalab/seurat}
#'
#' @export
#'
#' @concept palettes
#'

DiscretePalette_scCustomize <- function(
  num_colors,
  palette = NULL,
  shuffle_pal = FALSE,
  seed = 123
) {
  palette_list <- list(
    alphabet = as.vector(x = paletteer_d("pals::alphabet", 26)),
    alphabet2 = as.vector(x = paletteer_d("pals::alphabet2", 26)),
    glasbey = as.vector(x = paletteer_d("pals::glasbey", 32)),
    polychrome = as.vector(x = paletteer_d("pals::polychrome", 36)),
    stepped = as.vector(x = paletteer_d("pals::stepped", 24)),
    ditto_seq =  dittoColors(reps = 1, get.names = FALSE),
    varibow = varibow(n_colors = num_colors)
  )
  if (is.null(x = palette)) {
    cli_abort(message = c("Must specify a palette to return colors.",
                          "i" = "`palette` options are: {names(palette_list)}")
    )
  }
  palette_out <- palette_list[[palette]]
  if (num_colors > length(x = palette_out)) {
    warning("Not enough colors in specified palette.  ", '"', palette, '"', " only has ", length(x = palette_out), " colors.")
  }
  if (shuffle_pal) {
    set.seed(seed = seed)
    sample(palette_out[1:num_colors])
  } else {
    palette_out[1:num_colors]
  }
}


#' Single Color Palettes for Plotting
#'
#' Selects colors from modified versions of RColorBrewer single colors palettes
#'
#' @param pal_color color palette to select (Options are: 'reds', 'blues',
#' 'greens', 'purples', 'oranges', 'grays').
#' @param num_colors set number of colors (max = 7).
#' @param seed set seed for reproducibility (default: 123).
#'
#' @return A vector of colors
#'
#' @references
#' See RColorBrewer for more info on palettes
#' \url{https://CRAN.R-project.org/package=RColorBrewer}
#'
#' @export
#'
#' @concept palettes
#'

Single_Color_Palette <- function(pal_color,
                                 num_colors = NULL,
                                 seed_use = 123
) {
  # Check number of colors available
  if (is.null(x = num_colors)) {
    stop("No value provided to `num_colors`.")
  }
  if (num_colors > 7 || num_colors < 1) {
    stop("Value provided to `num_colors` (", num_colors, ") is greater than maximum number allowed (7).")
  }

  # Modified single palettes
  brewer_single_modified <- list(
    reds = c(
      "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15",
      "#67000D"
    ),
    blues = c(
      "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C",
      "#08306B"
    ),
    greens = c(
      "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C",
      "#00441B"
    ),
    purples = c(
      "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F",
      "#3F007D"
    ),
    oranges = c(
      "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603",
      "#7F2704"
      ),
    grays = c(
      "#D9D9D9", "#BDBDBD", "#969696", "#737373", "#525252", "#252525",
      "#000000"
    )
  )
  if (!pal_color %in% names(brewer_single_modified)) {
    stop("Palette name not found.  Please select one of following palette options:
         'reds', 'blues', 'greens', 'purples', or 'grays'")
  }
  set.seed(seed = seed_use)
  pal_use <- brewer_single_modified[[pal_color]]
  output_pal <- sample(pal_use, size = num_colors)
  return(output_pal)
}


#' Plot color palette in viewer
#'
#' Plots given color vector/palette in viewer to evaluate palette before plotting on data.
#'
#' @param palette a vector of colors (either named colors of hex codes).
#'
#' @import ggplot2
#'
#' @return Plot of all colors in supplied palette/vector
#'
#' @references
#' Adapted from colorway package `build_palette` internals (Licence: GPL-3).
#' \url{https://github.com/hypercompetent/colorway}.
#'
#' @export
#'
#' @concept palettes
#'
#' @examples
#' \dontrun{
#' PalettePlot(palette = varibow(n_colors = 36))
#' }
#'

PalettePlot <- function(palette = NULL) {
  # Check palette
  if (is.null(x = palette)) {
    stop("No value provided to `palette` parameter.")
  }

  # Generate data frame for plotting
  palette_data <- data.frame(x = 1:length(palette), y = 1, fill = palette)

  # Plot
  palette_plot <- ggplot(palette_data) +
    geom_tile(aes(x = x, y = y, fill = fill)) +
    geom_text(aes(x = x, y = y, label = x)) +
    scale_fill_identity() +
    theme_void()
  return(palette_plot)
}


#' Navy and Orange Dual Color Palette
#'
#' Shortcut to navy ornage color plot
#'
#' @param flip_order Whether to flip the order of colors.
#'
#' @return Navy orange palette
#'
#' @export
#'
#' @concept palettes
#'
#' @examples
#' \dontrun{
#' cols <- NavyAndOrange()
#' }
#'

NavyAndOrange <- function(
  flip_order = F
) {
  navy_orange <- c("navy", "orange")
  if (flip_order) {
    navy_orange <- rev(navy_orange)
  }
  return(navy_orange)
}


#' Four Color Palette (JCO)
#'
#' Shortcut to a specific JCO 4 color palette from ggsci package.
#'
#' @importFrom paletteer paletteer_d
#'
#' @return 4 color palette from the JCO ggsci palette
#'
#' @references
#' Selection of colors from the JCO palette from ggsci being called through paletteer.
#' See ggsci for more info on palettes
#' \url{https://cran.r-project.org/web/packages/ggsci/index.html}
#'
#' @export
#'
#' @concept palettes
#'
#' @examples
#' \dontrun{
#' cols <- JCO_Four()
#' }
#'

JCO_Four <- function(
) {
  jco <- as.vector(x = paletteer_d("ggsci::default_jco", 10))
  jco_four <- c(jco[1], jco[2], jco[4], jco[8])

  return(jco_four)
}


#' Dark2 Palette
#'
#' Shortcut to Dark2 color palette from RColorBrewer (8 Colors)
#'
#' @importFrom paletteer paletteer_d
#'
#' @return "Dark2" color palette (8 colors)
#'
#' @references
#' Dark2 palette from RColorBrewer being called through paletteer.
#' See RColorBrewer for more info on palettes
#' \url{https://CRAN.R-project.org/package=RColorBrewer}
#'
#' @export
#'
#' @concept palettes
#'
#' @examples
#' \dontrun{
#' cols <- Dark2_Pal()
#' }
#'

Dark2_Pal <- function(
) {
  dark2 <- as.vector(x = paletteer_d("RColorBrewer::Dark2", 8))

  return(dark2)
}


#' Hue_Pal
#'
#' Shortcut to hue_pal to return to ggplot2 defaults if user desires, from scales package.
#'
#' @param num_colors
#'
#' @importFrom scales hue_pal
#'
#' @return hue color palette (as many colors as desired)
#'
#' @export
#'
#' @concept palettes
#'
#' @examples
#' \dontrun{
#' cols <- Hue_Pal(num_colors = 8)
#' }
#'

Hue_Pal <- function(
  num_colors
) {
  pal <- hue_pal()(num_colors)

  return(pal)
}


#' Color Universal Design Short Palette
#'
#' Shortcut ta a modified 8 color palette based on Color Universal Design (CUD) colorblindness friendly palette.
#'
#' @return modified/reordered color palette (8 colors) based on ditto-seq
#'
#' @references palette is slightly modified version of the Color Universal Design (CUD) colorblindness
#' friendly palette (https://jfly.uni-koeln.de/color/).
#'
#' @export
#'
#' @concept palettes
#'
#' @examples
#' \dontrun{
#' cols <- ColorBlind_Pal()
#' }
#'

ColorBlind_Pal <- function(
) {
  color_blind_pal <- c("orange", "#0072B2","#009E73", "#CC79A7",
                       "#F0E442", "firebrick2", "#56B4E9", "gray55")

  return(color_blind_pal)
}
"firebrick1"


#' Color Palette Selection for scCustomize
#'
#' Function to return package default discrete palettes depending on number of groups plotted.
#'
#' @param num_groups number of groups to be plotted. If `ggplot_default_colors = FALSE` then by default:
#' \itemize{
#'       \item If number of levels plotted equal to 2 then colors will be `NavyAndOrange()`.
#'       \item If If number of levels plotted greater than 2 but less than or equal to 36 it will use "polychrome" from `DiscretePalette_scCustomize`.
#'       \item If greater than 36 will use "varibow" with shuffle = TRUE from `DiscretePalette_scCustomize`.
#'       }
#' @param ggplot_default_colors logical.  Whether to use default ggplot hue palette or not.
#' @param color_seed random seed to use for shuffling the "varibow" palette.
#'
#' @return vector of colors to use for plotting.
#'
#' @export
#'
#' @concept palettes
#'
#' @examples
#' \dontrun{
#' cols <- scCustomize_Palette(num_groups = 24, ggplot_default_colors = FALSE)
#' }
#'

scCustomize_Palette <- function(
  num_groups,
  ggplot_default_colors,
  color_seed = 123
) {
  # Set color palette depending on group length
  if (ggplot_default_colors) {
    colors_use <- Hue_Pal(num_colors = num_groups)
  } else {
    if (num_groups == 2) {
      colors_use <- NavyAndOrange()
    }
    if (num_groups > 2 && num_groups <= 36) {
      colors_use <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
    }
    if (num_groups > 36) {
      colors_use <- DiscretePalette_scCustomize(num_colors = num_groups, palette = "varibow", shuffle_pal = TRUE, seed = color_seed)
    }
  }
  return(colors_use)
}
