#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### VIRIDIS SHORTCUTS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Viridis Shortcuts
#'
#' Quick shortcuts to access viridis palettes
#'
#' @return A color palette for plotting
#'
#' @importFrom paletteer paletteer_c
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'
#' @examples
#' \dontrun{
#' FeaturePlot_scCustom(object = seurat_object, features = "Cx3cr1",
#' colors_use = viridis_plasma_dark_high, na_color = "lightgray")
#' }
#'

viridis_plasma_dark_high <- as.vector(x = paletteer_c(palette = "viridis::plasma", n = 250, direction = -1))

#' Viridis Shortcuts
#'
#' @importFrom paletteer paletteer_c
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'

viridis_plasma_light_high <- as.vector(x = paletteer_c(palette = "viridis::plasma", n = 250, direction = 1))

#' Viridis Shortcuts
#'
#' @importFrom paletteer paletteer_c
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'

viridis_inferno_dark_high <- as.vector(x = paletteer_c(palette = "viridis::inferno", n = 250, direction = -1))

#' Viridis Shortcuts
#'
#' @importFrom paletteer paletteer_c
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'

viridis_inferno_light_high <- as.vector(x = paletteer_c(palette = "viridis::inferno", n = 250, direction = 1))

#' Viridis Shortcuts
#'
#' @importFrom paletteer paletteer_c
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'

viridis_magma_dark_high <- as.vector(x = paletteer_c(palette = "viridis::magma", n = 250, direction = -1))

#' Viridis Shortcuts
#'
#' @importFrom paletteer paletteer_c
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'

viridis_magma_light_high <- as.vector(x = paletteer_c(palette = "viridis::magma", n = 250, direction = 1))

#' Viridis Shortcuts
#'
#' @importFrom paletteer paletteer_c
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'

viridis_dark_high <- as.vector(x = paletteer_c(palette = "viridis::viridis", n = 250, direction = -1))

#' Viridis Shortcuts
#'
#' @importFrom paletteer paletteer_c
#'
#' @export
#'
#' @concept palettes
#' @rdname viridis_shortcut
#'

viridis_light_high <- as.vector(x = paletteer_c(palette = "viridis::viridis", n = 250, direction = 1))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### PALETTE FUNCTIONS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Single Color Palettes for Plotting
#'
#' Selects colors from modified versions of RColorBrewer single colors palettes
#'
#' @param pal_color color palette to select (Options are: 'reds', 'blues',
#' 'greens', 'purples', 'oranges', 'grays').
#' @param num_colors set number of colors (max = 7).
#' @param seed_use set seed for reproducibility (default: 123).
#'
#' @return A vector of colors
#'
#' @references
#' See RColorBrewer for more info on palettes
#' \url{https://CRAN.R-project.org/package=RColorBrewer}
#'
#' @import cli
#'
#' @export
#'
#' @concept palettes
#'
#' @examples
#' pal <- Single_Color_Palette(pal_color = "reds", num_colors = 7)
#' PalettePlot(pal= pal)
#'

Single_Color_Palette <- function(pal_color,
                                 num_colors = NULL,
                                 seed_use = 123
) {
  # Check number of colors available
  if (is.null(x = num_colors)) {
    cli_abort(message = "No value provided to {.code num_colors}.")
  }
  if (num_colors > 7 || num_colors < 1) {
    cli_abort(message = c("Not enough colors.",
                          "i" = "Value provided to {.code num_colors} ({.field {num_colors}}) is greater than maximum number allowed ({.field 7}).")
    )
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
  if (!pal_color %in% names(x = brewer_single_modified)) {
    cli_abort(message = c("Paleete name not found.",
                          "i" = "Palette name not found.  Please select one of following palette options: {.field 'reds', 'blues', 'greens', 'purples', or 'grays'}")
    )
  }
  set.seed(seed = seed_use)
  pal_use <- brewer_single_modified[[pal_color]]
  output_pal <- sample(x = pal_use, size = num_colors)
  return(output_pal)
}


#' Navy and Orange Dual Color Palette
#'
#' Shortcut to navy orange color plot
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
#' cols <- NavyAndOrange()
#' PalettePlot(pal= cols)
#'

NavyAndOrange <- function(
  flip_order = FALSE
) {
  navy_orange <- c("navy", "orange")
  if (isTRUE(x = flip_order)) {
    navy_orange <- rev(x = navy_orange)
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
#' \url{https://CRAN.R-project.org/package=ggsci}
#'
#' @export
#'
#' @concept palettes
#'
#' @examples
#' cols <- JCO_Four()
#' PalettePlot(pal= cols)
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
#' cols <- Dark2_Pal()
#' PalettePlot(pal= cols)
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
#' @param num_colors number of colors to return in palette.
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
#' cols <- Hue_Pal(num_colors = 8)
#' PalettePlot(pal= cols)
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
#' friendly palette \url{https://jfly.uni-koeln.de/color/}.
#'
#' @export
#'
#' @concept palettes
#'
#' @examples
#' cols <- ColorBlind_Pal()
#' PalettePlot(pal = cols)
#'

ColorBlind_Pal <- function(
) {
  color_blind_pal <- c("orange", "#0072B2","#009E73", "#CC79A7",
                       "#F0E442", "firebrick2", "#56B4E9", "gray55")

  return(color_blind_pal)
}


#' Generate a rainbow palette with variation in saturation and value
#'
#' @param n_colors The number of colors to generate
#'
#' @importFrom grDevices rainbow
#'
#' @return a character vector of hex color values of length n_colors.
#'
#' @references Ported from colorway package for CRAN release.  See \url{https://github.com/hypercompetent/colorway/blob/master/R/palettes.R} (License: GPL-3).
#'
#' @keywords internal
#'
#' @noRd
#'

varibow_scCustom <- function(
  n_colors
) {
  sats <- rep_len(x = c(0.55,0.7,0.85,1), length.out = n_colors)
  vals <- rep_len(x = c(1,0.8,0.6), length.out = n_colors)
  rainbow(n_colors, s = sats, v = vals)
}


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
# #' @importFrom colorway varibow (now directly ported for CRAN compatibility)
#' @importFrom paletteer paletteer_d
#' @importFrom rlang is_installed
#'
#' @return A vector of colors
#'
#' @references
#' This function uses the paletteer package \url{https://github.com/EmilHvitfeldt/paletteer} to
#' provide simplified access to color palettes from many different R package sources while
#' minimizing scCustomize current and future dependencies.
#'
#' The following packages & palettes are called by this function (see individual packages for
#' palette references/citations):
#' \enumerate{
#'   \item pals (via paletteer) \url{https://CRAN.R-project.org/package=pals}
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
#' Function name and implementation modified from Seurat (License: GPL-3).
#' \url{https://github.com/satijalab/seurat}
#'
#' @export
#'
#' @concept palettes
#'
#' @examples
#' pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "varibow")
#' PalettePlot(pal= pal)
#'

DiscretePalette_scCustomize <- function(
  num_colors,
  palette = NULL,
  shuffle_pal = FALSE,
  seed = 123
) {
  if (is.null(x = palette)) {
    cli_abort(message = c("Must specify a palette to return colors.",
                          "i" = "{.code palette} options are: {.field {names(palette_list)}}")
    )
  }

  # dittoseq check
  if (palette == "ditto_seq") {
    dittoseq_check <- is_installed(pkg = "dittoSeq")
    if (isFALSE(x = dittoseq_check[1])) {
      cli_abort(message = c(
        "Please install the {.val dittoSeq} package to {.code palette = {symbol$dquote_left}ditto_seq{symbol$dquote_right}}",
        "i" = "This can be accomplished with the following commands:",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}BiocManager{symbol$dquote_right})`}",
        "{.field `BiocManager::install({symbol$dquote_left}dittoSeq{symbol$dquote_right})`}",
        "----------------------------------------"
      ))
    } else {
      palette_list <- list(
        alphabet = as.vector(x = paletteer_d("pals::alphabet", 26)),
        alphabet2 = as.vector(x = paletteer_d("pals::alphabet2", 26)),
        glasbey = as.vector(x = paletteer_d("pals::glasbey", 32)),
        polychrome = as.vector(x = paletteer_d("pals::polychrome", 36)),
        stepped = as.vector(x = paletteer_d("pals::stepped", 24)),
        ditto_seq =  dittoSeq::dittoColors(reps = 1, get.names = FALSE),
        varibow = varibow_scCustom(n_colors = num_colors)
      )
    }
  } else {
    palette_list <- list(
      alphabet = as.vector(x = paletteer_d("pals::alphabet", 26)),
      alphabet2 = as.vector(x = paletteer_d("pals::alphabet2", 26)),
      glasbey = as.vector(x = paletteer_d("pals::glasbey", 32)),
      polychrome = as.vector(x = paletteer_d("pals::polychrome", 36)),
      stepped = as.vector(x = paletteer_d("pals::stepped", 24)),
      varibow = varibow_scCustom(n_colors = num_colors)
    )
  }

  palette_out <- palette_list[[palette]]
  if (num_colors > length(x = palette_out)) {
    cli_abort(message = c("Not enough colors in specified palette.",
                          "*" = "{.val {palette}} only contains {.field {length(x = palette_out)}} colors.",
                          "i" = "Please adjust {.code num_colors} to be less than or equal to {.field {length(x = palette_out)}} or select a different {.code palette}.")
    )
  }
  if (isTRUE(x = shuffle_pal)) {
    set.seed(seed = seed)
    palette_out <- sample(x = palette_out[1:num_colors])
  } else {
    palette_out <- palette_out[1:num_colors]
  }
  return(palette_out)
}


#' Color Palette Selection for scCustomize
#'
#' Function to return package default discrete palettes depending on number of groups plotted.
#'
#' @param num_groups number of groups to be plotted. If `ggplot_default_colors = FALSE` then by default:
#' \itemize{
#'       \item If number of levels plotted equal to 2 then colors will be `NavyAndOrange()`.
#'       \item If number of levels plotted greater than 2 but less than or equal to 8 it will use `ColorBlind_Pal()`.
#'       \item If number of levels plotted greater than 2 but less than or equal to 36 it will use "polychrome" from `DiscretePalette_scCustomize()`.
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
#' cols <- scCustomize_Palette(num_groups = 24, ggplot_default_colors = FALSE)
#' PalettePlot(pal= cols)
#'

scCustomize_Palette <- function(
  num_groups,
  ggplot_default_colors = FALSE,
  color_seed = 123
) {
  # Set color palette depending on group length
  if (isTRUE(x = ggplot_default_colors)) {
    colors_use <- Hue_Pal(num_colors = num_groups)
  } else {
    if (num_groups == 1) {
      colors_use <- "dodgerblue"
    }
    if (num_groups == 2) {
      colors_use <- NavyAndOrange()
    }
    if (num_groups > 2 && num_groups <= 8) {
      colors_use <- ColorBlind_Pal()
    }
    if (num_groups > 8 && num_groups <= 36) {
      colors_use <- DiscretePalette_scCustomize(num_colors = num_groups, palette = "polychrome")
    }
    if (num_groups > 36) {
      colors_use <- DiscretePalette_scCustomize(num_colors = num_groups, palette = "varibow", shuffle_pal = TRUE, seed = color_seed)
    }
  }
  return(colors_use)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### PALETTE PLOTTING ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Plot color palette in viewer
#'
#' Plots given color vector/palette in viewer to evaluate palette before plotting on data.
#'
#' @param pal a vector of colors (either named colors of hex codes).
#' @param label_color_num logical, whether or not to numerically label the colors in output plot.
#' Default is TRUE is number of colors in `pal` is less than 75 and FALSE is greater than 75.
#'
#' @import cli
#' @import ggplot2
#'
#' @return Plot of all colors in supplied palette/vector
#'
#' @references
#' Adapted from colorway package `build_palette` internals (License: GPL-3).
#' \url{https://github.com/hypercompetent/colorway}.
#'
#' @export
#'
#' @concept palettes
#'
#' @examples
#' pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "varibow")
#' PalettePlot(pal = pal)
#'

PalettePlot <- function(
  pal = NULL,
  label_color_num = NULL
) {
  # Check palette
  if (is.null(x = pal)) {
    cli_abort(message = "No value provided to {.code palette} parameter.")
  }

  if (inherits(x = pal, what = "colors")) {
    pal <- as.vector(x = pal)
  }

  # Generate data frame for plotting
  palette_data <- data.frame(x = 1:length(pal), y = 1, fill = pal)

  # Decide color labeling
  if (is.null(x = label_color_num)) {
    if (length(x = pal) > 75) {
      label_color_num <- FALSE
    } else {
      label_color_num <- TRUE
    }
  }

  # Plot
  # Label plot
  if (isTRUE(x = label_color_num)) {
    palette_plot <- ggplot(palette_data) +
      geom_tile(aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["fill"]])) +
      geom_text(aes(x = .data[["x"]], y = .data[["y"]], label = .data[["x"]])) +
      scale_fill_identity() +
      theme_void()
  } else {
    palette_plot <- ggplot(palette_data) +
      geom_tile(aes(x = .data[["x"]], y = .data[["y"]], fill = .data[["fill"]])) +
      scale_fill_identity() +
      theme_void()
  }

  return(palette_plot)
}
