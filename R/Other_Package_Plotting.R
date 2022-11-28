#' Nebulosa Density Plot
#'
#' Allow for customization of Nebulosa plot_density.  Requires Nebulosa package from Bioconductor.
#'
#' @param seurat_object Seurat object name.
#' @param features Features to plot.
#' @param joint logical. Whether to return joint density plot. Default is FALSE.
#' @param viridis_palette default viridis palette to use (must be one of: "viridis", "magma", "cividis",
#' "inferno", "plasma").  Default is "magma".
#' @param custom_palette non-default color palette to be used in place of default viridis options.
#' @param pt.size Adjust point size for plotting.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default).
#' @param combine Create a single plot? If FALSE, a list with ggplot objects is returned.
#' @param ... Extra parameters passed to \code{\link[Nebulosa]{plot_density}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @import patchwork
# #' @importFrom Nebulosa plot_density
#' @importFrom SeuratObject DefaultDimReduc PackageCheck
#'
#' @export
#'
#' @concept other_seurat_plotting
#'
#' @examples
#' \dontrun{
#' Plot_Density_Custom(seurat_object = object, features = c("Cx3cr1", "Fcrls"))
#' }
#'

Plot_Density_Custom <- function(
  seurat_object,
  features,
  joint = FALSE,
  viridis_palette = "magma",
  custom_palette = NULL,
  pt.size = 1,
  reduction = NULL,
  combine = TRUE,
  ...
) {
  # Check Nebulosa installed
  Nebulosa_check <- PackageCheck("Nebulosa", error = FALSE)
  if (!Nebulosa_check[1]) {
    stop(
      "Please install the Nebulosa package to use 'Plot_Density_Custom'",
      "\nThis can be accomplished with the following commands: ",
      "\n----------------------------------------",
      "\ninstall.packages('BiocManager')",
      "\nBiocManager::install('Nebulosa')",
      "\n----------------------------------------",
      call. = FALSE
    )
  }

  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # check palettes
  if (!is.null(x = custom_palette) && viridis_palette != "magma") {
    cli_abort(message = c("Non-default values provided to both viridis_palette & custom_palette.",
                          "i" = "Please chose one non-default value.")
    )
  }

  # Extract default reduction
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)

  # Create plot list
  plot_list <- Nebulosa::plot_density(object = seurat_object,
                                      features = features,
                                      reduction = reduction,
                                      size = pt.size,
                                      combine = combine,
                                      pal = viridis_palette,
                                      joint = joint,
                                      ...)

  if (!is.null(x = custom_palette)) {
    suppressMessages(plot_list <- plot_list & scale_color_gradientn(colors = custom_palette))
    return(plot_list)
  }

  return(plot_list)
}


#' Nebulosa Joint Density Plot
#'
#' Return only the joint density plot from Nebulosa plot_density function.  Requires Nebulosa package from Bioconductor.
#'
#' @param seurat_object Seurat object name.
#' @param features Features to plot.
#' @param viridis_palette default viridis palette to use (must be one of: "viridis", "magma", "cividis",
#' "inferno", "plasma").  Default is "magma".
#' @param custom_palette non-default color palette to be used in place of default viridis options.
#' @param pt.size Adjust point size for plotting.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default).
#' @param ... Extra parameters passed to \code{\link[Nebulosa]{plot_density}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
# #' @importFrom Nebulosa plot_density
#' @importFrom SeuratObject DefaultDimReduc PackageCheck
#'
#' @export
#'
#' @concept other_seurat_plotting
#'
#' @examples
#' \dontrun{
#' Plot_Density_Joint_Only(seurat_object = object, features = c("Cx3cr1", "Fcrls"))
#' }
#'

Plot_Density_Joint_Only <- function(
  seurat_object,
  features,
  viridis_palette = "magma",
  custom_palette = NULL,
  pt.size = 1,
  reduction = NULL,
  ...
) {
  # Check Nebulosa installed
  Nebulosa_check <- PackageCheck("Nebulosa", error = FALSE)
  if (!Nebulosa_check[1]) {
    stop(
      "Please install the Nebulosa package to use 'Plot_Density_Joint_Only'",
      "\nThis can be accomplished with the following commands: ",
      "\n----------------------------------------",
      "\ninstall.packages('BiocManager')",
      "\nBiocManager::install('Nebulosa')",
      "\n----------------------------------------",
      call. = FALSE
    )
  }

  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check features
  if (length(x = features) < 2) {
    cli_abort(message = c("Less than 2 features provided.",
                          "i" = "Nebulosa joint density plots require two or more features.")
    )
  }

  # check palettes
  if (!is.null(x = custom_palette) && viridis_palette != "magma") {
    cli_abort(message = c("Non-default values provided to both viridis_palette & custom_palette.",
                          "i" = "Please chose one non-default value.")
    )
  }

  # Extract default reduction
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)

  # Create plot list
  plot_list <- Nebulosa::plot_density(object = seurat_object,
                                      features = features,
                                      reduction = reduction,
                                      size = pt.size,
                                      joint = TRUE,
                                      combine = FALSE,
                                      pal = viridis_palette,
                                      ...)

  # return the joint plot only
  plot <- plot_list[[length(x = plot_list)]]

  if (!is.null(x = custom_palette)) {
    suppressMessages(plot <- plot + scale_color_gradientn(colors = custom_palette))
    return(plot)
  }

  return(plot)
}
