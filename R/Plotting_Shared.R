#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### SHARED SEURAT & LIGER PLOTTING ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Factor Correlation Plot
#'
#' Plot positive correlations between gene loadings across `W` factor matrix in LIGER or
#' feature loadings in reduction slot of Seurat object.
#'
#' @param object LIGER or Seurat object.
#' @param reduction Seurat ONLY; name of dimensionality reduction containing NMF loadings.
#' @param colors_use Color palette to use for correlation values.
#' Default is `RColorBrewer::RdBu` if `positive_only = FALSE`.
#' If `positive_only = TRUE` the default is `viridis`.
#' Users can also supply vector of 3 colors (low, mid, high).
#' @param label logical, whether to add correlation values to plot result.
#' @param label_threshold threshold for adding correlation values if `label = TRUE`.  Default
#' is 0.5.
#' @param label_size size of correlation labels
#' @param plot_title Plot title.
#' @param plot_type Controls plotting full matrix, or just the upper or lower triangles.
#' Accepted values are: "full" (default), "upper", or "lower".
#' @param positive_only logical, whether to limit the plotted values to only positive
#' correlations (negative values set to 0); default is FALSE.
#' @param x_lab_rotate logical, whether to rotate the axes labels on the x-axis.  Default is TRUE.
#' @param cluster logical, whether to cluster the plot using `hclust` (default TRUE).  If FALSE
#' factors are listed in numerical order.
#' @param cluster_rect logical, whether to add rectangles around the clustered areas on plot,
#' default is FALSE. Uses `cutree` to create groups.
#' @param cluster_rect_num number of rectangles to add to the plot, default NULL.
#' Value is provided to `k` in `cutree`.
#' @param cluster_rect_col color to use for rectangles, default MULL (will set color automatically).
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom cowplot theme_cowplot
#' @importFrom dplyr arrange any_of
#' @importFrom magrittr "%>%"
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr drop_na pivot_longer
#'
#' @export
#'
#' @concept liger_plotting
#'
#' @examples
#' \dontrun{
#' Factor_Cor_Plot(object = obj)
#'}
#'

Factor_Cor_Plot <- function(
    object,
    reduction = NULL,
    colors_use = NULL,
    label = FALSE,
    label_threshold = 0.5,
    label_size = 5,
    plot_title = NULL,
    plot_type = "full",
    positive_only = FALSE,
    x_lab_rotate = TRUE,
    cluster = TRUE,
    cluster_rect = FALSE,
    cluster_rect_num = NULL,
    cluster_rect_col = NULL
) {
  # check plot type
  if (!plot_type %in% c("full", "lower", "upper")) {
    cli_abort(message = "{.code plot_type} must be one of {.field {glue_collapse_scCustom(input_string = c('full', 'lower,', 'upper'), and = FALSE)}}")
  }

  # check upper/lower and cluster parameters
  if (isTRUE(x = cluster) && plot_type %in% c("lower", "upper")) {
    cli::cli_abort(message = c("When setting {.code plot_type} to {.val upper} or {.val lower} the {.code cluster} parameter must be set to {.field FALSE}."))
  }

  # get data Seurat
  if (inherits(x = object, what = "Seurat")) {
    # check a reduction provided
    if (is.null(x = reduction)) {
      cli_abort(message = "Must supply name of reduction to {.code reduction} parameter when plotting from Seurat object")
    }

    # check reduction present
    if (!reduction %in% Reductions(object = object)) {
      cli_abort(message = "Provided reduction: {.field {reduction}} was not found in Seurat Object.")
    }

    # get cor data
    cor_mat <- Find_Factor_Cor(object = object, reduction = reduction)
  }

  # get data liger
  if (inherits(x = object, what = "liger")) {
    cor_mat <- Find_Factor_Cor(object = object)
  }

  # filter matrix by plot type
  if (plot_type == "upper") {
    plot_df <- upper_diag_cor_mat(cor_mat = cor_mat)
  }

  if (plot_type == "lower") {
    plot_df <- lower_diag_cor_mat(cor_mat = cor_mat)
  }

  if (plot_type == "full") {
    plot_df <- cor_mat
  }

  if (isTRUE(x = cluster)) {
    dist_mat <- stats::as.dist((1 - plot_df) / 2)
    hclust_res <- stats::hclust(dist_mat, method = "complete")

    plot_df <- plot_df[hclust_res$order, hclust_res$order]
  }

  # Reshape for plotting
  plot_df <- data.frame(plot_df) %>%
    rownames_to_column("rowname") %>%
    pivot_longer(cols = !any_of("rowname"), names_to = "Var", values_to = "corr") %>%
    drop_na()

  plot_df$rowname <- factor(plot_df$rowname, levels = rev(unique(plot_df$rowname)))

  if (isTRUE(x = cluster)) {
    plot_df$Var <- factor(plot_df$Var, levels = unique(plot_df$Var))
  }

  if (isTRUE(x = label)) {
    plot_df$label <- ifelse(plot_df$corr >= label_threshold, round(plot_df$corr, 2), NA)
    plot_df$label <- ifelse(plot_df$label == 1, NA, round(plot_df$label, 2))
  }

  factor_names <- levels(plot_df$rowname)

  # plot
  if (isTRUE(x = positive_only)) {
    colors_use <- colors_use %||% viridis_light_high
    cluster_rect_col <- cluster_rect_col %||% "white"

    plot <- ggplot(data = plot_df, mapping = aes(x = .data[["Var"]], y = .data[["rowname"]], fill = .data[["corr"]])) +
      theme_cowplot() +
      geom_tile() +
      scale_y_discrete(limits = factor_names, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_gradientn(colours = colors_use, limits = c(0,1), na.value = colors_use[1]) +
      xlab("") +
      ylab("")
  } else {
    colors_use <- colors_use %||% paletteer::paletteer_d("RColorBrewer::RdBu")
    cluster_rect_col <- cluster_rect_col %||% "black"

    plot <- ggplot(data = plot_df, mapping = aes(x = .data[["Var"]], y = .data[["rowname"]], fill = .data[["corr"]])) +
      theme_cowplot() +
      geom_tile() +
      scale_y_discrete(limits = factor_names, expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_gradientn(colours = colors_use, limits = c(-1,1), na.value = colors_use[1]) +
      xlab("") +
      ylab("")
  }

  # modify plot
  if (isTRUE(x = label)) {
    plot <- suppressMessages(plot + geom_text(aes(label=label), size = label_size))
  }

  if (!is.null(x = plot_title)) {
    plot <- plot + ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5))
  }

  if (isTRUE(x = x_lab_rotate)) {
    plot <- plot + RotatedAxis()
  }

  if (isTRUE(x = cluster_rect)) {
    rect_list <- create_factor_hclust_rect(cor_mat = cor_mat, num_rect = cluster_rect_num, num_factors = length(x = factor_names))

    plot <- plot + annotate(geom = "rect", xmin = rect_list[[1]][,1], xmax = rect_list[[1]][,2], ymin = rect_list[[2]][,1], ymax = rect_list[[2]][,2], fill = NA, color = cluster_rect_col)
  }

  # return plot
  return(plot)
}
