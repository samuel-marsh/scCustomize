#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### LIGER INTERNAL UTILS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Extract default dimensionality reduction
#'
#' Extract name of the default dimensionlity reduction for liger object.
#'
#' @param liger_object LIGER object name.
#'
#' @return name of default dimensionality reduction
#'
#' @import cli
#'
#' @noRd
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # return dimensionality reduction name
#' dim_reduc_name <- Default_DimReduc_LIGER(liger_object = obj)
#' }
#'

Default_DimReduc_LIGER <- function(
  liger_object
) {
  if (length(x = liger_object@dimReds) > 0) {
    default_reduc <- liger_object@uns$defaultDimRed

    return(default_reduc)
  } else {
    cli_abort(message = "No dimensionality reduction present.")
  }
}


#' Extract default clustering
#'
#' Extract name of the default clustering
#'
#' @param liger_object LIGER object name.
#'
#' @return name of default clustering
#'
#' @import cli
#'
#' @noRd
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # return dimensionality reduction name
#' dim_reduc_name <- LIGER_Default_Cluster_Name(liger_object = obj)
#' }
#'

LIGER_Default_Cluster_Name <- function(
  liger_object
) {
  if (is.null(x = rliger::defaultCluster(x = liger_object))) {
    cli_abort(message = "No default cell identity/cluster present in object.")
  } else {
    default_cluster_name <- liger_object@uns$defaultCluster
    return(default_cluster_name)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### LIGER PLOTTING UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' DimPlot LIGER Version
#'
#' Standard and modified version of LIGER's plotByDatasetAndCluster
#'
#' @param object Name of LIGER object.  Need to perform clustering before calling this function.
#' @param clusters Another clustering to use for coloring second plot (must have same names as
#' clusters slot) (default NULL).
#' @param shuffle Randomly shuffle points so that points from same dataset are not plotted one after
#' the other (default TRUE).
#' @param shuffle_seed Random seed for reproducibility of point shuffling (default 1).
#' @param redorder.idents logical whether to reorder the datasets from default order before plotting (default FALSE).
#' @param new.order new dataset factor order for plotting.  must set reorder.idents = TRUE.
#' @param group.by meta data varibale to group plots by
#' @param split.by meta data variable to splot plots by
#'
#' @return A data.frame with information for plotting
#'
#' @importFrom utils packageVersion
#'
#' @references This function is encompasses the first part of the LIGER function plotByDatasetAndCluster.
#' However, this function is modified to allow plotting other meta data variables.  In this case the function
#' just returns the data.frame needed for plotting rather than plots themselves.
#' \url{https://github.com/welch-lab/liger}. (License: GPL-3).
#'
#' @noRd
#'
#' @concept liger_plotting_util
#'

Generate_Plotting_df_LIGER <- function(
  object,
  clusters = NULL,
  shuffle = TRUE,
  shuffle_seed = 1,
  reorder.idents = FALSE,
  new.order = NULL,
  group.by = "dataset",
  split.by = NULL
) {
  tsne_df <- data.frame(object@tsne.coords)
  colnames(x = tsne_df) <- c("tsne1", "tsne2")
  tsne_df[[group.by]] <- object@cell.data[[group.by]]
  if (!is.null(x = split.by)) {
    tsne_df[[split.by]] <- object@cell.data[[split.by]]
  }

  if (isTRUE(x = reorder.idents)) {
    tsne_df[[group.by]]  <- factor(x = tsne_df[[group.by]], levels = new.order)
  }
  c_names <- names(x = object@clusters)
  if (is.null(x = clusters)) {
    # if clusters have not been set yet
    if (length(x = object@clusters) == 0) {
      clusters <- rep(1, nrow(x = object@tsne.coords))
      names(x = clusters) <- c_names <- rownames(x = object@tsne.coords)
    } else {
      clusters <- object@clusters
      c_names <- names(x = object@clusters)
    }
  }
  tsne_df[["Cluster"]] <- clusters[c_names]

  if (isTRUE(x = shuffle)) {
    set.seed(shuffle_seed)
    idx <- sample(x = 1:nrow(x = tsne_df))
    tsne_df <- tsne_df[idx, ]
  }
  return(tsne_df)
}


Generate_Plotting_df_LIGER2 <- function(
  object,
  reduction = NULL,
  clusters = NULL,
  shuffle = TRUE,
  shuffle_seed = 1,
  reorder.idents = FALSE,
  new.order = NULL,
  group.by = "dataset",
  split.by = NULL
) {
  # Set reduction if null
  if (!is.null(x = reduction)) {
    Embeddings(object = object, reduction = reduction, check_only = TRUE)
  } else {
    reduction <- reduction %||% Default_DimReduc_LIGER(liger_object = object)
  }

  reduc_df <- data.frame(Embeddings(object = object, reduction = reduction))
  reduc_df[[group.by]] <- object@cellMeta[[group.by]]
  if (!is.null(x = split.by)) {
    reduc_df[[split.by]] <- object@cellMeta[[split.by]]
  }

  if (isTRUE(x = reorder.idents)) {
    reduc_df[[group.by]]  <- factor(x = reduc_df[[group.by]], levels = new.order)
  }
  cluster_col <- LIGER_Default_Cluster_Name(liger_object = object)
  if (is.null(x = clusters)) {
    # if clusters have not been set yet
    if (length(x = object@cellMeta[[cluster_col]]) == 0) {
      clusters <- rep(1, nrow(x = reduc_df))
      names(x = clusters) <- rownames(x = reduc_df)
    } else {
      clusters <- object@cellMeta[[cluster_col]]
    }
  }
  reduc_df[["Cluster"]] <- clusters

  if (isTRUE(x = shuffle)) {
    set.seed(shuffle_seed)
    idx <- sample(x = 1:nrow(reduc_df))
    reduc_df <- reduc_df[idx, ]
  }
  return(reduc_df)
}



#' LIGER plot by cluster.
#'
#' Modified version of LIGER's plotByDatasetAndCluster just for plotting clusters.
#'
#' @param liger_object Name of LIGER object.  Need to perform clustering before calling this function.
#' @param colors_use colors to use for plotting by cluster.  By default if number of levels plotted is
#' less than or equal to 36 it will use "polychrome" and if greater than 36 will use "varibow" with
#' shuffle = TRUE both from \code{\link{DiscretePalette_scCustomize}}.
#' @param group.by Variable to be plotted.  If `NULL` will plot clusters from `liger@clusters` slot.
#' If `combination = TRUE` will plot both clusters and meta data variable.
#' @param split.by meta data variable to split plots by (i.e. "dataset").
#' @param title plot title.
#' @param pt_size Adjust point size for plotting.
#' @param reduction_label What to label the x and y axes of resulting plots.  LIGER does not store
#' name of technique and therefore needs to be set manually.  Default is "UMAP".
#' @param num_columns Number of columns to plot by if `split.by` is not NULL.
#' @param shuffle logical. Whether to randomly shuffle the order of points. This can be useful for
#' crowded plots if points of interest are being buried. (Default is TRUE).
#' @param shuffle_seed Sets the seed if randomly shuffling the order of points.
#' @param legend.size what to set legend size to.
#' @param label logical.  Whether or not to label the clusters.  Default is TRUE.
#' @param label_size size of cluster labels.
#' @param label_repel logical.  Whether to repel cluster labels from each other if plotting by
#' cluster (if `group.by = NULL` or `group.by = "cluster`).  Default is FALSE.
#' @param label_box logical.  Whether to put a box around the label text (uses `geom_text` vs `geom_label`).
#' Default is FALSE.
#' @param label_color Color to use for cluster labels.  Default is "black".
#' @param redorder.idents logical. should the idents plotted by reordered.  Default is FALSE.
#' @param new.order What should the new ident order be if `reorder.idents = TRUE`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#'
#' @return A ggplot/patchwork object
#'
#' @import ggplot2
#' @importFrom cowplot theme_cowplot
#' @importFrom dplyr summarize
#' @importFrom ggrepel geom_text_repel geom_label_repel
#' @importFrom patchwork wrap_plots
#' @importFrom scattermore geom_scattermore
#' @importFrom stats median
#' @importFrom utils packageVersion
#'
#' @references This function is encompasses part of the LIGER function plotByDatasetAndCluster.
#' However, this function is modified to just return cluster plots based on `Generate_Plotting_df_LIGER`.
#' \url{https://github.com/welch-lab/liger}. (Licence: GPL-3).
#'
#' @noRd
#'
#' @concept liger_plotting_util
#'

Plot_By_Cluster_LIGER <- function(
  liger_object,
  colors_use = NULL,
  group.by = "dataset",
  split.by = NULL,
  title = NULL,
  pt_size = NULL,
  reduction_label = "UMAP",
  num_columns = NULL,
  shuffle = TRUE,
  shuffle_seed = 1,
  legend.size = 5,
  label = TRUE,
  label_size = NA,
  label_repel = FALSE,
  label_box = FALSE,
  label_color = "black",
  reorder.idents = FALSE,
  new.order = NULL,
  raster = NULL,
  raster.dpi = c(512, 512),
  ggplot_default_colors = FALSE,
  color_seed = 123
) {
  # Check dimreduc present
  if (length(x = liger_object@tsne.coords) == 0) {
    cli_abort(message = "No dimensionality reduction coordinates found.")
  }

  # Create plotting data.frame
  tsne_df <- Generate_Plotting_df_LIGER(object = liger_object, group.by = group.by, split.by = split.by, reorder.idents = reorder.idents, shuffle = shuffle, shuffle_seed = shuffle_seed)

  if (!is.null(x = split.by)) {
    list_of_splits <- unique(x = tsne_df[[split.by]])
  }

  # Get length of meta data feature
  if (!is.null(x = split.by) && !is.null(x = num_columns)) {
    split.by_length <- length(x = list_of_splits)

    # Calculate number of rows for selected number of columns
    num_rows <- ceiling(x = split.by_length / num_columns)

    # Check column and row compatibility
    if (num_columns > split.by_length) {
      cli_abort(message = c("The number of columns specified is greater than the number of meta data variables.",
                            "*" = "{.field {split.by}} only contains: {.field {split.by_length}} variables.",
                            "i" = "Please adjust {.code num_columns} to be less than or equal to: {.field {split.by_length}}.")
      )
    }
  }

  centers <- tsne_df %>%
    group_by(.data[["Cluster"]]) %>%
    summarize(tsne1 = median(x = .data[["tsne1"]]), tsne2 = median(x = .data[["tsne2"]]))

  cluster_length <- length(x = unique(x = liger_object@clusters))

  # set default plot colors
  if (is.null(x = colors_use)) {
    colors_use <- scCustomize_Palette(num_groups = cluster_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
  }

  # Create accurate axis labels
  x_axis_label <- paste0(reduction_label, "_1")
  y_axis_label <- paste0(reduction_label, "_2")

  # plot
  if (isTRUE(x = raster)) {
    if (!is.null(x = split.by)) {
      p2 <- lapply(1:length(x = list_of_splits), function(x) {
        p2 <- ggplot(subset(tsne_df, tsne_df[[split.by]] %in% list_of_splits[x]), aes(x = .data[["tsne1"]], y = .data[["tsne2"]], color = .data[["Cluster"]])) +
          theme_cowplot() +
          geom_scattermore(pointsize = pt_size, pixels = raster.dpi) +
          guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(list_of_splits[x]) +
          scale_color_manual(values = colors_use) +
          theme(legend.position = "right",
                axis.text = element_text(size = rel(0.95)),
                plot.title = element_text(hjust = 0.5)) +
          guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
          xlab(x_axis_label) +
          ylab(y_axis_label)

        if (isTRUE(x = label_box)) {
          geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes(label = .data[["Cluster"]], fill = .data[["Cluster"]]), size = label_size,
            show.legend = FALSE, color = label_color
          ) + scale_fill_manual(values = colors_use)
        } else if (isTRUE(x = label)) {
          geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes(label = .data[["Cluster"]]), size = label_size, color = label_color,
            show.legend = FALSE
          )
        } else {
          p2 <- p2
        }
      })
    } else {
      p2 <- ggplot(tsne_df, aes(x = .data[["tsne1"]], y = .data[["tsne2"]], color = .data[["Cluster"]])) +
        theme_cowplot() +
        geom_scattermore(pointsize = pt_size, pixels = raster.dpi) +
        guides(color = guide_legend(override.aes = list(size = legend.size))) +
        scale_color_manual(values = colors_use) +
        theme(legend.position = "right",
              axis.text = element_text(size = rel(0.95)),
              plot.title = element_text(hjust = 0.5)) +
        guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
        xlab(x_axis_label) +
        ylab(y_axis_label)

      if (isTRUE(x = label_box)) {
        geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes(label = .data[["Cluster"]], fill = .data[["Cluster"]]), size = label_size,
          show.legend = FALSE, color = label_color
        ) + scale_fill_manual(values = colors_use)
      } else if (isTRUE(x = label)) {
        geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes(label = .data[["Cluster"]]), size = label_size, color = label_color,
          show.legend = FALSE
        )
      } else {
        p2 <- p2
      }

    }
  } else {
    if (!is.null(x = split.by)) {
      p2 <- lapply(1:length(x = list_of_splits), function(x) {
        p2 <- ggplot(subset(tsne_df, tsne_df[[split.by]] %in% list_of_splits[x]), aes(x = .data[["tsne1"]], y = .data[["tsne2"]], color = .data[["Cluster"]])) +
          theme_cowplot() +
          geom_point(size = pt_size) +
          guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(list_of_splits[x]) +
          scale_color_manual(values = colors_use) +
          theme(legend.position = "right",
                axis.text = element_text(size = rel(0.95)),
                plot.title = element_text(hjust = 0.5)) +
          guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
          xlab(x_axis_label) +
          ylab(y_axis_label)

        if (isTRUE(x = label_box)) {
          geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes(label = .data[["Cluster"]], fill = .data[["Cluster"]]), size = label_size,
            show.legend = FALSE, color = label_color
          ) + scale_fill_manual(values = colors_use)
        } else if (isTRUE(x = label)) {
          geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes(label = .data[["Cluster"]]), size = label_size, color = label_color,
            show.legend = FALSE
          )
        } else {
          p2 <- p2
        }
      })
    } else {
      p2 <- ggplot(tsne_df, aes(x = .data[["tsne1"]], y = .data[["tsne2"]], color = .data[["Cluster"]])) +
        theme_cowplot() +
        geom_point(size = pt_size) +
        guides(color = guide_legend(override.aes = list(size = legend.size))) +
        scale_color_manual(values = colors_use) +
        theme(legend.position = "right",
              axis.text = element_text(size = rel(0.95)),
              plot.title = element_text(hjust = 0.5)) +
        guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
        xlab(x_axis_label) +
        ylab(y_axis_label)

      if (isTRUE(x = label_box)) {
        geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes(label = .data[["Cluster"]], fill = .data[["Cluster"]]), size = label_size,
          show.legend = FALSE, color = label_color
        ) + scale_fill_manual(values = colors_use)
      } else if (isTRUE(x = label)) {
        geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes(label = .data[["Cluster"]]), size = label_size, color = label_color,
          show.legend = FALSE
        )
      } else {
        p2 <- p2
      }
    }
  }
  if (!is.null(x = split.by) && !is.null(x = num_columns)) {
    p2 <- wrap_plots(p2) + plot_layout(nrow = num_rows, ncol = num_columns, guides = "collect")
    return(p2)
  }
  if (!is.null(x = split.by) && is.null(x = num_columns)) {
    p2 <- wrap_plots(p2) + plot_layout(guides = "collect")
    return(p2)
  } else {
    return(p2)
  }
}


Plot_By_Cluster_LIGER2 <- function(
    liger_object,
    colors_use = NULL,
    group.by = "dataset",
    split.by = NULL,
    title = NULL,
    pt_size = NULL,
    reduction = NULL,
    num_columns = NULL,
    shuffle = TRUE,
    shuffle_seed = 1,
    legend.size = 5,
    label = TRUE,
    label_size = NA,
    label_repel = FALSE,
    label_box = FALSE,
    label_color = "black",
    reorder.idents = FALSE,
    new.order = NULL,
    raster = NULL,
    raster.dpi = c(512, 512),
    ggplot_default_colors = FALSE,
    color_seed = 123
) {
  # Set reduction
  reduction <- reduction %||% Default_DimReduc_LIGER(liger_object = liger_object)

  # Create plotting data.frame
  reduc_df <- Generate_Plotting_df_LIGER2(object = liger_object, group.by = group.by, split.by = split.by, reorder.idents = reorder.idents, shuffle = shuffle, shuffle_seed = shuffle_seed, reduction = reduction)

  if (!is.null(x = split.by)) {
    list_of_splits <- unique(x = reduc_df[[split.by]])
  }

  # Get length of meta data feature
  if (!is.null(x = split.by) && !is.null(x = num_columns)) {
    split.by_length <- length(x = list_of_splits)

    # Calculate number of rows for selected number of columns
    num_rows <- ceiling(x = split.by_length / num_columns)

    # Check column and row compatibility
    if (num_columns > split.by_length) {
      cli_abort(message = c("The number of columns specified is greater than the number of meta data variables.",
                            "*" = "{.field {split.by}} only contains: {.field {split.by_length}} variables.",
                            "i" = "Please adjust {.code num_columns} to be less than or equal to: {.field {split.by_length}}.")
      )
    }
  }

  # Create accurate axis labels
  x_axis_label <- names(x = reduc_df)[1]
  y_axis_label <- names(x = reduc_df)[2]

  centers <- reduc_df %>%
    group_by(.data[["Cluster"]]) %>%
    summarize(dr1 = median(x = .data[[x_axis_label]]), dr2 = median(x = .data[[y_axis_label]]))

  colnames(x = centers) <- c("Cluster", x_axis_label, y_axis_label)

  cluster_length <- length(x = unique(x = rliger::defaultCluster(x = liger_object)))

  if (is.null(x = colors_use)) {
    colors_use <- scCustomize_Palette(num_groups = cluster_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
  }

  # plot
  if (isTRUE(x = raster)) {
    if (!is.null(x = split.by)) {
      p2 <- lapply(1:length(x = list_of_splits), function(x) {
        p2 <- ggplot(data = subset(reduc_df, reduc_df[[split.by]] %in% list_of_splits[x]), aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = .data[["Cluster"]])) +
          theme_cowplot() +
          geom_scattermore(pointsize = pt_size, pixels = raster.dpi) +
          guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(list_of_splits[x]) +
          scale_color_manual(values = colors_use) +
          theme(legend.position = "right",
                axis.text = element_text(size = rel(0.95)),
                plot.title = element_text(hjust = 0.5)) +
          guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
          xlab(x_axis_label) +
          ylab(y_axis_label)

        if (isTRUE(x = label_box)) {
          geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes(label = .data[["Cluster"]], fill = .data[["Cluster"]]), size = label_size,
            show.legend = FALSE, color = label_color
          ) + scale_fill_manual(values = colors_use)
        } else if (isTRUE(x = label)) {
          geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes(label = .data[["Cluster"]]), size = label_size, color = label_color,
            show.legend = FALSE
          )
        } else {
          p2 <- p2
        }
      })
    } else {
      p2 <- ggplot(data = reduc_df, aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = .data[["Cluster"]])) +
        theme_cowplot() +
        geom_scattermore(pointsize = pt_size, pixels = raster.dpi) +
        guides(color = guide_legend(override.aes = list(size = legend.size))) +
        scale_color_manual(values = colors_use) +
        theme(legend.position = "right",
              axis.text = element_text(size = rel(0.95)),
              plot.title = element_text(hjust = 0.5)) +
        guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
        xlab(x_axis_label) +
        ylab(y_axis_label)

      if (isTRUE(x = label_box)) {
        geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes(label = .data[["Cluster"]], fill = .data[["Cluster"]]), size = label_size,
          show.legend = FALSE, color = label_color
        ) + scale_fill_manual(values = colors_use)
      } else if (isTRUE(x = label)) {
        geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes(label = .data[["Cluster"]]), size = label_size, color = label_color,
          show.legend = FALSE
        )
      } else {
        p2 <- p2
      }
    }
  } else {
    if (!is.null(x = split.by)) {
      p2 <- lapply(1:length(x = list_of_splits), function(x) {
        p2 <- ggplot(data = subset(reduc_df, reduc_df[[split.by]] %in% list_of_splits[x]), aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = .data[["Cluster"]])) +
          theme_cowplot() +
          geom_point(size = pt_size) +
          guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(list_of_splits[x]) +
          scale_color_manual(values = colors_use) +
          theme(legend.position = "right",
                axis.text = element_text(size = rel(0.95)),
                plot.title = element_text(hjust = 0.5)) +
          guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
          xlab(x_axis_label) +
          ylab(y_axis_label)

        if (isTRUE(x = label_box)) {
          geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes(label = .data[["Cluster"]], fill = .data[["Cluster"]]), size = label_size,
            show.legend = FALSE, color = label_color
          ) + scale_fill_manual(values = colors_use)
        } else if (isTRUE(x = label)) {
          geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
          p2 <- p2 + geom.use(
            data = centers,
            mapping = aes(label = .data[["Cluster"]]), size = label_size, color = label_color,
            show.legend = FALSE
          )
        } else {
          p2 <- p2
        }
      })
    } else {
      p2 <- ggplot(data = reduc_df, aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = .data[["Cluster"]])) +
        theme_cowplot() +
        geom_point(size = pt_size) +
        guides(color = guide_legend(override.aes = list(size = legend.size))) +
        scale_color_manual(values = colors_use) +
        theme(legend.position = "right",
              axis.text = element_text(size = rel(0.95)),
              plot.title = element_text(hjust = 0.5)) +
        guides(col = guide_legend(title = "", override.aes = list(size = 4)))

      if (isTRUE(x = label_box)) {
        geom.use <- ifelse(test = label_repel, yes = geom_label_repel, no = geom_label)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes(label = .data[["Cluster"]], fill = .data[["Cluster"]]), size = label_size,
          show.legend = FALSE, color = label_color
        ) + scale_fill_manual(values = colors_use)
      } else if (isTRUE(x = label)) {
        geom.use <- ifelse(test = label_repel, yes = geom_text_repel, no = geom_text)
        p2 <- p2 + geom.use(
          data = centers,
          mapping = aes(label = .data[["Cluster"]]), size = label_size, color = label_color,
          show.legend = FALSE
        )
      } else {
        p2 <- p2
      }
    }
  }
  if (!is.null(x = split.by) && !is.null(x = num_columns)) {
    p2 <- wrap_plots(p2) + plot_layout(nrow = num_rows, ncol = num_columns, guides = "collect")
    return(p2)
  }
  if (!is.null(x = split.by) && is.null(x = num_columns)) {
    p2 <- wrap_plots(p2) + plot_layout(guides = "collect")
    return(p2)
  } else {
    return(p2)
  }
}


#' LIGER plot by meta variables.
#'
#' Modified version of LIGER's plotByDatasetAndCluster just for plotting meta variables.
#'
#' @param liger_object Name of LIGER object.  Need to perform clustering before calling this function.
#' @param colors_use colors to use for plotting by cluster.  By default if number of levels plotted is
#' less than or equal to 36 it will use "polychrome" and if greater than 36 will use "varibow" with
#' shuffle = TRUE both from \code{\link{DiscretePalette_scCustomize}}.
#' @param group.by Variable to be plotted.  If `NULL` will plot clusters from `liger@clusters` slot.
#' If `combination = TRUE` will plot both clusters and meta data variable.
#' @param split.by meta data variable to split plots by (i.e. "dataset").
#' @param title plot title.
#' @param pt_size Adjust point size for plotting.
#' @param reduction_label What to label the x and y axes of resulting plots.  LIGER does not store name
#' of technique and therefore needs to be set manually.  Default is "UMAP".
#' @param num_columns Number of columns to plot by if `split.by` is not NULL.
#' @param shuffle logical. Whether to randomly shuffle the order of points. This can be useful for
#' crowded plots if points of interest are being buried. (Default is TRUE).
#' @param shuffle_seed Sets the seed if randomly shuffling the order of points.
#' @param legend.size what to set legend size to.
#' @param redorder.idents logical. should the idents plotted by reordered.  Default is FALSE.
#' @param new.order What should the new ident order be if `reorder.idents = TRUE`.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#'
#' @return A ggplot/patchwork object
#'
#' @import ggplot2
#' @importFrom cowplot theme_cowplot
#' @importFrom patchwork wrap_plots
#' @importFrom rlang sym "!!"
#' @importFrom scattermore geom_scattermore
#' @importFrom utils packageVersion
#'
#' @references This function is encompasses part of the LIGER function plotByDatasetAndCluster.
#' However, this function is modified to just return cluster plots based on `Generate_Plotting_df_LIGER`.
#' \url{https://github.com/welch-lab/liger}. (Licence: GPL-3).
#'
#' @noRd
#'
#' @concept liger_plotting_util
#'

Plot_By_Meta_LIGER <- function(
    liger_object,
    colors_use = NULL,
    group.by = "dataset",
    split.by = NULL,
    title = NULL,
    pt_size = NULL,
    reduction_label = "UMAP",
    num_columns = NULL,
    shuffle = TRUE,
    shuffle_seed = 1,
    legend.size = 3,
    reorder.idents = FALSE,
    new.order = NULL,
    raster = NULL,
    raster.dpi = c(512, 512),
    ggplot_default_colors = FALSE,
    color_seed = 123
) {
  # Check dimreduc present
  if (length(x = liger_object@tsne.coords) == 0) {
    cli_abort(message = "No dimensionality reduction coordinates found.")
  }

  tsne_df <- Generate_Plotting_df_LIGER(object = liger_object, group.by = group.by, split.by = split.by, reorder.idents = reorder.idents, shuffle = shuffle, shuffle_seed = shuffle_seed)

  if (!is.null(x = split.by)) {
    list_of_splits <- unique(x = tsne_df[[split.by]])
  }

  # Get length of meta data feature
  if (!is.null(x = split.by) && !is.null(x = num_columns)) {
    split.by_length <- length(x = list_of_splits)

    # Calculate number of rows for selected number of columns
    num_rows <- ceiling(x = split.by_length / num_columns)

    # Check column and row compatibility
    if (num_columns > split.by_length) {
      cli_abort(message = c("The number of columns specified is greater than the number of meta data variables.",
                            "*" = "{.field {split.by}} only contains: {.field {split.by_length}} variables.",
                            "i" = "Please adjust {.code num_columns} to be less than or equal to: {.field {split.by_length}}.")
      )
    }
  }

  meta_length <- length(x = unique(x = liger_object@cell.data[[group.by]]))

  # set default plot colors
  if (is.null(x = colors_use)) {
    colors_use <- scCustomize_Palette(num_groups = meta_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
  }

  # Create accurate axis labels
  x_axis_label <- paste0(reduction_label, "_1")
  y_axis_label <- paste0(reduction_label, "_2")

  group.by <- sym(x = group.by)

  if (isTRUE(x = raster)) {
    if (!is.null(x = split.by)) {
      p1 <- lapply(1:length(x = list_of_splits), function(x) {
        ggplot(subset(tsne_df, tsne_df[[split.by]] %in% list_of_splits[x]), aes(x = .data[["tsne1"]], y = .data[["tsne2"]], color = !!group.by)) +
          theme_cowplot() +
          geom_scattermore(pointsize = pt_size, pixels = raster.dpi) +
          guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(list_of_splits[x]) +
          scale_color_manual(values = colors_use) +
          theme(legend.position = "right",
                axis.text = element_text(size = rel(0.95)),
                plot.title = element_text(hjust = 0.5)) +
          guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
          xlab(x_axis_label) +
          ylab(y_axis_label)
      })
    } else {
      p1 <- ggplot(tsne_df, aes(x = .data[["tsne1"]], y = .data[["tsne2"]], color = !!group.by)) +
        theme_cowplot() +
        geom_scattermore(pointsize = pt_size, pixels = raster.dpi) +
        guides(color = guide_legend(override.aes = list(size = legend.size))) +
        scale_color_manual(values = colors_use) +
        theme(legend.position = "right",
              axis.text = element_text(size = rel(0.95)),
              plot.title = element_text(hjust = 0.5)) +
        guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
        xlab(x_axis_label) +
        ylab(y_axis_label)

    }
  } else {
    if (!is.null(x = split.by)) {
      p1 <- lapply(1:length(x = list_of_splits), function(x) {
        ggplot(subset(tsne_df, tsne_df[[split.by]] %in% list_of_splits[x]), aes(x = .data[["tsne1"]], y = .data[["tsne2"]], color = !!group.by)) +
          theme_cowplot() +
          geom_point(size = pt_size) +
          guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(list_of_splits[x]) +
          scale_color_manual(values = colors_use) +
          theme(legend.position = "right",
                axis.text = element_text(size = rel(0.95)),
                plot.title = element_text(hjust = 0.5)) +
          guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
          xlab(x_axis_label) +
          ylab(y_axis_label)
      })
    } else {
      p1 <- ggplot(tsne_df, aes(x = .data[["tsne1"]], y = .data[["tsne2"]], color = !!group.by)) +
        theme_cowplot() +
        geom_point(size = pt_size) +
        guides(color = guide_legend(override.aes = list(size = legend.size))) +
        scale_color_manual(values = colors_use) +
        theme(legend.position = "right",
              axis.text = element_text(size = rel(0.95)),
              plot.title = element_text(hjust = 0.5)) +
        guides(col = guide_legend(title = "", override.aes = list(size = 4))) +
        xlab(x_axis_label) +
        ylab(y_axis_label)
    }
  }
  if (!is.null(x = split.by) && !is.null(x = num_columns)) {
    p1 <- wrap_plots(p1) + plot_layout(nrow = num_rows, ncol = num_columns)
    return(p1)
  }
  if (!is.null(x = split.by) && is.null(x = num_columns)) {
    p1 <- wrap_plots(p1)
    return(p1)
  } else {
    return(p1)
  }
}


Plot_By_Meta_LIGER2 <- function(
    liger_object,
    colors_use = NULL,
    group.by = "dataset",
    split.by = NULL,
    title = NULL,
    pt_size = NULL,
    reduction = NULL,
    num_columns = NULL,
    shuffle = TRUE,
    shuffle_seed = 1,
    legend.size = 3,
    reorder.idents = FALSE,
    new.order = NULL,
    raster = NULL,
    raster.dpi = c(512, 512),
    ggplot_default_colors = FALSE,
    color_seed = 123
) {
  # Set reduction
  reduction <- reduction %||% Default_DimReduc_LIGER(liger_object = liger_object)

  reduc_df <- Generate_Plotting_df_LIGER2(object = liger_object, group.by = group.by, split.by = split.by, reorder.idents = reorder.idents, shuffle = shuffle, shuffle_seed = shuffle_seed, reduction = reduction)

  if (!is.null(x = split.by)) {
    list_of_splits <- unique(x = reduc_df[[split.by]])
  }

  # Get length of meta data feature
  if (!is.null(x = split.by) && !is.null(x = num_columns)) {
    split.by_length <- length(x = list_of_splits)

    # Calculate number of rows for selected number of columns
    num_rows <- ceiling(x = split.by_length / num_columns)

    # Check column and row compatibility
    if (num_columns > split.by_length) {
      cli_abort(message = c("The number of columns specified is greater than the number of meta data variables.",
                            "*" = "{.field {split.by}} only contains: {.field {split.by_length}} variables.",
                            "i" = "Please adjust {.code num_columns} to be less than or equal to: {.field {split.by_length}}.")
      )
    }
  }

  meta_length <- length(x = unique(x = liger_object@cellMeta[[group.by]]))

  # set default plot colors
  if (is.null(x = colors_use)) {
    colors_use <- scCustomize_Palette(num_groups = meta_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
  }

  # Create accurate axis labels
  x_axis_label <- names(x = reduc_df)[1]
  y_axis_label <- names(x = reduc_df)[2]

  group.by <- sym(x = group.by)

  if (isTRUE(x = raster)) {
    if (!is.null(x = split.by)) {
      p1 <- lapply(1:length(x = list_of_splits), function(x) {
        ggplot(subset(reduc_df, reduc_df[[split.by]] %in% list_of_splits[x]), aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = !!group.by)) +
          theme_cowplot() +
          geom_scattermore(pointsize = pt_size, pixels = raster.dpi) +
          guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(list_of_splits[x]) +
          scale_color_manual(values = colors_use) +
          theme(legend.position = "right",
                axis.text = element_text(size = rel(0.95)),
                plot.title = element_text(hjust = 0.5)) +
          guides(col = guide_legend(title = "", override.aes = list(size = 4)))
      })
    } else {
      p1 <- ggplot(reduc_df, aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = !!group.by)) +
        theme_cowplot() +
        geom_scattermore(pointsize = pt_size, pixels = raster.dpi) +
        guides(color = guide_legend(override.aes = list(size = legend.size))) +
        scale_color_manual(values = colors_use) +
        theme(legend.position = "right",
              axis.text = element_text(size = rel(0.95)),
              plot.title = element_text(hjust = 0.5)) +
        guides(col = guide_legend(title = "", override.aes = list(size = 4)))

    }
  } else {
    if (!is.null(x = split.by)) {
      p1 <- lapply(1:length(x = list_of_splits), function(x) {
        ggplot(subset(reduc_df, reduc_df[[split.by]] %in% list_of_splits[x]), aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = !!group.by)) +
          theme_cowplot() +
          geom_point(size = pt_size) +
          guides(color = guide_legend(override.aes = list(size = legend.size))) +
          ggtitle(list_of_splits[x]) +
          scale_color_manual(values = colors_use) +
          theme(legend.position = "right",
                axis.text = element_text(size = rel(0.95)),
                plot.title = element_text(hjust = 0.5)) +
          guides(col = guide_legend(title = "", override.aes = list(size = 4)))
      })
    } else {
      p1 <- ggplot(reduc_df, aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = !!group.by)) +
        theme_cowplot() +
        geom_point(size = pt_size) +
        guides(color = guide_legend(override.aes = list(size = legend.size))) +
        scale_color_manual(values = colors_use) +
        theme(legend.position = "right",
              axis.text = element_text(size = rel(0.95)),
              plot.title = element_text(hjust = 0.5)) +
        guides(col = guide_legend(title = "", override.aes = list(size = 4)))
    }
  }
  if (!is.null(x = split.by) && !is.null(x = num_columns)) {
    p1 <- wrap_plots(p1) + plot_layout(nrow = num_rows, ncol = num_columns)
    return(p1)
  }
  if (!is.null(x = split.by) && is.null(x = num_columns)) {
    p1 <- wrap_plots(p1)
    return(p1)
  } else {
    return(p1)
  }
}


#' Customized version of plotFactors
#'
#' Modified and optimized version of `plotFactors` function from LIGER package.
#'
#' @param liger_object \code{liger} liger_object.  Need to perform clustering and factorization before calling this function
#' @param num_genes Number of genes to display for each factor (Default 8).
#' @param colors_use_factors colors to use for plotting factor loadings  By default datasets will be
#' plotted using "varibow" with shuffle = TRUE from both from \code{\link{DiscretePalette_scCustomize}}.
#' @param colors_use_dimreduc colors to use for plotting factor loadings on dimensionality reduction
#' coordinates (tSNE/UMAP).  Default is c('lemonchiffon', 'red'),
#' @param pt.size_factors Adjust point size for plotting in the factor plots.
#' @param pt.size_dimreduc Adjust point size for plotting in dimensionality reduction plots.
#' @param reduction Name of dimensionality reduction to use for plotting.
#' @param plot_legend logical, whether to plot the legend on factor loading plots, default is TRUE.
#' Helpful if number of datasets is large to avoid crowding the plot with legend.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param order logical. Whether to plot higher loading cells on top of cells with lower loading values in the
#' dimensionality reduction plots (Default = FALSE).
#' @param plot_dimreduc logical.  Whether to plot factor loadings on dimensionality reduction coordinates.  Default is TRUE.
#' @param save_plots logical.  Whether to save plots.  Default is TRUE
#' @param file_path directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix to append after sample name.
#' @param return_plots logical. Whether or not to return plots to the environment.  (Default is FALSE)
#' @param cells.highlight Names of specific cells to highlight in plot (black) (default NULL).
#' @param ggplot_default_colors logical.  If `colors_use_factors = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "varibow" palette.
#' @param color_seed random seed for the palette shuffle if `colors_use_factors = NULL`.  Default = 123.
#'
#' @return A list of ggplot/patchwork objects and/or PDF file.
#'
#' @import cli
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom patchwork wrap_plots
#' @importFrom scattermore geom_scattermore
#'
#' @noRd
#'
#' @concept liger_plotting
#'
#' @author Velina Kozareva (Original code for modified function), Sam Marsh (Added/modified functionality)
#' @references Based on `plotFactors` functionality from original LIGER package.
#'
#' @examples
#' \dontrun{
#' plotFactors_scCustom(liger_object = liger_obj, return_plots = FALSE, plot_dimreduc = TRUE,
#' raster = FALSE, save_plots = TRUE)
#' }
#'

plotFactors_liger2_scCustom <- function(
    liger_object,
    num_genes = 8,
    colors_use_factors = NULL,
    colors_use_dimreduc = c("lemonchiffon", "red"),
    pt.size_factors = 1,
    pt.size_dimreduc = 1,
    reduction = "UMAP",
    plot_legend = TRUE,
    raster = TRUE,
    raster.dpi = c(512, 512),
    order = FALSE,
    plot_dimreduc = TRUE,
    save_plots = TRUE,
    file_path = NULL,
    file_name = NULL,
    return_plots = FALSE,
    cells.highlight = NULL,
    ggplot_default_colors = FALSE,
    color_seed = 123
) {
  # if returning and saving
  if (isTRUE(x = save_plots)) {
    # Check file path is valid
    if (!is.null(x = file_path) && file_path != "") {
      if (!dir.exists(paths = file_path)) {
        cli_abort(message = "Provided {.code file_path}: {.val {file_path}} does not exist.")
      }
    }

    # Set file_path before path check if current dir specified as opposed to leaving set to NULL
    if (is.null(x = file_path)) {
      file_path <- ""
    }

    # Check if file name provided
    file_ext <- grep(x = file_name, pattern = ".pdf$", ignore.case = TRUE)
    if (length(x = file_ext) == 0) {
      file_name <- file_name
    } else {
      file_name <- gsub(pattern = ".pdf", replacement = "", x = file_name, ignore.case = TRUE)
    }

    if (is.null(x = file_name)) {
      cli_abort(message = c("No file name provided.",
                            "i" = "Please provide a file name using {.code file_name}.")
      )
    }
  }

  # Extract dataset number
  num_datasets <- length(x = liger_object@datasets)

  # Default Colors for Factor Plots
  if (is.null(x = colors_use_factors)) {
    if (isTRUE(x = ggplot_default_colors)) {
      colors_use_factors <- Hue_Pal(num_colors = num_datasets)
    } else {
      colors_use_factors <- DiscretePalette_scCustomize(num_colors = num_datasets, palette = "varibow", shuffle_pal = TRUE, seed = color_seed)
    }
  }

  # Check valid number of colors for tsne/UMAP
  if (length(x = colors_use_dimreduc) < 2) {
    cli_abort(message = c("Less than two values provided to {.code colors_use_dimreduc}.",
                          "i" = "Must provided either two colors to use for creating a gradient or a larger color gradient.")
    )
  }

  # Get Data and Plot Factors
  k <- ncol(x = liger_object@H.norm)
  if (is.null(x = k)) {
    cli_abort(message = "{.code quantileNorm} must be run before plotting factors.")
  }

  cli_inform(message = "{.field Generating plots}")
  pb <- txtProgressBar(min = 0, max = k, style = 3)
  W <- liger_object@W
  rownames(x = W) <- rownames(x = liger_object@datasets[[1]]@scaleData)
  Hs_norm <- liger_object@H.norm
  dataset_names <- names(liger_object@datasets)
  H_raw_list <- lapply(1:num_datasets, function(x){
    H_raw <- t(liger_object@datasets[[x]]@H)
  })
  H_raw = do.call(rbind, H_raw_list)
  # Create accurate axis labels
  reduc_check <- Embeddings(object = liger_object, reduction = reduction, check_only = TRUE)

  x_axis_label <- paste0(reduction, "_1")
  y_axis_label <- paste0(reduction, "_2")
  plot_list = list()
  tsne_list = list()
  for (i in 1:k) {
    top_genes.W <- rownames(x = W)[order(W[, i], decreasing = T)[1:num_genes]]
    top_genes.W.string <- paste0(top_genes.W, collapse = ", ")
    factor_textstring <- paste0("Factor", i)
    plot_title1 <- paste(factor_textstring, "\n", top_genes.W.string, "\n")
    h_df = data.frame(x = 1:nrow(Hs_norm), h_norm = Hs_norm[, i],
                      h_raw = H_raw[, i], dataset = liger_object@cellMeta$dataset,
                      highlight = FALSE)
    if (isTRUE(x = raster)) {
      top <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_raw"]], col = .data[["dataset"]])) +
        geom_scattermore(pointsize = pt.size_factors, pixels = raster.dpi) +
        labs(x = "Cell", y = "Raw H Score") +
        ggtitle(plot_title1) +
        theme(legend.position = "none") +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        top <- top + NoLegend()
      }

      bottom <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_norm"]], col = .data[["dataset"]])) +
        geom_scattermore(pointsize = pt.size_factors, pixels = raster.dpi) +
        labs(x = "Cell", y = "H_norm Score") +
        theme(legend.position = "top",
              legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2))) +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        bottom <- bottom + NoLegend()
      }

    } else {
      top <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_raw"]], col = .data[["dataset"]])) +
        geom_point(size = pt.size_factors) +
        labs(x = "Cell", y = "Raw H Score") +
        ggtitle(plot_title1) +
        theme(legend.position = "none") +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        top <- top + NoLegend()
      }

      bottom <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_norm"]], col = .data[["dataset"]])) +
        geom_point(size = pt.size_factors) +
        labs(x = "Cell", y = "H_norm Score") +
        theme(legend.position = "top",
              legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2))) +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        bottom <- bottom + NoLegend()
      }

    }

    if (!is.null(cells.highlight)) {
      h_df[cells.highlight, "highlight"] = TRUE
      if (isTRUE(x = raster)) {
        top <- top + geom_scattermore(data = subset(h_df, .data[["highlight"]] == TRUE),
                                      aes(.data[["x"]], .data[["h_raw"]]),
                                      col = "black",
                                      pointsize = pt.size_factors,
                                      pixels = raster.dpi)
        bottom <- bottom + geom_scattermore(data = subset(h_df, .data[["highlight"]] == TRUE),
                                            aes(.data[["x"]], .data[["h_norm"]]),
                                            col = "black",
                                            pointsize = pt.size_factors,
                                            pixels = raster.dpi)
      } else {
        top <- top + geom_point(data = subset(h_df, .data[["highlight"]] == TRUE),
                                aes(.data[["x"]], .data[["h_raw"]]),
                                col = "black",
                                size = pt.size_factors)
        bottom <- bottom + geom_point(data = subset(h_df, .data[["highlight"]] == TRUE),
                                      aes(.data[["x"]], .data[["h_norm"]]),
                                      col = "black",
                                      size = pt.size_factors)
      }
    }
    full <- wrap_plots(top, bottom, ncol = 1)
    plot_list[[i]] = full

    # plot tSNE/UMAP
    if (isTRUE(x = plot_dimreduc)) {
      tsne_df <- data.frame(Hs_norm[, i], Embeddings(object = liger_object, reduction = reduction))
      factorlab <- paste0("Factor", i)
      colnames(x = tsne_df) <- c(factorlab, x_axis_label, y_axis_label)

      if (isTRUE(x = order)) {
        tsne_df <- tsne_df[order(tsne_df[,1], decreasing = FALSE),]
      }

      if (isTRUE(x = raster)) {
        p1 <- ggplot(tsne_df, aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = .data[[factorlab]])) +
          geom_scattermore(pointsize = pt.size_dimreduc, pixels = raster.dpi) +
          ggtitle(label = paste("Factor", i)) +
          theme(legend.position = "none") +
          xlab(x_axis_label) +
          ylab(y_axis_label) +
          if (length(x = colors_use_dimreduc) == 2) {
            scale_color_gradient(low = colors_use_dimreduc[1], high = colors_use_dimreduc[2])
          } else {
            scale_color_gradientn(colours = colors_use_dimreduc)
          }
      } else {
        p1 <- ggplot(tsne_df, aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = .data[[factorlab]])) +
          geom_point(size = pt.size_dimreduc) +
          ggtitle(label = paste("Factor", i)) +
          theme(legend.position = "none") +
          xlab(x_axis_label) +
          ylab(y_axis_label) +
          if (length(x = colors_use_dimreduc) == 2) {
            scale_color_gradient(low = colors_use_dimreduc[1], high = colors_use_dimreduc[2])
          } else {
            scale_color_gradientn(colours = colors_use_dimreduc)
          }
      }

      tsne_list[[i]] = p1
    }
    setTxtProgressBar(pb, i)
  }

  # save plots
  if (isTRUE(x = save_plots)) {
    cli_inform(message = "{.field Saving plots to file}")
    pdf(paste(file_path, file_name, ".pdf", sep=""))
    pb <- txtProgressBar(min = 0, max = length(x = 1:k), style = 3, file = stderr())
    for (i in 1:k) {
      if (isTRUE(x = plot_dimreduc)) {
        print(plot_list[[i]])
        print(tsne_list[[i]])
        setTxtProgressBar(pb = pb, value = i)
      } else {
        print(plot_list[[i]])
        setTxtProgressBar(pb = pb, value = i)
      }
    }
    close(con = pb)
    dev.off()
  }

  # return plots
  if (isTRUE(x = return_plots)) {
    return(list(factor_plots = plot_list,
                dimreduc_plots = tsne_list))
  }
}


#' Customized version of plotFactors
#'
#' Modified and optimized version of `plotFactors` function from LIGER package.
#'
#' @param liger_object \code{liger} liger_object.  Need to perform clustering and factorization before calling this function
#' @param num_genes Number of genes to display for each factor (Default 8).
#' @param colors_use_factors colors to use for plotting factor loadings  By default datasets will be
#' plotted using "varibow" with shuffle = TRUE from both from \code{\link{DiscretePalette_scCustomize}}.
#' @param colors_use_dimreduc colors to use for plotting factor loadings on dimensionality reduction
#' coordinates (tSNE/UMAP).  Default is c('lemonchiffon', 'red'),
#' @param pt.size_factors Adjust point size for plotting in the factor plots.
#' @param pt.size_dimreduc Adjust point size for plotting in dimensionality reduction plots.
#' @param reduction_label What to label the x and y axes of resulting plots.  LIGER does not store name of
#' technique and therefore needs to be set manually.  Default is "UMAP".
#' @param plot_legend logical, whether to plot the legend on factor loading plots, default is TRUE.
#' Helpful if number of datasets is large to avoid crowding the plot with legend.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param order logical. Whether to plot higher loading cells on top of cells with lower loading values in the
#' dimensionality reduction plots (Default = FALSE).
#' @param plot_dimreduc logical.  Whether to plot factor loadings on dimensionality reduction coordinates.  Default is TRUE.
#' @param save_plots logical.  Whether to save plots.  Default is TRUE
#' @param file_path directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix to append after sample name.
#' @param return_plots logical. Whether or not to return plots to the environment.  (Default is FALSE)
#' @param cells.highlight Names of specific cells to highlight in plot (black) (default NULL).
#' @param reorder_datasets New order to plot datasets in for the factor plots if different from current
#' factor level order in cell.data slot.
#' @param ggplot_default_colors logical.  If `colors_use_factors = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "varibow" palette.
#' @param color_seed random seed for the palette shuffle if `colors_use_factors = NULL`.  Default = 123.
#'
#' @return A list of ggplot/patchwork objects and/or PDF file.
#'
#' @import cli
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom patchwork wrap_plots
#' @importFrom scattermore geom_scattermore
#'
#' @noRd
#'
#' @concept liger_plotting
#'
#' @author Velina Kozareva (Original code for modified function), Sam Marsh (Added/modified functionality)
#' @references Based on `plotFactors` functionality from original LIGER package.
#'
#' @examples
#' \dontrun{
#' plotFactors_scCustom(liger_object = liger_obj, return_plots = FALSE, plot_dimreduc = TRUE,
#' raster = FALSE, save_plots = TRUE)
#' }
#'

plotFactors_liger_scCustom <- function(
    liger_object,
    num_genes = 8,
    colors_use_factors = NULL,
    colors_use_dimreduc = c('lemonchiffon', 'red'),
    pt.size_factors = 1,
    pt.size_dimreduc = 1,
    reduction_label = "UMAP",
    plot_legend = TRUE,
    raster = TRUE,
    raster.dpi = c(512, 512),
    order = FALSE,
    plot_dimreduc = TRUE,
    save_plots = TRUE,
    file_path = NULL,
    file_name = NULL,
    return_plots = FALSE,
    cells.highlight = NULL,
    reorder_datasets = NULL,
    ggplot_default_colors = FALSE,
    color_seed = 123
) {
  # if returning and saving
  if (isTRUE(x = save_plots)) {

    # Check file path is valid
    if (!is.null(x = file_path) && file_path != "") {
      if (!dir.exists(paths = file_path)) {
        cli_abort(message = "Provided {.code file_path}: {.val {file_path}} does not exist.")
      }
    }

    # Set file_path before path check if current dir specified as opposed to leaving set to NULL
    if (is.null(x = file_path)) {
      file_path <- ""
    }

    # Check if file name provided
    file_ext <- grep(x = file_name, pattern = ".pdf$", ignore.case = TRUE)
    if (length(x = file_ext) == 0) {
      file_name <- file_name
    } else {
      file_name <- gsub(pattern = ".pdf", replacement = "", x = file_name, ignore.case = TRUE)
    }

    if (is.null(x = file_name)) {
      cli_abort(message = c("No file name provided.",
                            "i" = "Please provide a file name using {.code file_name}.")
      )
    }
  }

  if (!is.null(x = reorder_datasets)) {
    # Check new order contains same dataset names and number of datasets
    if (length(x = levels(x = liger_object@cell.data$dataset)) != length(x = reorder_datasets)) {
      cli_abort(message = c("Error reordering datasets (number mismatch).",
                            "i" = "The number of datasets provided to {.code reorder_datasets} ({.field {length(x = reorder_datasets)}}) does not match number of datasets in LIGER object ({.field {length(x = levels(x = levels(liger_object@cell.data$dataset)))}}).")
      )
    } else {
      if (!all(levels(x = liger_object@cell.data$dataset) %in% reorder_datasets)) {
        cli_abort(message = c("Error reordering datasets (name mismatch).",
                              "*" = "Dataset names provided to {.code reorder_datasets} do not match names of datasets in LIGER object.",
                              "i" = "Please check spelling.")
        )
      } else {
        liger_object@cell.data$dataset <- factor(x = liger_object@cell.data$dataset, levels = reorder_datasets)
      }
    }
  }

  # Create accurate axis labels
  x_axis_label <- paste0(reduction_label, "_1")
  y_axis_label <- paste0(reduction_label, "_2")

  # Extract dataset number
  num_datasets <- length(x = liger_object@scale.data)

  # Default Colors for Factor Plots
  if (is.null(x = colors_use_factors)) {
    if (isTRUE(x = ggplot_default_colors)) {
      colors_use_factors <- Hue_Pal(num_colors = num_datasets)
    } else {
      colors_use_factors <- DiscretePalette_scCustomize(num_colors = num_datasets, palette = "varibow", shuffle_pal = TRUE, seed = color_seed)
    }
  }

  # Check valid number of colors for tsne/UMAP
  if (length(x = colors_use_dimreduc) < 2) {
    cli_abort(message = c("Less than two values provided to {.code colors_use_dimreduc}.",
                          "i" = "Must provided either two colors to use for creating a gradient or a larger color gradient.")
    )
  }

  # Add one time dim label warning
  if (getOption(x = "scCustomize_warn_LIGER_dim_labels_plotFactors", default = TRUE)) {
    cli_inform(message = c("",
                           "NOTE: {.field plotFactors_scCustom} uses the {.code reduction_label} parameter to set axis labels",
                           "on the dimensionality reduction plots.",
                           "By default this is set to {.val UMAP}.",
                           "Please take note of this parameter as LIGER objects do not store the name",
                           "of reduction technique used and therefore this needs to be set manually.",
                           "",
                           "-----This message will be shown once per session.-----"))
    options(scCustomize_warn_LIGER_dim_labels_plotFactors = FALSE)
  }

  # Get Data and Plot Factors
  cli_inform(message = "{.field Generating plots}")
  k <- ncol(x = liger_object@H.norm)
  pb <- txtProgressBar(min = 0, max = k, style = 3)
  W <- t(x = liger_object@W)
  rownames(x = W) <- colnames(x = liger_object@scale.data[[1]])
  Hs_norm <- liger_object@H.norm
  H_raw = do.call(rbind, liger_object@H)
  plot_list = list()
  tsne_list = list()
  for (i in 1:k) {
    top_genes.W <- rownames(x = W)[order(W[, i], decreasing = T)[1:num_genes]]
    top_genes.W.string <- paste0(top_genes.W, collapse = ", ")
    factor_textstring <- paste0("Factor", i)
    plot_title1 <- paste(factor_textstring, "\n", top_genes.W.string, "\n")
    h_df = data.frame(x = 1:nrow(Hs_norm), h_norm = Hs_norm[, i],
                      h_raw = H_raw[, i], dataset = liger_object@cell.data$dataset,
                      highlight = FALSE)
    if (isTRUE(x = raster)) {
      top <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_raw"]], col = .data[["dataset"]])) +
        geom_scattermore(pointsize = pt.size_factors, pixels = raster.dpi) +
        labs(x = "Cell", y = "Raw H Score") +
        ggtitle(plot_title1) +
        theme(legend.position = "none") +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        top <- top + NoLegend()
      }

      bottom <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_norm"]], col = .data[["dataset"]])) +
        geom_scattermore(pointsize = pt.size_factors, pixels = raster.dpi) +
        labs(x = "Cell", y = "H_norm Score") +
        theme(legend.position = "top",
              legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2))) +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        bottom <- bottom + NoLegend()
      }

    } else {
      top <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_raw"]], col = .data[["dataset"]])) +
        geom_point(size = pt.size_factors) +
        labs(x = "Cell", y = "Raw H Score") +
        ggtitle(plot_title1) +
        theme(legend.position = "none") +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        top <- top + NoLegend()
      }

      bottom <- ggplot(h_df, aes(x = .data[["x"]], y=.data[["h_norm"]], col = .data[["dataset"]])) +
        geom_point(size = pt.size_factors) +
        labs(x = "Cell", y = "H_norm Score") +
        theme(legend.position = "top",
              legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2))) +
        scale_color_manual(values = colors_use_factors)

      if (isFALSE(x = plot_legend)) {
        bottom <- bottom + NoLegend()
      }

    }

    if (!is.null(cells.highlight)) {
      h_df[cells.highlight, "highlight"] = TRUE
      if (isTRUE(x = raster)) {
        top <- top + geom_scattermore(data = subset(h_df, .data[["highlight"]] == TRUE),
                                      aes(.data[["x"]], .data[["h_raw"]]),
                                      col = "black",
                                      pointsize = pt.size_factors,
                                      pixels = raster.dpi)
        bottom <- bottom + geom_scattermore(data = subset(h_df, .data[["highlight"]] == TRUE),
                                            aes(.data[["x"]], .data[["h_norm"]]),
                                            col = "black",
                                            pointsize = pt.size_factors,
                                            pixels = raster.dpi)
      } else {
        top <- top + geom_point(data = subset(h_df, .data[["highlight"]] == TRUE),
                                aes(.data[["x"]], .data[["h_raw"]]),
                                col = "black",
                                size = pt.size_factors)
        bottom <- bottom + geom_point(data = subset(h_df, .data[["highlight"]] == TRUE),
                                      aes(.data[["x"]], .data[["h_norm"]]),
                                      col = "black",
                                      size = pt.size_factors)
      }
    }
    full <- wrap_plots(top, bottom, ncol = 1)
    plot_list[[i]] = full

    # plot tSNE/UMAP
    if (isTRUE(x = plot_dimreduc)) {
      tsne_df <- data.frame(Hs_norm[, i], liger_object@tsne.coords)
      factorlab <- paste0("Factor", i)
      colnames(x = tsne_df) <- c(factorlab, x_axis_label, y_axis_label)

      if (isTRUE(x = order)) {
        tsne_df <- tsne_df[order(tsne_df[,1], decreasing = FALSE),]
      }

      if (isTRUE(x = raster)) {
        p1 <- ggplot(tsne_df, aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = .data[[factorlab]])) +
          geom_scattermore(pointsize = pt.size_dimreduc, pixels = raster.dpi) +
          ggtitle(label = paste("Factor", i)) +
          theme(legend.position = "none") +
          xlab(x_axis_label) +
          ylab(y_axis_label) +
          if (length(x = colors_use_dimreduc) == 2) {
            scale_color_gradient(low = colors_use_dimreduc[1], high = colors_use_dimreduc[2])
          } else {
            scale_color_gradientn(colours = colors_use_dimreduc)
          }
      } else {
        p1 <- ggplot(tsne_df, aes(x = .data[[x_axis_label]], y = .data[[y_axis_label]], color = .data[[factorlab]])) +
          geom_point(size = pt.size_dimreduc) +
          ggtitle(label = paste("Factor", i)) +
          theme(legend.position = "none") +
          xlab(x_axis_label) +
          ylab(y_axis_label) +
          if (length(x = colors_use_dimreduc) == 2) {
            scale_color_gradient(low = colors_use_dimreduc[1], high = colors_use_dimreduc[2])
          } else {
            scale_color_gradientn(colours = colors_use_dimreduc)
          }
      }

      tsne_list[[i]] = p1
    }
    setTxtProgressBar(pb, i)
  }

  # save plots
  if (isTRUE(x = save_plots)) {
    cli_inform(message = "{.field Saving plots to file}")
    pdf(paste(file_path, file_name, ".pdf", sep=""))
    pb <- txtProgressBar(min = 0, max = length(x = 1:k), style = 3, file = stderr())
    for (i in 1:k) {
      if (isTRUE(x = plot_dimreduc)) {
        print(plot_list[[i]])
        print(tsne_list[[i]])
        setTxtProgressBar(pb = pb, value = i)
      } else {
        print(plot_list[[i]])
        setTxtProgressBar(pb = pb, value = i)
      }
    }
    close(con = pb)
    dev.off()
  }

  # return plots
  if (isTRUE(x = return_plots)) {
    return(list(factor_plots = plot_list,
                dimreduc_plots = tsne_list))
  }
}


#' DimPlot LIGER Version
#'
#' Standard and modified version of LIGER's plotByDatasetAndCluster
#'
#' @param liger_object \code{liger} liger_object.  Need to perform clustering before calling this function
#' @param group.by Variable to be plotted.  If `NULL` will plot clusters from `liger@clusters` slot.
#' If `combination = TRUE` will plot both clusters and meta data variable.
#' @param split.by Variable to split plots by.
#' @param colors_use_cluster colors to use for plotting by clusters.  By default if number of levels plotted is
#' less than or equal to 36 will use "polychrome" and if greater than 36 will use "varibow" with shuffle = TRUE
#' both from \code{\link{DiscretePalette_scCustomize}}.
#' @param colors_use_meta colors to use for plotting by meta data (cell.data) variable.  By default if number
#' of levels plotted is less than or equal to 36 it will use "polychrome" and if greater than 36 will use
#' "varibow" with shuffle = TRUE both from DiscretePalette_scCustomize.
#' @param pt_size Adjust point size for plotting.
#' @param shuffle logical. Whether to randomly shuffle the order of points. This can be useful for crowded plots
#' if points of interest are being buried. (Default is TRUE).
#' @param shuffle_seed Sets the seed if randomly shuffling the order of points.
#' @param reduction_label What to label the x and y axes of resulting plots.  LIGER does not store name of
#' technique and therefore needs to be set manually.  Default is "UMAP".
#' @param aspect_ratio Control the aspect ratio (y:x axes ratio length).  Must be numeric value;
#' Default is NULL.
#' @param label logical.  Whether or not to label the clusters.  ONLY applies to plotting by cluster.  Default is TRUE.
#' @param label_size size of cluster labels.
#' @param label_repel logical.  Whether to repel cluster labels from each other if plotting by
#' cluster (if `group.by = NULL` or `group.by = "cluster`).  Default is FALSE.
#' @param label_box logical.  Whether to put a box around the label text (uses `geom_text` vs `geom_label`).
#' Default is FALSE.
#' @param label_color Color to use for cluster labels.  Default is "black".
#' @param combination logical, whether to return patchwork displaying both plots side by side.  (Default is FALSE).
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param num_columns Number of columns in plot layout.  Only valid if `split.by != NULL`.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#'
#' @return A ggplot/patchwork object
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots
#' @importFrom utils packageVersion
#'
#' @noRd
#'
#' @concept liger_plotting
#'
#' @examples
#' \dontrun{
#' LIGER_DimPlot(liger_object = obj_name, reduction_label = "UMAP")
#' }
#'

LIGER_DimPlot <- function(
    liger_object,
    group.by = NULL,
    split.by = NULL,
    colors_use_cluster = NULL,
    colors_use_meta = NULL,
    pt_size = NULL,
    shuffle = TRUE,
    shuffle_seed = 1,
    reduction_label = "UMAP",
    aspect_ratio = NULL,
    label = TRUE,
    label_size = NA,
    label_repel = FALSE,
    label_box = FALSE,
    label_color = "black",
    combination = FALSE,
    raster = NULL,
    raster.dpi = c(512, 512),
    num_columns = NULL,
    ggplot_default_colors = FALSE,
    color_seed = 123
) {
  # Check LIGER
  Is_LIGER(liger_object = liger_object)

  # Set group.by defaults
  if (isFALSE(x = combination) && is.null(x = group.by)) {
    group.by <- "cluster"
  }

  if (isTRUE(x = combination) && is.null(x = group.by)) {
    group.by <- "dataset"
  }

  # Group by cluster options
  cluster_options <- c("cluster", "Cluster", "clusters", "Clusters")
  if (group.by %in% cluster_options) {
    group.by <- "cluster"
  }

  # Check group.by parameter
  if (!group.by == "cluster")
    group_by_var <- Meta_Present(object = liger_object, meta_col_names = group.by, print_msg = FALSE, omit_warn = FALSE)[[1]]

  if (!is.null(x = split.by)) {
    group_by_var <- Meta_Present(object = liger_object, meta_col_names = split.by, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  if (packageVersion(pkg = "rliger") < "2.0.0") {
    # Add one time dim label warning
    if (getOption(x = "scCustomize_warn_LIGER_dim_labels", default = TRUE)) {
      cli_inform(message = c("",
                             "NOTE: {.field DimPlot_LIGER} uses the {.code reduction_label} parameter to set axis labels ",
                             "on the plot.",
                             "By default this is set to {.val UMAP}.",
                             "Please take note of this parameter as LIGER objects do not store the name",
                             "of reduction technique used and therefore this needs to be set manually.",
                             "",
                             "-----This message will be shown once per session.-----"))
      options(scCustomize_warn_LIGER_dim_labels = FALSE)
    }
  }

  # cells in object
  cells_total <- Cells(x = liger_object)

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = cells_total) > 2e5)

  if (isTRUE(x = raster) && (length(x = cells_total) > 2e5) && getOption(x = "scCustomize_warn_raster_LIGER", default = TRUE)) {
    cli_inform(message = c("",
                           "Rasterizing points since number of points exceeds 200,000.",
                           "To disable this behavior set {.code raster = FALSE}",
                           "",
                           "-----This message will be shown once per session.-----"))
    options(scCustomize_warn_raster_LIGER = FALSE)
  }

  # Add point size
  if (is.null(x = pt_size)) {
    # modified version of the AutoPointSize() function from Seurat
    pt_size <- AutoPointSize_scCustom(data = cells_total, raster = raster)
  }

  # plot combination plot
  if (isTRUE(x = combination)) {
    p1 <- Plot_By_Cluster_LIGER(liger_object = liger_object,
                                colors_use = colors_use_cluster,
                                split.by = split.by,
                                pt_size = pt_size,
                                reduction_label = reduction_label,
                                shuffle = shuffle,
                                raster = raster,
                                raster.dpi = raster.dpi,
                                ggplot_default_colors = ggplot_default_colors,
                                num_columns = num_columns,
                                shuffle_seed = shuffle_seed,
                                label_size = label_size,
                                label_repel = label_repel,
                                label_box = label_box,
                                label_color = label_color,
                                label = label,
                                color_seed = color_seed)

    p2 <- Plot_By_Meta_LIGER(liger_object = liger_object,
                             colors_use = colors_use_meta,
                             group.by = group.by,
                             pt_size = pt_size,
                             reduction_label = reduction_label,
                             num_columns = num_columns,
                             shuffle = shuffle,
                             raster = raster,
                             raster.dpi = raster.dpi,
                             ggplot_default_colors = ggplot_default_colors,
                             split.by = split.by,
                             color_seed = color_seed,
                             shuffle_seed = shuffle_seed)

    p3 <- wrap_plots(p1 + p2)

    # Aspect ratio changes
    if (!is.null(x = aspect_ratio)) {
      if (!is.numeric(x = aspect_ratio)) {
        cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
      }
      p3 <- p3 & theme(aspect.ratio = aspect_ratio)
    }

    return(p3)
  }

  # Plot by cluster
  if (group.by == "cluster") {
    p1 <- Plot_By_Cluster_LIGER(liger_object = liger_object,
                                colors_use = colors_use_cluster,
                                split.by = split.by,
                                pt_size = pt_size,
                                reduction_label = reduction_label,
                                shuffle = shuffle,
                                raster = raster,
                                raster.dpi = raster.dpi,
                                ggplot_default_colors = ggplot_default_colors,
                                num_columns = num_columns,
                                shuffle_seed = shuffle_seed,
                                label_size = label_size,
                                label_repel = label_repel,
                                label_box = label_box,
                                label_color = label_color,
                                label = label,
                                color_seed = color_seed)
    # Aspect ratio changes
    if (!is.null(x = aspect_ratio)) {
      if (!is.numeric(x = aspect_ratio)) {
        cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
      }
      p1 <- p1 & theme(aspect.ratio = aspect_ratio)
    }

    return(p1)
  }

  # Plot by Meta
  if (group.by != "cluster") {
    p2 <- Plot_By_Meta_LIGER(liger_object = liger_object,
                             colors_use = colors_use_meta,
                             group.by = group.by,
                             pt_size = pt_size,
                             reduction_label = reduction_label,
                             num_columns = num_columns,
                             shuffle = shuffle,
                             raster = raster,
                             raster.dpi = raster.dpi,
                             ggplot_default_colors = ggplot_default_colors,
                             split.by = split.by,
                             shuffle_seed = shuffle_seed,
                             color_seed = color_seed)
    # Aspect ratio changes
    if (!is.null(x = aspect_ratio)) {
      if (!is.numeric(x = aspect_ratio)) {
        cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
      }
      p2 <- p2 & theme(aspect.ratio = aspect_ratio)
    }

    return(p2)
  }
}



#' DimPlot LIGER Version
#'
#' Standard and modified version of LIGER's plotByDatasetAndCluster
#'
#' @param liger_object \code{liger} liger_object.  Need to perform clustering before calling this function
#' @param group.by Variable to be plotted.  If `NULL` will plot clusters from `liger@clusters` slot.
#' If `combination = TRUE` will plot both clusters and meta data variable.
#' @param split.by Variable to split plots by.
#' @param colors_use_cluster colors to use for plotting by clusters.  By default if number of levels plotted is
#' less than or equal to 36 will use "polychrome" and if greater than 36 will use "varibow" with shuffle = TRUE
#' both from \code{\link{DiscretePalette_scCustomize}}.
#' @param colors_use_meta colors to use for plotting by meta data (cell.data) variable.  By default if number
#' of levels plotted is less than or equal to 36 it will use "polychrome" and if greater than 36 will use
#' "varibow" with shuffle = TRUE both from DiscretePalette_scCustomize.
#' @param pt_size Adjust point size for plotting.
#' @param shuffle logical. Whether to randomly shuffle the order of points. This can be useful for crowded plots
#' if points of interest are being buried. (Default is TRUE).
#' @param shuffle_seed Sets the seed if randomly shuffling the order of points.
#' @param reduction specify reduction to use when plotting.  Default is current object
#' default reduction.
#' @param aspect_ratio Control the aspect ratio (y:x axes ratio length).  Must be numeric value;
#' Default is NULL.
#' @param label logical.  Whether or not to label the clusters.  ONLY applies to plotting by cluster.  Default is TRUE.
#' @param label_size size of cluster labels.
#' @param label_repel logical.  Whether to repel cluster labels from each other if plotting by
#' cluster (if `group.by = NULL` or `group.by = "cluster`).  Default is FALSE.
#' @param label_box logical.  Whether to put a box around the label text (uses `geom_text` vs `geom_label`).
#' Default is FALSE.
#' @param label_color Color to use for cluster labels.  Default is "black".
#' @param combination logical, whether to return patchwork displaying both plots side by side.  (Default is FALSE).
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore().
#' Default is c(512, 512).
#' @param num_columns Number of columns in plot layout.  Only valid if `split.by != NULL`.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#'
#' @return A ggplot/patchwork object
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots
#' @importFrom utils packageVersion
#'
#' @noRd
#'
#' @concept liger_plotting
#'
#' @examples
#' \dontrun{
#' LIGER2_DimPlot(liger_object = obj_name, reduction_label = "UMAP")
#' }
#'

LIGER2_DimPlot <- function(
    liger_object,
    group.by = NULL,
    split.by = NULL,
    colors_use_cluster = NULL,
    colors_use_meta = NULL,
    pt_size = NULL,
    shuffle = TRUE,
    shuffle_seed = 1,
    reduction = reduction,
    aspect_ratio = NULL,
    label = TRUE,
    label_size = NA,
    label_repel = FALSE,
    label_box = FALSE,
    label_color = "black",
    combination = FALSE,
    raster = NULL,
    raster.dpi = c(512, 512),
    num_columns = NULL,
    ggplot_default_colors = FALSE,
    color_seed = 123
) {
  # Check LIGER
  Is_LIGER(liger_object = liger_object)

  # Set group.by defaults
  if (isFALSE(x = combination) && is.null(x = group.by)) {
    group.by <- "cluster"
  }

  if (isTRUE(x = combination) && is.null(x = group.by)) {
    group.by <- "dataset"
  }

  # Group by cluster options
  cluster_options <- c("cluster", "Cluster", "clusters", "Clusters")
  if (group.by %in% cluster_options) {
    group.by <- "cluster"
  }

  # Check group.by parameter
  if (!group.by == "cluster")
    group_by_var <- Meta_Present(object = liger_object, meta_col_names = group.by, print_msg = FALSE, omit_warn = FALSE)[[1]]

  if (!is.null(x = split.by)) {
    group_by_var <- Meta_Present(object = liger_object, meta_col_names = split.by, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  # cells in object
  cells_total <- Cells(x = liger_object)

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = cells_total) > 2e5)

  if (isTRUE(x = raster) && (length(x = cells_total) > 2e5) && getOption(x = "scCustomize_warn_raster_LIGER", default = TRUE)) {
    cli_inform(message = c("",
                           "Rasterizing points since number of points exceeds 200,000.",
                           "To disable this behavior set {.code raster = FALSE}",
                           "",
                           "-----This message will be shown once per session.-----"))
    options(scCustomize_warn_raster_LIGER = FALSE)
  }

  # Add point size
  if (is.null(x = pt_size)) {
    # modified version of the AutoPointSize() function from Seurat
    pt_size <- AutoPointSize_scCustom(data = cells_total, raster = raster)
  }

  # plot combination plot
  if (isTRUE(x = combination)) {
    p1 <- Plot_By_Cluster_LIGER2(liger_object = liger_object,
                                colors_use = colors_use_cluster,
                                split.by = split.by,
                                pt_size = pt_size,
                                reduction = reduction,
                                shuffle = shuffle,
                                raster = raster,
                                raster.dpi = raster.dpi,
                                ggplot_default_colors = ggplot_default_colors,
                                num_columns = num_columns,
                                shuffle_seed = shuffle_seed,
                                label_size = label_size,
                                label_repel = label_repel,
                                label_box = label_box,
                                label_color = label_color,
                                label = label,
                                color_seed = color_seed)

    p2 <- Plot_By_Meta_LIGER2(liger_object = liger_object,
                             colors_use = colors_use_meta,
                             group.by = group.by,
                             pt_size = pt_size,
                             reduction = reduction,
                             num_columns = num_columns,
                             shuffle = shuffle,
                             raster = raster,
                             raster.dpi = raster.dpi,
                             ggplot_default_colors = ggplot_default_colors,
                             split.by = split.by,
                             color_seed = color_seed,
                             shuffle_seed = shuffle_seed)

    p3 <- wrap_plots(p1 + p2)

    # Aspect ratio changes
    if (!is.null(x = aspect_ratio)) {
      if (!is.numeric(x = aspect_ratio)) {
        cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
      }
      p3 <- p3 & theme(aspect.ratio = aspect_ratio)
    }

    return(p3)
  }

  # Plot by cluster
  if (group.by == "cluster") {
    p1 <- Plot_By_Cluster_LIGER2(liger_object = liger_object,
                                colors_use = colors_use_cluster,
                                split.by = split.by,
                                pt_size = pt_size,
                                reduction = reduction,
                                shuffle = shuffle,
                                raster = raster,
                                raster.dpi = raster.dpi,
                                ggplot_default_colors = ggplot_default_colors,
                                num_columns = num_columns,
                                shuffle_seed = shuffle_seed,
                                label_size = label_size,
                                label_repel = label_repel,
                                label_box = label_box,
                                label_color = label_color,
                                label = label,
                                color_seed = color_seed)
    # Aspect ratio changes
    if (!is.null(x = aspect_ratio)) {
      if (!is.numeric(x = aspect_ratio)) {
        cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
      }
      p1 <- p1 & theme(aspect.ratio = aspect_ratio)
    }

    return(p1)
  }

  # Plot by Meta
  if (group.by != "cluster") {
    p2 <- Plot_By_Meta_LIGER2(liger_object = liger_object,
                             colors_use = colors_use_meta,
                             group.by = group.by,
                             pt_size = pt_size,
                             reduction = reduction,
                             num_columns = num_columns,
                             shuffle = shuffle,
                             raster = raster,
                             raster.dpi = raster.dpi,
                             ggplot_default_colors = ggplot_default_colors,
                             split.by = split.by,
                             shuffle_seed = shuffle_seed,
                             color_seed = color_seed)
    # Aspect ratio changes
    if (!is.null(x = aspect_ratio)) {
      if (!is.numeric(x = aspect_ratio)) {
        cli_abort(message = "{.code aspect_ratio} must be a {.field numeric} value.")
      }
      p2 <- p2 & theme(aspect.ratio = aspect_ratio)
    }

    return(p2)
  }
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### QC UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Add MSigDB Gene Lists Percentages
#'
#' Adds percentage of counts from 3 hallmark MSigDB hallmark gene sets: "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
#' "HALLMARK_APOPTOSIS", and "HALLMARK_DNA_REPAIR".
#'
#' @param liger_object object name.
#' @param species Species of origin for given Object.  Only accepted species are: mouse, human,
#' zebrafish, rat, drosophila, or rhesus macaque (name or abbreviation)
#' @param oxphos_name name to use for the new meta.data column containing percent MSigDB Hallmark oxidative
#' phosphorylation counts. Default is "percent_oxphos".
#' @param apop_name name to use for the new meta.data column containing percent MSigDB Hallmark apoptosis counts.
#' Default is "percent_apop".
#' @param dna_repair_name name to use for the new meta.data column containing percent MSigDB Hallmark DNA repair counts.
#' Default is "percent_oxphos".
#' @param ensembl_ids logical, whether feature names in the object are gene names or
#' ensembl IDs (default is FALSE; set TRUE if feature names are ensembl IDs).
#' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
#' function will abort if columns with any one of the names provided to `mito_name` `ribo_name` or
#' `mito_ribo_name` is present in meta.data slot.
#'
#' @return liger object
#'
#' @import cli
#'
#' @keywords internal
#'
#' @noRd
#'

Add_MSigDB_LIGER <- function(
    liger_object,
    species,
    oxphos_name = "percent_oxphos",
    apop_name = "percent_apop",
    dna_repair_name = "percent_dna_repair",
    ensembl_ids = FALSE,
    overwrite = FALSE
) {
  # Accepted species names
  accepted_names <- list(
    Mouse_Options = c("Mouse", "mouse", "Ms", "ms", "Mm", "mm"),
    Human_Options = c("Human", "human", "Hu", "hu", "Hs", "hs"),
    Marmoset_Options = c("Marmoset", "marmoset", "CJ", "Cj", "cj", NA),
    Zebrafish_Options = c("Zebrafish", "zebrafish", "DR", "Dr", "dr", NA),
    Rat_Options = c("Rat", "rat", "RN", "Rn", "rn", NA),
    Drosophila_Options = c("Drosophila", "drosophila", "DM", "Dm", "dm", NA),
    Macaque_Options = c("Macaque", "macaque", "Rhesus", "macaca", "mmulatta", NA),
    Chicken_Options = c("Chicken", "chicken", "Gallus", "gallus", "Gg", "gg")
  )

  if (!species %in% unlist(x = accepted_names)) {
    cli_inform(message = "The supplied species ({.field {species}}) is not currently supported.")
  }

  # Check liger
  Is_LIGER(liger_object = liger_object)

  # Check name collision
  if (any(duplicated(x = c(oxphos_name, apop_name, dna_repair_name)))) {
    cli_abort(message = "One or more of values provided to {.code oxphos_name}, {.code apop_name}, {.code dna_repair_name} are identical.")
  }

  # Overwrite check
  meta_names <- colnames(x = Fetch_Meta(object = liger_object))

  if (oxphos_name %in% meta_names || apop_name %in% meta_names || dna_repair_name %in% meta_names) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Columns with {.val {oxphos_name}} and/or {.val {apop_name}} already present in meta data.",
                            "i" = "*To run function and overwrite columns set parameter {.code overwrite = TRUE} or change respective {.code oxphos_name}, {.code apop_name}, and/or {.code dna_repair_name}*")
      )
    }
    cli_inform(message = c("Columns with {.val {oxphos_name}} and/or {.val {apop_name}} already present in meta data.",
                           "i" = "Overwriting those columns as {.code overwrite = TRUE.}")
    )
  }

  # Retrieve gene lists
  if (isFALSE(x = ensembl_ids)) {
    msigdb_gene_list <- Retrieve_MSigDB_Lists(species = species)
  } else {
    msigdb_gene_list <- Retrieve_MSigDB_Ensembl_Lists(species = species)
  }

  # Check features are present in object
  all_features <- Features(x = liger_object, by_dataset = FALSE)

  oxphos_found <- intersect(x = msigdb_gene_list[["oxphos"]], y = all_features)
  apop_found <- intersect(x = msigdb_gene_list[["apop"]], y = all_features)
  dna_repair_found <- intersect(x = msigdb_gene_list[["dna_repair"]], y = all_features)

  # Add meta data columns
  if (oxphos_found > 0) {
    if (packageVersion(pkg = "rliger") > "1.0.1") {
      object <- rliger::runGeneralQC(object = object, mito = FALSE, ribo = FALSE, hemo = FALSE, features = list(oxphos_name = oxphos_found), verbose = FALSE)
    } else {
      percent_oxphos <- unlist(lapply(object@raw.data, function(x) {
        (Matrix::colSums(x[oxphos_found, ]) / Matrix::colSums(x))*100}))
      object@cell.data[ , oxphos_name] <- percent_oxphos
    }
  }

  if (apop_found > 0) {
    if (packageVersion(pkg = "rliger") > "1.0.1") {
      object <- rliger::runGeneralQC(object = object, mito = FALSE, ribo = FALSE, hemo = FALSE, features = list(apop_name = apop_found), verbose = FALSE)
    } else {
      percent_apop <- unlist(lapply(object@raw.data, function(x) {
        (Matrix::colSums(x[apop_found, ]) / Matrix::colSums(x))*100}))
      object@cell.data[ , apop_name] <- percent_apop
    }
  }

  if (dna_repair_found > 0) {
    if (packageVersion(pkg = "rliger") > "1.0.1") {
      object <- rliger::runGeneralQC(object = object, mito = FALSE, ribo = FALSE, hemo = FALSE, features = list(dna_repair_name = dna_repair_found), verbose = FALSE)
    } else {
      percent_dna_repair <- unlist(lapply(object@raw.data, function(x) {
        (Matrix::colSums(x[dna_repair_found, ]) / Matrix::colSums(x))*100}))
      object@cell.data[, dna_repair_name] <- percent_dna_repair
    }
  }

  # return final object
  return(liger_object)
}


#' Add IEG Gene List Percentages
#'
#' Adds percentage of counts from IEG genes from mouse and human.
#'
#' @param liger_object object name.
#' @param species Species of origin for given Seurat Object.  Only accepted species are: mouse, human (name or abbreviation).
#' @param ieg_name name to use for the new meta.data column containing percent IEG gene counts. Default is "percent_ieg".
#' @param ensembl_ids logical, whether feature names in the object are gene names or
#' ensembl IDs (default is FALSE; set TRUE if feature names are ensembl IDs).
#' @param overwrite Logical.  Whether to overwrite existing meta data columns.  Default is FALSE meaning that
#' function will abort if columns with the name provided to `ieg_name` is present in meta data slot.
#'
#' @return liger object
#'
#' @import cli
#'
#' @keywords internal
#'
#' @noRd
#'

Add_IEG_LIGER <- function(
    liger_object,
    species,
    ieg_name = "percent_ieg",
    ensembl_ids = FALSE,
    overwrite = FALSE
) {
  # Accepted species names
  accepted_names <- list(
    Mouse_Options = c("Mouse", "mouse", "Ms", "ms", "Mm", "mm"),
    Human_Options = c("Human", "human", "Hu", "hu", "Hs", "hs"),
    Marmoset_Options = c("Marmoset", "marmoset", "CJ", "Cj", "cj", NA),
    Zebrafish_Options = c("Zebrafish", "zebrafish", "DR", "Dr", "dr", NA),
    Rat_Options = c("Rat", "rat", "RN", "Rn", "rn", NA),
    Drosophila_Options = c("Drosophila", "drosophila", "DM", "Dm", "dm", NA),
    Macaque_Options = c("Macaque", "macaque", "Rhesus", "macaca", "mmulatta", NA),
    Chicken_Options = c("Chicken", "chicken", "Gallus", "gallus", "Gg", "Gg")
  )

  if (!species %in% unlist(x = accepted_names)) {
    cli_inform(message = "The supplied species ({.field {species}}) is not currently supported.")
  }

  # Check Seurat
  Is_LIGER(liger_object = liger_object)

  # Overwrite check
  meta_names <- colnames(x = Fetch_Meta(object = liger_object))

  if (ieg_name %in% meta_names) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Column with {.val {ieg_name}} already present in meta data.",
                            "i" = "*To run function and overwrite column set parameter {.code overwrite = TRUE} or change respective {.code ieg_name}*")
      )
    }
    cli_inform(message = c("Column with {.val {ieg_name}} already present in meta data.",
                           "i" = "Overwriting those column as {.code overwrite = TRUE.}")
    )
  }

  # Retrieve gene lists
  if (isFALSE(x = ensembl_ids)) {
    ieg_gene_list <- Retrieve_IEG_Lists(species = species)
  } else {
    ieg_gene_list <- Retrieve_IEG_Ensembl_Lists(species = species)
  }

  all_features <- Features(x = liger_object, by_dataset = FALSE)

  ieg_found <- intersect(x = ieg_gene_list[["ieg"]], y = all_features)

  # Add ieg column
  if (length(x = ieg_found) > 0) {
    if (packageVersion(pkg = "rliger") > "1.0.1") {
      object <- rliger::runGeneralQC(object = object, mito = FALSE, ribo = FALSE, hemo = FALSE, features = list(ieg_name = ieg_found), verbose = FALSE)
    } else {
      percent_ieg <- unlist(lapply(object@raw.data, function(x) {
        (Matrix::colSums(x[ieg_found, ]) / Matrix::colSums(x))*100}))
      object@cell.data[, ieg_name] <- percent_ieg
    }
  }

  # return final object
  return(liger_object)
}
