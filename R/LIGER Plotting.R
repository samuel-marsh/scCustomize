#' DimPlot LIGER Version
#'
#' Standard and modified version of LIGER's plotByDatasetAndCluster
#'
#' @param liger_object \code{liger} liger_object.  Need to perform clustering before calling this function
#' @param group_by Variable to be plotted.  If `NULL` will plot clusters from `liger@clusters` slot.
#' If `combination = TRUE` will plot both clusters and meta data variable.
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
#' @param label logical.  Whether or not to label the clusters.  ONLY applies to plotting by cluster.  Default is TRUE.
#' @param label_size size of cluster labels.
#' @param label_repel logical.  Whether to repel cluster labels from each other if plotting by
#' cluster (if `group_by = NULL` or `group_by = "cluster`).  Default is FALSE.
#' @param label_box logical.  Whether to put a box around the label text (uses `geom_text` vs `geom_label`).
#' Default is FALSE.
#' @param label_color Color to use for cluster labels.  Default is "black".
#' @param combination logical, whether to return patchwork displaying both plots side by side.  (Default is FALSE).
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
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
#'
#' @export
#'
#' @concept liger_plotting
#'
#' @examples
#' \dontrun{
#' DimPlot_LIGER(liger_object = obj_name, reduction_label = "UMAP")
#' }
#'

DimPlot_LIGER <- function(
  liger_object,
  group_by = NULL,
  split_by = NULL,
  colors_use_cluster = NULL,
  colors_use_meta = NULL,
  pt_size = NULL,
  shuffle = TRUE,
  shuffle_seed = 1,
  reduction_label = "UMAP",
  label = TRUE,
  label_size = NA,
  label_repel = FALSE,
  label_box = FALSE,
  label_color = "black",
  combination = FALSE,
  raster = NULL,
  num_columns = NULL,
  ggplot_default_colors = FALSE,
  color_seed = 123
) {
  # Check LIGER
  Is_LIGER(liger_object = liger_object)

  # Set group_by defaults
  if (!combination && is.null(x = group_by)) {
    group_by <- "cluster"
  }

  if (combination && is.null(x = group_by)) {
    group_by <- "dataset"
  }

  # Group by cluster options
  cluster_options <- c("cluster", "Cluster", "clusters", "Clusters")
  if (group_by %in% cluster_options) {
    group_by <- "cluster"
  }

  # Check group_by parameter
  if (!group_by == "cluster")
    group_by_var <- Meta_Present_LIGER(liger_object = liger_object, meta_col_names = group_by, print_msg = FALSE)

  if (!is.null(x = split_by)) {
    group_by_var <- Meta_Present_LIGER(liger_object = liger_object, meta_col_names = split_by, print_msg = FALSE)
  }

  # Add one time dim label warning
  if (getOption(x = 'scCustomize_warn_LIGER_dim_labels', default = TRUE)) {
    message(
      "NOTE: DimPlot_LIGER uses the `reduction_label` parameter to set axis labels \n",
      "on the plot.\n",
      "By default this is set to UMAP.\n",
      "Please take note of this parameter as LIGER objects do not store the name\n",
      "of reduction technique used and therefore this needs to be set manually. \n
       \nThis message will be shown once per session.\n"
    )
    options(scCustomize_warn_LIGER_dim_labels = FALSE)
  }

  # Add raster check for scCustomize
  raster <- raster %||% (nrow(x = liger_object@cell.data) > 2e5)

  if (raster && (nrow(x = liger_object@cell.data) > 2e5) && getOption(x = 'scCustomize_warn_raster_LIGER', default = TRUE)) {
    message("Rasterizing points since number of points exceeds 200,000.",
            "\nTo disable this behavior set `raster=FALSE`
            \nThis message will be shown once per session.\n"
    )
    options(scCustomize_warn_raster_LIGER = FALSE)
  }

  # Add point size
  if (is.null(x = pt_size)) {
    cells_total <- nrow(x = liger_object@cell.data)
    # modified version of the AutoPointSize() function from Seurat
    pt_size <- AutoPointSize_scCustom(data = cells_total, raster = raster)
  }

  # Create accurate axis labels
  x_axis_label <- paste0(reduction_label, "_1")
  y_axis_label <- paste0(reduction_label, "_2")

  # plot combination plot
  if (combination) {
    p1 <- Plot_By_Cluster_LIGER(liger_object = liger_object,
                                colors_use = colors_use_cluster,
                                split_by = split_by,
                                pt_size = pt_size,
                                reduction_label = reduction_label,
                                shuffle = shuffle,
                                raster = raster,
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
                             group_by = group_by,
                             pt_size = pt_size,
                             reduction_label = reduction_label,
                             num_columns = num_columns,
                             shuffle = shuffle,
                             raster = raster,
                             ggplot_default_colors = ggplot_default_colors,
                             split_by = split_by,
                             color_seed = color_seed,
                             shuffle_seed = shuffle_seed)

    p3 <- wrap_plots(p1 + p2)
    return(p3)
  }

  # Plot by cluster
  if (group_by == "cluster") {
    p1 <- Plot_By_Cluster_LIGER(liger_object = liger_object,
                                colors_use = colors_use_cluster,
                                split_by = split_by,
                                pt_size = pt_size,
                                reduction_label = reduction_label,
                                shuffle = shuffle,
                                raster = raster,
                                ggplot_default_colors = ggplot_default_colors,
                                num_columns = num_columns,
                                shuffle_seed = shuffle_seed,
                                label_size = label_size,
                                label_repel = label_repel,
                                label_box = label_box,
                                label_color = label_color,
                                label = label,
                                color_seed = color_seed)
    return(p1)
  }

  # Plot by Meta
  if (group_by != "cluster") {
    p2 <- Plot_By_Meta_LIGER(liger_object = liger_object,
                             colors_use = colors_use_meta,
                             group_by = group_by,
                             pt_size = pt_size,
                             reduction_label = reduction_label,
                             num_columns = num_columns,
                             shuffle = shuffle,
                             raster = raster,
                             ggplot_default_colors = ggplot_default_colors,
                             split_by = split_by,
                             shuffle_seed = shuffle_seed,
                             color_seed = color_seed)
    return(p2)
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
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
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
#' @importFrom patchwork wrap_plots
#' @importFrom scattermore geom_scattermore
#'
#' @export
#'
#' @concept liger_plotting
#'
#' @author Velina Kozareva (Original code for modified function), Sam Marsh (Added/modified functionality)
#' @references Based on `plotFactors` functionality from original LIGER package.
#'
#' @examples
#' \dontrun{
#' plotFactors_scCustom(liger_object = liger_obj, return_plots = FALSE, plot_dimreduc = TRUE, raster = FALSE,
#' save_plots = TRUE)
#' }
#'

plotFactors_scCustom <- function(
  liger_object,
  num_genes = 8,
  colors_use_factors = NULL,
  colors_use_dimreduc = c('lemonchiffon', 'red'),
  pt.size_factors = 1,
  pt.size_dimreduc = 1,
  reduction_label = "UMAP",
  raster = TRUE,
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
  # Check LIGER
  Is_LIGER(liger_object = liger_object)

  # if returning and saving
  if (save_plots) {
    # Set file_path before path check if current dir specified as opposed to leaving set to NULL
    if (file_path == "") {
      file_path <- NULL
    }

    # Check file path is valid
    if (!is.null(x = file_path)) {
      if (!dir.exists(paths = file_path)) {
        cli_abort(message = "Provided `file_path`: '{file_path}' does not exist.")
      }
    }

    # Check if file name provided
    if (is.null(x = file_name)) {
      cli_abort(message = c("No file name provided.",
                            "i" = "Please provide a file name using `file_name`.")
      )
    }
  }

  if (!is.null(x = reorder_datasets)) {
    # Check new order contains same dataset names and number of datasets
    if (length(x = levels(x = liger_object@cell.data$dataset)) != length(x = reorder_datasets)) {
      cli_abort(message = c("Error reordering datasets (number mismatch).",
                            "i" = "The number of datasets provided to 'reorder_datasets' ({length(x = reorder_datasets)}) does not match number of datasets in LIGER object ({length(x = levels(x = levels(liger_object@cell.data$dataset)))}).")
      )
    } else {
      if (!all(levels(liger_object@cell.data$dataset) %in% reorder_datasets)) {
        cli_abort(message = c("Error reordering datasets (name mismatch).",
                              "*" = "Dataset names provided to 'reorder_datasets' do not match names of datasets in LIGER object.",
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
    if (ggplot_default_colors) {
      colors_use_factors <- Hue_Pal(num_colors = num_datasets)
    } else {
      colors_use_factors <- DiscretePalette_scCustomize(num_colors = num_datasets, palette = "varibow", shuffle_pal = TRUE, seed = color_seed)
    }
  }

  # Check valid number of colors for tsne/UMAP
  if (length(x = colors_use_dimreduc) < 2) {
    cli_abort(message = c("Less than two values provided to `colors_use_dimreduc`.",
                          "i" = "Must provided either two colors to use for creating a gradient or a larger color gradient.")
    )
  }

  # Add one time dim label warning
  if (getOption(x = 'scCustomize_warn_LIGER_dim_labels_plotFactors', default = TRUE)) {
    message(
      "NOTE: plotFactors_scCustom uses the `reduction_label` parameter to set axis labels \n",
      "on the plot.\n",
      "By default this is set to UMAP.\n",
      "Please take note of this parameter as LIGER objects do not store the name\n",
      "of reduction technique used and therefore this needs to be set manually. \n
       \nThis message will be shown once per session.\n"
    )
    options(scCustomize_warn_LIGER_dim_labels_plotFactors = FALSE)
  }

  # Get Data and Plot Factors
  cli_inform(message = "Generating plots")
  k <- ncol(liger_object@H.norm)
  pb <- txtProgressBar(min = 0, max = k, style = 3)
  W <- t(liger_object@W)
  rownames(W) <- colnames(liger_object@scale.data[[1]])
  Hs_norm <- liger_object@H.norm
  H_raw = do.call(rbind, liger_object@H)
  plot_list = list()
  tsne_list = list()
  for (i in 1:k) {
    top_genes.W <- rownames(W)[order(W[, i], decreasing = T)[1:num_genes]]
    top_genes.W.string <- paste0(top_genes.W, collapse = ", ")
    factor_textstring <- paste0("Factor", i)
    plot_title1 <- paste(factor_textstring, "\n", top_genes.W.string, "\n")
    h_df = data.frame(x = 1:nrow(Hs_norm), h_norm = Hs_norm[, i],
                      h_raw = H_raw[, i], dataset = liger_object@cell.data$dataset,
                      highlight = FALSE)
    if (raster) {
      top <- ggplot(h_df, aes(x = x, y=h_raw, col = dataset)) +
        scattermore::geom_scattermore(pointsize = pt.size_factors) +
        labs(x = 'Cell', y = 'Raw H Score') +
        ggtitle(plot_title1) +
        theme(legend.position = 'none') +
        scale_color_manual(values = colors_use_factors)

      bottom <- ggplot(h_df, aes(x = x, y=h_norm, col = dataset)) +
        scattermore::geom_scattermore(pointsize = pt.size_factors) +
        labs(x = 'Cell', y = 'H_norm Score') +
        theme(legend.position = 'top',
              legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2))) +
        scale_color_manual(values = colors_use_factors)


    } else {
      top <- ggplot(h_df, aes(x = x, y=h_raw, col = dataset)) +
        geom_point(size = pt.size_factors) +
        labs(x = 'Cell', y = 'Raw H Score') +
        ggtitle(plot_title1) +
        theme(legend.position = 'none') +
        scale_color_manual(values = colors_use_factors)

      bottom <- ggplot(h_df, aes(x = x, y=h_norm, col = dataset)) +
        geom_point(size = pt.size_factors) +
        labs(x = 'Cell', y = 'H_norm Score') +
        theme(legend.position = 'top',
              legend.title = element_blank()) +
        guides(colour = guide_legend(override.aes = list(size = 2))) +
        scale_color_manual(values = colors_use_factors)

    }

    if (!is.null(cells.highlight)) {
      h_df[cells.highlight, 'highlight'] = TRUE
      if (raster) {
        top <- top + scattermore::geom_scattermore(data = subset(h_df, highlight == TRUE),
                                                   aes(x, h_raw),
                                                   col = "black",
                                                   pointsize = pt.size_factors)
        bottom <- bottom + scattermore::geom_scattermore(data = subset(h_df, highlight == TRUE),
                                                         aes(x, h_norm),
                                                         col = "black",
                                                         pointsize = pt.size_factors)
      } else {
        top <- top + geom_point(data = subset(h_df, highlight == TRUE),
                                aes(x, h_raw),
                                col = "black",
                                size = pt.size_factors)
        bottom <- bottom + geom_point(data = subset(h_df, highlight == TRUE),
                                      aes(x, h_norm),
                                      col = "black",
                                      size = pt.size_factors)
      }
    }
    full <- wrap_plots(top, bottom, ncol = 1)
    plot_list[[i]] = full

    # plot tSNE/UMAP
    if (plot_dimreduc) {
      tsne_df <- data.frame(Hs_norm[, i], liger_object@tsne.coords)
      factorlab <- paste0("Factor", i)
      colnames(tsne_df) <- c(factorlab, x_axis_label, y_axis_label)

      if (order) {
        tsne_df <- tsne_df[order(tsne_df[,1], decreasing = FALSE),]
      }

      if (raster) {
        p1 <- ggplot(tsne_df, aes_string(x = x_axis_label, y = y_axis_label, color = factorlab)) +
          scattermore::geom_scattermore(pointsize = pt.size_dimreduc) +
          ggtitle(label = paste('Factor', i)) +
          theme(legend.position = 'none') +
          xlab(x_axis_label) +
          ylab(y_axis_label) +
          if (length(x = colors_use_dimreduc) == 2) {
            scale_color_gradient(low = colors_use_dimreduc[1], high = colors_use_dimreduc[2])
          } else {
            scale_color_gradientn(colours = colors_use_dimreduc)
          }
      } else {
        p1 <- ggplot(tsne_df, aes_string(x = x_axis_label, y = y_axis_label, color = factorlab)) +
          geom_point(size = pt.size_dimreduc) +
          ggtitle(label = paste('Factor', i)) +
          theme(legend.position = 'none') +
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
  if (save_plots) {
    cli_inform(message = "\nSaving plots to file")
    pdf(paste(file_name, ".pdf", sep=""))
    pb <- txtProgressBar(min = 0, max = length(x = 1:(2 * k)), style = 3, file = stderr())
    for (i in 1:k) {
      if (plot_dimreduc) {
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
  if (return_plots) {
    return(list(factor_plots = plot_list,
                dimreduc_plots = tsne_list))
  }
}
