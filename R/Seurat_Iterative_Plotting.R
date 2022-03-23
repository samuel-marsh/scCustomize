#' Iterate PC Loading Plots
#'
#' Plot PC Heatmaps and Dim Loadings for exploratory analysis
#'
#' @param seurat_object Seurat object name.
#' @param dims_plot number of PCs to plot (integer).  Default is all dims present in PCA.
#' @param file_path directory file path to save file.
#' @param name_prefix prefix for file name (optional).
#' @param file_name suffix for file name.  Default is "PC_Loading_Plots".
#' @param return_plots Whether to return the plot list (Default is FALSE).  Must assign to environment
#' to save plot list.
#'
#' @return A list of plots outputted as pdf
#'
#' @import patchwork
#' @import ggplot2
#' @importFrom pbapply pblapply pboptions
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @seealso \code{\link[Seurat]{PCHeatmap}} and \code{\link[Seurat]{VizDimLoadings}}
#'
#' @export
#'
#' @concept iterative_plotting
#'
#' @examples
#' \dontrun{
#' Iterate_PC_Loading_Plots(seurat_object = seurat, dims_plot = 25, file_path = "plots/")
#' }
#'

Iterate_PC_Loading_Plots <- function(
  seurat_object,
  dims_plot = NULL,
  file_path = NULL,
  name_prefix = NULL,
  file_name = "PC_Loading_Plots",
  return_plots = FALSE
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path)) {
    if (!dir.exists(paths = file_path)) {
      stop("Provided `file_path`: ", '"', file_path, '"', " does not exist.")
    }
  }

  # Set dims to plot if not specified
  num_pc_present <- length(seurat_object@reductions$pca@stdev)
  if (is.null(x = dims_plot)) {
    dims_plot <- num_pc_present
  }

  # Check pca present
  reduc_present <- names(seurat_object@reductions)
  if (!"pca" %in% reduc_present) {
    stop("Cannot find 'pca' in this Seurat Object.")
  }
  # Check dims present in object
  if (dims_plot > num_pc_present) {
    stop("The number of PCs specified to `dims_plot` (", dims_plot, ") is greater than number of PCs present in Seurat Object (", num_pc_present, ").")
  }
  # Check file path is valid
  if (!dir.exists(paths = file_path)) {
    stop("Provided `file_path`: ", '"', file_path, '"', " does not exist.")
  }

  dims_list <- 1:dims_plot
  # Create list of all plots
  message("Generating plots")
  pboptions(char = "=")
  all_plots <- pblapply(dims_list, function(x) {
    PC_Plotting(seurat_object = seurat_object, dim_number = x)
  })
  # Save plots
  message("Saving plots to file")
  pdf(file = paste(file_path, name_prefix, file_name, ".pdf", sep=""), height = 11, width = 8.5)
  pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
  for (i in 1:length(all_plots)) {
    print(all_plots[[i]])
    setTxtProgressBar(pb = pb, value = i)
  }
  close(con = pb)
  dev.off()
  if (return_plots) {
    return(all_plots)
  }
}


#' Iterate DimPlot By Sample
#'
#' Iterate Dimplot by orig.ident column from Seurat object metadata
#'
#' @param seurat_object Seurat object name.
#' @param file_path/prefix directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix to append after sample name.
#' @param file_type File type to save output as.  Must be one of following: ".pdf", ".png", ".tiff", ".jpeg", or ".svg".
#' @param single_pdf saves all plots to single PDF file (default = FALSE).  `file_type`` must be .pdf
#' @param color color scheme to use.
#' @param dpi dpi for image saving.
#' @param reduction Dimensionality Reduction to use (default is object default).
#' @param dims Dimensions to plot.
#' @param pt.size Adjust point size for plotting.
#' @param ... Extra parameters passed to \code{\link[Seurat]{DimPlot}}.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom pbapply pblapply pboptions
#' @importFrom Seurat DimPlot
#' @importFrom SeuratObject DefaultDimReduc
#' @importFrom stringr str_detect
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @concept iterative_plotting
#'
#' @examples
#' \dontrun{
#' Iterate_DimPlot_bySample(seurat_object = object, file_path = "plots/", file_name = "tsne",
#' file_type = ".jpg", dpi = 600, color = "black")
#' }
#'

Iterate_DimPlot_bySample <- function(
  seurat_object,
  file_path = NULL,
  file_name = NULL,
  file_type = NULL,
  single_pdf = FALSE,
  dpi = 600,
  color = "black",
  reduction = NULL,
  dims = c(1, 2),
  pt.size = NULL,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path)) {
    if (!dir.exists(paths = file_path)) {
      stop("Provided `file_path`: ", '"', file_path, '"', " does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name)) {
    stop("No file name provided.  Please provide a file name using `file_name`.")
  }

  # Set file type for single pdf option
  if (single_pdf && is.null(x = file_type)) {
    file_type <- ".pdf"
  }
  if (single_pdf && !is.null(x = file_type) && str_detect(file_type, ".pdf") == FALSE) {
    message("WARNING: non-PDF 'file_type' specified but 'single_pdf = TRUE' selected.  Changing file_type to .pdf for output.")
    file_type <- ".pdf"
  }

  # Check file_type parameter
  file_type_options <- c(".pdf", ".png", ".tiff", ".jpeg", ".svg")
  if (is.null(x = file_type)) {
    stop("'file_type' not specified must be one of the following: '.pdf', '.png', '.tiff', '.jpeg', '.svg'")
  }
  if (!file_type %in% file_type_options) {
    stop("'file_type' must be one of the following: '.pdf', '.png', '.tiff', '.jpeg', '.svg'")
  }

  # Extract reduction coordinates
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)
  cells <- colnames(x = seurat_object)
  reduc_coordinates <- Embeddings(object = seurat_object[[reduction]])[cells, dims]
  reduc_coordinates <- as.data.frame(x = reduc_coordinates)
  x_axis <- c(min(reduc_coordinates[, 1]),
              max(reduc_coordinates[, 1]))
  y_axis <- c(min(reduc_coordinates[, 2]),
              max(reduc_coordinates[, 2]))

  # Extract orig.ident
  column_list <- as.character(unique(seurat_object@meta.data$orig.ident))

  # Create list of cells per sample
  cells_per_sample <- lapply(column_list, function(sample) {
    row.names(seurat_object@meta.data)[which(seurat_object@meta.data$orig.ident == sample)]
  })

  # Single PDF option
  if (single_pdf == TRUE) {
    message("Generating plots")
    pboptions(char = "=")
    all_plots <- pblapply(cells_per_sample,function(cells) {DimPlot(object = seurat_object, cells = cells, group.by = "orig.ident", cols = color, reduction = reduction, pt.size = pt.size, ...) +
        xlim(x_axis) +
        ylim(y_axis)})
    message("Saving plots to file")
    pdf(paste(file_path, file_name, file_type, sep=""))
    pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
    for (i in 1:length(all_plots)) {
      print(all_plots[[i]])
      setTxtProgressBar(pb = pb, value = i)
    }
    close(con = pb)
    dev.off()
  }
  else{
    # Code for non-PDF figure
    if (str_detect(file_type, ".pdf") == FALSE) {
      message("Generating plots and saving plots to file")
      pb <- txtProgressBar(min = 0, max = length(cells_per_sample), style = 3, file = stderr())
      for (i in 1:length(cells_per_sample)) {
        DimPlot(object = seurat_object, cells = cells_per_sample[[i]], group.by = "orig.ident", cols = color, reduction = reduction, pt.size = pt.size, ...) +
          xlim(x_axis) +
          ylim(y_axis)
        suppressMessages(ggsave(filename = paste(file_path, column_list[[i]], file_name, file_type, sep=""), dpi = dpi))
        setTxtProgressBar(pb = pb, value = i)
        }
      close(con = pb)
      }
    # Code for PDF Version
    if (str_detect(file_type, ".pdf") == TRUE) {
      message("Generating plots and saving plots to file")
      pb <- txtProgressBar(min = 0, max = length(cells_per_sample), style = 3, file = stderr())
      for (i in 1:length(cells_per_sample)) {
        DimPlot(object = seurat_object, cells = cells_per_sample[[i]], group.by = "orig.ident", cols = color, reduction = reduction, pt.size = pt.size, ...) +
          xlim(x_axis) +
          ylim(y_axis)
        suppressMessages(ggsave(filename = paste(file_path, column_list[[i]], file_name, file_type, sep=""), useDingbats = FALSE))
        setTxtProgressBar(pb = pb, value = i)
        }
      close(con = pb)
    }
  }
}


#' Iterate Cluster Highlight Plot
#'
#' Iterate the create plots with cluster of interest highlighted across all cluster (active.idents)
#' in given Seurat Object
#'
#' @param seurat_object Seurat object name.
#' @param highlight_color Color to highlight cells (default "navy").  Can provide either single color to use for
#' all clusters/plots or a vector of colors equal to the number of clusters to use (in order) for the clusters/plots.
#' @param background_color non-highlighted cell colors.
#' @param pt.size point size for both highlighted cluster and background.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default).
#' @param file_path directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix to append after sample name.
#' @param file_type File type to save output as.  Must be one of following: ".pdf", ".png", ".tiff", ".jpeg", or ".svg".
#' @param single_pdf saves all plots to single PDF file (default = FALSE).  `file_type`` must be .pdf
#' @param dpi dpi for image saving.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param ... Extra parameters passed to\code{\link[Seurat]{DimPlot}}.
#'
#' @return Saved plots
#'
#' @import ggplot2
#' @importFrom pbapply pbmapply pboptions
#' @importFrom Seurat DimPlot
#' @importFrom SeuratObject DefaultDimReduc
#' @importFrom stringr str_detect
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @concept iterative_plotting
#'
#' @examples
#' \dontrun{
#' Iterate_Cluster_Highlight_Plot(seurat_object = object, highlight_color = "navy", background_color = "lightgray",
#' file_path = "path/", file_name = "name", file_type = "pdf", single_pdf = TRUE)
#' }
#'

Iterate_Cluster_Highlight_Plot <- function(
  seurat_object,
  highlight_color = "navy",
  background_color = "lightgray",
  pt.size = NULL,
  reduction = NULL,
  file_path = NULL,
  file_name = NULL,
  file_type = NULL,
  single_pdf = FALSE,
  dpi = 600,
  raster = NULL,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path)) {
    if (!dir.exists(paths = file_path)) {
      stop("Provided `file_path`: ", '"', file_path, '"', " does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name)) {
    stop("No file name provided.  Please provide a file name using `file_name`.")
  }

  # Set file type for single pdf option
  if (single_pdf && is.null(x = file_type)) {
    file_type <- ".pdf"
  }
  if (single_pdf && !is.null(x = file_type) && str_detect(file_type, ".pdf") == FALSE) {
    message("WARNING: non-PDF 'file_type' specified but 'single_pdf = TRUE' selected.  Changing file_type to .pdf for output.")
    file_type <- ".pdf"
  }

  # Check file_type parameter
  file_type_options <- c(".pdf", ".png", ".tiff", ".jpeg", ".svg")
  if (is.null(x = file_type)) {
    stop("'file_type' not specified must be one of the following: '.pdf', '.png', '.tiff', '.jpeg', '.svg'")
  }
  if (!file_type %in% file_type_options) {
    stop("'file_type' must be one of the following: '.pdf', '.png', '.tiff', '.jpeg', '.svg'")
  }

  # Extract default reduction
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = colnames(x = seurat_object)) > 2e5)

  # Get number of clusters/identities
  list_idents <- levels(x = seurat_object@active.ident)
  num_idents <- length(x = list_idents)

  # Modify ident names if needed for saving plots
  list_idents_save <- gsub(pattern = "/", replacement = "-", x = list_idents)

  # Get color palette
  if (length(x = highlight_color) == 1) {
    highlight_color <- rep(highlight_color, num_idents)
  }
  if (length(x = highlight_color) > 1) {
    if (length(x = highlight_color) != num_idents) {
      stop("Number of colors provided to `highlight_color` (", length(highlight_color), ") is not equal to the number of clusters (", num_idents, ").")
    } else {
      highlight_color <- highlight_color
    }
  }

  # Single PDF option
  if (single_pdf == TRUE) {
    message("Generating plots")
    pboptions(char = "=")
    all_plots <- pbmapply(function(x, arg1) {
      cells_to_highlight <- CellsByIdentities(seurat_object, idents = x)
      suppressMessages(DimPlot(object = seurat_object,
              cells.highlight = cells_to_highlight,
              cols.highlight = arg1,
              cols = background_color,
              sizes.highlight = pt.size,
              pt.size = pt.size,
              order = TRUE,
              reduction = reduction,
              raster = raster,
              ...))
    }, list_idents, highlight_color, SIMPLIFY = FALSE)
    message("Saving plots to file")
    pdf(paste(file_path, file_name, file_type, sep=""))
    pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
    for (i in 1:length(all_plots)) {
      print(all_plots[[i]])
      setTxtProgressBar(pb = pb, value = i)
    }
    close(con = pb)
    dev.off()
  }
  else {
    # Code for non-PDF figure
    if (str_detect(file_type, ".pdf") == FALSE) {
      message("Generating plots and saving plots to file")
      pb <- txtProgressBar(min = 0, max = num_idents, style = 3, file = stderr())
      cells_to_highlight <- CellsByIdentities(seurat_object)
      for (i in 1:length(cells_to_highlight)) {
        suppressMessages(DimPlot(object = seurat_object,
                cells.highlight = cells_to_highlight[[i]],
                cols.highlight = highlight_color[i],
                cols = background_color,
                sizes.highlight = pt.size,
                pt.size = pt.size,
                order = TRUE,
                reduction = reduction,
                raster = raster,
                ...))
        suppressMessages(ggsave(filename = paste(file_path, list_idents_save[i], "_", file_name, file_type, sep=""), dpi = dpi))
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
    # Code for PDF Version
    if (str_detect(file_type, ".pdf") == TRUE) {
      message("Generating plots and saving plots to file")
      pb <- txtProgressBar(min = 0, max = num_idents, style = 3, file = stderr())
      cells_to_highlight <- CellsByIdentities(seurat_object)
      for (i in 1:length(cells_to_highlight)) {
        suppressMessages(DimPlot(object = seurat_object,
                cells.highlight = cells_to_highlight[[i]],
                cols.highlight = highlight_color[i],
                cols = background_color,
                sizes.highlight = pt.size,
                pt.size = pt.size,
                order = TRUE,
                reduction = reduction,
                raster = raster,
                ...))
        suppressMessages(ggsave(filename = paste(file_path, list_idents_save[[i]], "_", file_name, file_type, sep=""), useDingbats = FALSE))
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
  }
}


#' Iterate Meta Highlight Plot
#'
#' Iterate the create plots with meta data variable of interest highlighted.
#'
#' @param seurat_object Seurat object name.
#' @param meta_data_column Name of the column in `seurat_object@meta.data` slot to pull value from for highlighting.
#' @param new_meta_order The order in which to plot each level within `meta_data_column` if `single_PDF` is TRUE.
#' @param meta_data_sort logical.  Whether or not to sort and relevel the levels in `meta_data_column` if
#' `single_PDF` is TRUE.  Default is TRUE.
#' @param highlight_color Color to highlight cells (default "navy").  Can provide either single color to use for
#' all clusters/plots or a vector of colors equal to the number of clusters to use (in order) for the clusters/plots.
#' @param background_color non-highlighted cell colors.
#' @param pt.size point size for both highlighted cluster and background.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default).
#' @param file_path directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix to append after sample name.
#' @param file_type File type to save output as.  Must be one of following: ".pdf", ".png", ".tiff", ".jpeg", or ".svg".
#' @param single_pdf saves all plots to single PDF file (default = FALSE).  `file_type`` must be .pdf
#' @param dpi dpi for image saving.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param ... Extra parameters passed to\code{\link[Seurat]{DimPlot}}.
#'
#' @return Saved plots
#'
#' @import ggplot2
#' @importFrom forcats fct_relevel
#' @importFrom pbapply pbmapply pboptions
#' @importFrom Seurat DimPlot
#' @importFrom SeuratObject DefaultDimReduc
#' @importFrom stringr str_detect
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @concept iterative_plotting
#'
#' @examples
#' \dontrun{
#' Iterate_Meta_Highlight_Plot(seurat_object = object, meta_data_column = "sample_id", highlight_color = "navy",
#' background_color = "lightgray", file_path = "path/", file_name = "name", file_type = "pdf", single_pdf = TRUE)
#' }
#'

Iterate_Meta_Highlight_Plot <- function(
  seurat_object,
  meta_data_column,
  new_meta_order = NULL,
  meta_data_sort = TRUE,
  highlight_color = "navy",
  background_color = "lightgray",
  pt.size = NULL,
  reduction = NULL,
  file_path = NULL,
  file_name = NULL,
  file_type = NULL,
  single_pdf = FALSE,
  dpi = 600,
  raster = NULL,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check meta data
  meta_data_column <- Meta_Present(seurat_object = seurat_object, meta_col_names = meta_data_column, omit_warn = FALSE, print_msg = FALSE)[[1]]

  # stop if none found
  if (length(x = meta_data_column) == 0) {
    stop("The 'meta_data_column': ", meta_data_column, " was not found in object meta.data slot.")
  }

  # Check that meta data is factor or character
  accepted_meta_types <- c("factor", "character", "logical")

  if (!class(x = seurat_object@meta.data[[meta_data_column]]) %in% accepted_meta_types) {
    stop("The 'meta_data_column': ", meta_data_column, " is of class: ", '"', class(x = seurat_object@meta.data[[meta_data_column]]), '"', " only meta data variables of classes: factor, character, or logical can be used with Meta_Highlight_Plot().")
  }

  # Change active ident for plotting
  Idents(seurat_object) <- meta_data_column

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path)) {
    if (!dir.exists(paths = file_path)) {
      stop("Provided `file_path`: ", '"', file_path, '"', " does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name)) {
    stop("No file name provided.  Please provide a file name using `file_name`.")
  }

  # Set file type for single pdf option
  if (single_pdf && is.null(x = file_type)) {
    file_type <- ".pdf"
  }
  if (single_pdf && !is.null(x = file_type) && str_detect(file_type, ".pdf") == FALSE) {
    message("WARNING: non-PDF 'file_type' specified but 'single_pdf = TRUE' selected.  Changing file_type to .pdf for output.")
    file_type <- ".pdf"
  }

  # Check file_type parameter
  file_type_options <- c(".pdf", ".png", ".tiff", ".jpeg", ".svg")
  if (is.null(x = file_type)) {
    stop("'file_type' not specified must be one of the following: '.pdf', '.png', '.tiff', '.jpeg', '.svg'")
  }
  if (!file_type %in% file_type_options) {
    stop("'file_type' must be one of the following: '.pdf', '.png', '.tiff', '.jpeg', '.svg'")
  }

  # Extract default reduction
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = colnames(x = seurat_object)) > 2e5)

  # Relevel idents for plotting to sorted order
  if (single_pdf && is.null(x = new_meta_order) && meta_data_sort) {
    Idents(seurat_object) <- fct_relevel(Idents(seurat_object), sort)
  }

  # Relevel idents to custom order
  if (single_pdf && !is.null(x = new_meta_order)) {
    if (length(x = new_meta_order) != length(x = levels(x = seurat_object@active.ident))) {
      stop("The length of 'new_meta_order' (", length(x = new_meta_order), ") does not equal the number of levels in ", meta_data_column, " (", length(x = levels(x = seurat_object@active.ident)), ").")
    }
    Idents(seurat_object) <- factor(Idents(seurat_object), levels = new_meta_order)
  }

  # Get number of clusters/identities
  list_idents <- levels(x = seurat_object@active.ident)
  num_idents <- length(x = list_idents)

  # Modify ident names if needed for saving plots
  list_idents_save <- gsub(pattern = "/", replacement = "-", x = list_idents)

  # Get color palette
  if (length(x = highlight_color) == 1) {
    highlight_color <- rep(highlight_color, num_idents)
  }
  if (length(x = highlight_color) > 1) {
    if (length(x = highlight_color) != num_idents) {
      stop("Number of colors provided to `highlight_color` (", length(highlight_color), ") is not equal to the number of clusters (", num_idents, ").")
    } else {
      highlight_color <- highlight_color
    }
  }

  # Single PDF option
  if (single_pdf == TRUE) {
    message("Generating plots")
    pboptions(char = "=")
    all_plots <- pbmapply(function(x, arg1) {
      cells_to_highlight <- CellsByIdentities(seurat_object, idents = x)
      suppressMessages(DimPlot(object = seurat_object,
                               cells.highlight = cells_to_highlight,
                               cols.highlight = arg1,
                               cols = background_color,
                               sizes.highlight = pt.size,
                               pt.size = pt.size,
                               order = TRUE,
                               reduction = reduction,
                               raster = raster,
                               ...))
    }, list_idents, highlight_color, SIMPLIFY = FALSE)
    message("Saving plots to file")
    pdf(paste(file_path, file_name, file_type, sep=""))
    pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
    for (i in 1:length(all_plots)) {
      print(all_plots[[i]])
      setTxtProgressBar(pb = pb, value = i)
    }
    close(con = pb)
    dev.off()
  }
  else {
    # Code for non-PDF figure
    if (str_detect(file_type, ".pdf") == FALSE) {
      message("Generating plots and saving plots to file")
      pb <- txtProgressBar(min = 0, max = num_idents, style = 3, file = stderr())
      cells_to_highlight <- CellsByIdentities(seurat_object)
      for (i in 1:length(cells_to_highlight)) {
        suppressMessages(DimPlot(object = seurat_object,
                                 cells.highlight = cells_to_highlight[[i]],
                                 cols.highlight = highlight_color[i],
                                 cols = background_color,
                                 sizes.highlight = pt.size,
                                 pt.size = pt.size,
                                 order = TRUE,
                                 reduction = reduction,
                                 raster = raster,
                                 ...))
        suppressMessages(ggsave(filename = paste(file_path, list_idents_save[i], "_", file_name, file_type, sep=""), dpi = dpi))
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
    # Code for PDF Version
    if (str_detect(file_type, ".pdf") == TRUE) {
      message("Generating plots and saving plots to file")
      pb <- txtProgressBar(min = 0, max = num_idents, style = 3, file = stderr())
      cells_to_highlight <- CellsByIdentities(seurat_object)
      for (i in 1:length(cells_to_highlight)) {
        suppressMessages(DimPlot(object = seurat_object,
                                 cells.highlight = cells_to_highlight[[i]],
                                 cols.highlight = highlight_color[i],
                                 cols = background_color,
                                 sizes.highlight = pt.size,
                                 pt.size = pt.size,
                                 order = TRUE,
                                 reduction = reduction,
                                 raster = raster,
                                 ...))
        suppressMessages(ggsave(filename = paste(file_path, list_idents_save[[i]], "_", file_name, file_type, sep=""), useDingbats = FALSE))
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
  }
}


#' Iterative Plotting of Gene Lists using Custom Featureplots
#'
#' Create and Save plots for Gene list with Single Command
#'
#' @param seurat_object Seurat object name.
#' @param gene_list vector of genes to plot.  If a named vector is provided then the names for each gene
#'  will be incorporated into plot title if `single_pdf = TRUE` or into file name if `FALSE`.
#' @param colors_use color scheme to use.
#' @param na_color color for non-expressed cells.
#' @param na_cutoff Value to use as minimum expression cutoff.  To set no cutoff set to `NA`.
#' @param split.by Variable in `@meta.data` to split the plot by.
#' @param order whether to move positive cells to the top (default = TRUE).
#' @param return_plots logical. Whether to return plots to list instead of saving them to file(s).  Default is FALSE.
#' @param file_path/prefix directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix and file extension.
#' @param file_type File type to save output as.  Must be one of following: ".pdf", ".png", ".tiff", ".jpeg", or ".svg".
#' @param single_pdf saves all plots to single PDF file (default = FALSE).  `file_type`` must be .pdf.
#' @param dpi dpi for image saving.
#' @param pt.size Adjust point size for plotting.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default).
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param alpha_exp new alpha level to apply to expressing cell color palette (`colors_use`).  Must be
#' value between 0-1.
#' @param alpha_na_exp new alpha level to apply to non-expressing cell color palette (`na_color`).  Must be
#' value between 0-1.
#' @param ... Extra parameters passed to \code{\link[Seurat]{FeaturePlot}}.
#'
#' @import ggplot2
#' @importFrom pbapply pblapply pboptions
#' @importFrom Seurat FeaturePlot
#' @importFrom SeuratObject DefaultDimReduc
#' @importFrom stringr str_detect
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @concept iterative_plotting
#'
#' @examples
#' \dontrun{
#' Iterate_FeaturePlot_scCustom(seurat_object = object, gene_list = DEG_list, colors_use = viridis_plasma_dark_high,
#' na_color = "lightgray", file_path = "plots/", file_name = "tsne", file_type = ".jpg", dpi = 600)
#' }
#'

Iterate_FeaturePlot_scCustom <- function(
  seurat_object,
  gene_list,
  colors_use = viridis_plasma_dark_high,
  na_color = "lightgray",
  na_cutoff = 0.000000001,
  split.by = NULL,
  order = TRUE,
  return_plots = FALSE,
  file_path = NULL,
  file_name = NULL,
  file_type = NULL,
  single_pdf = FALSE,
  dpi = 600,
  pt.size = NULL,
  reduction = NULL,
  raster = NULL,
  alpha_exp = NULL,
  alpha_na_exp = NULL,
  ...
) {
  # temp turn off message call from FeaturePlot_scCustomize
  op <- options(scCustomize_warn_na_cutoff = FALSE)
  on.exit(options(op))

  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = colnames(x = seurat_object)) > 2e5)

  # Return plot check
  if (return_plots) {
    if (!is.null(x = file_type) | !is.null(x = file_path) | !is.null(x = file_name) | single_pdf) {
      stop("Cannot return plots to list and save plots to file with single function call.  If saving plots please set 'return_plots = FALSE'.  If returning plots please leave 'file_type', 'file_path', 'file_name' and 'single_pdf' at their default settings.")
    }
  }

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path) && !return_plots) {
    if (!dir.exists(paths = file_path)) {
      stop("Provided `file_path`: ", '"', file_path, '"', " does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name) && !return_plots) {
    stop("No file name provided.  Please provide a file name using `file_name`.")
  }

  # Extract default reduction
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)

  # Set file type for single pdf option
  if (single_pdf && is.null(x = file_type)) {
    file_type <- ".pdf"
  }
  if (single_pdf && !is.null(x = file_type) && str_detect(file_type, ".pdf") == FALSE) {
    message("WARNING: non-PDF 'file_type' specified but 'single_pdf = TRUE' selected.  Changing file_type to .pdf for output.")
    file_type <- ".pdf"
  }

  # Check file_type parameter
  file_type_options <- c(".pdf", ".png", ".tiff", ".jpeg", ".svg")
  if (is.null(x = file_type) && !return_plots) {
    stop("'file_type' not specified must be one of the following: '.pdf', '.png', '.tiff', '.jpeg', '.svg'")
  }
  if (!file_type %in% file_type_options && !return_plots) {
    stop("'file_type' must be one of the following: '.pdf', '.png', '.tiff', '.jpeg', '.svg'")
  }

  # Check whether features are present in object
  gene_list <- Gene_Present(data = seurat_object, gene_list = gene_list, print_msg = FALSE, case_check = TRUE)[[1]]

  # Modify Cluster Labels names if needed for saving plots
  if (!is.null(x = names(gene_list)) && !single_pdf) {
    names_vec_mod <- gsub(pattern = "/", replacement = "-", x = names(x = gene_list))
    names(gene_list) <- names_vec_mod
  }

  # Return plots instead of saving them
  if (return_plots) {
    message("Generating plots")
    pboptions(char = "=")
    all_plots <- pblapply(gene_list,function(gene) {FeaturePlot_scCustom(seurat_object = seurat_object, features = gene, colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, split.by = split.by, order = order, pt.size = pt.size, reduction = reduction, raster = raster, alpha_exp = alpha_exp, alpha_na_exp = alpha_na_exp, ...)})
    return(all_plots)
  }

  # Single PDF option
  if (single_pdf == TRUE) {
    message("Generating plots")
    pboptions(char = "=")
    all_plots <- pblapply(gene_list,function(gene) {FeaturePlot_scCustom(seurat_object = seurat_object, features = gene, colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, split.by = split.by, order = order, pt.size = pt.size, reduction = reduction, raster = raster, alpha_exp = alpha_exp, alpha_na_exp = alpha_na_exp,...)})
    message("Saving plots to file")
    # save plots with cluster annotation
    if (!is.null(x = names(x = gene_list)) && is.null(x = split.by)) {
      pdf(paste(file_path, file_name, file_type, sep=""))
      pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
      for (i in 1:length(all_plots)) {
        print(all_plots[[i]] + ggtitle((paste0(gene_list[i], "_", names(x = gene_list)[i]))))
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
      dev.off()
    } else {
      # Save plots without cluster annotation
      pdf(paste(file_path, file_name, file_type, sep=""))
      pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
      for (i in 1:length(all_plots)) {
        print(all_plots[[i]])
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
      dev.off()
    }
  }
  else {
    if (str_detect(file_type, ".pdf") == FALSE) {
      message("Generating plots and saving plots to file")
      pb <- txtProgressBar(min = 0, max = length(gene_list), style = 3, file = stderr())
      for (i in 1:length(gene_list)) {
        FeaturePlot_scCustom(seurat_object = seurat_object, features = gene_list[i], colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, split.by = split.by, order = order, pt.size = pt.size, reduction = reduction, raster = raster, alpha_exp = alpha_exp, alpha_na_exp = alpha_na_exp, ...)
        if (!is.null(x = names(x = gene_list))) {
          suppressMessages(ggsave(filename = paste(file_path, gene_list[i], "_", names(x = gene_list)[i], "_", file_name, file_type, sep=""), dpi = dpi))
        } else {
          suppressMessages(ggsave(filename = paste(file_path, gene_list[i], "_", file_name, file_type, sep=""), dpi = dpi))
        }
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
    if (str_detect(file_type, ".pdf") == TRUE) {
      message("Generating plots and saving plots to file")
      pb <- txtProgressBar(min = 0, max = length(gene_list), style = 3, file = stderr())
      for (i in 1:length(gene_list)) {
        FeaturePlot_scCustom(seurat_object = seurat_object, features = gene_list[i], colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, split.by = split.by, order = order, pt.size = pt.size, reduction = reduction, raster = raster, alpha_exp = alpha_exp, alpha_na_exp = alpha_na_exp, ...)
        if (!is.null(x = names(x = gene_list))) {
          suppressMessages(ggsave(filename = paste(file_path, gene_list[i], "_", names(x = gene_list)[i], "_", file_name, file_type, sep=""), useDingbats = FALSE))
        } else {
          suppressMessages(ggsave(filename = paste(file_path, gene_list[i], "_", file_name, file_type, sep=""), useDingbats = FALSE))
        }
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
  }

  # One warning rastering
  if (!raster && single_pdf && getOption(x = 'scCustomize_warn_raster_iterative', default = TRUE)) {
    message(
      "NOTE: The raster parameter is currently set to FALSE and `single_pdf = TRUE`.\n",
      "Saving large numbers of plots in vector form can result in very large\n",
      "file sizes. Suggest setting `raster = TRUE` when plotting large numbers \n",
      "of features in single output file.
      \nThis message will be shown once per session.\n"
    )
    options(scCustomize_warn_raster_iterative = FALSE)
  }
}


#' Iterative Plotting of Gene Lists using VlnPlot_scCustom
#'
#' Create and Save plots for Gene list with Single Command
#'
#' @param seurat_object Seurat object name.
#' @param gene_list list of genes to plot.
#' @param colors_use color palette to use for plotting.  By default if number of levels plotted is less than
#' or equal to 36 it will use "polychrome" and if greater than 36 will use "varibow" with shuffle = TRUE
#' both from `DiscretePalette_scCustomize`.
#' @param pt.size point size for plotting.
#' @param group.by Name of one or more metadata columns to group (color) plot by (for example, orig.ident);
#' default is the current active.ident of the object.
#' @param split.by Feature to split plots by (i.e. "orig.ident").
#' @param file_path/prefix directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix and file extension.
#' @param file_type File type to save output as.  Must be one of following: ".pdf", ".png", ".tiff", ".jpeg", or ".svg".
#' @param single_pdf saves all plots to single PDF file (default = FALSE).  `file_type`` must be .pdf.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 100,000 total points plotted (# Cells x # of features).
#' @param dpi dpi for image saving.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param ... Extra parameters passed to \code{\link[Seurat]{VlnPlot}}.
#'
#' @import ggplot2
#' @importFrom pbapply pblapply pboptions
#' @importFrom Seurat VlnPlot
#' @importFrom stringr str_detect
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @concept iterative_plotting
#'
#' @examples
#' \dontrun{
#' Iterate_VlnPlot_scCustom(seurat_object = object, gene_list = DEG_list, colors = color_list, file_path = "plots/",
#' file_name = "_vln", file_type = ".jpg", dpi = 600)
#' }
#'

Iterate_VlnPlot_scCustom <- function(
  seurat_object,
  gene_list,
  colors_use = NULL,
  pt.size = NULL,
  group.by = NULL,
  split.by = NULL,
  file_path = NULL,
  file_name = NULL,
  file_type = NULL,
  single_pdf = FALSE,
  raster = NULL,
  dpi = 600,
  ggplot_default_colors = FALSE,
  color_seed = 123,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Add pt.size check
  pt.size <- pt.size %||% AutoPointSize_scCustom(data = seurat_object)

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path)) {
    if (!dir.exists(paths = file_path)) {
      stop("Provided `file_path`: ", '"', file_path, '"', " does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name)) {
    stop("No file name provided.  Please provide a file name using `file_name`.")
  }

  # Set file type for single pdf option
  if (single_pdf && is.null(x = file_type)) {
    file_type <- ".pdf"
  }
  if (single_pdf && !is.null(x = file_type) && str_detect(file_type, ".pdf") == FALSE) {
    message("WARNING: non-PDF 'file_type' specified but 'single_pdf = TRUE' selected.  Changing file_type to .pdf for output.")
    file_type <- ".pdf"
  }

  # Check file_type parameter
  file_type_options <- c(".pdf", ".png", ".tiff", ".jpeg", ".svg")
  if (is.null(x = file_type)) {
    stop("'file_type' not specified must be one of the following: '.pdf', '.png', '.tiff', '.jpeg', '.svg'")
  }
  if (!file_type %in% file_type_options) {
    stop("'file_type' must be one of the following: '.pdf', '.png', '.tiff', '.jpeg', '.svg'")
  }

  # check for single pdf file type
  if (single_pdf && str_detect(file_type, ".pdf") == FALSE) {
    stop("File type set to non-PDF type but 'single_pdf = TRUE' also selected.  Confirm output desired and re-run function.")
  }

  # Check whether features are present in object
  gene_list <- Gene_Present(data = seurat_object, gene_list = gene_list, print_msg = FALSE, case_check = TRUE)[[1]]

  # Set default color palette based on number of levels being plotted
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }

  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use) && ggplot_default_colors) {
    stop("Cannot provide both custom palette to `colors_use` and specify `ggplot_default_colors = TRUE`.")
  }
  if (is.null(x = colors_use)) {
    # set default plot colors
    if (is.null(x = colors_use)) {
      colors_use <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
    }
  }

  # Add one time raster warning
  if (single_pdf && pt.size != 0 && getOption(x = 'scCustomize_warn_vln_raster_iterative', default = TRUE)) {
    message(
      "NOTE: `single_pdf = TRUE` and `pt.size` > 0, so all points are plotted.\n",
      "Seurat::VlnPlot does not currently support raster plotting and therefore\n",
      "Saving large numbers of plots in vector form can result in very large\n",
      "file sizes. Suggest setting `pt.size = 0` or splitting in small batches\n",
      "when plotting large numbers of features in single output file.
      \nThis message will be shown once per session.\n"
    )
    options(scCustomize_warn_vln_raster_iterative = FALSE)
  }

  # Single PDF option
  if (single_pdf == TRUE) {
    message("Generating plots")
    pboptions(char = "=")
    all_plots <- pblapply(gene_list,function(gene) {VlnPlot_scCustom(seurat_object = seurat_object, features = gene, colors_use = colors_use, pt.size = pt.size, group.by = group.by, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, split.by = split.by, ...)})
    message("Saving plots to file")
    pdf(paste(file_path, file_name, file_type, sep=""))
    pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
    for (i in 1:length(all_plots)) {
      print(all_plots[[i]])
      setTxtProgressBar(pb = pb, value = i)
    }
    close(con = pb)
    dev.off()
  }
  else {
    if (str_detect(file_type, ".pdf") == FALSE) {
      message("Generating plots and saving plots to file")
      pb <- txtProgressBar(min = 0, max = length(gene_list), style = 3, file = stderr())
      for (i in 1:length(gene_list)) {
        VlnPlot_scCustom(seurat_object = seurat_object, features = gene_list[i], colors_use = colors_use, pt.size = pt.size, group.by = group.by, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, split.by = split.by, ...)
        suppressMessages(ggsave(filename = paste(file_path, gene_list[i], file_name, file_type, sep=""), dpi = dpi))
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
    if (str_detect(file_type, ".pdf") == TRUE) {
      message("Generating plots and saving plots to file")
      pb <- txtProgressBar(min = 0, max = length(gene_list), style = 3, file = stderr())
      for (i in 1:length(gene_list)) {
        VlnPlot_scCustom(seurat_object = seurat_object, features = gene_list[i], colors_use = colors_use, pt.size = pt.size, group.by = group.by, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, split.by = split.by, ...)
        suppressMessages(ggsave(filename = paste(file_path, gene_list[i], file_name, file_type, sep=""), useDingbats = FALSE))
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
  }
}


#' Iterative Plotting of Gene Lists using Custom Density Plots
#'
#' Create and save plots for gene list with single command.  Requires Nebulosa package from Bioconductor.
#'
#' @param seurat_object Seurat object name.
#' @param gene_list vector of genes to plot.  If a named vector is provided then the names for each gene
#' will be incorporated into plot title if `single_pdf = TRUE` or into file name if `FALSE`.
#' @param viridis_palette color scheme to use.
#' @param custom_palette color for non-expressed cells.
#' #' @param pt.size Adjust point size for plotting.
#' @param file_path/prefix directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix and file extension.
#' @param file_type File type to save output as.  Must be one of following: ".pdf", ".png", ".tiff", ".jpeg", or ".svg".
#' @param single_pdf saves all plots to single PDF file (default = FALSE).  `file_type`` must be .pdf.
#' @param dpi dpi for image saving.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default)
#' #' @param joint NULL.  This function only supports `joint = FALSE`.  Leave as NULL to generate plots.  To iterate joint plots see function: `Iterate_Plot_Density_Joint`.
#' @param combine Create a single plot? If FALSE, a list with ggplot objects is returned.
#' @param ... Extra parameters passed to \code{\link[Nebulosa]{plot_density}}.
#'
#' @import ggplot2
#' @importFrom pbapply pblapply pboptions
#' @importFrom SeuratObject DefaultDimReduc
#' @importFrom stringr str_detect
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @concept iterative_plotting
#'
#' @examples
#' \dontrun{
#' Iterate_Plot_Density_Custom(seurat_object = object, gene_list = DEG_list, viridis_palette = "magma",
#' file_path = "plots/", file_name = "_density_plots", file_type = ".jpg", dpi = 600)
#' }
#'

Iterate_Plot_Density_Custom <- function(
  seurat_object,
  gene_list,
  viridis_palette = "magma",
  custom_palette = NULL,
  pt.size = 1,
  file_path = NULL,
  file_name = NULL,
  file_type = NULL,
  single_pdf = FALSE,
  dpi = 600,
  reduction = NULL,
  combine = TRUE,
  joint = FALSE,
  ...
) {
  # Check Nebulosa installed
  Nebulosa_check <- PackageCheck("Nebulosa", error = FALSE)
  if (!Nebulosa_check[1]) {
    stop(
      "Please install the Nebulosa package to use 'Iterate_Plot_Density_Custom'",
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

  # joint check
  if (!is.null(x = joint)) {
    stop("`Iterate_Plot_Density_Custom` only supports 'joint = FALSE'.  Leave as NULL to generate plots.  To iterate joint plots see function: 'Iterate_Plot_Density_Joint'.")
  }

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path)) {
    if (!dir.exists(paths = file_path)) {
      stop("Provided `file_path`: ", '"', file_path, '"', " does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name)) {
    stop("No file name provided.  Please provide a file name using `file_name`.")
  }

  # Extract default reduction
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)

  # Set file type for single pdf option
  if (single_pdf && is.null(x = file_type)) {
    file_type <- ".pdf"
  }
  if (single_pdf && !is.null(x = file_type) && str_detect(file_type, ".pdf") == FALSE) {
    message("WARNING: non-PDF 'file_type' specified but 'single_pdf = TRUE' selected.  Changing file_type to .pdf for output.")
    file_type <- ".pdf"
  }

  # Check file_type parameter
  file_type_options <- c(".pdf", ".png", ".tiff", ".jpeg", ".svg")
  if (is.null(x = file_type)) {
    stop("'file_type' not specified must be one of the following: '.pdf', '.png', '.tiff', '.jpeg', '.svg'")
  }
  if (!file_type %in% file_type_options) {
    stop("'file_type' must be one of the following: '.pdf', '.png', '.tiff', '.jpeg', '.svg'")
  }

  # Check whether features are present in object
  gene_list <- Gene_Present(data = seurat_object, gene_list = gene_list, print_msg = FALSE, case_check = TRUE)[[1]]

  # check palettes
  if (!is.null(x = custom_palette) && viridis_palette != "magma") {
    stop("Non-default values provided to both viridis_palette & custom_palette.  Please chose one non-default value.")
  }

  # Modify Cluster Labels names if needed for saving plots
  if (!is.null(x = names(gene_list)) && !single_pdf) {
    names_vec_mod <- gsub(pattern = "/", replacement = "-", x = names(x = gene_list))
    names(gene_list) <- names_vec_mod
  }

  # Single PDF option
  if (single_pdf == TRUE) {
    message("Generating plots")
    pboptions(char = "=")
    all_plots <- pblapply(gene_list,function(gene) {
      Plot_Density_Custom(seurat_object = seurat_object, features = gene, joint = FALSE, viridis_palette = viridis_palette, custom_palette = custom_palette, pt.size = pt.size, reduction = reduction, ...)})
    message("Saving plots to file")
    # save plots with cluster annotation
    if (!is.null(x = names(x = gene_list))) {
      pdf(paste(file_path, file_name, file_type, sep=""))
      pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
      for (i in 1:length(all_plots)) {
        print(all_plots[[i]] + ggtitle((paste0(gene_list[i], "_", names(x = gene_list)[i]))))
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
      dev.off()
    } else {
      # Save plots without cluster annotation
      pdf(paste(file_path, file_name, file_type, sep=""))
      pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
      for (i in 1:length(all_plots)) {
        print(all_plots[[i]])
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
      dev.off()
    }
  }
  else {
    if (str_detect(file_type, ".pdf") == FALSE) {
      message("Generating plots and saving plots to file")
      pb <- txtProgressBar(min = 0, max = length(gene_list), style = 3, file = stderr())
      for (i in 1:length(gene_list)) {
        Plot_Density_Custom(seurat_object = seurat_object, features = gene_list[i], joint = FALSE, viridis_palette = viridis_palette, custom_palette = custom_palette, pt.size = pt.size, reduction = reduction, ...)
        if (!is.null(x = names(x = gene_list))) {
          suppressMessages(ggsave(filename = paste(file_path, gene_list[i], "_", names(x = gene_list)[i], "_", file_name, file_type, sep=""), dpi = dpi))
        } else {
          suppressMessages(ggsave(filename = paste(file_path, gene_list[i], "_", file_name, file_type, sep=""), dpi = dpi))
        }
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
    if (str_detect(file_type, ".pdf") == TRUE) {
      message("Generating plots and saving plots to file")
      pb <- txtProgressBar(min = 0, max = length(gene_list), style = 3, file = stderr())
      for (i in 1:length(gene_list)) {
        Plot_Density_Custom(seurat_object = seurat_object, features = gene_list[i], joint = FALSE, viridis_palette = viridis_palette, custom_palette = custom_palette, pt.size = pt.size, reduction = reduction, ...)
        if (!is.null(x = names(x = gene_list))) {
          suppressMessages(ggsave(filename = paste(file_path, gene_list[i], "_", names(x = gene_list)[i], "_", file_name, file_type, sep=""), useDingbats = FALSE))
        } else {
          suppressMessages(ggsave(filename = paste(file_path, gene_list[i], "_", file_name, file_type, sep=""), useDingbats = FALSE))
        }
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
  }
}


#' Iterative Plotting of Gene Lists using Custom Joint Density Plots
#'
#' Create and save plots for gene list with single command.  Requires Nebulosa package from Bioconductor.
#'
#' @param seurat_object Seurat object name.
#' @param gene_list a list of vectors of genes to plot jointly.  Each entry in the list will be plotted
#' for the joint density.  All entries in list must be greater than 2 features.  If a named list is provided
#' then the names for each list entry will be incorporated into plot title if `single_pdf = TRUE` or
#' into file name if `FALSE`.
#' @param viridis_palette color scheme to use.
#' @param custom_palette color for non-expressed cells.
#' #' @param pt.size Adjust point size for plotting.
#' @param file_path/prefix directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix and file extension.
#' @param file_type File type to save output as.  Must be one of following: ".pdf", ".png", ".tiff", ".jpeg", or ".svg".
#' @param single_pdf saves all plots to single PDF file (default = FALSE).  `file_type`` must be .pdf.
#' @param dpi dpi for image saving.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default)
#' #' @param joint NULL.  This function only supports `joint = FALSE`.  Leave as NULL to generate plots.  To iterate joint plots see function: `Iterate_Plot_Density_Joint`.
#' @param combine Create a single plot? If FALSE, a list with ggplot objects is returned.
#' @param ... Extra parameters passed to \code{\link[Nebulosa]{plot_density}}.
#'
#' @import ggplot2
#' @importFrom pbapply pblapply pboptions
#' @importFrom purrr discard keep
#' @importFrom SeuratObject DefaultDimReduc
#' @importFrom stringr str_detect
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @concept iterative_plotting
#'
#' @examples
#' \dontrun{
#' Iterate_Plot_Density_Joint(seurat_object = object, gene_list = DEG_list, viridis_palette = "magma",
#' file_path = "plots/", file_name = "joint_plots", file_type = ".jpg", dpi = 600)
#' }
#'

Iterate_Plot_Density_Joint <- function(
  seurat_object,
  gene_list,
  viridis_palette = "magma",
  custom_palette = NULL,
  pt.size = 1,
  file_path = NULL,
  file_name = NULL,
  file_type = NULL,
  single_pdf = FALSE,
  dpi = 600,
  reduction = NULL,
  combine = TRUE,
  joint = NULL,
  ...
) {
  # Check Nebulosa installed
  Nebulosa_check <- PackageCheck("Nebulosa", error = FALSE)
  if (!Nebulosa_check[1]) {
    stop(
      "Please install the Nebulosa package to use 'Iterate_Plot_Density_Joint'",
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

  # Check gene list is in list form
  if (class(x = gene_list) != "list") {
    stop("For 'Iterate_Plot_Density_Joint' the 'gene_list' must be of class() 'list', please reformat from current class(): ", '"', class(x = gene_list), '"', ".")
  }

  # joint check
  if (!is.null(x = joint)) {
    stop("`Iterate_Plot_Density_Joint` only supports 'joint = TRUE'.  Leave as NULL to generate plots.  To iterate individual plots see function: 'Iterate_Plot_Density_Custom'.")
  }

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path)) {
    if (!dir.exists(paths = file_path)) {
      stop("Provided `file_path`: ", '"', file_path, '"', " does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name)) {
    stop("No file name provided.  Please provide a file name using `file_name`.")
  }

  # Extract default reduction
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)

  # Set file type for single pdf option
  if (single_pdf && is.null(x = file_type)) {
    file_type <- ".pdf"
  }
  if (single_pdf && !is.null(x = file_type) && str_detect(file_type, ".pdf") == FALSE) {
    message("WARNING: non-PDF 'file_type' specified but 'single_pdf = TRUE' selected.  Changing file_type to .pdf for output.")
    file_type <- ".pdf"
  }

  # Check file_type parameter
  file_type_options <- c(".pdf", ".png", ".tiff", ".jpeg", ".svg")
  if (is.null(x = file_type)) {
    stop("'file_type' not specified must be one of the following: '.pdf', '.png', '.tiff', '.jpeg', '.svg'")
  }
  if (!file_type %in% file_type_options) {
    stop("'file_type' must be one of the following: '.pdf', '.png', '.tiff', '.jpeg', '.svg'")
  }

  # Check whether features are present in object
  checked_gene_list <- lapply(1:length(gene_list), function(x){
    genes <- Gene_Present(data = seurat_object, gene_list = gene_list[[x]], print_msg = FALSE, case_check = TRUE, return_none = TRUE)[[1]]
  })

  if (!is.null(x = names(x = gene_list))) {
    names(checked_gene_list) <- names(x = gene_list)
  }

  # remove any empty entries in list
  checked_gene_list <- discard(checked_gene_list,  ~length(.x) == 0)

  # Check for lists less than 2
  bad_gene_lists <- keep(checked_gene_list,  ~length(.x) < 2)
  if (length(x = bad_gene_lists) > 0) {
    warning("A total of ", length(x = bad_gene_lists), " list entries from `gene_list contain less than two features and were excluded from plotting.")
  }

  # Create final good gene list
  final_gene_list <- keep(checked_gene_list,  ~length(.x) > 1)

  # check palettes
  if (!is.null(x = custom_palette) && viridis_palette != "magma") {
    stop("Non-default values provided to both viridis_palette & custom_palette.  Please chose one non-default value.")
  }

  # Modify Cluster Labels names if needed for saving plots
  if (!is.null(x = names(gene_list)) && !single_pdf) {
    names_vec_mod <- gsub(pattern = "/", replacement = "-", x = names(x = gene_list))
    names(gene_list) <- names_vec_mod
  }

  # Single PDF option
  if (single_pdf == TRUE) {
    message("Generating plots")
    pboptions(char = "=")
    all_plots <- pblapply(1:length(final_gene_list),function(i) {
      plot <- Plot_Density_Joint_Only(seurat_object = seurat_object, features = final_gene_list[[i]], viridis_palette = viridis_palette, custom_palette = custom_palette, pt.size = pt.size, reduction = reduction, ...)})
    message("Saving plots to file")
    # save plots with cluster annotation
    if (!is.null(x = names(x = final_gene_list))) {
      pdf(paste(file_path, file_name, file_type, sep=""))
      pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
      for (i in 1:length(all_plots)) {
        print(all_plots[[i]] + ggtitle((paste0(paste(final_gene_list[[i]], collapse = "_"), "_", names(x = final_gene_list)[i]))))
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
      dev.off()
    } else {
      # Save plots without cluster annotation
      pdf(paste(file_path, file_name, file_type, sep=""))
      pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
      for (i in 1:length(all_plots)) {
        print(all_plots[[i]])
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
      dev.off()
    }
  }
  else {
    if (str_detect(file_type, ".pdf") == FALSE) {
      message("Generating plots and saving plots to file")
      pb <- txtProgressBar(min = 0, max = length(final_gene_list), style = 3, file = stderr())
      for (i in 1:length(final_gene_list)) {
        Plot_Density_Joint_Only(seurat_object = seurat_object, features = final_gene_list[[i]], viridis_palette = viridis_palette, custom_palette = custom_palette, pt.size = pt.size, reduction = reduction, ...)
        if (!is.null(x = names(x = final_gene_list))) {
          suppressMessages(ggsave(filename = paste(file_path, paste(final_gene_list[[i]], collapse = "_"), "_", names(x = final_gene_list)[i], "_", file_name, file_type, sep=""), dpi = dpi))
        } else {
          suppressMessages(ggsave(filename = paste(file_path, paste(final_gene_list[[i]], collapse = "_"), "_", file_name, file_type, sep=""), dpi = dpi))
        }
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
    if (str_detect(file_type, ".pdf") == TRUE) {
      message("Generating plots and saving plots to file")
      pb <- txtProgressBar(min = 0, max = length(final_gene_list), style = 3, file = stderr())
      for (i in 1:length(final_gene_list)) {
        Plot_Density_Joint_Only(seurat_object = seurat_object, features = final_gene_list[[i]], viridis_palette = viridis_palette, custom_palette = custom_palette, pt.size = pt.size, reduction = reduction, ...)
        if (!is.null(x = names(x = final_gene_list))) {
          suppressMessages(ggsave(filename = paste(file_path, paste(final_gene_list[[i]], collapse = "_"), "_", names(x = final_gene_list)[i], "_", file_name, file_type, sep=""), useDingbats = FALSE))
        } else {
          suppressMessages(ggsave(filename = paste(file_path, paste(final_gene_list[[i]], collapse = "_"), "_", file_name, file_type, sep=""), useDingbats = FALSE))
        }
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
  }
}
