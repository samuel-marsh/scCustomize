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
#' @import cli
#' @import patchwork
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
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
      cli_abort(message = "Provided {.code file_path}: {symbol$dquote_left}{.field {file_path}}{symbol$dquote_right} does not exist.")
    }
  }

  # Set dims to plot if not specified
  num_pc_present <- length(seurat_object@reductions$pca@stdev)
  if (is.null(x = dims_plot)) {
    dims_plot <- num_pc_present
  }

  # Check pca present
  reduc_present <- names(x = seurat_object@reductions)
  if (!"pca" %in% reduc_present) {
    cli_abort(message = "Cannot find reduction 'pca' in this Seurat Object.")
  }
  # Check dims present in object
  if (dims_plot > num_pc_present) {
    cli_abort(message = "The number of PCs specified to {.code dims_plot} ({.field {dims_plot}}) is greater than number of PCs present in Seurat Object ({.field {num_pc_present}})")
  }

  dims_list <- 1:dims_plot
  # Create list of all plots
  cli_inform(message = "Generating plots")
  pboptions(char = "=")
  all_plots <- pblapply(dims_list, function(x) {
    PC_Plotting(seurat_object = seurat_object, dim_number = x)
  })
  # Save plots
  cli_inform(message = "Saving plots to file")
  pdf(file = paste(file_path, name_prefix, file_name, ".pdf", sep=""), height = 11, width = 8.5)
  pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
  for (i in 1:length(all_plots)) {
    print(all_plots[[i]])
    setTxtProgressBar(pb = pb, value = i)
  }
  close(con = pb)
  dev.off()
  if (isTRUE(x = return_plots)) {
    return(all_plots)
  }
}


#' Iterate DimPlot By Sample
#'
#' Iterate DimPlot by orig.ident column from Seurat object metadata
#'
#' @param seurat_object Seurat object name.
#' @param sample_column name of meta.data column containing sample names/ids (default is "orig.ident").
#' @param file_path directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix to append after sample name.
#' @param file_type File type to save output as.  Must be one of following: ".pdf", ".png", ".tiff", ".jpeg", or ".svg".
#' @param single_pdf saves all plots to single PDF file (default = FALSE).  `file_type`` must be .pdf
#' @param color color scheme to use.
#' @param no_legend logical, whether or not to include plot legend, default is TRUE.
#' @param title_prefix Value that should be used for plot title prefix if `no_legend = TRUE`.
#' If NULL the value of `meta_data_column` will be used.  Default is NULL.
#' @param title_prefix Value that should be used for plot title prefix if `no_legend = TRUE`.
#' If NULL the value of `meta_data_column` will be used.  Default is NULL.
#' @param dpi dpi for image saving.
#' @param reduction Dimensionality Reduction to use (default is object default).
#' @param dims Dimensions to plot.
#' @param pt.size Adjust point size for plotting.
#' @param raster Convert points to raster format.  Default is NULL which will rasterize by default if
#' greater than 200,000 cells.
#' @param ... Extra parameters passed to \code{\link[Seurat]{DimPlot}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
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
  sample_column = "orig.ident",
  file_path = NULL,
  file_name = NULL,
  file_type = NULL,
  single_pdf = FALSE,
  output_width = NULL,
  output_height = NULL,
  dpi = 600,
  color = "black",
  no_legend = TRUE,
  title_prefix = NULL,
  reduction = NULL,
  dims = c(1, 2),
  pt.size = NULL,
  raster = NULL,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Harmonize pt.size across all plots
  pt.size <- pt.size %||% AutoPointSize_scCustom(data = seurat_object)

  # Check meta.data column if not orig.ident
  if (sample_column != "orig.ident") {
    # Check meta data
    sample_column <- Meta_Present(object = seurat_object, meta_col_names = sample_column, omit_warn = FALSE, print_msg = FALSE)[[1]]

    # stop if none found
    if (length(x = sample_column) == 0) {
      cli_abort(message = c("No meta.data column found.",
                            "i" = "Column {.field {sample_column}} was not found in the meta.data slot."))
    }
  }

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path)) {
    if (!dir.exists(paths = file_path)) {
      cli_abort(message = "Provided {.code file_path}: {symbol$dquote_left}{.field {file_path}}{symbol$dquote_right} does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name)) {
    cli_abort(message = "No file name provided.  Please provide a file name using {.code file_name}.")
  }

  # Set file type for single pdf option
  if (isTRUE(x = single_pdf) && is.null(x = file_type)) {
    file_type <- ".pdf"
  }
  if (isTRUE(x = single_pdf) && !is.null(x = file_type) && str_detect(file_type, ".pdf") == FALSE) {
    cli_inform(message = "WARNING: non-PDF {.code file_type} specified but {.code single_pdf = TRUE} selected.  Changing file_type to {.val .pdf} for output.")
    file_type <- ".pdf"
  }

  # Check file_type parameter
  file_type_options <- c(".pdf", ".png", ".tiff", ".jpeg", ".svg")
  if (is.null(x = file_type)) {
    cli_abort(message = c("{.code file_type} not specified.",
                          "*" = "Must specify output file type format from the following:",
                          "i" = "{.field {glue_collapse_scCustom(input_string = file_type_options, and = TRUE)}}"))
  }
  if (!file_type %in% file_type_options) {
    cli_abort(message = "{.code file_type} must be one of the following: {.field {glue_collapse_scCustom(input_string = file_type_options, and = TRUE)}}")
  }

  # Extract reduction coordinates
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)
  cells <- Cells(x = seurat_object)
  reduc_coordinates <- Embeddings(object = seurat_object[[reduction]])[cells, dims]
  reduc_coordinates <- as.data.frame(x = reduc_coordinates)
  x_axis <- c(min(reduc_coordinates[, 1]),
              max(reduc_coordinates[, 1]))
  y_axis <- c(min(reduc_coordinates[, 2]),
              max(reduc_coordinates[, 2]))

  # Extract sample id column
  column_list <- as.character(x = unique(x = seurat_object@meta.data[[sample_column]]))
  num_idents <- length(x = column_list)

  # Create plot titles if needed.
  if (!is.null(x = title_prefix) && isFALSE(x = no_legend)) {
    cli_warn(message = "{.code title_prefix} was omitted as {.code no_legend = FALSE}.")
  }

  if (is.null(x = title_prefix) && isTRUE(x = no_legend)) {
    plot_title <- lapply(1:num_idents, function(z) {
      paste0(sample_column, ": ", column_list[z])
    })
  } else {
    plot_title <- lapply(1:num_idents, function(z) {
      paste0(title_prefix, ": ", column_list[z])
    })
  }

  if (!is.null(x = title_prefix) && length(x = title_prefix) != 1 && isTRUE(x = no_legend)) {
    cli_abort(message = "{.field `title_prefix`} must be vector of length 1.")
  }

  # Create list of cells per sample
  cells_per_sample <- lapply(column_list, function(sample) {
    row.names(x = seurat_object@meta.data)[which(x = seurat_object@meta.data[[sample_column]] == sample)]
  })

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)


  # Single PDF option
  if (isTRUE(x = single_pdf)) {
    cli_inform(message = "{.field Generating plots}")
    pboptions(char = "=")
    all_plots <- pblapply(1:num_idents,function(x) {
      if (isTRUE(x = no_legend)) {
        DimPlot(object = seurat_object, cells = cells_per_sample[[x]], group.by = sample_column, cols = color, reduction = reduction, pt.size = pt.size, raster = raster, ...) +
          xlim(x_axis) +
          ylim(y_axis) +
          NoLegend() +
          ggtitle(plot_title[x]) +
          CenterTitle()
      } else {
        DimPlot(object = seurat_object, cells = cells_per_sample[[x]], group.by = sample_column, cols = color, reduction = reduction, pt.size = pt.size, raster = raster, ...) +
          xlim(x_axis) +
          ylim(y_axis)
      }
      })
    cli_inform(message = "{.field Saving plots to file}")
    pdf(paste(file_path, file_name, file_type, sep=""), width = output_width, height = output_height)
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
      cli_inform(message = "{.field Generating plots and saving plots to file}")
      pb <- txtProgressBar(min = 0, max = length(cells_per_sample), style = 3, file = stderr())
      for (i in 1:length(cells_per_sample)) {
        if (isTRUE(x = no_legend)) {
          DimPlot(object = seurat_object, cells = cells_per_sample[[i]], group.by = sample_column, cols = color, reduction = reduction, pt.size = pt.size, raster = raster, ...) +
            xlim(x_axis) +
            ylim(y_axis) +
            NoLegend() +
            ggtitle(plot_title[i]) +
            CenterTitle()
        } else {
          DimPlot(object = seurat_object, cells = cells_per_sample[[i]], group.by = sample_column, cols = color, reduction = reduction, pt.size = pt.size, raster = raster, ...) +
            xlim(x_axis) +
            ylim(y_axis)
        }
        suppressMessages(ggsave(filename = paste(file_path, column_list[[i]], file_name, file_type, sep=""), dpi = dpi))
        setTxtProgressBar(pb = pb, value = i)
        }
      close(con = pb)
      }
    # Code for PDF Version
    if (str_detect(file_type, ".pdf") == TRUE) {
      cli_inform(message = "{.field Generating plots and saving plots to file}")
      pb <- txtProgressBar(min = 0, max = length(cells_per_sample), style = 3, file = stderr())
      for (i in 1:length(cells_per_sample)) {
        if (isTRUE(x = no_legend)) {
          DimPlot(object = seurat_object, cells = cells_per_sample[[i]], group.by = sample_column, cols = color, reduction = reduction, pt.size = pt.size, raster = raster, ...) +
            xlim(x_axis) +
            ylim(y_axis) +
            NoLegend() +
            ggtitle(plot_title[i]) +
            CenterTitle()
        } else {
          DimPlot(object = seurat_object, cells = cells_per_sample[[i]], group.by = sample_column, cols = color, reduction = reduction, pt.size = pt.size, raster = raster, ...) +
            xlim(x_axis) +
            ylim(y_axis)
        }
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
#' @import cli
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
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
#' Iterate_Cluster_Highlight_Plot(seurat_object = object, highlight_color = "navy",
#' background_color = "lightgray", file_path = "path/", file_name = "name", file_type = "pdf",
#' single_pdf = TRUE)
#' }
#'

Iterate_Cluster_Highlight_Plot <- function(
    seurat_object,
    highlight_color = "dodgerblue",
    background_color = "lightgray",
    pt.size = NULL,
    reduction = NULL,
    file_path = NULL,
    file_name = NULL,
    file_type = NULL,
    single_pdf = FALSE,
    output_width = NULL,
    output_height = NULL,
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
      cli_abort(message = "Provided {.code file_path}: {symbol$dquote_left}{.field {file_path}}{symbol$dquote_right} does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name)) {
    cli_abort(message = "No file name provided.  Please provide a file name using {.code file_name}.")
  }

  # Set file type for single pdf option
  if (isTRUE(x = single_pdf) && is.null(x = file_type)) {
    file_type <- ".pdf"
  }
  if (isTRUE(x = single_pdf) && !is.null(x = file_type) && str_detect(file_type, ".pdf") == FALSE) {
    cli_inform(message = "WARNING: non-PDF {.code file_type} specified but {.code single_pdf = TRUE} selected.  Changing file_type to {.val .pdf} for output.")
    file_type <- ".pdf"
  }

  # Check file_type parameter
  file_type_options <- c(".pdf", ".png", ".tiff", ".jpeg", ".svg")
  if (is.null(x = file_type)) {
    cli_abort(message = c("{.code file_type} not specified.",
                          "*" = "Must specify output file type format from the following:",
                          "i" = "{.field {glue_collapse_scCustom(input_string = file_type_options, and = TRUE)}}"))
  }

  if (!file_type %in% file_type_options) {
    cli_abort(message = "{.code file_type} must be one of the following: {.field {glue_collapse_scCustom(input_string = file_type_options, and = TRUE)}}")
  }

  # Extract default reduction
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)

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
      cli_abort(message = "Number of colors provided to {.code highlight_color} ({.field {length(highlight_color)}}) is not equal to the number of clusters (({.field {num_idents}}).")
    } else {
      highlight_color <- highlight_color
    }
  }

  # Single PDF option
  if (isTRUE(x = single_pdf)) {
    cli_inform(message = "{.field Generating plots}")
    pboptions(char = "=")
    all_plots <- pblapply(1:num_idents, function(x) {
      suppressMessages(Cluster_Highlight_Plot(seurat_object = seurat_object,
                                              cluster_name = list_idents[x],
                                              highlight_color = highlight_color[x],
                                              background_color = background_color,
                                              pt.size = pt.size,
                                              reduction = reduction,
                                              raster = raster,
                                              ...))
    })
    cli_inform(message = "{.field Saving plots to file}")
    pdf(paste(file_path, file_name, file_type, sep=""), width = output_width, height = output_height)
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
      cli_inform(message = "{.field Generating plots and saving plots to file}")
      pb <- txtProgressBar(min = 0, max = num_idents, style = 3, file = stderr())
      for (i in 1:num_idents) {
        suppressMessages(Cluster_Highlight_Plot(seurat_object = seurat_object,
                                                cluster_name = list_idents[i],
                                                highlight_color = highlight_color[i],
                                                background_color = background_color,
                                                pt.size = pt.size,
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
      cli_inform(message = "{.field Generating plots and saving plots to file}")
      pb <- txtProgressBar(min = 0, max = num_idents, style = 3, file = stderr())
      for (i in 1:num_idents) {
        suppressMessages(Cluster_Highlight_Plot(seurat_object = seurat_object,
                                                cluster_name = list_idents[i],
                                                highlight_color = highlight_color[i],
                                                background_color = background_color,
                                                pt.size = pt.size,
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
#' @param no_legend logical, whether or not to remove plot legend and move to plot title.  Default is FALSE.
#' @param title_prefix Value that should be used for plot title prefix if `no_legend = TRUE`.
#' If NULL the value of `meta_data_column` will be used.  Default is NULL.
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
#' @import cli
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom forcats fct_relevel
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
#' Iterate_Meta_Highlight_Plot(seurat_object = object, meta_data_column = "sample_id",
#' highlight_color = "navy", background_color = "lightgray", file_path = "path/",
#' file_name = "name", file_type = "pdf", single_pdf = TRUE)
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
  no_legend = FALSE,
  title_prefix = NULL,
  reduction = NULL,
  file_path = NULL,
  file_name = NULL,
  file_type = NULL,
  single_pdf = FALSE,
  output_width = NULL,
  output_height = NULL,
  dpi = 600,
  raster = NULL,
  ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check meta data
  meta_data_column <- Meta_Present(object = seurat_object, meta_col_names = meta_data_column, omit_warn = FALSE, print_msg = FALSE)[[1]]

  # stop if none found
  if (length(x = meta_data_column) == 0) {
    cli_abort(message = "None of values provided to {.code meta_data_column} were found in the meta.data slot.")
  }

  # Check that meta data is factor or character
  accepted_meta_types <- c("factor", "character", "logical")

  if (!class(x = seurat_object@meta.data[[meta_data_column]]) %in% accepted_meta_types) {
    cli_abort(message = c("The {.code meta_data_column}: {.val {meta_data_column}} is of class: {.field {class(x = seurat_object@meta.data[[meta_data_column]])}}",
                          "i" = "Only meta data variables of classes: {.field factor, character, or logical} can be used with {.code Meta_Highlight_Plot()}."))
  }

  # Change active ident for plotting
  Idents(object = seurat_object) <- meta_data_column

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path)) {
    if (!dir.exists(paths = file_path)) {
      cli_abort(message = "Provided {.code file_path}: {symbol$dquote_left}{.field {file_path}}{symbol$dquote_right} does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name)) {
    cli_abort(message = "No file name provided.  Please provide a file name using {.code file_name}.")
  }

  # Set file type for single pdf option
  if (isTRUE(x = single_pdf) && is.null(x = file_type)) {
    file_type <- ".pdf"
  }
  if (isTRUE(x = single_pdf) && !is.null(x = file_type) && str_detect(file_type, ".pdf") == FALSE) {
    cli_inform(message = "WARNING: non-PDF {.code file_type} specified but {.code single_pdf = TRUE} selected.  Changing file_type to {.val .pdf} for output.")
    file_type <- ".pdf"
  }

  # Check file_type parameter
  file_type_options <- c(".pdf", ".png", ".tiff", ".jpeg", ".svg")
  if (is.null(x = file_type)) {
    cli_abort(message = c("{.code file_type} not specified.",
                          "*" = "Must specify output file type format from the following:",
                          "i" = "{.field {glue_collapse_scCustom(input_string = file_type_options, and = TRUE)}}"))
  }
  if (!file_type %in% file_type_options) {
    cli_abort(message = "{.code file_type} must be one of the following: {.field {glue_collapse_scCustom(input_string = file_type_options, and = TRUE)}}")
  }

  # Extract default reduction
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)

  # Add raster check for scCustomize
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)

  # Relevel idents for plotting to sorted order
  if (isTRUE(x = single_pdf) && is.null(x = new_meta_order) && meta_data_sort) {
    Idents(object = seurat_object) <- fct_relevel(Idents(object = seurat_object), sort)
  }

  # Relevel idents to custom order
  if (isTRUE(x = single_pdf) && !is.null(x = new_meta_order)) {
    if (length(x = new_meta_order) != length(x = levels(x = seurat_object@active.ident))) {
      cli_abort(message = c("The length of 'new_meta_order' ({.field {length(x = new_meta_order)}}) does not equal the number of levels in {.code meta_data_column}: {.val {meta_data_column}} ({.field {length(x = levels(x = seurat_object@active.ident))}})"))
    }
    Idents(object = seurat_object) <- factor(Idents(object = seurat_object), levels = new_meta_order)
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
      cli_abort(message = "Number of colors provided to {.code highlight_color} ({.field {length(highlight_color)}}) is not equal to the number of clusters (({.field {num_idents}}).")
    } else {
      highlight_color <- highlight_color
    }
  }

  # Create plot titles if needed.
  if (!is.null(x = title_prefix) && isFALSE(x = no_legend)) {
    cli_warn(message = "{.code title_prefix} was omitted as {.code no_legend = FALSE}.")
  }

  if (is.null(x = title_prefix) && isTRUE(x = no_legend)) {
    plot_title <- lapply(1:num_idents, function(z) {
      paste0(meta_data_column, ": ", list_idents[z])
    })
  } else {
    plot_title <- lapply(1:num_idents, function(z) {
      paste0(title_prefix, ": ", list_idents[z])
    })
  }

  if (!is.null(x = title_prefix) && length(x = title_prefix) != 1 && isTRUE(x = no_legend)) {
    cli_abort(message = "{.field `title_prefix`} must be vector of length 1.")
  }

  # Single PDF option
  if (isTRUE(x = single_pdf)) {
    cli_inform(message = "{.field Generating plots}")
    pboptions(char = "=")
    all_plots <- pblapply(1:num_idents, function(x) {
      if (isTRUE(x = no_legend)) {
        suppressMessages(Meta_Highlight_Plot(seurat_object = seurat_object,
                                             meta_data_column = meta_data_column,
                                             meta_data_highlight = list_idents[x],
                                             highlight_color = highlight_color[x],
                                             background_color = background_color,
                                             pt.size = pt.size,
                                             reduction = reduction,
                                             raster = raster,
                                             ...) +
                           NoLegend() +
                           ggtitle(plot_title[x]) +
                           CenterTitle())
      } else {
        suppressMessages(Meta_Highlight_Plot(seurat_object = seurat_object,
                                             meta_data_column = meta_data_column,
                                             meta_data_highlight = list_idents[x],
                                             highlight_color = highlight_color[x],
                                             background_color = background_color,
                                             pt.size = pt.size,
                                             reduction = reduction,
                                             raster = raster,
                                             ...))
      }

    })
    cli_inform(message = "{.field Saving plots to file}")
    pdf(paste(file_path, file_name, file_type, sep=""), width = output_width, height = output_height)
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
      cli_inform(message = "{.field Generating plots and saving plots to file}")
      pb <- txtProgressBar(min = 0, max = num_idents, style = 3, file = stderr())
      for (i in 1:num_idents) {
        if (isTRUE(x = no_legend)) {
          suppressMessages(Meta_Highlight_Plot(seurat_object = seurat_object,
                                               meta_data_column = meta_data_column,
                                               meta_data_highlight = list_idents[i],
                                               highlight_color = highlight_color[i],
                                               background_color = background_color,
                                               pt.size = pt.size,
                                               reduction = reduction,
                                               raster = raster,
                                               ...) +
                             NoLegend() +
                             ggtitle(plot_title[i]) +
                             CenterTitle())
        } else {
          suppressMessages(Meta_Highlight_Plot(seurat_object = seurat_object,
                                               meta_data_column = meta_data_column,
                                               meta_data_highlight = list_idents[i],
                                               highlight_color = highlight_color[i],
                                               background_color = background_color,
                                               pt.size = pt.size,
                                               reduction = reduction,
                                               raster = raster,
                                               ...))
        }

        suppressMessages(ggsave(filename = paste(file_path, list_idents_save[i], "_", file_name, file_type, sep=""), dpi = dpi))
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
    # Code for PDF Version
    if (str_detect(file_type, ".pdf") == TRUE) {
      cli_inform(message = "{.field Generating plots and saving plots to file}")
      pb <- txtProgressBar(min = 0, max = num_idents, style = 3, file = stderr())
      for (i in 1:num_idents) {
        if (isTRUE(x = no_legend)) {
          suppressMessages(Meta_Highlight_Plot(seurat_object = seurat_object,
                                               meta_data_column = meta_data_column,
                                               meta_data_highlight = list_idents[i],
                                               highlight_color = highlight_color[i],
                                               background_color = background_color,
                                               pt.size = pt.size,
                                               reduction = reduction,
                                               raster = raster,
                                               ...) +
                             NoLegend() +
                             ggtitle(plot_title[i]) +
                             CenterTitle())
        } else {
          suppressMessages(Meta_Highlight_Plot(seurat_object = seurat_object,
                                               meta_data_column = meta_data_column,
                                               meta_data_highlight = list_idents[i],
                                               highlight_color = highlight_color[i],
                                               background_color = background_color,
                                               pt.size = pt.size,
                                               reduction = reduction,
                                               raster = raster,
                                               ...))
        }

        suppressMessages(ggsave(filename = paste(file_path, list_idents_save[[i]], "_", file_name, file_type, sep=""), useDingbats = FALSE))
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
  }
}


#' Iterative Plotting of Gene Lists using Custom FeaturePlots
#'
#' Create and Save plots for Gene list with Single Command
#'
#' @param seurat_object Seurat object name.
#' @param features vector of features to plot.  If a named vector is provided then the names for each gene
#'  will be incorporated into plot title if `single_pdf = TRUE` or into file name if `FALSE`.
#' @param colors_use color scheme to use.
#' @param na_color color for non-expressed cells.
#' @param na_cutoff Value to use as minimum expression cutoff.  To set no cutoff set to `NA`.
#' @param split.by Variable in `@meta.data` to split the plot by.
#' @param order whether to move positive cells to the top (default = TRUE).
#' @param return_plots logical. Whether to return plots to list instead of saving them to file(s).  Default is FALSE.
#' @param file_path directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix and file extension.
#' @param file_type File type to save output as.  Must be one of following: ".pdf", ".png", ".tiff", ".jpeg", or ".svg".
#' @param single_pdf saves all plots to single PDF file (default = FALSE).
#' @param features_per_page numeric, number of features to plot on single page if `single_pdf = TRUE`.  Default is 1.
#' @param num_columns Number of columns in plot layout (only applicable if `single_pdf = TRUE` AND
#' `features_per_page` > 1).
#' @param landscape logical, when plotting multiple features per page in single PDF whether to use landscape or portrait
#' page dimensions (default is TRUE).
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
#' @import cli
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom pbapply pblapply pboptions
#' @importFrom Seurat FeaturePlot
#' @importFrom SeuratObject DefaultDimReduc
#' @importFrom stringr str_detect
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @return Saved plots
#'
#' @export
#'
#' @concept iterative_plotting
#'
#' @examples
#' \dontrun{
#' Iterate_FeaturePlot_scCustom(seurat_object = object, features = DEG_list,
#' colors_use = viridis_plasma_dark_high, na_color = "lightgray", file_path = "plots/",
#' file_name = "tsne", file_type = ".jpg", dpi = 600)
#' }
#'

Iterate_FeaturePlot_scCustom <- function(
  seurat_object,
  features,
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
  output_width = NULL,
  output_height = NULL,
  features_per_page = 1,
  num_columns = NULL,
  landscape = TRUE,
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
  raster <- raster %||% (length(x = Cells(x = seurat_object)) > 2e5)

  # Return plot check
  if (isTRUE(x = return_plots)) {
    if (!is.null(x = file_type) | !is.null(x = file_path) | !is.null(x = file_name) | isTRUE(x = single_pdf)) {
      cli_abort(message = c("Cannot return plots to list and save plots to file with single function call.",
                            "i" = "If {.field saving plots} please set {.code return_plots = FALSE}.",
                            "i" = "If {.field returning plots} please leave {.code file_type}, {.code file_path}, {.code file_name} and {.code single_pdf} at their default settings."))
    }
  }

  # Check num_columns validity
  if (!is.null(x = num_columns) && !isTRUE(x = single_pdf)) {
    cli_warn(message = c("{.code num_columns} is only valid when {.code single_pdf = TRUE}",
                         "i" = "Setting {.num_columns = NULL}"))
    num_columns <- NULL
  }

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path) && isFALSE(x = return_plots)) {
    if (!dir.exists(paths = file_path)) {
      cli_abort(message = "Provided {.code file_path}: {symbol$dquote_left}{.field {file_path}}{symbol$dquote_right} does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name) && isFALSE(return_plots)) {
    cli_abort(message = "No file name provided.  Please provide a file name using {.code file_name}.")
  }

  # multi-plot checks
  if (isFALSE(x = single_pdf) && features_per_page != 1) {
    cli_warn(message = "{.code features_per_page} only applicable when {.code single_pdf = TRUE}.")
  }

  if (isFALSE(x = is.numeric(x = features_per_page))) {
    cli_abort(message = "{.code features_per_page} must be numeric value.")
  }

  if (isTRUE(x = is.numeric(x = features_per_page)) && isFALSE(x = check_whole_num(x = features_per_page))) {
    cli_abort(message = "{.code features_per_page} must be whole numeric value.")
  }

  # Extract default reduction
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)

  # Set file type for single pdf option
  if (isTRUE(x = single_pdf) && is.null(x = file_type)) {
    file_type <- ".pdf"
  }
  if (isTRUE(x = single_pdf) && !is.null(x = file_type) && str_detect(file_type, ".pdf") == FALSE) {
    cli_inform(message = "WARNING: non-PDF {.code file_type} specified but {.code single_pdf = TRUE} selected.  Changing file_type to {.val .pdf} for output.")
    file_type <- ".pdf"
  }

  # Check file_type parameter
  file_type_options <- c(".pdf", ".png", ".tiff", ".jpeg", ".svg")
  if (is.null(x = file_type) && isFALSE(x = return_plots)) {
    cli_abort(message = c("{.code file_type} not specified.",
                          "*" = "Must specify output file type format from the following:",
                          "i" = "{.field {glue_collapse_scCustom(input_string = file_type_options, and = TRUE)}}"))
  }

  if (!file_type %in% file_type_options && isFALSE(x = return_plots)) {
    cli_abort(message = "{.code file_type} must be one of the following: {.field {glue_collapse_scCustom(input_string = file_type_options, and = TRUE)}}")
  }

  # Check whether features are present in object (dependent on whether vector is named)
  if (is.null(x = names(x = features))) {
    all_found_features <- Feature_PreCheck(object = seurat_object, features = features)
  } else {
    all_found_features <- features
  }

  if (any(features %in% colnames(seurat_object@meta.data)) && any(features %in% rownames(seurat_object))) {
    cli_warn(message = c("Some of the {.code features} provided are from both assay features and meta.data",
                         "*" = "This could cause problems in plot output due to differences in {.field na_cutoff} parameter.",
                         "i" = "Suggest splitting {.code features} and running {.field Iterate_FeaturePlot_scCustom} once for each feature list."))
  }

  # Modify Cluster Labels names if needed for saving plots
  if (!is.null(x = names(x = all_found_features)) && isFALSE(x = single_pdf)) {
    names_vec_mod <- gsub(pattern = "/", replacement = "-", x = names(x = all_found_features))
    names(x = all_found_features) <- names_vec_mod
  }

  # Return plots instead of saving them
  if (isTRUE(x = return_plots)) {
    cli_inform(message = "{.field Generating plots}")
    pboptions(char = "=")
    all_plots <- pblapply(all_found_features,function(gene) {FeaturePlot_scCustom(seurat_object = seurat_object, features = gene, colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, split.by = split.by, order = order, pt.size = pt.size, reduction = reduction, raster = raster, alpha_exp = alpha_exp, alpha_na_exp = alpha_na_exp, ...)})
    return(all_plots)
  }

  # Single PDF option
  if (isTRUE(x = single_pdf)) {
    # plot if one fearture per page
    if (features_per_page == 1) {
      cli_inform(message = "{.field Generating plots}")
      pboptions(char = "=")
      all_plots <- pblapply(all_found_features,function(gene) {FeaturePlot_scCustom(seurat_object = seurat_object, features = gene, colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, split.by = split.by, order = order, pt.size = pt.size, reduction = reduction, raster = raster, alpha_exp = alpha_exp, alpha_na_exp = alpha_na_exp,...)})
      cli_inform(message = "{.field Saving plots to file}")
      # save plots with cluster annotation
      if (!is.null(x = names(x = all_found_features)) && is.null(x = split.by)) {
        pdf(paste(file_path, file_name, file_type, sep=""), width = output_width, height = output_height)
        pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
        for (i in 1:length(all_plots)) {
          print(all_plots[[i]] + ggtitle((paste0(all_found_features[i], "_", names(x = all_found_features)[i]))))
          setTxtProgressBar(pb = pb, value = i)
        }
        close(con = pb)
        dev.off()
      } else {
        # Save plots without cluster annotation
        pdf(paste(file_path, file_name, file_type, sep=""), width = output_width, height = output_height)
        pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
        for (i in 1:length(all_plots)) {
          print(all_plots[[i]])
          setTxtProgressBar(pb = pb, value = i)
        }
        close(con = pb)
        dev.off()
      }
    } else {
      # for plotting multiple features per page

      # split features by
      features_split <- Split_Vector(x = all_found_features, chunk_size = features_per_page, verbose = FALSE)

      cli_inform(message = "{.field Generating plots}")
      pboptions(char = "=")
      all_plots <- pblapply(features_split, function(z) {FeaturePlot_scCustom(seurat_object = seurat_object, features = z, colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, split.by = split.by, order = order, pt.size = pt.size, reduction = reduction, raster = raster, alpha_exp = alpha_exp, alpha_na_exp = alpha_na_exp, num_columns = num_columns, ...)})



      cli_inform(message = "{.field Saving plots to file}")
      if (isTRUE(x = landscape)) {
        # save plots with cluster annotation
        if (!is.null(x = names(x = all_found_features)) && is.null(x = split.by)) {
          pdf(paste(file_path, file_name, file_type, sep=""), width = 22, height = 17)
          pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())

          list_names <- lapply(1:length(x = features_split), function(k){
            feat_name <- features_split[[k]]
            clu_name <- names(x = features_split[[k]])
            new_names <- paste0(feat_name, "_", clu_name)
          })

          all_plots <- lapply(1:length(x = all_plots), function (j){
            plot_split <- all_plots[[j]]
            for (k in 1:length(x = list_names[[j]])) {
              plot_split[[k]][["labels"]][["title"]] <- list_names[[j]][k]
            }
            return(plot_split)
          })

          for (i in 1:length(x = all_plots)) {
            print(all_plots[[i]])
            setTxtProgressBar(pb = pb, value = i)
          }
          close(con = pb)
          dev.off()
        } else {
          # Save plots without cluster annotation
          pdf(paste(file_path, file_name, file_type, sep=""), width = 22, height = 17)
          pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
          for (i in 1:length(all_plots)) {
            print(all_plots[[i]])
            setTxtProgressBar(pb = pb, value = i)
          }
          close(con = pb)
          dev.off()
        }
      } else {
        if (!is.null(x = names(x = all_found_features)) && is.null(x = split.by)) {
          pdf(paste(file_path, file_name, file_type, sep=""), width = 17, height = 22)
          pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())

          list_names <- lapply(1:length(x = features_split), function(k){
            feat_name <- features_split[[k]]
            clu_name <- names(x = features_split[[k]])
            new_names <- paste0(feat_name, "_", clu_name)
          })

          all_plots <- lapply(1:length(x = all_plots), function (j){
            plot_split <- all_plots[[j]]
            for (k in 1:length(x = list_names[[j]])) {
              plot_split[[k]][["labels"]][["title"]] <- list_names[[j]][k]
            }
            return(plot_split)
          })

          for (i in 1:length(x = all_plots)) {
            print(all_plots[[i]])
            setTxtProgressBar(pb = pb, value = i)
          }
          close(con = pb)
          dev.off()
        } else {
          # Save plots without cluster annotation
          pdf(paste(file_path, file_name, file_type, sep=""), width = 17, height = 22)
          pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
          for (i in 1:length(all_plots)) {
            print(all_plots[[i]])
            setTxtProgressBar(pb = pb, value = i)
          }
          close(con = pb)
          dev.off()
        }
      }
    }
  }
  else {
    if (str_detect(file_type, ".pdf") == FALSE) {
      cli_inform(message = "{.field Generating plots and saving plots to file}")
      pb <- txtProgressBar(min = 0, max = length(all_found_features), style = 3, file = stderr())
      for (i in 1:length(all_found_features)) {
        FeaturePlot_scCustom(seurat_object = seurat_object, features = all_found_features[i], colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, split.by = split.by, order = order, pt.size = pt.size, reduction = reduction, raster = raster, alpha_exp = alpha_exp, alpha_na_exp = alpha_na_exp, ...)
        if (!is.null(x = names(x = all_found_features))) {
          suppressMessages(ggsave(filename = paste(file_path, all_found_features[i], "_", names(x = all_found_features)[i], "_", file_name, file_type, sep=""), dpi = dpi))
        } else {
          suppressMessages(ggsave(filename = paste(file_path, all_found_features[i], "_", file_name, file_type, sep=""), dpi = dpi))
        }
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
    if (str_detect(file_type, ".pdf") == TRUE) {
      cli_inform(message = "{.field Generating plots and saving plots to file}")
      pb <- txtProgressBar(min = 0, max = length(all_found_features), style = 3, file = stderr())
      for (i in 1:length(all_found_features)) {
        FeaturePlot_scCustom(seurat_object = seurat_object, features = all_found_features[i], colors_use = colors_use, na_color = na_color, na_cutoff = na_cutoff, split.by = split.by, order = order, pt.size = pt.size, reduction = reduction, raster = raster, alpha_exp = alpha_exp, alpha_na_exp = alpha_na_exp, ...)
        if (!is.null(x = names(x = all_found_features))) {
          suppressMessages(ggsave(filename = paste(file_path, all_found_features[i], "_", names(x = all_found_features)[i], "_", file_name, file_type, sep=""), useDingbats = FALSE))
        } else {
          suppressMessages(ggsave(filename = paste(file_path, all_found_features[i], "_", file_name, file_type, sep=""), useDingbats = FALSE))
        }
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
  }

  # One warning rastering
  if (isFALSE(x = raster) && isTRUE(x = single_pdf) && getOption(x = 'scCustomize_warn_raster_iterative', default = TRUE)) {
    cli_inform(message = c("",
                           "NOTE: {.code single_pdf = TRUE} and {.code raster = FALSE},",
                           "Saving large numbers of plots in vector form can result in very large file sizes.",
                           "Suggest setting {.code raster = TRUE} when plotting large numbers of features in",
                           "single output file.",
                           "",
                           "-----This message will be shown once per session.-----"))
    options(scCustomize_warn_raster_iterative = FALSE)
  }
}


#' Iterative Plotting of Gene Lists using VlnPlot_scCustom
#'
#' Create and Save plots for Gene list with Single Command
#'
#' @param seurat_object Seurat object name.
#' @param features vector of features to plot.
#' @param colors_use color palette to use for plotting.  By default if number of levels plotted is less than
#' or equal to 36 it will use "polychrome" and if greater than 36 will use "varibow" with shuffle = TRUE
#' both from `DiscretePalette_scCustomize`.
#' @param pt.size point size for plotting.
#' @param group.by Name of one or more metadata columns to group (color) plot by (for example, orig.ident);
#' default is the current active.ident of the object.
#' @param split.by Feature to split plots by (i.e. "orig.ident").
#' @param file_path directory file path and/or file name prefix.  Defaults to current wd.
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
#' @import cli
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom pbapply pblapply pboptions
#' @importFrom Seurat VlnPlot
#' @importFrom stringr str_detect
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @return Saved plots
#'
#' @export
#'
#' @concept iterative_plotting
#'
#' @examples
#' \dontrun{
#' Iterate_VlnPlot_scCustom(seurat_object = object, features = DEG_list, colors = color_list,
#' file_path = "plots/", file_name = "_vln", file_type = ".jpg", dpi = 600)
#' }
#'

Iterate_VlnPlot_scCustom <- function(
  seurat_object,
  features,
  colors_use = NULL,
  pt.size = NULL,
  group.by = NULL,
  split.by = NULL,
  file_path = NULL,
  file_name = NULL,
  file_type = NULL,
  single_pdf = FALSE,
  output_width = NULL,
  output_height = NULL,
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
      cli_abort(message = "Provided {.code file_path}: {symbol$dquote_left}{.field {file_path}}{symbol$dquote_right} does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name)) {
    cli_abort(message = "No file name provided.  Please provide a file name using {.code file_name}.")
  }

  # Set file type for single pdf option
  if (isTRUE(x = single_pdf) && is.null(x = file_type)) {
    file_type <- ".pdf"
  }
  if (isTRUE(x = single_pdf) && !is.null(x = file_type) && str_detect(file_type, ".pdf") == FALSE) {
    cli_inform(message = "WARNING: non-PDF {.code file_type} specified but {.code single_pdf = TRUE} selected.  Changing file_type to {.val .pdf} for output.")
    file_type <- ".pdf"
  }

  # Check file_type parameter
  file_type_options <- c(".pdf", ".png", ".tiff", ".jpeg", ".svg")
  if (is.null(x = file_type)) {
    cli_abort(message = c("{.code file_type} not specified.",
                          "*" = "Must specify output file type format from the following:",
                          "i" = "{.field {glue_collapse_scCustom(input_string = file_type_options, and = TRUE)}}"))
  }
  if (!file_type %in% file_type_options) {
    cli_abort(message = "{.code file_type} must be one of the following: {.field {glue_collapse_scCustom(input_string = file_type_options, and = TRUE)}}")
  }

  # Check whether features are present in object (dependent on whether vector is named)
  if (is.null(x = names(x = features))) {
    all_found_features <- Feature_PreCheck(object = seurat_object, features = features)
  } else {
    all_found_features <- features
  }

  # Set default color palette based on number of levels being plotted
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }

  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use) && isTRUE(x = ggplot_default_colors)) {
    cli_abort(message = "Cannot provide both custom palette to {.code colors_use} and specify {.code ggplot_default_colors = TRUE}.")
  }
  if (is.null(x = colors_use)) {
    # set default plot colors
    if (is.null(x = colors_use)) {
      colors_use <- scCustomize_Palette(num_groups = group_by_length, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed)
    }
  }

  # Add one time raster warning
  if (isTRUE(x = single_pdf) && pt.size != 0 && isFALSE(raster) && getOption(x = 'scCustomize_warn_vln_raster_iterative', default = TRUE)) {
    cli_inform(message = c("",
                           "NOTE: {.code single_pdf = TRUE} and {.code pt.size} > 0 and {.code raster = FALSE},",
                           "so all points are plotted.",
                           "Saving large numbers of plots in vector form can result in very large",
                           "file sizes. Suggest setting {.code pt.size = 0} or {.code raster = TRUE}\n",
                           "when plotting large numbers of features in single output file.",
                           "",
                           "-----This message will be shown once per session.-----"))
    options(scCustomize_warn_vln_raster_iterative = FALSE)
  }

  # Single PDF option
  if (isTRUE(x = single_pdf)) {
    cli_inform(message = "{.field Generating plots}")
    pboptions(char = "=")
    all_plots <- pblapply(all_found_features,function(gene) {VlnPlot_scCustom(seurat_object = seurat_object, features = gene, colors_use = colors_use, pt.size = pt.size, group.by = group.by, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, split.by = split.by, ...)})
    cli_inform(message = "{.field Saving plots to file}")
    pdf(paste(file_path, file_name, file_type, sep=""), width = output_width, height = output_height)
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
      cli_inform(message = "{.field Generating plots and saving plots to file}")
      pb <- txtProgressBar(min = 0, max = length(x = all_found_features), style = 3, file = stderr())
      for (i in 1:length(x = all_found_features)) {
        VlnPlot_scCustom(seurat_object = seurat_object, features = all_found_features[i], colors_use = colors_use, pt.size = pt.size, group.by = group.by, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, split.by = split.by, ...)
        suppressMessages(ggsave(filename = paste(file_path, all_found_features[i], file_name, file_type, sep=""), dpi = dpi))
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
    }
    if (str_detect(file_type, ".pdf") == TRUE) {
      cli_inform(message = "{.field Generating plots and saving plots to file}")
      pb <- txtProgressBar(min = 0, max = length(x = all_found_features), style = 3, file = stderr())
      for (i in 1:length(x = all_found_features)) {
        VlnPlot_scCustom(seurat_object = seurat_object, features = all_found_features[i], colors_use = colors_use, pt.size = pt.size, group.by = group.by, raster = raster, ggplot_default_colors = ggplot_default_colors, color_seed = color_seed, split.by = split.by, ...)
        suppressMessages(ggsave(filename = paste(file_path, all_found_features[i], file_name, file_type, sep=""), useDingbats = FALSE))
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
#' @param pt.size Adjust point size for plotting.
#' @param file_path directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix and file extension.
#' @param file_type File type to save output as.  Must be one of following: ".pdf", ".png", ".tiff", ".jpeg", or ".svg".
#' @param single_pdf saves all plots to single PDF file (default = FALSE).  `file_type`` must be .pdf.
#' @param dpi dpi for image saving.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default)
#' @param joint NULL.  This function only supports `joint = FALSE`.  Leave as NULL to generate plots.  To iterate joint plots see function: `Iterate_Plot_Density_Joint`.
#' @param combine Create a single plot? If FALSE, a list with ggplot objects is returned.
#' @param ... Extra parameters passed to \code{\link[Nebulosa]{plot_density}}.
#'
#' @import cli
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom pbapply pblapply pboptions
#' @importFrom rlang is_installed
#' @importFrom SeuratObject DefaultDimReduc
#' @importFrom stringr str_detect
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @return Saved plots
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
  output_width = NULL,
  output_height = NULL,
  dpi = 600,
  reduction = NULL,
  combine = TRUE,
  joint = FALSE,
  ...
) {
  # Check Nebulosa installed
  Nebulosa_check <- is_installed(pkg = "Nebulosa")
  if (isFALSE(x = Nebulosa_check)) {
    cli_abort(message = c(
      "Please install the {.val Nebulosa} package to use {.code Iterate_Plot_Density_Custom}",
      "i" = "This can be accomplished with the following commands: ",
      "----------------------------------------",
      "{.field `install.packages({symbol$dquote_left}BiocManager{symbol$dquote_right})`}",
      "{.field `BiocManager::install({symbol$dquote_left}Nebulosa{symbol$dquote_right})`}",
      "----------------------------------------"
    ))
  }

  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # joint check
  if (!is.null(x = joint)) {
    cli_abort(message = c("{.code Iterate_Plot_Density_Custom} only supports {.code joint = FALSE}.",
                          "i" = "Leave as NULL to generate plots.",
                          "i" = "To iterate joint plots see function: {.code Iterate_Plot_Density_Joint}."))
  }

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path)) {
    if (!dir.exists(paths = file_path)) {
      cli_abort(message = "Provided {.code file_path}: {symbol$dquote_left}{.field {file_path}}{symbol$dquote_right} does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name)) {
    cli_abort(message = "No file name provided.  Please provide a file name using {.code file_name}.")
  }

  # Extract default reduction
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)

  # Set file type for single pdf option
  if (isTRUE(x = single_pdf) && is.null(x = file_type)) {
    file_type <- ".pdf"
  }
  if (isTRUE(x = single_pdf) && !is.null(x = file_type) && str_detect(file_type, ".pdf") == FALSE) {
    cli_inform(message = "WARNING: non-PDF {.code file_type} specified but {.code single_pdf = TRUE} selected.  Changing file_type to {.val .pdf} for output.")
    file_type <- ".pdf"
  }

  # Check file_type parameter
  file_type_options <- c(".pdf", ".png", ".tiff", ".jpeg", ".svg")
  if (is.null(x = file_type)) {
    cli_abort(message = c("{.code file_type} not specified.",
                          "*" = "Must specify output file type format from the following:",
                          "i" = "{.field {glue_collapse_scCustom(input_string = file_type_options, and = TRUE)}}"))
  }
  if (!file_type %in% file_type_options) {
    cli_abort(message = "{.code file_type} must be one of the following: {.field {glue_collapse_scCustom(input_string = file_type_options, and = TRUE)}}")
  }

  # Check whether features are present in object
  gene_list <- Feature_Present(data = seurat_object, features = gene_list, print_msg = FALSE, case_check = TRUE)[[1]]

  # check palettes
  if (!is.null(x = custom_palette) && viridis_palette != "magma") {
    cli_abort(message = c("Non-default values provided to both {.code viridis_palette} & {.code custom_palette}.",
                          "i" = "Please chose one non-default value."))
  }

  # Modify Cluster Labels names if needed for saving plots
  if (!is.null(x = names(x = gene_list)) && isFALSE(x = single_pdf)) {
    names_vec_mod <- gsub(pattern = "/", replacement = "-", x = names(x = gene_list))
    names(x = gene_list) <- names_vec_mod
  }

  # Single PDF option
  if (isTRUE(x = single_pdf)) {
    cli_inform(message = "{.field Generating plots}")
    pboptions(char = "=")
    all_plots <- pblapply(gene_list,function(gene) {
      Plot_Density_Custom(seurat_object = seurat_object, features = gene, joint = FALSE, viridis_palette = viridis_palette, custom_palette = custom_palette, pt.size = pt.size, reduction = reduction, ...)})
    cli_inform(message = "{.field Saving plots to file}")
    # save plots with cluster annotation
    if (!is.null(x = names(x = gene_list))) {
      pdf(paste(file_path, file_name, file_type, sep=""), width = output_width, height = output_height)
      pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
      for (i in 1:length(all_plots)) {
        print(all_plots[[i]] + ggtitle((paste0(gene_list[i], "_", names(x = gene_list)[i]))))
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
      dev.off()
    } else {
      # Save plots without cluster annotation
      pdf(paste(file_path, file_name, file_type, sep=""), width = output_width, height = output_height)
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
      cli_inform(message = "{.field Generating plots and saving plots to file}")
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
      cli_inform(message = "{.field Generating plots and saving plots to file}")
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
#' @param pt.size Adjust point size for plotting.
#' @param file_path directory file path and/or file name prefix.  Defaults to current wd.
#' @param file_name name suffix and file extension.
#' @param file_type File type to save output as.  Must be one of following: ".pdf", ".png", ".tiff", ".jpeg", or ".svg".
#' @param single_pdf saves all plots to single PDF file (default = FALSE).  `file_type`` must be .pdf.
#' @param dpi dpi for image saving.
#' @param reduction Dimensionality Reduction to use (if NULL then defaults to Object default)
#' @param joint NULL.  This function only supports `joint = FALSE`.  Leave as NULL to generate plots.  To iterate joint plots see function: `Iterate_Plot_Density_Joint`.
#' @param combine Create a single plot? If FALSE, a list with ggplot objects is returned.
#' @param ... Extra parameters passed to \code{\link[Nebulosa]{plot_density}}.
#'
#' @import cli
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom pbapply pblapply pboptions
#' @importFrom purrr discard keep
#' @importFrom rlang is_installed
#' @importFrom SeuratObject DefaultDimReduc
#' @importFrom stringr str_detect
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @return Saved plots
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
  output_width = NULL,
  output_height = NULL,
  dpi = 600,
  reduction = NULL,
  combine = TRUE,
  joint = NULL,
  ...
) {
  # Check Nebulosa installed
  Nebulosa_check <- is_installed("Nebulosa")
  if (isFALSE(Nebulosa_check)) {
    cli_abort(message = c(
      "Please install the {.val Nebulosa} package to use {.code Iterate_Plot_Density_Joint}",
      "i" = "This can be accomplished with the following commands: ",
      "----------------------------------------",
      "{.field `install.packages({symbol$dquote_left}BiocManager{symbol$dquote_right})`}",
      "{.field `BiocManager::install({symbol$dquote_left}Nebulosa{symbol$dquote_right})`}",
      "----------------------------------------"
    ))
  }

  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check gene list is in list form
  if (!inherits(x = gene_list, what = "list")) {
    cli_abort(message = c("For {.code Iterate_Plot_Density_Joint} the {.code gene_list} must be of class(): {.val list}",
                          "i" = "Please reformat from current class(): {.val {class(x = gene_list)}}"
    ))
  }

  # joint check
  if (!is.null(x = joint)) {
    cli_abort(message = c("{.code Iterate_Plot_Density_Joint} only supports {.code joint = FALSE}.",
                          "i" = "Leave as NULL to generate plots.",
                          "i" = "To iterate joint plots see function: {.code Iterate_Plot_Density_Custom}."))
  }

  # Set file_path before path check if current dir specified as opposed to leaving set to NULL
  if (!is.null(x = file_path) && file_path == "") {
    file_path <- NULL
  }

  # Check file path is valid
  if (!is.null(x = file_path)) {
    if (!dir.exists(paths = file_path)) {
      cli_abort(message = "Provided {.code file_path}: {symbol$dquote_left}{.field {file_path}}{symbol$dquote_right} does not exist.")
    }
  }

  # Check if file name provided
  if (is.null(x = file_name)) {
    cli_abort(message = "No file name provided.  Please provide a file name using {.code file_name}.")
  }

  # Extract default reduction
  reduction <- reduction %||% DefaultDimReduc(object = seurat_object)

  # Set file type for single pdf option
  if (isTRUE(x = single_pdf) && is.null(x = file_type)) {
    file_type <- ".pdf"
  }
  if (isTRUE(x = single_pdf) && !is.null(x = file_type) && str_detect(file_type, ".pdf") == FALSE) {
    cli_inform(message = "WARNING: non-PDF {.code file_type} specified but {.code single_pdf = TRUE} selected.  Changing file_type to {.val .pdf} for output.")
    file_type <- ".pdf"
  }

  # Check file_type parameter
  file_type_options <- c(".pdf", ".png", ".tiff", ".jpeg", ".svg")
  if (is.null(x = file_type)) {
    cli_abort(message = c("{.code file_type} not specified.",
                          "*" = "Must specify output file type format from the following:",
                          "i" = "{.field {glue_collapse_scCustom(input_string = file_type_options, and = TRUE)}}"))
  }
  if (!file_type %in% file_type_options) {
    cli_abort(message = "{.code file_type} must be one of the following: {.field {glue_collapse_scCustom(input_string = file_type_options, and = TRUE)}}")
  }

  # Check whether features are present in object
  checked_gene_list <- lapply(1:length(gene_list), function(x){
    genes <- Feature_Present(data = seurat_object, features = gene_list[[x]], print_msg = FALSE, case_check = TRUE, return_none = TRUE)[[1]]
  })

  if (!is.null(x = names(x = gene_list))) {
    names(x = checked_gene_list) <- names(x = gene_list)
  }

  # remove any empty entries in list
  checked_gene_list <- discard(checked_gene_list,  ~length(.x) == 0)

  # Check for lists less than 2
  bad_gene_lists <- keep(checked_gene_list,  ~length(.x) < 2)
  if (length(x = bad_gene_lists) > 0) {
    cli_warn(message = "A total of {.field {length(x = bad_gene_lists)}} list entries from {.code gene_list} contain less than two features and were excluded from plotting.")
  }

  # Create final good gene list
  final_gene_list <- keep(checked_gene_list,  ~length(.x) > 1)

  # check palettes
  if (!is.null(x = custom_palette) && viridis_palette != "magma") {
    cli_abort(message = c("Non-default values provided to both {.code viridis_palette} & {.code custom_palette}.",
                          "i" = "Please chose one non-default value."))
  }

  # Modify Cluster Labels names if needed for saving plots
  if (!is.null(x = names(x = gene_list)) && isFALSE(x = single_pdf)) {
    names_vec_mod <- gsub(pattern = "/", replacement = "-", x = names(x = gene_list))
    names(x = gene_list) <- names_vec_mod
  }

  # Single PDF option
  if (isTRUE(x = single_pdf)) {
    cli_inform(message = "{.field Generating plots}")
    pboptions(char = "=")
    all_plots <- pblapply(1:length(final_gene_list),function(i) {
      plot <- Plot_Density_Joint_Only(seurat_object = seurat_object, features = final_gene_list[[i]], viridis_palette = viridis_palette, custom_palette = custom_palette, pt.size = pt.size, reduction = reduction, ...)})
    cli_inform(message = "{.field Saving plots to file}")
    # save plots with cluster annotation
    if (!is.null(x = names(x = final_gene_list))) {
      pdf(paste(file_path, file_name, file_type, sep=""), width = output_width, height = output_height)
      pb <- txtProgressBar(min = 0, max = length(all_plots), style = 3, file = stderr())
      for (i in 1:length(all_plots)) {
        print(all_plots[[i]] + ggtitle((paste0(paste(final_gene_list[[i]], collapse = "_"), "_", names(x = final_gene_list)[i]))))
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
      dev.off()
    } else {
      # Save plots without cluster annotation
      pdf(paste(file_path, file_name, file_type, sep=""), width = output_width, height = output_height)
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
      cli_inform(message = "{.field Generating plots and saving plots to file}")
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
      cli_inform(message = "{.field Generating plots and saving plots to file}")
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
