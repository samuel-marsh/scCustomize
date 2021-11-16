#' Plot Median Genes per Cell per Sample
#'
#' Plot of median genes per cell per sample grouped by desired meta data variable.
#'
#' @param seurat_object Seurat object name.
#' @param sample_col Specify which column in meta.data specifies sample ID (i.e. orig.ident).
#' @param group_by Column in meta.data slot to group results by (i.e. "Treatment").
#' @param colors_use List of colors or color palette to use.  Only applicable if `group_by` is not NULL.
#' @param plot_title Plot title.
#' @param y_axis_label Label for y axis.
#' @param x_axis_label Label for x axis.
#' @param legend_title Label for plot legend.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom dplyr select slice left_join
#' @importFrom magrittr "%>%"
#'
#' @export
#'
#' @concept stats_plotting
#'
#' @examples
#' \dontrun{
#' Plot_Median_Genes(seurat_object = obj, sample_col = "orig.ident",  group_by = "Treatment")
#' }
#'

Plot_Median_Genes <- function(
  seurat_object,
  sample_col = "orig.ident",
  group_by = NULL,
  colors_use = NULL,
  plot_title = "Median Genes/Cell per Sample",
  y_axis_label = "Median Genes",
  x_axis_label = NULL,
  legend_title = NULL,
  x_lab_rotate = TRUE,
  color_seed = 123
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check group by is valid
  group_by <- Meta_Present(seurat_object = seurat_object, meta_col_names = group_by, print_msg = FALSE)[[1]]

  # Check sample_col is valid
  sample_col <- Meta_Present(seurat_object = seurat_object, meta_col_names = sample_col, print_msg = FALSE)[[1]]

  # Calculate medians and merge with meta.data
  medians <- Median_Stats(seurat_object = seurat_object, group_by_var = sample_col, median_var = "nFeature_RNA", default_var = FALSE) %>%
    slice(-n()) %>%
    droplevels()

  meta <- seurat_object@meta.data

  if (!is.null(x = group_by)) {
    meta <- meta %>%
      select(.data[[sample_col]], .data[[group_by]])
  } else {
    meta <- meta %>%
      select(.data[[sample_col]])
  }

  meta[[sample_col]] <- factor(meta[[sample_col]], ordered = FALSE)

  meta <- data.frame(meta[!duplicated(meta[,sample_col]),])

  if (is.null(x = group_by)) {
    colnames(meta) <- sample_col
  }

  merged <- suppressMessages(left_join(medians, meta))

  # Check colors_use
  if (!is.null(x = group_by)) {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group_by]]))
    if (is.null(x = colors_use)) {
      if (group_by_length <= 8) {
        colors_use <- Dark2_Pal()
      } else if (group_by_length < 36) {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "polychrome")
      } else {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "varibow", seed = color_seed)
      }
    }
  }

  # Generate base plot
  if (is.null(x = group_by)) {
    merged$samples_plotting <- "Samples"

    plot <- ggplot(merged, aes(x = samples_plotting, y = Median_nFeature_RNA)) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("") +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(data = merged, mapping = aes(x = .data[[group_by]], y = Median_nFeature_RNA, fill = .data[[group_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      theme_ggprism_mod() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("")
  }

  # Modify base plot
  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (!is.null(x = x_axis_label)) {
    plot <- plot + xlab(x_axis_label)
  }

  if (!is.null(x = legend_title)) {
    plot <- plot + labs(fill = legend_title)
  }

  # Return plot
  return(plot)
}


#' Plot Median UMIs per Cell per Sample
#'
#' Plot of median UMIs per cell per sample grouped by desired meta data variable.
#'
#' @param seurat_object Seurat object name.
#' @param sample_col Specify which column in meta.data specifies sample ID (i.e. orig.ident).
#' @param group_by Column in meta.data slot to group results by (i.e. "Treatment").
#' @param colors_use List of colors or color palette to use.  Only applicable if `group_by` is not NULL.
#' @param plot_title Plot title.
#' @param y_axis_label Label for y axis.
#' @param x_axis_label Label for x axis.
#' @param legend_title Label for plot legend.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom dplyr select slice left_join
#' @importFrom magrittr "%>%"
#'
#' @export
#'
#' @concept stats_plotting
#'
#' @examples
#' \dontrun{
#' Plot_Median_UMIs(seurat_object = obj, sample_col = "orig.ident",  group_by = "Treatment")
#' }
#'

Plot_Median_UMIs <- function(
  seurat_object,
  sample_col = "orig.ident",
  group_by = NULL,
  colors_use = NULL,
  plot_title = "Median UMIs/Cell per Sample",
  y_axis_label = "Median UMIs",
  x_axis_label = NULL,
  legend_title = NULL,
  x_lab_rotate = TRUE,
  color_seed = 123
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check group by is valid
  group_by <- Meta_Present(seurat_object = seurat_object, meta_col_names = group_by, print_msg = FALSE)[[1]]

  # Check sample_col is valid
  sample_col <- Meta_Present(seurat_object = seurat_object, meta_col_names = sample_col, print_msg = FALSE)[[1]]

  # Calculate medians and merge with meta.data
  medians <- Median_Stats(seurat_object = seurat_object, group_by_var = sample_col, median_var = "nCount_RNA", default_var = FALSE) %>%
    slice(-n()) %>%
    droplevels()

  meta <- seurat_object@meta.data

  if (!is.null(x = group_by)) {
    meta <- meta %>%
      select(.data[[sample_col]], .data[[group_by]])
  } else {
    meta <- meta %>%
      select(.data[[sample_col]])
  }

  meta[[sample_col]] <- factor(meta[[sample_col]], ordered = FALSE)

  meta <- data.frame(meta[!duplicated(meta[,sample_col]),])

  if (is.null(x = group_by)) {
    colnames(meta) <- sample_col
  }

  merged <- suppressMessages(left_join(medians, meta))

  # Check colors_use
  if (!is.null(x = group_by)) {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group_by]]))
    if (is.null(x = colors_use)) {
      if (group_by_length <= 8) {
        colors_use <- Dark2_Pal()
      } else if (group_by_length < 36) {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "polychrome")
      } else {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "varibow", seed = color_seed)
      }
    }
  }

  # Generate base plot
  if (is.null(x = group_by)) {
    merged$samples_plotting <- "Samples"

    plot <- ggplot(merged, aes(x = samples_plotting, y = Median_nCount_RNA)) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("") +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(data = merged, mapping = aes(x = .data[[group_by]], y = Median_nCount_RNA, fill = .data[[group_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      theme_ggprism_mod() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("")
  }

  # Modify base plot
  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (!is.null(x = x_axis_label)) {
    plot <- plot + xlab(x_axis_label)
  }

  if (!is.null(x = legend_title)) {
    plot <- plot + labs(fill = legend_title)
  }

  # Return plot
  return(plot)
}


#' Plot Median Percent Mito per Cell per Sample
#'
#' Plot of median percent mito per cell per sample grouped by desired meta data variable.
#'
#' @param seurat_object Seurat object name.
#' @param sample_col Specify which column in meta.data specifies sample ID (i.e. orig.ident).
#' @param group_by Column in meta.data slot to group results by (i.e. "Treatment").
#' @param colors_use List of colors or color palette to use.  Only applicable if `group_by` is not NULL.
#' @param plot_title Plot title.
#' @param y_axis_label Label for y axis.
#' @param x_axis_label Label for x axis.
#' @param legend_title Label for plot legend.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom dplyr select slice left_join
#' @importFrom magrittr "%>%"
#'
#' @export
#'
#' @concept stats_plotting
#'
#' @examples
#' \dontrun{
#' Plot_Median_Mito(seurat_object = obj, sample_col = "orig.ident",  group_by = "Treatment")
#' }
#'

Plot_Median_Mito <- function(
  seurat_object,
  sample_col = "orig.ident",
  group_by = NULL,
  colors_use = NULL,
  plot_title = "Median % Mito per Sample",
  y_axis_label = "Percent Mitochondrial Reads",
  x_axis_label = NULL,
  legend_title = NULL,
  x_lab_rotate = TRUE,
  color_seed = 123
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check group by is valid
  group_by <- Meta_Present(seurat_object = seurat_object, meta_col_names = group_by, print_msg = FALSE)[[1]]

  # Check sample_col is valid
  sample_col <- Meta_Present(seurat_object = seurat_object, meta_col_names = sample_col, print_msg = FALSE)[[1]]

  # Calculate medians and merge with meta.data
  medians <- Median_Stats(seurat_object = seurat_object, group_by_var = sample_col, median_var = "percent_mito", default_var = FALSE) %>%
    slice(-n()) %>%
    droplevels()

  meta <- seurat_object@meta.data

  if (!is.null(x = group_by)) {
    meta <- meta %>%
      select(.data[[sample_col]], .data[[group_by]])
  } else {
    meta <- meta %>%
      select(.data[[sample_col]])
  }

  meta[[sample_col]] <- factor(meta[[sample_col]], ordered = FALSE)

  meta <- data.frame(meta[!duplicated(meta[,sample_col]),])

  if (is.null(x = group_by)) {
    colnames(meta) <- sample_col
  }

  merged <- suppressMessages(left_join(medians, meta))

  # Check colors_use
  if (!is.null(x = group_by)) {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group_by]]))
    if (is.null(x = colors_use)) {
      if (group_by_length <= 8) {
        colors_use <- Dark2_Pal()
      } else if (group_by_length < 36) {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "polychrome")
      } else {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "varibow", seed = color_seed)
      }
    }
  }

  # Generate base plot
  if (is.null(x = group_by)) {
    merged$samples_plotting <- "Samples"

    plot <- ggplot(merged, aes(x = samples_plotting, y = Median_percent_mito)) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("") +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(data = merged, mapping = aes(x = .data[[group_by]], y = Median_percent_mito, fill = .data[[group_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      theme_ggprism_mod() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("")
  }

  # Modify base plot
  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (!is.null(x = x_axis_label)) {
    plot <- plot + xlab(x_axis_label)
  }

  if (!is.null(x = legend_title)) {
    plot <- plot + labs(fill = legend_title)
  }

  # Return plot
  return(plot)
}


#' Plot Median other variable per Cell per Sample
#'
#' Plot of median other variable per cell per sample grouped by desired meta data variable.
#'
#' @param seurat_object Seurat object name.
#' @param median_var Variable in meta.data slot to calculate and plot median values for.
#' @param sample_col Specify which column in meta.data specifies sample ID (i.e. orig.ident).
#' @param group_by Column in meta.data slot to group results by (i.e. "Treatment").
#' @param colors_use List of colors or color palette to use.  Only applicable if `group_by` is not NULL.
#' @param plot_title Plot title.
#' @param y_axis_label Label for y axis.
#' @param x_axis_label Label for x axis.
#' @param legend_title Label for plot legend.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom dplyr select slice left_join
#' @importFrom magrittr "%>%"
#'
#' @export
#'
#' @concept stats_plotting
#'
#' @examples
#' \dontrun{
#' Plot_Median_Other(seurat_object = obj, median_var = "module_score", sample_col = "orig.ident", group_by = "Treatment")
#' }
#'

Plot_Median_Other <- function(
  seurat_object,
  median_var,
  sample_col = "orig.ident",
  group_by = NULL,
  colors_use = NULL,
  plot_title = NULL,
  y_axis_label = NULL,
  x_axis_label = NULL,
  legend_title = NULL,
  x_lab_rotate = TRUE,
  color_seed = 123
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Set plot and y axis labels if NULL
  if (is.null(x = plot_title)) {
    plot_title <- paste0("Median ", median_var, " per Sample")
  }

  if (is.null(x = y_axis_label)) {
    y_axis_label <- paste0("Median ", median_var)
  }

  # Check group by is valid
  group_by <- Meta_Present(seurat_object = seurat_object, meta_col_names = group_by, print_msg = FALSE)[[1]]

  # Check sample_col is valid
  sample_col <- Meta_Present(seurat_object = seurat_object, meta_col_names = sample_col, print_msg = FALSE)[[1]]

  # Calculate medians and merge with meta.data
  medians <- Median_Stats(seurat_object = seurat_object, group_by_var = sample_col, median_var = median_var, default_var = FALSE) %>%
    slice(-n()) %>%
    droplevels()

  meta <- seurat_object@meta.data

  if (!is.null(x = group_by)) {
    meta <- meta %>%
      select(.data[[sample_col]], .data[[group_by]])
  } else {
    meta <- meta %>%
      select(.data[[sample_col]])
  }

  meta[[sample_col]] <- factor(meta[[sample_col]], ordered = FALSE)

  meta <- data.frame(meta[!duplicated(meta[,sample_col]),])

  if (is.null(x = group_by)) {
    colnames(meta) <- sample_col
  }

  merged <- suppressMessages(left_join(medians, meta))

  # Check colors_use
  if (!is.null(x = group_by)) {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group_by]]))
    if (is.null(x = colors_use)) {
      if (group_by_length <= 8) {
        colors_use <- Dark2_Pal()
      } else if (group_by_length < 36) {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "polychrome")
      } else {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "varibow", seed = color_seed)
      }
    }
  }

  # Generate base plot
  if (is.null(x = group_by)) {
    merged$samples_plotting <- "Samples"

    plot <- ggplot(merged, aes(x = samples_plotting, y = .data[[paste0("Median_", median_var)]])) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("") +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(data = merged, mapping = aes(x = .data[[group_by]], y = .data[[paste0("Median_", median_var)]], fill = .data[[group_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      theme_ggprism_mod() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("")
  }

  # Modify base plot
  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (!is.null(x = x_axis_label)) {
    plot <- plot + xlab(x_axis_label)
  }

  if (!is.null(x = legend_title)) {
    plot <- plot + labs(fill = legend_title)
  }

  # Return plot
  return(plot)
}


#' Plot Number of Cells/Nuclei per Sample
#'
#' Plot of total cell or nuclei number per sample grouped by another meta data variable.
#'
#' @param seurat_object Seurat object name.
#' @param sample_col Specify which column in meta.data specifies sample ID (i.e. orig.ident).
#' @param group_by Column in meta.data slot to group results by (i.e. "Treatment").
#' @param colors_use List of colors or color palette to use.
#' @param plot_title Plot title.
#' @param y_axis_label Label for y axis.
#' @param x_axis_label Label for x axis.
#' @param legend_title Label for plot legend.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom dplyr select slice left_join rename
#' @importFrom magrittr "%>%"
#'
#' @export
#'
#' @concept stats_plotting
#'
#' @examples
#' \dontrun{
#' Plot_Cells_per_Sample(seurat_object = obj, sample_col = "orig.ident", group_by = "Treatment")
#' }
#'

Plot_Cells_per_Sample <- function(
  seurat_object,
  sample_col = "orig.ident",
  group_by = NULL,
  colors_use = NULL,
  plot_title = "Cells/Nuclei per Sample",
  y_axis_label = "Number of Cells",
  x_axis_label = NULL,
  legend_title = NULL,
  x_lab_rotate = TRUE,
  color_seed = 123
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check grouping variable is present
  if (is.null(x = group_by)) {
    stop("Must provided meta data variable to `group_by` in order to plot data.")
  }

  # Check group by is valid
  group_by <- Meta_Present(seurat_object = seurat_object, meta_col_names = group_by, print_msg = FALSE)[[1]]

  # Check sample_col is valid
  sample_col <- Meta_Present(seurat_object = seurat_object, meta_col_names = sample_col, print_msg = FALSE)[[1]]

  # Calculate total cells and merge with meta.data
  total_cells <- table(seurat_object@meta.data[[sample_col]]) %>%
    data.frame() %>%
    rename(!!sample_col := Var1, Number_of_Cells = Freq)

  meta <- seurat_object@meta.data

  meta <- meta %>%
    select(.data[[sample_col]], .data[[group_by]])

  meta[[sample_col]] <- factor(meta[[sample_col]], ordered = FALSE)

  meta <- data.frame(meta[!duplicated(meta[,sample_col]),])

  merged <- suppressMessages(left_join(total_cells, meta))

  # Check colors_use
  if (!is.null(x = group_by)) {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group_by]]))
    if (is.null(x = colors_use)) {
      if (group_by_length <= 8) {
        colors_use <- Dark2_Pal()
      } else if (group_by_length < 36) {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "polychrome")
      } else {
        colors_use <- DiscretePalette_scCustomize(num_colors = group_by_length, palette = "varibow", seed = color_seed)
      }
    }
  }

  # Generate base plot
  plot <- ggplot(data = merged, mapping = aes(x = .data[[group_by]], y = Number_of_Cells, fill = .data[[group_by]])) +
    geom_boxplot(fill = "white") +
    geom_dotplot(binaxis ='y', stackdir = 'center') +
    scale_fill_manual(values = colors_use) +
    theme_ggprism_mod() +
    ggtitle(plot_title) +
    ylab(y_axis_label) +
    xlab("")

  # Modify base plot
  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (!is.null(x = x_axis_label)) {
    plot <- plot + xlab(x_axis_label)
  }

  if (!is.null(x = legend_title)) {
    plot <- plot + labs(fill = legend_title)
  }

  # Return plot
  return(plot)
}
