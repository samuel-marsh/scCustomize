#' Plot Median Genes per Cell per Sample
#'
#' Plot of median genes per cell per sample grouped by desired meta data variable.
#'
#' @param seurat_object Seurat object name.
#' @param sample_col Specify which column in meta.data specifies sample ID (i.e. orig.ident).
#' @param group_by Column in meta.data slot to group results by (i.e. "Treatment").
#' @param colors_use List of colors or color palette to use.  Only applicable if `group_by` is not NULL.
#' @param dot_size size of the dots plotted if `group_by` is not NULL.  Default is 1.
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
#' @importFrom dplyr n select slice left_join any_of
#' @importFrom magrittr "%>%"
#'
#' @export
#'
#' @concept stats_plotting
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' # Create example groups
#' pbmc_small$sample_id <- sample(c("sample1", "sample2"), size = ncol(pbmc_small), replace = TRUE)
#'
#' # Plot
#' Plot_Median_Genes(seurat_object = pbmc_small, sample_col = "orig.ident",  group_by = "sample_id")
#'}
#'

Plot_Median_Genes <- function(
  seurat_object,
  sample_col = "orig.ident",
  group_by = NULL,
  colors_use = NULL,
  dot_size = 1,
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
  group_by <- Meta_Present(object = seurat_object, meta_col_names = group_by, print_msg = FALSE)[[1]]

  # Check sample_col is valid
  sample_col <- Meta_Present(object = seurat_object, meta_col_names = sample_col, print_msg = FALSE)[[1]]

  # Calculate medians and merge with meta.data
  medians <- Median_Stats(seurat_object = seurat_object, group_by_var = sample_col, median_var = "nFeature_RNA", default_var = FALSE) %>%
    slice(-n()) %>%
    droplevels()

  meta <- Fetch_Meta(object = seurat_object)

  if (!is.null(x = group_by)) {
    meta <- meta %>%
      select(any_of(c(sample_col, group_by)))
  } else {
    meta <- meta %>%
      select(any_of(sample_col))
  }

  meta[[sample_col]] <- factor(meta[[sample_col]], ordered = FALSE)

  meta <- data.frame(meta[!duplicated(x = meta[,sample_col]),])

  if (is.null(x = group_by)) {
    colnames(x = meta) <- sample_col
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

    plot <- ggplot(merged, aes(x = .data[["samples_plotting"]], y = .data[["Median_nFeature_RNA"]])) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("") +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(data = merged, mapping = aes(x = .data[[group_by]], y = .data[["Median_nFeature_RNA"]], fill = .data[[group_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center', dotsize = dot_size) +
      scale_fill_manual(values = colors_use) +
      theme_ggprism_mod() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("")
  }

  # Modify base plot
  if (isTRUE(x = x_lab_rotate)) {
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
#' @param dot_size size of the dots plotted if `group_by` is not NULL.  Default is 1.
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
#' @importFrom dplyr n select slice left_join any_of
#' @importFrom magrittr "%>%"
#'
#' @export
#'
#' @concept stats_plotting
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' # Create example groups
#' pbmc_small$sample_id <- sample(c("sample1", "sample2"), size = ncol(pbmc_small), replace = TRUE)
#'
#' # Plot
#' Plot_Median_UMIs(seurat_object = pbmc_small, sample_col = "orig.ident",  group_by = "sample_id")
#' }
#'

Plot_Median_UMIs <- function(
  seurat_object,
  sample_col = "orig.ident",
  group_by = NULL,
  colors_use = NULL,
  dot_size = 1,
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
  group_by <- Meta_Present(object = seurat_object, meta_col_names = group_by, print_msg = FALSE)[[1]]

  # Check sample_col is valid
  sample_col <- Meta_Present(object = seurat_object, meta_col_names = sample_col, print_msg = FALSE)[[1]]

  # Calculate medians and merge with meta.data
  medians <- Median_Stats(seurat_object = seurat_object, group_by_var = sample_col, median_var = "nCount_RNA", default_var = FALSE) %>%
    slice(-n()) %>%
    droplevels()

  meta <- Fetch_Meta(object = seurat_object)

  if (!is.null(x = group_by)) {
    meta <- meta %>%
      select(any_of(c(sample_col, group_by)))
  } else {
    meta <- meta %>%
      select(any_of(sample_col))
  }

  meta[[sample_col]] <- factor(meta[[sample_col]], ordered = FALSE)

  meta <- data.frame(meta[!duplicated(meta[,sample_col]),])

  if (is.null(x = group_by)) {
    colnames(x = meta) <- sample_col
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

    plot <- ggplot(merged, aes(x = .data[["samples_plotting"]], y = .data[["Median_nCount_RNA"]])) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("") +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(data = merged, mapping = aes(x = .data[[group_by]], y = .data[["Median_nCount_RNA"]], fill = .data[[group_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center', dotsize = dot_size) +
      scale_fill_manual(values = colors_use) +
      theme_ggprism_mod() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("")
  }

  # Modify base plot
  if (isTRUE(x = x_lab_rotate)) {
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
#' @param dot_size size of the dots plotted if `group_by` is not NULL.  Default is 1.
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
#' @importFrom dplyr n select slice left_join any_of
#' @importFrom magrittr "%>%"
#'
#' @export
#'
#' @concept stats_plotting
#'
#' @examples
#' \dontrun{
#' # Add mito
#' obj <- Add_Mito_Ribo_Seurat(seurat_object = obj, species = "human")
#'
#' # Plot
#' Plot_Median_Mito(seurat_object = obj, sample_col = "orig.ident",  group_by = "sample_id")
#'}
#'

Plot_Median_Mito <- function(
  seurat_object,
  sample_col = "orig.ident",
  group_by = NULL,
  colors_use = NULL,
  dot_size = 1,
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
  group_by <- Meta_Present(object = seurat_object, meta_col_names = group_by, print_msg = FALSE)[[1]]

  # Check sample_col is valid
  sample_col <- Meta_Present(object = seurat_object, meta_col_names = sample_col, print_msg = FALSE)[[1]]

  # Calculate medians and merge with meta.data
  medians <- Median_Stats(seurat_object = seurat_object, group_by_var = sample_col, median_var = "percent_mito", default_var = FALSE) %>%
    slice(-n()) %>%
    droplevels()

  meta <- Fetch_Meta(object = seurat_object)

  if (!is.null(x = group_by)) {
    meta <- meta %>%
      select(any_of(c(sample_col, group_by)))
  } else {
    meta <- meta %>%
      select(any_of(sample_col))
  }

  meta[[sample_col]] <- factor(meta[[sample_col]], ordered = FALSE)

  meta <- data.frame(meta[!duplicated(meta[,sample_col]),])

  if (is.null(x = group_by)) {
    colnames(x = meta) <- sample_col
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

    plot <- ggplot(merged, aes(x = .data[["samples_plotting"]], y = .data[["Median_percent_mito"]])) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("") +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(data = merged, mapping = aes(x = .data[[group_by]], y = .data[["Median_percent_mito"]], fill = .data[[group_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center', dotsize = dot_size) +
      scale_fill_manual(values = colors_use) +
      theme_ggprism_mod() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("")
  }

  # Modify base plot
  if (isTRUE(x = x_lab_rotate)) {
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
#' @param dot_size size of the dots plotted if `group_by` is not NULL.  Default is 1.
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
#' @importFrom dplyr n select slice left_join any_of
#' @importFrom magrittr "%>%"
#'
#' @export
#'
#' @concept stats_plotting
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' cd_features <- list(c('CD79B', 'CD79A', 'CD19', 'CD180', 'CD200', 'CD3D', 'CD2','CD3E',
#' 'CD7','CD8A', 'CD14', 'CD1C', 'CD68', 'CD9', 'CD247'))
#'
#' pbmc_small <- AddModuleScore(object = pbmc_small, features = cd_features, ctrl = 5,
#' name = 'CD_Features')
#'
#' Plot_Median_Other(seurat_object = pbmc_small, median_var = "CD_Features1",
#' sample_col = "orig.ident", group_by = "Treatment")
#' }
#'

Plot_Median_Other <- function(
  seurat_object,
  median_var,
  sample_col = "orig.ident",
  group_by = NULL,
  colors_use = NULL,
  dot_size = 1,
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
  group_by <- Meta_Present(object = seurat_object, meta_col_names = group_by, print_msg = FALSE)[[1]]

  # Check sample_col is valid
  sample_col <- Meta_Present(object = seurat_object, meta_col_names = sample_col, print_msg = FALSE)[[1]]

  # Calculate medians and merge with meta.data
  medians <- Median_Stats(seurat_object = seurat_object, group_by_var = sample_col, median_var = median_var, default_var = FALSE) %>%
    slice(-n()) %>%
    droplevels()

  meta <- Fetch_Meta(object = seurat_object)

  if (!is.null(x = group_by)) {
    meta <- meta %>%
      select(any_of(c(sample_col, group_by)))
  } else {
    meta <- meta %>%
      select(any_of(sample_col))
  }

  meta[[sample_col]] <- factor(meta[[sample_col]], ordered = FALSE)

  meta <- data.frame(meta[!duplicated(meta[,sample_col]),])

  if (is.null(x = group_by)) {
    colnames(x = meta) <- sample_col
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

    plot <- ggplot(merged, aes(x = .data[["samples_plotting"]], y = .data[[paste0("Median_", median_var)]])) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("") +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(data = merged, mapping = aes(x = .data[[group_by]], y = .data[[paste0("Median_", median_var)]], fill = .data[[group_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center', dotsize = dot_size) +
      scale_fill_manual(values = colors_use) +
      theme_ggprism_mod() +
      ggtitle(plot_title) +
      ylab(y_axis_label) +
      xlab("")
  }

  # Modify base plot
  if (isTRUE(x = x_lab_rotate)) {
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
#' @param dot_size size of the dots plotted if `group_by` is not NULL.  Default is 1.
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
#' @import cli
#' @import ggplot2
#' @importFrom rlang "%||%" ":="
#' @importFrom dplyr select slice left_join rename all_of
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
  dot_size = 1,
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
    cli_abort(message = "Must provide meta data variable to {.code group_by} in order to plot data.")
  }

  # Check group by is valid
  group_by <- Meta_Present(object = seurat_object, meta_col_names = group_by, print_msg = FALSE)[[1]]

  # Check sample_col is valid
  sample_col <- Meta_Present(object = seurat_object, meta_col_names = sample_col, print_msg = FALSE)[[1]]

  # Calculate total cells and merge with meta.data
  total_cells <- table(seurat_object@meta.data[[sample_col]]) %>%
    data.frame() %>%
    rename(!!sample_col := all_of("Var1"), Number_of_Cells = all_of("Freq"))

  meta <- Fetch_Meta(object = seurat_object)

  meta <- meta %>%
    select(all_of(c(sample_col, group_by)))

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
  plot <- ggplot(data = merged, mapping = aes(x = .data[[group_by]], y = .data[["Number_of_Cells"]], fill = .data[[group_by]])) +
    geom_boxplot(fill = "white") +
    geom_dotplot(binaxis ='y', stackdir = 'center', dotsize = dot_size) +
    scale_fill_manual(values = colors_use) +
    theme_ggprism_mod() +
    ggtitle(plot_title) +
    ylab(y_axis_label) +
    xlab("")

  # Modify base plot
  if (isTRUE(x = x_lab_rotate)) {
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
#' @param feature_diff_df name of data.frame created using \code{\link[scCustomize]{CellBender_Feature_Diff}}.
#' @param pct_diff_threshold threshold to use for feature plotting.  Resulting plot will only contain
#' features which exhibit percent change >= value.  Default is 25.
#' @param num_features Number of features to plot.  Will ignore `pct_diff_threshold` and return
#' plot with specified number of features.  Default is NULL.
#' @param label logical, whether or not to label the features that have largest percent difference
#' between raw and CellBender counts (Default is TRUE).
#' @param num_labels Number of features to label if `label = TRUE`, (default is 20).
#' @param min_count_label Minimum number of raw counts per feature necessary to be included in
#' plot labels (default is 1)
#' @param repel logical, whether to use geom_text_repel to create a nicely-repelled labels; this is
#' slow when a lot of points are being plotted. If using repel, set xnudge and ynudge to 0, (Default is TRUE).
#' @param custom_labels A custom set of features to label instead of the features most different between
#' raw and CellBender counts.
#' @param plot_line logical, whether to plot diagonal line with slope = 1 (Default is TRUE).
#' @param plot_title Plot title.
#' @param x_axis_label Label for x axis.
#' @param y_axis_label Label for y axis.
#' @param xnudge Amount to nudge X and Y coordinates of labels by.
#' @param ynudge Amount to nudge X and Y coordinates of labels by.
#' @param max.overlaps passed to \code{\link[ggrepel]{geom_text_repel}}, exclude text labels that
#' overlap too many things. Defaults to 100.
#' @param label_color Color to use for text labels.
#' @param fontface font face to use for text labels (“plain”, “bold”, “italic”, “bold.italic”) (Default is "bold").
#' @param label_size text size for feature labels (passed to \code{\link[ggrepel]{geom_text_repel}}).
#' @param bg.color color to use for shadow/outline of text labels (passed to \code{\link[ggrepel]{geom_text_repel}}) (Default is white).
#' @param bg.r radius to use for shadow/outline of text labels (passed to \code{\link[ggrepel]{geom_text_repel}}) (Default is 0.15).
#' @param ... Extra parameters passed to \code{\link[ggrepel]{geom_text_repel}} through
#' \code{\link[Seurat]{LabelPoints}}.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom cowplot theme_cowplot
#' @importFrom dplyr filter
#' @importFrom magrittr "%>%"
#' @importFrom tidyr drop_na
#'
#' @export
#'
#' @concept stats_plotting
#'
#' @examples
#' \dontrun{
#' # get cell bender differences data.frame
#' cb_stats <- CellBender_Feature_Diff(seurat_object - obj, raw_assay = "RAW",
#' cell_bender_assay = "RNA")
#'
#' # plot
#' CellBender_Diff_Plot(feature_diff_df = cb_stats, pct_diff_threshold = 25)
#' }
#'

CellBender_Diff_Plot <- function(
  feature_diff_df,
  pct_diff_threshold = 25,
  num_features = NULL,
  label = TRUE,
  num_labels = 20,
  min_count_label = 1,
  repel = TRUE,
  custom_labels = NULL,
  plot_line = TRUE,
  plot_title = "Raw Counts vs. Cell Bender Counts",
  x_axis_label = "Raw Data Counts",
  y_axis_label = "Cell Bender Counts",
  xnudge = 0,
  ynudge = 0,
  max.overlaps = 100,
  label_color = "dodgerblue",
  fontface = "bold",
  label_size = 3.88,
  bg.color = "white",
  bg.r = 0.15,
  ...
) {
  # Remove unshared features
  feature_diff_df_filtered <- feature_diff_df %>%
    drop_na(all_of(c("Raw_Counts", "CellBender_Counts")))

  diff_features <- symdiff(x = rownames(x = feature_diff_df), y = rownames(x = feature_diff_df_filtered))

  if (length(x = diff_features > 0)) {
    cli_warn(message = c("The following features are not present in both assays and were omitted:",
                         "*" = "{.field diff_features}}")
    )
  }

  num_features_total <- nrow(x = feature_diff_df_filtered)

  # Check how to filter data.frame
  if (!is.null(x = pct_diff_threshold) && !is.null(x = num_features)) {
    cli_abort(message = c("{.code pct_diff_threshold} and {.code num_features} cannot both have values.",
                          "i" = "Set undesired parameter to NULL."))
  }

  # Filter plot
  if (!is.null(x = pct_diff_threshold)) {
    feature_diff_df_filtered <- feature_diff_df_filtered %>%
      filter(.data[["Pct_Diff"]] >= pct_diff_threshold)
  } else {
    feature_diff_df_filtered <- feature_diff_df_filtered[1:num_features, ]
  }

  num_features_plotted <- nrow(x = feature_diff_df_filtered)

  # Extract max plotted value
  axis_lim <- max(feature_diff_df_filtered$Raw_Counts)

  # Make plot
  plot <- ggplot(feature_diff_df_filtered, aes(x = .data[["Raw_Counts"]], y = .data[["CellBender_Counts"]])) +
    geom_point() +
    scale_x_log10(limits = c(1, axis_lim)) +
    scale_y_log10(limits = c(1, axis_lim)) +
    ylab(y_axis_label) +
    xlab(x_axis_label) +
    theme_cowplot() +
    if (!is.null(x = pct_diff_threshold)) {
      ggtitle(plot_title, subtitle = paste0("Plotting features which exhibit difference of ", pct_diff_threshold, "% or greater (", num_features_plotted, "/", num_features_total, ")." ))
    } else {
      ggtitle(plot_title, subtitle = paste0("Plotting ", num_features_plotted, "/", num_features_total, " features." ))
    }

  # Label points
  if (isTRUE(x = label)) {
    if (is.null(x = custom_labels)) {
      # Subset the labels based on min count threshold
      labels_use <- feature_diff_df_filtered %>%
        filter(.data[["Raw_Counts"]] >= min_count_label) %>%
        rownames()

      # Return message of features not found
      if (length(x = labels_use) == 0) {
        cli_warn(message = c("No features met the labeling criteria.",
                             "i" = "Try adjusting {.field min_count_label} and/or {.field pct_diff_threshold}.")
        )

        plot <- plot
      } else {
        plot <- LabelPoints(plot = plot, points = labels_use[1:num_labels], repel = repel, xnudge = xnudge, ynudge = ynudge, max.overlaps = max.overlaps, color = label_color, fontface = fontface, size = label_size, bg.color = bg.color, bg.r = bg.r, ...)
      }
    } else {
      # check for features
      features_list <- Feature_Present(data = feature_diff_df_filtered, features = custom_labels, omit_warn = FALSE, print_msg = FALSE, case_check_msg = FALSE, return_none = TRUE)

      all_not_found_features <- features_list[[2]]

      all_found_features <- features_list[[1]]

      # Stop if no features found
      if (length(x = all_found_features) < 1) {
        cli_abort(message = c("None of features in {.code custom_labels} were found in plot data.",
                              "i" = "Check both raw data and adjust {.code pct_diff_threshold} if needed.")
        )
      }

      # Return message of features not found
      if (length(x = all_not_found_features) > 0) {
        op <- options(warn = 1)
        on.exit(options(op))
        cli_warn(message = c("The following features in {.code custom_labels} were omitted as they were not found:",
                             "*" = "{.field {glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}}",
                             "i" = "Check both raw data and adjust {.code pct_diff_threshold} if needed.")
        )
      }
      # plot with custom labels
      plot <- LabelPoints(plot = plot, points = all_found_features, repel = repel, xnudge = xnudge, ynudge = ynudge, max.overlaps = max.overlaps, color = label_color, fontface = fontface, size = label_size, bg.color = bg.color, bg.r = bg.r, ...)
    }
  }

  if (isTRUE(x = plot_line)) {
    plot <- plot + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
  }

  # return plot
  return(plot)
}
