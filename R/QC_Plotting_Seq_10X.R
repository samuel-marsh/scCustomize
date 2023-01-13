#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### 10X SEQ QC ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' QC Plots Sequencing metrics
#'
#' Plot the mean number of reads per cell
#'
#' @param metrics_dataframe data.frame  contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#'  less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom SeuratObject PackageCheck
#' @importFrom utils combn
#'
#' @export
#'
#' @concept seq_qc_plotting_basic
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_Reads_per_Cell(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_Reads_per_Cell <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  significance = FALSE,
  ...
) {
  if (!plot_by %in% colnames(x = metrics_dataframe)) {
    cli_abort(message = "{.val {plot_by}} is not a column in the provided {.code metrics_dataframe}.")
  }

  # Change plot_by to character vector to make significance functions show all comparisons
  if (inherits(x = metrics_dataframe[[plot_by]], what = "factor")) {
    stats_dataframe <- metrics_dataframe
    stats_dataframe[[plot_by]] <- as.character(stats_dataframe[[plot_by]])
  } else {
    stats_dataframe <- metrics_dataframe
  }

  # Create color palette if null and check valid if provided
  length_plotby <- length(x = unique(x = metrics_dataframe[[plot_by]]))

  if (is.null(x = colors_use) && !plot_by == "sample_id") {
    if (length_plotby <= 8) {
      colors_use <- Dark2_Pal()
    } else {
      colors_use <- DiscretePalette_scCustomize(num_colors = length_plotby, palette = "polychrome")
    }
  } else {
    if (length(x = colors_use) < length_plotby && !plot_by == "sample_id") {
      cli_abort(message = c("Not enough colors provided.",
                            "i" = "The number of colors supplied: {.field {length(x = colors_use)}}, is less than the number of groups in {.val {plot_by}} column: {.field {length_plotby}}.")
      )
    } else {
      colors_use <- colors_use
    }
  }

  if (plot_by == "sample_id") {
    metrics_dataframe$samples_plotting <- "Samples"

    plot <- ggplot(metrics_dataframe, aes(x = .data[["samples_plotting"]], y = .data[["Mean_Reads_per_Cell"]])) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(face = "bold", size = 14),
            plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
      ggtitle("Mean Reads per Cell per Sample") +
      ylab('Mean Reads per Cell') +
      xlab("") +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(metrics_dataframe, aes(x=.data[[plot_by]], y = .data[["Mean_Reads_per_Cell"]], fill = .data[[plot_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(face = "bold", size = 14),
            plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
      scale_fill_manual(values = colors_use) +
      ggtitle("Mean Reads per Cell per Sample") +
      ylab('Mean Reads per Cell') +
      xlab(plot_by) +
      theme_ggprism_mod()
  }

  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (significance) {
    ggpubr_check <- PackageCheck("ggpubr", error = FALSE)
    if (!ggpubr_check[1]) {
      cli_abort(message = c(
        "Please install the {.val ggpubr} package to calculate/plot significance values.",
        "i" = "This can be accomplished with the following commands: ",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}ggpubr{symbol$dquote_right})`}",
        "----------------------------------------"
      ))
    }

    if (length(x = unique(x = stats_dataframe[[plot_by]])) < 2) {
      cli_abort(message = "Cannot calculate statistics when {.val {plot_by}} column contains less than 2 groups.")
    }
    groups <- unique(x = stats_dataframe[[plot_by]])

    comparisons <- combn(groups, 2)
    comparisons <- data.frame(comparisons, stringsAsFactors = FALSE)
    comparisons <- lapply(1:length(x = colnames(x = comparisons)), function(x){
      indiv_comp <- as.character(x = comparisons[[x]])
    })
    plot <- plot + ggpubr::stat_compare_means(comparisons = comparisons, ...)
  }

  return(plot)
}


#' QC Plots Sequencing metrics
#'
#' Plot the number of cells per sample
#'
#' @param metrics_dataframe data.frame  contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#'  less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom SeuratObject PackageCheck
#' @importFrom utils combn
#'
#' @export
#'
#' @concept seq_qc_plotting_basic
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_Number_Cells(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_Number_Cells <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  significance = FALSE,
  ...
) {
  if (!plot_by %in% colnames(x = metrics_dataframe)) {
    cli_abort(message = "{.val {plot_by}} is not a column in the provided {.code metrics_dataframe}.")
  }

  # Change plot_by to character vector to make significance functions show all comparisons
  if (inherits(x = metrics_dataframe[[plot_by]], what = "factor")) {
    stats_dataframe <- metrics_dataframe
    stats_dataframe[[plot_by]] <- as.character(stats_dataframe[[plot_by]])
  } else {
    stats_dataframe <- metrics_dataframe
  }

  # Create color palette if null and check valid if provided
  length_plotby <- length(x = unique(x = metrics_dataframe[[plot_by]]))

  if (is.null(x = colors_use) && !plot_by == "sample_id") {
    if (length_plotby <= 8) {
      colors_use <- Dark2_Pal()
    } else {
      colors_use <- DiscretePalette_scCustomize(num_colors = length_plotby, palette = "polychrome")
    }
  } else {
    if (length(x = colors_use) < length_plotby && !plot_by == "sample_id") {
      cli_abort(message = c("Not enough colors provided.",
                            "i" = "The number of colors supplied: {.field {length(x = colors_use)}}, is less than the number of groups in {.val {plot_by}} column: {.field {length_plotby}}.")
      )
    } else {
      colors_use <- colors_use
    }
  }

  if (plot_by == "sample_id") {
    metrics_dataframe$samples_plotting <- "Samples"

    plot <- ggplot(metrics_dataframe, aes(x = .data[["samples_plotting"]], y = .data[["Estimated_Number_of_Cells"]])) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(face = "bold", size = 14),
            plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
      ggtitle("Cells per Sample") +
      ylab('Cells') +
      xlab("") +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(metrics_dataframe, aes(x=.data[[plot_by]], y = .data[["Estimated_Number_of_Cells"]], fill = .data[[plot_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(face = "bold", size = 14),
            plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
      scale_fill_manual(values = colors_use) +
      ggtitle("Cells per Sample") +
      ylab('Cells') +
      xlab(plot_by) +
      theme_ggprism_mod()
  }

  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (significance) {
    ggpubr_check <- PackageCheck("ggpubr", error = FALSE)
    if (!ggpubr_check[1]) {
      cli_abort(message = c(
        "Please install the {.val ggpubr} package to calculate/plot significance values.",
        "i" = "This can be accomplished with the following commands: ",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}ggpubr{symbol$dquote_right})`}",
        "----------------------------------------"
      ))
    }

    if (length(x = unique(x = stats_dataframe[[plot_by]])) < 2) {
      cli_abort(message = "Cannot calculate statistics when {.val {plot_by}} column contains less than 2 groups.")
    }
    groups <- unique(x = stats_dataframe[[plot_by]])

    comparisons <- combn(groups, 2)
    comparisons <- data.frame(comparisons, stringsAsFactors = FALSE)
    comparisons <- lapply(1:length(x = colnames(x = comparisons)), function(x){
      indiv_comp <- as.character(x = comparisons[[x]])
    })
    plot <- plot + ggpubr::stat_compare_means(comparisons = comparisons, ...)
  }

  return(plot)
}


#' QC Plots Sequencing metrics
#'
#' Plot the median genes per cell per sample
#'
#' @param metrics_dataframe data.frame  contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#' less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom SeuratObject PackageCheck
#' @importFrom utils combn
#'
#' @export
#'
#' @concept seq_qc_plotting_basic
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_Genes(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_Genes <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  significance = FALSE,
  ...
) {
  if (!plot_by %in% colnames(x = metrics_dataframe)) {
    cli_abort(message = " is not a column in the provided `metrics_dataframe`.")
  }

  # Change plot_by to character vector to make significance functions show all comparisons
  if (inherits(x = metrics_dataframe[[plot_by]], what = "factor")) {
    stats_dataframe <- metrics_dataframe
    stats_dataframe[[plot_by]] <- as.character(stats_dataframe[[plot_by]])
  } else {
    stats_dataframe <- metrics_dataframe
  }

  # Create color palette if null and check valid if provided
  length_plotby <- length(x = unique(x = metrics_dataframe[[plot_by]]))

  if (is.null(x = colors_use) && !plot_by == "sample_id") {
    if (length_plotby <= 8) {
      colors_use <- Dark2_Pal()
    } else {
      colors_use <- DiscretePalette_scCustomize(num_colors = length_plotby, palette = "polychrome")
    }
  } else {
    if (length(x = colors_use) < length_plotby && !plot_by == "sample_id") {
      cli_abort(message = c("Not enough colors provided.",
                            "i" = "The number of colors supplied: {.field {length(x = colors_use)}}, is less than the number of groups in {.val {plot_by}} column: {.field {length_plotby}}.")
      )
    } else {
      colors_use <- colors_use
    }
  }

  if (plot_by == "sample_id") {
    metrics_dataframe$samples_plotting <- "Samples"

    plot <- ggplot(metrics_dataframe, aes(x = .data[["samples_plotting"]], y = .data[["Median_Genes_per_Cell"]])) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle("Median Genes per Cell") +
      ylab('Median Genes') +
      xlab("") +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(metrics_dataframe, aes(x=.data[[plot_by]], y = .data[["Median_Genes_per_Cell"]], fill = .data[[plot_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      ggtitle("Median Genes per Cell") +
      ylab('Median Genes') +
      xlab(plot_by) +
      theme_ggprism_mod()
  }

  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (significance) {
    ggpubr_check <- PackageCheck("ggpubr", error = FALSE)
    if (!ggpubr_check[1]) {
      cli_abort(message = c(
        "Please install the {.val ggpubr} package to calculate/plot significance values.",
        "i" = "This can be accomplished with the following commands: ",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}ggpubr{symbol$dquote_right})`}",
        "----------------------------------------"
      ))
    }

    if (length(x = unique(x = stats_dataframe[[plot_by]])) < 2) {
      cli_abort(message = "Cannot calculate statistics when {.val {plot_by}} column contains less than 2 groups.")
    }
    groups <- unique(x = stats_dataframe[[plot_by]])

    comparisons <- combn(groups, 2)
    comparisons <- data.frame(comparisons, stringsAsFactors = FALSE)
    comparisons <- lapply(1:length(x = colnames(x = comparisons)), function(x){
      indiv_comp <- as.character(x = comparisons[[x]])
    })
    plot <- plot + ggpubr::stat_compare_means(comparisons = comparisons, ...)
  }

  return(plot)
}


#' QC Plots Sequencing metrics
#'
#' Plot the median UMIs per cell per sample
#'
#' @param metrics_dataframe data.frame contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#' less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom SeuratObject PackageCheck
#' @importFrom utils combn
#'
#' @export
#'
#' @concept seq_qc_plotting_basic
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_UMIs(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_UMIs <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  significance = FALSE,
  ...
) {
  if (!plot_by %in% colnames(x = metrics_dataframe)) {
    cli_abort(message = "{.val {plot_by}} is not a column in the provided {.code metrics_dataframe}.")
  }

  # Change plot_by to character vector to make significance functions show all comparisons
  if (inherits(x = metrics_dataframe[[plot_by]], what = "factor")) {
    stats_dataframe <- metrics_dataframe
    stats_dataframe[[plot_by]] <- as.character(stats_dataframe[[plot_by]])
  } else {
    stats_dataframe <- metrics_dataframe
  }

  # Create color palette if null and check valid if provided
  length_plotby <- length(x = unique(x = metrics_dataframe[[plot_by]]))

  if (is.null(x = colors_use) && !plot_by == "sample_id") {
    if (length_plotby <= 8) {
      colors_use <- Dark2_Pal()
    } else {
      colors_use <- DiscretePalette_scCustomize(num_colors = length_plotby, palette = "polychrome")
    }
  } else {
    if (length(x = colors_use) < length_plotby && !plot_by == "sample_id") {
      cli_abort(message = c("Not enough colors provided.",
                            "i" = "The number of colors supplied: {.field {length(x = colors_use)}}, is less than the number of groups in {.val {plot_by}} column: {.field {length_plotby}}.")
      )
    } else {
      colors_use <- colors_use
    }
  }

  if (plot_by == "sample_id") {
    metrics_dataframe$samples_plotting <- "Samples"

    plot <- ggplot(metrics_dataframe, aes(x = .data[["samples_plotting"]], y = .data[["Median_UMI_Counts_per_Cell"]])) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle("Median UMIs per Cell") +
      ylab('Median UMIs') +
      xlab("") +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(metrics_dataframe, aes(x=.data[[plot_by]], y = .data[["Median_UMI_Counts_per_Cell"]], fill = .data[[plot_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      ggtitle("Median UMIs per Cell") +
      ylab('Median UMIs') +
      xlab(plot_by) +
      theme_ggprism_mod()
  }

  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (significance) {
    ggpubr_check <- PackageCheck("ggpubr", error = FALSE)
    if (!ggpubr_check[1]) {
      cli_abort(message = c(
        "Please install the {.val ggpubr} package to calculate/plot significance values.",
        "i" = "This can be accomplished with the following commands: ",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}ggpubr{symbol$dquote_right})`}",
        "----------------------------------------"
      ))
    }

    if (length(x = unique(x = stats_dataframe[[plot_by]])) < 2) {
      cli_abort(message = "Cannot calculate statistics when {.val {plot_by}} column contains less than 2 groups.")
    }
    groups <- unique(x = stats_dataframe[[plot_by]])

    comparisons <- combn(groups, 2)
    comparisons <- data.frame(comparisons, stringsAsFactors = FALSE)
    comparisons <- lapply(1:length(x = colnames(x = comparisons)), function(x){
      indiv_comp <- as.character(x = comparisons[[x]])
    })
    plot <- plot + ggpubr::stat_compare_means(comparisons = comparisons, ...)
  }

  return(plot)
}


#' QC Plots Sequencing metrics
#'
#' Plot the total genes detected per sample
#'
#' @param metrics_dataframe data.frame  contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#' less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom SeuratObject PackageCheck
#' @importFrom utils combn
#'
#' @export
#'
#' @concept seq_qc_plotting_basic
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_Total_Genes(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_Total_Genes <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  significance = FALSE,
  ...
) {
  if (!plot_by %in% colnames(x = metrics_dataframe)) {
    cli_abort(message = "{.val {plot_by}} is not a column in the provided {.code metrics_dataframe}.")
  }

  # Change plot_by to character vector to make significance functions show all comparisons
  if (inherits(x = metrics_dataframe[[plot_by]], what = "factor")) {
    stats_dataframe <- metrics_dataframe
    stats_dataframe[[plot_by]] <- as.character(stats_dataframe[[plot_by]])
  } else {
    stats_dataframe <- metrics_dataframe
  }

  # Create color palette if null and check valid if provided
  length_plotby <- length(x = unique(x = metrics_dataframe[[plot_by]]))

  if (is.null(x = colors_use) && !plot_by == "sample_id") {
    if (length_plotby <= 8) {
      colors_use <- Dark2_Pal()
    } else {
      colors_use <- DiscretePalette_scCustomize(num_colors = length_plotby, palette = "polychrome")
    }
  } else {
    if (length(x = colors_use) < length_plotby && !plot_by == "sample_id") {
      cli_abort(message = c("Not enough colors provided.",
                            "i" = "The number of colors supplied: {.field {length(x = colors_use)}}, is less than the number of groups in {.val {plot_by}} column: {.field {length_plotby}}.")
      )
    } else {
      colors_use <- colors_use
    }
  }

  if (plot_by == "sample_id") {
    metrics_dataframe$samples_plotting <- "Samples"

    plot <- ggplot(metrics_dataframe, aes(x = .data[["samples_plotting"]], y = .data[["Total_Genes_Detected"]])) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle("Total Genes Detected per Sample") +
      ylab('Total Genes') +
      xlab("") +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(metrics_dataframe, aes(x=.data[[plot_by]], y = .data[["Total_Genes_Detected"]], fill = .data[[plot_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      ggtitle("Total Genes Detected per Sample") +
      ylab('Total Genes') +
      xlab(plot_by) +
      theme_ggprism_mod()
  }

  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (significance) {
    ggpubr_check <- PackageCheck("ggpubr", error = FALSE)
    if (!ggpubr_check[1]) {
      cli_abort(message = c(
        "Please install the {.val ggpubr} package to calculate/plot significance values.",
        "i" = "This can be accomplished with the following commands: ",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}ggpubr{symbol$dquote_right})`}",
        "----------------------------------------"
        ))
    }

    if (length(x = unique(x = stats_dataframe[[plot_by]])) < 2) {
      cli_abort(message = "Cannot calculate statistics when {.val {plot_by}} column contains less than 2 groups.")
    }
    groups <- unique(x = stats_dataframe[[plot_by]])

    comparisons <- combn(groups, 2)
    comparisons <- data.frame(comparisons, stringsAsFactors = FALSE)
    comparisons <- lapply(1:length(x = colnames(x = comparisons)), function(x){
      indiv_comp <- as.character(x = comparisons[[x]])
    })
    plot <- plot + ggpubr::stat_compare_means(comparisons = comparisons, ...)
  }

  return(plot)
}


#' QC Plots Sequencing metrics
#'
#' Plot the sequencing saturation percentage per sample
#'
#' @param metrics_dataframe data.frame  contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#' less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom scales label_percent
#' @importFrom SeuratObject PackageCheck
#' @importFrom utils combn
#'
#' @export
#'
#' @concept seq_qc_plotting_basic
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_Saturation(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_Saturation <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  significance = FALSE,
  ...
) {
  if (!plot_by %in% colnames(x = metrics_dataframe)) {
    cli_abort(message = "{.val {plot_by}} is not a column in the provided {.code metrics_dataframe}.")
  }

  # Change plot_by to character vector to make significance functions show all comparisons
  if (inherits(x = metrics_dataframe[[plot_by]], what = "factor")) {
    stats_dataframe <- metrics_dataframe
    stats_dataframe[[plot_by]] <- as.character(stats_dataframe[[plot_by]])
  } else {
    stats_dataframe <- metrics_dataframe
  }

  # Create color palette if null and check valid if provided
  length_plotby <- length(x = unique(x = metrics_dataframe[[plot_by]]))

  if (is.null(x = colors_use) && !plot_by == "sample_id") {
    if (length_plotby <= 8) {
      colors_use <- Dark2_Pal()
    } else {
      colors_use <- DiscretePalette_scCustomize(num_colors = length_plotby, palette = "polychrome")
    }
  } else {
    if (length(x = colors_use) < length_plotby && !plot_by == "sample_id") {
      cli_abort(message = c("Not enough colors provided.",
                            "i" = "The number of colors supplied: {.field {length(x = colors_use)}}, is less than the number of groups in {.val {plot_by}} column: {.field {length_plotby}}.")
      )
    } else {
      colors_use <- colors_use
    }
  }

  # Modify dataframe
  metrics_dataframe[,"Sequencing_Saturation"] <- as.numeric(gsub("%", "", metrics_dataframe[,"Sequencing_Saturation"]))


  if (plot_by == "sample_id") {
    metrics_dataframe$samples_plotting <- "Samples"

    plot <- ggplot(metrics_dataframe, aes(x = .data[["samples_plotting"]], y = .data[["Sequencing_Saturation"]]), color = .data[["samples_plotting"]]) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle("Sequencing Saturation") +
      ylab('Sequencing Saturation Percent') +
      xlab("") +
      scale_y_continuous(labels = label_percent(accuracy = 0.01, scale = 1)) +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(metrics_dataframe, aes(x=.data[[plot_by]], y = .data[["Sequencing_Saturation"]], fill = .data[[plot_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      ggtitle("Sequencing Saturation") +
      ylab('Sequencing Saturation Percent') +
      xlab(plot_by)+
      scale_y_continuous(labels = label_percent(accuracy = 1, scale = 1)) +
      theme_ggprism_mod()
  }

  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (significance) {
    ggpubr_check <- PackageCheck("ggpubr", error = FALSE)
    if (!ggpubr_check[1]) {
      cli_abort(message = c(
        "Please install the {.val ggpubr} package to calculate/plot significance values.",
        "i" = "This can be accomplished with the following commands: ",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}ggpubr{symbol$dquote_right})`}",
        "----------------------------------------"
      ))
    }

    if (length(x = unique(x = stats_dataframe[[plot_by]])) < 2) {
      cli_abort(message = "Cannot calculate statistics when {.val {plot_by}} column contains less than 2 groups.")
    }
    groups <- unique(x = stats_dataframe[[plot_by]])

    comparisons <- combn(groups, 2)
    comparisons <- data.frame(comparisons, stringsAsFactors = FALSE)
    comparisons <- lapply(1:length(x = colnames(x = comparisons)), function(x){
      indiv_comp <- as.character(x = comparisons[[x]])
    })
    plot <- plot + ggpubr::stat_compare_means(comparisons = comparisons, ...)
  }

  return(plot)
}


#' QC Plots Sequencing metrics
#'
#' Plot the fraction of reads in cells per sample
#'
#' @param metrics_dataframe data.frame  contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#' less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom scales label_percent
#' @importFrom SeuratObject PackageCheck
#' @importFrom utils combn
#'
#' @export
#'
#' @concept seq_qc_plotting_basic
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_Reads_in_Cells(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_Reads_in_Cells <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  significance = FALSE,
  ...
) {
  if (!plot_by %in% colnames(x = metrics_dataframe)) {
    cli_abort(message = "{.val {plot_by}} is not a column in the provided {.code metrics_dataframe}.")
  }

  # Change plot_by to character vector to make significance functions show all comparisons
  if (inherits(x = metrics_dataframe[[plot_by]], what = "factor")) {
    stats_dataframe <- metrics_dataframe
    stats_dataframe[[plot_by]] <- as.character(stats_dataframe[[plot_by]])
  } else {
    stats_dataframe <- metrics_dataframe
  }

  # Create color palette if null and check valid if provided
  length_plotby <- length(x = unique(x = metrics_dataframe[[plot_by]]))

  if (is.null(x = colors_use) && !plot_by == "sample_id") {
    if (length_plotby <= 8) {
      colors_use <- Dark2_Pal()
    } else {
      colors_use <- DiscretePalette_scCustomize(num_colors = length_plotby, palette = "polychrome")
    }
  } else {
    if (length(x = colors_use) < length_plotby && !plot_by == "sample_id") {
      cli_abort(message = c("Not enough colors provided.",
                            "i" = "The number of colors supplied: {.field {length(x = colors_use)}}, is less than the number of groups in {.val {plot_by}} column: {.field {length_plotby}}.")
      )
    } else {
      colors_use <- colors_use
    }
  }

  # Modify dataframe
  metrics_dataframe[,"Fraction_Reads_in_Cells"] <- as.numeric(gsub("%", "", metrics_dataframe[,"Fraction_Reads_in_Cells"]))


  if (plot_by == "sample_id") {
    metrics_dataframe$samples_plotting <- "Samples"

    plot <- ggplot(metrics_dataframe, aes(x = .data[["samples_plotting"]], y = .data[["Fraction_Reads_in_Cells"]]), color = .data[["samples_plotting"]]) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle("Fraction of Reads in Cells per Sample") +
      ylab('Fraction of Reads in Cells') +
      xlab("") +
      scale_y_continuous(labels = label_percent(accuracy = 0.01, scale = 1)) +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(metrics_dataframe, aes(x=.data[[plot_by]], y = .data[["Fraction_Reads_in_Cells"]], fill = .data[[plot_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      ggtitle("Fraction of Reads in Cells per Sample") +
      ylab('Fraction of Reads in Cells') +
      xlab(plot_by)+
      scale_y_continuous(labels = label_percent(accuracy = 1, scale = 1)) +
      theme_ggprism_mod()
  }

  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (significance) {
    ggpubr_check <- PackageCheck("ggpubr", error = FALSE)
    if (!ggpubr_check[1]) {
      cli_abort(message = c(
        "Please install the {.val ggpubr} package to calculate/plot significance values.",
        "i" = "This can be accomplished with the following commands: ",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}ggpubr{symbol$dquote_right})`}",
        "----------------------------------------"
      ))
    }

    if (length(x = unique(x = stats_dataframe[[plot_by]])) < 2) {
      cli_abort(message = "Cannot calculate statistics when {.val {plot_by}} column contains less than 2 groups.")
    }
    groups <- unique(x = stats_dataframe[[plot_by]])

    comparisons <- combn(groups, 2)
    comparisons <- data.frame(comparisons, stringsAsFactors = FALSE)
    comparisons <- lapply(1:length(x = colnames(x = comparisons)), function(x){
      indiv_comp <- as.character(x = comparisons[[x]])
    })
    plot <- plot + ggpubr::stat_compare_means(comparisons = comparisons, ...)
  }

  return(plot)
}


#' QC Plots Sequencing metrics (Alignment)
#'
#' Plot the fraction of reads confidently mapped to transcriptome
#'
#' @param metrics_dataframe data.frame  contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#' less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom scales label_percent
#' @importFrom SeuratObject PackageCheck
#' @importFrom utils combn
#'
#' @export
#'
#' @concept seq_qc_plotting_alignment
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_Transcriptome(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_Transcriptome <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  significance = FALSE,
  ...
) {
  if (!plot_by %in% colnames(x = metrics_dataframe)) {
    cli_abort(message = "{.val {plot_by}} is not a column in the provided {.code metrics_dataframe}.")
  }

  # Change plot_by to character vector to make significance functions show all comparisons
  if (inherits(x = metrics_dataframe[[plot_by]], what = "factor")) {
    stats_dataframe <- metrics_dataframe
    stats_dataframe[[plot_by]] <- as.character(stats_dataframe[[plot_by]])
  } else {
    stats_dataframe <- metrics_dataframe
  }

  # Create color palette if null and check valid if provided
  length_plotby <- length(x = unique(x = metrics_dataframe[[plot_by]]))

  if (is.null(x = colors_use) && !plot_by == "sample_id") {
    if (length_plotby <= 8) {
      colors_use <- Dark2_Pal()
    } else {
      colors_use <- DiscretePalette_scCustomize(num_colors = length_plotby, palette = "polychrome")
    }
  } else {
    if (length(x = colors_use) < length_plotby && !plot_by == "sample_id") {
      cli_abort(message = c("Not enough colors provided.",
                            "i" = "The number of colors supplied: {.field {length(x = colors_use)}}, is less than the number of groups in {.val {plot_by}} column: {.field {length_plotby}}.")
      )
    } else {
      colors_use <- colors_use
    }
  }

  # Modify dataframe
  metrics_dataframe[,"Reads_Mapped_Confidently_to_Transcriptome"] <- as.numeric(gsub("%", "", metrics_dataframe[,"Reads_Mapped_Confidently_to_Transcriptome"]))


  if (plot_by == "sample_id") {
    metrics_dataframe$samples_plotting <- "Samples"

    plot <- ggplot(metrics_dataframe, aes(x = .data[["samples_plotting"]], y = .data[["Reads_Mapped_Confidently_to_Transcriptome"]]), color = .data[["samples_plotting"]]) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle("Percent of Reads Confidently Mapped to Transcriptome") +
      ylab('Percent of Reads ') +
      xlab("") +
      scale_y_continuous(labels = label_percent(accuracy = 0.01, scale = 1)) +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(metrics_dataframe, aes(x=.data[[plot_by]], y = .data[["Reads_Mapped_Confidently_to_Transcriptome"]], fill = .data[[plot_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      ggtitle("Percent of Reads Confidently Mapped to Transcriptome") +
      ylab('Percent of Reads') +
      xlab(plot_by) +
      scale_y_continuous(labels = label_percent(accuracy = 1, scale = 1)) +
      theme_ggprism_mod()
  }
  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (significance) {
    ggpubr_check <- PackageCheck("ggpubr", error = FALSE)
    if (!ggpubr_check[1]) {
      cli_abort(message = c(
        "Please install the {.val ggpubr} package to calculate/plot significance values.",
        "i" = "This can be accomplished with the following commands: ",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}ggpubr{symbol$dquote_right})`}",
        "----------------------------------------"
      ))
    }

    if (length(x = unique(x = stats_dataframe[[plot_by]])) < 2) {
      cli_abort(message = "Cannot calculate statistics when {.val {plot_by}} column contains less than 2 groups.")
    }
    groups <- unique(x = stats_dataframe[[plot_by]])

    comparisons <- combn(groups, 2)
    comparisons <- data.frame(comparisons, stringsAsFactors = FALSE)
    comparisons <- lapply(1:length(x = colnames(x = comparisons)), function(x){
      indiv_comp <- as.character(x = comparisons[[x]])
    })
    plot <- plot + ggpubr::stat_compare_means(comparisons = comparisons, ...)
  }

  return(plot)
}


#' QC Plots Sequencing metrics (Alignment)
#'
#' Plot the fraction of reads confidently mapped to genome
#'
#' @param metrics_dataframe data.frame  contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#' less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom scales label_percent
#' @importFrom SeuratObject PackageCheck
#' @importFrom utils combn
#'
#' @export
#'
#' @concept seq_qc_plotting_alignment
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_Genome(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_Genome <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  significance = FALSE,
  ...
) {
  if (!plot_by %in% colnames(x = metrics_dataframe)) {
    cli_abort(message = "{.val {plot_by}} is not a column in the provided {.code metrics_dataframe}.")
  }

  # Change plot_by to character vector to make significance functions show all comparisons
  if (inherits(x = metrics_dataframe[[plot_by]], what = "factor")) {
    stats_dataframe <- metrics_dataframe
    stats_dataframe[[plot_by]] <- as.character(stats_dataframe[[plot_by]])
  } else {
    stats_dataframe <- metrics_dataframe
  }

  # Create color palette if null and check valid if provided
  length_plotby <- length(x = unique(x = metrics_dataframe[[plot_by]]))

  if (is.null(x = colors_use) && !plot_by == "sample_id") {
    if (length_plotby <= 8) {
      colors_use <- Dark2_Pal()
    } else {
      colors_use <- DiscretePalette_scCustomize(num_colors = length_plotby, palette = "polychrome")
    }
  } else {
    if (length(x = colors_use) < length_plotby && !plot_by == "sample_id") {
      cli_abort(message = c("Not enough colors provided.",
                            "i" = "The number of colors supplied: {.field {length(x = colors_use)}}, is less than the number of groups in {.val {plot_by}} column: {.field {length_plotby}}.")
      )
    } else {
      colors_use <- colors_use
    }
  }

  # Modify dataframe
  metrics_dataframe[,"Reads_Mapped_Confidently_to_Genome"] <- as.numeric(gsub("%", "", metrics_dataframe[,"Reads_Mapped_Confidently_to_Genome"]))


  if (plot_by == "sample_id") {
    metrics_dataframe$samples_plotting <- "Samples"

    plot <- ggplot(metrics_dataframe, aes(x = .data[["samples_plotting"]], y = .data[["Reads_Mapped_Confidently_to_Genome"]]), color = .data[["samples_plotting"]]) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle("Percent of Reads Confidently Mapped to Genome") +
      ylab('Percent of Reads ') +
      xlab("") +
      scale_y_continuous(labels = label_percent(accuracy = 0.01, scale = 1)) +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(metrics_dataframe, aes(x=.data[[plot_by]], y = .data[["Reads_Mapped_Confidently_to_Genome"]], fill = .data[[plot_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      ggtitle("Percent of Reads Confidently Mapped to Genome") +
      ylab('Percent of Reads') +
      xlab(plot_by) +
      scale_y_continuous(labels = label_percent(accuracy = 1, scale = 1)) +
      theme_ggprism_mod()
  }
  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (significance) {
    ggpubr_check <- PackageCheck("ggpubr", error = FALSE)
    if (!ggpubr_check[1]) {
      cli_abort(message = c(
        "Please install the {.val ggpubr} package to calculate/plot significance values.",
        "i" = "This can be accomplished with the following commands: ",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}ggpubr{symbol$dquote_right})`}",
        "----------------------------------------"
      ))
    }

    if (length(x = unique(x = stats_dataframe[[plot_by]])) < 2) {
      cli_abort(message = "Cannot calculate statistics when {.val {plot_by}} column contains less than 2 groups.")
    }
    groups <- unique(x = stats_dataframe[[plot_by]])

    comparisons <- combn(groups, 2)
    comparisons <- data.frame(comparisons, stringsAsFactors = FALSE)
    comparisons <- lapply(1:length(x = colnames(x = comparisons)), function(x){
      indiv_comp <- as.character(x = comparisons[[x]])
    })
    plot <- plot + ggpubr::stat_compare_means(comparisons = comparisons, ...)
  }

  return(plot)
}


#' QC Plots Sequencing metrics (Alignment)
#'
#' Plot the fraction of reads confidently mapped to intergenic regions
#'
#' @param metrics_dataframe data.frame  contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#' less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom scales label_percent
#' @importFrom SeuratObject PackageCheck
#' @importFrom utils combn
#'
#' @export
#'
#' @concept seq_qc_plotting_alignment
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_Intergeneic(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_Intergenic <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  significance = FALSE,
  ...
) {
  if (!plot_by %in% colnames(x = metrics_dataframe)) {
    cli_abort(message = "{.val {plot_by}} is not a column in the provided {.code metrics_dataframe}.")
  }

  # Change plot_by to character vector to make significance functions show all comparisons
  if (inherits(x = metrics_dataframe[[plot_by]], what = "factor")) {
    stats_dataframe <- metrics_dataframe
    stats_dataframe[[plot_by]] <- as.character(stats_dataframe[[plot_by]])
  } else {
    stats_dataframe <- metrics_dataframe
  }

  # Create color palette if null and check valid if provided
  length_plotby <- length(x = unique(x = metrics_dataframe[[plot_by]]))

  if (is.null(x = colors_use) && !plot_by == "sample_id") {
    if (length_plotby <= 8) {
      colors_use <- Dark2_Pal()
    } else {
      colors_use <- DiscretePalette_scCustomize(num_colors = length_plotby, palette = "polychrome")
    }
  } else {
    if (length(x = colors_use) < length_plotby && !plot_by == "sample_id") {
      cli_abort(message = c("Not enough colors provided.",
                            "i" = "The number of colors supplied: {.field {length(x = colors_use)}}, is less than the number of groups in {.val {plot_by}} column: {.field {length_plotby}}.")
      )
    } else {
      colors_use <- colors_use
    }
  }

  # Modify dataframe
  metrics_dataframe[,"Reads_Mapped_Confidently_to_Intergenic_Regions"] <- as.numeric(gsub("%", "", metrics_dataframe[,"Reads_Mapped_Confidently_to_Intergenic_Regions"]))


  if (plot_by == "sample_id") {
    metrics_dataframe$samples_plotting <- "Samples"

    plot <- ggplot(metrics_dataframe, aes(x = .data[["samples_plotting"]], y = .data[["Reads_Mapped_Confidently_to_Intergenic_Regions"]]), color = .data[["samples_plotting"]]) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle("Percent of Reads Confidently Mapped to Intergenic Regions") +
      ylab('Percent of Reads ') +
      xlab("") +
      scale_y_continuous(labels = label_percent(accuracy = 0.01, scale = 1)) +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(metrics_dataframe, aes(x=.data[[plot_by]], y = .data[["Reads_Mapped_Confidently_to_Intergenic_Regions"]], fill = .data[[plot_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      ggtitle("Percent of Reads Confidently Mapped to Intergenic Regions") +
      ylab('Percent of Reads') +
      xlab(plot_by) +
      scale_y_continuous(labels = label_percent(accuracy = 1, scale = 1)) +
      theme_ggprism_mod()
  }
  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (significance) {
    ggpubr_check <- PackageCheck("ggpubr", error = FALSE)
    if (!ggpubr_check[1]) {
      cli_abort(message = c(
        "Please install the {.val ggpubr} package to calculate/plot significance values.",
        "i" = "This can be accomplished with the following commands: ",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}ggpubr{symbol$dquote_right})`}",
        "----------------------------------------"
      ))
    }

    if (length(x = unique(x = stats_dataframe[[plot_by]])) < 2) {
      cli_abort(message = "Cannot calculate statistics when {.val {plot_by}} column contains less than 2 groups.")
    }
    groups <- unique(x = stats_dataframe[[plot_by]])

    comparisons <- combn(groups, 2)
    comparisons <- data.frame(comparisons, stringsAsFactors = FALSE)
    comparisons <- lapply(1:length(x = colnames(x = comparisons)), function(x){
      indiv_comp <- as.character(x = comparisons[[x]])
    })
    plot <- plot + ggpubr::stat_compare_means(comparisons = comparisons, ...)
  }

  return(plot)
}


#' QC Plots Sequencing metrics (Alignment)
#'
#' Plot the fraction of reads confidently mapped to intronic regions
#'
#' @param metrics_dataframe data.frame  contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#' less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom scales label_percent
#' @importFrom SeuratObject PackageCheck
#' @importFrom utils combn
#'
#' @export
#'
#' @concept seq_qc_plotting_alignment
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_Intronic(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_Intronic <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  significance = FALSE,
  ...
) {
  if (!plot_by %in% colnames(x = metrics_dataframe)) {
    cli_abort(message = "{.val {plot_by}} is not a column in the provided {.code metrics_dataframe}.")
  }

  # Change plot_by to character vector to make significance functions show all comparisons
  if (inherits(x = metrics_dataframe[[plot_by]], what = "factor")) {
    stats_dataframe <- metrics_dataframe
    stats_dataframe[[plot_by]] <- as.character(stats_dataframe[[plot_by]])
  } else {
    stats_dataframe <- metrics_dataframe
  }

  # Create color palette if null and check valid if provided
  length_plotby <- length(x = unique(x = metrics_dataframe[[plot_by]]))

  if (is.null(x = colors_use) && !plot_by == "sample_id") {
    if (length_plotby <= 8) {
      colors_use <- Dark2_Pal()
    } else {
      colors_use <- DiscretePalette_scCustomize(num_colors = length_plotby, palette = "polychrome")
    }
  } else {
    if (length(x = colors_use) < length_plotby && !plot_by == "sample_id") {
      cli_abort(message = c("Not enough colors provided.",
                            "i" = "The number of colors supplied: {.field {length(x = colors_use)}}, is less than the number of groups in {.val {plot_by}} column: {.field {length_plotby}}.")
      )
    } else {
      colors_use <- colors_use
    }
  }

  # Modify dataframe
  metrics_dataframe[,"Reads_Mapped_Confidently_to_Intronic_Regions"] <- as.numeric(gsub("%", "", metrics_dataframe[,"Reads_Mapped_Confidently_to_Intronic_Regions"]))


  if (plot_by == "sample_id") {
    metrics_dataframe$samples_plotting <- "Samples"

    plot <- ggplot(metrics_dataframe, aes(x = .data[["samples_plotting"]], y = .data[["Reads_Mapped_Confidently_to_Intronic_Regions"]]), color = .data[["samples_plotting"]]) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle("Percent of Reads Confidently Mapped to Intronic Regions") +
      ylab('Percent of Reads ') +
      xlab("") +
      scale_y_continuous(labels = label_percent(accuracy = 0.01, scale = 1)) +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(metrics_dataframe, aes(x=.data[[plot_by]], y = .data[["Reads_Mapped_Confidently_to_Intronic_Regions"]], fill = .data[[plot_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      ggtitle("Percent of Reads Confidently Mapped to Intronic Regions") +
      ylab('Percent of Reads') +
      xlab(plot_by) +
      scale_y_continuous(labels = label_percent(accuracy = 1, scale = 1)) +
      theme_ggprism_mod()
  }
  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (significance) {
    ggpubr_check <- PackageCheck("ggpubr", error = FALSE)
    if (!ggpubr_check[1]) {
      cli_abort(message = c(
        "Please install the {.val ggpubr} package to calculate/plot significance values.",
        "i" = "This can be accomplished with the following commands: ",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}ggpubr{symbol$dquote_right})`}",
        "----------------------------------------"
      ))
    }

    if (length(x = unique(x = stats_dataframe[[plot_by]])) < 2) {
      cli_abort(message = "Cannot calculate statistics when {.val {plot_by}} column contains less than 2 groups.")
    }
    groups <- unique(x = stats_dataframe[[plot_by]])

    comparisons <- combn(groups, 2)
    comparisons <- data.frame(comparisons, stringsAsFactors = FALSE)
    comparisons <- lapply(1:length(x = colnames(x = comparisons)), function(x){
      indiv_comp <- as.character(x = comparisons[[x]])
    })
    plot <- plot + ggpubr::stat_compare_means(comparisons = comparisons, ...)
  }

  return(plot)
}


#' QC Plots Sequencing metrics (Alignment)
#'
#' Plot the fraction of reads confidently mapped to Exonic regions
#'
#' @param metrics_dataframe data.frame  contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#' less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom scales label_percent
#' @importFrom SeuratObject PackageCheck
#' @importFrom utils combn
#'
#' @export
#'
#' @concept seq_qc_plotting_alignment
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_Exonic(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_Exonic <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  significance = FALSE,
  ...
) {
  if (!plot_by %in% colnames(x = metrics_dataframe)) {
    cli_abort(message = "{.val {plot_by}} is not a column in the provided {.code metrics_dataframe}.")
  }

  # Change plot_by to character vector to make significance functions show all comparisons
  if (inherits(x = metrics_dataframe[[plot_by]], what = "factor")) {
    stats_dataframe <- metrics_dataframe
    stats_dataframe[[plot_by]] <- as.character(stats_dataframe[[plot_by]])
  } else {
    stats_dataframe <- metrics_dataframe
  }

  # Create color palette if null and check valid if provided
  length_plotby <- length(x = unique(x = metrics_dataframe[[plot_by]]))

  if (is.null(x = colors_use) && !plot_by == "sample_id") {
    if (length_plotby <= 8) {
      colors_use <- Dark2_Pal()
    } else {
      colors_use <- DiscretePalette_scCustomize(num_colors = length_plotby, palette = "polychrome")
    }
  } else {
    if (length(x = colors_use) < length_plotby && !plot_by == "sample_id") {
      cli_abort(message = c("Not enough colors provided.",
                            "i" = "The number of colors supplied: {.field {length(x = colors_use)}}, is less than the number of groups in {.val {plot_by}} column: {.field {length_plotby}}.")
      )
    } else {
      colors_use <- colors_use
    }
  }

  # Modify dataframe
  metrics_dataframe[,"Reads_Mapped_Confidently_to_Exonic_Regions"] <- as.numeric(gsub("%", "", metrics_dataframe[,"Reads_Mapped_Confidently_to_Exonic_Regions"]))


  if (plot_by == "sample_id") {
    metrics_dataframe$samples_plotting <- "Samples"

    plot <- ggplot(metrics_dataframe, aes(x = .data[["samples_plotting"]], y = .data[["Reads_Mapped_Confidently_to_Exonic_Regions"]]), color = .data[["samples_plotting"]]) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle("Percent of Reads Confidently Mapped to Exonic Regions") +
      ylab('Percent of Reads ') +
      xlab("") +
      scale_y_continuous(labels = label_percent(accuracy = 0.01, scale = 1)) +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(metrics_dataframe, aes(x=.data[[plot_by]], y = .data[["Reads_Mapped_Confidently_to_Exonic_Regions"]], fill = .data[[plot_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      ggtitle("Percent of Reads Confidently Mapped to Exonic Regions") +
      ylab('Percent of Reads') +
      xlab(plot_by) +
      scale_y_continuous(labels = label_percent(accuracy = 1, scale = 1)) +
      theme_ggprism_mod()
  }
  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (significance) {
    ggpubr_check <- PackageCheck("ggpubr", error = FALSE)
    if (!ggpubr_check[1]) {
      cli_abort(message = c(
        "Please install the {.val ggpubr} package to calculate/plot significance values.",
        "i" = "This can be accomplished with the following commands: ",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}ggpubr{symbol$dquote_right})`}",
        "----------------------------------------"
      ))
    }

    if (length(x = unique(x = stats_dataframe[[plot_by]])) < 2) {
      cli_abort(message = "Cannot calculate statistics when {.val {plot_by}} column contains less than 2 groups.")
    }
    groups <- unique(x = stats_dataframe[[plot_by]])

    comparisons <- combn(groups, 2)
    comparisons <- data.frame(comparisons, stringsAsFactors = FALSE)
    comparisons <- lapply(1:length(x = colnames(x = comparisons)), function(x){
      indiv_comp <- as.character(x = comparisons[[x]])
    })
    plot <- plot + ggpubr::stat_compare_means(comparisons = comparisons, ...)
  }

  return(plot)
}


#' QC Plots Sequencing metrics (Alignment)
#'
#' Plot the fraction of reads mapped Antisense to Gene
#'
#' @param metrics_dataframe data.frame  contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#' less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import cli
#' @import ggplot2
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom scales label_percent
#' @importFrom SeuratObject PackageCheck
#' @importFrom utils combn
#'
#' @export
#'
#' @concept seq_qc_plotting_alignment
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_Antisense(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_Antisense <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  significance = FALSE,
  ...
) {
  if (!plot_by %in% colnames(x = metrics_dataframe)) {
    cli_abort(message = "{.val {plot_by}} is not a column in the provided {.code metrics_dataframe}.")
  }

  # Change plot_by to character vector to make significance functions show all comparisons
  if (inherits(x = metrics_dataframe[[plot_by]], what = "factor")) {
    stats_dataframe <- metrics_dataframe
    stats_dataframe[[plot_by]] <- as.character(stats_dataframe[[plot_by]])
  } else {
    stats_dataframe <- metrics_dataframe
  }

  # Create color palette if null and check valid if provided
  length_plotby <- length(x = unique(x = metrics_dataframe[[plot_by]]))

  if (is.null(x = colors_use) && !plot_by == "sample_id") {
    if (length_plotby <= 8) {
      colors_use <- Dark2_Pal()
    } else {
      colors_use <- DiscretePalette_scCustomize(num_colors = length_plotby, palette = "polychrome")
    }
  } else {
    if (length(x = colors_use) < length_plotby && !plot_by == "sample_id") {
      cli_abort(message = c("Not enough colors provided.",
                            "i" = "The number of colors supplied: {.field {length(x = colors_use)}}, is less than the number of groups in {.val {plot_by}} column: {.field {length_plotby}}.")
      )
    } else {
      colors_use <- colors_use
    }
  }

  # Modify dataframe
  metrics_dataframe[,"Reads_Mapped_Antisense_to_Gene"] <- as.numeric(gsub("%", "", metrics_dataframe[,"Reads_Mapped_Antisense_to_Gene"]))


  if (plot_by == "sample_id") {
    metrics_dataframe$samples_plotting <- "Samples"

    plot <- ggplot(metrics_dataframe, aes(x = .data[["samples_plotting"]], y = .data[["Reads_Mapped_Antisense_to_Gene"]]), color = .data[["samples_plotting"]]) +
      geom_boxplot(fill = "white", outlier.color = NA) +
      geom_quasirandom() +
      ggtitle("Percent of Reads Confidently Mapped to Antisense to Gene") +
      ylab('Percent of Reads ') +
      xlab("") +
      scale_y_continuous(labels = label_percent(accuracy = 0.01, scale = 1)) +
      theme_ggprism_mod()
  } else {
    plot <- ggplot(metrics_dataframe, aes(x=.data[[plot_by]], y = .data[["Reads_Mapped_Antisense_to_Gene"]], fill = .data[[plot_by]])) +
      geom_boxplot(fill = "white") +
      geom_dotplot(binaxis ='y', stackdir = 'center') +
      scale_fill_manual(values = colors_use) +
      ggtitle("Percent of Reads Confidently Mapped to Antisense to Gene") +
      ylab('Percent of Reads') +
      xlab(plot_by) +
      scale_y_continuous(labels = label_percent(accuracy = 1, scale = 1)) +
      theme_ggprism_mod()
  }
  if (x_lab_rotate) {
    plot <- plot + theme_ggprism_mod(axis_text_angle = 45)
  }

  if (significance) {
    ggpubr_check <- PackageCheck("ggpubr", error = FALSE)
    if (!ggpubr_check[1]) {
      cli_abort(message = c(
        "Please install the {.val ggpubr} package to calculate/plot significance values.",
        "i" = "This can be accomplished with the following commands: ",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}ggpubr{symbol$dquote_right})`}",
        "----------------------------------------"
      ))
    }

    if (length(x = unique(x = stats_dataframe[[plot_by]])) < 2) {
      cli_abort(message = "Cannot calculate statistics when {.val {plot_by}} column contains less than 2 groups.")
    }
    groups <- unique(x = stats_dataframe[[plot_by]])

    comparisons <- combn(groups, 2)
    comparisons <- data.frame(comparisons, stringsAsFactors = FALSE)
    comparisons <- lapply(1:length(x = colnames(x = comparisons)), function(x){
      indiv_comp <- as.character(x = comparisons[[x]])
    })
    plot <- plot + ggpubr::stat_compare_means(comparisons = comparisons, ...)
  }

  return(plot)
}


#' QC Plots Sequencing metrics (Layout)
#'
#' Plot a combined plot of the basic QC metrics from sequencing output.
#'
#' @param metrics_dataframe data.frame  contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#' less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param patchwork_title Title to use for the patchworked plot output.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom SeuratObject PackageCheck
#' @importFrom stringr str_wrap
#'
#' @export
#'
#' @concept seq_qc_plotting_layout
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_Basic_Combined(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_Basic_Combined <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  patchwork_title = "Sequencing QC Plots: Basic Cell Metrics",
  significance = FALSE,
  ...
) {
  # Create rotated axis value
  if (x_lab_rotate) {
    axis_angle <- 45
  } else {
    axis_angle <- 0
  }

  # Create Plots & modify for plotting together
  p1 <- Seq_QC_Plot_Number_Cells(metrics_dataframe = metrics_dataframe, plot_by = plot_by, colors_use = colors_use, significance = significance, ...)
  p1 <- p1 +
    labs(title = str_wrap(p1$labels$title, 18)) +
    theme_ggprism_mod(base_size = 10, axis_text_angle = axis_angle)

  p2 <- Seq_QC_Plot_Reads_per_Cell(metrics_dataframe = metrics_dataframe, plot_by = plot_by, colors_use = colors_use, significance = significance, ...)
  p2 <- p2 + labs(title = str_wrap(p2$labels$title, 18)) +
    theme_ggprism_mod(base_size = 10, axis_text_angle = axis_angle)

  p3 <- Seq_QC_Plot_Genes(metrics_dataframe = metrics_dataframe, plot_by = plot_by, colors_use = colors_use, significance = significance, ...)
  p3 <- p3 + labs(title = str_wrap(p3$labels$title, 18)) +
    theme_ggprism_mod(base_size = 10, axis_text_angle = axis_angle)

  p4 <- Seq_QC_Plot_UMIs(metrics_dataframe = metrics_dataframe, plot_by = plot_by, colors_use = colors_use, significance = significance, ...)
  p4 <- p4 + labs(title = str_wrap(p4$labels$title, 18)) +
    theme_ggprism_mod(base_size = 10, axis_text_angle = axis_angle)

  p5 <- Seq_QC_Plot_Total_Genes(metrics_dataframe = metrics_dataframe, plot_by = plot_by, colors_use = colors_use, significance = significance, ...)
  p5 <- p5 + labs(title = str_wrap(p5$labels$title, 18)) +
    theme_ggprism_mod(base_size = 10, axis_text_angle = axis_angle)

  p6 <- Seq_QC_Plot_Saturation(metrics_dataframe = metrics_dataframe, plot_by = plot_by, colors_use = colors_use, significance = significance, ...)
  p6 <- p6 + labs(title = str_wrap(p6$labels$title, 18)) +
    theme_ggprism_mod(base_size = 10, axis_text_angle = axis_angle)

  p7 <- Seq_QC_Plot_Reads_in_Cells(metrics_dataframe = metrics_dataframe, plot_by = plot_by, colors_use = colors_use, significance = significance, ...)
  p7 <- p7 + labs(title = str_wrap(p7$labels$title, 18)) +
    theme_ggprism_mod(base_size = 10, axis_text_angle = axis_angle)

  p8 <- Seq_QC_Plot_Transcriptome(metrics_dataframe = metrics_dataframe, plot_by = plot_by, colors_use = colors_use, significance = significance, ...)
  p8 <- p8 + labs(title = str_wrap(p8$labels$title, 18)) +
    theme_ggprism_mod(base_size = 10, axis_text_angle = axis_angle)

  # Assemble plots and unifying legends
  plot <- (p1 | p3 | p5 | p7) /
    (p2 | p4 | p6 | p8)
  plot <- plot + plot_layout(guides = 'collect') + plot_annotation(title = patchwork_title, theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5))))

  # Print plots
  suppressMessages(print(plot))
}


#' QC Plots Sequencing metrics (Alignment) (Layout)
#'
#' Plot a combined plot of the Alignment QC metrics from sequencing output.
#'
#' @param metrics_dataframe data.frame  contain Cell Ranger QC Metrics (see \code{\link{Read_Metrics_10X}}).
#' @param plot_by Grouping factor for the plot.  Default is to plot as single group with single point per sample.
#' @param colors_use colors to use for plot if plotting by group.  Defaults to RColorBrewer Dark2 palette if
#' less than 8 groups and `DiscretePalette_scCustomize(palette = "polychrome")` if more than 8.
#' @param x_lab_rotate logical.  Whether to rotate the axes labels on the x-axis.  Default is FALSE.
#' @param patchwork_title Title to use for the patchworked plot output.
#' @param significance logical.  Whether to calculate and plot p-value comparisons when plotting by
#' grouping factor.  Default is FALSE.
#' @param ... Other variables to pass to `ggpubr::stat_compare_means` when doing significance testing.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom SeuratObject PackageCheck
#' @importFrom stringr str_wrap
#'
#' @export
#'
#' @concept seq_qc_plotting_layout
#'
#' @examples
#' \dontrun{
#' Seq_QC_Plot_Alignment_Combined(metrics_dataframe = metrics)
#' }
#'

Seq_QC_Plot_Alignment_Combined <- function(
  metrics_dataframe,
  plot_by = "sample_id",
  colors_use = NULL,
  x_lab_rotate = FALSE,
  patchwork_title = "Sequencing QC Plots: Read Alignment Metrics",
  significance = FALSE,
  ...
) {
  # Create rotated axis value
  if (x_lab_rotate) {
    axis_angle <- 45
  } else {
    axis_angle <- 0
  }

  # Create Plots & modify for plotting together
  p1 <- Seq_QC_Plot_Genome(metrics_dataframe = metrics_dataframe, plot_by = plot_by, colors_use = colors_use, significance = significance, ...)
  p1 <- p1 +
    labs(title = str_wrap(p1$labels$title, 18)) +
    theme_ggprism_mod(base_size = 10, axis_text_angle = axis_angle)

  p2 <- Seq_QC_Plot_Intergenic(metrics_dataframe = metrics_dataframe, plot_by = plot_by, colors_use = colors_use, significance = significance, ...)
  p2 <- p2 + labs(title = str_wrap(p2$labels$title, 18)) +
    theme_ggprism_mod(base_size = 10, axis_text_angle = axis_angle)

  p3 <- Seq_QC_Plot_Transcriptome(metrics_dataframe = metrics_dataframe, plot_by = plot_by, colors_use = colors_use, significance = significance, ...)
  p3 <- p3 + labs(title = str_wrap(p3$labels$title, 18)) +
    theme_ggprism_mod(base_size = 10, axis_text_angle = axis_angle)

  p4 <- Seq_QC_Plot_Exonic(metrics_dataframe = metrics_dataframe, plot_by = plot_by, colors_use = colors_use, significance = significance, ...)
  p4 <- p4 + labs(title = str_wrap(p4$labels$title, 18)) +
    theme_ggprism_mod(base_size = 10, axis_text_angle = axis_angle)

  p5 <- Seq_QC_Plot_Intronic(metrics_dataframe = metrics_dataframe, plot_by = plot_by, colors_use = colors_use, significance = significance, ...)
  p5 <- p5 + labs(title = str_wrap(p5$labels$title, 18)) +
    theme_ggprism_mod(base_size = 10, axis_text_angle = axis_angle)

  p6 <- Seq_QC_Plot_Antisense(metrics_dataframe = metrics_dataframe, plot_by = plot_by, colors_use = colors_use, significance = significance, ...)
  p6 <- p6 + labs(title = str_wrap(p6$labels$title, 18)) +
    theme_ggprism_mod(base_size = 10, axis_text_angle = axis_angle)

  # Assemble plots and unifying legends
  plot <- (p1 | p3 | p5) /
    (p2 | p4 | p6)
  plot <- plot + plot_layout(guides = 'collect') + plot_annotation(title = patchwork_title, theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = rel(1.5))))

  # Print plots
  suppressMessages(print(plot))
}
