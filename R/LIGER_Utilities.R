#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### DATA ACCESS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname Fetch_Meta
#' @importFrom methods slot
#' @export
#' @concept liger_object_util
#' @method Fetch_Meta liger

Fetch_Meta.liger <- function(
    object,
    ...
) {
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    object_meta <- rliger::cellMeta(x = object, as.data.frame = TRUE)
  } else {
    object_meta <- object_meta <- slot(object = object, name = "cell.data")
  }

  # return meta
  return(object_meta)
}


#' Extract Features from LIGER Object
#'
#' Extract all unique features from LIGER object
#'
#' @param liger_object LIGER object name.
#' @param by_dataset logical, whether to return list with vector of features for each dataset in
#' LIGER object or to return single vector of unique features across all datasets in object
#' (default is FALSE; return vector of unique features)
#'
#' @return vector or list depending on `by_dataset` parameter
#'
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # return single vector of all unique features
#' all_features <- LIGER_Features(liger_object = object, by_dataset = FALSE)
#'
#' # return list of vectors containing features from each individual dataset in object
#' dataset_features <- LIGER_Features(liger_object = object, by_dataset = TRUE)
#' }
#'

LIGER_Features <- function(
    liger_object,
    by_dataset = FALSE
) {
  # check liger
  Is_LIGER(liger_object = liger_object)

  # liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    # Extract features
    features_by_dataset <- lapply(1:length(x = liger_object@datasets), function(x) {
      rownames(x = liger_object@datasets[[x]]@featureMeta)
    })
  } else {
    # Extract features
    features_by_dataset <- lapply(1:length(x = liger_object@raw.data), function(x) {
      rownames(x = liger_object@raw.data[[x]])
    })
  }

  # Return features
  if (isFALSE(x = by_dataset)) {
    features <- unique(x = unlist(x = features_by_dataset))
    return(features)
  } else {
    return(features_by_dataset)
  }
}


#' Extract Cells from LIGER Object
#'
#' Extract all cell barcodes from LIGER object
#'
#' @param liger_object LIGER object name.
#' @param by_dataset logical, whether to return list with vector of cell barcodes for each
#' dataset in LIGER object or to return single vector of cell barcodes across all
#' datasets in object (default is FALSE; return vector of cells)
#'
#' @return vector or list depending on `by_dataset` parameter
#'
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # return single vector of all cells
#' all_features <- LIGER_Cells(liger_object = object, by_dataset = FALSE)
#'
#' # return list of vectors containing cells from each individual dataset in object
#' dataset_features <- LIGER_Cells(liger_object = object, by_dataset = TRUE)
#' }
#'

LIGER_Cells <- function(
    liger_object,
    by_dataset = FALSE
) {
  # check liger
  Is_LIGER(liger_object = liger_object)

  # liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    # Extract features
    cells_by_dataset <- lapply(1:length(x = liger_object@datasets), function(x) {
      colnames(x = liger_object@datasets[[x]])
    })
    names(cells_by_dataset) <- names(liger_object@datasets)
  } else {
    # Extract features
    cells_by_dataset <- lapply(1:length(x = liger_object@raw.data), function(x) {
      colnames(x = liger_object@raw.data[[x]])
    })
  }

  # Return features
  if (isFALSE(x = by_dataset)) {
    cells <- unlist(x = cells_by_dataset)
    return(cells)
  } else {
    return(cells_by_dataset)
  }
}


#' Extract Cells by identity
#'
#' Extract all cell barcodes by identity from LIGER object
#'
#' @param liger_object LIGER object name.
#' @param group.by name of meta data column to use, default is current default clustering.
#' @param by_dataset logical, whether to return list with entries for cell barcodes for each
#' identity in `group.by`
#' or to return list of lists (1 entry per dataset and each ident within the dataset)
#' (default is FALSE; return list)
#'
#' @return list or list of lists depending on `by_dataset` parameter
#'
#' @import cli
#' @importFrom dplyr filter select all_of
#' @importFrom magrittr "%>%"
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # return single vector of all cells
#' cells_by_idents <- LIGER_Cells_by_Identities(liger_object = object, by_dataset = FALSE)
#'
#' # return list of vectors containing cells from each individual dataset in object
#' cells_by_idents_by_dataset <- LIGER_Cells_by_Identities(liger_object = object, by_dataset = TRUE)
#' }
#'

LIGER_Cells_by_Identities <- function(
    liger_object,
    group.by = NULL,
    by_dataset = FALSE
) {
  # Check new liger object
  if (!"cellMeta" %in% slotNames(liger_object)) {
    cli_abort(message = "This function is only for objects created with rliger >= v2.0.0")
  }

  # check group.by is valid
  if (!is.null(x = group.by)) {
    Meta_Present(object = liger_object, meta_col_names = group.by, print_msg = FALSE)
  }

  # set group.by if not set
  group.by <- group.by %||% LIGER_Default_Cluster(liger_object = liger_object)

  # Check cluster df
  cell_df <- Fetch_Meta(object = liger_object) %>%
    select(all_of(c(group.by, "dataset")))

  if (inherits(x = cell_df[[group.by]], what = "factor")) {
    ident_levels <- levels(x = cell_df[[group.by]])
  } else {
    ident_levels <- unique(x = cell_df[[group.by]])
  }

  # Get cells for object overall
  if (isFALSE(x = by_dataset)) {
    cells_list <- lapply(ident_levels, function(x) {
      cells <- cell_df %>%
        filter(.data[[group.by]] == x) %>%
        rownames()
    })

    names(cells_list) <- ident_levels
  } else {
    # Get cells by cluster by dataset
    dataset_names <- names(x = rliger::datasets(x = liger_object))
    cells_list <- lapply(1:length(x = dataset_names), function(x) {
      sample_cells_df <- cell_df %>%
        filter(.data[["dataset"]] == dataset_names[x])

      sample_cells <- lapply(ident_levels, function(y) {
        sample_cells_df %>%
          filter(.data[[group.by]] == y) %>%
          rownames()
      })
      names(sample_cells) <- ident_levels

      return(sample_cells)
    })
    names(cells_list) <- dataset_names
  }

  return(cells_list)
}


#' Subset LIGER object
#'
#' Subset LIGER object by cluster or other meta data variable.
#'
#' @param liger_object LIGER object name.
#' @param cluster Name(s) of cluster to subset from object.
#' @param cluster_col name of `@cellMeta` column containing cluster names, default is "leiden_cluster".
#' @param ident variable within `ident_col` to use in sub-setting object.
#' @param ident_col column in `@cellMeta` that contains values provided to `ident`.
#' @param invert logical, whether to subset the inverse of the clusters or idents provided, default is FALSE.
#'
#' @return liger object
#'
#' @import cli
#' @importFrom dplyr pull filter
#' @importFrom magrittr "%>%"
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # subset clusters 3 and 5
#' sub_liger <- subset_liger(liger_object = liger_object, cluster = c(3, 5))
#'
#' # subset control samples from column "Treatment"
#' sub_liger <- subset_liger(liger_object = liger_object, ident = "control",
#' ident_col = "Treatment")
#'
#' # subset control samples from column "Treatment" in clusters 3 and 5
#' sub_liger <- subset_liger(liger_object = liger_object, ident = "control",
#' ident_col = "Treatment", cluster = c(3, 5))
#'
#' # Remove cluster 9
#' sub_liger <- subset_liger(liger_object = liger_object, cluster = 9, invert = TRUE)
#' }
#'

Subset_LIGER <- function(
    liger_object,
    cluster = NULL,
    cluster_col = "leiden_cluster",
    ident = NULL,
    ident_col = NULL,
    invert = FALSE
) {
  # Check new liger object
  if (!"cellMeta" %in% slotNames(liger_object)) {
    cli_abort(message = "This function is only for objects created with rliger >= v2.0.0")
  }

  # Check value provided
  if (is.null(x = cluster) && is.null(x = ident)) {
    cli_abort(message = "No values provided to subset object")
  }

  # Check meta present
  if (!is.null(x = ident_col)) {
    ident_col <- Meta_Present(object = liger_object, meta_col_names = ident_col, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  # Check meta present
  if (!is.null(x = cluster_col)) {
    cluster_col <- Meta_Present(object = liger_object, meta_col_names = cluster_col, print_msg = FALSE, omit_warn = FALSE)[[1]]
  }

  # pull meta data
  meta <- Fetch_Meta(object = liger_object)

  # check subset value ok
  if (!is.null(x = ident)) {
    ident_values <- meta %>%
      pull(.data[[ident_col]]) %>%
      unique()

    if (!all(ident %in% ident_values)) {
      cli_abort(message = "One or more of provided ident values ({.field {ident}}) were not found in provided ident_col ({.field {ident_col}})")
    }
  }

  # check sub set value ok
  if (!is.null(x = cluster)) {
    cluster_values <- meta %>%
      pull(.data[[cluster_col]]) %>%
      unique()

    if (!all(cluster %in% cluster_values)) {
      cli_abort(message = "One or more of provided cluster values ({.field {cluster}}) were not found in provided cluster_col ({.field {cluster_col}})")
    }
  }

  # filter just by cluster
  if (!is.null(x = cluster) && is.null(x = ident)) {
    cells_filter <- meta %>%
      filter(.data[[cluster_col]] %in% cluster) %>%
      rownames()
  }

  # filter just by ident
  if (!is.null(x = ident) && is.null(x = cluster)) {
    cells_filter <- meta %>%
      filter(.data[[ident_col]] %in% ident) %>%
      rownames()
  }

  # Filter by ident and cluster
  if (!is.null(x = ident) && !is.null(x = cluster)) {
    cells_filter <- meta %>%
      filter(.data[[ident_col]] %in% ident & .data[[cluster_col]] %in% cluster) %>%
      rownames()
  }

  # invert filtering
  if (isTRUE(x = invert)) {
    # get vector of call cells
    all_cells <- LIGER_Cells(liger_object = liger_object)

    # setdiff to get inverse
    cells_filter <- setdiff(x = all_cells, y = cells_filter)
  }

  sub_obj <- rliger::subsetLiger(object = liger_object, cellIdx = cells_filter)

  return(sub_obj)
}


#' Extract top loading genes for LIGER factor
#'
#' Extract vector to the top loading genes for specified LIGER iNMF factor
#'
#' @param liger_object LIGER object name.
#' @param liger_factor LIGER factor number to pull genes from.
#' @param num_genes number of top loading genes to return as vector.
#'
#' @return A LIGER Object
#'
#' @import cli
#' @importFrom utils packageVersion
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' top_genes_factor10 <- Top_Genes_Factor(liger_object = object, num_genes = 10)
#' }
#'

Top_Genes_Factor <- function(
    liger_object,
    liger_factor,
    num_genes = 10
) {
  # LIGER object check
  Is_LIGER(liger_object = liger_object)

  # check number of factors present
  if (!liger_factor %in% 1:dim(x = liger_object@W)[[1]]) {
    cli_abort(message = c("{.code liger_factor} provided: {.field {liger_factor}} not found",
                          "i" = "{.code liger_object} only contains {.field {dim(x = liger_object@W)[[1]]}} factors.")
    )
  }

  # temp liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    W <- liger_object@W
    rownames(x = W) <- rownames(x = liger_object@datasets[[1]]@scaleData)
    top_genes <- rownames(x = W)[order(W[, liger_factor], decreasing = TRUE)[1:num_genes]]
    return(top_genes)
  } else {
    # Extract genes
    W <- t(x = liger_object@W)
    rownames(x = W) <- colnames(x = liger_object@scale.data[[1]])
    top_genes <- rownames(x = W)[order(W[, liger_factor], decreasing = TRUE)[1:num_genes]]
    return(top_genes)
  }
}





#' Extract dimensionality reduction coordinates from Liger object
#'
#' Extract data.frame containing dimensionality reduction coordinates from new format of
#' Liger objects
#'
#' @param liger_object LIGER object name.
#' @param reduction name of dimensionality reduction stored in cellMeta slot.  Default is
#' NULL, which will use liger object's default reduction.
#' @param check_only logical, return `TRUE` if valid reduction is present.
#'
#' @return dimensionality reduction coordinates in 2 column format
#'
#' @import cli
#' @importFrom methods slotNames
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # return dimensionality reduction coordinates
#' umap_coords <- LIGER_DimReduc(liger_object = object)
#'
#' # return logical to see if reduction is present
#' reduc_present <- LIGER_DimReduc(liger_object = object, reduction = "umap",
#' check_only = TRUE)
#' }
#'

LIGER_DimReduc <- function(
    liger_object,
    reduction = NULL,
    check_only = FALSE
) {
  # Check new liger object
  if (!"cellMeta" %in% slotNames(liger_object)) {
    cli_abort(message = "This function is only for objects created with rliger >= v2.0.0")
  }

  # reduction to use
  reduction_use <- reduction %||% Default_DimReduc_LIGER(liger_object = liger_object)

  # check reduction in cellMeta
  if (reduction_use %in% names(x = rliger::dimReds(x = liger_object))) {
      if (isTRUE(x = check_only)) {
        return(TRUE)
      }
      # get coords
      reduc_coords <- rliger::dimReds(x = liger_object)[[reduction_use]]
  } else {
    cli_abort("The reduction {.field {reduction_use}} is not present in dimReds slot.")
  }

  # return coords
  return(reduc_coords)
}


#' Find Factor Correlations
#'
#' Calculate correlations between gene loadings for all factors in liger object.
#'
#' @param liger_object LIGER object name.
#'
#' @return correlation matrix
#'
#' @import cli
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' factor_correlations <- Find_Factor_Cor(liger_object = object)
#' }
#'

Find_Factor_Cor <- function(
    liger_object
) {
  Is_LIGER(liger_object = liger_object)

  # Check new liger object
  if (!"cellMeta" %in% slotNames(liger_object)) {
    cli_abort(message = "This function is only for objects created with rliger >= v2.0.0")
  }

  # Get loadings
  factor_loadings <- data.frame(rliger::getMatrix(x = liger_object, slot = "W"))

  # Rename is zero padding
  colnames(x = factor_loadings) <- paste0("Factor_", seq_zeros(seq_length = ncol(x = factor_loadings), num_zeros = 1))

  # Correlation
  cor_mat <- cor(x = factor_loadings)

  return(cor_mat)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### QC UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param species Species of origin for given Seurat Object.  If mouse, human, marmoset, zebrafish, rat,
#' drosophila, rhesus macaque, or chicken (name or abbreviation) are provided the function will automatically
#' generate patterns and features.
#' @param add_mito_ribo logical, whether to add percentage of counts belonging to mitochondrial/ribosomal
#' genes to object (Default is TRUE).
#' @param add_complexity logical, whether to add Cell Complexity to object (Default is TRUE).
#' @param add_top_pct logical, whether to add Top Gene Percentages to object (Default is TRUE).
#' @param add_MSigDB logical, whether to add percentages of counts belonging to genes from of mSigDB hallmark
#' gene lists: "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_APOPTOSIS", and "HALLMARK_DNA_REPAIR" to
#' object (Default is TRUE).
#' @param add_IEG logical, whether to add percentage of counts belonging to IEG genes to object (Default is TRUE).
#' @param add_hemo logical, whether to add percentage of counts belonging to homoglobin genes to object (Default is TRUE).
#' @param mito_name name to use for the new meta.data column containing percent mitochondrial counts.
#' Default is "percent_mito".
#' @param ribo_name name to use for the new meta.data column containing percent ribosomal counts.
#' Default is "percent_ribo".
#' @param mito_ribo_name name to use for the new meta.data column containing percent
#' mitochondrial+ribosomal counts.  Default is "percent_mito_ribo".
#' @param complexity_name name to use for new meta data column for `Add_Cell_Complexity`.
#' Default is "log10GenesPerUMI".
#' @param top_pct_name name to use for new meta data column for `Add_Top_Gene_Pct`.
#' Default is "percent_topXX", where XX is equal to the value provided to `num_top_genes`.
#' @param oxphos_name name to use for new meta data column for percentage of MSigDB oxidative phosphorylation
#' counts.  Default is "percent_oxphos".
#' @param apop_name name to use for new meta data column for percentage of MSigDB apoptosis counts.
#' Default is "percent_apop".
#' @param dna_repair_name name to use for new meta data column for percentage of MSigDB DNA repair
#' counts.  Default is "percent_dna_repair"..
#' @param ieg_name name to use for new meta data column for percentage of IEG counts.  Default is "percent_ieg".
#' @param hemo_name name to use for the new meta.data column containing percent hemoglobin counts.
#' Default is "percent_mito".
#' @param mito_pattern A regex pattern to match features against for mitochondrial genes (will set automatically if
#' species is mouse or human; marmoset features list saved separately).
#' @param ribo_pattern A regex pattern to match features against for ribosomal genes
#' (will set automatically if species is in default list).
#' @param hemo_pattern A regex pattern to match features against for hemoglobin genes
#' (will set automatically if species is in default list).
#' @param mito_features A list of mitochondrial gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param ribo_features A list of ribosomal gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param hemo_features A list of hemoglobin gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param ensembl_ids logical, whether feature names in the object are gene names or
#' ensembl IDs (default is FALSE; set TRUE if feature names are ensembl IDs).
#' @param num_top_genes An integer vector specifying the size(s) of the top set of high-abundance genes.
#' Used to compute the percentage of library size occupied by the most highly expressed genes in each cell.
#' @param overwrite Logical.  Whether to overwrite existing an meta.data column.  Default is FALSE meaning that
#' function will abort if column with name provided to `meta_col_name` is present in meta.data slot.
#'
#' @import cli
#'
#' @return A liger Object
#'
#' @method Add_Cell_QC_Metrics liger
#'
#' @export
#' @rdname Add_Cell_QC_Metrics
#'
#'
#' @concept qc_util
#'
#' @examples
#' \dontrun{
#' obj <- Add_Cell_QC_Metrics(object = obj, species = "Human")
#'}
#'

Add_Cell_QC_Metrics.liger <- function(
    object,
    add_mito_ribo = TRUE,
    add_complexity = TRUE,
    add_top_pct = TRUE,
    add_MSigDB = TRUE,
    add_IEG = TRUE,
    add_hemo = TRUE,
    add_cell_cycle = TRUE,
    species,
    mito_name = "percent_mito",
    ribo_name = "percent_ribo",
    mito_ribo_name = "percent_mito_ribo",
    complexity_name = "log10GenesPerUMI",
    top_pct_name = NULL,
    oxphos_name = "percent_oxphos",
    apop_name = "percent_apop",
    dna_repair_name = "percent_dna_repair",
    ieg_name = "percent_ieg",
    hemo_name = "percent_hemo",
    mito_pattern = NULL,
    ribo_pattern = NULL,
    hemo_pattern = NULL,
    mito_features = NULL,
    ribo_features = NULL,
    hemo_features = NULL,
    ensembl_ids = FALSE,
    num_top_genes = 50,
    assay = NULL,
    overwrite = FALSE,
    ...
) {
  # Accepted species names
  accepted_names <- data.frame(
    Mouse_Options = c("Mouse", "mouse", "Ms", "ms", "Mm", "mm"),
    Human_Options = c("Human", "human", "Hu", "hu", "Hs", "hs"),
    Marmoset_Options = c("Marmoset", "marmoset", "CJ", "Cj", "cj", NA),
    Zebrafish_Options = c("Zebrafish", "zebrafish", "DR", "Dr", "dr", NA),
    Rat_Options = c("Rat", "rat", "RN", "Rn", "rn", NA),
    Drosophila_Options = c("Drosophila", "drosophila", "DM", "Dm", "dm", NA),
    Macaque_Options = c("Macaque", "macaque", "Rhesus", "macaca", "mmulatta", NA)
  )

  # Species Spelling Options
  mouse_options <- accepted_names$Mouse_Options
  human_options <- accepted_names$Human_Options
  marmoset_options <- accepted_names$Marmoset_Options
  zebrafish_options <- accepted_names$Zebrafish_Options
  rat_options <- accepted_names$Rat_Options
  drosophila_options <- accepted_names$Drosophila_Options
  macaque_options <- accepted_names$Macaque_Options

  # Add mito/ribo
  if (isTRUE(x = add_mito_ribo)) {
    cli_inform(message = c("*" = "Adding {.field Mito/Ribo Percentages} to meta.data."))
    liger_object <- Add_Mito_Ribo(object = liger_object, species = species, mito_name = mito_name, ribo_name = ribo_name, mito_ribo_name = mito_ribo_name, mito_pattern = mito_pattern, ribo_pattern = ribo_pattern, mito_features = mito_features, ribo_features = ribo_features, ensembl_ids = ensembl_ids, overwrite = overwrite)
  }

  # Add complexity
  if (isTRUE(x = add_complexity)) {
    cli_inform(message = c("*" = "Adding {.field Cell Complexity #1 (log10GenesPerUMI)} to meta.data."))
    liger_object <- Add_Cell_Complexity(object = liger_object, meta_col_name = complexity_name, overwrite = overwrite)
  }

  # Add top gene expression percent
  if (isTRUE(x = add_top_pct)) {
    cli_inform(message = c("*" = "Adding {.field Cell Complexity #2 (Top {num_top_genes} Percentages)} to meta.data."))
    liger_object <- Add_Top_Gene_Pct(object = liger_object, num_top_genes = num_top_genes, meta_col_name = top_pct_name, overwrite = overwrite)
  }

  # Add MSigDB
  if (isTRUE(x = add_MSigDB)) {
    if (species %in% marmoset_options) {
      cli_warn(message = c("{.val Marmoset} is not currently a part of MSigDB gene list database.",
                           "i" = "No columns will be added to object meta.data"))
    } else {
      cli_inform(message = c("*" = "Adding {.field MSigDB Oxidative Phosphorylation, Apoptosis, and DNA Repair Percentages} to meta.data."))
      liger_object <- Add_MSigDB_LIGER(liger_object = liger_object, species = species, oxphos_name = oxphos_name, apop_name = apop_name, dna_repair_name = dna_repair_name, overwrite = overwrite)
    }
  }

  # Add IEG
  if (isTRUE(x = add_IEG)) {
    if (species %in% c(marmoset_options, rat_options, zebrafish_options, macaque_options, drosophila_options)) {
      cli_warn(message = c("{.val Rat, Marmoset, Macaque, Zebrafish, and Drosophila} are not currently supported.",
                           "i" = "No column will be added to object meta.data"))
    } else {
      cli_inform(message = c("*" = "Adding {.field IEG Percentages} to meta.data."))
      liger_object <- Add_IEG_LIGER(liger_object = liger_object, species = species, ieg_name = ieg_name, overwrite = overwrite)
    }
  }

  # Add hemo
  if (isTRUE(x = add_hemo)) {
    cli_inform(message = c("*" = "Adding {.field Hemoglobin Percentages} to meta.data."))
    liger_object <- Add_Hemo(object = liger_object, species = species, hemo_name = hemo_name, hemo_pattern = hemo_pattern, hemo_features = hemo_features, overwrite = overwrite)
  }

  # return object
  return(liger_object)
}



#' Add Mito and Ribo percentages
#'
#' @param species Species of origin for given Object.  If mouse, human, marmoset, zebrafish, rat,
#' drosophila, rhesus macaque, or chicken (name or abbreviation) are provided the function will automatically
#' generate mito_pattern and ribo_pattern values.
#' @param mito_name name to use for the new meta.data column containing percent mitochondrial counts.
#' Default is "percent_mito".
#' @param ribo_name name to use for the new meta.data column containing percent ribosomal counts.
#' Default is "percent_ribo".
#' @param mito_ribo_name name to use for the new meta.data column containing percent mitochondrial+ribosomal
#' counts.  Default is "percent_mito_ribo".
#' @param mito_pattern A regex pattern to match features against for mitochondrial genes (will set automatically
#' if species is mouse or human; marmoset features list saved separately).
#' @param ribo_pattern A regex pattern to match features against for ribosomal genes (will set automatically
#' if species is mouse, human, or marmoset).
#' @param mito_features A list of mitochondrial gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param ribo_features A list of ribosomal gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param ensembl_ids logical, whether feature names in the object are gene names or
#' ensembl IDs (default is FALSE; set TRUE if feature names are ensembl IDs).
#' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
#' function will abort if columns with any one of the names provided to `mito_name` `ribo_name` or `mito_ribo_name`
#' is present in meta.data slot.
#' @param list_species_names returns list of all accepted values to use for default species names which
#' contain internal regex/feature lists (human, mouse, marmoset, zebrafish, rat, drosophila, and
#' rhesus macaque).  Default is FALSE.
#'
#' @import cli
#' @importFrom dplyr mutate select intersect
#' @importFrom magrittr "%>%"
#' @importFrom rlang ":="
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom utils packageVersion
#'
#' @method Add_Mito_Ribo liger
#'
#' @export
#' @rdname Add_Mito_Ribo
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # Liger
#' liger_object <- Add_Mito_Ribo(object = liger_object, species = "human")
#' }
#'

Add_Mito_Ribo.liger <- function(
  object,
  species,
  mito_name = "percent_mito",
  ribo_name = "percent_ribo",
  mito_ribo_name = "percent_mito_ribo",
  mito_pattern = NULL,
  ribo_pattern = NULL,
  mito_features = NULL,
  ribo_features = NULL,
  ensembl_ids = FALSE,
  overwrite = FALSE,
  list_species_names = FALSE,
  ...
) {
  # Accepted species names
  accepted_names <- data.frame(
    Mouse_Options = c("Mouse", "mouse", "Ms", "ms", "Mm", "mm"),
    Human_Options = c("Human", "human", "Hu", "hu", "Hs", "hs"),
    Marmoset_Options = c("Marmoset", "marmoset", "CJ", "Cj", "cj", NA),
    Zebrafish_Options = c("Zebrafish", "zebrafish", "DR", "Dr", "dr", NA),
    Rat_Options = c("Rat", "rat", "RN", "Rn", "rn", NA),
    Drosophila_Options = c("Drosophila", "drosophila", "DM", "Dm", "dm", NA),
    Macaque_Options = c("Macaque", "macaque", "Rhesus", "macaca", "mmulatta", NA),
    Chicken_Options = c("Chicken", "chicken", "Gallus", "gallus", "Gg", "gg")
  )

  # Return list of accepted default species name options
  if (isTRUE(x = list_species_names)) {
    return(accepted_names)
    stop_quietly()
  }

  # LIGER object check
  Is_LIGER(liger_object = object)

  # Check name collision
  if (any(duplicated(x = c(mito_name, ribo_name, mito_ribo_name)))) {
    cli_abort(message = "One or more of values provided to {.code mito_name}, {.code ribo_name}, {.code mito_ribo_name} are identical.")
  }

  # Overwrite check
  meta_names <- colnames(x = Fetch_Meta(object = object))

  if (mito_name %in% meta_names || ribo_name %in% meta_names || mito_ribo_name %in% meta_names) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Columns with {.val {mito_name}} and/or {.val {ribo_name}} already present in meta data.",
                            "i" = "*To run function and overwrite columns set parameter {.code overwrite = TRUE} or change respective {.code mito_name}, {.code ribo_name}, and/or {.code mito_ribo_name}.*")
      )
    }
    cli_inform(message = c("Columns with {.val {mito_name}} and/or {.val {ribo_name}} already present in meta data.",
                           "i" = "Overwriting those columns as {.code overwrite = TRUE}.")
    )
  }

  # Checks species
  if (is.null(x = species)) {
    cli_abort(message = c("No species name or abbreivation was provided to {.code species} parameter.",
                          "i" = "If not using default species please set {.code species = other}.")
    )
  }

  # Species Spelling Options
  mouse_options <- accepted_names$Mouse_Options
  human_options <- accepted_names$Human_Options
  marmoset_options <- accepted_names$Marmoset_Options
  zebrafish_options <- accepted_names$Zebrafish_Options
  rat_options <- accepted_names$Rat_Options
  drosophila_options <- accepted_names$Drosophila_Options
  macaque_options <- accepted_names$Macaque_Options
  chicken_options <- accepted_names$Chicken_Options

  # Check ensembl vs patterns
  if (isTRUE(x = ensembl_ids) && species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options, chicken_options) && any(!is.null(x = mito_pattern), !is.null(x = ribo_pattern), !is.null(x = mito_features), !is.null(x = ribo_features))) {
    cli_warn(message = c("When using a default species and setting {.code ensembl_ids = TRUE} provided patterns or features are ignored.",
                         "*" = "Supplied {.code mito_pattern}, {.code ribo_pattern}, {.code mito_features}, {.code ribo_features} will be disregarded.")
    )
  }

  # Assign mito/ribo pattern to stored species
  if (species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options, chicken_options) && any(!is.null(x = mito_pattern), !is.null(x = ribo_pattern))) {
    cli_warn(message = c("Pattern expressions for included species are set by default.",
                         "*" = "Supplied {.code mito_pattern} and {.code ribo_pattern} will be disregarded.",
                         "i" = "To override defaults please supply a feature list for mito and/or ribo genes.")
    )
  }

  # default patterns or features
  if (species %in% mouse_options) {
    mito_pattern <- "^mt-"
    ribo_pattern <- "^Rp[sl]"
  }
  if (species %in% human_options) {
    mito_pattern <- "^MT-"
    ribo_pattern <- "^RP[SL]"
  }
  if (species %in% c(marmoset_options, macaque_options, chicken_options)) {
    mito_features <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6")
    ribo_pattern <- "^RP[SL]"
  }
  if (species %in% zebrafish_options) {
    mito_pattern <- "^mt-"
    ribo_pattern <- "^rp[sl]"
  }
  if (species %in% rat_options) {
    mito_pattern <- "^Mt-"
    ribo_pattern <- "^Rp[sl]"
  }
  if (species %in% drosophila_options) {
    mito_pattern <- "^mt:"
    ribo_pattern <- "^Rp[SL]"
  }

  # Check that values are provided for mito and ribo
  if (is.null(x = mito_pattern) && is.null(x = mito_features) && is.null(x = ribo_pattern) && is.null(x = ribo_features)) {
    cli_abort(message = c("No features or patterns provided for mito/ribo genes.",
                          "i" = "Please provide a default species name or pattern/features."))
  }

  # Retrieve ensembl ids if TRUE
  if (isTRUE(x = ensembl_ids)) {
    mito_features <- Retrieve_Ensembl_Mito(species = species)
    ribo_features <- Retrieve_Ensembl_Ribo(species = species)
  }

  all_features <- LIGER_Features(liger_object = object)

  # get features from patterns
  mito_features <- mito_features %||% grep(pattern = mito_pattern, x = all_features, value = TRUE)

  ribo_features <- ribo_features %||% grep(pattern = ribo_pattern, x = all_features, value = TRUE)

  # Check features are present in object
  length_mito_features <- length(x = intersect(x = mito_features, y = all_features))

  length_ribo_features <- length(x = intersect(x = ribo_features, y = all_features))

  # Check length of mito and ribo features found in object
  if (length_mito_features < 1 && length_ribo_features < 1) {
    cli_abort(message = c("No Mito or Ribo features found in object using patterns/feature list provided.",
                          "i" = "Please check pattern/feature list and/or gene names in object.")
    )
  }
  if (length_mito_features < 1) {
    cli_warn(message = c("No Mito features found in object using pattern/feature list provided.",
                         "i" = "No column will be added to meta.data.")
    )
  }
  if (length_ribo_features < 1) {
    cli_warn(message = c("No Ribo features found in object using pattern/feature list provided.",
                         "i" = "No column will be added to meta.data.")
    )
  }

  # Add mito and ribo percent
  if (length_mito_features > 0) {
    good_mito <- mito_features[mito_features %in% all_features]

    if (packageVersion(pkg = 'rliger') > "1.0.1") {
      object <- rliger::runGeneralQC(object = object, mito = FALSE, ribo = FALSE, hemo = FALSE, features = list(mito_name = good_mito), verbose = FALSE)
    } else {
      percent_mito <- unlist(lapply(object@raw.data, function(x) {
        (Matrix::colSums(x[good_mito, ])/Matrix::colSums(x))*100}))
      object@cell.data[ , mito_name] <- percent_mito
    }
  }

  if (length_ribo_features > 0){
    good_ribo <- ribo_features[ribo_features %in% all_features]

    if (packageVersion(pkg = 'rliger') > "1.0.1") {
      object <- rliger::runGeneralQC(object = object, mito = FALSE, ribo = FALSE, hemo = FALSE, features = list(ribo_name = good_ribo), verbose = FALSE)
    } else {
      percent_ribo <- unlist(lapply(object@raw.data, function(x) {
        (Matrix::colSums(x[good_ribo, ])/Matrix::colSums(x))*100}))
      object@cell.data[ , ribo_name] <- percent_ribo
    }
  }

  # Create combined mito ribo column if both present
  if (length_mito_features > 0 && length_ribo_features > 0) {
    if (packageVersion(pkg = 'rliger') > "1.0.1") {
      object@cellMeta[[mito_ribo_name]] <- object@cellMeta[[mito_name]] + object@cellMeta[[ribo_name]]
    } else {
      object_meta <- Fetch_Meta(object = object) %>%
        rownames_to_column("barcodes")

      object_meta <- object_meta %>%
        mutate({{mito_ribo_name}} := .data[[mito_name]] + .data[[ribo_name]])


      object@cell.data[ , mito_ribo_name] <- object_meta[[mito_ribo_name]]
    }
  }

  # return object
  return(object)
}


#' Add Cell Complexity Value
#'
#' @param meta_col_name name to use for new meta data column.  Default is "log10GenesPerUMI".
#' @param overwrite Logical.  Whether to overwrite existing an meta.data column.  Default is FALSE meaning that
#' function will abort if column with name provided to `meta_col_name` is present in meta.data slot.
#'
#' @import cli
#' @importFrom utils packageVersion
#'
#' @method Add_Cell_Complexity liger
#'
#' @export
#' @rdname Add_Cell_Complexity
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' # Liger
#' liger_object <- Add_Cell_Complexity(object = liger_object)
#' }
#'

Add_Cell_Complexity.liger <- function(
  object,
  meta_col_name = "log10GenesPerUMI",
  overwrite = FALSE,
  ...
) {
  # temp liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    cli_abort(message = c("Liger functionality is currently restricted to rliger v1.0.1 or lower.",
                          "i" = "Functionality with rliger v2+ is currently in development."))
  }

  # Check liger
  Is_LIGER(liger_object = object)

  # Check columns for overwrite
  meta_names <- colnames(x = Fetch_Meta(object = object))

  if (meta_col_name %in% meta_names) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Column {.val {meta_col_name}} already present in meta data.",
                            "i" = "*To run function and overwrite column, set parameter {.code overwrite = TRUE} or change respective {.code meta_col_name}*.")
      )
    }
    cli_inform(message = c("Column {.val {meta_col_name}} already present in meta data slot",
                           "i" = "Overwriting those columns as `overwrite = TRUE`.")
    )
  }

  # Add score
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    object@cellMeta[[meta_col_name]] <- log10(object@cellMeta$nGene) / log10(object@cellMeta$nUMI)
  } else {
    object@cell.data[ , meta_col_name] <- log10(object@cell.data$nGene) / log10(object@cell.data$nUMI)
  }

  #return object
  return(object)
}


#' @param species Species of origin for given Object.  If mouse, human, marmoset, zebrafish, rat,
#' drosophila, rhesus macaque, or chicken (name or abbreviation) are provided the function will automatically
#' generate hemo_pattern values.
#' @param hemo_name name to use for the new meta.data column containing percent hemoglobin counts.
#' Default is "percent_hemo".
#' @param hemo_pattern A regex pattern to match features against for hemoglobin genes (will set automatically if
#' species is mouse or human; marmoset features list saved separately).
#' @param hemo_features A list of hemoglobin gene names to be used instead of using regex pattern.
#' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
#' function will abort if columns with any one of the names provided to `hemo_name` is
#' present in meta.data slot.
#' @param list_species_names returns list of all accepted values to use for default species names which
#' contain internal regex/feature lists (human, mouse, marmoset, zebrafish, rat, drosophila, and
#' rhesus macaque).  Default is FALSE.
#'
#' @import cli
#' @importFrom magrittr "%>%"
#'
#' @method Add_Hemo liger
#'
#' @export
#' @rdname Add_Hemo
#'
#' @concept qc_util
#'
#' @examples
#' \dontrun{
#' # Liger
#' liger_object <- Add_Hemo(object = liger_object, species = "human")
#'}
#'

Add_Hemo.liger <- function(
    object,
    species,
    hemo_name = "percent_hemo",
    hemo_pattern = NULL,
    hemo_features = NULL,
    overwrite = FALSE,
    list_species_names = FALSE,
    ...
) {
  # Accepted species names
  accepted_names <- data.frame(
    Mouse_Options = c("Mouse", "mouse", "Ms", "ms", "Mm", "mm"),
    Human_Options = c("Human", "human", "Hu", "hu", "Hs", "hs"),
    Marmoset_Options = c("Marmoset", "marmoset", "CJ", "Cj", "cj", NA),
    Zebrafish_Options = c("Zebrafish", "zebrafish", "DR", "Dr", "dr", NA),
    Rat_Options = c("Rat", "rat", "RN", "Rn", "rn", NA),
    Drosophila_Options = c("Drosophila", "drosophila", "DM", "Dm", "dm", NA),
    Macaque_Options = c("Macaque", "macaque", "Rhesus", "macaca", "mmulatta", NA),
    Chicken_Options = c("Chicken", "chicken", "Gallus", "gallus", "Gg", "gg")
  )

  # Return list of accepted default species name options
  if (isTRUE(x = list_species_names)) {
    return(accepted_names)
    stop_quietly()
  }

  # Check liger
  Is_LIGER(liger_object = object)

  # Overwrite check
  meta_names <- colnames(x = Fetch_Meta(object = object))

  if (hemo_name %in% meta_names) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Columns with {.val {hemo_name}} already present in meta data.",
                            "i" = "*To run function and overwrite columns set parameter {.code overwrite = TRUE} or change {.code hemo_name}.*")
      )
    }
    cli_inform(message = c("Columns with {.val {hemo_name}} already present in meta data.",
                           "i" = "Overwriting those columns as {.code overwrite = TRUE}.")
    )
  }

  # Checks species
  if (is.null(x = species)) {
    cli_abort(message = c("No species name or abbreivation was provided to {.code species} parameter.",
                          "i" = "If not using default species please set {.code species = other}.")
    )
  }

  # Species Spelling Options
  mouse_options <- accepted_names$Mouse_Options
  human_options <- accepted_names$Human_Options
  marmoset_options <- accepted_names$Marmoset_Options
  zebrafish_options <- accepted_names$Zebrafish_Options
  rat_options <- accepted_names$Rat_Options
  drosophila_options <- accepted_names$Drosophila_Options
  macaque_options <- accepted_names$Macaque_Options
  chicken_options <- accepted_names$Chicken_Options

  # Assign hemo pattern to stored species
  if (species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options, macaque_options, chicken_options) && any(!is.null(x = hemo_pattern))) {
    cli_warn(message = c("Pattern expressions for included species are set by default.",
                         "*" = "Supplied {.code hemo_pattern} and {.code hemo_pattern} will be disregarded.",
                         "i" = "To override defaults please supply a feature list for hemo genes.")
    )
  }

  if (species %in% mouse_options) {
    species_use <- "Mouse"
    hemo_pattern <- "^Hb[^(P)]"
  }
  if (species %in% human_options) {
    species_use <- "Human"
    hemo_pattern <- "^HB[^(P)]"
  }
  if (species %in% c(marmoset_options, macaque_options)) {
    species_use <- "Marmoset/Macaque"
    hemo_pattern <- "^^HB[^(P)]"
  }
  if (species %in% zebrafish_options) {
    species_use <- "Zebrafish"
    hemo_pattern <- "^hb[^(P)]"
  }
  if (species %in% rat_options) {
    species_use <- "Rat"
    hemo_pattern <- "^Hb[^(P)]"
  }
  if (species %in% drosophila_options) {
    species_use <- "Drosophila"
    hemo_pattern <- "^glob"
  }
  if (species %in% chicken_options) {
    species_use <- "Chicken"
    hemo_pattern <- "^HB[^(P)]"
  }

  # Check that values are provided for mito and ribo
  if (is.null(x = hemo_pattern) && is.null(x = hemo_features)) {
    cli_abort(message = c("No features or patterns provided for hemo genes.",
                          "i" = "Please provide a default species name or pattern/features."))
  }

  # get all features
  all_features <- LIGER_Features(liger_object = object, by_dataset = FALSE)

  # get features from patterns
  hemo_features <- hemo_features %||% grep(pattern = hemo_pattern, x = all_features, value = TRUE)

  # Check features are present in object
  length_hemo_features <- length(x = intersect(x = hemo_features, y = all_features))

  # Check length of hemo features found in object
  if (length_hemo_features < 1) {
    cli_warn(message = c("No hemoglobin features found in object using pattern/feature list provided.",
                         "i" = "No column will be added to meta.data.")
    )
  }

  # Add hemo column
  cli_inform(message = "Adding Percent Hemoglobin for {.field {species_use}} using gene symbol pattern: {.val {hemo_pattern}}.")
  if (length_hemo_features > 0) {
    good_hemo <- hemo_features[hemo_features %in% all_features]

    if (packageVersion(pkg = 'rliger') > "1.0.1") {
      object <- rliger::runGeneralQC(object = object, mito = FALSE, ribo = FALSE, hemo = FALSE, features = list(hemo_name = good_hemo), verbose = FALSE)
    } else {
      percent_hemo <- unlist(lapply(object@raw.data, function(x) {
        (Matrix::colSums(x[good_hemo, ])/Matrix::colSums(x))*100}))
      object@cell.data[ , hemo_name] <- percent_hemo
    }
  }

  # return final object
  return(object)
}


#' @param num_top_genes An integer vector specifying the size(s) of the top set of high-abundance genes.
#' Used to compute the percentage of library size occupied by the most highly expressed genes in each cell.
#' @param meta_col_name name to use for new meta data column.  Default is "percent_topXX", where XX is
#' equal to the value provided to `num_top_genes`.
#' @param overwrite Logical.  Whether to overwrite existing an meta.data column.  Default is FALSE meaning that
#' function will abort if column with name provided to `meta_col_name` is present in meta.data slot.
#' @param verbose logical, whether to print messages with status updates, default is TRUE.
#'
#' @import cli
#' @importFrom dplyr select all_of bind_rows
#' @importFrom magrittr "%>%"
#' @importFrom rlang is_installed
#'
#' @return A liger Object
#'
#' @method Add_Top_Gene_Pct liger
#'
#' @export
#' @rdname Add_Top_Gene_Pct
#'
#' @concept qc_util
#'
#' @examples
#' \dontrun{
#' liger_object <- Add_Top_Gene_Pct(object = liger_object, num_top_genes = 50)
#' }
#'

Add_Top_Gene_Pct.liger <- function(
    object,
    num_top_genes = 50,
    meta_col_name = NULL,
    overwrite = FALSE,
    verbose = TRUE,
    ...
){
  # Check for scuttle first
  scuttle_check <- is_installed(pkg = "scuttle")
  if (isFALSE(x = scuttle_check)) {
    cli_abort(message = c(
      "Please install the {.val scuttle} package to calculate/add top {num_top_genes} genes percentage.",
      "i" = "This can be accomplished with the following commands: ",
      "----------------------------------------",
      "{.field `install.packages({symbol$dquote_left}BiocManager{symbol$dquote_right})`}",
      "{.field `BiocManager::install({symbol$dquote_left}scuttle{symbol$dquote_right})`}",
      "----------------------------------------"
    ))
  }

  # Check Liger
  Is_LIGER(liger_object = object)

  # Set colnames
  scuttle_colname <- paste0("percent.top_", num_top_genes)
  if (is.null(x = meta_col_name)) {
    meta_col_name <- paste0("percent_top", num_top_genes)
  }

  # Check columns for overwrite
  meta_names <- colnames(x = Fetch_Meta(object = object))

  if (meta_col_name %in% meta_names) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Column {.val {meta_col_name}} already present in meta data.",
                            "i" = "*To run function and overwrite column, set parameter {.code overwrite = TRUE} or change respective {.code meta_col_name}*.")
      )
    }
    cli_inform(message = c("Column {.val {meta_col_name}} already present in meta data.",
                           "i" = "Overwriting those columns as {.code overwrite = TRUE}.")
    )
  }

  # Get number of datasets
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    num_datasets <- length(x = object@datasets)
  } else {
    num_datasets <- length(x = object@raw.data)
  }

  # Extract matrix
  if (isTRUE(x = verbose)) {
    cli_inform(message = "Calculating percent expressing top {num_top_genes} across all datasets.")
  }

  # apply over all datasets
  res_list <- lapply(1:num_datasets, function(x) {
    if (packageVersion(pkg = 'rliger') > "1.0.1") {
      dataset_mat <- rliger::getMatrix(x = object, slot = "rawData")[[x]]
    } else {
      dataset_mat <- object@raw.data[[x]]
    }

    # run scuttle
    dataset_res <- as.data.frame(scuttle::perCellQCMetrics(x = dataset_mat, percent.top = num_top_genes))
    # select results column
    dataset_res <- dataset_res %>%
      select(all_of(scuttle_colname))
  })

  # combine results
  if (isTRUE(x = verbose)) {
    cli_inform(message = "Combining data from all datasets.")
  }
  res <- bind_rows(res_list)

  # Add to object and return
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    object@cellMeta[[meta_col_name]] <- res
  } else {
    object@cell.data[ , meta_col_name] <- res
  }

  # return object
  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### ANALYSIS UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Perform variable gene selection over whole dataset
#'
#' Performs variable gene selection for LIGER object across the entire object instead of by
#' dataset and then taking union.
#'
#' @param liger_object LIGER object name.
#' @param num_genes Number of genes to find. Optimizes the value of `var.thresh`  to get
#' this number of genes, (Default is NULL).
#' @param var.thresh Variance threshold. Main threshold used to identify variable genes.
#' Genes with expression variance greater than threshold (relative to mean) are selected.
#' (higher threshold -> fewer selected genes).
#' @param alpha.thresh Alpha threshold. Controls upper bound for expected mean gene
#' expression (lower threshold -> higher upper bound). (default 0.99)
#' @param tol Tolerance to use for optimization if num.genes values passed in (default 0.0001).
#' @param do.plot Display log plot of gene variance vs. gene expression. Selected genes are
#' plotted in green. (Default FALSE)
#' @param pt.size Point size for plot.
#' @param chunk size of chunks in hdf5 file. (Default 1000)
#'
#' @return A LIGER Object with variable genes in correct slot.
#'
#' @import cli
#' @importFrom utils packageVersion
#'
#' @references Matching function parameter text descriptions are taken from `rliger::selectGenes`
#' which is called by this function after creating new temporary object/dataset.
#' \url{https://github.com/welch-lab/liger}. (License: GPL-3).
#'
#' @export
#'
#' @concept liger_object_util
#'
#' @examples
#' \dontrun{
#' liger_obj <- Variable_Features_ALL_LIGER(liger_object = liger_obj, num_genes = 2000)
#' }
#'

Variable_Features_ALL_LIGER <- function(
  liger_object,
  num_genes = NULL,
  var.thresh = 0.3,
  alpha.thresh = 0.99,
  tol = 0.0001,
  do.plot = FALSE,
  pt.size = 0.3,
  chunk=1000
) {
  # temp liger version check
  if (packageVersion(pkg = 'rliger') > "1.0.1") {
    cli_abort(message = c("Liger functionality is currently restricted to rliger v1.0.1 or lower.",
                          "i" = "Functionality with rliger v2+ is currently in development."))
  }

  Is_LIGER(liger_object = liger_object)

  raw_data <- liger_object@raw.data

  cli_inform(message = "Creating temporary object with combined data.")

  temp_liger <- rliger::createLiger(raw.data = list("dataset" = Merge_Sparse_Data_All(raw_data)), remove.missing = FALSE)

  rm(raw_data)
  gc()

  cli_inform(message = "Normalizing and identifying variable features.")

  temp_liger <- rliger::normalize(object = temp_liger)
  temp_liger <- rliger::selectGenes(object = temp_liger, var.thresh = var.thresh, do.plot = do.plot, num.genes = num_genes, tol = tol, alpha.thresh = alpha.thresh, cex.use = pt.size, chunk = chunk)
  var_genes <- temp_liger@var.genes

  rm(temp_liger)
  gc()

  liger_object@var.genes <- var_genes
  return(liger_object)
}
