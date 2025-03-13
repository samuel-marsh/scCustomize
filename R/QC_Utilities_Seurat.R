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
#' @param add_IEG_module_score logical, whether to add module score belonging to IEG genes to object (Default is TRUE).
#' @param add_hemo logical, whether to add percentage of counts belonging to homoglobin genes to object (Default is TRUE).
#' @param add_cell_cycle logical, whether to addcell cycle scores and phase based on
#' \code{\link[Seurat]{CellCycleScoring}}.  Only applicable if `species = "human"`.  (Default is TRUE).
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
#' @param ieg_module_name name to use for new meta data column for module score of IEGs.  Default is "ieg_score".
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
#' @param assay assay to use in calculation.  Default is "RNA".  *Note* This should only be changed if
#' storing corrected and uncorrected assays in same object (e.g. outputs of both Cell Ranger and Cell Bender).
#' @param list_species_names returns list of all accepted values to use for default species names which
#' contain internal regex/feature lists (human, mouse, marmoset, zebrafish, rat, drosophila, rhesus macaque, and
#' chicken).  Default is FALSE.
#' @param overwrite Logical.  Whether to overwrite existing an meta.data column.  Default is FALSE meaning that
#' function will abort if column with name provided to `meta_col_name` is present in meta.data slot.
#'
#' @import cli
#' @importFrom SeuratObject Layers
#'
#' @return A Seurat Object
#'
#' @method Add_Cell_QC_Metrics Seurat
#'
#' @export
#' @rdname Add_Cell_QC_Metrics
#'
#' @concept qc_util
#'
#' @examples
#' \dontrun{
#' obj <- Add_Cell_QC_Metrics(object = obj, species = "Human")
#'}
#'

Add_Cell_QC_Metrics.Seurat <- function(
    object,
    species,
    add_mito_ribo = TRUE,
    add_complexity = TRUE,
    add_top_pct = TRUE,
    add_MSigDB = TRUE,
    add_IEG = TRUE,
    add_IEG_module_score = TRUE,
    add_hemo = TRUE,
    add_cell_cycle = TRUE,
    mito_name = "percent_mito",
    ribo_name = "percent_ribo",
    mito_ribo_name = "percent_mito_ribo",
    complexity_name = "log10GenesPerUMI",
    top_pct_name = NULL,
    oxphos_name = "percent_oxphos",
    apop_name = "percent_apop",
    dna_repair_name = "percent_dna_repair",
    ieg_name = "percent_ieg",
    ieg_module_name = "ieg_score",
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
    list_species_names = FALSE,
    overwrite = FALSE,
    ...
) {
  # Set assay
  assay <- assay %||% DefaultAssay(object = object)

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

  # Species Spelling Options
  mouse_options <- accepted_names$Mouse_Options
  human_options <- accepted_names$Human_Options
  marmoset_options <- accepted_names$Marmoset_Options
  zebrafish_options <- accepted_names$Zebrafish_Options
  rat_options <- accepted_names$Rat_Options
  drosophila_options <- accepted_names$Drosophila_Options
  macaque_options <- accepted_names$Macaque_Options
  chicken_options <- accepted_names$Chicken_Options

  # Add mito/ribo
  if (isTRUE(x = add_mito_ribo)) {
    cli_inform(message = c("*" = "Adding {.field Mito/Ribo Percentages} to meta.data."))
    object <- Add_Mito_Ribo(object = object, species = species, mito_name = mito_name, ribo_name = ribo_name, mito_ribo_name = mito_ribo_name, mito_pattern = mito_pattern, ribo_pattern = ribo_pattern, mito_features = mito_features, ribo_features = ribo_features, ensembl_ids = ensembl_ids, assay = assay, overwrite = overwrite)
  }

  # Add complexity
  if (isTRUE(x = add_complexity)) {
    cli_inform(message = c("*" = "Adding {.field Cell Complexity #1 (log10GenesPerUMI)} to meta.data."))
    object <- Add_Cell_Complexity(object = object, meta_col_name = complexity_name, assay = assay, overwrite = overwrite)
  }

  # Add top gene expression percent
  if (isTRUE(x = add_top_pct)) {
    cli_inform(message = c("*" = "Adding {.field Cell Complexity #2 (Top {num_top_genes} Percentages)} to meta.data."))
    object <- Add_Top_Gene_Pct(object = object, num_top_genes = num_top_genes, meta_col_name = top_pct_name, assay = assay, overwrite = overwrite)
  }

  # Add MSigDB
  if (isTRUE(x = add_MSigDB)) {
    if (species %in% marmoset_options) {
      cli_warn(message = c("{.val Marmoset} is not currently a part of MSigDB gene list database.",
                           "i" = "No columns will be added to object meta.data"))
    } else {
      cli_inform(message = c("*" = "Adding {.field MSigDB Oxidative Phosphorylation, Apoptosis, and DNA Repair Percentages} to meta.data."))
      object <- Add_MSigDB_Seurat(seurat_object = object, species = species, oxphos_name = oxphos_name, apop_name = apop_name, dna_repair_name = dna_repair_name, assay = assay, overwrite = overwrite, ensembl_ids = ensembl_ids)
    }
  }

  # Add IEG
  if (isTRUE(x = add_IEG)) {
    if (species %in% c(marmoset_options, rat_options, zebrafish_options, macaque_options, drosophila_options, chicken_options)) {
      cli_warn(message = c("{.val Rat, Marmoset, Macaque, Zebrafish, Drosophila, Chicken} are not currently supported.",
                           "i" = "No column will be added to object meta.data"))
    } else {
      if (isTRUE(x = add_IEG_module_score)) {
        cli_inform(message = c("*" = "Adding {.field IEG Percentages} & {.field IEG Module Score} to meta.data."))
      } else {
        cli_inform(message = c("*" = "Adding {.field IEG Percentages} to meta.data."))
      }
      object <- Add_IEG_Seurat(seurat_object = object, species = species, ieg_name = ieg_name, assay = assay, overwrite = overwrite, ensembl_ids = ensembl_ids, ieg_module_score = add_IEG_module_score, ieg_module_name = ieg_module_name)
    }
  }

  # Add hemo
  if (isTRUE(x = add_hemo)) {
    cli_inform(message = c("*" = "Adding {.field Hemoglobin Percentages} to meta.data."))
    object <- Add_Hemo(object = object, species = species, hemo_name = hemo_name, hemo_pattern = hemo_pattern, hemo_features = hemo_features, assay = assay, overwrite = overwrite)
  }

  # Add cell cycle
  if (isTRUE(x = add_cell_cycle)) {
    if (!species %in% human_options) {
      cli_warn(message = c("x" = "Cell Cycle Scoring is only supported for human in this function.",
                           "i" = "To add score for other species, use {.code Seurat::CellCycleScoring} function separately with correct species cell cycle gene list."
      ))
    } else {
      cli_inform(message = c("*" = "Adding {.field Cell Cycle Scoring} to meta.data."))
      if (isFALSE(x = Check_Normalized(object = object, assay = assay, error = FALSE))) {
        cli_inform(message = c("Layer with normalized data not present.",
                               "i" = "Normalizing Data."))
        object <- NormalizeData(object = object)
      }

      # Overwrite check
      if ("S.Score" %in% colnames(x = object@meta.data) || "G2M.Score" %in% colnames(x = object@meta.data) || "Phase" %in% colnames(x = object@meta.data)) {
        if (isFALSE(x = overwrite)) {
          cli_abort(message = c("Columns with {.val S.Score}, {.val G2M.Score} and/or {.val Phase} already present in meta.data slot.",
                                "i" = "*To run function and overwrite columns set parameter {.code overwrite = TRUE}*")
          )
        }
        cli_inform(message = c("Columns with {.val S.Score}, {.val G2M.Score} and/or {.val Phase} already present in meta.data slot.",
                               "i" = "Overwriting those columns as {.code overwrite = TRUE.}")
        )
      }

      # Add Cell Cycle Scoring
      cli_inform(message = "Calculating {.field Cell Cycle Scores}.")
      object <- CellCycleScoring(object = object, s.features = Seurat::cc.genes.updated.2019$s.genes, g2m.features = Seurat::cc.genes.updated.2019$g2m.genes)
    }
  }

  # Log Command
  object <- LogSeuratCommand(object = object)

  # return object
  return(object)
}


#' @param species Species of origin for given Seurat Object.  If mouse, human, marmoset, zebrafish, rat,
#' drosophila, rhesus macaque, or chicken (name or abbreviation) are provided the function will automatically
#' generate mito_pattern and ribo_pattern values.
#' @param mito_name name to use for the new meta.data column containing percent mitochondrial counts.
#' Default is "percent_mito".
#' @param ribo_name name to use for the new meta.data column containing percent ribosomal counts.
#' Default is "percent_ribo".
#' @param mito_ribo_name name to use for the new meta.data column containing percent
#' mitochondrial+ribosomal counts.  Default is "percent_mito_ribo".
#' @param mito_pattern A regex pattern to match features against for mitochondrial genes (will set automatically if
#' species is mouse, human, zebrafish, rat, drosophila, rhesus macaque, or chicken;
#' marmoset features list saved separately).
#' @param ribo_pattern A regex pattern to match features against for ribosomal genes
#' (will set automatically if species is mouse, human, marmoset, zebrafish, rat,
#' drosophila, rhesus macaque, or chicken).
#' @param mito_features A list of mitochondrial gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param ribo_features A list of ribosomal gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param ensembl_ids logical, whether feature names in the object are gene names or
#' ensembl IDs (default is FALSE; set TRUE if feature names are ensembl IDs).
#' @param assay Assay to use (default is the current object default assay).
#' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
#' function will abort if columns with any one of the names provided to `mito_name` `ribo_name` or
#' `mito_ribo_name` is present in meta.data slot.
#' @param list_species_names returns list of all accepted values to use for default species names which
#' contain internal regex/feature lists (human, mouse, marmoset, zebrafish, rat, drosophila, rhesus macaque, and
#' chicken).  Default is FALSE.
#' @param species_prefix the species prefix in front of gene symbols in object if providing two species for
#' multi-species aligned dataset.
#'
#' @import cli
#' @importFrom dplyr mutate select intersect all_of
#' @importFrom magrittr "%>%"
#' @importFrom rlang ":="
#' @importFrom Seurat PercentageFeatureSet AddMetaData
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @method Add_Mito_Ribo Seurat
#'
#' @export
#' @rdname Add_Mito_Ribo
#'
#' @concept qc_util
#'
#' @examples
#' \dontrun{
#' # Seurat
#' seurat_object <- Add_Mito_Ribo(object = seurat_object, species = "human")
#'}
#'

Add_Mito_Ribo.Seurat <- function(
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
    assay = NULL,
    overwrite = FALSE,
    list_species_names = FALSE,
    species_prefix = NULL,
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

  # Check Seurat
  Is_Seurat(seurat_object = object)

  # Check name collision
  if (any(duplicated(x = c(mito_name, ribo_name, mito_ribo_name)))) {
    cli_abort(message = "One or more of values provided to {.code mito_name}, {.code ribo_name}, {.code mito_ribo_name} are identical.")
  }

  # Overwrite check
  if (mito_name %in% colnames(x = object@meta.data) || ribo_name %in% colnames(x = object@meta.data) || mito_ribo_name %in% colnames(x = object@meta.data)) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Columns with {.val {mito_name}} and/or {.val {ribo_name}} already present in meta.data slot.",
                            "i" = "*To run function and overwrite columns set parameter {.code overwrite = TRUE} or change respective {.code mito_name}, {.code ribo_name}, and/or {.code mito_ribo_name}*")
      )
    }
    cli_inform(message = c("Columns with {.val {mito_name}} and/or {.val {ribo_name}} already present in meta.data slot.",
                           "i" = "Overwriting those columns as {.code overwrite = TRUE.}")
    )
  }

  # Checks species
  if (is.null(x = species)) {
    cli_abort(message = c("No species name or abbreivation was provided to {.code species} parameter.",
                          "i" = "If not using default species please set {.code species = other}.")
    )
  }

  # Dual species checks
  if (length(x = species) > 1 && length(x = species) != length(x = species_prefix)) {
    cli_abort(message = "The length of {.code species} must be equal to length of {.code species_prefix}.")
  }

  # Set default assay
  assay <- assay %||% DefaultAssay(object = object)

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
  if (isTRUE(x = ensembl_ids) && all(species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options, chicken_options) && any(!is.null(x = mito_pattern)), !is.null(x = ribo_pattern), !is.null(x = mito_features), !is.null(x = ribo_features))) {
    cli_warn(message = c("When using a default species and setting {.code ensembl_ids = TRUE} provided patterns or features are ignored.",
                         "*" = "Supplied {.code mito_pattern}, {.code ribo_pattern}, {.code mito_features}, {.code ribo_features} will be disregarded.")
    )
  }

  # Assign mito/ribo pattern to stored species
  if (all(species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options, chicken_options)) && any(!is.null(x = mito_pattern), !is.null(x = ribo_pattern))) {
    cli_warn(message = c("Pattern expressions for included species are set by default.",
                         "*" = "Supplied {.code mito_pattern} and {.code ribo_pattern} will be disregarded.",
                         "i" = "To override defaults please supply a feature list for mito and/or ribo genes.")
    )
  }

  if (length(x = species) == 1) {
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

    mito_features <- mito_features %||% grep(pattern = mito_pattern, x = rownames(x = object[[assay]]), value = TRUE)

    ribo_features <- ribo_features %||% grep(pattern = ribo_pattern, x = rownames(x = object[[assay]]), value = TRUE)

    # Check features are present in object
    length_mito_features <- length(x = intersect(x = mito_features, y = rownames(x = object[[assay]])))

    length_ribo_features <- length(x = intersect(x = ribo_features, y = rownames(x = object[[assay]])))

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
  } else {
    # get dual species gene lists
    mito_features <- Retrieve_Dual_Mito_Features(object = object, species = species, species_prefix = species_prefix, assay = assay)

    ribo_features <- Retrieve_Dual_Ribo_Features(object = object, species = species, species_prefix = species_prefix, assay = assay)

    # Check features are present in object
    length_mito_features <- length(x = intersect(x = mito_features, y = rownames(x = object[[assay]])))

    length_ribo_features <- length(x = intersect(x = ribo_features, y = rownames(x = object[[assay]])))
  }

  # Add mito and ribo columns
  cli_inform(message = "Adding Percent Mitochondrial genes for {.field {species}} using gene symbol pattern: {.val {mito_pattern}}.")
  if (length_mito_features > 0) {
    good_mito <- mito_features[mito_features %in% rownames(x = object)]
    object[[mito_name]] <- PercentageFeatureSet(object = object, features = good_mito, assay = assay)
  }
  if (length_ribo_features > 0) {
    cli_inform(message = "Adding Percent Ribosomal genes for {.field {species}} using gene symbol pattern: {.val {ribo_pattern}}.")
    good_ribo <- ribo_features[ribo_features %in% rownames(x = object)]
    object[[ribo_name]] <- PercentageFeatureSet(object = object, features = good_ribo, assay = assay)
  }

  # Create combined mito ribo column if both present
  if (length_mito_features > 0 && length_ribo_features > 0) {
    cli_inform(message = "Adding Percent Mito+Ribo by adding Mito & Ribo percentages.")
    object_meta <- Fetch_Meta(object = object) %>%
      rownames_to_column("temp_barcodes")

    object_meta <- object_meta %>%
      mutate({{mito_ribo_name}} := .data[[mito_name]] + .data[[ribo_name]])

    object_meta <- object_meta %>%
      select(all_of(c("temp_barcodes", mito_ribo_name))) %>%
      column_to_rownames("temp_barcodes")

    object <- AddMetaData(object = object, metadata = object_meta)
  }

  # Log Command
  object <- LogSeuratCommand(object = object)

  # return final object
  return(object)
}


#' @param species Species of origin for given Seurat Object.  If mouse, human, marmoset, zebrafish, rat,
#' drosophila, rhesus macaque, or chicken (name or abbreviation) are provided the function will automatically
#' generate hemo_pattern values.
#' @param hemo_name name to use for the new meta.data column containing percent hemoglobin counts.
#' Default is "percent_hemo".
#' @param hemo_pattern A regex pattern to match features against for hemoglobin genes (will set automatically if
#' species is mouse or human; marmoset features list saved separately).
#' @param hemo_features A list of hemoglobin gene names to be used instead of using regex pattern.
#' @param ensembl_ids logical, whether feature names in the object are gene names or
#' ensembl IDs (default is FALSE; set TRUE if feature names are ensembl IDs).
#' @param assay Assay to use (default is the current object default assay).
#' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
#' function will abort if columns with any one of the names provided to `hemo_name` is
#' present in meta.data slot.
#' @param list_species_names returns list of all accepted values to use for default species names which
#' contain internal regex/feature lists (human, mouse, marmoset, zebrafish, rat, drosophila, and
#' rhesus macaque).  Default is FALSE.
#'
#' @import cli
#' @importFrom dplyr mutate select intersect all_of
#' @importFrom magrittr "%>%"
#' @importFrom rlang ":="
#' @importFrom Seurat PercentageFeatureSet AddMetaData
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @method Add_Hemo Seurat
#'
#' @export
#' @rdname Add_Hemo
#'
#' @concept qc_util
#'
#' @examples
#' \dontrun{
#' # Seurat
#' seurat_object <- Add_Hemo(object = seurat_object, species = "human")
#'}
#'

Add_Hemo.Seurat <- function(
    object,
    species,
    hemo_name = "percent_hemo",
    hemo_pattern = NULL,
    hemo_features = NULL,
    ensembl_ids = FALSE,
    assay = NULL,
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

  # Check Seurat
  Is_Seurat(seurat_object = object)

  # Overwrite check
  if (hemo_name %in% colnames(x = object@meta.data)) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Columns with {.val {hemo_name}} already present in meta.data slot.",
                            "i" = "*To run function and overwrite columns set parameter {.code overwrite = TRUE} or change {.code hemo_name}*")
      )
    }
    cli_inform(message = c("Columns with {.val {hemo_name}} already present in meta.data slot.",
                           "i" = "Overwriting column as {.code overwrite = TRUE.}")
    )
  }

  # Checks species
  if (is.null(x = species)) {
    cli_abort(message = c("No species name or abbreivation was provided to {.code species} parameter.",
                          "i" = "If not using default species please set {.code species = other}.")
    )
  }

  # Set default assay
  assay <- assay %||% DefaultAssay(object = object)

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
    hemo_pattern <- "^HB[^(P)]"
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
    cli_abort(message = c("No features or patterns provided for hemoglobin genes.",
                          "i" = "Please provide a default species name or pattern/features."))
  }

  # Retrieve ensembl ids if TRUE
  if (isTRUE(x = ensembl_ids)) {
    hemo_features <- Retrieve_Ensembl_Hemo(species = species)
  }

  hemo_features <- hemo_features %||% grep(pattern = hemo_pattern, x = rownames(x = object[[assay]]), value = TRUE)

  # Check features are present in object
  length_hemo_features <- length(x = intersect(x = hemo_features, y = rownames(x = object[[assay]])))

  # Check length of hemo features found in object
  if (length_hemo_features < 1) {
    cli_warn(message = c("No hemoglobin features found in object using pattern/feature list provided.",
                         "i" = "No column will be added to meta.data.")
    )
  }

  # Add hemo columns
  cli_inform(message = "Adding Percent Hemoglobin for {.field {species_use}} using gene symbol pattern: {.val {hemo_pattern}}.")
  if (length_hemo_features > 0) {
    good_hemo <- hemo_features[hemo_features %in% rownames(x = object)]
    object[[hemo_name]] <- PercentageFeatureSet(object = object, features = good_hemo, assay = assay)
  }

  # Log Command
  object <- LogSeuratCommand(object = object)

  # return final object
  return(object)
}


#' Add Cell Complexity Value
#'
#' @param meta_col_name name to use for new meta data column.  Default is "log10GenesPerUMI".
#' @param assay assay to use in calculation.  Default is "RNA".  *Note* This should only be changed if
#' storing corrected and uncorrected assays in same object (e.g. outputs of both Cell Ranger and Cell Bender).
#' @param overwrite Logical.  Whether to overwrite existing an meta.data column.  Default is FALSE meaning that
#' function will abort if column with name provided to `meta_col_name` is present in meta.data slot.
#'
#' @import cli
#'
#' @method Add_Cell_Complexity Seurat
#'
#' @export
#' @rdname Add_Cell_Complexity
#'
#' @concept qc_util
#'
#' @examples
#' # Seurat
#' library(Seurat)
#' pbmc_small <- Add_Cell_Complexity(object = pbmc_small)
#'

Add_Cell_Complexity.Seurat <- function(
    object,
    meta_col_name = "log10GenesPerUMI",
    assay = "RNA",
    overwrite = FALSE,
    ...
) {
  # Check Seurat
  Is_Seurat(seurat_object = object)

  # Add assay warning message
  if (assay != "RNA") {
    cli_warn(message = "Assay is set to value other than 'RNA'. This should only be done in rare instances.  See documentation for more info ({.code ?Add_Cell_Complexity_Seurat}).",
             .frequency = "once",
             .frequency_id = "assay_warn")
  }

  # Check columns for overwrite
  if (meta_col_name %in% colnames(x = object@meta.data)) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Column {.val {meta_col_name}} already present in meta.data slot.",
                            "i" = "*To run function and overwrite column, set parameter {.code overwrite = TRUE} or change respective {.code meta_col_name}*.")
      )
    }
    cli_inform(message = c("Column {.val {meta_col_name}} already present in meta.data slot",
                           "i" = "Overwriting those columns as {.code overwrite = TRUE}.")
    )
  }

  # variable names
  feature_name <- paste0("nFeature_", assay)
  count_name <- paste0("nCount_", assay)

  # Add score
  object[[meta_col_name]] <- log10(object[[feature_name]]) / log10(object[[count_name]])

  # Log Command
  object <- LogSeuratCommand(object = object)

  #return object
  return(object)
}


#' @param num_top_genes An integer vector specifying the size(s) of the top set of high-abundance genes.
#' Used to compute the percentage of library size occupied by the most highly expressed genes in each cell.
#' @param meta_col_name name to use for new meta data column.  Default is "percent_topXX", where XX is
#' equal to the value provided to `num_top_genes`.
#' @param assay assay to use in calculation.  Default is "RNA".  *Note* This should only be changed if
#' storing corrected and uncorrected assays in same object (e.g. outputs of both Cell Ranger and Cell Bender).
#' @param overwrite Logical.  Whether to overwrite existing an meta.data column.  Default is FALSE meaning that
#' function will abort if column with name provided to `meta_col_name` is present in meta.data slot.
#' @param verbose logical, whether to print messages with status updates, default is TRUE.
#'
#' @import cli
#' @importFrom dplyr select all_of bind_rows
#' @importFrom magrittr "%>%"
#' @importFrom rlang is_installed
#' @importFrom SeuratObject LayerData
#'
#' @return A Seurat Object
#'
#' @method Add_Top_Gene_Pct Seurat
#'
#' @export
#' @rdname Add_Top_Gene_Pct
#'
#' @concept qc_util
#'
#' @references This function uses scuttle package (license: GPL-3) to calculate the percent of expression
#' coming from top XX genes in each cell.  Parameter description for `num_top_genes` also from scuttle.
#' If using this function in analysis, in addition to citing scCustomize, please cite scuttle:
#' McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater: pre-processing, quality control,
#' normalisation and visualisation of single-cell RNA-seq data in R.” Bioinformatics, 33, 1179-1186.
#' \url{doi:10.1093/bioinformatics/btw777}.
#' @seealso \url{https://bioconductor.org/packages/release/bioc/html/scuttle.html}
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' pbmc_small <- Add_Top_Gene_Pct(seurat_object = pbmc_small, num_top_genes = 50)
#' }
#'

Add_Top_Gene_Pct.Seurat <- function(
    object,
    num_top_genes = 50,
    meta_col_name = NULL,
    assay = "RNA",
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

  # Check Seurat
  Is_Seurat(seurat_object = object)

  # Add assay warning message
  if (assay != "RNA") {
    cli_warn(message = "Assay is set to value other than 'RNA'. This should only be done in rare instances.  See documentation for more info ({.code ?Add_Top_Gene_Pct_Seurat}).",
             .frequency = "once",
             .frequency_id = "assay_warn")
  }

  # Set colnames
  scuttle_colname <- paste0("percent.top_", num_top_genes)
  if (is.null(x = meta_col_name)) {
    meta_col_name <- paste0("percent_top", num_top_genes)
  }

  # Check columns for overwrite
  if (meta_col_name %in% colnames(x = object@meta.data)) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Column {.val {meta_col_name}} already present in meta.data slot.",
                            "i" = "*To run function and overwrite column, set parameter {.code overwrite = TRUE} or change respective {.code meta_col_name}*.")
      )
    }
    cli_inform(message = c("Column {.val {meta_col_name}} already present in meta.data slot",
                           "i" = "Overwriting those columns as {.code overwrite = TRUE}.")
    )
  }

  count_layers_present <- Layers(object = object, search = "counts")

  # Extract matrix
  if (length(x = count_layers_present) == 1) {
    if (isTRUE(x = verbose)) {
      cli_inform(message = "Calculating percent expressing top {num_top_genes} for layer: {.field {count_layers_present}}")
    }

    count_mat <- LayerData(object = object, assay = assay, layer = "counts")

    # calculate
    res <- as.data.frame(scuttle::perCellQCMetrics(x = count_mat, percent.top = num_top_genes))

    # select percent column
    res <- res %>%
      select(all_of(scuttle_colname))
  }


  if (length(x = count_layers_present) > 1) {
    res_list <- lapply(1:length(x = count_layers_present), function(x) {
      if (isTRUE(x = verbose)) {
        cli_inform(message = "Calculating percent expressing top {num_top_genes} for layer: {.field {count_layers_present[x]}}")
      }

      # Get layer data
      layer_count <- LayerData(object = object, assay = assay, layer = count_layers_present[x])

      # run scuttle
      layer_res <- as.data.frame(scuttle::perCellQCMetrics(x = layer_count, percent.top = num_top_genes))
      # select results column
      layer_res <- layer_res %>%
        select(all_of(scuttle_colname))
    })

    # combine results
    if (isTRUE(x = verbose)) {
      cli_inform(message = "Combining data from: {.field {count_layers_present}}")
    }
    res <- bind_rows(res_list)
  }

  # Add to object and return
  object <- AddMetaData(object = object, metadata = res, col.name = meta_col_name)

  # Log Command
  object <- LogSeuratCommand(object = object)

  return(object)
}


#' @param species Species of origin for given Seurat Object.  Only accepted species are: mouse, human (name or abbreviation).
#' @param sample_col column name in meta.data that contains sample ID information.
#' @param malat1_threshold_name name to use for the new meta.data column containing percent IEG gene counts.
#' Default is set dependent on species gene symbol.
#' @param ensembl_ids logical, whether feature names in the object are gene names or
#' ensembl IDs (default is FALSE; set TRUE if feature names are ensembl IDs).
#' @param assay Assay to use (default is the current object default assay).
#' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
#' function will abort if columns with the name provided to `malat1_threshold_name` is present in meta.data slot.
#' @param print_plots logical, should plots be printed to output when running function (default is NULL).
#' Will automatically set to FALSE if performing across samples or TRUE if performing across whole object.
#' @param save_plots logical, whether or not to save plots to pdf (default is FALSE).
#' @param save_plot_path path to save location for plots (default is NULL; current working directory).
#' @param save_plot_name name for pdf file containing plots.
#' @param plot_width the width (in inches) for output page size.  Default is 11.
#' @param plot_height the height (in inches) for output page size.  Default is 8.
#' @param whole_object logical, whether to perform calculation on whole object (default is FALSE).
#' Should be only be run if object contains single sample.
#' @param homolog_name feature name for MALAT1 homolog in non-default species (if annotated).
#' @param bw The "bandwidth" value when plotting the density function to the MALAT1 distribution;
#' default is bw = 0.1, but this parameter should be lowered (e.g. to 0.01) if you run the function and
#' the line that's produced doesn't look like it's tracing the shape of the histogram accurately (this will
#' make the line less "stiff" and more fitted to the data)
#' @param lwd The "line width" fed to the abline function which adds the vertical red line to the output plots;
#' default is 2, and it can be increased or decreased depending on the user's plotting preferences
#' @param breaks The number of bins used for plotting the histogram of normalized MALAT1 values; default is 100
#' @param chosen_min The minimum MALAT1 value cutoff above which a MALAT1 peak in the density function should
#' be found. This value is necessary to determine which peak in the density function fitted to the MALAT1
#' distribution is likely representative of what we would expect to find in real cells. This is because
#' some samples may have large numbers of cells or empty droplets with lower than expected normalized MALAT1 values,
#' and therefore have a peak close to or at zero. Ideally, "chosen_min" would be manually chosen after looking at
#' a histogram of MALAT1 values, and be the normalized MALAT1 value that cuts out all of the cells that look like
#' they stray from the expected distribution (a unimodal distribution above zero). The default value is 1 as this
#' works well in many test cases, but different types of normalization may make the user want to change this
#' parameter (e.g. Seurat's original normalization function generates different results to their SCT function)
#' which may change the MALAT1 distribution). Increase or decrease chosen_min depending on where your MALAT1 peak is located.
#' @param smooth The "smoothing parameter" fed into the "smooth.spline" function that adjusts the trade-off between
#' the smoothness of the line fitting the histogram, and how closely it fits the histogram; the default is 1,
#' and can be lowered if it looks like the line is underfitting the data, and raised in the case of overfitting.
#' The ideal scenario is for the line to trace the histogram in a way where the only inflection point(s) are between
#' major peaks, e.g. separating the group of poor-quality cells or empty droplets with lower normalized MALAT1
#' expression from higher-quality cells with higher normalized MALAT1 expression.
#' @param abs_min The absolute lowest value allowed as the MALAT1 threshold. This parameter increases the
#' robustness of the function if working with an outlier data distribution (e.g. an entire sample is poor
#' quality so there is a unimodal MALAT1 distribution that is very low but above zero, but also many
#' values close to zero) and prevents a resulting MALAT1 threshold of zero. In the case where a calculated
#' MALAT1 value is zero, the function will return 0.3 by default.
#' @param rough_max A rough value for the location of a MALAT1 peak if a peak is not found. This is possible
#' if there are so few cells with higher MALAT1 values, that a distribution fitted to the data finds no local maxima.
#' For example, if a sample only has poor-quality cells such that all have near-zero MALAT1 expression,
#' the fitted function may look similar to a positive quadratic function which has no local maxima.
#' In this case, the function searches for the closest MALAT1 value to the default value, 2, to use in place of
#' a real local maximum.
#'
#' @import cli
#' @import ggplot2
#' @import pbapply
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @method Add_MALAT1_Threshold Seurat
#'
#' @references This function incorporates a threshold calculation and procedure as described in
#' Clarke & Bader (2024). bioRxiv \url{doi.org/10.1101/2024.07.14.603469}.  Please cite this preprint
#' whenever using this function.
#'
#' @author Zoe Clark (original function and manuscript) & Samuel Marsh (wrappers and updates for inclusion in package)
#'
#' @return Seurat object with added meta.data column
#'
#' @export
#' @rdname Add_MALAT1_Threshold
#'
#' @concept qc_util
#'
#' @examples
#' \dontrun{
#' object <- Add_MALAT1_Threshold(object = object, species = "Human")
#' }
#'

Add_MALAT1_Threshold.Seurat <- function(
    object,
    species,
    sample_col = NULL,
    malat1_threshold_name = NULL,
    ensembl_ids = FALSE,
    assay = NULL,
    overwrite = FALSE,
    print_plots = NULL,
    save_plots = FALSE,
    save_plot_path = NULL,
    save_plot_name = NULL,
    plot_width = 11,
    plot_height = 8,
    whole_object = FALSE,
    homolog_name = NULL,
    bw = 0.1,
    lwd = 2,
    breaks = 100,
    chosen_min = 1,
    smooth = 1,
    abs_min = 0.3,
    rough_max = 2
) {
  # Check Seurat
  Is_Seurat(seurat_object = object)

  # Check for sample column
  if (is.null(x = sample_col) && isFALSE(x = whole_object)) {
    cli_abort(message = c("No sample column provided to {.code sample_col} parameter.",
                          "i" = "Please provide name of column in meta.data that contains sample IDs to use for MALAT1 thresholding.")
    )
  }

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

  # Species Spelling Options
  mouse_options <- accepted_names$Mouse_Options
  human_options <- accepted_names$Human_Options
  marmoset_options <- accepted_names$Marmoset_Options
  zebrafish_options <- accepted_names$Zebrafish_Options
  rat_options <- accepted_names$Rat_Options
  drosophila_options <- accepted_names$Drosophila_Options
  macaque_options <- accepted_names$Macaque_Options
  chicken_options <- accepted_names$Chicken_Options

  if (species %in% c(marmoset_options, zebrafish_options, rat_options, drosophila_options, macaque_options, chicken_options)) {
    cli_abort(message = "Rat, Marmoset, Macaque, Zebrafish, Drosophila, and Chicken are not currently supported by default.  Please supply MALAT1 homolog/ortholog feature name.")
  }

  # Set meta data column name
  if (species %in% mouse_options) {
    if (is.null(x = malat1_threshold_name)) {
      malat1_threshold_name <- "Malat1_Threshold"
    }
  }

  if (species %in% human_options) {
    if (is.null(x = malat1_threshold_name)) {
      malat1_threshold_name <- "MALAT1_Threshold"
    }
  }

  # Overwrite check
  if (malat1_threshold_name %in% colnames(x = object@meta.data)) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Column with {.val {malat1_threshold_name}} already present in meta.data slot.",
                            "i" = "*To run function and overwrite column set parameter {.code overwrite = TRUE} or change respective {.code malat1_threshold_name}*")
      )
    }
    cli_inform(message = c("Column with {.val {malat1_threshold_name}} already present in meta.data slot.",
                           "i" = "Overwriting those column as {.code overwrite = TRUE.}")
    )
  }

  # Checks species
  if (is.null(x = species)) {
    cli_abort(message = c("No species name or abbreivation was provided to {.code species} parameter.",
                          "i" = "If not using default species please set {.code species = other}.")
    )
  }

  # check normalized
  Check_Normalized(object = object, assay = assay, error = TRUE)

  if (isTRUE(x = save_plots)) {
    # Set file_path before path check if current dir specified as opposed to leaving set to NULL
    if (!is.null(x = save_plot_path) && save_plot_path == "") {
      save_plot_path <- NULL
    }

    # Check file path is valid
    if (!is.null(x = save_plot_path)) {
      if (!dir.exists(paths = save_plot_path)) {
        cli_abort(message = "Provided {.code save_plot_path}: {symbol$dquote_left}{.field {save_plot_path}}{symbol$dquote_right} does not exist.")
      }
    }

    # Check if file name provided
    if (is.null(x = save_plot_name)) {
      cli_abort(message = "No file name provided.  Please provide a file name using {.code save_plot_name}.")
    }

    # append .pdf if needed
    file_ext <- grep(x = save_plot_name, pattern = ".pdf$")
    if (length(x = file_ext) == 0) {
      cli_abort(message = "{.code save_plot_name} must end with file extension '.pdf'.")
    }

    # check non-default output values
    if (!is.numeric(x = plot_width)) {
      cli_abort(message = "The value provided to {.code plot_width} ({.field {plot_width}}) is not numeric.")
    }

    if (!is.numeric(x = plot_height)) {
      cli_abort(message = "The value provided to {.code plot_height} ({.field {plot_height}}) is not numeric.")
    }
  }

  # Set default assay
  assay <- assay %||% DefaultAssay(object = object)

  # Retrieve gene lists
  if (isFALSE(x = ensembl_ids)) {
    if (species %in% mouse_options) {
      malat_id <- "Malat1"
    }
    if (species %in% human_options) {
      malat_id <- "MALAT1"
    }
  } else {
    malat_id <- Retrieve_MALAT1_Ensembl_Lists(species = species)
  }

  # check malat1 present
  malat_id <- Feature_PreCheck(object = object, features = malat_id)

  # Get data
  cli_inform(message = c("*" = "Adding MALAT1 Threshold for {.field {species}} using gene id: {.val {malat_id}}."))
  cli_inform(message = c("i" = "{col_cyan('Please cite')} {.field Clarke & Bader (2024). doi.org/10.1101/2024.07.14.603469} {col_cyan('when using MALAT1 thresholding function.')}"))

  # split data and run by sample or whole object
  if (isTRUE(x = whole_object)) {
    malat_norm_data <- as.numeric(x = FetchData(object, vars = malat_id, layer = "data")[,1])

    # set plot params
    plot_title <- NULL
    print_plots <- print_plots %||% TRUE

    # run threshold function
    res <- define_malat1_threshold(counts = malat_norm_data, bw = bw, lwd = lwd, breaks = breaks, chosen_min = chosen_min, smooth = smooth, abs_min = abs_min, rough_max = rough_max, print_plots = print_plots, return_plots = TRUE, plot_title = plot_title)

    # split results
    threshold <- res[[1]]

    # save plots
    if (isTRUE(x = save_plots)) {
      plots <- res[[2]]
      ggsave(plots, filename = paste(save_plot_path, save_plot_name, sep=""), width = replace_null(parameter = plot_width), height = replace_null(parameter = plot_height))
    }

    malat1_threshold <- malat_norm_data > threshold
    object[[malat1_threshold_name]] <- malat1_threshold
    object[[malat1_threshold_name]] <- factor(object[[malat1_threshold_name]][,1], levels = c("TRUE","FALSE"))
  } else {
    Idents(object = object) <- sample_col
    cells_by_sample <- CellsByIdentities(object = object)

    sample_col_names <- names(x = cells_by_sample)

    # calculate threshold
    cli_inform(message = "Calculating thresholds across {.field {length(x = sample_col_names)}} samples from meta.data column {.field {sample_col}}.")
    threshold_all <- pblapply(1:length(x = sample_col_names), function(x) {
      malat_norm_data <- as.numeric(x = FetchData(object, vars = malat_id, layer = "data", cells = cells_by_sample[[x]])[,1])

      # run threshold function
      res <- define_malat1_threshold(counts = malat_norm_data, bw = bw, lwd = lwd, breaks = breaks, chosen_min = chosen_min, smooth = smooth, abs_min = abs_min, rough_max = rough_max, print_plots = print_plots, return_plots = TRUE, plot_title = sample_col_names[x])

      threshold <- res[[1]]

      malat1_threshold <- malat_norm_data > threshold
      malat1_threshold <- data.frame(malat1_threshold_name = malat1_threshold)
      rownames(malat1_threshold) <- cells_by_sample[[x]]

      res_list <- list("thresholds" = malat1_threshold,
                       "plots" = res[[2]])
    })

    # save plots
    if (isTRUE(x = save_plots)) {
      cli_inform(message = "{.field Saving plots to file}")

      # Extract plots into a list
      plots_list <- lapply(threshold_all, function(res) res$plots)

      # save plots
      pdf(file = paste(save_plot_path, save_plot_name, sep=""), width = plot_width, height = plot_height)
      pb <- txtProgressBar(min = 0, max = length(x = plots_list), style = 3, file = stderr())
      for (i in 1:length(x = plots_list)) {
        print(plots_list[[i]])
        setTxtProgressBar(pb = pb, value = i)
      }
      close(con = pb)
      dev.off()
    }

    # Combine results and add to object
    cli_inform("Adding results to object as {.val {malat1_threshold_name}}.")
    # Extract thresholds and bind them into a single data frame
    thresholds_list <- lapply(threshold_all, function(res) res$thresholds)
    thresholds_df <- bind_rows(thresholds_list)

    object[[malat1_threshold_name]] <- thresholds_df
    object[[malat1_threshold_name]] <- factor(object[[malat1_threshold_name]][,1], levels = c("TRUE","FALSE"))
  }

  object <- LogSeuratCommand(object = object)

  return(object)
}


#' Add exAM Gene List Module Scores
#'
#' Adds module scores from exAM genes from mouse and human.
#'
#' @param seurat_object object name.
#' @param species Species of origin for given Seurat Object.  Only accepted species are: mouse, human (name or abbreviation).
#' @param exam_module_name name to use for the new meta.data column containing module scores.
#' @param ensembl_ids logical, whether feature names in the object are gene names or
#' ensembl IDs (default is FALSE; set TRUE if feature names are ensembl IDs).
#' @param assay Assay to use (default is the current object default assay).
#' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
#' function will abort if columns with the name provided to `exam_module_name` is present in meta.data slot.
#' @param exclude_unfound logical, whether to exclude features not presne tin current object (default is FALSE).
#' @param seed seed for reproducibility (default is 1).
#'
#' @return Seurat object
#'
#' @import cli
#'
#' @references Gene list is from: SI Table 22 Marsh et al., 2022 (Nature Neuroscience) from \doi{10.1038/s41593-022-01022-8}.
#' See data-raw directory for scripts used to create gene list.
#'
#' @export
#'
#' @concept qc_util
#'
#' @examples
#' \dontrun{
#' # Seurat
#' seurat_object <- exAM_Scoring(seurat_object = seurat_object, species = "human")
#'}
#'

exAM_Scoring <- function(
    seurat_object,
    species,
    exam_module_name = NULL,
    method = "Seurat",
    ensembl_ids = FALSE,
    assay = NULL,
    overwrite = FALSE,
    exclude_unfound = FALSE,
    seed = 1
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

  # Species Spelling Options
  mouse_options <- accepted_names$Mouse_Options
  human_options <- accepted_names$Human_Options
  marmoset_options <- accepted_names$Marmoset_Options
  zebrafish_options <- accepted_names$Zebrafish_Options
  rat_options <- accepted_names$Rat_Options
  drosophila_options <- accepted_names$Drosophila_Options
  macaque_options <- accepted_names$Macaque_Options
  chicken_options <- accepted_names$Chicken_Options

  if (!species %in% unlist(x = accepted_names)) {
    cli_inform(message = "The supplied species ({.field {species}}) is not currently supported.")
  }

  if (species %in% c(marmoset_options, rat_options, zebrafish_options, macaque_options, drosophila_options, chicken_options)) {
    cli_abort(message = c("{.val Rat, Marmoset, Macaque, Zebrafish, Drosophila, Chicken} are not currently supported.",
                         "i" = "No column will be added to object meta.data"))
  }

  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  if (method == "UCell") {
    # Check Nebulosa installed
    UCell_check <- is_installed(pkg = "UCell")
    if (isFALSE(x = UCell_check)) {
      cli_abort(message = c(
        "Please install the {.val UCell} package to use {.code method = {symbol$dquote_left}UCell{symbol$dquote_right}}",
        "i" = "This can be accomplished with the following commands: ",
        "----------------------------------------",
        "{.field `install.packages({symbol$dquote_left}BiocManager{symbol$dquote_right})`}",
        "{.field `BiocManager::install({symbol$dquote_left}UCell{symbol$dquote_right})`}",
        "----------------------------------------"
      ))
    }
  }

  # Overwrite check
  if (exam_module_name %in% colnames(x = seurat_object@meta.data)) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Column with {.val {exam_module_name}} already present in meta.data slot.",
                            "i" = "*To run function and overwrite column set parameter {.code overwrite = TRUE} or change respective {.code exam_module_name}*")
      )
    }
    cli_inform(message = c("Column with {.val {exam_module_name}} already present in meta.data slot.",
                           "i" = "Overwriting those column as {.code overwrite = TRUE.}")
    )
  }

  # Set default assay
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Retrieve gene lists
  if (isFALSE(x = ensembl_ids)) {
    exAM_gene_list <- Retrieve_exAM_Lists(species = species)
  } else {
    exAM_gene_list <- Retrieve_exAM_Ensembl_Lists(species = species)
  }

  if (isTRUE(x = exclude_unfound)) {
    # check features present
    exAM_found <- Feature_PreCheck(object = seurat_object, features = exAM_gene_list[["exAM_union"]])

    if (species %in% human) {
      exAM_found2 <- Feature_PreCheck(object = seurat_object, features = exAM_gene_list[["exAM_micro"]])
    }
  } else {
    # check features present
    exAM_found <- exAM_gene_list[["exAM_union"]]

    if (species %in% human) {
      exAM_found2 <- exAM_gene_list[["exAM_micro"]]
    }
  }

  # set module score names
  if (is.null(x = exam_module_name)) {
    if (species %in% human_options) {
      exam_module_name <- c("exAM_Union_Score", "exAM_Microglia_Score")
    } else {
      exam_module_name <- c("exAM_Union_Score")
    }
  } else {
    if (species %in% human_options && length(x = exam_module_name) != 2) {
      cli_abort("{.code exam_module_name} must be length 2 when {.code species = {symbol$dquote_left}human{symbol$dquote_right}}.")
    }
  }

  # Add mito and ribo columns
  if (length(x = exAM_found) > 0) {
    cli_inform(message = "Adding module score for exAM union gene list as {.field {symbol$dquote_left}{exam_module_name[1]}{symbol$dquote_right}}.")
    seurat_object <- AddModuleScore(object = seurat_object, features = list(exAM_found), name = exam_module_name[1], search = search, seed = seed)
    colnames(seurat_object@meta.data) <- gsub(pattern = paste0(exam_module_name[1], "1"), replacement = exam_module_name[1], x = colnames(seurat_object@meta.data))
  }

  if (species %in% human_options) {
    if (length(x = exAM_found2) > 0) {
      cli_inform(message = "Adding module score for exAM union gene list as {.field {symbol$dquote_left}{exam_module_name[2]}{symbol$dquote_right}}.")
      seurat_object <- AddModuleScore(object = seurat_object, features = list(exAM_found2), name = exam_module_name[2], search = search, seed = seed)
      colnames(seurat_object@meta.data) <- gsub(pattern = paste0(exam_module_name[2], "1"), replacement = exam_module_name[2], x = colnames(seurat_object@meta.data))
    }
  }

  # Log Command
  seurat_object <- LogSeuratCommand(object = seurat_object)

  # return final object
  return(seurat_object)
}


#' Calculate and add differences post-cell bender analysis
#'
#' Calculate the difference in features and UMIs per cell when both cell bender and raw assays are present.
#'
#' @param seurat_object object name.
#' @param raw_assay_name name of the assay containing the raw data.
#' @param cell_bender_assay_name name of the assay containing the Cell Bender'ed data.
#'
#' @importFrom dplyr mutate
#' @importFrom magrittr "%>%"
#'
#' @return Seurat object with 2 new columns in the meta.data slot.
#'
#' @export
#'
#' @concept qc_util
#'
#' @examples
#' \dontrun{
#' object <- Add_CellBender_Diff(seurat_object = obj, raw_assay_name = "RAW",
#' cell_bender_assay_name = "RNA")
#' }
#'

Add_CellBender_Diff <- function(
    seurat_object,
    raw_assay_name,
    cell_bender_assay_name
) {
  # Is Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check assays present
  assays_not_found <- Assay_Present(seurat_object = seurat_object, assay_list = c(raw_assay_name, cell_bender_assay_name), print_msg = FALSE, omit_warn = TRUE)[[2]]

  if (!is.null(x = assays_not_found)) {
    stop_quietly()
  }

  # pull meta
  meta_seurat <- seurat_object@meta.data

  # Set variable names
  raw_nFeature <- paste0("nFeature_", raw_assay_name)
  raw_nCount <- paste0("nCount_", raw_assay_name)

  cb_nFeature <- paste0("nFeature_", cell_bender_assay_name)
  cb_nCount <- paste0("nCount_", cell_bender_assay_name)

  # Mutate meta data
  meta_modified <- meta_seurat %>%
    mutate(nFeature_Diff = .data[[raw_nFeature]] - .data[[cb_nFeature]],
           nCount_Diff = .data[[raw_nCount]] - .data[[cb_nCount]])

  # Remove already in object meta data
  meta_modified <- Meta_Remove_Seurat(meta_data = meta_modified, seurat_object = seurat_object)

  # Add back to Seurat Object
  seurat_object <- AddMetaData(object = seurat_object, metadata = meta_modified)

  # Log Command
  seurat_object <- LogSeuratCommand(object = seurat_object)

  return(seurat_object)
}

