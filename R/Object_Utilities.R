#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GENERAL OBJECT UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Merge a list of Seurat Objects
#'
#' Enables easy merge of a list of Seurat Objects.  See  See \code{\link[SeuratObject]{merge}} for more information,
#'
#' @param list_seurat list composed of multiple Seurat Objects.
#' @param add.cell.ids A character vector of equal length to the number of objects in `list_seurat`.
#' Appends the corresponding values to the start of each objects' cell names.  See \code{\link[SeuratObject]{merge}}.
#' @param merge.data Merge the data slots instead of just merging the counts (which requires renormalization).
#' This is recommended if the same normalization approach was applied to all objects.
#' See \code{\link[SeuratObject]{merge}}.
#' @param project Project name for the Seurat object. See \code{\link[SeuratObject]{merge}}.
#'
#' @import cli
#' @importFrom purrr reduce
#'
#' @return A Seurat Object
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' \dontrun{
#' object_list <- list(obj1, obj2, obj3, ...)
#' merged_object <- Merge_Seurat_List(list_seurat = object_list)
#' }
#'

Merge_Seurat_List <- function(
    list_seurat,
    add.cell.ids = NULL,
    merge.data = TRUE,
    project = "SeuratProject"
) {
  # Check list_seurat is list
  if (!inherits(x = list_seurat, what = "list")) {
    cli_abort(message = "{.code list_seurat} must be environmental variable of class {.val list}")
  }

  # Check list_seurat is only composed of Seurat objects
  for (i in 1:length(x = list_seurat)) {
    if (!inherits(x = list_seurat[[i]], what = "Seurat")) {
      cli_abort("One or more of entries in {.code list_seurat} are not objects of class {.val Seurat}")
    }
  }

  # Check all barcodes are unique to begin with
  duplicated_barcodes <- list_seurat %>%
    lapply(colnames) %>%
    unlist() %>%
    duplicated() %>%
    any()

  if (isTRUE(x = duplicated_barcodes) && is.null(x = add.cell.ids)) {
    cli_abort(message = c("There are overlapping cell barcodes present in the input objects",
                          "i" = "Please rename cells or provide prefixes to {.code add.cell.ids} parameter to make unique.")
    )
  }

  # Check right number of suffix/prefix ids are provided
  if (!is.null(x = add.cell.ids) && length(x = add.cell.ids) != length(x = list_seurat)) {
    cli_abort(message = "The number of prefixes in {.code add.cell.ids} must be equal to the number of objects supplied to {.code list_seurat}.")
  }

  # Rename cells if provided
  list_seurat <- lapply(1:length(x = list_seurat), function(x) {
    list_seurat[[x]] <- RenameCells(object = list_seurat[[x]], add.cell.id = add.cell.ids[x])
  })

  # Merge objects
  merged_object <- reduce(list_seurat, function(x, y) {
    merge(x = x, y = y, merge.data = merge.data, project = project)
  })
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### QC UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Add Mito and Ribo percentages
#'
#' Add Mito, Ribo, & Mito+Ribo percentages to meta.data slot of Seurat Object
#'
#' @param seurat_object object name.
#' @param species Species of origin for given Seurat Object.  If mouse, human, marmoset, zebrafish, rat,
#' drosophila, or rhesus macaque (name or abbreviation) are provided the function will automatically
#' generate mito_pattern and ribo_pattern values.
#' @param mito_name name to use for the new meta.data column containing percent mitochondrial counts.
#' Default is "percent_mito".
#' @param ribo_name name to use for the new meta.data column containing percent ribosomal counts.
#' Default is "percent_ribo".
#' @param mito_ribo_name name to use for the new meta.data column containing percent
#' mitochondrial+ribosomal counts.  Default is "percent_mito_ribo".
#' @param mito_pattern A regex pattern to match features against for mitochondrial genes (will set automatically if
#' species is mouse or human; marmoset features list saved separately).
#' @param ribo_pattern A regex pattern to match features against for ribosomal genes
#' (will set automatically if species is mouse, human, or marmoset).
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
#' contain internal regex/feature lists (human, mouse, marmoset, zebrafish, rat, drosophila, and
#' rhesus macaque).  Default is FALSE.
#'
#' @import cli
#' @importFrom dplyr mutate select intersect all_of
#' @importFrom magrittr "%>%"
#' @importFrom Seurat PercentageFeatureSet AddMetaData
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @return A Seurat Object
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' \dontrun{
#' obj <- Add_Mito_Ribo_Seurat(seurat_object = obj, species = "human")
#'}
#'

Add_Mito_Ribo_Seurat <- function(
  seurat_object,
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
  list_species_names = FALSE
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

  # Return list of accepted default species name options
  if (isTRUE(x = list_species_names)) {
    return(accepted_names)
    stop_quietly()
  }

  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check name collision
  if (any(duplicated(x = c(mito_name, ribo_name, mito_ribo_name)))) {
    cli_abort(message = "One or more of values provided to {.code mito_name}, {.code ribo_name}, {.code mito_ribo_name} are identical.")
  }

  # Overwrite check
  if (mito_name %in% colnames(x = seurat_object@meta.data) || ribo_name %in% colnames(x = seurat_object@meta.data) || mito_ribo_name %in% colnames(x = seurat_object@meta.data)) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Columns with {.val {mito_name}} and/or {.val {ribo_name}} already present in meta.data slot.",
                            "i" = "*To run function and overwrite columns set parameter {.code overwrite = TRUE} or change respective {.code mito_name}, {.code ribo_name}, and/or {.code mito_ribo_name}*")
      )
    }
    cli_inform(message = c("Columns with {.val {mito_name}} and/or {.val {ribo_name}} already present in meta.data slot.",
                           "i" = "Overwriting those columns as .code {overwrite = TRUE.}")
    )
  }

  # Checks species
  if (is.null(x = species)) {
    cli_abort(message = c("No species name or abbreivation was provided to {.code species} parameter.",
                          "i" = "If not using default species please set {.code species = other}.")
    )
  }

  # Set default assay
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Species Spelling Options
  mouse_options <- accepted_names$Mouse_Options
  human_options <- accepted_names$Human_Options
  marmoset_options <- accepted_names$Marmoset_Options
  zebrafish_options <- accepted_names$Zebrafish_Options
  rat_options <- accepted_names$Rat_Options
  drosophila_options <- accepted_names$Drosophila_Options
  macaque_options <- accepted_names$Macaque_Options

  # Check ensembl vs patterns
  if (isTRUE(x = ensembl_ids) && species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options) && any(!is.null(x = mito_pattern), !is.null(x = ribo_pattern), !is.null(x = mito_features), !is.null(x = ribo_features))) {
    cli_warn(message = c("When using a default species and setting {.code ensembl_ids = TRUE} provided patterns or features are ignored.",
                         "*" = "Supplied {.code mito_pattern}, {.code ribo_pattern}, {.code mito_features}, {.code ribo_features} will be disregarded.")
    )
  }

  # Assign mito/ribo pattern to stored species
  if (species %in% c(mouse_options, human_options, marmoset_options, zebrafish_options, rat_options, drosophila_options) && any(!is.null(x = mito_pattern), !is.null(x = ribo_pattern))) {
    cli_warn(message = c("Pattern expressions for included species are set by default.",
                         "*" = "Supplied {.code mito_pattern} and {.code ribo_pattern} will be disregarded.",
                         "i" = "To override defaults please supply a feature list for mito and/or ribo genes.")
    )
  }

  if (species %in% mouse_options) {
    mito_pattern <- "^mt-"
    ribo_pattern <- "^Rp[sl]"
  }
  if (species %in% human_options) {
    mito_pattern <- "^MT-"
    ribo_pattern <- "^RP[SL]"
  }
  if (species %in% c(marmoset_options, macaque_options)) {
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

  mito_features <- mito_features %||% grep(pattern = mito_pattern, x = rownames(x = seurat_object[[assay]]), value = TRUE)

  ribo_features <- ribo_features %||% grep(pattern = ribo_pattern, x = rownames(x = seurat_object[[assay]]), value = TRUE)

  # Check features are present in object
  length_mito_features <- length(x = intersect(x = mito_features, y = rownames(x = seurat_object[[assay]])))

  length_ribo_features <- length(x = intersect(x = ribo_features, y = rownames(x = seurat_object[[assay]])))

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

  # Add mito and ribo columns
  if (length_mito_features > 0) {
    good_mito <- mito_features[mito_features %in% rownames(x = seurat_object)]
    seurat_object[[mito_name]] <- PercentageFeatureSet(object = seurat_object, features = good_mito, assay = assay)
  }
  if (length_ribo_features > 0) {
    good_ribo <- ribo_features[ribo_features %in% rownames(x = seurat_object)]
    seurat_object[[ribo_name]] <- PercentageFeatureSet(object = seurat_object, features = good_ribo, assay = assay)
  }

  # Create combined mito ribo column if both present
  if (length_mito_features > 0 && length_ribo_features > 0) {
    object_meta <- Fetch_Meta(object = seurat_object) %>%
      rownames_to_column("barcodes")

    object_meta <- object_meta %>%
      mutate({{mito_ribo_name}} := .data[[mito_name]] + .data[[ribo_name]])

    object_meta <- object_meta %>%
      select(all_of(c("barcodes", mito_ribo_name))) %>%
      column_to_rownames("barcodes")

    seurat_object <- AddMetaData(object = seurat_object, metadata = object_meta)
  }

  # return final object
  return(seurat_object)
}


#' Add Cell Complexity Value
#'
#' Add measure of cell complexity/novelty (log10PerUMI) for data QC.
#'
#' @param seurat_object object name.
#' @param meta_col_name name to use for new meta data column.  Default is "log10GenesPerUMI".
#' @param assay assay to use in calculation.  Default is "RNA".  *Note* This should only be changed if
#' storing corrected and uncorrected assays in same object (e.g. outputs of both Cell Ranger and Cell Bender).
#' @param overwrite Logical.  Whether to overwrite existing an meta.data column.  Default is FALSE meaning that
#' function will abort if column with name provided to `meta_col_name` is present in meta.data slot.
#'
#' @import cli
#'
#' @return A Seurat Object
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' library(Seurat)
#' pbmc_small <- Add_Cell_Complexity_Seurat(seurat_object = pbmc_small)
#'

Add_Cell_Complexity_Seurat <- function(
  seurat_object,
  meta_col_name = "log10GenesPerUMI",
  assay = "RNA",
  overwrite = FALSE
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Add assay warning message
  if (assay != "RNA") {
    cli_warn(message = "Assay is set to value other than 'RNA'. This should only be done in rare instances.  See documentation for more info ({.code ?Add_Cell_Complexity_Seurat}).",
             .frequency = "once",
             .frequency_id = "assay_warn")
  }

  # Check columns for overwrite
  if (meta_col_name %in% colnames(x = seurat_object@meta.data)) {
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
  seurat_object[[meta_col_name]] <- log10(seurat_object[[feature_name]]) / log10(seurat_object[[count_name]])

  #return object
  return(seurat_object)
}


#' Add Percent of High Abundance Genes
#'
#' Add the percentage of counts occupied by the top XX most highly expressed genes in each cell.
#'
#' @param seurat_object object name.
#' @param num_top_genes An integer vector specifying the size(s) of the top set of high-abundance genes.
#' Used to compute the percentage of library size occupied by the most highly expressed genes in each cell.
#' @param meta_col_name name to use for new meta data column.  Default is "percent_topXX", where XX is
#' equal to the value provided to `num_top_genes`.
#' @param assay assay to use in calculation.  Default is "RNA".  *Note* This should only be changed if
#' storing corrected and uncorrected assays in same object (e.g. outputs of both Cell Ranger and Cell Bender).
#' @param overwrite Logical.  Whether to overwrite existing an meta.data column.  Default is FALSE meaning that
#' function will abort if column with name provided to `meta_col_name` is present in meta.data slot.
#'
#' @import cli
#' @importFrom dplyr select all_of
#' @importFrom magrittr "%>%"
#' @importFrom rlang is_installed
#'
#' @return A Seurat Object
#'
#' @export
#'
#' @concept object_util
#'
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
#' library(Seurat)
#' pbmc_small <- Add_Top_Gene_Pct_Seurat(seurat_object = pbmc_small, num_top_genes = 50)
#'

Add_Top_Gene_Pct_Seurat <- function(
    seurat_object,
    num_top_genes = 50,
    meta_col_name = NULL,
    assay = "RNA",
    overwrite = FALSE
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
  Is_Seurat(seurat_object = seurat_object)

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
  if (meta_col_name %in% colnames(x = seurat_object@meta.data)) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Column {.val {meta_col_name}} already present in meta.data slot.",
                            "i" = "*To run function and overwrite column, set parameter {.code overwrite = TRUE} or change respective {.code meta_col_name}*.")
      )
    }
    cli_inform(message = c("Column {.val {meta_col_name}} already present in meta.data slot",
                           "i" = "Overwriting those columns as {.code overwrite = TRUE}.")
    )
  }

  # Extract matrix
  count_mat <- LayerData(object = seurat_object, assay = assay)

  # calculate
  res <- as.data.frame(scuttle::perCellQCMetrics(x = count_mat, percent.top = num_top_genes))

  # select percent column
  res <- res %>%
    select(all_of(scuttle_colname))

  # Add to object and return
  seurat_object <- AddMetaData(object = seurat_object, metadata = res, col.name = meta_col_name)

  return(seurat_object)
}


#' Add Multiple Cell Quality Control Values with Single Function
#'
#' Add Mito/Ribo %, Cell Complexity (log10GenesPerUMI), Top Gene Percent with single function call
#'
#' @param seurat_object object name.
#' @param add_mito_ribo logical, whether to add percentage of counts belonging to mitochondrial/ribosomal
#' genes to object (Default is TRUE).
#' @param add_complexity logical, whether to add Cell Complexity to object (Default is TRUE).
#' @param add_top_pct logical, whether to add Top Gene Percentages to object (Default is TRUE).
#' @param add_MSigDB logical, whether to add percentages of counts belonging to genes from of mSigDB hallmark
#' gene lists: "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_APOPTOSIS", and "HALLMARK_DNA_REPAIR" to
#' object (Default is TRUE).
#' @param add_IEG logical, whether to add percentage of counts belonging to IEG genes to object (Default is TRUE).
#' @param add_cell_cycle logical, whether to addcell cycle scores and phase based on
#' \code{\link[Seurat]{CellCycleScoring}}.  Only applicable if `species = "human"`.  (Default is TRUE).
#' @param species Species of origin for given Seurat Object.  If mouse, human, marmoset, zebrafish, rat,
#' drosophila, or rhesus macaque (name or abbreviation) are provided the function will automatically
#' generate mito_pattern and ribo_pattern values.
#' @param mito_name name to use for the new meta.data column containing percent mitochondrial counts.
#' Default is "percent_mito".
#' @param ribo_name name to use for the new meta.data column containing percent ribosomal counts.
#' Default is "percent_ribo".
#' @param mito_ribo_name name to use for the new meta.data column containing percent
#' mitochondrial+ribosomal counts.  Default is "percent_mito_ribo".
#' @param complexity_name name to use for new meta data column for `Add_Cell_Complexity_Seurat`.
#' Default is "log10GenesPerUMI".
#' @param top_pct_name name to use for new meta data column for `Add_Top_Gene_Pct_Seurat`.
#' Default is "percent_topXX", where XX is equal to the value provided to `num_top_genes`.
#' @param oxphos_name name to use for new meta data column for percentage of MSigDB oxidative phosphorylation
#' counts.  Default is "percent_oxphos".
#' @param apop_name name to use for new meta data column for percentage of MSigDB apoptosis counts.
#' Default is "percent_apop".
#' @param dna_repair_name name to use for new meta data column for percentage of MSigDB DNA repair
#' counts.  Default is "percent_dna_repair"..
#' @param ieg_name name to use for new meta data column for percentage of IEG counts.  Default is "percent_ieg".
#' @param mito_pattern A regex pattern to match features against for mitochondrial genes (will set automatically if
#' species is mouse or human; marmoset features list saved separately).
#' @param ribo_pattern A regex pattern to match features against for ribosomal genes
#' (will set automatically if species is mouse, human, or marmoset).
#' @param mito_features A list of mitochondrial gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param ribo_features A list of ribosomal gene names to be used instead of using regex pattern.
#' Will override regex pattern if both are present (including default saved regex patterns).
#' @param ensembl_ids logical, whether feature names in the object are gene names or
#' ensembl IDs (default is FALSE; set TRUE if feature names are ensembl IDs).
#' @param num_top_genes An integer vector specifying the size(s) of the top set of high-abundance genes.
#' Used to compute the percentage of library size occupied by the most highly expressed genes in each cell.
#' @param assay assay to use in calculation.  Default is "RNA".  *Note* This should only be changed if
#' storing corrected and uncorrected assays in same object (e.g. outputs of both Cell Ranger and Cell Bender).
#' @param overwrite Logical.  Whether to overwrite existing an meta.data column.  Default is FALSE meaning that
#' function will abort if column with name provided to `meta_col_name` is present in meta.data slot.
#'
#' @import cli
#'
#' @return A Seurat Object
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' \dontrun{
#' obj <- Add_Cell_QC_Metrics(seurat_object = obj, species = "Human")
#'}
#'

Add_Cell_QC_Metrics <- function(
    seurat_object,
    add_mito_ribo = TRUE,
    add_complexity = TRUE,
    add_top_pct = TRUE,
    add_MSigDB = TRUE,
    add_IEG = TRUE,
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
    mito_pattern = NULL,
    ribo_pattern = NULL,
    mito_features = NULL,
    ribo_features = NULL,
    ensembl_ids = FALSE,
    num_top_genes = 50,
    assay = NULL,
    overwrite = FALSE
) {
  # Set assay
  assay <- assay %||% DefaultAssay(object = seurat_object)

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
    cli_inform(message = "Adding {.field Mito/Ribo Percentages} to meta.data.")
    seurat_object <- Add_Mito_Ribo_Seurat(seurat_object = seurat_object, species = species, mito_name = mito_name, ribo_name = ribo_name, mito_ribo_name = mito_ribo_name, mito_pattern = mito_pattern, ribo_pattern = ribo_pattern, mito_features = mito_features, ribo_features = ribo_features, ensembl_ids = ensembl_ids, assay = assay, overwrite = overwrite)
  }

  # Add complexity
  if (isTRUE(x = add_complexity)) {
    cli_inform(message = "Adding {.field Cell Complexity #1 (log10GenesPerUMI)} to meta.data.")
    seurat_object <- Add_Cell_Complexity_Seurat(seurat_object = seurat_object, meta_col_name = complexity_name, assay = assay, overwrite = overwrite)
  }

  # Add top gene expression percent
  if (isTRUE(x = add_top_pct)) {
    cli_inform(message = "Adding {.field Cell Complexity #2 (Top {num_top_genes} Percentages)} to meta.data.")
    seurat_object <- Add_Top_Gene_Pct_Seurat(seurat_object = seurat_object, num_top_genes = num_top_genes, meta_col_name = top_pct_name, assay = assay, overwrite = overwrite)
  }

  # Add MSigDB
  if (isTRUE(x = add_MSigDB)) {
    if (species %in% marmoset_options) {
      cli_warn(message = c("{.val Marmoset} is not currently a part of MSigDB gene list database.",
                           "i" = "No columns will be added to object meta.data"))
    } else {
      cli_inform(message = "Adding {.field MSigDB Oxidative Phosphorylation, Apoptosis, and DNA Repair Percentages.} to meta.data.")
      seurat_object <- Add_MSigDB_Seurat(seurat_object = seurat_object, species = species, oxphos_name = oxphos_name, apop_name = apop_name, dna_repair_name = dna_repair_name, assay = assay, overwrite = overwrite)
    }
  }

  # Add IEG
  if (isTRUE(x = add_IEG)) {
    if (species %in% c(marmoset_options, rat_options, zebrafish_options, macaque_options, drosophila_options)) {
      cli_warn(message = c("{.val Rat, Marmoset, Macaque, Zebrafish, and Drosophila} are not currently supported.",
                           "i" = "No column will be added to object meta.data"))
    } else {
      cli_inform(message = "Adding {.field MSigDB Oxidative Phosphorylation, Apoptosis, and DNA Repair Percentages.} to meta.data.")
      seurat_object <- Add_IEG_Seurat(seurat_object = seurat_object, species = species, ieg_name = ieg_name, assay = assay, overwrite = overwrite)
    }
  }

  if (isTRUE(x = add_cell_cycle)) {
    if (!species %in% human_options) {
      cli_abort(message = c("Cell Cycle Scoring is only supported for human in this function.",
                            "i" = "To add score for other species supply cell cycle gene list of `CellCycleScoring` function."
                ))
    } else {
      if (length(grep(x = Layers(object = seurat_object), pattern = "data", value = T)) == 0) {
        cli_inform(message = c("Layer with normalized data not present.",
                               "i" = "Normalizing Data."))
        seurat_object <- NormalizeData(object = seurat_object)
      }

      # Overwrite check
      if ("S.Score" %in% colnames(x = seurat_object@meta.data) || "G2M.Score" %in% colnames(x = seurat_object@meta.data) || "Phase" %in% colnames(x = seurat_object@meta.data)) {
        if (!overwrite) {
          cli_abort(message = c("Columns with {.val S.Score}, {.val G2M.Score} and/or {.val Phase} already present in meta.data slot.",
                                "i" = "*To run function and overwrite columns set parameter {.code overwrite = TRUE}*")
          )
        }
        cli_inform(message = c("Columns with {.val S.Score}, {.val G2M.Score} and/or {.val Phase} already present in meta.data slot.",
                               "i" = "Overwriting those columns as .code {overwrite = TRUE.}")
        )
      }

      # Add Cell Cycle Scoring
      cli_inform(message = "Adding {.field Cell Cycle Scoring to meta.data.} to meta.data.")
      seurat_object <- CellCycleScoring(object = seurat_object, s.features = Seurat::cc.genes.updated.2019$s.genes, g2m.features = Seurat::cc.genes.updated.2019$g2m.genes)
    }
  }

  # return object
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
#' @concept object_util
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

  return(seurat_object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### META DATA UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Remove meta data columns containing Seurat Defaults
#'
#' Remove any columns from new meta_data data.frame in preparation for adding back to Seurat Object
#'
#' @param meta_data data.frame containing meta data.
#' @param seurat_object object name.
#' @param barcodes_to_rownames logical, are barcodes present as column and should they be moved to
#' rownames (to be compatible with `Seurat::AddMetaData`).  Default is FALSE.
#' @param barcodes_colname name of barcodes column in meta_data.  Required if `barcodes_to_rownames = TRUE`.
#'
#' @import cli
#' @importFrom dplyr select all_of
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#'
#' @return data.frame with only new columns.
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' \dontrun{
#' new_meta <- Meta_Remove_Seurat(meta_data = meta_data_df, seurat_object = object)
#' object <- AddMetaData(object = object, metadata = new_meta)
#' }
#'

Meta_Remove_Seurat <- function(
  meta_data,
  seurat_object,
  barcodes_to_rownames = FALSE,
  barcodes_colname = "barcodes"
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Get existing meta.data column names
  existing_names <- colnames(x = seurat_object@meta.data)

  # Filter meta_data to remove already existing columns
  meta_data_filtered <- meta_data %>%
    select(-all_of(x = existing_names))

  if (isTRUE(x = barcodes_to_rownames)) {
    # Check barcodes colname exists
    if (!barcodes_colname %in% colnames(x = meta_data)) {
      cli_abort(message = "{.code barcodes_colname}: {.val {barcodes_colname}} was not present in the column names of meta_data data.frame provided.")
    }
    # Move barcodes to rownames
    meta_data_filtered <- meta_data_filtered %>%
      column_to_rownames(barcodes_colname)
    # return data.frame
    return(meta_data_filtered)
  }

  # Return filtered data.frame
  return(meta_data_filtered)
}


#' Add Sample Level Meta Data
#'
#' Add meta data from ample level data.frame/tibble to cell level seurat `@meta.data` slot
#'
#' @param seurat_object object name.
#' @param meta_data data.frame/tibble containing meta data or path to file to read.  Must be formatted as
#' either data.frame or tibble.
#' @param join_by_seurat name of the column in `seurat_object@meta.data` that contains matching
#' variables to `join_by_meta` in `meta_data.`
#' @param join_by_meta name of the column in `meta_data` that contains matching
#' variables to `join_by_seurat` in `seurat_object@meta.data`.
#' @param na_ok logical, is it ok to add NA values to `seurat_object@meta.data`.  Default is FALSE.
#' Be very careful if setting TRUE because if there is error in join operation it may result in all
#' `@meta.data` values being replaced with NA.
#' @param overwrite logical, if there are shared columns between `seurat_object@meta.data` and `meta_data`
#' should the current `seurat_object@meta.data` columns be overwritten.  Default is FALSE.  This parameter
#' excludes values provided to `join_by_seurat` and `join_by_meta`.
#'
#' @import cli
#' @importFrom data.table fread
#' @importFrom dplyr select left_join all_of
#' @importFrom magrittr "%>%"
#' @importFrom stats setNames
#' @importFrom tibble column_to_rownames rownames_to_column
#'
#' @return Seurat object with new `@meta.data` columns
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' \dontrun{
#' # meta_data present in environment
#' sample_level_meta <- data.frame(...)
#' obj <- Add_Sample_Meta(seurat_object = obj, meta_data = sample_level_meta,
#' join_by_seurat = "orig.ident", join_by_meta = "sample_ID")
#'
#' # from meta data file
#' obj <- Add_Sample_Meta(seurat_object = obj, meta_data = "meta_data/sample_level_meta.csv",
#' join_by_seurat = "orig.ident", join_by_meta = "sample_ID")
#' }
#'

Add_Sample_Meta <- function(
  seurat_object,
  meta_data,
  join_by_seurat,
  join_by_meta,
  na_ok = FALSE,
  overwrite = FALSE
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check variable vs. file path
  if (isTRUE(x = exists(x = deparse(expr = substitute(expr = meta_data))))) {
    meta_data <- meta_data
  } else {
    if (file.exists(meta_data)) {
      meta_data <- fread(file = meta_data, data.table = FALSE)
    } else {
      cli_abort(message = c("Could not find {.code meta_data} {.field {deparse(expr = substitute(expr = meta_data))}}.",
                            "*" = "If providing environmental variable please check {.code meta_data} name.",
                            "i" = "If providing path to file please check path is correct.")
      )
    }
  }

  # Check meta data structure
  if (!inherits(x = meta_data, what = "data.frame")) {
    cli_abort(message = c("{.code meta_data} not in correct format",
                          "*" = "{.code meta_data} must be a data.frame or tibble.",
                          "i" = "Change format and re-run function.")
    )
  }

  # Check NA in meta data
  if (anyNA(x = meta_data)) {
    cli_abort(message = c("{.code meta_data} contains NA values.",
                          "i" = "If you would like NA values added to Seurat meta data please set {.code na_ok = TRUE}.")
    )
  }

  # Check join variables exist
  if (!join_by_seurat %in% colnames(x = seurat_object@meta.data)) {
    cli_abort(message = "The column {.val {join_by_seurat}} was not found in object @meta.data slot."
    )
  }

  if (!join_by_meta %in% colnames(x = meta_data)) {
    cli_abort(message = "The column {.val {join_by_meta}} was not found in supplied `meta_data`."
    )
  }

  # Check if any duplicate column names
  dup_columns <- colnames(x = meta_data)[colnames(x = meta_data) %in% colnames(x = seurat_object@meta.data)]

  if (length(x = dup_columns) > 0) {
    dup_columns <- dup_columns[!dup_columns %in% c(join_by_seurat, join_by_meta)]

    if (any(dup_columns %in% colnames(x = seurat_object@meta.data)) && !overwrite) {
      cli_abort(message = c(" Duplicate {.code meta_data} contains column names in object @meta.data.",
                            "i" = "{.code meta_data} and object@meta.data both contain columns: {.field {glue_collapse_scCustom(input_string = dup_columns)}}.",
                            "*" = "To overwrite existing object @meta.data columns with those in {.code meta_data} set {.code overwrite = TRUE}.")
      )
    }
  }

  # Pull meta data
  meta_seurat <- seurat_object@meta.data %>%
    rownames_to_column("barcodes")

  # remove
  if (isTRUE(x = overwrite)) {
    meta_seurat <- meta_seurat %>%
      select(-all_of(x = dup_columns))
  } else {
    meta_seurat <- meta_seurat
  }

  # join
  meta_merged <- left_join(x = meta_seurat, y = meta_data, by = setNames(join_by_meta, join_by_seurat))

  # Remove existing Seurat meta
  if (length(x = dup_columns) > 0 && isTRUE(x = overwrite)) {
    meta_merged <- meta_merged %>%
      column_to_rownames("barcodes")
  } else {
    meta_merged <- Meta_Remove_Seurat(meta_data = meta_merged, seurat_object = seurat_object) %>%
      column_to_rownames("barcodes")
  }

  # check NA
  if (anyNA(x = meta_merged) && !na_ok) {
    cli_abort(message = c("NAs found in new seurat meta.data.",
                          "*" = "No new meta data added to Seurat object.",
                          "i" = "Check to make sure all levels of joining factors are present in both sets of meta data.")
    )
  }

  # Add meta data
  seurat_object <- AddMetaData(object = seurat_object, metadata = meta_merged)

  # Return object
  return(seurat_object)
}


#' Extract sample level meta.data
#'
#' Returns a by identity meta.data data.frame with one row per sample.  Useful for downstream
#' quick view of sample breakdown, meta data table creation, and/or use in pseudobulk analysis
#'
#' @param object Seurat object
#' @param sample_name meta.data column to use as sample.  Output data.frame will contain one row per
#' level or unique value in this variable.
#' @param variables_include `@meta.data` columns to keep in final data.frame.  All other columns will
#' be discarded.  Default is NULL.
#' @param variables_exclude columns to discard in final data.frame.  Many cell level columns are
#' irrelevant at the sample level (e.g., nFeature_RNA, percent_mito).
#' \itemize{
#' \item Default parameter value is `NULL` but internally will set to discard nFeature_ASSAY(s),
#' nCount_ASSAY(s), percent_mito, percent_ribo, percent_mito_ribo, and log10GenesPerUMI.
#' \item If sample level median values are desired for these type of variables the output of this
#' function can be joined with output of \code{\link[scCustomize]{Median_Stats}}.
#' \item Set parameter to `include_all = TRUE` to prevent any columns from being excluded.
#' }
#' @param include_all logical, whether or not to include all object meta data columns in output data.frame.
#' Default is FALSE.
#'
#' @return Returns a data.frame with one row per `sample_name`.
#'
#' @import cli
#' @importFrom dplyr any_of grouped_df select slice
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' library(Seurat)
#' pbmc_small$batch <- sample(c("batch1", "batch2"), size = ncol(pbmc_small), replace = TRUE)
#'
#' sample_meta <- Extract_Sample_Meta(object = pbmc_small, sample_name = "orig.ident")
#'
#' # Only return specific columns from meta data (orig.ident and batch)
#' sample_meta <- Extract_Sample_Meta(object = pbmc_small, sample_name = "orig.ident",
#' variables_include = "batch")
#'
#' # Return all columns from meta data
#' sample_meta <- Extract_Sample_Meta(object = pbmc_small, sample_name = "orig.ident",
#' include_all = TRUE)
#'

Extract_Sample_Meta <- function(
  object,
  sample_name = "orig.ident",
  variables_include = NULL,
  variables_exclude = NULL,
  include_all = FALSE
) {
  # Pull meta data
  meta_df <- Fetch_Meta(object = object)

  # Check sample name parameter is present
  if (!sample_name %in% colnames(x = meta_df)) {
    cli_abort(message = "The {.code sample_name} parameter: {.val {sample_name}} was not found in object meta.data")
  }

  if (!is.null(x = variables_include) && !is.null(x = variables_exclude) && include_all) {
    cli_abort(message = "If {.code include_all = TRUE} then {.code variables_include} and {.code variables_exclude} must be NULL.")
  }

  # Generate nCount and nFeature variable vectors for exclusion
  if (is.null(x = variables_exclude)) {
    nFeature_cols <- grep(x = colnames(x = object@meta.data), pattern = "^nFeature", value = TRUE)

    nCount_cols <- grep(x = colnames(x = object@meta.data), pattern = "^nCount", value = TRUE)

    combined_exclude <- c(nFeature_cols, nCount_cols, "percent_mito", "percent_ribo", "percent_mito_ribo", "log10GenesPerUMI")

    variables_exclude <- Meta_Present(seurat_object = object, meta_col_names = combined_exclude, omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)[[1]]
  }

  # Ensure include exclude are unique
  if (!is.null(x = variables_include) && !is.null(x = variables_exclude)) {
    variable_intersect <- intersect(x = variables_include, y = variables_exclude)
    if (length(x = variable_intersect > 0)) {
      cli_abort(message = c("Following variable is present in both {.code variables_include} and {.code variables_exclude} parameters.",
                            "i" = "{.field {variable_intersect}}")
      )
    }
  }

  # Check variables include/exclude are present
  if (!is.null(x = variables_include)) {
    include_meta_list <- Meta_Present(seurat_object = object, meta_col_names = variables_include, omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)
  } else {
    include_meta_list <- NULL
  }

  if (!is.null(x = variables_exclude)) {
    exclude_meta_list <- Meta_Present(seurat_object = object, meta_col_names = variables_exclude, omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)
  } else {
    exclude_meta_list <- NULL
  }

  all_not_found_features <- c(include_meta_list[[2]], exclude_meta_list[[2]])

  all_found_features <- c(include_meta_list[[1]], exclude_meta_list[[1]])

  # Stop if no features found
  if (length(x = all_found_features) < 1) {
    cli_abort(message = c("No features were found.",
                          "*" = "The following are not present in object meta data:",
                          "i" = "{.field {glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}}")
    )
  }

  # Return message of features not found
  if (length(x = all_not_found_features) > 0) {
    op <- options(warn = 1)
    on.exit(options(op))
    cli_warn(message = c("The following meta data variables were omitted as they were not found:",
                         "i" = "{.field {glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}}")
    )
  }

  # Create by sample data.frame
  sample_meta_df <- meta_df %>%
    grouped_df(vars = sample_name) %>%
    slice(1)

  # remove rownames
  rownames(x = sample_meta_df) <- NULL

  # Filter data.frame
  if (isTRUE(x = include_all)) {
    sample_meta_df_filtered <- sample_meta_df
  } else {
    if (length(x = include_meta_list[[1]]) > 0) {
      sample_meta_df_filtered <- sample_meta_df %>%
        select(any_of(c(include_meta_list[[1]], sample_name)))
      if (length(x = exclude_meta_list[[1]]) > 0) {
        sample_meta_df_filtered <- sample_meta_df_filtered %>%
          select(-any_of(exclude_meta_list[[1]]))
      }
    } else {
      if (length(x = exclude_meta_list[[1]]) > 0) {
        sample_meta_df_filtered <- sample_meta_df %>%
          select(-any_of(exclude_meta_list[[1]]))
      }
    }
  }

  # return data
  return(sample_meta_df_filtered)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### MISC SLOT UTILITIES ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Store misc data in Seurat object
#'
#' Wrapper function save variety of data types to the `object@misc` slot of Seurat object.
#'
#' @param seurat_object object name.
#' @param data_to_store data to be stored in `@misc` slot.  Can be single piece of data or list.
#' If list of data see `list_as_list` parameter for control over data storage.
#' @param data_name name to give the entry in `@misc` slot.  Must be of equal length of the number
#' of data items being stored.
#' @param list_as_list logical.  If `data_to_store` is a list, this dictates whether to store in `@misc` slot
#' as list (TRUE) or whether to store each entry in the list separately (FALSE).  Default is FALSE.
#' @param overwrite Logical.  Whether to overwrite existing items with the same name.  Default is FALSE, meaning
#' that function will abort if item with `data_name` is present in misc slot.
#'
#' @return Seurat Object with new entries in the `@misc` slot.
#'
#' @import cli
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' library(Seurat)
#' clu_pal <- c("red", "green", "blue")
#'
#' pbmc_small <- Store_Misc_Info_Seurat(seurat_object = pbmc_small, data_to_store = clu_pal,
#' data_name = "rd1_colors")
#'

Store_Misc_Info_Seurat <- function(
  seurat_object,
  data_to_store,
  data_name,
  list_as_list = FALSE,
  overwrite = FALSE
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Check if name already present
  misc_present <- names(x = seurat_object@misc)
  if (data_name %in% misc_present) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Item(s) named: {.val {data_name}} already present in @misc slot.",
                            "i" = "*To run function and overwrite items set parameter {.code overwrite = TRUE} or change {.code data_name}*")
      )
    } else {
      cli_inform(message = c("Items named {.val {data_name}} already present in @misc slot.",
                             "i" = "Overwriting those items as overwrite = TRUE.")
      )
    }
  }

  # Commenting our for now.  Not sure why you would need this check...
  # # Check length of data
  # if (length(x = data_to_store) != 1 && class(x = data_to_store) != "list") {
  #   stop("'data_to_store' must be either single vector/string or a single list of items.")
  # }

  # Check class of data
  if (inherits(x = data_to_store, what = "list")) {
    if (isTRUE(x = list_as_list)) {
      # Check length of name
      if (length(x = data_name) != 1) {
        cli_abort(message = "When storing as list the length {.code data_name} must be {.field 1 (one)}.")
      }

      # Add data
      seurat_object@misc[[data_name]] <- data_to_store
      cli_inform(message = c("Seurat Object now contains the following items in @misc slot: ",
                             "i" = "{.field {paste(shQuote(names(x = seurat_object@misc)), collapse=", ")}}")
      )
      return(seurat_object)
    }

    # length of list
    data_list_length <- length(x = data_to_store)

    if (length(x = data_name) != data_list_length) {
      cli_abort(message = "The lengths of {.code data_to_store} ({.field {data_list_length}}) and {.code data_name} ({.field {length(x = data_name)}}) must be equal.")
    }

    # Add data
    for (i in 1:data_list_length) {
      seurat_object@misc[[data_name[i]]] <- data_to_store[[i]]
    }
    cli_inform(message = c("Seurat Object now contains the following items in @misc slot: ",
                           "i" = "{.field {paste(shQuote(names(x = seurat_object@misc)), collapse=", ")}}")
    )
    return(seurat_object)
  } else {
    # Check length of name
    if (length(x = data_name) != 1) {
      cli_abort(message = "When storing a string/vector the length {.code data_name} must be {.field 1 (one)}.")
    }

    # Add data
    seurat_object@misc[[data_name]] <- data_to_store
    misc_names <- shQuote(string = names(x = seurat_object@misc))
    cli_inform(message = c("Seurat Object now contains the following items in @misc slot: ",
                           "i" = "{.field {glue_collapse_scCustom(input_string = misc_names, and = TRUE)}}")
    )
    return(seurat_object)
  }
}


#' Store color palette in Seurat object
#'
#' Wrapper function around `Store_Misc_Info_Seurat` to store color palettes.
#'
#' @param seurat_object object name.
#' @param palette vector or list of vectors containing color palettes to store.  If list of palettes
#' see `list_as_list` parameter for control over data storage.
#' @param palette_name name to give the palette(s) in `@misc` slot.  Must be of equal length to the number
#' of data items being stored.
#' @param list_as_list logical.  If `data_to_store` is a list, this dictates whether to store in `@misc` slot
#' as list (TRUE) or whether to store each entry in the list separately (FALSE).  Default is FALSE.
#' @param overwrite Logical.  Whether to overwrite existing items with the same name.  Default is FALSE, meaning
#' that function will abort if item with `data_name` is present in misc slot.
#'
#' @return Seurat Object with new entries in the `@misc` slot.
#'
#' @export
#'
#' @concept object_util
#'
#' @examples
#' library(Seurat)
#' clu_pal <- c("red", "green", "blue")
#'
#' pbmc_small <- Store_Misc_Info_Seurat(seurat_object = pbmc_small, data_to_store = clu_pal,
#' data_name = "rd1_colors")
#'
#'

Store_Palette_Seurat <- function(
  seurat_object,
  palette,
  palette_name,
  list_as_list = FALSE,
  overwrite = FALSE
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  seurat_object <- Store_Misc_Info_Seurat(seurat_object = seurat_object, data_to_store = palette, data_name = palette_name, list_as_list = list_as_list, overwrite = overwrite)
  return(seurat_object)
}
