#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### Operators ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Set a default value if an object is NOT null
#'
#' @param lhs An object to set if it's NOT null
#' @param rhs The value to provide if x is NOT null
#'
#' @return lhs if lhs is null, else rhs
#'
#' @author Hadley Wickham
#' @references \url{https://adv-r.hadley.nz/functions.html#missing-arguments}
#'
#' @noRd
#'

`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}


#' Set a default value if an object is null
#'
#' @param lhs An object to set if it's null
#' @param rhs The value to provide if x is null
#'
#' @return rhs if lhs is null, else lhs
#'
#' @author Hadley Wickham
#' @references \url{https://adv-r.hadley.nz/functions.html#missing-arguments}
#'
#' @noRd
#'

`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### Object Checks ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Check Seurat Object
#'
#' Checks if object is of class: Seurat and returns error message if not.
#'
#' @param seurat_object Seurat object name.
#'
#' @import cli
#'
#' @return error if not Seurat object.
#'
#' @noRd
#'

Is_Seurat <- function(
  seurat_object
) {
  if (!inherits(what = "Seurat", x = seurat_object)) {
    cli_abort(message = "{.code seurat_object} provided is not an object of class: Seurat.")
  }
}


#' Check LIGER Object
#'
#' Checks if object is of class: liger and returns error message if not.
#'
#' @param liger_object liger object name.
#'
#' @import cli
#'
#' @return error is not LIGER object
#'
#' @noRd
#'

Is_LIGER <- function(
  liger_object
) {
  if (class(x = liger_object)[[1]] != "liger") {
    cli_abort(message = "{.code liger_object} provided is not an object of class: liger")
  }
}


#' Check for assays present in object
#
#' Checks Seurat object for the presence of assays against list of specified assays
#'
#' @param seurat_object Seurat object name.
#' @param assay_list vector of genes to check.
#' @param print_msg logical. Whether message should be printed if all features are found.  Default is TRUE.
#' @param omit_warn logical. Whether to print message about features that are not found in current object.
#'
#' @import cli
#'
#' @return List of found vs not found assays.
#'
#' @noRd
#'

Assay_Present <- function(
  seurat_object,
  assay_list,
  print_msg = TRUE,
  omit_warn = TRUE
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # get all features
  possible_assays <- Assays(object = seurat_object)

  # If any features not found
  if (any(!assay_list %in% possible_assays)) {
    bad_assays <- assay_list[!assay_list %in% possible_assays]
    found_assays <- assay_list[assay_list %in% possible_assays]
    if (length(x = found_assays) == 0) {
      cli_abort(message = "No requested assays found.")
    }

    # Return message of assays not found
    if (length(x = bad_assays) > 0 && omit_warn) {
      cli_warn(message = c("The following assays were omitted as they were not found:",
                           "i" = "{.field {glue_collapse_scCustom(input_string = bad_assays, and = TRUE)}}.")
      )
    }

    # Combine into list and return
    assay_list <- list(
      found_assays = found_assays,
      bad_assays = bad_assays
    )
    return(assay_list)
  }

  # Print all found message if TRUE
  if (print_msg) {
    cli_inform(message = "All assays present.")
  }

  # Return full input gene list.
  # Combine into list and return
  assay_list <- list(
    found_assays = assay_list,
    bad_assays = NULL
  )
  return(assay_list)
}


#' Check whether assay is V5
#
#' Checks Seurat object to verify whether it is composed of "Assay" or "Assay5" slots.
#'
#' @param seurat_object Seurat object name.
#' @param assay name of assay to check, default is NULL.
#'
#' @return TRUE if seurat_object contains "Assay5" class.
#'
#' @noRd
#'

Assay5_Check <- function(
    seurat_object,
    assay = NULL
){
  assay <- assay %||% DefaultAssay(object = seurat_object)

  if (inherits(x = seurat_object@assays[[assay]], what = "Assay")) {
    return(FALSE)
  }
  if (inherits(x = seurat_object@assays[[assay]], what = "Assay5")) {
    return(TRUE)
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### WARN/ERROR MESSAGING ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Stop function without error message
#'
#' Modifies R options within the function call only to hide the error message from `stop` while
#' keeping global options preserved outside of function.
#'
#' @return stops function without error message
#'
#' @author Stibu
#' @references \url{https://stackoverflow.com/a/42945293/15568251}
#' @details \url{https://creativecommons.org/licenses/by-sa/3.0/}
#'
#' @noRd
#'

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}


#' Custom glue collapse
#
#' Customized glue_collapse that is based on the number of items in input list
#'
#' @param input_string input string to be collapsed.
#' @param and logical.  Whether to use "and" or "or" in the collapsed string
#'
#' @return collapsed string
#'
#' @importFrom glue glue_collapse
#'
#' @noRd
#'

glue_collapse_scCustom <- function(
  input_string,
  and = TRUE
) {
  # Check length of input string
  input_length <- length(x = input_string)

  # set last seperator
  if (and) {
    last_sep <- " and "
  } else {
    last_sep <- " or "
  }

  if (input_length < 3) {
    glue_collapse(x = input_string, sep = ", ", last = last_sep)
  } else {
    glue_collapse(x = input_string, sep = ", ", last = paste0(",", last_sep))
  }
}


#' Perform Feature and Meta Checks before plotting
#'
#' Wraps the `Gene_Present`, `Meta_Present`, `Reduction_Loading_Present`, and `Case_Check` into
#' single function to perform feature checks before plotting.
#'
#' @param object Seurat object
#' @param features vector of features and/or meta data variables to plot.
#' @param assay Assay to use (default is the current object default assay).
#'
#' @return vector of features and/or meta data that were found in object.
#'
#' @noRd
#'
#' @keywords internal
#'

Feature_PreCheck <- function(
    object,
    features,
    assay = NULL
) {
  # set assay (if null set to active assay)
  assay <- assay %||% DefaultAssay(object = object)

  # Check features and meta to determine which features present
  features_list <- Gene_Present(data = object, gene_list = features, omit_warn = FALSE, print_msg = FALSE, case_check_msg = FALSE, return_none = TRUE, seurat_assay = assay)

  meta_list <- Meta_Present(seurat_object = object, meta_col_names = features_list[[2]], omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)

  reduction_list <- Reduction_Loading_Present(seurat_object = object, reduction_names = meta_list[[2]], omit_warn = FALSE, print_msg = FALSE, return_none = TRUE)

  all_not_found_features <- reduction_list[[2]]

  all_found_features <- c(features_list[[1]], meta_list[[1]], reduction_list[[1]])

  # Stop if no features found
  if (length(x = all_found_features) < 1) {
    cli_abort(message = c("No features were found.",
                          "*" = "The following are not present in object:",
                          "i" = "{.field {glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}}")
    )
  }

  # Return message of features not found
  if (length(x = all_not_found_features) > 0) {
    op <- options(warn = 1)
    on.exit(options(op))
    cli_warn(message = c("The following features were omitted as they were not found:",
                         "i" = "{.field {glue_collapse_scCustom(input_string = all_not_found_features, and = TRUE)}}")
    )
  }

  # Check feature case and message if found
  Case_Check(seurat_object = object, gene_list = all_not_found_features, case_check_msg = TRUE, return_features = FALSE)

  # return all found features
  return(all_found_features)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### QC HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Ensembl Mito IDs
#'
#' Retrieves Ensembl IDs for mitochondrial genes
#'
#' @param species species to retrieve IDs.
#'
#' @return vector of Ensembl Gene IDs
#'
#' @import cli
#'
#' @keywords internal
#'
#' @noRd
#'

Retrieve_Ensembl_Mito <- function(
    species
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

  if (species %in% marmoset_options) {
    cli_abort(message = "Marmoset mitochondrial genome is not part of current Ensembl build.")
  }

  if (species %in% mouse_options) {
    mito_ensembl <- ensembl_mito_id$Mus_musculus_mito_ensembl
  }
  if (species %in% human_options) {
    mito_ensembl <- ensembl_mito_id$Homo_sapiens_mito_ensembl
  }
  if (species %in% zebrafish_options) {
    mito_ensembl <- ensembl_mito_id$Danio_rerio_mito_ensembl
  }
  if (species %in% rat_options) {
    mito_ensembl <- ensembl_mito_id$Rattus_norvegicus_mito_ensembl
  }
  if (species %in% drosophila_options) {
    mito_ensembl <- ensembl_mito_id$Drosophila_melanogaster_mito_ensembl
  }
  if (species %in% macaque_options) {
    mito_ensembl <- ensembl_mito_id$Macaca_mulatta_mito_ensembl
  }

  return(mito_ensembl)
}


#' Ensembl Ribo IDs
#'
#' Retrieves Ensembl IDs for ribsomal genes
#'
#' @param species species to retrieve IDs.
#'
#' @return vector of Ensembl Gene IDs
#'
#' @import cli
#'
#' @keywords internal
#'
#' @noRd
#'

Retrieve_Ensembl_Ribo <- function(
    species
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

  if (species %in% mouse_options) {
    ribo_ensembl <- ensembl_ribo_id$Mus_musculus_ribo_ensembl
  }
  if (species %in% human_options) {
    ribo_ensembl <- ensembl_ribo_id$Homo_sapiens_ribo_ensembl
  }
  if (species %in% zebrafish_options) {
    ribo_ensembl <- ensembl_ribo_id$Callithrix_jacchus_ribo_ensembl
  }
  if (species %in% zebrafish_options) {
    ribo_ensembl <- ensembl_ribo_id$Danio_rerio_ribo_ensembl
  }
  if (species %in% rat_options) {
    ribo_ensembl <- ensembl_ribo_id$Rattus_norvegicus_ribo_ensembl
  }
  if (species %in% drosophila_options) {
    ribo_ensembl <- ensembl_ribo_id$Drosophila_melanogaster_ribo_ensembl
  }
  if (species %in% macaque_options) {
    ribo_ensembl <- ensembl_ribo_id$Macaca_mulatta_ribo_ensembl
  }

  return(ribo_ensembl)
}


#' Retrieve MSigDB Gene Lists
#'
#' Retrieves species specifc gene lists for MSigDB QC Hallmark lists: "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
#' "HALLMARK_APOPTOSIS", and "HALLMARK_DNA_REPAIR".
#'
#' @param species species to retrieve IDs.
#'
#' @return list of 3 sets of gene_symbols
#'
#' @import cli
#'
#' @keywords internal
#'
#' @noRd
#'

 Retrieve_MSigDB_Lists <- function(
     species
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

   if (species %in% marmoset_options) {
     cli_abort(message = "Marmoset is not currently a part of MSigDB gene list database.")
   }

   # set prefix
   if (species %in% mouse_options) {
     prefix <- "Mus_musculus_"
   }
   if (species %in% human_options) {
     prefix <- "Homo_sapiens_"
   }
   if (species %in% zebrafish_options) {
     prefix <- "Dario_rerio_"
   }
   if (species %in% rat_options) {
     prefix <- "Rattus_norvegicus_"
   }
   if (species %in% drosophila_options) {
     prefix <- "Drosophila_melanogaster_"
   }
   if (species %in% macaque_options) {
     prefix <- "Macaca_mulatta_"
   }

   # set list names
   oxphos <- paste0(prefix, "msigdb_oxphos")
   apop <- paste0(prefix, "msigdb_apop")
   dna_repair <- paste0(prefix, "msigdb_dna_repair")

   # pull lists
   qc_gene_list <- list(
     oxphos <- msigdb_qc_gene_list[[oxphos]],
     apop <- msigdb_qc_gene_list[[apop]],
     dna_repair <- msigdb_qc_gene_list[[dna_repair]]
   )

   return(qc_gene_list)
 }


 #' Add MSigDB Gene Lists Percentages
 #'
 #' Adds percentage of counts from 3 hallmark MSigDB hallmark gene sets: "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
 #' "HALLMARK_APOPTOSIS", and "HALLMARK_DNA_REPAIR".
 #'
 #' @param seurat_object object name.
 #' @param species Species of origin for given Seurat Object.  Only accepted species are: mouse, human,
 #' zebrafish, rat, drosophila, or rhesus macaque (name or abbreviation)
 #' @param oxphos_name name to use for the new meta.data column containing percent MSigDB Hallmark oxidative
 #' phosphorylation counts. Default is "percent_oxphos".
 #' @param apop_name name to use for the new meta.data column containing percent MSigDB Hallmark apoptosis counts.
 #' Default is "percent_apop".
 #' @param dna_repair_name name to use for the new meta.data column containing percent MSigDB Hallmark DNA repair counts.
 #' Default is "percent_oxphos".
 #' @param assay Assay to use (default is the current object default assay).
 #' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
 #' function will abort if columns with any one of the names provided to `mito_name` `ribo_name` or
 #' `mito_ribo_name` is present in meta.data slot.
 #'
 #' @return list of 3 sets of gene_symbols
 #'
 #' @import cli
 #'
 #' @keywords internal
 #'
 #' @noRd
 #'


 Add_MSigDB_Seurat <- function(
     seurat_object,
     species,
     oxphos_name = "percent_oxphos",
     apop_name = "percent_apop",
     dna_repair_name = "percent_dna_repair",
     assay = NULL,
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
     Macaque_Options = c("Macaque", "macaque", "Rhesus", "macaca", "mmulatta", NA)
   )

   if (!species %in% unlist(x = accepted_names)) {
     cli::cli_inform(message = "The supplied species ({.field {species}}) is not currently supported.")
   }

   # Check Seurat
   Is_Seurat(seurat_object = seurat_object)

   # Check name collision
   if (any(duplicated(x = c(oxphos_name, apop_name, dna_repair_name)))) {
     cli_abort(message = "One or more of values provided to {.code oxphos_name}, {.code apop_name}, {.code dna_repair_name} are identical.")
   }

   # Overwrite check
   if (oxphos_name %in% colnames(x = seurat_object@meta.data) || apop_name %in% colnames(x = seurat_object@meta.data) || dna_repair_name %in% colnames(x = seurat_object@meta.data)) {
     if (!overwrite) {
       cli_abort(message = c("Columns with {.val {oxphos_name}} and/or {.val {apop_name}} already present in meta.data slot.",
                             "i" = "*To run function and overwrite columns set parameter {.code overwrite = TRUE} or change respective {.code oxphos_name}, {.code apop_name}, and/or {.code dna_repair_name}*")
       )
     }
     cli_inform(message = c("Columns with {.val {oxphos_name}} and/or {.val {apop_name}} already present in meta.data slot.",
                            "i" = "Overwriting those columns as .code {overwrite = TRUE.}")
     )
   }

   # Set default assay
   assay <- assay %||% DefaultAssay(object = seurat_object)

   # Retrieve gene lists
   msigdb_gene_list <- Retrieve_MSigDB_Lists(species = species)

   oxphos_found <- Feature_PreCheck(object = seurat_object, features = msigdb_gene_list[["oxphos"]])
   apop_found <- Feature_PreCheck(object = seurat_object, features = msigdb_gene_list[["apop"]])
   dna_repair_found <- Feature_PreCheck(object = seurat_object, features = msigdb_gene_list[["dna_repair"]])

   # Add mito and ribo columns
   if (length(x = oxphos_found) > 0) {
     seurat_object[[oxphos_name]] <- PercentageFeatureSet(object = seurat_object, features = oxphos_found, assay = assay)
   }
   if (length(x = apop_found) > 0) {
     seurat_object[[apop_name]] <- PercentageFeatureSet(object = seurat_object, features = apop_found, assay = assay)
   }
   if (length(x = dna_repair_found) > 0) {
     seurat_object[[dna_repair_name]] <- PercentageFeatureSet(object = seurat_object, features = dna_repair_found, assay = assay)
   }

   # return final object
   return(seurat_object)
 }


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GENERAL HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Calculate the percentage of a vector above some threshold
#'
#' @param x Vector of values
#' @param threshold Threshold to use when calculating percentage
#'
#' @return Returns the percentage of `x` values above the given threshold
#'
#' @author Satija Lab & all co-Authors of Seurat Package
#' @references See Utilities.R in source code of Seurat \url{https://github.com/satijalab/seurat/blob/master/R/utilities.R}  (Licence: GPL-3).
#'
#' @note modified from Seurat version to return a percentage instead of proportion/decimal as part of `Percent_Expressing` function.  To be replaced following next Seurat version update.
#'
#' @keywords internal
#'
#' @noRd
#'

PercentAbove_Seurat <- function(x, threshold) {
  return((length(x = x[x > threshold]) / length(x = x))*100)
}


#' Extract delimiter information from a string.
#'
#' Parses a string (usually a cell name) and extracts fields based on a delimiter
#'
#' @param string String to parse.
#' @param field Integer(s) indicating which field(s) to extract. Can be a vector multiple numbers.
#' @param delim Delimiter to use, set to underscore by default.
#'
#' @return A new string, that parses out the requested fields, and (if multiple), rejoins them with the same delimiter
#'
#' @references See Utilities.R in source code of Seurat \url{https://github.com/satijalab/seurat/blob/master/R/utilities.R}  (Licence: GPL-3).
#'
#' @keywords internal
#'
#' @noRd
#'

ExtractField <- function(string, field = 1, delim = "_") {
  fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
  if (length(x = fields) == 1) {
    return(strsplit(x = string, split = delim)[[1]][field])
  }
  return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}


#' Calculate the middle value between two numbers
#'
#' @param min Lower value.
#' @param max Higher value.
#'
#' @return Returns number in middle of two provided values.
#'
#' @references Code to calculate middle value from: adapted from: \url{https://stackoverflow.com/a/54147509/15568251}.
#' Renamed and wrapped into function by Samuel Marsh.
#'
#' @keywords internal
#'
#' @noRd
#'

Middle_Number <- function(
  min,
  max
) {
  min_max <- c(min, max)
  middle <- min_max[-length(min_max)] + diff(min_max) / 2
  return(middle)
}


#' Symmetrical setdiff
#'
#' tests for differences between two vectors symmetrically.
#'
#' @param x first vector to test
#' @param y second vector to test
#'
#' @return vector differences x vs. y and y vs. x
#'
#' @references Function name and code from R-bloggers post:
#' \url{https://www.r-bloggers.com/2013/06/symmetric-set-differences-in-r/}
#'
#' @keywords internal
#'
#' @noRd
#'

symdiff <- function(
  x,
  y
) {
  setdiff(x = union(x = x, y = y), intersect(x = x, y = y))
}


#' Read Gene Expression Statistics from 10X Cell Ranger Count
#'
#' Get data.frame with all metrics from the Cell Ranger `count` analysis (present in web_summary.html)
#'
#' @param base_path path to the parent directory which contains all of the subdirectories of interest.
#' @param secondary_path path from the parent directory to count "outs/" folder which contains the
#' "metrics_summary.csv" file.
#' @param default_10X logical (default TRUE) sets the secondary path variable to the default 10X directory structure.
#' @param lib_list a list of sample names (matching directory names) to import.  If `NULL` will read
#' in all samples in parent directory.
#' @param lib_names a set of sample names to use for each sample.  If `NULL` will set names to the
#' directory name of each sample.
#'
#' @return A data frame with sample metrics produced by Cell Ranger `count` pipeline.
#'
#' @import cli
#' @import pbapply
#' @importFrom dplyr bind_rows setdiff
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' count_metrics <- Metrics_Count_GEX(base_path = base_path, lib_list = lib_list, secondary_path = secondary_path, lib_names = lib_names)
#' }
#'

Metrics_Count_GEX <- function(
    lib_list,
    base_path,
    secondary_path,
    lib_names
){
  cli_inform(message = "Reading {.field Gene Expression} Metrics")
  raw_data_list <- pblapply(1:length(x = lib_list), function(x) {
    if (is.null(x = secondary_path)) {
      file_path <- file.path(base_path, lib_list[x])
    } else {
      file_path <- file.path(base_path, lib_list[x], secondary_path)
    }

    raw_data <- read.csv(file = file.path(file_path, "metrics_summary.csv"), stringsAsFactors = F)
    # Change format of numeric columns to due commas in data csv output.
    column_numbers <- grep(pattern = ",", x = raw_data[1, ])
    raw_data[,c(column_numbers)] <- lapply(raw_data[,c(column_numbers)],function(x){as.numeric(gsub(",", "", x))})


    column_numbers_pct <- grep(pattern = "%", x = raw_data[1, ])
    all_columns <- 1:ncol(x = raw_data)

    column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

    raw_data[,c(column_numbers_numeric)] <- lapply(raw_data[,c(column_numbers_numeric)],function(x){as.numeric(x)})

    return(raw_data)
  })

  # Name the list items
  if (is.null(x = lib_names)) {
    names(x = raw_data_list) <- lib_list
  } else {
    names(x = raw_data_list) <- lib_names
  }

  # Combine the list and add sample_id column
  full_data <- bind_rows(raw_data_list, .id = "sample_id")

  # Change column nams to use "_" separator instead of "." for readability
  colnames(x = full_data) <- gsub(pattern = "\\.", replacement = "_", x = colnames(x = full_data))

  rownames(x = full_data) <- full_data$sample_id

  return(full_data)

}


#' Read Gene Expression Statistics from 10X Cell Ranger multi
#'
#' Get data.frame with all gene expression metrics from the Cell Ranger `multi` analysis (present in web_summary.html)
#'
#' @param base_path path to the parent directory which contains all of the subdirectories of interest.
#' @param secondary_path path from the parent directory to count "outs/" folder which contains the
#' "metrics_summary.csv" file.
#' @param default_10X logical (default TRUE) sets the secondary path variable to the default 10X directory structure.
#' @param lib_list a list of sample names (matching directory names) to import.  If `NULL` will read
#' in all samples in parent directory.
#' @param lib_names a set of sample names to use for each sample.  If `NULL` will set names to the
#' directory name of each sample.
#'
#' @return A data frame with sample gene expression metrics produced by Cell Ranger `multi` pipeline.
#'
#' @import cli
#' @import pbapply
#' @importFrom dplyr all_of bind_rows filter rename select setdiff
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' count_multi_metrics <- Metrics_Multi_GEX(base_path = base_path, lib_list = lib_list, secondary_path = secondary_path, lib_names = lib_names)
#' }
#'

Metrics_Multi_GEX <- function(
    lib_list,
    base_path,
    secondary_path,
    lib_names
){
  cli_inform(message = "Reading {.field Gene Expression} Metrics")

  raw_data_list <- pblapply(1:length(x = lib_list), function(x) {
    if (is.null(x = secondary_path)) {
      file_path <- file.path(base_path, lib_list[x])
    } else {
      file_path <- file.path(base_path, lib_list[x], secondary_path, lib_list[x])
    }

    raw_data <- read.csv(file = file.path(file_path, "metrics_summary.csv"), stringsAsFactors = F)

    # Change format to column based and select relevant metrics
    GEX_metrics <- raw_data %>%
      filter(Grouped.By == "Physical library ID" & Library.Type == "Gene Expression") %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    GEX_metrics2 <- raw_data %>%
      filter(Metric.Name %in% c(c("Median UMI counts per cell", "Median genes per cell", "Median reads per cell", "Total genes detected"))) %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    raw_data_gex <- cbind(GEX_metrics, GEX_metrics2)

    # Change format of numeric columns to due commas in data csv output.
    column_numbers <- grep(pattern = ",", x = raw_data_gex[1, ])
    raw_data_gex[,c(column_numbers)] <- lapply(raw_data_gex[,c(column_numbers)],function(x){as.numeric(gsub(",", "", x))})

    # Rename multi columns to match names from count
    names_to_replace <- c(Reads.Mapped.to.Genome = "Mapped.to.genome",
                          Reads.Mapped.Confidently.to.Genome = "Confidently.mapped.to.genome",
                          Reads.Mapped.Confidently.to.Intergenic.Regions = "Confidently.mapped.to.intergenic.regions",
                          Reads.Mapped.Confidently.to.Intronic.Regions = "Confidently.mapped.to.intronic.regions",
                          Reads.Mapped.Confidently.to.Exonic.Regions = "Confidently.mapped.to.exonic.regions",
                          Reads.Mapped.Confidently.to.Transcriptome = "Confidently.mapped.to.transcriptome",
                          Reads.Mapped.Antisense.to.Gene = "Confidently.mapped.antisense",
                          Fraction.Reads.in.Cells = "Confidently.mapped.reads.in.cells",
                          Estimated.Number.of.Cells = "Estimated.number.of.cells",
                          Mean.Reads.per.Cell = "Mean.reads.per.cell",
                          Median.Genes.per.Cell = "Median.genes.per.cell",
                          Number.of.Reads = "Number.of.reads",
                          Valid.Barcodes = "Valid.barcodes",
                          Sequencing.Saturation = "Sequencing.saturation",
                          Total.Genes.Detected = "Total.genes.detected",
                          Median.UMI.Counts.per.Cell = "Median.UMI.counts.per.cell")

    raw_data_gex <- raw_data_gex %>%
      rename(all_of(names_to_replace))

    column_numbers_pct <- grep(pattern = "%", x = raw_data_gex[1, ])
    all_columns <- 1:ncol(x = raw_data_gex)

    column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

    raw_data_gex[,c(column_numbers_numeric)] <- lapply(raw_data_gex[,c(column_numbers_numeric)],function(x){as.numeric(x)})

    return(raw_data_gex)
  })

  # Name the list items
  if (is.null(x = lib_names)) {
    names(x = raw_data_list) <- lib_list
  } else {
    names(x = raw_data_list) <- lib_names
  }

  # Combine the list and add sample_id column
  full_data <- bind_rows(raw_data_list, .id = "sample_id")

  # Change column nams to use "_" separator instead of "." for readability
  colnames(x = full_data) <- gsub(pattern = "\\.", replacement = "_", x = colnames(x = full_data))

  rownames(x = full_data) <- full_data$sample_id

  return(full_data)
}


#' Read VDJ T Statistics from 10X Cell Ranger multi
#'
#' Get data.frame with all VDJ T metrics from the Cell Ranger `multi` analysis (present in web_summary.html)
#'
#' @param base_path path to the parent directory which contains all of the subdirectories of interest.
#' @param secondary_path path from the parent directory to count "outs/" folder which contains the
#' "metrics_summary.csv" file.
#' @param default_10X logical (default TRUE) sets the secondary path variable to the default 10X directory structure.
#' @param lib_list a list of sample names (matching directory names) to import.  If `NULL` will read
#' in all samples in parent directory.
#' @param lib_names a set of sample names to use for each sample.  If `NULL` will set names to the
#' directory name of each sample.
#'
#' @return A data frame with sample VDJ T metrics produced by Cell Ranger `multi` pipeline.
#'
#' @import cli
#' @import pbapply
#' @importFrom dplyr all_of bind_rows filter select setdiff
#' @importFrom magrittr "%>%"
#' @importFrom tibble column_to_rownames
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' vdj_multi_metrics <- Metrics_Multi_VDJT(base_path = base_path, lib_list = lib_list, secondary_path = secondary_path, lib_names = lib_names)
#' }
#'

Metrics_Multi_VDJT <- function(
    lib_list,
    base_path,
    secondary_path,
    lib_names
){
  cli_inform(message = "Reading {.field VDJ T} Metrics")

  raw_data_list <- pblapply(1:length(x = lib_list), function(x) {
    if (is.null(x = secondary_path)) {
      file_path <- file.path(base_path, lib_list[x])
    } else {
      file_path <- file.path(base_path, lib_list[x], secondary_path, lib_list[x])
    }

    raw_data <- read.csv(file = file.path(file_path, "metrics_summary.csv"), stringsAsFactors = F)

    VDJ_T_Metrics <- raw_data %>%
      filter(Grouped.By == "Physical library ID" & Library.Type == "VDJ T") %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    VDJ_T_Metrics2 <- raw_data %>%
      filter(Metric.Name %in% c("Cells with productive TRA contig", "Cells with productive TRB contig", "Cells with productive V-J spanning (TRA, TRB) pair", "Cells with productive V-J spanning pair", "Median TRA UMIs per Cell", "Median TRB UMIs per Cell", "Number of cells with productive V-J spanning pair", "Paired clonotype diversity")
      ) %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    raw_data_vdjt <- cbind(VDJ_T_Metrics, VDJ_T_Metrics2)

    column_numbers <- grep(pattern = ",", x = raw_data_vdjt[1, ])
    raw_data_vdjt[,c(column_numbers)] <- lapply(raw_data_vdjt[,c(column_numbers)],function(x){as.numeric(gsub(",", "", x))})

    column_numbers_pct <- grep(pattern = "%", x = raw_data_vdjt[1, ])
    all_columns <- 1:ncol(raw_data_vdjt)

    column_numbers_numeric <- setdiff(x = all_columns, y = column_numbers_pct)

    raw_data_vdjt[,c(column_numbers_numeric)] <- lapply(raw_data_vdjt[,c(column_numbers_numeric)],function(x){as.numeric(x)})

    return(raw_data_vdjt)
  })

  # Name the list items
  if (is.null(x = lib_names)) {
    names(x = raw_data_list) <- lib_list
  } else {
    names(x = raw_data_list) <- lib_names
  }

  # test_return <- lapply(1:length(test_return), function(i) {
  #   test_return[[i]]$Estimated.number.of.cells <- as.numeric(test_return[[i]]$Estimated.number.of.cells)
  # })

  # Combine the list and add sample_id column
  full_data <- bind_rows(raw_data_list, .id = "sample_id")

  # Change column nams to use "_" separator instead of "." for readability
  colnames(x = full_data) <- gsub(pattern = "\\.", replacement = "_", x = colnames(x = full_data))

  rownames(x = full_data) <- full_data$sample_id

  return(full_data)
}

