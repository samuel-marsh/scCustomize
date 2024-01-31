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
    if (length(x = bad_assays) > 0 && isTRUE(x = omit_warn)) {
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
  if (isTRUE(x = print_msg)) {
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
  if (isTRUE(x = and)) {
    last_sep <- " and "
  } else {
    last_sep <- " or "
  }

  if (input_length <= 3) {
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


#' Ask yes/no question to proceed
#'
#' Asks the user to answer yes/no question and returns logical value depending on
#' the answer.
#'
#' @return logical
#'
#' @references function modified from function in devtools R package (License: MIT) \url{https://github.com/r-lib/devtools}.
#' @details \url{https://github.com/r-lib/devtools/blob/9f27cc3e6335e74d6f51ed331509ebda56747901/R/release.R#L147-L156}.
#'
#' @import cli
#'
#' @noRd
#'

yesno <- function(msg, .envir = parent.frame()) {
  yeses <- c("Yes")
  nos <- c("No")

  cli_inform(message = msg, .envir = .envir)
  qs <- c("Yes", "No")
  rand <- sample(length(qs))

  utils::menu(qs[rand]) != which(rand == 1)
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
#' Retrieves Ensembl IDs for ribosomal genes
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
#' Retrieves species specific gene lists for MSigDB QC Hallmark lists: "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
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
     oxphos = msigdb_qc_gene_list[[oxphos]],
     apop = msigdb_qc_gene_list[[apop]],
     dna_repair = msigdb_qc_gene_list[[dna_repair]]
   )

   return(qc_gene_list)
 }


#' Retrieve IEG Gene Lists
#'
#' Retrieves species specific IEG gene lists
#'
#' @param species species to retrieve IDs.
#'
#' @return list of 2 sets of gene_symbols
#'
#' @import cli
#'
#' @keywords internal
#'
#' @noRd
#'

 Retrieve_IEG_Lists <- function(
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

   if (species %in% c(marmoset_options, zebrafish_options, rat_options, drosophila_options, macaque_options)) {
     cli_abort(message = "Rat, Marmoset, Macaque, Zebrafish, and Drosophila are not currently supported.")
   }

   # set prefix
   if (species %in% mouse_options) {
     prefix <- "Mus_musculus_"
   }
   if (species %in% human_options) {
     prefix <- "Homo_sapiens_"
   }

   # set list names
   ieg <- paste0(prefix, "IEG")

   # pull lists
   qc_gene_list <- list(
     ieg = ieg_gene_list[[ieg]]
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
#' @return Seurat object
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
     cli_inform(message = "The supplied species ({.field {species}}) is not currently supported.")
   }

   # Check Seurat
   Is_Seurat(seurat_object = seurat_object)

   # Check name collision
   if (any(duplicated(x = c(oxphos_name, apop_name, dna_repair_name)))) {
     cli_abort(message = "One or more of values provided to {.code oxphos_name}, {.code apop_name}, {.code dna_repair_name} are identical.")
   }

   # Overwrite check
   if (oxphos_name %in% colnames(x = seurat_object@meta.data) || apop_name %in% colnames(x = seurat_object@meta.data) || dna_repair_name %in% colnames(x = seurat_object@meta.data)) {
     if (isFALSE(x = overwrite)) {
       cli_abort(message = c("Columns with {.val {oxphos_name}} and/or {.val {apop_name}} already present in meta.data slot.",
                             "i" = "*To run function and overwrite columns set parameter {.code overwrite = TRUE} or change respective {.code oxphos_name}, {.code apop_name}, and/or {.code dna_repair_name}*")
       )
     }
     cli_inform(message = c("Columns with {.val {oxphos_name}} and/or {.val {apop_name}} already present in meta.data slot.",
                            "i" = "Overwriting those columns as {.code overwrite = TRUE.}")
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

   # Log Command
   seurat_object <- LogSeuratCommand(object = seurat_object)

   # return final object
   return(seurat_object)
 }



#' Add IEG Gene List Percentages
#'
#' Adds percentage of counts from IEG genes from mouse and human.
#'
#' @param seurat_object object name.
#' @param species Species of origin for given Seurat Object.  Only accepted species are: mouse, human (name or abbreviation).
#' @param ieg_name name to use for the new meta.data column containing percent IEG gene counts. Default is "percent_ieg".
#' @param assay Assay to use (default is the current object default assay).
#' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
#' function will abort if columns with the name provided to `ieg_name` is present in meta.data slot.
#'
#' @return Seurat object
#'
#' @import cli
#'
#' @keywords internal
#'
#' @noRd
#'


 Add_IEG_Seurat <- function(
    seurat_object,
    species,
    ieg_name = "percent_ieg",
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
     cli_inform(message = "The supplied species ({.field {species}}) is not currently supported.")
   }

   # Check Seurat
   Is_Seurat(seurat_object = seurat_object)

   # Overwrite check
   if (ieg_name %in% colnames(x = seurat_object@meta.data)) {
     if (isFALSE(x = overwrite)) {
       cli_abort(message = c("Column with {.val {ieg_name}} already present in meta.data slot.",
                             "i" = "*To run function and overwrite column set parameter {.code overwrite = TRUE} or change respective {.code ieg_name}*")
       )
     }
     cli_inform(message = c("Column with {.val {ieg_name}} already present in meta.data slot.",
                            "i" = "Overwriting those column as {.code overwrite = TRUE.}")
     )
   }

   # Set default assay
   assay <- assay %||% DefaultAssay(object = seurat_object)

   # Retrieve gene lists
   ieg_gene_list <- Retrieve_IEG_Lists(species = species)

   ieg_found <- Feature_PreCheck(object = seurat_object, features = ieg_gene_list[["ieg"]])

   # Add mito and ribo columns
   if (length(x = ieg_found) > 0) {
     seurat_object[[ieg_name]] <- PercentageFeatureSet(object = seurat_object, features = ieg_found, assay = assay)
   }

   # Log Command
   seurat_object <- LogSeuratCommand(object = seurat_object)

   # return final object
   return(seurat_object)
 }


#' Return default QC features
#'
#' Returns default QC features full names when provided with shortcut name.
#'
#' @param seurat_object object name.
#' @param features vector of features to check against defaults.
#' @param print_defaults return the potential accepted default values.
#'
#' @return list of found and not found features
#'
#' @import cli
#'
#' @keywords internal
#'
#' @noRd
#'

 Return_QC_Defaults <- function(
    seurat_object,
    features,
    print_defaults = FALSE
 ) {
   # default values
   feature_defaults <- list(
     feature = c("features", "Features", "genes", "Genes"),
     UMIs = c("counts", "Counts", "umis", "umi", "UMI", "UMIs", "UMIS"),
     mito = c("mito", "Mito"),
     ribo = c("ribo", "Ribo"),
     mito_ribo = c("mito_ribo", "Mito_Ribo"),
     complexity = c("complexity", "Complexity"),
     top_pct = c("top_pct", "Top_Pct"),
     IEG = c("ieg", "IEG"),
     OXPHOS = c("oxphos", "OXPHOS"),
     APOP = c("apop", "Apop"),
     DNA_Repair = c("dna_repair", "DNA_Repair")
   )

   # if print is TRUE
   if (isTRUE(x = print_defaults)) {
     cli_inform(message = c("Accepted default values are:",
                            "{.field {glue_collapse_scCustom(input_string = unlist(feature_defaults), and = TRUE)}}"))
     stop_quietly()
   }

   # Assign values
   if (any(features %in% feature_defaults[[1]])) {
     default1 <- "nFeature_RNA"
   } else {
     default1 <- NULL
   }
   if (any(features %in% feature_defaults[[2]])) {
     default2 <- "nCount_RNA"
   } else {
     default2 <- NULL
   }
   if (any(features %in% feature_defaults[[3]])) {
     default3 <- "percent_mito"
   } else {
     default3 <- NULL
   }
   if (any(features %in% feature_defaults[[4]])) {
     default4 <- "percent_ribo"
   } else {
     default4 <- NULL
   }
   if (any(features %in% feature_defaults[[5]])) {
     default5 <- "percent_mito_ribo"
   } else {
     default5 <- NULL
   }
   if (any(features %in% feature_defaults[[6]])) {
     default6 <- "log10GenesPerUMI"
   } else {
     default6 <- NULL
   }
   if (any(features %in% feature_defaults[[7]])) {
     default7 <- grep(pattern = "percent_top", x = colnames(x = seurat_object@meta.data), value = TRUE)
   } else {
     default7 <- NULL
   }
   if (any(features %in% feature_defaults[[8]])) {
     default8 <- "percent_ieg"
   } else {
     default8 <- NULL
   }
   if (any(features %in% feature_defaults[[9]])) {
     default9 <- "percent_oxphos"
   } else {
     default9 <- NULL
   }
   if (any(features %in% feature_defaults[[10]])) {
     default10 <- "percent_apop"
   } else {
     default10 <- NULL
   }
   if (any(features %in% feature_defaults[[11]])) {
     default11 <- "percent_dna_repair"
   } else {
     default11 <- NULL
   }

   # All found defaults
   all_found_defaults <- c(default1, default2, default3, default4, default5, default6, default7, default8, default9, default10, default11)

   # get not found features
   not_found_defaults <- features[!features %in% unlist(feature_defaults)]

   # create return list
   feat_list <- list(
     found_defaults = all_found_defaults,
     not_found_defaults = not_found_defaults
   )

   # return feature list
   return(feat_list)
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


#' Whole number check
#'
#' Checks whether a number is whole
#'
#' @param x number to check
#'
#' @return NULL or error message if number is not whole
#'
#' @keywords internal
#'
#' @noRd
#'

check_whole_num <- function(
    x
) {
  round_x <- round(x = x)

  res <-  identical(x = x, y = round_x)

  return(res)
}


#' Remove single value columns of data.frame
#'
#' Checks all columns within data.frame and returns data.frame minus columns that have the same value in all rows.
#'
#' @param df data.frame to filter
#'
#' @references Code used in function has been slightly modified from `sceasy:::.regularise_df` function of
#' sceasy package \url{https://github.com/cellgeni/sceasy} (License: GPL-3).
#' Code modified to match scCustomize & tidyverse style, add error checks, and
#' add cli formatted messages.
#'
#' @import cli
#' @importFrom dplyr select all_of
#' @importFrom magrittr "%>%"
#'
#' @return data.frame
#'
#' @keywords internal
#'
#' @noRd
#'

drop_single_value_cols <- function(
    df
) {
  if (!inherits(what = "data.frame", x = df)) {
    cli_abort(message = "{.code df} must of be of class data.frame.")
  }

  single_val_columns <- sapply(df, function(x) {
    length(x = unique(x = x)) == 1
    })

  col_names_single <- df %>%
    select(which(single_val_columns)) %>%
    colnames()

  if (length(x = col_names_single) > 0) {
    cli_inform(message = c("The following columns were removed as they contain identical values for all rows:",
                           "i" = "{.field {col_names_single}}"))
  }

  # filter df
  df_filtered <- df %>%
    select(-(all_of(col_names_single)))

  # return df
  return(df_filtered)
}


#' Check valid color
#'
#' Checks if input values are valid colors representations in R.
#'
#' @param colors vector of color(s) to check
#'
#' @references Code for function \url{https://stackoverflow.com/a/13290832/15568251}.
#' Renamed by Samuel Marsh.
#'
#' @importFrom grDevices col2rgb
#'
#' @return logical named vector
#'
#' @keywords internal
#'
#' @noRd
#'

Is_Color <- function(
    color
) {
  sapply(colors, function(X) {
    tryCatch(is.matrix(x = col2rgb(col = X)),
             error = function(e) FALSE)
  })
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### METRICS HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
      filter(.data[["Grouped.By"]] == "Physical library ID" & .data[["Library.Type"]] == "Gene Expression") %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    GEX_metrics2 <- raw_data %>%
      filter(.data[["Metric.Name"]] %in% c(c("Median UMI counts per cell", "Median genes per cell", "Median reads per cell", "Total genes detected"))) %>%
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
      filter(.data[["Grouped.By"]]== "Physical library ID" & .data[["Library.Type"]] == "VDJ T") %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    VDJ_T_Metrics2 <- raw_data %>%
      filter(.data[["Metric.Name"]] %in% c("Cells with productive TRA contig", "Cells with productive TRB contig", "Cells with productive V-J spanning (TRA, TRB) pair", "Cells with productive V-J spanning pair", "Median TRA UMIs per Cell", "Median TRB UMIs per Cell", "Number of cells with productive V-J spanning pair", "Paired clonotype diversity")
      ) %>%
      select(all_of(c("Metric.Name", "Metric.Value"))) %>%
      column_to_rownames("Metric.Name") %>%
      t() %>%
      data.frame()

    raw_data_vdjt <- cbind(VDJ_T_Metrics, VDJ_T_Metrics2)

    column_numbers <- grep(pattern = ",", x = raw_data_vdjt[1, ])
    raw_data_vdjt[,c(column_numbers)] <- lapply(raw_data_vdjt[,c(column_numbers)],function(x){as.numeric(gsub(",", "", x))})

    column_numbers_pct <- grep(pattern = "%", x = raw_data_vdjt[1, ])
    all_columns <- 1:ncol(x = raw_data_vdjt)

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GENE NAME/FILE CACHE HELPERS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' BiocFileCache Interface
#'
#' Internal function to manage/call BiocFileCache.
#'
#' @return cache
#'
#' @references \url{https://bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html#cache-to-manage-package-data}
#'
#' @noRd
#'

.get_bioc_cache <- function(
) {
    cache <- tools::R_user_dir(package = "scCustomize", which="cache")
    BiocFileCache::BiocFileCache(cache)
  }


#' Download HGNC Dataset
#'
#' Internal function to download and cache the latest version of HGNC dataset for use with renaming genes.
#'
#' @param update logical, whether to manually override update parameters and download new data.
#'
#' @import cli
#'
#' @return path to data cache
#'
#' @references \url{https://bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html}
#'
#' @noRd
#'


download_hgnc_data <- function(
    update = NULL
) {
  # Get cache
  bfc <- .get_bioc_cache()

  # URL from https://www.genenames.org/download/statistics-and-files/
  hgnc_ftp_url <- "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt"

  # bfc <- BiocFileCache::BiocFileCache(hgnc_ftp_url)

  rid <- BiocFileCache::bfcquery(bfc, hgnc_ftp_url, "fpath")$rid
  if (!length(rid)) {               # not in cache, add but do not download
    rid <- names(BiocFileCache::bfcadd(bfc, hgnc_ftp_url, download = FALSE))
  }

  if (isTRUE(x = update)) {
    update <- update
  } else {
    update <- BiocFileCache::bfcneedsupdate(bfc, rid)  # TRUE if newly added or stale
  }

  # download & process
  if (!isFALSE(x = update)) {
    cli_inform(message = "Downloading HGNC data from: {.field {hgnc_ftp_url}}")
    BiocFileCache::bfcdownload(bfc, rid, ask = FALSE, FUN = process_hgnc_data)
  }

  rpath <- BiocFileCache::bfcrpath(bfc, rids=rid)    # path to processed result

  return(rpath)
}


#' Process HGNC Dataset
#'
#' Internal function process/filter and save HGNC dataset during cache process
#'
#' @param from input (cache location).
#' @param to output (cached data).
#'
#' @importFrom dplyr mutate select filter any_of contains
#' @importFrom magrittr "%>%"
#' @importFrom tidyr separate_wider_delim pivot_longer
#'
#' @return path to data cache
#'
#' @references \url{https://bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html}
#'
#' @noRd
#'

process_hgnc_data <- function(
    from,
    to
) {
  # read in data
  hgnc_full_data <- data.table::fread(file = from, data.table = F)

  # filter data: Approved Genes > select relevant categories
  hgnc_filtered_data <- hgnc_full_data %>%
    filter(.data[["status"]] == "Approved") %>%
    select(any_of(c("hgnc_id", "symbol", "status", "alias_symbol", "prev_symbol", "date_symbol_changed", "entrez_id", "ensembl_gene_id")))

  # Select needed for renaming > split prev symbol column by number of additional columns needed > pivot wider without NAs > mutate
  hgnc_long_data <- hgnc_filtered_data %>%
    select(any_of(c("symbol", "prev_symbol"))) %>%
    separate_wider_delim(cols = "prev_symbol", delim = "|", names_sep = "_", names = NULL, too_few = "align_start") %>%
    pivot_longer(cols = contains("_symbol"),
                 names_to = "column",
                 values_to = "prev_symbol",
                 values_drop_na = TRUE) %>%
    mutate("prev_symbol" = ifelse(.data[["prev_symbol"]] %in% "", .data[["symbol"]], .data[["prev_symbol"]]))

  # save processed data
  saveRDS(hgnc_long_data, file = to)
  TRUE
}
