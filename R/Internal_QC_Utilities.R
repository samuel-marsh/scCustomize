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
  if (species %in% chicken_options) {
    mito_ensembl <- ensembl_mito_id$Gallus_gallus_mito_ensembl
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

  if (species %in% mouse_options) {
    ribo_ensembl <- ensembl_ribo_id$Mus_musculus_ribo_ensembl
  }
  if (species %in% human_options) {
    ribo_ensembl <- ensembl_ribo_id$Homo_sapiens_ribo_ensembl
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
  if (species %in% chicken_options) {
    ribo_ensembl <- ensembl_ribo_id$Gallus_gallus_ribo_ensembl
  }

  return(ribo_ensembl)
}


#' Ensembl Hemo IDs
#'
#' Retrieves Ensembl IDs for hemoglobin genes
#'
#' @param species species to retrieve IDs.
#'
#' @return vector of Ensembl Gene IDs
#'
#' @keywords internal
#'
#' @noRd
#'

Retrieve_Ensembl_Hemo <- function(
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

  if (species %in% mouse_options) {
    hemo_ensembl <- ensembl_hemo_id$Mus_musculus_hemo_ensembl
  }
  if (species %in% human_options) {
    hemo_ensembl <- ensembl_hemo_id$Homo_sapiens_hemo_ensembl
  }
  if (species %in% zebrafish_options) {
    hemo_ensembl <- ensembl_hemo_id$Danio_rerio_hemo_ensembl
  }
  if (species %in% rat_options) {
    hemo_ensembl <- ensembl_hemo_id$Rattus_norvegicus_hemo_ensembl
  }
  if (species %in% drosophila_options) {
    hemo_ensembl <- ensembl_hemo_id$Drosophila_melanogaster_hemo_ensembl
  }
  if (species %in% macaque_options) {
    hemo_ensembl <- ensembl_hemo_id$Macaca_mulatta_hemo_ensembl
  }
  if (species %in% chicken_options) {
    hemo_ensembl <- ensembl_hemo_id$Gallus_gallus_hemo_ensembl
  }

  return(hemo_ensembl)
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
  if (species %in% chicken_options) {
    prefix <- "Gallus_gallus_"
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


#' Retrieve MSigDB Ensembl Lists
#'
#' Retrieves species specific gene lists (ensembl IDs) for MSigDB QC Hallmark lists: "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
#' "HALLMARK_APOPTOSIS", and "HALLMARK_DNA_REPAIR".
#'
#' @param species species to retrieve IDs.
#'
#' @return list of 3 sets of ensembl IDs
#'
#' @keywords internal
#'
#' @noRd
#'

Retrieve_MSigDB_Ensembl_Lists <- function(
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
  if (species %in% chicken_options) {
    prefix <- "Gallus_gallus_"
  }

  # set list names
  oxphos <- paste0(prefix, "msigdb_oxphos")
  apop <- paste0(prefix, "msigdb_apop")
  dna_repair <- paste0(prefix, "msigdb_dna_repair")

  # pull lists
  qc_gene_list <- list(
    oxphos = msigdb_qc_ensembl_list[[oxphos]],
    apop = msigdb_qc_ensembl_list[[apop]],
    dna_repair = msigdb_qc_ensembl_list[[dna_repair]]
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
    cli_abort(message = "Rat, Marmoset, Macaque, Zebrafish, Drosophila, and Chicken are not currently supported.")
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


#' Retrieve IEG Gene Lists (Ensembl)
#'
#' Retrieves species specific IEG gene lists with ensembl IDs
#'
#' @param species species to retrieve IDs.
#'
#' @return list of 2 sets of ensembl IDs
#'
#' @keywords internal
#'
#' @noRd
#'

Retrieve_IEG_Ensembl_Lists <- function(
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
  ieg <- paste0(prefix, "IEG_ensembl")

  # pull lists
  qc_gene_list <- list(
    ieg = ensembl_ieg_list[[ieg]]
  )

  return(qc_gene_list)
}


#' Retrieve exAM Gene Lists
#'
#' Retrieves species specific exAM gene lists
#'
#' @param species species to retrieve IDs.
#'
#' @return list of 2 sets of gene_symbols
#'
#' @keywords internal
#'
#' @noRd
#'

Retrieve_exAM_Lists <- function(
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
    cli_abort(message = "Rat, Marmoset, Macaque, Zebrafish, Drosophila, and Chicken are not currently supported.")
  }

  # set prefix
  if (species %in% mouse_options) {
    prefix <- "Mus_musculus_"
  }
  if (species %in% human_options) {
    prefix <- "Homo_sapiens_"
  }

  if (species %in% human_options) {
    # set list names
    exAM1 <- paste0(prefix, "exAM_union")
    exAM2 <- paste0(prefix, "exAM_micro")

    # pull lists
    qc_gene_list <- list(
      "exAM_union" = exAM_gene_list[[exAM1]],
      "exAM_micro" = exAM_gene_list[[exAM2]]
    )
  } else {
    # set list names
    exAM1 <- paste0(prefix, "exAM_union")

    # pull lists
    qc_gene_list <- list(
      "exAM_union" = exAM_gene_list[[exAM1]]
    )
  }

  return(qc_gene_list)
}


#' Retrieve exAM Gene Lists (Ensembl)
#'
#' Retrieves species specific exAM gene lists with ensembl IDs
#'
#' @param species species to retrieve IDs.
#'
#' @return list of 2 sets of ensembl IDs
#'
#' @keywords internal
#'
#' @noRd
#'

Retrieve_exAM_Ensembl_Lists <- function(
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
    cli_abort(message = "Rat, Marmoset, Macaque, Zebrafish, and Drosophila are not currently supported.")
  }

  # set prefix
  if (species %in% mouse_options) {
    prefix <- "Mus_musculus_"
  }
  if (species %in% human_options) {
    prefix <- "Homo_sapiens_"
  }

  if (species %in% human_options) {
    # set list names
    exAM1 <- paste0(prefix, "exAM_union_ensembl")
    exAM2 <- paste0(prefix, "exAM_micro_ensembl")

    # pull lists
    qc_gene_list <- list(
      "exAM_union" = exAM_gene_list[[exAM1]],
      "exAM_micro" = exAM_gene_list[[exAM2]]
    )
  } else {
    # set list names
    exAM1 <- paste0(prefix, "exAM_union_ensembl")

    # pull lists
    qc_gene_list <- list(
      "exAM_union" = exAM_gene_list[[exAM1]]
    )
  }

  return(qc_gene_list)
}


#' Ensembl lncRNA IDs
#'
#' Retrieves Ensembl IDs for lncRNA genes
#'
#' @param species species to retrieve IDs.
#'
#' @return vector of Ensembl Gene IDs
#'
#' @keywords internal
#'
#' @noRd
#'

Retrieve_Ensembl_lncRNA <- function(
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

  if (species %in% mouse_options) {
    lncRNA_ensembl <- ensembl_lncRNA_id$Mus_musculus_lncRNA_ensembl
  }
  if (species %in% human_options) {
    lncRNA_ensembl <- ensembl_lncRNA_id$Homo_sapiens_lncRNA_ensembl
  }
  if (species %in% marmoset_options) {
    lncRNA_ensembl <- ensembl_lncRNA_id$Callithrix_jacchus_lncRNA_ensembl
  }
  if (species %in% zebrafish_options) {
    lncRNA_ensembl <- ensembl_lncRNA_id$Danio_rerio_lncRNA_ensembl
  }
  if (species %in% rat_options) {
    lncRNA_ensembl <- ensembl_lncRNA_id$Rattus_norvegicus_lncRNA_ensembl
  }
  if (species %in% macaque_options) {
    lncRNA_ensembl <- ensembl_lncRNA_id$Macaca_mulatta_lncRNA_ensembl
  }
  if (species %in% chicken_options) {
    lncRNA_ensembl <- ensembl_lncRNA_id$Gallus_gallus_lncRNA_ensembl
  }

  return(lncRNA_ensembl)
}


#' lncRNA IDs
#'
#' Retrieves gene symbol IDs for lncRNA genes
#'
#' @param species species to retrieve IDs.
#'
#' @return vector of Gene IDs
#'
#' @keywords internal
#'
#' @noRd
#'

Retrieve_lncRNA <- function(
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

  if (species %in% mouse_options) {
    lncRNA <- lncRNA_gene_list$Mus_musculus_lncRNA
  }
  if (species %in% human_options) {
    lncRNA <- lncRNA_gene_list$Homo_sapiens_lncRNA
  }
  if (species %in% zebrafish_options) {
    lncRNA <- lncRNA_gene_list$Danio_rerio_lncRNA
  }
  if (species %in% rat_options) {
    lncRNA <- lncRNA_gene_list$Rattus_norvegicus_lncRNA
  }
  if (species %in% macaque_options) {
    lncRNA <- lncRNA_gene_list$Macaca_mulatta_lncRNA
  }
  if (species %in% chicken_options) {
    lncRNA <- lncRNA_gene_list$Gallus_gallus_lncRNA
  }

  return(lncRNA)
}


#' Retrieve dual species gene lists mitochondrial
#'
#' Returns vector of all mitochondrial genes across all species in dataset.
#'
#' @param seurat_object object name.
#' @param species Species of origin for given Seurat Object.  Only accepted species are: mouse, human,
#' zebrafish, rat, drosophila, rhesus macaque, or chicken (name or abbreviation).
#' @param species_prefix the species prefix in front of gene symbols in object.
#' @param assay Assay to use (default is the current object default assay).
#'
#' @return vector of gene ids
#'
#' @keywords internal
#'
#' @noRd
#'

Retrieve_Dual_Mito_Features <- function(
    object,
    species,
    species_prefix,
    assay
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

  # Species Spelling Options
  mouse_options <- accepted_names$Mouse_Options
  human_options <- accepted_names$Human_Options
  marmoset_options <- accepted_names$Marmoset_Options
  zebrafish_options <- accepted_names$Zebrafish_Options
  rat_options <- accepted_names$Rat_Options
  drosophila_options <- accepted_names$Drosophila_Options
  macaque_options <- accepted_names$Macaque_Options
  chicken_options <- accepted_names$Chicken_Options


  mito_features_list <- lapply(1:length(x = species), function(x) {
    if (species[x] %in% mouse_options) {
      mito_pattern <- "mt-"
    }
    if (species[x] %in% human_options) {
      mito_pattern <- "MT-"
    }
    if (species[x] %in% c(marmoset_options, macaque_options, chicken_options)) {
      mito_features <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6")
    }
    if (species[x] %in% zebrafish_options) {
      mito_pattern <- "mt-"
    }
    if (species[x] %in% rat_options) {
      mito_pattern <- "Mt-"
    }
    if (species[x] %in% drosophila_options) {
      mito_pattern <- "mt:"
    }

    mito_pattern <- paste0("^", species_prefix[x], mito_pattern)

    mito_features <- grep(pattern = mito_pattern, x = rownames(x = object[[assay]]), value = TRUE)

    # Check features are present in object
    length_mito_features <- length(x = intersect(x = mito_features, y = rownames(x = object[[assay]])))

    # Check length of mito and ribo features found in object
    if (length_mito_features < 1) {
      cli_warn(message = c("No Mito features found in object using pattern/feature list provided.",
                           "i" = "No column will be added to meta.data.")
      )
    }
    mito_features
  })

  # combine the lists
  full_list <- unlist(x = mito_features_list)

  return(full_list)
}


#' Retrieve dual species gene lists ribosomal
#'
#' Returns vector of all ribosomal genes across all species in dataset.
#'
#' @param seurat_object object name.
#' @param species Species of origin for given Seurat Object.  Only accepted species are: mouse, human,
#' zebrafish, rat, drosophila, rhesus macaque, or chicken (name or abbreviation).
#' @param species_prefix the species prefix in front of gene symbols in object.
#' @param assay Assay to use (default is the current object default assay).
#'
#' @return vector of gene ids
#'
#' @keywords internal
#'
#' @noRd
#'

Retrieve_Dual_Ribo_Features <- function(
    object,
    species,
    species_prefix,
    assay = NULL
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

  # Species Spelling Options
  mouse_options <- accepted_names$Mouse_Options
  human_options <- accepted_names$Human_Options
  marmoset_options <- accepted_names$Marmoset_Options
  zebrafish_options <- accepted_names$Zebrafish_Options
  rat_options <- accepted_names$Rat_Options
  drosophila_options <- accepted_names$Drosophila_Options
  macaque_options <- accepted_names$Macaque_Options
  chicken_options <- accepted_names$Chicken_Options


  ribo_features_list <- lapply(1:length(x = species), function(x) {
    if (species[x] %in% mouse_options) {
      ribo_pattern <- "Rp[sl]"
    }
    if (species[x] %in% human_options) {
      ribo_pattern <- "RP[SL]"
    }
    if (species[x] %in% c(marmoset_options, macaque_options, chicken_options)) {
      ribo_pattern <- "RP[SL]"
    }
    if (species[x] %in% zebrafish_options) {
      ribo_pattern <- "rp[sl]"
    }
    if (species[x] %in% rat_options) {
      ribo_pattern <- "Rp[sl]"
    }
    if (species[x] %in% drosophila_options) {
      ribo_pattern <- "Rp[SL]"
    }

    ribo_pattern <- paste0("^", species_prefix[x], ribo_pattern)

    ribo_features <- grep(pattern = ribo_pattern, x = rownames(x = object[[assay]]), value = TRUE)

    # Check features are present in object
    length_ribo_features <- length(x = intersect(x = ribo_features, y = rownames(x = object[[assay]])))

    # Check length of mito and ribo features found in object
    if (length_ribo_features < 1) {
      cli_warn(message = c("No Ribo features found in object using pattern/feature list provided.",
                           "i" = "No column will be added to meta.data.")
      )
    }

    ribo_features
  })

  # combine the lists
  full_list <- unlist(x = ribo_features_list)

  return(full_list)
}


#' Retrieve MALAT1 Ensembl Gene IDs
#'
#' Retrieves species specific MALAT1 gene ensembl IDs
#'
#' @param species species to retrieve IDs.
#'
#' @return list of 2 sets of ensembl IDs
#'
#' @keywords internal
#'
#' @noRd
#'

Retrieve_MALAT1_Ensembl_Lists <- function(
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
    cli_abort(message = "Rat, Marmoset, Macaque, Zebrafish, Drosophila, and Chicken are not currently supported.")
  }

  # set prefix
  if (species %in% mouse_options) {
    prefix <- "Mus_musculus_"
  }
  if (species %in% human_options) {
    prefix <- "Homo_sapiens_"
  }

  # set list names
  malat <- paste0(prefix, "MALAT1_ensembl")

  # pull lists
  qc_gene_list <- list(
    malat1 = ensembl_malat1_list[[malat]]
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
#' zebrafish, rat, drosophila, rhesus macaque, or chicken (name or abbreviation)
#' @param oxphos_name name to use for the new meta.data column containing percent MSigDB Hallmark oxidative
#' phosphorylation counts. Default is "percent_oxphos".
#' @param apop_name name to use for the new meta.data column containing percent MSigDB Hallmark apoptosis counts.
#' Default is "percent_apop".
#' @param dna_repair_name name to use for the new meta.data column containing percent MSigDB Hallmark DNA repair counts.
#' Default is "percent_oxphos".
#' @param ensembl_ids logical, whether feature names in the object are gene names or
#' ensembl IDs (default is FALSE; set TRUE if feature names are ensembl IDs).
#' @param assay Assay to use (default is the current object default assay).
#' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
#' function will abort if columns with any one of the names provided to `mito_name` `ribo_name` or
#' `mito_ribo_name` is present in meta.data slot.
#'
#' @return Seurat object
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
    ensembl_ids = FALSE,
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
    Macaque_Options = c("Macaque", "macaque", "Rhesus", "macaca", "mmulatta", NA),
    Chicken_Options = c("Chicken", "chicken", "Gallus", "gallus", "Gg", "gg")
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
  if (isFALSE(x = ensembl_ids)) {
    msigdb_gene_list <- Retrieve_MSigDB_Lists(species = species)
  } else {
    msigdb_gene_list <- Retrieve_MSigDB_Ensembl_Lists(species = species)
  }


  oxphos_found <- Feature_PreCheck(object = seurat_object, features = msigdb_gene_list[["oxphos"]])
  apop_found <- Feature_PreCheck(object = seurat_object, features = msigdb_gene_list[["apop"]])
  dna_repair_found <- Feature_PreCheck(object = seurat_object, features = msigdb_gene_list[["dna_repair"]])

  # Add meta data columns
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
#' @param ensembl_ids logical, whether feature names in the object are gene names or
#' ensembl IDs (default is FALSE; set TRUE if feature names are ensembl IDs).
#' @param assay Assay to use (default is the current object default assay).
#' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
#' function will abort if columns with the name provided to `ieg_name` is present in meta.data slot.
#'
#' @return Seurat object
#'
#' @keywords internal
#'
#' @noRd
#'

Add_IEG_Seurat <- function(
    seurat_object,
    species,
    ieg_name = "percent_ieg",
    ieg_module_score = TRUE,
    ieg_module_name = "ieg_score",
    ensembl_ids = FALSE,
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
    Macaque_Options = c("Macaque", "macaque", "Rhesus", "macaca", "mmulatta", NA),
    Chicken_Options = c("Chicken", "chicken", "Gallus", "gallus", "Gg", "gg")
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
  if (isFALSE(x = ensembl_ids)) {
    ieg_gene_list <- Retrieve_IEG_Lists(species = species)
  } else {
    ieg_gene_list <- Retrieve_IEG_Ensembl_Lists(species = species)
  }

  ieg_found <- Feature_PreCheck(object = seurat_object, features = ieg_gene_list[["ieg"]])

  # Add mito and ribo columns
  if (length(x = ieg_found) > 0) {
    seurat_object[[ieg_name]] <- PercentageFeatureSet(object = seurat_object, features = ieg_found, assay = assay)

    if (isTRUE(x = ieg_module_score)) {
      # Overwrite check
      if (ieg_module_name %in% colnames(x = seurat_object@meta.data)) {
        if (isFALSE(x = overwrite)) {
          cli_abort(message = c("Column with {.val {ieg_module_name}} already present in meta.data slot.",
                                "i" = "*To run function and overwrite column set parameter {.code overwrite = TRUE} or change respective {.code ieg_module_name}*")
          )
        }
        cli_inform(message = c("Column with {.val {ieg_module_name}} already present in meta.data slot.",
                               "i" = "Overwriting those column as {.code overwrite = TRUE.}")
        )
      }

      # check normalized
      if (isFALSE(x = Check_Normalized(object = seurat_object, assay = assay, error = FALSE))) {
        cli_inform(message = c("x" = "Layer with normalized data not present.",
                               "i" = "Normalizing Data."))
        seurat_object <- NormalizeData(object = seurat_object)
      }

      seurat_object <- AddModuleScore(object = seurat_object, features = list(ieg_found), search = FALSE, name = ieg_module_name)
    }
  }

  # Log Command
  seurat_object <- LogSeuratCommand(object = seurat_object)

  # return final object
  return(seurat_object)
}


#' Add lncRNA Gene List Percentages
#'
#' Adds percentage of counts from lncRNA genes from mouse and human.
#'
#' @param seurat_object object name.
#' @param species Species of origin for given Seurat Object.
#' @param lncRNA_name name to use for the new meta.data column containing percent lncRNA gene counts. Default is "percent_lncRNA".
#' @param assay Assay to use (default is the current object default assay).
#' @param overwrite Logical.  Whether to overwrite existing meta.data columns.  Default is FALSE meaning that
#' function will abort if columns with the name provided to `lncRNA_name` is present in meta.data slot.
#'
#' @return Seurat object
#'
#' @keywords internal
#'
#' @noRd
#'

Add_lncRNA_Seurat <- function(
    seurat_object,
    species,
    lncRNA_name = "percent_lncRNA",
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
    Macaque_Options = c("Macaque", "macaque", "Rhesus", "macaca", "mmulatta", NA),
    Chicken_Options = c("Chicken", "chicken", "Gallus", "gallus", "Gg", "gg")
  )

  if (!species %in% unlist(x = accepted_names)) {
    cli_inform(message = "The supplied species ({.field {species}}) is not currently supported.")
  }

  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)

  # Overwrite check
  if (lncRNA_name %in% colnames(x = seurat_object@meta.data)) {
    if (isFALSE(x = overwrite)) {
      cli_abort(message = c("Column with {.val {lncRNA_name}} already present in meta.data slot.",
                            "i" = "*To run function and overwrite column set parameter {.code overwrite = TRUE} or change respective {.code lncRNA_name}*")
      )
    }
    cli_inform(message = c("Column with {.val {lncRNA_name}} already present in meta.data slot.",
                           "i" = "Overwriting those column as {.code overwrite = TRUE.}")
    )
  }

  # Set default assay
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Retrieve gene lists
  complete_lnc_list <- c(Retrieve_lncRNA(species = species), Retrieve_Ensembl_lncRNA(species = species))

  lncRNA_found <- Feature_PreCheck(object = seurat_object, features = complete_lnc_list, print_missing = FALSE)

  # Add mito and ribo columns
  if (length(x = lncRNA_found) > 0) {
    seurat_object[[lncRNA_name]] <- PercentageFeatureSet(object = seurat_object, features = lncRNA_found, assay = assay)
  }

  # Log Command
  seurat_object <- LogSeuratCommand(object = seurat_object)

  # return final object
  return(seurat_object)
}


#' Define MALAT1 QC Threshold
#'
#' This function is a threshold calculation and procedure as described in
#' Clarke & Bader (2024). bioRxiv \doi{10.1101/2024.07.14.603469}.  Please cite this preprint
#' whenever using this function.  The original code and function can be found here:
#' \url{https://github.com/BaderLab/MALAT1_threshold}
#'
#' @param counts A vector of normalized MALAT1 counts; this is ideally from a single unintegrated sample consisting
#' of multiple cell types
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
#' @param print_plots logical, whether to print the plots while running the function (default is TRUE).
#' @param return_plots logical, whether to return plots as part of the results of the function (default is FALSE).
#' @param plot_title title to add to plot patchwork, default is NULL.
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom cowplot theme_cowplot
#'
#' @noRd
#'
#' @references This function incorporates a threshold calculation and procedure as described in
#' Clarke & Bader (2024). bioRxiv \doi{10.1101/2024.07.14.603469}.  Please cite this preprint
#' whenever using this function.
#'
#' @author Zoe Clark
#'
#' @return threshold value to use in
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' object <- define_malat1_threshold(object = normalized_counts)
#' }
#'

define_malat1_threshold <- function(
    counts,
    bw = 0.1,
    lwd = 2,
    breaks = 100,
    chosen_min = 1,
    smooth = 1,
    abs_min = 0.3,
    rough_max = 2,
    print_plots = TRUE,
    return_plots = FALSE,
    plot_title = NULL
) {
  tryCatch({
    # Calculate the density values
    density_data <- density(counts, bw = bw,
                            from = min(counts),
                            to = max(counts))
    # Fit a smooth line to the data
    fit <- smooth.spline(density_data$x, density_data$y, spar = smooth)
    # Predict the first derivative
    first_derivative <- predict(fit, density_data$x, deriv = 1)

    # Find the local maxima
    local_maxima <- density_data$x[which(diff(sign(first_derivative$y)) == -2)]
    if(length(local_maxima) == 0) {
      # Find x-val closest to rough_max
      local_maxima <- density_data$x[which(abs(density_data$x - rough_max) == min(abs(density_data$x - rough_max)))]
    }

    # Create data frames for plotting
    density_df <- data.frame(x = density_data$x, y = density_data$y)
    fit_df <- data.frame(x = density_data$x, y = predict(fit, density_data$x)$y)
    local_maxima_df <- data.frame(x = local_maxima, y = predict(fit, local_maxima)$y)

    # Plot the density with the local maxima using ggplot2
    p1 <- ggplot() +
      geom_line(data = density_df, aes(x = .data[["x"]], y = .data[["y"]]), color = "black", na.rm = TRUE) +
      geom_line(data = fit_df, aes(x = .data[["x"]], y = .data[["y"]]), color = "blue", na.rm = TRUE) +
      geom_point(data = local_maxima_df, aes(x = .data[["x"]], y = .data[["y"]]), color = "red", size = 3, na.rm = TRUE) +
      ggtitle("Density Plot with Local Maxima") +
      theme_cowplot() +
      ylim(-0.05, NA)


    # Find the local minima
    local_minima <- density_data$x[which(diff(sign(first_derivative$y)) == 2)]
    if(length(local_minima) == 0) {
      local_minima <- abs_min
    }

    # Create data frame for local minima
    local_minima_df <- data.frame(x = local_minima, y = predict(fit, local_minima)$y)

    # Plot the density with the local minima using ggplot2
    p2 <- ggplot() +
      geom_line(data = density_df, aes(x = .data[["x"]], y = .data[["y"]]), color = "black", na.rm = TRUE) +
      geom_line(data = fit_df, aes(x = .data[["x"]], y = .data[["y"]]), color = "blue", na.rm = TRUE) +
      geom_point(data = local_minima_df, aes(x = .data[["x"]], y = .data[["y"]]), color = "red", size = 3, na.rm = TRUE) +
      ggtitle("Density Plot with Local Minima") +
      theme_cowplot() +
      ylim(-0.05, NA)

    # Find biggest local maximum greater than desired value (default 1)
    biggest_y <- max(density_data$y[density_data$x %in% local_maxima[local_maxima > chosen_min]])
    maxi <- density_data$x[density_data$y == biggest_y]

    # Find local minimum closest to the left of that
    local_minima <- local_minima[local_minima < maxi]
    if(length(local_minima) < abs_min) {
      local_minima <- abs_min
    }
    mini <- local_minima[(maxi - local_minima) == min(maxi - local_minima)]

    # Calculate delta to get range of values to isolate for quadratic calculation
    delta <- maxi - mini

    # Subset dataframe
    df <- as.data.frame(cbind(density_data$x, density_data$y))
    colnames(df) <- c("x","y")
    subset_df <- df[df[,1] >= (maxi - delta) & df[,1] <= (maxi + delta), ]

    # Fit a quadratic model (y ~ x + I(x^2))
    quad_model <- lm(y ~ poly(x, 2, raw = TRUE), data = subset_df)

    p3 <- ggplot(df, aes(x = .data[["x"]], y = .data[["y"]])) +
      geom_line(color = "black") +
      geom_point(data = subset_df, aes(x = .data[["x"]], y = .data[["y"]]), color = "blue", size = 2, na.rm = TRUE) +
      stat_function(fun = function(x) predict(quad_model, newdata = data.frame(x = x)), color = "red", linewidth = 1, na.rm = TRUE) +
      labs(
        title = "Density Plot with Quadratic Fit",
        x = "Normalized MALAT1 expression",
        y = "Density value"
      ) +
      theme_cowplot() +
      ylim(-0.05, NA)

    # Grab intercept
    # Extract coefficients from the summary
    coefficients <- summary(quad_model)$coefficients
    # Extract coefficients
    c <- coefficients[1, 1]
    b <- coefficients[2, 1]
    a <- coefficients[3, 1]
    # Calculate discriminant
    discriminant <- b^2 - 4 * a * c
    # Grab first intercept
    x_intercept1 <- (-b + sqrt(discriminant)) / (2 * a)
    if(x_intercept1 < 0) {
      x_intercept1 <- abs_min
    }

    p4 <- ggplot(data = data.frame(counts), aes(x = .data[["counts"]])) +
      geom_histogram(color = "black", bins = breaks, fill = "dodgerblue", na.rm = TRUE) +
      theme_cowplot() +
      geom_vline(xintercept = x_intercept1, color = "red", linewidth = lwd) +
      ggtitle("MALAT1")

    plots <- wrap_plots(p1, p2, p3, p4, ncol = 2) + plot_annotation(title = plot_title, theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = rel(2))), tag_levels = "A")

    if (isTRUE(x = print_plots)) {
      print(plots)
    }

    if (isTRUE(x = return_plots)) {
      res <- list("threshold" = x_intercept1,
                  "plots" = plots)
      return(res)
    } else {
      return(x_intercept1)
    }

  }, error = function(e) {
    # Code to execute if an error occurs
    ggplot(data = data.frame(counts), aes(x = .data[["counts"]])) +
      geom_histogram(color = "black", bins = breaks, fill = "dodgerblue") +
      theme_cowplot() +
      ggtitle("MALAT1")

    stop(" An error occurred: Please make sure you have use a vector of normalized counts as input. This may also indicate that you have no high MALAT1 peaks, meaning this particular sample may be poor quality. Please check your histogram of normalized MALAT1 counts to investigate: ", e$message)

  })
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
