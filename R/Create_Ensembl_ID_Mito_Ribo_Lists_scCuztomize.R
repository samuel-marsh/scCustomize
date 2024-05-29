# Code and functions to create the lists of ensembl IDs for mitochondrial and ribosomal genes


# Functions -----------------------------------------------------------------------------------
library(tidyverse)
library(AnnotationHub)

Create_Ensembl_Ribo_List <- function(
) {

  refreshHub(hubClass="AnnotationHub")

  species_list <- c("Mus musculus", "Homo sapiens", "Callithrix jacchus", "Danio rerio", "Rattus norvegicus", "Drosophila melanogaster", "Macaca mulatta", "Gallus gallus")
  ribo_pattern_list <- c("^Rp[sl]", "^RP[SL]", "^RP[SL]", "^rp[sl]", "^Rp[sl]", "^Rp[SL]", "^RP[SL]", "^RP[SL]")

  ah <- AnnotationHub()

  ribo_list <- lapply(1:length(x = species_list), function(x){

    cli::cli_inform("Retrieving ensembl ID for {species_list[x]}")
    # Access the Ensembl database for organism
    ahDb <- query(ah,
                  pattern = c(species_list[x], "EnsDb"),
                  ignore.case = TRUE)

    # Check versions of databases available
    ahDb %>%
      mcols()


    # Acquire the latest annotation files
    id <- ahDb %>%
      mcols() %>%
      rownames() %>%
      tail(n = 1)

    # Download the appropriate Ensembldb database
    edb <- ah[[id]]


    # Extract gene-level information from database
    annotations <- genes(edb,
                         return.type = "data.frame")


    # Select annotations of interest
    annotations <- annotations %>%
      dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)

    ribo_ids <- annotations %>%
      dplyr::filter(str_detect(string = gene_name, pattern = ribo_pattern_list[x]), gene_biotype != "LRG_gene") %>%
      dplyr::pull(gene_id)


    cli::cli_alert("Complete")
    return(ribo_ids)

  })

  names(ribo_list) <- paste0(gsub(pattern = " ", replacement = "_", x = species_list), "_ribo_ensembl")

  return(ribo_list)
}


Create_Ensembl_Mito_List <- function(
) {

  refreshHub(hubClass="AnnotationHub")

  species_list <- c("Mus musculus", "Homo sapiens", "Danio rerio", "Rattus norvegicus", "Drosophila melanogaster", "Macaca mulatta", "Gallus gallus")

  ah <- AnnotationHub()

  mito_list <- lapply(species_list, function(x){

    cli::cli_inform("Retrieving ensembl ID for {x}")
    # Access the Ensembl database for organism
    ahDb <- query(ah,
                  pattern = c(x, "EnsDb"),
                  ignore.case = TRUE)

    # Check versions of databases available
    ahDb %>%
      mcols()


    # Acquire the latest annotation files
    id <- ahDb %>%
      mcols() %>%
      rownames() %>%
      tail(n = 1)

    # Download the appropriate Ensembldb database
    edb <- ah[[id]]


    # Extract gene-level information from database
    annotations <- genes(edb,
                         return.type = "data.frame")


    # Select annotations of interest
    annotations <- annotations %>%
      dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)

    if (x == "Drosophila melanogaster") {
      mito_ids <- annotations %>%
        dplyr::filter(seq_name == "mitochondrion_genome") %>%
        dplyr::pull(gene_id)
    } else {
      mito_ids <- annotations %>%
        dplyr::filter(seq_name == "MT") %>%
        dplyr::pull(gene_id)
    }

    cli::cli_alert("Complete")
    return(mito_ids)


  })

  names(mito_list) <- paste0(gsub(pattern = " ", replacement = "_", x = species_list), "_mito_ensembl")

  return(mito_list)
}


# Create & Save Lists -------------------------------------------------------------------------
ensembl_mito_id <- Create_Ensembl_Mito_List()

save(ensembl_mito_id, file = "data/ensembl_mito_id.rda")

ensembl_ribo_id <- Create_Ensembl_Ribo_List()

save(ensembl_ribo_id, file = "data/ensembl_ribo_id.rda")



# Deprecated Code -----------------------------------------------------------------------------
# OLD MARMOSET VERSION
# library(biomaRt)
# mart <- useMart(biomart = "ensembl", dataset = "cjacchus_gene_ensembl")
# listDatasets(mart = mart)
# attributes <- c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand", "gene_biotype", "description")
#
# marmoset_ensembl_ids <- getBM(attributes = attributes, mart = mart)
#
# dplyr::filter(str_detect(string = gene_name, pattern = ribo_pattern_list[2]), gene_biotype != "LRG_gene")
#
# marmoset_ribo_ids <- marmoset_ensembl_ids %>%
#   filter(str_detect(string = external_gene_name, pattern = "^RP[SL]"), gene_biotype != "LRG_gene") %>%
#   dplyr::pull(ensembl_gene_id)
#
#
# ensembl_ribo_id[["Callithrix_jacchus_ribo_ensembl"]] <- marmoset_ribo_ids



