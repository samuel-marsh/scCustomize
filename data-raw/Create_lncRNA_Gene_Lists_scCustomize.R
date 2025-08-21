# Code and functions to create the lists of ensembl IDs and gene symbols for lncRNA


# Functions -----------------------------------------------------------------------------------
library(tidyverse)
library(AnnotationHub)

Create_Ensembl_lncRNA_List <- function(
) {

  refreshHub(hubClass="AnnotationHub")

  species_list <- c("Mus musculus", "Homo sapiens", "Callithrix jacchus", "Danio rerio", "Rattus norvegicus", "Macaca mulatta", "Gallus gallus")
      # "Drosophila melanogaster" doesn't have separate lncRNA biotype, only ncRNA.  Therefore not adding to list.

  ah <- AnnotationHub()

  lncRNA_ensembl_list <- lapply(1:length(x = species_list), function(x){

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

    if (isTRUE(species_list[x] %in% c("Mus musculus", "Homo sapiens", "Callithrix jacchus", "Rattus norvegicus", "Macaca mulatta", "Gallus gallus"))) {
      lncRNA_ensembl_ids <- annotations %>%
        dplyr::filter(gene_biotype == "lncRNA") %>%
        dplyr::pull(gene_id)
    }

    if (isTRUE(species_list[x] == "Danio rerio")) {
      lncRNA_ensembl_ids <- annotations %>%
        dplyr::filter(gene_biotype == "lincRNA") %>%
        dplyr::pull(gene_id)
    }

    cli::cli_alert("Complete")
    return(lncRNA_ensembl_ids)

  })

  names(lncRNA_ensembl_list) <- paste0(gsub(pattern = " ", replacement = "_", x = species_list), "_lncRNA_ensembl")

  return(lncRNA_ensembl_list)
}


Create_lncRNA_List <- function(
) {

  refreshHub(hubClass="AnnotationHub")

  species_list <- c("Mus musculus", "Homo sapiens", "Danio rerio", "Rattus norvegicus", "Macaca mulatta", "Gallus gallus")
      # "Callithrix jacchus" do not have any annotated symbol/names for lncRNA and therefore excluded.
      # "Drosophila melanogaster" doesn't have separate lncRNA biotype, only ncRNA.  Therefore not adding to list.

  ah <- AnnotationHub()

  lncRNA_list <- lapply(1:length(x = species_list), function(x){

    cli::cli_inform("Retrieving gene ID for {species_list[x]}")
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
      dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid, symbol)

    if (isTRUE(species_list[x] %in% c("Mus musculus", "Homo sapiens", "Rattus norvegicus", "Macaca mulatta", "Gallus gallus"))) {
      lncRNA_ids <- annotations %>%
        dplyr::filter(gene_biotype == "lncRNA") %>%
        dplyr::pull(gene_name) %>%
        stringi::stri_remove_empty_na() %>%
        unique()
    }

    if (isTRUE(species_list[x] == "Danio rerio")) {
      lncRNA_ids <- annotations %>%
        dplyr::filter(gene_biotype == "lincRNA") %>%
        dplyr::pull(gene_name) %>%
        stringi::stri_remove_empty_na() %>%
        unique()
    }


    cli::cli_alert("Complete")
    return(lncRNA_ids)

  })

  names(lncRNA_list) <- paste0(gsub(pattern = " ", replacement = "_", x = species_list), "_lncRNA")

  return(lncRNA_list)
}



# Create & Save Lists -------------------------------------------------------------------------
ensembl_lncRNA_id <- Create_Ensembl_lncRNA_List()

save(ensembl_lncRNA_id, file = "data/ensembl_lncRNA_id.rda")

lncRNA_gene_list <- Create_lncRNA_List()

save(lncRNA_gene_list, file = "data/lncRNA_gene_list.rda")

