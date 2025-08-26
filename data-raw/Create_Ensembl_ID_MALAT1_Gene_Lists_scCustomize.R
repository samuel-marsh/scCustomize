# Scripts to create MALAT1 ensembl ID gene list

library(dplyr)
library(AnnotationHub)


Create_Ensembl_MALAT1_List <- function(
) {

  refreshHub(hubClass="AnnotationHub")

  species_list <- c("Mus musculus", "Homo sapiens")
  malat_symbol <- c("Malat1", "MALAT1")

  ah <- AnnotationHub()

  malat_list <- lapply(1:length(x = species_list), function(x){

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

    malat_ids <- annotations %>%
      dplyr::filter(gene_name %in% malat_symbol, gene_biotype != "LRG_gene") %>%
      dplyr::pull(gene_id)


    cli::cli_alert("Complete")
    return(malat_ids)

  })

  names(malat_list) <- paste0(gsub(pattern = " ", replacement = "_", x = species_list), "_MALAT1_ensembl")

  return(malat_list)
}

ensembl_malat1_list <- Create_Ensembl_MALAT1_List()

save(ensembl_malat1_list, file = "data/ensembl_malat1_list.rda")
