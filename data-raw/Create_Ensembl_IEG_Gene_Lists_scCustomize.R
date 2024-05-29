# Scripts to create ensembl_id gene lists for IEGs

# Mouse gene list is from: SI Table 4 from \doi{10.1016/j.neuron.2017.09.026}.  Human
# gene list was compiled by first creating homologous gene list using biomaRt and then adding some manually curated
# homologs according to HGNC.  See data-raw directory for scripts used to create gene list.

library(dplyr)
library(AnnotationHub)


Create_Ensembl_IEG_List <- function(
) {

  refreshHub(hubClass="AnnotationHub")

  species_list <- c("Mus musculus", "Homo sapiens")
  ieg_symbol_list <- scCustomize:::ieg_gene_list

  ah <- AnnotationHub()

  ieg_list <- lapply(1:length(x = species_list), function(x){

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

    ieg_ids <- annotations %>%
      dplyr::filter(gene_name %in% ieg_symbol_list[[x]], gene_biotype != "LRG_gene") %>%
      dplyr::pull(gene_id)


    cli::cli_alert("Complete")
    return(ieg_ids)

  })

  names(ieg_list) <- paste0(gsub(pattern = " ", replacement = "_", x = species_list), "_IEG_ensembl")

  return(ieg_list)
}

ensembl_ieg_list <- Create_Ensembl_IEG_List()

save(ensembl_ieg_list, file = "data/ensembl_ieg_list.rda")

