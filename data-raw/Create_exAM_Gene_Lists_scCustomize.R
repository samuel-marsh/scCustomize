# Scripts to create symbol and ensembl_id gene lists for exAM

# Symbol list is from: SI Table 22 from Marsh et al., 2022 \doi{10.1038/s41593-022-01022-8}.


# Gene List & Symbol Updates ------------------------------------------------------------------

# Changed Hist2h2aa4 to H2ac19 as new symbol was not updating for mice (works for humans)
# Changed C2orf40 to Ecrg4 as new symbol was not updating for mice (works for humans)
exAM_mouse_union <- c("Rgs1", "H2ac19", "Junb", "Dusp1", "Hspa1a", "Hsp90aa1", "Fos", "Hspa1b", "Jun", "Zfp36", "Atf3", "Btg2", "Fosb", "Hist1h1d", "Hist1h2ac", "Hist1h4i", "Nfkbiz", "Klf2", "Ccl3", "Jund", "Nfkbid", "Gem", "Ccl4", "Ier5", "Txnip", "Hist1h2bc", "Hist1h1c", "Egr1", "Rhob", "Socs3", "Ecrg4", "Hist1h1e", "Folr1", "Serpine1")

exAM_mouse_union_update <- Updated_MGI_Symbols(input_data = exAM_mouse_union)
exAM_mouse_union_update <- exAM_mouse_union_update$Output_Features


exAM_human_union <- c("RGS1", "HIST2H2AA4", "JUNB", "DUSP1", "HSPA1A", "HSP90AA1", "FOS", "HSPA1B", "JUN", "ZFP36", "ATF3", "BTG2", "FOSB", "HIST1H1D", "HIST1H2AC", "HIST1H4I", "NFKBIZ", "KLF2", "CCL3", "JUND", "NFKBID", "GEM", "CCL4", "IER5", "TXNIP", "HIST1H2BC", "HIST1H1C", "EGR1", "RHOB", "SOCS3", "C2orf40", "HIST1H1E", "FOLR1", "SERPINE1")

exAM_human_union_update <- Updated_HGNC_Symbols(input_data = exAM_human_union)
exAM_human_union_update <- exAM_human_union_update$Output_Features


exAM_human_microglia_pm <- c("HSPA1A", "HSPA1B", "FOS", "HIST1H2AC", "UBC", "HSP90AA1", "DUSP1", "HIST1H2BG", "RGS1", "HIST1H2BF", "SLC2A3", "HIST2H2BE", "DDIT4", "HIST1H2BJ", "HSPB1", "HIST1H2BD", "HIST1H1D", "SRGN", "GLUL", "BCAS2", "HIST2H2AA4", "HSPH1", "DNAJB1", "BTG2", "HIST1H2BH", "ZFP36L1", "PPP1R15A", "HSPE1", "JUN", "GADD45B", "P4HA1", "HMOX1", "RHOB", "HIST1H3H", "ZFP36", "RGS16", "HIST2H2BF", "GPR183")

exAM_human_microglia_pm_update <- Updated_HGNC_Symbols(input_data = exAM_human_microglia_pm)
exAM_human_microglia_pm_update <- exAM_human_microglia_pm_update$Output_Features



exAM_gene_list <- list("Mus_musculus_exAM_union" = exAM_mouse_union_update,
                       "Homo_sapiens_exAM_union" = exAM_human_union_update,
                       "Homo_sapiens_exAM_micro" = exAM_human_microglia_pm_update)
save(exAM_gene_list, file = "data/exAM_gene_list.rda")


# Ensembl IDs ---------------------------------------------------------------------------------


library(dplyr)
library(AnnotationHub)


Create_Ensembl_exAM_union_List <- function(
) {

  refreshHub(hubClass="AnnotationHub")

  species_list <- c("Mus musculus", "Homo sapiens")
  exAM_union_symbol_list <- list("Mus_musculus_exAM" = exAM_mouse_union_update,
                                 "Homo_sapiens_exAM" = exAM_human_union_update)

  ah <- AnnotationHub()

  exAM_union_list <- lapply(1:length(x = species_list), function(x){

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
      dplyr::filter(gene_name %in% exAM_union_symbol_list[[x]], gene_biotype != "LRG_gene") %>%
      dplyr::pull(gene_id)


    cli::cli_alert("Complete")
    return(ieg_ids)

  })

  names(exAM_union_list) <- paste0(gsub(pattern = " ", replacement = "_", x = species_list), "_exAM_union_ensembl")

  return(exAM_union_list)
}

ensembl_exAM_union_list <- Create_Ensembl_exAM_union_List()



Create_Ensembl_exAM_micro_List <- function(
) {

  refreshHub(hubClass="AnnotationHub")

  species_list <- c("Homo sapiens")
  exAM_union_symbol_list <- list("Homo_sapiens_exAM_micro" = exAM_human_microglia_pm_update)

  ah <- AnnotationHub()

  exAM_union_list <- lapply(1:length(x = species_list), function(x){

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
      dplyr::filter(gene_name %in% exAM_union_symbol_list[[x]], gene_biotype != "LRG_gene") %>%
      dplyr::pull(gene_id)


    cli::cli_alert("Complete")
    return(ieg_ids)

  })

  names(exAM_union_list) <- paste0(gsub(pattern = " ", replacement = "_", x = species_list), "_exAM_micro_ensembl")

  return(exAM_union_list)
}

ensembl_exAM_micro_list <- Create_Ensembl_exAM_micro_List()

ensembl_exAM_list <- c(ensembl_exAM_union_list, ensembl_exAM_micro_list)



save(ensembl_exAM_list, file = "data/ensembl_exAM_list.rda")

