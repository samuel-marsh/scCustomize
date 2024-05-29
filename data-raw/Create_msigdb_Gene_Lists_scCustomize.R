# Instructions for creation of package msigdb gene lists


# Create Gene Symbol Lists --------------------------------------------------------------------
library(dplyr)
library(msigdbr)

msigdbr_species()

msig_dbr <- msigdbr(species = "Homo sapiens", category = "H")

msig_oxphos_direct <- msig_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  pull(gene_symbol) %>%
  unique()


msig_apoptosis_direct <- msig_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(gene_symbol) %>%
  unique()


msig_DNA_repair_direct <- msig_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_DNA_REPAIR") %>%
  pull(gene_symbol) %>%
  unique()




msig_mouse_dbr <- msigdbr(species = "Mus musculus", category = "H")

msig_mouse_oxphos_direct <- msig_mouse_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  pull(gene_symbol) %>%
  unique()


msig_mouse_apoptosis_direct <- msig_mouse_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(gene_symbol) %>%
  unique()


msig_mouse_DNA_repair_direct <- msig_mouse_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_DNA_REPAIR") %>%
  pull(gene_symbol) %>%
  unique()


# zebrafish
msig_zebra_dbr <- msigdbr(species = "Danio rerio", category = "H")

msig_zebra_oxphos_direct <- msig_zebra_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  pull(gene_symbol) %>%
  unique()


msig_zebra_apoptosis_direct <- msig_zebra_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(gene_symbol) %>%
  unique()


msig_zebra_DNA_repair_direct <- msig_zebra_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_DNA_REPAIR") %>%
  pull(gene_symbol) %>%
  unique()


# rat
msig_rat_dbr <- msigdbr(species = "Rattus norvegicus", category = "H")

msig_rat_oxphos_direct <- msig_rat_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  pull(gene_symbol) %>%
  unique()


msig_rat_apoptosis_direct <- msig_rat_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(gene_symbol) %>%
  unique()


msig_rat_DNA_repair_direct <- msig_rat_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_DNA_REPAIR") %>%
  pull(gene_symbol) %>%
  unique()



# fly
msig_fly_dbr <- msigdbr(species = "Drosophila melanogaster", category = "H")

msig_fly_oxphos_direct <- msig_fly_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  pull(gene_symbol) %>%
  unique()


msig_fly_apoptosis_direct <- msig_fly_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(gene_symbol) %>%
  unique()


msig_fly_DNA_repair_direct <- msig_fly_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_DNA_REPAIR") %>%
  pull(gene_symbol) %>%
  unique()


# macaque
msig_macaque_dbr <- msigdbr(species = "Macaca mulatta", category = "H")

msig_macaque_oxphos_direct <- msig_macaque_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  pull(gene_symbol) %>%
  unique()


msig_macaque_apoptosis_direct <- msig_macaque_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(gene_symbol) %>%
  unique()


msig_macaque_DNA_repair_direct <- msig_macaque_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_DNA_REPAIR") %>%
  pull(gene_symbol) %>%
  unique()


msig_chicken_dbr <- msigdbr(species = "Gallus gallus", category = "H")

msig_chicken_oxphos_direct <- msig_chicken_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  pull(gene_symbol) %>%
  unique()


msig_chicken_apoptosis_direct <- msig_chicken_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(gene_symbol) %>%
  unique()


msig_chicken_DNA_repair_direct <- msig_chicken_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_DNA_REPAIR") %>%
  pull(gene_symbol) %>%
  unique()



msigdb_qc_gene_list <- list(
  Homo_sapiens_msigdb_oxphos = msig_oxphos_direct,
  Homo_sapiens_msigdb_apop = msig_apoptosis_direct,
  Homo_sapiens_msigdb_dna_repair = msig_DNA_repair_direct,
  Mus_musculus_msigdb_oxphos = msig_mouse_oxphos_direct,
  Mus_musculus_msigdb_apop = msig_mouse_apoptosis_direct,
  Mus_musculus_msigdb_dna_repair = msig_mouse_DNA_repair_direct,
  Rattus_norvegicus_msigdb_oxphos = msig_rat_oxphos_direct,
  Rattus_norvegicus_msigdb_apop = msig_rat_apoptosis_direct,
  Rattus_norvegicus_msigdb_dna_repair = msig_rat_DNA_repair_direct,
  Drosophila_melanogaster_msigdb_oxphos = msig_fly_oxphos_direct,
  Drosophila_melanogaster_msigdb_apop = msig_fly_apoptosis_direct,
  Drosophila_melanogaster_msigdb_dna_repair = msig_fly_DNA_repair_direct,
  Dario_rerio_msigdb_oxphos = msig_zebra_oxphos_direct,
  Dario_rerio_msigdb_apop = msig_zebra_apoptosis_direct,
  Dario_rerio_msigdb_dna_repair = msig_zebra_DNA_repair_direct,
  Macaca_mulatta_msigdb_oxphos = msig_macaque_oxphos_direct,
  Macaca_mulatta_msigdb_apop = msig_macaque_apoptosis_direct,
  Macaca_mulatta_msigdb_dna_repair = msig_macaque_DNA_repair_direct,
  Gallus_gallus_msigdb_oxphos = msig_chicken_oxphos_direct,
  Gallus_gallus_msigdb_apop = msig_chicken_apoptosis_direct,
  Gallus_gallus_msigdb_dna_repair = msig_chicken_DNA_repair_direct
)

save(msigdb_qc_gene_list, file = "data/msigdb_qc_gene_list.rda")

# Create Ensembl ID Lists --------------------------------------------------------------------
library(dplyr)
library(msigdbr)

msigdbr_species()

msig_dbr <- msigdbr(species = "Homo sapiens", category = "H")

msig_oxphos_direct_ensembl <- msig_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  pull(ensembl_gene) %>%
  unique()


msig_apoptosis_direct_ensembl <- msig_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(ensembl_gene) %>%
  unique()


msig_DNA_repair_direct_ensembl <- msig_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_DNA_REPAIR") %>%
  pull(ensembl_gene) %>%
  unique()




msig_mouse_dbr <- msigdbr(species = "Mus musculus", category = "H")

msig_mouse_oxphos_direct_ensembl <- msig_mouse_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  pull(ensembl_gene) %>%
  unique()


msig_mouse_apoptosis_direct_ensembl <- msig_mouse_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(ensembl_gene) %>%
  unique()


msig_mouse_DNA_repair_direct_ensembl <- msig_mouse_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_DNA_REPAIR") %>%
  pull(ensembl_gene) %>%
  unique()


# zebrafish
msig_zebra_dbr <- msigdbr(species = "Danio rerio", category = "H")

msig_zebra_oxphos_direct_ensembl <- msig_zebra_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  pull(ensembl_gene) %>%
  unique()


msig_zebra_apoptosis_direct_ensembl <- msig_zebra_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(ensembl_gene) %>%
  unique()


msig_zebra_DNA_repair_direct_ensembl <- msig_zebra_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_DNA_REPAIR") %>%
  pull(ensembl_gene) %>%
  unique()


# rat
msig_rat_dbr <- msigdbr(species = "Rattus norvegicus", category = "H")

msig_rat_oxphos_direct_ensembl <- msig_rat_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  pull(ensembl_gene) %>%
  unique()


msig_rat_apoptosis_direct_ensembl <- msig_rat_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(ensembl_gene) %>%
  unique()


msig_rat_DNA_repair_direct_ensembl <- msig_rat_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_DNA_REPAIR") %>%
  pull(ensembl_gene) %>%
  unique()



# fly
msig_fly_dbr <- msigdbr(species = "Drosophila melanogaster", category = "H")

msig_fly_oxphos_direct_ensembl <- msig_fly_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  pull(ensembl_gene) %>%
  unique()


msig_fly_apoptosis_direct_ensembl <- msig_fly_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(ensembl_gene) %>%
  unique()


msig_fly_DNA_repair_direct_ensembl <- msig_fly_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_DNA_REPAIR") %>%
  pull(ensembl_gene) %>%
  unique()


# macaque
msig_macaque_dbr <- msigdbr(species = "Macaca mulatta", category = "H")

msig_macaque_oxphos_direct_ensembl <- msig_macaque_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  pull(ensembl_gene) %>%
  unique()


msig_macaque_apoptosis_direct_ensembl <- msig_macaque_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(ensembl_gene) %>%
  unique()


msig_macaque_DNA_repair_direct_ensembl <- msig_macaque_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_DNA_REPAIR") %>%
  pull(ensembl_gene) %>%
  unique()


msig_chicken_dbr <- msigdbr(species = "Gallus gallus", category = "H")

msig_chicken_oxphos_direct_ensembl <- msig_chicken_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>%
  pull(ensembl_gene) %>%
  unique()


msig_chicken_apoptosis_direct_ensembl <- msig_chicken_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_APOPTOSIS") %>%
  pull(ensembl_gene) %>%
  unique()


msig_chicken_DNA_repair_direct_ensembl <- msig_chicken_dbr %>%
  dplyr::filter(gs_name == "HALLMARK_DNA_REPAIR") %>%
  pull(ensembl_gene) %>%
  unique()



msigdb_qc_ensembl_list <- list(
  Homo_sapiens_msigdb_oxphos = msig_oxphos_direct_ensembl,
  Homo_sapiens_msigdb_apop = msig_apoptosis_direct_ensembl,
  Homo_sapiens_msigdb_dna_repair = msig_DNA_repair_direct_ensembl,
  Mus_musculus_msigdb_oxphos = msig_mouse_oxphos_direct_ensembl,
  Mus_musculus_msigdb_apop = msig_mouse_apoptosis_direct_ensembl,
  Mus_musculus_msigdb_dna_repair = msig_mouse_DNA_repair_direct_ensembl,
  Rattus_norvegicus_msigdb_oxphos = msig_rat_oxphos_direct_ensembl,
  Rattus_norvegicus_msigdb_apop = msig_rat_apoptosis_direct_ensembl,
  Rattus_norvegicus_msigdb_dna_repair = msig_rat_DNA_repair_direct_ensembl,
  Drosophila_melanogaster_msigdb_oxphos = msig_fly_oxphos_direct_ensembl,
  Drosophila_melanogaster_msigdb_apop = msig_fly_apoptosis_direct_ensembl,
  Drosophila_melanogaster_msigdb_dna_repair = msig_fly_DNA_repair_direct_ensembl,
  Dario_rerio_msigdb_oxphos = msig_zebra_oxphos_direct_ensembl,
  Dario_rerio_msigdb_apop = msig_zebra_apoptosis_direct_ensembl,
  Dario_rerio_msigdb_dna_repair = msig_zebra_DNA_repair_direct_ensembl,
  Macaca_mulatta_msigdb_oxphos = msig_macaque_oxphos_direct_ensembl,
  Macaca_mulatta_msigdb_apop = msig_macaque_apoptosis_direct_ensembl,
  Macaca_mulatta_msigdb_dna_repair = msig_macaque_DNA_repair_direct_ensembl,
  Gallus_gallus_msigdb_oxphos = msig_chicken_oxphos_direct_ensembl,
  Gallus_gallus_msigdb_apop = msig_chicken_apoptosis_direct_ensembl,
  Gallus_gallus_msigdb_dna_repair = msig_chicken_DNA_repair_direct_ensembl
)

save(msigdb_qc_ensembl_list, file = "data/msigdb_qc_ensembl_list.rda")
