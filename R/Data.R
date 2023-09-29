#' Ensembl Mito IDs
#'
#' A list of ensembl ids for mitochondrial genes (Ensembl version 105)
#'
#' @format A list of six vectors
#' \describe{
#'   \item{Mus_musculus_mito_ensembl}{Ensembl IDs for mouse mitochondrial genes}
#'   \item{Homo_sapiens_mito_ensembl}{Ensembl IDs for human mitochondrial genes}
#'   \item{Danio_rerio_mito_ensembl}{Ensembl IDs for zebrafish mitochondrial genes}
#'   \item{Rattus_norvegicus_mito_ensembl}{Ensembl IDs for rat mitochondrial genes}
#'   \item{Drosophila_melanogaster_mito_ensembl}{Ensembl IDs for fly mitochondrial genes}
#'   \item{Macaca_mulatta_mito_ensembl}{Ensembl IDs for macaque mitochondrial genes}
#'
#' }
#' @concept data
#'
"ensembl_mito_id"


#' Ensembl Ribo IDs
#'
#' A list of ensembl ids for ribosomal genes (Ensembl version 105)
#'
#' @format A list of seven vectors
#' \describe{
#'   \item{Mus_musculus_ribo_ensembl}{Ensembl IDs for mouse ribosomal genes}
#'   \item{Homo_sapiens_ribo_ensembl}{Ensembl IDs for human ribosomal genes}
#'   \item{Callithrix_jacchus_ribo_ensembl}{Ensembl IDs for marmoset ribosomal genes}
#'   \item{Danio_rerio_ribo_ensembl}{Ensembl IDs for zebrafish ribosomal genes}
#'   \item{Rattus_norvegicus_ribo_ensembl}{Ensembl IDs for rat ribosomal genes}
#'   \item{Drosophila_melanogaster_ribo_ensembl}{Ensembl IDs for fly ribosomal genes}
#'   \item{Macaca_mulatta_ribo_ensembl}{Ensembl IDs for macaque ribosomal genes}
#'
#' }
#' @concept data
#'
"ensembl_ribo_id"


#' QC Gene Lists
#'
#' A list gene symbols for qc percentages from MSigDB database.  The gene sets are from 3 MSigDB lists:
#' "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_APOPTOSIS", and "HALLMARK_DNA_REPAIR".
#'
#' @format A list of 18 vectors
#' \describe{
#'   \item{Homo_sapiens_msigdb_oxphos}{Genes in msigdb "HALLMARK_OXIDATIVE_PHOSPHORYLATION" list for human}
#'   \item{Homo_sapiens_msigdb_apop}{Genes in msigdb "HALLMARK_APOPTOSIS" list for human}
#'   \item{Homo_sapiens_msigdb_dna_repair}{Genes in msigdb "HALLMARK_DNA_REPAIR" list for human}
#'   \item{Mus_musculus_msigdb_oxphos}{Genes in msigdb "HALLMARK_OXIDATIVE_PHOSPHORYLATION" list for mouse}
#'   \item{Mus_musculus_msigdb_apop}{Genes in msigdb "HALLMARK_APOPTOSIS" list for mouse}
#'   \item{Mus_musculus_msigdb_dna_repair}{Genes in msigdb "HALLMARK_DNA_REPAIR" list for mouse}
#'   \item{Rattus_norvegicus_msigdb_oxphos}{Genes in msigdb "HALLMARK_OXIDATIVE_PHOSPHORYLATION" list for rat}
#'   \item{Rattus_norvegicus_msigdb_apop}{Genes in msigdb "HALLMARK_APOPTOSIS" list for rat}
#'   \item{Rattus_norvegicus_msigdb_dna_repair}{Genes in msigdb "HALLMARK_DNA_REPAIR" list for rat}
#'   \item{Drosophila_melanogaster_msigdb_oxphos}{Genes in msigdb "HALLMARK_OXIDATIVE_PHOSPHORYLATION" list for fly}
#'   \item{Drosophila_melanogaster_msigdb_apop}{Genes in msigdb "HALLMARK_APOPTOSIS" list for fly}
#'   \item{Drosophila_melanogaster_msigdb_dna_repair}{Genes in msigdb "HALLMARK_DNA_REPAIR" list for fly}
#'   \item{Dario_rerio_msigdb_oxphos}{Genes in msigdb "HALLMARK_OXIDATIVE_PHOSPHORYLATION" list for zebrafish}
#'   \item{Dario_rerio_msigdb_apop}{Genes in msigdb "HALLMARK_APOPTOSIS" list for zebrafish}
#'   \item{Dario_rerio_msigdb_dna_repair}{Genes in msigdb "HALLMARK_DNA_REPAIR" list for zebrafish}
#'   \item{Macaca_mulatta_msigdb_oxphos}{Genes in msigdb "HALLMARK_OXIDATIVE_PHOSPHORYLATION" list for macaque}
#'   \item{Macaca_mulatta_msigdb_apop}{Genes in msigdb "HALLMARK_APOPTOSIS" list for macaque}
#'   \item{Macaca_mulatta_msigdb_dna_repair}{Genes in msigdb "HALLMARK_DNA_REPAIR" list for macaque}
#'
#' }
#' @concept data
#'
#' @source MSigDB gene sets via msigdbr package \url{https://cran.r-project.org/package=msigdbr}
#'
"msigdb_qc_gene_list"


#' Immediate Early Gene (IEG) gene lists
#'
#' A list of gene symbols for immediate early genes
#'
#' @format A list of seven vectors
#' \describe{
#'   \item{Mus_musculus_IEGs}{Gene symbols for IEGs from source publication (see below)}
#'   \item{Homo_sapiens_ribo_ensembl}{Human gene symbols for homologous genes from mouse gene list}
#'
#' }
#' @concept data
#'
#' @source SI Table 4 from \url{https://doi.org/10.1016/j.neuron.2017.09.026}.
#'
"ieg_gene_list"
