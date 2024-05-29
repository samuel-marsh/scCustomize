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
#'   \item{Gallus_gallus_ribo_ensembl}{Ensembl IDs for chicken mitochondrial genes}
#' }
#' @concept data
#' @source See data-raw directory for scripts used to create gene list.
#'
#'
"ensembl_mito_id"


#' Ensembl Ribo IDs
#'
#' A list of ensembl ids for ribosomal genes (Ensembl version 105)
#'
#' @format A list of eight vectors
#' \describe{
#'   \item{Mus_musculus_ribo_ensembl}{Ensembl IDs for mouse ribosomal genes}
#'   \item{Homo_sapiens_ribo_ensembl}{Ensembl IDs for human ribosomal genes}
#'   \item{Callithrix_jacchus_ribo_ensembl}{Ensembl IDs for marmoset ribosomal genes}
#'   \item{Danio_rerio_ribo_ensembl}{Ensembl IDs for zebrafish ribosomal genes}
#'   \item{Rattus_norvegicus_ribo_ensembl}{Ensembl IDs for rat ribosomal genes}
#'   \item{Drosophila_melanogaster_ribo_ensembl}{Ensembl IDs for fly ribosomal genes}
#'   \item{Macaca_mulatta_ribo_ensembl}{Ensembl IDs for macaque ribosomal genes}
#'   \item{Gallus_gallus_ribo_ensembl}{Ensembl IDs for chicken ribosomal genes}
#' }
#' @concept data
#' @source See data-raw directory for scripts used to create gene list.
#'
"ensembl_ribo_id"


#' Ensembl Hemo IDs
#'
#' A list of ensembl ids for hemoglobin genes (Ensembl version 112)
#'
#' @format A list of six vectors
#' \describe{
#'   \item{Mus_musculus_hemo_ensembl}{Ensembl IDs for mouse hemoglobin genes}
#'   \item{Homo_sapiens_hemo_ensembl}{Ensembl IDs for human hemoglobin genes}
#'   \item{Danio_rerio_hemo_ensembl}{Ensembl IDs for zebrafish hemoglobin genes}
#'   \item{Rattus_norvegicus_hemo_ensembl}{Ensembl IDs for rat hemoglobin genes}
#'   \item{Drosophila_melanogaster_hemo_ensembl}{Ensembl IDs for fly hemoglobin genes}
#'   \item{Macaca_mulatta_hemo_ensembl}{Ensembl IDs for macaque hemoglobin genes}
#'   \item{Gallus_gallus_ribo_ensembl}{Ensembl IDs for chicken hemoglobin genes}
#' }
#' @concept data
#' @source See data-raw directory for scripts used to create gene list.
#'
#'
"ensembl_hemo_id"


#' QC Gene Lists
#'
#' Gene symbols for qc percentages from MSigDB database.  The gene sets are from 3 MSigDB lists:
#' "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_APOPTOSIS", and "HALLMARK_DNA_REPAIR".
#'
#' @format A list of 21 vectors
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
#'   \item{Gallus_gallus_msigdb_oxphos}{Genes in msigdb "HALLMARK_OXIDATIVE_PHOSPHORYLATION" list for chicken}
#'   \item{Gallus_gallus_msigdb_apop}{Genes in msigdb "HALLMARK_APOPTOSIS" list for chicken}
#'   \item{Gallus_gallus_msigdb_dna_repair}{Genes in msigdb "HALLMARK_DNA_REPAIR" list for chicken}
#' }
#' @concept data
#'
#' @source MSigDB gene sets (gene symbols) via msigdbr package \url{https://cran.r-project.org/package=msigdbr}.  See data-raw directory for scripts used to create gene list.
#'
"msigdb_qc_gene_list"


#' QC Gene Lists
#'
#' Ensembl IDs for qc percentages from MSigDB database.  The gene sets are from 3 MSigDB lists:
#' "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_APOPTOSIS", and "HALLMARK_DNA_REPAIR".
#'
#' @format A list of 21 vectors
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
#'   \item{Gallus_gallus_msigdb_oxphos}{Genes in msigdb "HALLMARK_OXIDATIVE_PHOSPHORYLATION" list for chicken}
#'   \item{Gallus_gallus_msigdb_apop}{Genes in msigdb "HALLMARK_APOPTOSIS" list for chicken}
#'   \item{Gallus_gallus_msigdb_dna_repair}{Genes in msigdb "HALLMARK_DNA_REPAIR" list for chicken}
#' }
#' @concept data
#'
#' @source MSigDB gene sets (ensembl IDs) via msigdbr package \url{https://cran.r-project.org/package=msigdbr}.  See data-raw directory for scripts used to create gene list.
#'
"msigdb_qc_ensembl_list"


#' Immediate Early Gene (IEG) gene lists
#'
#' Gene symbols for immediate early genes
#'
#' @format A list of seven vectors
#' \describe{
#'   \item{Mus_musculus_IEGs}{Gene symbols for IEGs from source publication (see below)}
#'   \item{Homo_sapiens_IEGs}{Human gene symbols for homologous genes from mouse gene list}
#'
#' }
#' @concept data
#'
#' @source Mouse gene list is from: SI Table 4 from \doi{10.1016/j.neuron.2017.09.026}.  Human
#' gene list was compiled by first creating homologous gene list using biomaRt and then adding some manually curated
#' homologs according to HGNC.  See data-raw directory for scripts used to create gene list.
#'
"ieg_gene_list"
