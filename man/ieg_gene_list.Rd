% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data.R
\docType{data}
\name{ieg_gene_list}
\alias{ieg_gene_list}
\title{Immediate Early Gene (IEG) gene lists}
\format{
A list of seven vectors
\describe{
\item{Mus_musculus_IEGs}{Gene symbols for IEGs from source publication (see below)}
\item{Homo_sapiens_IEGs}{Human gene symbols for homologous genes from mouse gene list}

}
}
\source{
Mouse gene list is from: SI Table 4 from \doi{10.1016/j.neuron.2017.09.026}.  Human
gene list was compiled by first creating homologous gene list using biomaRt and then adding some manually curated
homologs according to HGNC.
}
\usage{
ieg_gene_list
}
\description{
Gene symbols for immediate early genes
}
\concept{data}
\keyword{datasets}
