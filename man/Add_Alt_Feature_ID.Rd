% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Object_Utilities.R
\name{Add_Alt_Feature_ID}
\alias{Add_Alt_Feature_ID}
\title{Add Alternative Feature IDs}
\usage{
Add_Alt_Feature_ID(
  seurat_object,
  features_tsv_file = NULL,
  hdf5_file = NULL,
  assay = NULL
)
}
\arguments{
\item{seurat_object}{object name.}

\item{features_tsv_file}{output file from Cell Ranger used for creation of Seurat object.
(Either provide this of \code{hdf5_file})}

\item{hdf5_file}{output file from Cell Ranger used for creation of Seurat object.
(Either provide this of \code{features_tsv_file})}

\item{assay}{name of assay(s) to add the alternative features to.  Can specify "all"
to add to all assays.}
}
\value{
Seurat Object with new entries in the \code{obj@assays$ASSAY@meta.data} slot.
}
\description{
Add alternative feature ids to the assay level meta.data slot in Assay5 compatible object (Seurat V5.0.0 or greater)
}
\examples{
\dontrun{
# Using features.tsv.gz file
   # Either file from filtered or raw outputs can be used as they are identical.
obj <- Add_Alt_Feature_ID(seurat_object = obj,
features_tsv = "sample01/outs/filtered_feature_bc_matrix/features.tsv.gz", assay = "RNA")

#' # Using hdf5 file
   # Either filtered_feature_bc or raw_feature_bc can be used as the features slot is identical
   # Though it is faster to load filtered_feature_bc file due to droplet filtering
obj <- Add_Alt_Feature_ID(seurat_object = obj,
hdf5_file = "sample01/outs/outs/filtered_feature_bc_matrix.h5", assay = "RNA")
}

}
\concept{get_set_util}
