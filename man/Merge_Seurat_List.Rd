% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Object_Utilities.R
\name{Merge_Seurat_List}
\alias{Merge_Seurat_List}
\title{Merge a list of Seurat Objects}
\usage{
Merge_Seurat_List(
  list_seurat,
  add.cell.ids = NULL,
  merge.data = TRUE,
  project = "SeuratProject"
)
}
\arguments{
\item{list_seurat}{list composed of multiple Seurat Objects.}

\item{add.cell.ids}{A character vector of equal length to the number of objects in \code{list_seurat}.
Appends the corresponding values to the start of each objects' cell names.  See \code{\link[SeuratObject]{merge}}.}

\item{merge.data}{Merge the data slots instead of just merging the counts (which requires renormalization).
This is recommended if the same normalization approach was applied to all objects.
See \code{\link[SeuratObject]{merge}}.}

\item{project}{Project name for the Seurat object. See \code{\link[SeuratObject]{merge}}.}
}
\value{
A Seurat Object
}
\description{
Enables easy merge of a list of Seurat Objects.  See  See \code{\link[SeuratObject]{merge}} for more information,
}
\examples{
\dontrun{
object_list <- list(obj1, obj2, obj3, ...)
merged_object <- Merge_Seurat_List(list_seurat = object_list)
}

}
\concept{misc_util}
