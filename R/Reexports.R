#' @importFrom SeuratObject as.Seurat
#' @export
#' @note See \code{\link{as.Seurat.liger}} for scCustomize extension of this generic to converting Liger objects.
#'
#'
SeuratObject::as.Seurat

#' @importFrom SeuratObject WhichCells
#' @export
#' @note See \code{\link{WhichCells.liger}} for scCustomize extension of this generic to extract cell barcodes.
#'
#'
SeuratObject::WhichCells

#' @importFrom SeuratObject Cells
#' @export
#' @note See \code{\link{Cells.liger}} for scCustomize extension of this generic to extract cell barcodes.
#'
#'
SeuratObject::Cells

#' @importFrom SeuratObject Features
#' @export
#' @note See \code{\link{Features.liger}} for scCustomize extension of this generic to extract dataset features.
#'
#'
SeuratObject::Features

#' @importFrom SeuratObject Embeddings
#' @export
#' @note See \code{\link{Embeddings.liger}} for scCustomize extension of this generic to extract embeddings.
#'
#'
SeuratObject::Embeddings

#' @importFrom SeuratObject Idents
#' @export
#' @note See \code{\link{Idents.liger}} for scCustomize extension of this generic to extract cell identities.
#'
#'
SeuratObject::Idents

#' @importFrom SeuratObject Idents<-
#' @export
#' @note See \code{\link{Idents.liger}} for scCustomize extension of this generic to set cell identities.
#'
#'
SeuratObject::`Idents<-`
