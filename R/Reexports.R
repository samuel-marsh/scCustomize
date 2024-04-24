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
