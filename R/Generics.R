#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#################### GENERICS ####################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Convert objects to LIGER objects
#'
#' @param x An object to convert to class \code{LIGER}
#' @param ... Arguments passed to other methods
#'
#' @rdname as.LIGER
#' @export as.LIGER
#'

as.LIGER <- function(x, ...) {
  UseMethod(generic = 'as.LIGER', object = x)
}
