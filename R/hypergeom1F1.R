#' Title
#'
#' @param z
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples
hypergeom1F1 <- function(z, a, b) {
  sz <- dim(z)

  z <- c(z)

  unlist(lapply(z, function(z) cchg(z, a, b)))

}
