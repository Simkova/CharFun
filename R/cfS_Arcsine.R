#' @title Characteristic function of symmetric zero-mean Arcsine distribution
#'
#' @description
#' cfS_Arcsine(t) evaluates the characteristic function cf(t) of
#' the symmetric zero-mean Arcsine distribution on the interval
#' (-1,1) (U-shaped distribution with mean = 0 and variance = 1/2
#' \deqn{cfS_Arcsine(t) = besselj(0,t)}
#'
#' @family Characteristic Function
#' @family Continuous Probability distribution
#' @family Symetric Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Arcsine_distribution}
#'

#' @param t numerical values (number, vector...)
#' @return characteristic function cf(t) of the Arcsine distribution
#' @usage cfS_Arcsine(t)
#'
#' @example Examples/example_cfS_Arcsine.R
#'
#' @export
#'
cfS_Arcsine <- function(t) {
  szt <- dim(t)
  t <- c(t)

  cf <- unlist(lapply(t, function(t) tryCatch(BesselJ(t, 0), error = function(e) 0)))
  cf[t == 0] <- 1
  cf[]

  dim(cf) <- szt

  return(cf)
}
