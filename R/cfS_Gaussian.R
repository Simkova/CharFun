#' @title Characteristic function of Normal distribution N(0,1)
#'
#' @description
#' cfS_Gaussian(t) evaluates the characteristic function cf(t) of
#' the symmetric zero-mean standard Gaussian distribution (i.e. the standard
#' normal distribution with mean = 0 and variance = 1: N(0, 1))
#' cfS_Gaussian(t) = exp(-t^2/2)
#'
#' @family Continuous Probability distribution
#' @family Symetric Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Normal_distribution}
#'

#' @param t numerical values (number, vector...)
#' @return characteristic function cf(t) of the normal distribution N(0, 1)
#' @usage cfS_Gaussian(t)
#'
#' @example Examples/example_cfS_Gaussian.R
#'
#' @export
#'
cfS_Gaussian <- function(t) {
  szt <- dim(t)
  t <- c(t)

  cf <- exp(-t^2/2)
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
