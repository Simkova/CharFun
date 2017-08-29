#' @title Characteristic function of Rectangular distribution
#'
#' @description
#' cfS_Rectangular(t) evaluates the characteristic function cf(t) of
#' the symmetric zero-mean Rectangular distribution on the interval (-a, a)
#' (Rectangular distribution with mean = 0 and variance = 1/3*a^2)
#' cfS_Rectangular(t, a) = sin(at)/at
#'
#' @family Continuous Probability distribution
#' @family Symetric Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)}
#'
#' @param t numerical values (number, vector...)
#' @param a number, a > 0, default value a = 1
#'
#' @return characteristic function cf(t) of the Rectangular distribution on the interval (-a, a)
#'
#' @example R/Examples/example_cfS_Rectangular.R
#'
#' @export
#'
cfS_Rectangular <- function(t, a = 1) {
  szt <- dim(t)
  t <- c(t)

  cf <- pmin(1, sin(a*t)/(a*t))
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
