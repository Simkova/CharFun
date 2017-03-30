#' @title Characteristic function of Triangular distribution
#'
#' @description
#' cfS_Triangular(t, a) evaluates the characteristic function cf(t) of
#' the Triangular distribution on the interval (-a, a) with mode 0
#' (Triangular distribution with mean = 0 and variance = 1/18(2a^2 + a)
#' cfS_Triangula(t, a) = (2 - 2cos(at)) / (a^2 * t^2)
#'
#' @family Characteristic Function
#' @family Continuous Probability distribution
#' @family Symetric Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Triangular_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param a number, a > 0, default value a = 1
#' @return characteristic function cf(t) of the Triangular distribution on the interval (-a, a) with mode 0
#' @usage cfX_Triangular(t, a)
#'
#' @example Examples/example_cfS_Triangular.R
#'
#' @export

cfS_Triangular <- function(t, a = 1) {
  szt <- dim(t)
  t <- c(t)

  cf <- (2 - 2*cos(a*t)) / ((a*t)^2)
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
