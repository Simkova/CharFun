#' @title Characteristic function of Trapezoidal distribution
#'
#' @description
#' cfS_Trapezoidal(t, a, c) (a > 0, c > 0, \eqn{c \le a}) evaluates the characteristic function cf(t) of
#' the Trapezoidal distribution on the interval (-a, a) with mode on the interval (-c, c)
#' (Trapezoidal distribution with mean = 0 and variance = ???
#' cfS_Trapezoidal(t, a, c) = (sin(w*at)/(w*at))*(sin((1-w)*at)/((1-w)*at))
#'
#' @family Continuous Probability distribution
#' @family Symetric Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Trapezoidal_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param a number, a > 0, default value a = 1
#' @param c number, (\eqn{0 \le c \le a}), default value c = 1/3
#' @return characteristic function cf(t) of the Triangular distribution on the interval (a, b) with mode c
#' @usage cfX_Triangular(t, a, c)
#'
#' @example Examples/example_cfS_Trapezoidal.R
#'
#' @export
#'
cfS_Trapezoidal <- function(t, a = 1, c = 1/3) {
  szt <- dim(t)
  t <- c(t)

  w = (1+c/a)/2

  cf <- cfS_Rectangular(w*a*t)*cfS_Rectangular((1-w*a)*t)
    cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
