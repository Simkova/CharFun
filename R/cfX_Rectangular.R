#' @title Characteristic function of Rectangular distribution
#'
#' @description
#' cfX_Rectangular(t, a, b) evaluates the characteristic function cf(t) of
#' the Rectangular distribution on the interval (a, b)
#' (Rectangular distribution with mean = (a + b)/2 and variance = 1/12(b - a)^2)
#' cfX_Rectangular(t, a, b) = (exp(ibt) - exp(iat))/(i(b - a)t)
#'
#' @family Characteristic Function
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Normal_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param a number, default value a = -1
#' @param b number, default value b = 1
#' @return characteristic function cf(t) of the Rectangular distribution on the interval (a, b)
#' @usage cfX_Rectangular(t, a, b)
#'
#' @example Examples/example_cfX_Rectangular.R
#'
#' @export

cfX_Rectangular <- function(t, a = -1, b = 1) {
  szt <- dim(t)
  t <- c(t)

  cf <- (exp((0+1i)*b*t) - exp((0+1i)*a*t))/((0+1i)*(b - a)*t)
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
