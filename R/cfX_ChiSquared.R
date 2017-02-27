#' @title Characteristic function of Chi-Squared distribution
#'
#' @description
#' cfX_ChiSquared(t,df) evaluates the characteristic function cf(t) of
#' the CHI-SUQARED distribution with the parameter df (degrees of freedom, df > 0),
#' i.e.
#' cfX_ChiSquared(t, df) = (1 - 2it )^(-df/2) = cfX_GAMMA(t, df/2, 1/2)
#'
#' @family Characteristic Function
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Chi-squared_distribution}
#'

#' @param t numerical values (number, vector...)
#' @param df degrees of freedom, df > 0
#' @return characteristic function cf(t) of the CHI-SUQARED distribution
#' @usage cfX_ChiSquared(t, df)
#'
#' @example
#'
#' @export
#'
cfX_ChiSquared <- function(t, df) {
  szt <- dim(t)
  t <- c(t)

  cf <- (1 - 2*(0+1i)*t)^(-df/2)
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
