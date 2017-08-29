#' @title Characteristic function of Exponential distribution
#'
#' @description
#' cfX_Exponential(t, lambda) evaluates the characteristic function cf(t) of
#' the Exponential distribution with the parameter lambda (rate, lambda > 0)
#'
#' cfX_Exponential(t, lambda) = lambda / (lambda - it)
#'
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Exponential_distribution}
#'

#' @param t numerical values (number, vector...)
#' @param lambda rate, lambda > 0, default value lambda = 1
#' @return characteristic function cf(t) of the Exponential distribution
#'
#' @example R/Examples/example_cfX_Exponential.R
#'
#' @export

cfX_Exponential <- function(t, lambda = 1) {
  szt <- dim(t)
  t <- c(t)

  cf <- cfX_Gamma(t, 1, lambda)
  cf[t == 0] <- 1
  dim(cf) <- szt

  return(cf)
}

