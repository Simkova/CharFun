#' @title Characteristic function of Gamma distribution
#'
#' @description
#' cfX_Gamma(t, alpha, beta) evaluates the characteristic function cf(t) of
#' the Gamma distribution with the parameters alpha (shape, alpha > 0) and
#' beta (rate, beta > 0), i.e.
#'
#' cfX_Gamma(t, alpha, beta) = (1 - it/beta)^(-alpha)
#'
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Gamma_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param alpha shape, alpha > 0, default value alpha = 1
#' @param beta rate > 0, default value beta = 1
#' @return characteristic function cf(t) of the Gamma distribution
#' @usage cfX_Gamma(t, alpha, beta)
#'
#' @example Examples/example_cfX_Gamma.R
#'
#' @export
#'
cfX_Gamma <- function(t, alpha = 1, beta = 1) {
  szt <- dim(t)
  t <- c(t)

  cf <- (1 - (0+1i)*t/beta)^(-alpha)
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
