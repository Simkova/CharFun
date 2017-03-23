#' @title Characteristic function of Inverse Gamma distribution
#'
#' @description
#' cfX_InverseGamma(t,alpha,beta) evaluates the characteristic function cf(t) of
#' the Inverse Gamma distribution with the parameters alpha (shape, alpha > 0) and
#' beta (rate, beta > 0), i.e.
#'
#' cfX_InverseGamma(t, alpha, beta) = (1 - it/beta)^(-alpha)
#'
#' @family Characteristic Function
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Inverse-gamma_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param alpha shape, alpha > 0, default value alpha = 1
#' @param beta rate > 0, default value beta = 1
#' @return characteristic function cf(t) of the Inverse Gamma distribution
#' @usage cfX_InverseGamma(t, alpha, beta)
#'
#' @example
#'
#' @export
#'
cfX_InverseGamma <- function(t, alpha = 1, beta = 1) {
  szt <- dim(t)
  t <- c(t)

  cf <- unlist(lapply(t, function(t) tryCatch(BesselK(sqrt((0-4i)*beta*t), alpha), error = function(e) NA)))
  cf <- 2/gamma(alpha) * ((0-1i)*beta*t)^(alpha/2) * cf
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
