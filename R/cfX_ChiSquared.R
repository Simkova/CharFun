#' @title Characteristic function of Noncentral Chi-Squared distribution
#'
#' @description
#' cfX_ChiSquared(t, df, ncp = 0) evaluates the characteristic function cf(t) of
#' the Chi-Squared distribution with the parameter df (degrees of freedom, df > 0)
#' and npc (non-centrality, \eqn{npc \ge 0}), i.e.
#'
#' cfX_ChiSquared(t, df, npc) = exp((i*npc*t)/(1 - 2it)) / (1 - 2it)^(df/2)
#'
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution}
#'

#' @param t numerical values (number, vector...)
#' @param df degrees of freedom
#' @param npc non-centrality parameter, default value = 0
#' @return characteristic function cf(t) of the CHI-SUQARED distribution
#'
#' @example R/Examples/example_cfX_ChiSquared.R
#'
#' @export
#'
cfX_ChiSquared <- function(t, df, npc = 0) {
  szt <- dim(t)
  t <- c(t)

  cf <- exp(((0+1i)*npc*t)/(1-2*(0+1i)*t))/(1-2*(0+1i)*t)^(df/2)
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
