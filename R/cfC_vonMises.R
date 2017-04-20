#' @title Characteristic function of von Mises distribution
#'
#' @description
#' cfC_vonMises(t) evaluates the characteristic function cf(t) of
#' the von Mises distribution (circular normal distribution) with the
#' parameters mu in (-pi,pi) and kappa > 0 (mu and 1/kappa are analogous to
#' mu and sigma^2, the mean and variance in the normal distribution), on a
#' circle e.g. the interval (-pi,pi), i.e.
#' cf(t) = besseli(t,kappa)/besseli(0,kappa) .* exp(1i*t*mu).
#'
#' @family Circular Probability Distribution
#' @family Characteristic Function
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Von_Mises_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param mu in (-pi, pi)
#' @param kappa > 0
#'
#' @return characteristic function cf(t) of the von Mises distribution with the parameters mu and kappa > 0
#' @usage cfC_vonMises(t, mu, kappa)
#'
#' @example Examples/example_cfC_vonMises.R
#'
#' @export

cfC_vonMises <- function(t, mu = 0, kappa = 1) {
  szt <- dim(t)
  t <- c(t)



  cf <- unlist(lapply(t, function(t) tryCatch((BesselI(kappa, abs(t), TRUE) / BesselI(kappa, 0, TRUE)) * exp(1i*t*mu), error = function(e) 0)))
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
