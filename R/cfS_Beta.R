#' @title Characteristic function of Beta distribution
#'
#' @description
#' cfX_Beta(t,alpha, beta) evaluates the characteristic function cf(t) of
#' the symmetric zero-mean Beta distribution with shape parameter theta >0,
#' defined on the interval (-1,1), i.e. symmetric beta distribution with
#' zero mean and variance VAR = 1/(1+2*theta). The standard deviation is
#' given by STD = sqrt(1/(1+2*theta))), i.e.
#' cf(t) = cfS_Beta(t,theta)
#'       = gamma(1/2+theta) * (t/2)^(1/2-theta) * besselj(theta-1/2,t).
#' Special cases (for specific values of the shape parameter theta):
#'  1) theta = 1/2; Arcsine distribution on (-1,1):     cf(t) = besselj(0,t).
#'  2) theta = 1;   Rectangular distribution on (-1,1): cf(t) = sin(t)/t;
#'
#' @family Continuous Symetric distribution
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Beta_distribution}
#'

#' @param t numerical values (number, vector...)
#'
#' @param theta default value theta = 1
#'
#' @return characteristic function cf(t) of the Beta distribution
#'
#' @example Examples/example_cfX_Beta.R
#'
#' @export
#'
cfS_Beta <- function(t, theta = 1) {
  szt <- dim(t)
  t <- c(t)

  cf <- pmin(1, gamma(0.5+theta)*(0.5*t)^(0.5-theta) * unlist(lapply(t, function(x) tryCatch(BesselJ(x,theta-0.5), error = function(e) 0))))

  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
