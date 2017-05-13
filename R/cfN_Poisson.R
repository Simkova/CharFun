#' @title Characteristic function of Poisson distribution
#'
#' @description
#' cfN_Poisson(t, lambda) evaluates the characteristic function cf(t) of the
#' Poisson distribution with the rate parameter lambda > 0, i.e.
#' cfN_Poisson(t, lambda) = exp(lambda*(exp(1i*t)-1))
#'
#' cfN_Poisson(t, lambda, cfX) evaluates the compound characteristic function
#' cf(t) = cfN_Poisson(-1i*log(cfX(t)), lambda), where cfX is function
#' handle of the characteristic function cfX(t) of a continuous distribution
#' and/or random variable X.
#'
#' Note that such CF is characteristic function of the compound distribution,
#' i.e. distribution of the random variable Y = X_1 + ... + X_N, where X_i ~ F_X
#' are i.i.d. random variables with common CF cfX(t), and N ~ F_N is
#' independent RV with its CF given by cfN(t).
#'
#' @family Discrete Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Poisson_distribution}
#' \url{https://en.wikipedia.org/wiki/Compound_Poisson_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param lambda rate, lambda > 0, default value lambda = 1
#' @param cfX function
#'
#' @return characteristic function cf(t) of the Poisson distribution with rate lambda
#' @usage cfN_Poisson(t, lambda)
#' cfN_Poisson(t, lambda, cfX)
#'
#' @example Examples/example_cfN_Poisson.R
#'
#' @export
#'
cfN_Poisson <- function(t, lambda = 1, cfX) {
  szt <- dim(t)
  t <- c(t)

  if (missing(cfX)) {expit <- exp(1i*t)
  } else {expit = cfX(t)}

  cf <- exp(lambda * (expit-1))
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
