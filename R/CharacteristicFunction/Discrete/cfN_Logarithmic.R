#' @title Characteristic function of Logarithmic distribution
#'
#' @description
#' cfN_Logarithmic(t, p) evaluates the characteristic function cf(t) of the
#' Logarithmic distribution defined on non-negative integers n =0, 1, ...,
#' with the parameters p (success probability p in [0,1]), i.e.
#' cfN_Logarithmic(t, p) = log(1 - p * exp(1i*t)) / log(1 - p)
#'
#' cfN_Logarithmic(t,p,cfX) evaluates the compound characteristic function
#' cf(t) = cfN_Logarithmic(-1i*log(cfX(t)), p), where cfX is function
#' handle of the characteristic function cfX(t) of a continuous distribution
#' and/or random variable X.
#'
#' Note that such CF is characteristic function of the compound distribution,
#' i.e. distribution of the random variable Y = X_1 + ... + X_N, where X_i ~ F_X
#' are i.i.d. random variables with common CF cfX(t), and N ~ F_N is
#' independent RV with its CF given by cfN(t).
#'
#' @family Discrete Probability Distribution
#' @family Characteristic Function
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Logarithmic_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param p success probability, \eqn{0 \le p \le 1}, default value p = 1/2
#' @param cfX function
#'
#' @return characteristic function cf(t) of the Logarithmic distribution p success probability
#' @usage cfN_Logarithmic(t, p)
#' cfN_Logarithmic(t, p, cfX)
#'
#' @example
#'
#' @export
#'
cfN_Logarithmic <- function(t, p = 0.5, cfX) {
  szt <- dim(t)
  t <- c(t)

  if (missing(cfX)) {expit <- exp(1i*t)
  } else {expit = cfX(t)}

  cf <- log(1 - p * expit) / log(1 - p);
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
