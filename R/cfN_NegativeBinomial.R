#' @title Characteristic function of Negative-Binomial distribution
#'
#' @description
#' cfN_NegativeBinomial(t, r, p) evaluates the characteristic function cf(t) of the
#' Negative-Binomial distribution with the parameters r (number of failures until
#' the experiment is stopped, r in N) and p (success probability in each experiment, p in [0,1]), i.e.
#' cfN_NegativeBinomial(t, r, p) = p^r * (1 - (1-p) * e^(1i*t))^(-r)
#'
#' cfN_NegativeBinomial(t, r, p, cfX) evaluates the compound characteristic function
#' cf(t) = cfN_NegativeBinomial(-1i*log(cfX(t)), r, p), where cfX is function
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
#' \url{https://en.wikipedia.org/wiki/Negative_binomial_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param r number of trials
#' @param p success probability, \eqn{0 \le p \le 1}, default value p = 1/2
#' @param cfX function
#'
#' @return characteristic function cf(t) of the Negative-Binomial distribution with n trials and p success probability
#' @usage cfN_NegativeBinomial(t, n, p)
#' cfN_NegativeBinomial(t, n, p, cfX)
#'
#' @example Examples/example_cfN_NegativeBinomial.R
#'
#' @export
#'
#'
#'
cfN_NegativeBinomial <- function(t, r = 10, p = 0.5, cfX) {
  szt <- dim(t)
  t <- c(t)

  if (missing(cfX)) {expit <- exp(1i*t)
  } else {expit = cfX(t)}

  cf <- p^r * (1 - (1-p)*expit)^(-r)
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
