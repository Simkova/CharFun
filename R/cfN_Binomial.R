#' @title Characteristic function of Binomial distribution
#'
#' @description
#' cfN_Binomial(t, n, p) evaluates the characteristic function cf(t) of the
#' Binomial distribution with the parameters n (number of trials, n in N)
#' and p (success probability, p in [0,1]), i.e.
#' cfN_Binomial(t, n, p) = (1 - p + p*exp(1i*t))^n
#'
#' cfN_Binomial(t, n, p, cfX) evaluates the compound characteristic function
#' cf(t) = cfN_Binomial(-1i*log(cfX(t)), n, p), where cfX is function
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
#' \url{https://en.wikipedia.org/wiki/Binomial_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param n number of trials
#' @param p success probability, \eqn{0 \le p \le 1}, default value p = 1/2
#' @param cfX function
#'
#' @return characteristic function cf(t) of the Binomial distribution with n trials and p success probability
#'
#' @example Examples/example_cfN_Binomial.R
#'
#' @export

cfN_Binomial <- function(t, n = 10, p = 1/2, cfX) {
  szt <- dim(t)
  t <- c(t)

  if (missing(cfX)) {expit <- exp(1i*t)
    } else {expit = cfX(t)}

  cf <- (1 - p + p*expit)^n
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
