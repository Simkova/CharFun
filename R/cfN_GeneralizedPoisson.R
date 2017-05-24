#' @title Characteristic function of Generalized Poisson distribution
#'
#' @description
#' cfN_GeneralizedPoisson(t, a, p) evaluates the characteristic function cf(t) of the
#' Generalized Poisson with the parameter a (variable mean, a > 0),
#' and p (success probability, \eqn{0 \le p \le 1}) i.e.
#' cfN_GeneralizedPoisson(t, a, p) = exp(a*(sum_{j=1}^Inf ((p*j)^(j-1)*e^(-p*j)/j!)*e^(1i*t*j)-1))
#'
#' The Generalized-Poisson distribution is equivalent
#' with the Borel-Tanner distribution with parameters (p,m)
#'
#' cfN_GeneralizedPoisson(t, a, p, cfX) evaluates the compound characteristic function
#' cf(t) = cfN_GeneralizedPoisson(-1i*log(cfX(t)), a, p), where cfX is function
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
#' @seealso For more details
#'
#' @param t numerical values (number, vector...)
#' @param a variable mean, a > 0
#' @param p success probability, \eqn{0 \le p \le 1}, default value p = 1/2
#' @param cfX function
#'
#' @return characteristic function cf(t) of the Poisson distribution with rate lambda
#'
#' @example Examples/example_cfN_GeneralizedPoisson.R
#'
#' @export

cfN_GeneralizedPoisson <- function(t, a, p = 0.5, cfX) {
  szt <- dim(t)
  t <- c(t)

  if (missing(cfX)) {expit <- exp(1i*t)
  } else {expit = cfX(t)}

  exppt = exp(-p) * expit
  c     = exppt;
  cf    = c;
  N     = ceiling(-36.84/(1+log(p)-p))
  for (j in seq(1:N))
  {
    c  = p * ((1+j)/j)^(j-1) * exppt * c
    cf = cf + c
  }

  cf = exp(a * (cf-1));
  cf[t==0] = 1;

  dim(cf) <- szt

  return(cf)
}
