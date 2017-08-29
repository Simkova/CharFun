#' @title Characteristic function of Polya-Eggenberger
#'
#' @description
#' cfN_PolyaEggenberger(t, a, b, m) evaluates the characteristic function cf(t) of the
#' Polya-Eggenberger distribution with the parameters a (a real), b (b real), and m (m integer), i.e.
#' cfN_PolyaEggenberger(t, a, b, m) = 2F1(-m,a,a+b,1-e^(1i*t))
#' where 2F1 denotes the Gauss hypergeometric function.
#'
#' cfN_PolyaEggenberger(t, a, b, m) evaluates the compound characteristic function
#' cf(t) = cfN_PolyaEggenberge(-1i*log(cfX(t)), a, b, m), where cfX is function
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
#'
#' @param t numerical values (number, vector...)
#' @param a real number
#' @param b real number
#' @param m integer
#' @param cfX function
#'
#' @return characteristic function cf(t) of the Polya-Eggenberger distribution
#'
#' @example R/Examples/example_cfN_PolyaEggenberger.R
#'
#' @export
#'
cfN_PolyaEggenberger <- function(t, a, b, m, cfX) {
  szt <- dim(t)
  t <- c(t)

  if (missing(cfX)) {expit <- exp(1i*t)
  } else {expit = cfX(t)}

  const1 <- 1
  const2 <- rep.int(1, length(t))
  cf <- const2

  for (i in 0:(m-1)) {
    const1 <- const1 * (b + i)/(a + b + i)
    const2 <- (-m + i) * (a + i)/(-m - b + 1 + i)/(i + 1) * const2 * expit
    cf <- cf + const2
  }

  cf = cf * const1
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
