#' @title Characteristic function of Delaporte distribution
#'
#' @description
#' cfN_Delaporte(t, a, b, c) evaluates the characteristic function cf(t) of the
#' Delaporte distribution with the parameters a (parameter of variable mean, a > 0)
#' b (parameter of variable mean, b > 0 ), and c (fixed mean, c > 0), i.e.
#' cfN_Delaporte(t, a, b, c) = (b/(1+b))^a * (1-e^(1i*t)/(b+1))^(-a) * exp(-c*(1-e^(1i*t)))
#'
#' cfN_Delaporte(t, a, b, c, cfX) evaluates the compound characteristic function
#' cf(t) = cfN_Delaporte(-1i*log(cfX(t)), a, b, c), where cfX is function
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
#' \url{https://en.wikipedia.org/wiki/Delaporte_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param a variable mean, a > 0
#' @param b variable mean, b > 0
#' @param c fixed mean, c > 0
#' @param cfX function
#'
#' @return characteristic function cf(t) of the Delaporte distribution with a and b variable means and c fixed mean
#'
#' @example Examples/example_cfN_Delaporte.R
#'
#' @export
#'
#'


cfN_Delaporte <- function(t, a, b, c, cfX) {
  szt <- dim(t)
  t <- c(t)

  if (missing(cfX)) {expit <- exp(1i*t)
  } else {expit = cfX(t)}

  cf <- (b/(1+b))^a * (1-expit/(b+1))^(-a) * exp(-c*(1-expit));
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
