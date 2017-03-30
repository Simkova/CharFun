#' @title Characteristic function of Geometric distribution
#'
#' @description
#' cfN_Geometric(t, p, type, cfX) evaluates the characteristic function cf(t) of the
#' Geometric distribution.
#'
#'  The standard Geometric distribution (type = "standard" or "zero") is
#'  defined on non-negative integers k = 0,1,... .
#'
#'  The shifted Geometric distribution (type = "shifted")
#'  is defined on positive integers k = 1,2,... .
#'
#' Both types are parametrized by the success probability parameter p in [0,1]), i.e.
#' cfN_Geometric(t, p, "standard") = p / (1 - (1-p) * exp(1i*t)),
#' cfN_Geometric(t, p, "shifted")  = exp(1i*t) * (p / (1 - (1-p) * exp(1i*t))).
#'
#' cfN_Geometric(t, p, type, cfX) evaluates the compound characteristic function
#' cf(t) = Geometric(-1i*log(cfX(t)), p), where cfX is function
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
#' @family Discrete Distribution
#' @family Characteristic Function
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Geometric_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param p success probability, \eqn{0 \le p \le 1}, default value p = 1/2
#' @param type standard = 1, shifted = 2, default type = standard
#' @param cfX function
#'
#' @return characteristic function cf(t) of the Geometric distribution with p success probability
#' @usage cfN_Geometric(t, p, type)
#' cfN_Geometric(t, p, type, cfX)
#'
#' @example Examples/example_cfN_Geometric.R
#'
#' @export

cfN_Geometric <- function(t, p = 1/2, type = "standard", cfX) {
  szt <- dim(t)
  t <- c(t)

  if (missing(cfX)) {expit <- exp(1i*t)
  } else {expit = cfX(t)}

  cf <- switch(type,
         standard = 1,
         shifted = expit)

  cf <- cf * (p / (1 - (1-p) * expit))
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
