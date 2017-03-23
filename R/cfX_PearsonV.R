#' @title Characteristic function of Pearson type V distribution
#'
#' @description
#' cfX_PearsonV(t, alpha, beta) evaluates the characteristic function cf(t) of
#' the Pearson type V distribution with the parameters alpha (shape, alpha > 0) and
#' beta (scale, beta > 0), computed for real vector argument t, i.e.
#'
#' cfX_PearsonV(t, alpha, beta) = (2/gamma(alpha)) * (-1i*t/beta)^(alpha/2) * besselk(alpha,2*sqrt(-1i*t/beta)),
#' where besselk(a,z) denotes the modified Bessel function of the second.
#'
#' @family Characteristic Function
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Pearson_distribution}
#'

#' @param t numerical real values (number, vector...)
#' @param alpha shape, alpha > 0, default value alpha = 1
#' @param beta scale > 0, default value beta = 1
#' @return characteristic function cf(t) of the Gamma distribution
#' @usage cfX_PearsonV(t, alpha, beta)
#'
#' @example
#'
#' @export
#'
cfX_PearsonV <- function(t, alpha = 1, beta = 1) {
  szt <- dim(t)

  cf <- unlist(lapply(t, function(t) tryCatch(BesselK(2*sqrt((0-1i*t)/beta), nu = alpha), error = function(e) NA)))
  cf <- (2/gamma(alpha)) * ((0-1i)*t/beta)^(alpha/2) * cf

  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
