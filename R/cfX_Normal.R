#' @title Characteristic function of Normal distribution
#'
#' @description
#' cfX_Normal(t, mean, variance) evaluates the characteristic function cf(t) of
#' the Normal distribution with mean = mean and variance = variance: N(mean, variance))
#' cfX_Normal(t, mean, variance) = exp(imeant -1/2variance^2t^2)
#'
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Normal_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param mean real number, mean or expextation of the distribution, default value mean = 0
#' @param variance real number, standard deviation, variance > 0, default value variance = 1
#' @return characteristic function cf(t) of the normal distribution
#'
#' @example Examples/example_cfX_Normal.R
#'
#' @export
#'
cfX_Normal <- function(t, mean = 0, variance = 1) {
  szt <- dim(t)
  t <- c(t)

  cf <- exp(1i*mean*t -1/2*variance^2*t^2)
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
