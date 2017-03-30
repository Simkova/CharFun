#' @title Characteristic function of Student's t-distribution
#'
#' @description
#' cfS_StudentT(t, df) evaluates the characteristic function cf(t) of
#' the Student's t-distribution with the parameters df (degrees of freedom, df > 0)
#' computed for real vector argument t, i.e.
#'
#' cfS_StudentT(t,df) = besselK(df/2,abs(t)*sqrt(df),1) * exp(-abs(t)*sqrt(df)) *
#'                      (sqrt(df)*abs(t))^(df/2) / 2^(df/2-1)/gamma(df/2)
#'
#' @family Characteristic Function
#' @family Continuous Probability distribution
#' @family Symetric Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Student's_t-distribution}
#'
#' @param t numerical values (number, vector...)
#' @param df df degrees of freedom, df > 0
#'
#' @note bessel function nefunguje dobre
#'
#' @return characteristic function cf(t) of the Student's t-distribution
#'
#' @example Examples/example_cfS_StudentT.R
#'
#' @export
#'
cfS_StudentT <- function(t, df) {
  szt <- dim(t)
  t <- c(t)

  cf <- unlist(lapply(t, function(t) tryCatch(BesselK(abs(t)*sqrt(df), df/2, TRUE), error = function(e) NA, message("Error in Bessel function: nu must be integer."))))
  cf <- exp(-abs(t)*sqrt(df)) * (sqrt(df)*abs(t))^(df/2) / 2^(df/2-1)/gamma(df/2) * cf
  cf[t == 0] <- 1

  dim(cf) <- szt

  return(cf)
}
