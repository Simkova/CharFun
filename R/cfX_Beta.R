#' @title Characteristic function of Beta distribution
#'
#' @description
#' cfX_Beta(t,alpha, beta) evaluates the characteristic function cf(t) of
#' the Beta distribution with the parameter alpha (shape, alpha > 0) and beta (shape, beta > 0)
#' defined on the interval (0,1), i.e. beta distribution with the Mean = alpha / (alpha + beta)
#' and the Variance = (alpha*beta) / ((alpha+beta)^2*(alpha+beta+1)).
#' Then, the standard deviation is given by STD = sqrt(Variance)
#' i.e.
#' cf(t) = cfX_Beta(t,alpha,beta) = 1F1(alpha ,alpha + beta , i*t)
#' where 1F1(.;.;.) is the Confluent hypergeometric function.
#'
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Beta_distribution}
#'

#' @param t numerical values (number, vector...)
#' @param alpha shape, aplha > 0, default value aplha = 1
#' @param beta shape, beta > 0, default value beta = 1
#' @return characteristic function cf(t) of the Beta distribution
#'
#' @example Examples/example_cfX_Beta.R
#'
#' @export
#'
cfX_Beta <- function(t, alpha = 1, beta = 1) {
     szt <- dim(t)
     t <- c(t)

     cf <- hypergeom1F1(1i*t, alpha, alpha + beta)


     cf[t == 0] <- 1

     dim(cf) <- szt

     return(cf)
}
