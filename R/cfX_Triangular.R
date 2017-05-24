#' @title Characteristic function of Triangular distribution
#'
#' @description
#' cfX_Triangular(t, a, b, c) (\eqn{a \le c \le b}) evaluates the characteristic function cf(t) of
#' the Triangular distribution on the interval (a, b) with mode c
#' (Triangular distribution with mean = a + b + c)/3 and variance = 1/18(a^2 + b^2 + c^2 + - ab - ac - bc))
#' cfX_Triangula(t, a, b, c) = -2((b-c)exp(iat) - (b-a)exp(ict) + (c-a)exp(ibt))/((b-a)(c-a)(b-c)t^2)
#'
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Triangular_distribution}
#'
#' @param t numerical values (number, vector...)
#' @param a number, default value a = -1
#' @param b number, default value b = 1
#' @param c number, (\eqn{a \le c \le b}), default value c = 0
#' @return characteristic function cf(t) of the Triangular distribution on the interval (a, b) with mode c
#'
#' @example Examples/example_cfX_Triangular.R
#'
#' @export

cfX_Triangular <- function(t, a = -1, b = 1, c = 0) {
     szt <- dim(t)
     t <- c(t)

     cf <- -2*((b-c)*exp(1i*a*t) - (b-a)*exp(1i*c*t) + (c-a)*exp(1i*b*t))/((b-a)*(c-a)*(b-c)*t^2)
     cf[t == 0] <- 1

     dim(cf) <- szt

     return(cf)
}
