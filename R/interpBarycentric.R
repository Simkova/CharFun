#' @title
#' Barycentric Interpolation
#'
#' @description
#' \code{interpBarycentric} evaluates (interpolates) function values funNew at
#' given points xNew by barycentric interpolation from function values fun
#' given at chebpoints x.
#'
#' @family
#'
#' @seealso
#' \url{https://en.wikipedia.org/wiki/Barycentric_coordinate_system}
#'
#' @param x points in which is fun given (for more accuracy use chebpoints)
#' @param fun function values of fun given at points x
#' @param xNew point in which fun will be estimated
#' @param isChebPts default value TRUE
#'
#' @return function fun evaluated in points xNew
#'
#' @usage funNew <- interpBarycentric(x, fun, xNew)
#'
#' @export
interpBarycentric <- function(x, fun, xNew, isChebPts = TRUE) {
  if (missing(xNew)) {xNew <- seq(min(x), max(x))}

  x <- c(x)
  fun <- c(fun)

  sztNew <- dim(xNew)
  xNew <- c(xNew)
  nx <- length(x)
  nxNew <- length(xNew)
  funNew <- seq(0, 0, length.out = length(xNew))

  w <- (-1)^seq(0, length(x) - 1)
  w[1] <- w[1]/2
  w[length(x)] <- w[length(x)]/2

  for (i in seq(1, length(xNew))) {
    A <- 0
    B <- 0
    for (j in seq(1, length(x))) {
      if (xNew[i] == x[j]) {
        exactVal <- TRUE
        funNew[i] <- fun[j]
      } else {
        exactVal <- FALSE
        weight <- w[j] / (xNew[i] - x[j])
        A <- A + weight * fun[j]
        B <- B + weight
      }

      if (exactVal) {break}
      else {funNew[i] <- A / B}
    }

  }

  dim(xNew) <- sztNew
  dim(funNew) <- sztNew

  return(funNew)
}
