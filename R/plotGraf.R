#' @title Plot Graf f(x)
#'
#' @description
#' plotGraf(f, x, xmin, xmax) plots graf of complex function f in real arguments x or in range [xmin, xmax] with step dx
#'
#' @param f function
#' @param x numerical values (number or vector)
#' @param xmin number, default value xmin = -1
#' @param xmax number, default value xmax = 1
#' @param dx positive number, default value dx = 0.1
#' @param px numerical values (number or vector) in which are auxiliary line on axis x render default value \eqn{px = 0}
#' @param py numerical values (number or vector) in which are auxiliary line on axis y render default value \eqn{py = 0}
#' @param labelx rendering auxiliary line default value TRUE
#' @param labely rendering auxiliary line default value TRUE
#'
#' @details
#' function with optional parameters must be pass as lambda function f = function(x) f(x, optional parameters)
#'
#'
#' @returng Graf of function f in argument x or in range [xmin, xmax] with step dx
#' @export

plotGraf <- function(f, x, xmin, xmax, dx = 0.1, px = 0, py = 0, labelx = TRUE, labely = TRUE) {

  if (missing(x)) {
    if (missing(xmin) && missing(xmax)) {
      xmin = -1
      xmax = 1
    } else if (missing(xmin)) {
      xmin <- xmax - 2
    } else if (missing(xmin)) {
        xmax <- xmin + 2
    }
    x <- seq(xmin, xmax, by=dx)
  } else {
    if (missing(xmin)) {xmin <- min(x)}
    if (missing(xmax)) {xmax <- max(x)}
    x <- x[x >= xmin && x <= xmax]
  }

  real = Re(f(x))
  imag = Im(f(x))

  plot(x,real,
       main = "Function f(x)",
       xlab = "x",
       ylab = "f(x)",
       type = "l",
       col = "blue",
       xlim = c(xmin, xmax),
       ylim = c(min(real, imag), max(real, imag)))
  lines(x, imag,
        lty = 1,
        col = "red")

  if (labelx) {
    for (xi in c(px)) {
      lines(c(xi, xi), c(min(real, imag) - 1, max(real, imag) + 1),
            lty = 1,
            col = "gray")
    }
  }
  if (labely) {
    for (yi in c(py)) {
      lines(c(xmin, xmax), c(yi, yi),
            lty = 1,
            col = "gray")
    }
  }

  legend("bottomright",
         c("Re(f))","Im(f)"),
         fill=c("blue","red")
  )
}
