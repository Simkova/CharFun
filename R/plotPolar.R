#' Title
#'
#' @param angle angle
#' @param radius amplitude
#' @param title an overall title for the plot
#' @param type type of lines, dots
#' @param color color
#'
#' @return
#' @export
#'
#' @examples
#'
#'
plotPolar <- function(angle, radius, title = "", type = "l", color = "blue") {


  sgag
  x <- radius * cos(angle)
  y <- radius * sin(angle)

  m <- max(abs(radius*exp(1i*angle)))
  m <- ceiling(m)

  par(pty="s")

  plot(x = x, y = y,
       main = title,
       xlab = "",
       ylab = "",
       type = "l",
       col = "blue",
       xlim = c(-m, m),
       ylim = c(-m, m),
       asp = 1
  )

  n <- 12

  for (i in seq(0, n/2)) {
    a <- 2*pi/n
    lines(m*c(cos(i*a), -cos(i*a)), m*c(sin(i*a), -sin(i*a)),
          lty = 1,
          col = "gray")
  }

  a <- seq(0, 2*pi, length.out = 101)
  p <- 4
  m <- m/p

  for (i in seq(1:p)) {
    x <- m*i*cos(a)
    y <- m*i*sin(a)

    lines(
        x = x,
        y = y,
        type = "l",
        col = "gray"
      )
  }

  par(pty="m")
}
