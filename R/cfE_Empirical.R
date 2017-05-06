#' @title Characteristic function of Empirical distribution
#'
#' @description
#' cfE_Empirical(t, data, cfX) evaluates the characteristic function cf(t)
#' of the Empirical distribution, based on the observed data. In particular,
#'   cf(t) = cfE_Empirical(t, data)
#'         = \deqn{(1/n) * sum_{j=1}^n cf_Dirac(data(j)*t)},
#' where cfX is function handle of the characteristic function cfX(t) of the
#' random variable X (as e.g. another empirical CF based on observed data of X).
#'
#' @family Empirical Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Empirical_characteristic_function}
#'
#' @param t numerical values (number, vector...)
#' @param data set of observed data
#' @param cfX function
#'
#' @return characteristic function cf(t) of the Empirical distribution, based on the observed data
#' @usage cfE_Empirical(t, data)
#'cfE_Empirical(t, data, cfX)
#'
#' @example Examples/example_cfE_Empirical.R
#'
#' @export

cfE_Empirical <- function(t, data = 1, cfX) {
  weights <- 1 / length(data)
  data <- c(data)


  # Special treatment for mixtures with large number of variables

  szcoefs <- dim(t(data))
  szcoefs <- szcoefs[1] * szcoefs[2]

  szt <- dim(t)
  sz <- dim(t(t))[1] * dim(t(t))[2]

  szcLimit <- ceiling(1e3 / (sz / 2 ^ 16))
  idc <- (1:(trunc(szcoefs / szcLimit) + 1))

  # Characteristic function of a weighted mixture of Dirac variables

  t <- c(t)
  idx0 <- 1
  cf <- 0

  for (j in idc) {
    idx1 <- min(idc[j] * szcLimit, szcoefs)
    idx <- (idx0:idx1)
    idx0 <- idx1 + 1

    if (missing(cfX)) {
      aux = exp(1i * t(t(t)) %*% data[idx])
    } else {
      aux = t(apply(t(cfX(t)), 2, FUN = '^', data[idx]))
    }

    cf = cf + apply(weights * aux, 1, sum)
  }

  dim(cf) <- szt

  return(cf)
}
