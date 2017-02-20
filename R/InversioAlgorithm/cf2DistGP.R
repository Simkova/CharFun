#' @title
#' Evaluating CDF/PDF/QF from CF of a continous distribution F by using the Gil-Pelaez inversion formulae.
#'
#'
#' @description
#' \code{cf2DistGP(cf, x, prob, options)} evaluates the CDF/PDF/QF (Cumulative Distribution Function,
#' Probability Density Function, Quantile Function) from the Characteristic Function CF
#' by using the Gil-Pelaez inversion formulae: \cr
#' cdf(x) = 1/2 + (1/pi) * Integral_0^inf imag(exp(-1i*t*x)*cf(t)/t)*dt. \cr
#' pdf(x) = (1/pi) * Integral_0^inf real(exp(-1i*t*x)*cf(t))*dt.
#'
#' @family Algorithm
#'
#' @seealso
#'
#' @param cf      function handle for the characteristic function CF
#' @param x       vector of values from domain of the distribution F, if missing,
#'                cf2DistFFT automatically selects vector x from the domain
#' @param prob    vector of values from [0,1] for which the quantiles will be estimated
#' @param option  list with the following (default) parameters:
#' \itemize{
#'     \item option$isCompound   = FALSE   treat the compound distributions
#'     \item option$N            = 2^10    set N points used by FFT
#'     \item option$xMin         = -Inf    set the lower limit of X
#'     \item option$xMax         = Inf     set the upper limit of X
#'     \item option$xMean        = NULL    set the MEAN value of X
#'     \item option$xStd         = NULL    set the STD value of X
#'     \item option$dt           = NULL    set grid step dt = 2*pi/xRange
#'     \item option$T            = NULL    set upper limit of (0,T), T = N*dt
#'     \item option$SixSigmaRule = 6       set the rule for computing domain
#'     \item option$tolDiff      = 1e-4    tol for numerical differentiation
#'     \item option$isPlot       = TRUE    plot the graphs of PDF/CDF
#'
#'     \item option$DIST                   list with information for future evaluations.
#'                                         option$DIST is created automatically after first call:
#'     \itemize{
#'         \item option$DIST$xMean
#'         \item option$DIST$xMin
#'         \item option$DIST$xMax
#'         \item option$DIST$xRange
#'         \item option$DIST$cft
#'         \item option$DIST$N
#'         \item option$DIST$k
#'         \item option$DIST$t
#'         \item option$DIST$dt
#'         }
#'     }
#'
#' @param isCompound    treat the compound distributions, default FALSE
#' @param N             N points used by GP, default 2^10 (2^14 for compound distribution)
#' @param SixSigmaRule  set the rule for computing domain, default 6
#' @param xMin          set the lower limit of X
#' @param xMax          set the upper limit of X
#' @param isPlot        plot the graphs of PDF/CDF, default TRUE
#'
#' @return result - list with the following results values:
#'    \item{result$x}{<- x}
#'    \item{result$cdf}{<- cdf (Cumulative Distribution Function)}
#'    \item{result$pdf}{<- pdf (Probability Density Function)}
#'    \item{result$prob}{<- prob}
#'    \item{result$qf}{<- qf (Quantile Function)}
#'    \item{result$SixSigmaRule}{<- option$SixSigmaRule}
#'    \item{result$N}{<- option$N}
#'    \item{result$dt}{<- dt}
#'    \item{result$T}{<- t[(length(t))]}
#'    \item{result$PrecisionCrit}{<- PrecisionCrit}
#'    \item{result$myPrecisionCrit}{<- options$crit}
#'    \item{result$isPrecisionOK}{<- isPrecisionOK}
#'    \item{result$xMean}{<- xMean}
#'    \item{result$xRange}{<- xRange}
#'    \item{result$xSTD}{<- xStd}
#'    \item{result$xMin}{<- xMin}
#'    \item{result$xMAx}{<- xMax}
#'    \item{result$cf}{<- cf}
#'    \item{result$const}{<- const}
#'    \item{result$isCompound}{<- option$isCompound}
#'    \item{result$details$count}{<- count}
#'    \item{result$details$correction}{<- correction}
#'    \item{result$option}{<- option}
#'    \item{result$tictoc}{<- toc}
#'
#' @usage result <- cf2DistGP(cf)
#' result <- cf2DistGP(cf, x, prob, options)
#'
#' @note
#' If options.DIST is provided, then cf2DistGP evaluates CDF/PDF based on
#' this information only (it is useful, e.g., for subsequent evaluation of
#' the quantiles). options.DIST is created automatically after first call.
#' This is supposed to be much faster, bacause there is no need for further
#' evaluations of the characteristic function. In fact, in such case the
#' function handle of the CF is not required, i.e. in such case set cf = []
#'
#' The required integrals are evaluated approximately by using the simple
#' trapezoidal rule on the interval(0,T), where T = N * dt is a sufficienly
#' large integration upper limit in the frequency domain.
#'
#' If the optimum values of N and T are unknown, we suggest, as a simple
#' rule of thumb, to start with the application of the six-sigma-rule for
#' determining the value of dt = (2*pi)/(xMax-xMin), where xMax = xMean +
#' 6*xStd, and xMin = xMean - 6*xStd, see \code{[1]}.
#'
#' Please note that THIS (TRAPEZOIDAL) METHOD IS AN APPROXIMATE METHOD:
#' Frequently, with relatively low numerical precision of the results of the
#' calculated PDF/CDF/QF, but frequently more efficient and more precise
#' than comparable Monte Carlo methods.
#'
#' However, the numerical error (truncation error and/or the integration
#' error) could be and should be properly controled!
#'
#' CONTROLING THE PRECISION:
#' Simple criterion for controling numerical precision is as follows: Set N
#' and T = N*dt such that the value of the integrand function
#' Imag(e^(-1i*t*x) * cf(t)/t) is sufficiently small for all t > T, i.e.
#' PrecisionCrit = abs(cf(t)/t) <= tol,
#' for pre-selected small tolerance value, say tol = 10^-8. If this
#' criterion is not valid, the numerical precission of the result is
#' violated, and the method should be improved (e.g. by selecting larger N
#' or considering other more sofisticated algorithm - not considered here).
#' For more details consult the references below.
#'
#' @references
#' \enumerate{
#'     \item WITKOVSKY, V.: On the exact computation of the density and of
#'     the quantiles of linear combinations of t and F random
#'     variables. Journal of Statistical Planning and Inference 94
#'     (2001), 113.
#'     \item WITKOVSKY, V.: Matlab algorithm TDIST: The distribution of a
#'     linear combination of Students t random variables. In COMPSTAT
#'     2004 Symposium (2004), J. Antoch, Ed., Physica-Verlag/Springer
#'     2004, Heidelberg, Germany, pp. 19952002.
#'     \item WITKOVSKY, V.: WIMMER,G., DUBY, T. Logarithmic Lambert W x F
#'     random variables for the family of chi-squared distributions
#'     and their applications. Statistics & Probability Letters 96
#'     (2015), 223231.
#'     \item WITKOVSKY V. (2016). Numerical inversion of a characteristic
#'     function: An alternative tool to form the probability distribution of
#'     output quantity in linear measurement models. Acta IMEKO, 5(3), 32-44.
#'     \item WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
#'     distribution based on numerical inversion of the compound empirical
#'     characteristic function of frequency and severity. Preprint submitted
#'     to Insurance: Mathematics and Economics.
#'     \item WITKOVSKY V., WIMMER G., DUBY T. (2016). Computing the aggregate loss
#'     distribution based on numerical inversion of the compound empirical
#'     characteristic function of frequency and severity. Preprint submitted
#'     to Insurance: Mathematics and Economics.
#'     \item DUBY T., WIMMER G., WITKOVSKY V.(2016). MATLAB toolbox CRM for
#'     computing distributions of collective risk models. Preprint submitted
#'     to Journal of Statistical Software.
#'     }
#'
#' @export
#'
#'
cf2DistGP <- function(cf, x, prob, option, isCompound, N, SixSigmaRule, xMin, xMax, isPlot) {

# check and set the default input arguments
  if (missing(option)) option <- list()

  if (!missing(isCompound)) {
    option$isCompound = isCompound
  } else if (!"isCompound" %in% names(option)) option$isCompound = FALSE
  if (!missing(N)) {
    option$N = N
  } else if (!"N" %in% names(option)) {
    if (option$isCompound) option$N = 2^14 else option$N = 2^10}
  if (!missing(xMin)) {
    option$xMin = xMin
  } else if (!"xMin" %in% names(option)) {
    if (option$isCompound) option$xMin = 0 else option$xMin = -Inf}
  if (!missing(xMax)) {
    option$xMax = xMax
  } else if (!"xMax" %in% names(option)) {option$xMax = Inf}

#  if ("xMean" %in% names(argg)) option$xMean = argg$xMean
#  if ("xStd" %in% names(argg)) option$xStd = argg$xStd
#  if ("dt" %in% names(argg)) option$dt = argg$dt
#  if ("T" %in% names(argg)) option$T = argg$T
#  if (!"xMean" %in% names(option)) {option$xMean = NULL}
#  if (!"xStd" %in% names(option)) {option$xStd = NULL}
#  if (!"dt" %in% names(option)) {option$dt = NULL}
#  if (!"T" %in% names(option)) {option$T = NULL}

  if (!missing(SixSigmaRule)) {
    option$SixSigmaRule =SixSigmaRule
  } else if (!"SixSigmaRule" %in% names(option)) {
    if (option$isCompound) option$SixSigmaRule = 10 else option$SixSigmaRule = 6}

  if (!"tolDiff" %in% names(option)) {option$tolDiff = 1e-04}
  if (!"crit" %in% names(option)) {option$crit = 1e-12}
  if (!missing(isPlot)) {
    option$isPlot = isPlot
  } else if (!"isPlot" %in% names(option)) {option$isPlot = TRUE}

# Other options parameters

  if (!"qf0" %in% names(option)) {option$qf0 = (cf(1e-4)-cf(-1e-4))/(2e-4*1i)}
  if (!"maxiter" %in% names(option)) {option$maxiter = 1000}
  if (!"xN" %in% names(option)) {option$xN = 101}

  const <- Re(cf(1e+30))
  if (option$isCompound) {
    cfOld <- cf
    if (const > 1e-13) cf <- function(x) ((cf(x) - const)/(1 - const))
  }

  if ("DIST" %in% names(option)) {
    xMean <- option$DIST$xMean
    cft <- option$DIST$cft
    xMin <- option$DIST$xMin
    xMax <- option$DIST$xMax
    N <- option$DIST$N
    k <- option$DIST$k
    xRange <- option$DIST$xRange
    dt <- option$DIST$dt
    t <- option$DIST$t
    xStd <- NaN
  } else {
    N <- 2*option$N
    dt <- option$dt
    T <- option$T
    xMin <- option$xMin
    xMax <- option$xMax
    xMean <- option$xMean
    xStd <- option$xStd
    SixSigmaRule <- option$SixSigmaRule
    tolDiff <- option$tolDiff
    cft <- cf(tolDiff*(1:4))

    if (is.null(xMean)) {
      xMean <- Re((-cft[2] + 8*cft[1]-8*Conj(cft[1]) + Conj(cft[2]))/(1i*12*tolDiff))
    }
    if (is.null(xStd)) {
      xM2 <- Re(-(Conj(cft[4]) - 16*Conj(cft[3]) + 64*Conj(cft[2]) + 16*Conj(cft[1]) - 130 + 16*Conj(cft[1]) + 64*cft[2] - 16*cft[3] + cft[4]) / (144*tolDiff^2))
      xStd <- sqrt(xM2 - xMean^2)
    }
    if (is.finite(xMin) && is.finite(xMax)) {
      xRange <- xMax - xMin
    } else if (!is.null(T)) {
      xRange <- 2*pi / (T/N)
      if (is.finite(xMax)) {
        xMin = xMax - xRange
      } else if (is.finite(xMin)) {
        xMax = xMin + xRange
      } else {
        xMin = xMean - xRange/2
        xMax = xMean + xRange/2
      }
    } else if (!is.null(dt)) {
        xRange <- 2*pi / dt
        if (is.finite(xMax)) {
          xMin = xMax - xRange
        } else if (is.finite(xMin)) {
          xMax = xMin + xRange
        } else {
          xMin = xMean - xRange/2
          xMax = xMean + xRange/2
        }
    } else {
      if (is.finite(xMin)) {
        xMax <- xMean + SixSigmaRule * xStd
      } else if (is.finite(xMax)) {
        xMin <- xMean - SixSigmaRule * xStd
      } else {
        xMin <- xMean - SixSigmaRule * xStd
        xMax <- xMean + SixSigmaRule * xStd
      }
      xRange <- xMax - xMin
    }

    xRange  <- xRange
    dt      <- 2*pi / xRange
    t       <- (1:N) * dt
    cft     <- cf(t)
    cft[N]    <- cft[N]/2

    option$DIST$xMean   <- xMean
    option$DIST$cft     <- cft
    option$DIST$xMin    <- xMin
    option$DIST$xMax    <- xMax
    option$DIST$N       <- N
    option$DIST$xRange  <- xRange
    option$DIST$dt      <- dt
    option$DIST$t       <- t
  }


# ALGORITHM ---------------------------------------------------------------

  if (missing(x)) x <- seq(xMin, xMax, length.out = option$xN)

  #PRIDAT WARNING?

# Evaluate the required functions

  szx <- dim(x)
  x <- c(x)
  E <- exp((-1i)*x%*%t(t))

# CDF estimate computed by using the simple trapezoidal quadrature rule

  cdf <- (xMean - x)/2 + Im(E %*% (cft / t))
  cdf <- 0.5 - (cdf %*% dt) / pi
  cdf <-pmax(0,pmin(1,Re(cdf)))
  dim(cdf) <- szx

# PDF estimate computed by using the simple trapezoidal quadrature rule

  pdf <- 0.5 + Re(E %*% cft)
  pdf <- (pdf %*% dt) / pi
  pdf <- pmax(0, pdf)
  dim(pdf) <- szx
  dim(x) <- szx

# REMARK:
# Note that, exp(-1i*x_i*0) = cos(x_i*0) + 1i*sin(x_i*0) = 1. Moreover,
# cf(0) = 1 and lim_{t -> 0} cf(t)/t = E(X) - x. Hence, the leading term of
# the trapezoidal rule for computing the CDF integral is CDFfun_1 = (xMean- x)/2,
# and PDFfun_1 = 1/2 for the PDF integral, respectively.

# Reset the transformed CF, PDF, and CDF to the original values

  if (option$isCompound) {
    cf  <- cfOld
    cdf <- const + cdf * (1 - const)
    pdf = pdf * (1-const)
    pdf[x == 0] = Inf
    pdf[x == xMax] = NA
  }

# Calculate the precision criterion PrecisionCrit = abs(cf(t)/t) <= tol,
# PrecisionCrit should be small for t > T, smaller than tolerance
# options$crit

  PrecisionCrit = abs(cft[length(cft)]/t[length(t)])
  isPrecisionOK = (PrecisionCrit <= option$crit)

# QF evaluated by the Newton-Raphson iterative scheme

  if (!missing(prob)) {
    isPlot = option$isPlot
    option$isPlot = FALSE
    szp <- dim(prob)
    maxiter <- option$maxiter
    crit <- option$crit
    qf <- Re(option$qf0)
    criterion <- TRUE
    count <- 0
    res <- cf2DistGP(cf, qf, option = option)
    cdfQ <- res$cdf
    pdfQ <- res$pdf


    while (criterion) {
      count <- count + 1
      correction = (cdfQ - prob) / pdfQ
      qf = pmax(xMin, qf - correction)

      res <- cf2DistGP(function(x) cf(x), x = qf, option = option)
      cdfQ <- res$cdf
      pdfQ <- res$pdf

      criterion <- any(abs(correction) > crit * abs(qf)) && max(abs(correction)) > crit && count < maxiter
    }

    dim(qf) <- szp
    dim(prob) <- szp
    option$isPlot <- isPlot

  } else {
    qf <- c()
    count = c()
    correction = c()
    prob = c()
  }

  result <- list(
    "x"                  = x,
    "cdf"                = cdf,
    "pdf"                = pdf,
    "prob"               = prob,
    "qf"                 = qf,
    "SixSigmaRule"       = option$SixSigmaRule,
    "N"                  = N,
    "dt"                 = dt,
    "T"                  = t[length(t)],
    "PrecisionCrit"      = PrecisionCrit,
    "myPrecisionCrit"    = option$crit,
    "isPrecisionOK"      = isPrecisionOK,
    "xRange"             = xRange,
    "xMean"              = xMean,
    "xStd"               = xStd,
    "xMin"               = xMin,
    "xMax"               = xMax,
    "cf"                 = cf,
    "const"              = const,
    "isCompound"         = option$isCompound,
    "details$count"      = count,
    "details$correction" = correction,
    "option"             = option
  )

# PLOT the PDF / CDF

  if (length(x)) option$isPot = FALSE

  if (option$isPlot) {
    plot(x = x, y = pdf,
         main="PDF Specified by the Characteristic Function CF",
         xlab="x",
         ylab="pdf",
         type="l",
         col="blue")

    plot(x = x, y = cdf,
         main="CDF Specified by the Characteristic Function CF",
         xlab="x",
         ylab="cdf",
         type="l",
         col="red")
  }

  return(result)

}
