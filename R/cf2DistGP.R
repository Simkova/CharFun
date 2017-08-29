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
#' @importFrom stats runif
#' @importFrom graphics plot grid
#'
#' @seealso For more details see:
#' \url{https://arxiv.org/pdf/1701.08299.pdf}
#'
#' @param cf      function handle for the characteristic function CF
#' @param x       vector of values from domain of the distribution F, if missing,
#'                cf2DistFFT automatically selects vector x from the domain
#' @param prob    vector of values from [0,1] for which the quantiles will be estimated
#' @param option  list with the following (default) parameters:
#' \itemize{
#'     \item option$isCompound   = FALSE   treat the compound distributions
#'     \item option$isCircular   = FALSE   treat the circular distributions
#'     \item option$isInterp     = FALSE   create and use the interpolant functions for PDF/CDF/RND
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
#' @param isCircular    treat the circular distributions, default FALSE
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
#'
#' @example R/Examples/example_cf2DistGP.R
#' @export
#'
#'
cf2DistGP <-
  function(cf,
           x,
           prob,
           option,
           isCompound,
           isCircular,
           N,
           SixSigmaRule,
           xMin,
           xMax,
           isPlot) {
    # check and set the default input arguments
    if (missing(option)) {
      option <- list()
    }

    if (!missing(isCompound)) {
      option$isCompound = isCompound
    } else if (!"isCompound" %in% names(option)) {
      option$isCompound = FALSE
    }
    if (!missing(isCircular)) {
      option$isCircular = isCircular
    } else if (!"isCircular" %in% names(option)) {
      option$isCircular = FALSE
    }
    if (!missing(N)) {
      option$N = N
    } else if (!"N" %in% names(option)) {
      if (option$isCompound) {
        option$N = 2 ^ 14
      }  else {
          option$N = 2 ^ 10
      }
    }
    if (!missing(xMin)) {
      option$xMin = xMin
    } else if (!"xMin" %in% names(option)) {
      if (option$isCompound) {
        option$xMin = 0
      }  else {
          option$xMin = -Inf
      }
    }
    if (!missing(xMax)) {
      option$xMax = xMax
    } else if (!"xMax" %in% names(option)) {
      option$xMax = Inf
    }
    if (!missing(SixSigmaRule)) {
      option$SixSigmaRule = SixSigmaRule
    } else if (!"SixSigmaRule" %in% names(option)) {
      if (option$isCompound) {
        option$SixSigmaRule = 10
      }  else {
          option$SixSigmaRule = 6
      }
    }

    if (!"tolDiff" %in% names(option)) {
      option$tolDiff = 1e-04
    }
    if (!"crit" %in% names(option)) {
      option$crit = 1e-12
    }
    if (!missing(isPlot)) {
      option$isPlot = isPlot
    } else if (!"isPlot" %in% names(option)) {
      option$isPlot = TRUE
    }

    # Other options parameters


    if (!"qf0" %in% names(option)) {
      option$qf0 = Re((cf(1e-4) - cf(-1e-4)) / (2e-4 * 1i))
    }
    if (!"maxiter" %in% names(option)) {
      option$maxiter = 1000
    }
    if (!"xN" %in% names(option)) {
      option$xN = 101
    }
    if (!"chebyPts" %in% names(option)) {
      option$chebyPts = 2^8
    }
    if (!"CorrectCDF" %in% names(option)) {
      if (option$isCircular) {
        option$CorrectCDF = TRUE
      } else {
        option$CorrectCDF = FALSE
      }
    }
    if (!"isInterp" %in% names(option)) {
      option$isInterp = FALSE
    }

    # First, set a special treatment if the real value of CF at infinity (large value)
    # is positive, i.e. const = real(cf(Inf)) > 0. In this case the
    # compound CDF has jump at 0 of size equal to this value, i.e. cdf(0) =
    # const, and pdf(0) = Inf. In order to simplify the calculations, here we
    # calculate PDF and CDF of a distribution given by transformed CF, i.e.
    # cf_new(t) = (cf(t)-const) / (1-const); which is converging to 0 at Inf,
    # i.e. cf_new(Inf) = 0. Using the transformed CF requires subsequent
    # recalculation of the computed CDF and PDF, in order to get the originaly
    # required values: Set pdf_original(0) =  Inf & pdf_original(x) = pdf_new(x) * (1-const),
    # for x > 0. Set cdf_original(x) =  const + cdf_new(x) * (1-const).


    const <- abs(Re(cf(1e+30)))
    if (option$isCompound) {
      cfOld <- cf
      if (const > 1e-13) {
        cf <- function(x)
          ((cfOld(x) - const) / (1 - const))
      }
    }

    if ("DIST" %in% names(option)) {
      # Set values from the last evaluation.

      xMean <- option$DIST$xMean
      cft <- option$DIST$cft
      xMin <- option$DIST$xMin
      xMax <- option$DIST$xMax
      N <- option$DIST$N
      k <- option$DIST$k
      xRange <- xMax - xMin
      dt <- 2 * pi / xRange
      t <- (1:N) * dt
      xStd <- NaN
    } else {
      N <- 2 * option$N
      dt <- option$dt
      T <- option$T
      xMin <- option$xMin
      xMax <- option$xMax
      xMean <- option$xMean
      xStd <- option$xStd
      SixSigmaRule <- option$SixSigmaRule
      tolDiff <- option$tolDiff
      cft <- cf(tolDiff * (1:4))


      # Evaluate values from input

      if (is.null(xMean)) {
        if (option$isCircular) {
          xMean <- Arg(cf(1))
        } else {
          xMean <-
            (-3*Im(cft[4]) + 32*Im(cft[3]) - 168*Im(cft[2]) + 672*Im(cft[1])) / (420 * tolDiff)
        }
      }
      if (is.null(xStd)) {
        if (option$isCircular) {
          xStd <- sqrt(-2 * log(abs(cf(1))))
        } else {
          xM2   = (-9*Re(cft[4]) + 128 * Re(cft[3]) - 1008 * Re(cft[2]) + 8064 * Re(cft[1]) - 7175) / (2520 * tolDiff^2)
          xStd  = sqrt(-xM2 - xMean^2)
        }
      }
      if (is.finite(xMin) && is.finite(xMax)) {
        xRange <- xMax - xMin
      } else if (!is.null(T)) {
        xRange <- 2 * pi / (T / N)
        if (is.finite(xMax)) {
          xMin = xMax - xRange
        } else if (is.finite(xMin)) {
          xMax = xMin + xRange
        } else {
          xMin = xMean - xRange / 2
          xMax = xMean + xRange / 2
        }
      } else if (!is.null(dt)) {
        xRange <- 2 * pi / dt
        if (is.finite(xMax)) {
          xMin = xMax - xRange
        } else if (is.finite(xMin)) {
          xMax = xMin + xRange
        } else {
          xMin = xMean - xRange / 2
          xMax = xMean + xRange / 2
        }
      } else {
        if (option$isCircular) {
          xMin <- -pi
          xMax <- pi
        } else {
          if (is.finite(xMin)) {
            xMax <- xMean + SixSigmaRule * xStd
          } else if (is.finite(xMax)) {
            xMin <- xMean - SixSigmaRule * xStd
          } else {
            xMin <- xMean - SixSigmaRule * xStd
            xMax <- xMean + SixSigmaRule * xStd
          }
        }
        xRange <- xMax - xMin
      }

      xRange  <- xRange
      dt      <- 2 * pi / xRange
      t       <- (1:N) * dt
      cft     <- cf(t)
      cft[N]    <- cft[N] / 2

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

    if (missing(x)) {
      x <- seq(xMin, xMax, length.out = option$xN)
    }

    if (option$isInterp) {
      xOrg <- x
      #Chebysev points
      x <- (xMax - xMin) * (-cos(pi * (0:option$chebyPts) / option$chebyPts) + 1) / 2 + xMin
    } else {
      xOrg <- c()
    }

    #WARNING: OUT of range

    if (any(x < xMin || any(x > xMax))) {
      warning(
        "CharFun: cf2DistGP" ,
        "x out of range (the used support): [xMin, xMax] = [",
        xMin,
        ", ",
        xMax,
        "]!"
      )
    }

    # Evaluate the required functions

    szx <- dim(x)
    x <- c(x)
    E <- exp((-1i) * x %*% t(t))

    # CDF estimate computed by using the simple trapezoidal quadrature rule

    cdf <- (xMean - x) / 2 + Im(E %*% (cft / t))
    cdf <- 0.5 - (cdf %*% dt) / pi

    # Correct the CDF (if the computed result is out of (0,1))
    # This is useful for circular distributions over intervals of length 2*pi,
    # as e.g. the von Mises distribution

    corrCDF <- 0
    if (option$CorrectCDF) {
      if (min(cdf) < 0) {
        corrCDF <- min(cdf)
        cdf <- cdf - corrCDF
      }
      if (max(cdf) > 1) {
        corrCDF <- max(cdf)-1
        cdf <- cdf - corrCDF
      }
    }

    dim(cdf) <- szx

    # PDF estimate computed by using the simple trapezoidal quadrature rule

    pdf <- 0.5 + Re(E %*% cft)
    pdf <- (pdf %*% dt) / pi
    pdf[pdf < 0] <- 0
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
      pdf = pdf * (1 - const)
      pdf[x == 0] = Inf
      pdf[x == xMax] = NA
    }

    # Calculate the precision criterion PrecisionCrit = abs(cf(t)/t) <= tol,
    # PrecisionCrit should be small for t > T, smaller than tolerance
    # options$crit

    PrecisionCrit = abs(cft[length(cft)] / t[length(t)])
    isPrecisionOK = (PrecisionCrit <= option$crit)

    # QF evaluated by the Newton-Raphson iterative scheme

    if (!missing(prob)) {
      isPlot = option$isPlot
      option$isPlot = FALSE
      szp <- dim(prob)
      maxiter <- option$maxiter
      crit <- option$crit
      qf <- option$qf0
      criterion <- TRUE
      count <- 0
      res <- cf2DistGP(cf, qf, option = option)
      cdfQ <- res$cdf
      pdfQ <- res$pdf


      while (criterion) {
        count <- count + 1
        correction <- (cdfQ - corrCDF - prob) / pdfQ
        qf <- pmax(xMin, pmin(xMax, qf - correction))

        res <- cf2DistGP(function(x)
          cf(x), x = qf, option = option)
        cdfQ <- res$cdf
        pdfQ <- res$pdf

        criterion <- any(abs(correction) > crit * abs(qf)) &&
          max(abs(correction)) > crit &&
          count < maxiter
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

    if (option$isInterp) {
      id <- is.finite(pdf)
      PDF <- function(xNew) pmax(0, interpBarycentric(x[id], pdf[id],xNew))

      id <- is.finite(cdf)
      CDF <- function(xNew) pmax(0, pmin(1, interpBarycentric(x[id]), cdf[id], xNew))

      QF <- function(prob) interpBarycentric(cdf[id], x[id], prob)

      RND <- function(n) QF(runif(n))
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

    if (option$isInterp) {
      result$PDF <- PDF
      result$CDF <- CDF
      result$QF <- QF
      result$RND <- RND
    }

    # PLOT the PDF / CDF

    if (length(x))
      option$isPot = FALSE

    if (option$isPlot) {
      plot(
        x = x,
        y = pdf,
        main = "PDF Specified by the Characteristic Function",
        xlab = "x",
        ylab = "pdf",
        type = "l",
        col = "blue"
      )
      grid()

      plot(
        x = x,
        y = cdf,
        main = "CDF Specified by the Characteristic Function",
        xlab = "x",
        ylab = "cdf",
        type = "l",
        col = "blue"
      )
      grid()
    }

    return(result)

  }
