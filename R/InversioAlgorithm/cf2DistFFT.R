#' @title
#' Evaluating CDF/PDF/QF from CF of a continous distribution F by using the FFT algorithm
#'
#' @description
#' \code{cf2DistFFT(cf, x, prob, options)} evaluates the CDF/PDF/QF (Cumulative Distribution Function,
#' Probability Density Function, Quantile Function) from the Characteristic Function CF
#' of a (continuous) distribution F by using the Fast Fourier Transform (FFT) algorithm.
#'
#' The algorithm cf2DistFFT evaluates the approximate values CDF(x), PDF(x),
#' and/or the quantiles QF(prob) for given x and prob, by interpolation from
#' the PDF-estimate computed by the numerical inversion of the given
#' characteristic function CF by using the FFT algorithm.
#'
#' @family Algorithm
#'
#' @seealso
#'
#' @param cf      function handle for the characteristic function CF
#' @param x       vector of values from domain of the distribution F, if missing,
#'                cf2DistFFT automatically selects vector x from the domain
#' @param prob    vector of values from [0,1] for which the quantiles will be estimated,
#'                default vector prob = c(0.9,0.95,0.975,0.99,0.995,0.999)
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
#' @param N             N points used by FFT, default 2^10 (2^14 for compound distribution)
#' @param SixSigmaRule  set the rule for computing domain, default 6
#' @param xMin          set the lower limit of X
#' @param xMax          set the upper limit of X
#' @param isPlot        plot the graphs of PDF/CDF, default TRUE
#'
#' @return result - list with the following results values:
#'    \item{result$x}{<- x}
#'    \item{result$cdf}{<- cdf}
#'    \item{result$pdf}{<- pdf}
#'    \item{result$prob}{<- prob}
#'    \item{result$qf}{<- qf}
#'    \item{result$xFFT}{<- xFFT}
#'    \item{result$pdfFFT}{<- pdfFFT}
#'    \item{result$cdfFFT}{<- cdfFFT}
#'    \item{result$SixSigmaRule}{<- option$SixSigmaRule}
#'    \item{result$N}{<- option$N}
#'    \item{result$dt}{<- dt}
#'    \item{result$T}{<- t[(length(t))]}
#'    \item{result$PrecisionCrit}{<- PrecisionCrit}
#'    \item{result$myPrecisionCrit}{<- options$crit}
#'    \item{result$isPrecisionOK}{<- isPrecisionOK}
#'    \item{result$xMean}{<- xMean}
#'    \item{result$xSTD}{<- xStd}
#'    \item{result$xMin}{<- xMin}
#'    \item{result$xMAx}{<- xMax}
#'    \item{result$cf}{<- cf}
#'    \item{result$option}{<- option}
#'    \item{result$tictoc}{<- proces time}
#'
#' @usage result <- cf2DistFFT(cf)
#' result <- cf2DistFFT(cf, x, prob, options)
#'
#' @note
#' The outputs of the algorithm cf2DistFFT are approximate values! The
#' precission of the presented results depends on several different
#' factors:
#' \itemize{
#'     \item  application of the FFT algorithm for numerical inversion of the CF
#'     \item selected number of points used by the FFT algorithm (by default
#'    option$N = 2^10),
#'     \item estimated/calculated domain [A,B] of the distribution F, used with the
#'    FFT algorithm. Optimally, [A,B] covers large part of the
#'    distribution domain, say more than 99%. However, the default
#'    automatic procedure for selection of the domain [A,B] may fail. It
#'    is based on the 'SixSigmaRule': \cr
#'    A = MEAN - SixSigmaRule * STD, and B = MEAN + SixSigmaRule * STD. \cr
#'    Alternatively, change the option$SixSigmaRule to different value, say 12,
#'    or use the option$xMin and option$xMax to set manually the values of A and B.
#'    }
#'
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
cf2DistFFT <- function(cf, x, prob, option, isCompound, N, SixSigmaRule, xMin, xMax, isPlot) {

# check and set the default input arguments
  if (missing(option)) {option <- list()}

  if (!missing(isCompound)) {option$isCompound = isCompound
  } else if (!"isCompound" %in% names(option)) {option$isCompound = FALSE}
  if (!missing(N)) {option$N = N
  } else if (!"N" %in% names(option)) {
    if (option$isCompound) {option$N = 2 ^ 14
    } else {option$N = 2 ^ 10}
  }
  if (!missing(xMin)) {option$xMin = xMin
  } else if (!"xMin" %in% names(option)) {
    if (option$isCompound) {option$xMin = 0
    } else {option$xMin = -Inf}
  }
  if (!missing(xMax)) {option$xMax = xMax
  } else if (!"xMax" %in% names(option)) {option$xMax = Inf}
  if (!missing(SixSigmaRule)) {
    option$SixSigmaRule =SixSigmaRule
  } else if (!"SixSigmaRule" %in% names(option)) {
    if (option$isCompound) {option$SixSigmaRule = 10} else {option$SixSigmaRule = 6}}
  if (!"tolDiff" %in% names(option)) {option$tolDiff = 1e-04}
  if (!"crit" %in% names(option)) {option$crit = 1e-12}
  if (!missing(isPlot)) {option$isPlot = isPlot
  } else if (!"isPlot" %in% names(option)) {option$isPlot = TRUE}

# Other options parameters

  if (!"isPlotFFT" %in% names(option)) {option$isPlotFFT = FALSE}
  if (!"xN" %in% names(option)) {option$xN = 101}

# GET/SET the DEFAULT parameters and the OPTION

# Set a special treatment if the real value of CF at infinity (large value) is positive,
# i.e. const = real(cf(Inf)) > 0. In this case the compound CDF has jump at 0 of size equal to this value,
# i.e. cdf(0) = const, and pdf(0) = Inf. In order to simplify the calculations,
# here we calculate PDF and CDF of a distribution given by transformed CF, i.e.
# cf_new(t) = (cf(t)-const) / (1-const); which is converging to 0 at Inf,
# i.e. cf_new(Inf) = 0. Using the transformed CF requires subsequent
# recalculation of the computed CDF and PDF, in order to get the originaly required values:
# Set pdf_original(0) =  Inf & pdf_original(x) = pdf_new(x) * (1-const),
# for x > 0. Set cdf_original(x) =  const + cdf_new(x) * (1-const).

  const <- Re(cf(1e+30))
  if (option$isCompound) {
    cfOld <- cf
    if (const > 1e-13) {
      cf <- function(x) ((cf(x) - const) / (1 - const))
    }
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
    xStd <- c()
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

    if (is.null(xMean)) {
      xMean <-
        Re((-cft[2] + 8 * cft[1] - 8 * Conj(cft[1]) + Conj(cft[2])) / (1i * 12 *
                                                                         tolDiff))
    }
    if (is.null(xStd)) {
      xM2 <-
        Re(
          -(
            Conj(cft[4]) - 16 * Conj(cft[3]) + 64 * Conj(cft[2]) + 16 * Conj(cft[1]) - 130 + 16 *
              Conj(cft[1]) + 64 * cft[2] - 16 * cft[3] + cft[4]
          ) / (144 * tolDiff ^ 2)
        )
      xStd <- sqrt(xM2 - xMean ^ 2)
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
    dt      <- 2 * pi / xRange
    k       <- (0:(N - 1))
    t       <- (k - N / 2 + 0.5) * dt
    cft     <- cf(t[((N / 2) + 1):length(t)])
    cft     <- c(Conj(cft[c(length(cft):1)]), cft)
    cft[1]  <- cft[1] / 2
    cft[length(cft)]  <- cft[length(cft)] / 2

    option$DIST$xMean   <- xMean
    option$DIST$xMin    <- xMin
    option$DIST$xMax    <- xMax
    option$DIST$xRange  <- xRange
    option$DIST$cft     <- cft
    option$DIST$N       <- N
    option$DIST$k       <- k
    option$DIST$dt      <- dt
    option$DIST$t       <- t
  }

# ALGORITHM ---------------------------------------------------------------

  A      <- xMin
  B      <- xMax
  dx     <- (B - A) / N
  c      <- ((-1 + 0i) ^ (A * (N - 1) / (B - A))) / (B - A)
  C      <- c * (-1 + 0i) ^ ((1 - 1 / N) * k)
  D      <- (-1 + 0i) ^ (-2 * (A / (B - A)) * k)
  pdfFFT <- pmax(0, Re(C * fft(D * cft)))
  cdfFFT <- pmin(1, pmax(0, 0.5 + Re(1i * C * fft(D * cft / t))))
  xFFT   <- A + k * dx

# Reset the transformed CF, PDF, and CDF to the original values

  if (option$isCompound) {
    cf <- cfOld
    cdfFFT <- const + cdfFFT * (1 - const)
    pdfFFT <- pdfFFT * (1 - const)
    pdfFFT[x == 0] <- Inf
  }

# Calculate the precision criterion PrecisionCrit = abs(cf(t)/t) <= tol,
# PrecisionCrit should be small for t > T, smaller than tolerance options.crit

  PrecisionCrit = abs(cft[length(cft)]/t[length(t)])
  isPrecisionOK = (PrecisionCrit<=option$crit)

  if (missing(prob)) prob = c(0.9,0.95,0.975,0.99,0.995,0.999)

  cdfu <- unique(cdfFFT)
  id   <- !duplicated(cdfFFT)
  xxU  <- xFFT[id]
  szp  <- dim(prob)

  cdfuxxU <-
    matrix(c(cdfu, xxU),
           nrow  = 2,
           ncol  = length(cdfu),
           byrow = TRUE)
  cdfuxxU <- cdfuxxU[, order(cdfuxxU[1, ], cdfuxxU[2, ])]

  qfFun <- function(prob)
    pchip(xi =  (cdfuxxU[1,]),
          yi =  c(cdfuxxU[2,]) + dx / 2,
          x =  prob)
  qf        <- qfFun(prob)
  dim(qf)   <- szp
  dim(prob) <- szp

# INTERPOLATE CDF required values: CDF(x)

  if (missing(x)) x <- {seq(xMin, xMax, length.out = option$xN)}
  szx <- dim(x)
  x <- c(x)
  cdfFun <- function(x) approx(x = xxU, y = cdfu,xout = x)
  cdf <- matrix(unlist(cdfFun(x)), nrow = 2, byrow = TRUE)[2, ]
  dim(cdf) <- szx
  dim(x) <- szx

# TRY INTERPOLATE PDF required values: PDF(x)

  pdfFun = function(x) approx(xFFT,pdfFFT,x);

  pdf <- try(matrix(unlist(pdfFun(x)), nrow = 2, byrow = TRUE)[2,])*seq(from = 1,to = 1,length.out = length(x))
  pdf <-
    tryCatch(
      expr = matrix(unlist(pdfFun(x)), nrow = 2, byrow = TRUE)[2, ]
    ) * seq(from = 1, to = 1, length.out = length(x))

  dim(x) <- szx

# RESULT

  result <- list(
    "x"                  <- x,
    "cdf"                <- cdf,
    "pdf"                <- pdf,
    "prob"               <- prob,
    "qf"                 <- qf,
    "xFFT"               <- xFFT,
    "pdfFFT"             <- pdfFFT,
    "cdfFFT"             <- cdfFFT,
    "SixSigmaRule"       <- option$SixSigmaRule,
    "N"                  <- N,
    "dt"                 <- dt,
    "T"                  <- t[length(t)],
    "PrecisionCrit"      <- PrecisionCrit,
    "myPrecisionCrit"    <- option$crit,
    "isPrecisionOK"      <- isPrecisionOK,
    "xMean"              <- xMean,
    "xStd"               <- xStd,
    "xMin"               <- xMin,
    "xMax"               <- xMax,
    "cf"                 <- cf,
    "isCompound"         <- option$isCompound,
    "option"             <- option,
    "tictoc"             <- proc.time()
  )

  if (length(x) == 1) option$isPlot = FALSE

  if (option$isPlotFFT) {
    x <- xFFT
    pdf <- pdfFFT
    cdf <- cdfFFT
  }

# PLOT THE PDF/CDF, if required

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
