% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cfN_NegativeBinomial.R
\name{cfN_NegativeBinomial}
\alias{cfN_NegativeBinomial}
\title{Characteristic function of Negative-Binomial distribution}
\usage{
cfN_NegativeBinomial(t, r = 10, p = 0.5, cfX)
}
\arguments{
\item{t}{numerical values (number, vector...)}

\item{r}{number of trials}

\item{p}{success probability, \eqn{0 \le p \le 1}, default value p = 1/2}

\item{cfX}{function}
}
\value{
characteristic function cf(t) of the Negative-Binomial distribution with n trials and p success probability
}
\description{
cfN_NegativeBinomial(t, r, p) evaluates the characteristic function cf(t) of the
Negative-Binomial distribution with the parameters r (number of failures until
the experiment is stopped, r in N) and p (success probability in each experiment, p in [0,1]), i.e.
cfN_NegativeBinomial(t, r, p) = p^r * (1 - (1-p) * e^(1i*t))^(-r)

cfN_NegativeBinomial(t, r, p, cfX) evaluates the compound characteristic function
cf(t) = cfN_NegativeBinomial(-1i*log(cfX(t)), r, p), where cfX is function
handle of the characteristic function cfX(t) of a continuous distribution
and/or random variable X.

Note that such CF is characteristic function of the compound distribution,
i.e. distribution of the random variable Y = X_1 + ... + X_N, where X_i ~ F_X
are i.i.d. random variables with common CF cfX(t), and N ~ F_N is
independent RV with its CF given by cfN(t).
}
\examples{
## EXAMPLE1 (CF of the Negative Binomial distribution with r = 5, p = 0.3)
r <- 5
p <- 0.3
t <- seq(-15, 15, length.out = 1001)
plotGraf(function(t)
  cfN_NegativeBinomial(t, r, p), t,
  title = "CF of the Negative Binomial distribution with r = 5, p = 0.3")

## EXAMPLE2 (PDF/CDF of the compound NegativeBinomial-Exponential distribution)
r <- 5
p <- 0.3
lambda <- 5
cfX <- function(t)
  cfX_Exponential(t, lambda)
cf <- function(t)
  cfN_NegativeBinomial(t, r, p, cfX)
x <- seq(0, 10, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
result <- cf2DistGP(cf, x, prob, isCompound = TRUE)
}
\seealso{
For more details see WIKIPEDIA:
\url{https://en.wikipedia.org/wiki/Negative_binomial_distribution}

Other Discrete Probability Distribution: \code{\link{cfN_Binomial}},
  \code{\link{cfN_Delaporte}},
  \code{\link{cfN_GeneralizedPoisson}},
  \code{\link{cfN_Geometric}},
  \code{\link{cfN_Logarithmic}}, \code{\link{cfN_Poisson}},
  \code{\link{cfN_PolyaEggenberger}}
}

