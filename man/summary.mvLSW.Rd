% File src/library/base/man/summary.mvLSW.Rd
\name{summary.mvLSW}
\alias{summary.mvLSW}
\alias{print.mvLSW}
\title{Print a Summary of mvLSW Object}
\description{
  Prints a summary of the information contained within a \code{mvLSW} 
  classed object.
}
\usage{
  \method{summary}{mvLSW}(object, ...)
  \method{print}{mvLSW}(x, ...)
}

\arguments{
  \item{object, x}{A \code{mvLSW} object.}
  \item{...}{Additional arguments.}
}

\details{
  The command prints to screen a summary of the information 
  contained within a \code{mvLSW} object. Information printed 
  includes: dimensions, wavelet function, the smoothing regime 
  applied, smoothing kernel(s), generalized cross-validation 
  gamma deviance criteria score, application of the bias correction 
  and minimum eigenvalue from across all spectral matrices.
}

\value{
  This command returns nothing, only prints a summary to the console.
}

\references{
  Taylor, S.A.C., Park, T.A. and Eckley, I. (2019) Multivariate 
  locally stationary wavelet analysis with the mvLSW R package.
  \emph{Journal of statistical software} \strong{90}(11) pp. 1--16,
  doi: 10.18637/jss.v090.i11.
}

\seealso{
  \code{\link{mvEWS}}, \code{\link{as.mvLSW}}.
}

\examples{
## Generate a bivariate time series
set.seed(100)
X <- matrix(rnorm(2 * 2^8), ncol = 2)
X[1:2^7, 2] <- 3 * (X[1:2^7, 2] + 0.95 * X[1:2^7, 1])
X[-(1:2^7), 2] <- X[-(1:2^7), 2] - 0.95 * X[-(1:2^7), 1]
X[-(1:2^7), 1] <- X[-(1:2^7), 1] * 4
X <- as.ts(X)

## Haar wavelet, apply same smoothing to all levels & optimize
EWS <- mvEWS(X, kernel.name = "daniell", kernel.param = 20,
  optimize = TRUE)
summary(EWS)
print(EWS)
plot(EWS, style = 2, info = 1)
}

\concept{summary.mvLSW}
\concept{print.mvLSW}
