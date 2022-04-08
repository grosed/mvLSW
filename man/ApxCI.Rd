% File src/library/base/man/ApxCI.Rd
\name{ApxCI}
\alias{ApxCI}
\title{Evaluate the Approximate Confidence Interval of a mvEWS Estimate}
\description{
  Evaluate the approximate confidence interval of a multivariate 
  evolutionary wavelet spectrum.
}
\usage{
  ApxCI(object, var = NULL, alpha = 0.05, ...)
}

\arguments{
  \item{object}{A \code{mvLSW} object containing the multivariate 
    evolutionary wavelet spectrum estimate.}
  \item{var}{A \code{mvLSW} object containing the variance estimate of 
    the wavelet spectrum. If this is \code{NULL} (default) then the 
    variance is estimates by calling the \code{varEWS} and using 
    \code{object}.}
  \item{alpha}{Type I error, a single numerical value within (0,0.5].}
  \item{...}{Additional arguments to be passed to the 
    \code{\link{varEWS}} command.}
}

\details{
  The command evaluates the approximate Gaussian confidence 
  intervals for the elements of the mvEWS estimate.
}

\value{
  Invisibly returns a list containing two \code{mvLSW} classed
  objects with names "L" and "U" that respectively identify the 
  lower and upper interval estimates.
}

\references{
  Taylor, S.A.C., Park, T.A. and Eckley, I. (2019) Multivariate 
  locally stationary wavelet analysis with the mvLSW R package.
  \emph{Journal of statistical software} \strong{90}(11) pp. 1--16,
  doi: 10.18637/jss.v090.i11. 

  Park, T. (2014) Wavelet Methods for Multivariate Nonstationary 
  Time Series, PhD thesis, Lancaster University, pp. 91-111.
}

\seealso{
  \code{\link{mvEWS}}, \code{\link{as.mvLSW}}, \code{\link{varEWS}}.
}

\examples{
## Define evolutionary wavelet spectrum, structure only on level 2
Spec <- array(0, dim = c(3, 3, 8, 256))
Spec[1, 1, 2, ] <- 10
Spec[2, 2, 2, ] <- c(rep(5, 64), rep(0.6, 64), rep(5, 128))
Spec[3, 3, 2, ] <- c(rep(2, 128), rep(8, 128))
Spec[2, 1, 2, ] <- Spec[1, 2, 2, ] <- punif(1:256, 65, 192)
Spec[3, 1, 2, ] <- Spec[1, 3, 2, ] <- c(rep(-1, 128), rep(5, 128))
Spec[3, 2, 2, ] <- Spec[2, 3, 2, ] <- -0.5
EWS <- as.mvLSW(x = Spec, filter.number = 1, family = "DaubExPhase", 
  min.eig.val = NA)

## Sample time series and estimate the EWS.
set.seed(10)
X <- rmvLSW(Spectrum = EWS)
EWS_X <- mvEWS(X, kernel.name = "daniell", kernel.param = 20)

## Evaluate asymptotic spectral variance 
SpecVar <- varEWS(EWS_X)

## Plot Estimate & 95% confidence interval
CI <- ApxCI(object = EWS_X, var = SpecVar, alpha = 0.05)
plot(x = EWS_X, style = 2, info = 2, Interval = CI)

}

\keyword{ApxCI}
