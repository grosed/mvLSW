% File src/library/base/man/varEWS.Rd
\name{varEWS}
\alias{varEWS}
\title{Asymptotic Variance of the mvEWS Estimate}
\description{
  Calculates the asymptotic variance of a multivariate 
  evolutionary wavelet spectrum estimate.
}
\usage{
  varEWS(object, ACWIP = NULL, verbose = FALSE)
}

\arguments{
  \item{object}{A \code{mvLSW} object containing the multivariate 
    evolutionary wavelet spectrum. Matrices must be positive definite, 
	i.e. information item \code{min.eig.val} must be greater than zero.}
  \item{ACWIP}{4D array containing the wavelet autocorrelation
    inner product functions. Set to \code{NULL} by default and 
	therefore evaluated within the command based on the 
	information supplied by \code{object}.}
  \item{verbose}{Logical. Controls the printing of messages whist 
    the computations progress. Set as \code{FALSE} as default.}
}

\details{
  The \code{varEWS} commands evaluate the asymptotic variance of a 
  multivariate evolutionary wavelet spectrum (mvEWS) estimate. Note, 
  the variance is only applicable when the mvEWS is smoothed 
  consistently across all levels with list item \code{smooth.type="all"}.
  This can be written in terms of the smoothed 
  periodogram relating to the bias correction of the mvEWS estimate, 
  where \eqn{A_{j,k}} is the inner product matrix of the wavelet 
  autocorrelation function:
  
  \deqn{Var( \hat{S}^{(p,q)}_{j,k} ) 
    = \sum_{l_1,l_2=1}^{J} (A^{-1})_{j,l_1} (A^{-1})_{j,l_2} 
	Cov( \tilde{I}^{(p,q)}_{l_1,k}, \tilde{I}^{(p,q)}_{l_1,k})}
  
  The covariance between elements of the smoothed periodogram can also 
  be expressed in terms of the raw wavelet periodogram:
  
  \deqn{Cov( \tilde{I}^{(p,q)}_{l_1,k}, \tilde{I}^{(p,q)}_{l_1,k}) 
    = \sum_{m_1,m_2} W_{m_1} W_{m_2} Cov( I^{(p,q)}_{l_1,m_1}, 
	  I^{(p,q)}_{l_2,m_2} )}
  
  The weights \eqn{W_i}, for integer i, define the smoothing kernel function 
  that is evaluated by the \code{kernel} command. Note that \eqn{W_i = W_{-i}}
  and \eqn{\sum_i W_i = 1}.
    
  The final step is to derive the covariance of the raw periodogram. This has 
  a long derivation, which can be concisely calculated by:

  \deqn{Cov( I^{(p,q)}_{j,k}, I^{(p,q)}_{l,m} ) 
    = E(p,j,k,q,l,m)^2 + E(p,j,k,p,l,m)E(q,j,k,q,l,m)}
  where	
  \deqn{E(p,j,k,q,l,m) = \sum_{h=1}^{J} A^{k-m}_{j,l,h} S^{(p,q)}_h((k+m)/2T) }
  
  Here, \eqn{A^{\lambda}_{j,l,h}} defines the autocorrelation 
  wavelet inner product function and \eqn{S^{(p,q)}_{j}(k/T)} 
  is the true spectrum of the process between channels p & q, 
  level j and location k. The true spectrum is not always available 
  and so this may be substituted with the smoothed and bias corrected 
  mvEWS estimate. For practical purposes, if k+m is odd then the 
  average between the available spectrum values at neighbouring 
  locations are substituted.
  
  For efficiency purpose, if the \code{varEWS} command is going to 
  be called multiple times then it is highly recommended that the 
  autocorrelation wavelet inner product should be evaluated beforehand 
  by \code{AutoCorrIP} and supplied via the \code{ACWIP} argument.
}

\value{
  Invisibly returns a \code{mvLSW} object containing the asymptotic variance 
    of the multivariate evolutionary wavelet spectrum.
}

\references{
  Taylor, S.A.C., Park, T.A. and Eckley, I. (2019) Multivariate 
  locally stationary wavelet analysis with the mvLSW R package.
  \emph{Journal of statistical software} \strong{90}(11) pp. 1--16,
  doi: 10.18637/jss.v090.i11.

%  Park, T., Eckley, I. and Ombao, H.C. (2014) Estimating 
%  time-evolving partial coherence between signals via multivariate 
%  locally stationary wavelet processes. \emph{Signal Processing, 
%  IEEE Transactions on} \strong{62}(20) pp. 5240-5250.
  Park, T. (2014) Wavelet Methods for Multivariate Nonstationary 
  Time Series, PhD thesis, Lancaster University, pp. 91-111.
}

\seealso{
  \code{ipndacw}, \code{\link{AutoCorrIP}}, 
  \code{\link{as.mvLSW}}, \code{\link{mvEWS}}
}

\examples{
## Define evolutionary wavelet spectrum, structure only on level 2
Spec <- array(0, dim=c(3, 3, 8, 256))
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

\keyword{varEWS}
