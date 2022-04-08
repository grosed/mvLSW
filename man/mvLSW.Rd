% File src/library/base/man/mvLSW.Rd
\name{mvLSW}
\alias{mvLSW}
\title{Multivariate, Locally Stationary Wavelet Process Estimation}
\description{
  The mvLSW package provides an implementation of the multivariate locally 
  stationary time series modelling approach proposed by Park, Eckley and 
  Ombao (2014).
  
  The approach extends the locally stationary wavelet time series work of 
  Nason, von Sachs and Kroisandt (2000) to a multivariate setting, introducing 
  wavelet-based measures of local coherence and local partial coherence. The 
  package implements the estimation scheme by Park et al. (2014) for such 
  processes. Note that mvLSW should be used in conjunction with the 
  \code{wavethresh} package developed by Nason (2016).
}

\details{
  Package: mvLSW
  
  Type: Package
  
  Version: 1.2.3
  
  Date: 2019-08-05
  
  License: GPL(>=3)
}

\author{
  Simon Taylor, <s.taylor2@lancaster.ac.uk>
}

\references{
  Taylor, S.A.C., Park, T.A. and Eckley, I. (2019) Multivariate 
  locally stationary wavelet analysis with the mvLSW R package.
  \emph{Journal of statistical software} \strong{90}(11) pp. 1--16,
  doi: 10.18637/jss.v090.i11. 

  Park, T.A., Eckley, I. and Ombao, H.C. (2014) Estimating 
  time-evolving partial coherence between signals via multivariate 
  locally stationary wavelet processes \emph{IEEE Transactions on 
  Signal Processing} \strong{62}(20), pp. 5240--5250.
  
  Nason, G.P., von Sachs, R. and Kroisandt, G. (2000) Wavelet 
  processes and adaptive estimation of the evolutionary wavelet 
  spectrum \emph{Journal of the Royal Statistical Society B} 
  \strong{62}, pp. 271--292.

  Nason, G. (2016) wavethresh: Wavelets Statistics and Transforms. 
  R package version 4.6.8. 
  
  https://CRAN.R-project.org/package=wavethresh
}

\seealso{
  \code{mvEWS}, \code{coherence}, \code{rmvLSW},
  \code{\link{as.mvLSW}}, \code{\link{mvEWS}}
}

\examples{
#
# See examples in individual help pages
#
}

\keyword{varEWS}
