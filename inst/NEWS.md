# Changelog for mvLSW R package

## Version 1.1 (24/02/2017)   

- Support time series classes "zoo" and "xts" 
- Add S3 print() and simulate() commands
- Remove internal only commands from NAMESPACE

## Version 1.2 (04/08/2017)

- Correct error when tol=-Inf in mvEWS() to not apply regularisation
- Corrections and general edits to manual pages
- Convert Asymp_Quantile to ApxCI, make appropriate edits to plot.mvLSW

## Version 1.2.1 (18/08/2017)
- Allow tol=NA in mvEWS() to specify not to apply positiveness correction                 
- Correct error in RawPeriodogram when defining attributes of output
- Added Vignette to package

## Version 1.2.3 (05/08/2019) 
- Added JSS paper citation

## Version 1.2.4 (02/09/2011)
- Fixed some minor CRAN check issues
- Changed maintainership