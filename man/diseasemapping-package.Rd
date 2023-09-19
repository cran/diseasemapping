\name{diseasemapping-package}
\alias{diseasemapping-package}
\alias{diseasemapping}
\docType{package}
\title{
Disease Mapping
}
\description{
Functions for calculating observed and expected counts by region, and manipulating posterior 
samples from Bayesian models produced by glmmBUGS.
}


\author{
Patrick Brown
  

}

\examples{
# creating SMR's
data('kentucky')
kentucky = terra::unwrap(kentucky)


kentucky2 = getSMR(kentucky, larynxRates, larynx, 
		regionCode="County")

if(require('mapmisc', quietly=TRUE)) {
	mycol = colourScale(kentucky2$SMR, breaks=9, 
	dec=-log10(0.5), style='equal', transform='sqrt')
	plot(kentucky2, col=mycol$plot)
	legendBreaks('topleft', mycol)
}

if(require("INLA", quietly=TRUE)) {


kBYM = bym(observed ~ offset(logExpected) + poverty, 
	data= kentucky2, 
	prior = list(sd=c(0.1, 0.5), propSpatial=c(0.5, 0.5)),
	num.threads=2, blas.num.threads=2
)

kBYM$par$summary

}

}


