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
\details{

}

\author{
Patrick Brown
  

}

\examples{
# creating SMR's
data(popdata)
data(casedata)

model = getRates(casedata, popdata, ~age*sex, breaks=seq(0, 90, by=10) )
ontario = getSMR(popdata,model,  casedata)
spplot(ontario, 'SMR')


# an example of a spatial model with glmmBUGS

\dontrun{
library(glmmBUGS)
# run the model


# merge results back in to popdata
}


# running the same thing with INLA



}
