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
data(popdata)
data(casedata)

model = getRates(casedata, popdata, ~age*sex, breaks=seq(0, 90, by=10) )
ontario = getSMR(popdata,model,  casedata,regionCode="CSDUID", regionCodeCases="CSD2006")
\dontrun{
spplot(ontario, 'SMR')
}


# an example of a spatial model with glmmBUGS

\dontrun{
# run the model
library(spdep)
popDataAdjMat = poly2nb(ontario,row.names=as.character(ontario[["CSDUID"]]))

library(glmmBUGS)
forBugs = glmmBUGS(formula=observed + logExpected ~ 1,
  effects="CSDUID", family="poisson", spatial=popDataAdjMat,
  data=ontario@data)
startingValues = forBugs$startingValues
source("getInits.R")
library(R2WinBUGS)
ontarioResult = bugs(forBugs$ragged, getInits, parameters.to.save = names(getInits()),
    model.file="model.bug", n.chain=3, n.iter=100, n.burnin=10, n.thin=2,
      program="winbugs", debug=TRUE)

data(ontarioResult)
ontarioParams = restoreParams(ontarioResult, forBugs$ragged)
ontarioSummary = summaryChain(ontarioParams)

# merge results back in to popdata
ontario = mergeBugsData(ontario, ontarioSummary)
}


# running the same thing with INLA



}
