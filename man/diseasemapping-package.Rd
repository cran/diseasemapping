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
data(kentucky)

\dontrun{
larynxRates= cancerRates("USA", year=1998:2002,site="Larynx")
}

kentucky2 = getSMR(kentucky, larynxRates, larynx, 
		regionCode="County")

\dontrun{
library(spdep)
spplot(kentucky, "SMR")
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
\dontrun{
library(INLA)
# get rid of regions with no neighbours
ontario2 = ontario[! as.character(ontario$CSDUID) %in% 
c("3510005" ,"3501007", "3537001", "3551031", "3560065", "3560062"),]
popDataAdjMat2 = poly2nb(ontario2,
	row.names=as.character(ontario2[["CSDUID"]]))
nb2INLA("nb.graph",popDataAdjMat2)

ontario2$CSDUID = as.character(ontario2$CSDUID)
 
prior.iid=prior.besag=c(.42,.00015)
formula.bym = observed ~  f(CSDUID, 
    model = "bym", graph = "nb.graph", values=CSDUID ,
    param = c(prior.iid, prior.besag))
 
result1.bym = inla(formula.bym,family="poisson",data=ontario2@data,
	offset=ontario2@data$logExpected,
    verbose=FALSE, keep = TRUE, 
    control.predictor=list(compute=TRUE))
    

tomerge = result1.bym$summary.random$CSDUID
rownames(tomerge) = tomerge$ID
 
ontario2@data =   cbind(ontario2@data, 
	tomerge[as.character(ontario2$CSDUID),] )   
	
spplot(ontario2, "mean")


bivariate = inla.hyperpar(result1.bym)$joint.hyper
colnames(bivariate) = c("iid","spatial","logDens")
bivariate = SpatialPixelsDataFrame(bivariate[,1:2], bivariate)

	
}

}
