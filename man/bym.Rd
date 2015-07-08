\name{bym-methods}
\docType{methods}
\alias{bym}
\alias{bym-methods}
\alias{bym,formula,ANY,ANY,missing-method}
\alias{bym,formula,ANY,missing,missing-method}
\alias{bym,formula,SpatialPolygonsDataFrame,NULL,character-method}
\alias{bym,formula,SpatialPolygonsDataFrame,missing,character-method}
\alias{bym,formula,SpatialPolygonsDataFrame,nb,character-method}
\alias{bym,formula,data.frame,nb,character-method}
\alias{nbToInlaGraph}


\title{
Fit the BYM model
}
\description{
Uses inla to fit a Besag, York and Mollie disease mapping model
}
\usage{
\S4method{bym}{formula,ANY,ANY,missing}(formula,data,adjMat,region.id,...)
\S4method{bym}{formula,ANY,missing,missing}(formula,data,adjMat,region.id,...)
\S4method{bym}{formula,SpatialPolygonsDataFrame,NULL,character}(formula, data, adjMat, region.id, ...)
\S4method{bym}{formula,SpatialPolygonsDataFrame,missing,character}(formula, data, adjMat, region.id, ...)
\S4method{bym}{formula,SpatialPolygonsDataFrame,nb,character}(formula,data,adjMat,region.id,...)
\S4method{bym}{formula,data.frame,nb,character}(
formula,data,adjMat,region.id,
priorCI=list(sdSpatial=c(0.01,2),sdIndep=c(0.01,2)),
family="poisson",formula.fitted=formula,...)
nbToInlaGraph(adjMat, graphFile="graph.dat")
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{formula}{model formula, defaults to intercept-only model suitable for
 output from  \code{\link[diseasemapping]{getSMR}} if \code{data} is a
 SpatialPointsDataFrame. }
   \item{data}{
The observations and covariates for the model, can be output from 
 \code{\link[diseasemapping]{getSMR}}.
}
\item{adjMat}{An object of class \code{nb} containing the adjacency matrix.  If not supplied
it will be computed from \code{data}, but is required if \code{data} is a SpatialPolygonDataFrame}
\item{graphFile}{name of file to save adjacency matrix to.}
\item{region.id}{Variable in \code{data} giving identifiers for the spatial regions.  
If not supplied, row numbers will be used.}
\item{priorCI}{prior 95 pct intervals for the standard deviation of the spatial and
independent random effects.
}
\item{family}{distribution of the observations, defaults to \code{poisson}}
\item{formula.fitted}{formula to use to compute the fitted values, defaults to the model 
formula but may, for example, exclude individual-level covariates.}
  \item{\dots}{
Additional arguments passed to  \code{inla} 
}
}
 
\value{
A list containing
\item{inla}{results from the call to \code{inla}.  Two additional elements are added: marginals.bym for the marginal distributions of the spatial random
effects, and marginals.fitted.bym for the marginals of the fitted values.  }
\item{data}{A data frame or SpatialPolygonsDataFrame containing posterior means
and quantiles of the spatial random effect and fitted values.}
\item{parameters}{summaries of posterior distributions, see \code{\link[geostatsp]{glgm}}}
}

\author{
Patrick Brown
}


\seealso{
\url{http://r-inla.org}, \code{\link[geostatsp]{glgm}}, \code{\link[diseasemapping]{getSMR}} 
}
\examples{

if( require("spdep", quietly=TRUE)) {
data('kentucky')

if(FALSE) {
	# must have an internet connection to do the following
	larynxRates= cancerRates("USA", year=1998:2002,site="Larynx")
	# get rid of under 10's
	larynxRates = larynxRates[-grep("_(0|5)$",names(larynxRates))]
	dput(larynxRates)
} else {
larynxRates = structure(c(0, 0, 0, 0, 1e-06, 6e-06, 2.3e-05, 4.5e-05, 9.9e-05, 
0.000163, 0.000243, 0.000299, 0.000343, 0.000308, 0.000291, 0.000217, 
0,0, 0, 0, 1e-06, 1e-06, 3e-06, 8e-06, 1.3e-05, 2.3e-05, 3.5e-05, 
5.8e-05, 6.8e-05, 7.5e-05, 5.5e-05, 4.1e-05, 3e-05), .Names = c("M_10", 
"M_15", "M_20", "M_25", "M_30", "M_35", "M_40", "M_45", "M_50", 
"M_55", "M_60", "M_65", "M_70", "M_75", "M_80", "M_85", "F_0", "F_10", 
"F_15", "F_20", "F_25", "F_30", "F_35", "F_40", "F_45", "F_50", 
"F_55", "F_60", "F_65", "F_70", "F_75", "F_80", "F_85"))

}

kentucky = getSMR(kentucky, larynxRates, larynx,
 regionCode="County")

kBYM = bym(observed ~ offset(logExpected) + poverty, kentucky, 
	 priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)),
	 control.mode=list(theta=c(3.52, 3.35),restart=TRUE))

	kBYM$par$summary
	
	if(require('geostatsp', quietly=TRUE))
		kBYM$data$exc1 = excProb(kBYM$inla$marginals.fitted.bym, log(1.2))
}


\dontrun{

library('mapmisc')

kmap = openmap(kBYM$data,)
map.new(kBYM$data)
plot(kmap,add=TRUE)
thecol = colourScale(kBYM$data$fitted.exp, breaks=5, dec=1,
opacity = 0.7)
plot(kBYM$data, col=thecol$plot,add=TRUE)
legendBreaks("topleft", thecol)


}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Bayesian inference }
\keyword{ disease mapping }% __ONLY ONE__ keyword per line
