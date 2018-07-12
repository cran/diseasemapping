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
}

\arguments{

\item{formula}{model formula, defaults to intercept-only model suitable for
 output from  \code{\link[diseasemapping]{getSMR}} if \code{data} is a
 \code{SpatialPolygonsDataFrame}. }
   \item{data}{
The observations and covariates for the model, can be output from 
 \code{\link[diseasemapping]{getSMR}}.
}
\item{adjMat}{An object of class \code{nb} containing the adjacency matrix.  If not supplied
it will be computed from \code{data}, but is required if \code{data} is a \code{SpatialPolygonDataFrame}}
\item{region.id}{Variable in \code{data} giving identifiers for the spatial regions.  
If not supplied, row numbers will be used.}
\item{priorCI}{named list of vectors specifying priors, see Details}
\item{family}{distribution of the observations, defaults to \code{poisson}}
\item{formula.fitted}{formula to use to compute the fitted values, defaults to the model 
formula but may, for example, exclude individual-level covariates.}
  \item{\dots}{
Additional arguments passed to  
\Sexpr[results=rd]{c(
  '\\\\code{inla} in the \\\\code{INLA} package', 
  '\\\\code{\\\\link[INLA]{inla}}'
  )[1+requireNamespace('INLA', quietly=TRUE)]
  }, such as 
\Sexpr[results=rd]{c(
  '\\\\code{control.inla}', 
  '\\\\code{\\\\link[INLA]{control.inla}}'
  )[1+requireNamespace('INLA', quietly=TRUE)]}
}
}

\details{

The Besag, York and Mollie model for Poisson distributed case counts is:

\deqn{Y_i \sim Poisson(O_i \lambda_i)}{Y_i~Poisson(O_i \lambda_i)}
\deqn{\log(\mu_i) = X_i \beta + U_i}{log(\mu_i) = X_i \beta + U_i}
\deqn{U_i \sim BYM(\sigma_1^2 , \sigma_2^2)}{U_i ~ BYM(\sigma_1^2 , \sigma_2^2)}

\itemize{
\item{\eqn{Y_i} is the response variable for region \eqn{i}, on the left side of the \code{formula} argument.}
\item{\eqn{O_i} is the 'baseline' expected count, which is specified 
in \code{formula} on the log scale with \eqn{\log(O_i)}{log(O_i)} an \code{offset} variable.}
\item{\eqn{X_i} are covariates, on the right side of \code{formula}}
\item{\eqn{U_i} is a spatial random effect, with a spatially structured variance parameter \eqn{\sigma_1^2}
and a spatially independent variance \eqn{\sigma_2^2}.}
}

The \code{priorCI} argument can be 
a \code{list} containing elements named \code{sdSpatial} and \code{sdIndep}, each being a vector of length 2 with 
2.5pct and 97.5pct quantiles for the prior distributions of the standard deviations
\eqn{\sigma_1} and \eqn{\sigma_2} respectively. 
Gamma prior distributions for the precision parameters 
\eqn{1/\sigma_1^2} and \eqn{1/\sigma_2^2}
yielding quantiles specified
for the standard deviations are computed, and used with the \code{model="bym"} option to 
\Sexpr[results=rd]{c(
  '\\\\code{f}', 
  '\\\\code{\\\\link[INLA]{f}}'
  )[1+requireNamespace('INLA', quietly=TRUE)]
  }.  

The other possible format for \code{priorCI} is to have elements named \code{sd} and \code{propSpatial}, which 
specifies \code{model="bym2"} should be used with penalized complexity priors.  
The \code{sd} element gives a prior for the marginal standard deviation
\eqn{\sigma_0 =\sqrt{\sigma_1^2+\sigma_2^2}}{\sigma_0 = sqrt(\sigma_1^2+\sigma_2^2)}.  
This prior is approximately exponential, and \code{priorCI$sd = c(1, 0.01)} specifies a 
 prior probability \eqn{pr(\sigma_0 > 1) = 0.01}.  
 The \code{propSpatial} element gives the prior for the ratio
 \eqn{\phi = \sigma_1/\sigma_0}.  Having \code{priorCI$propSpatial = c(0.5, 0.9)} implies 
 \eqn{pr(\phi < 0.5) = 0.9}.
}
 
\value{
A list containing
\item{inla}{results from the call to 
\Sexpr[results=rd]{c(
  '\\\\code{inla}', 
  '\\\\code{\\\\link[INLA]{inla}}'
  )[1+requireNamespace('INLA', quietly=TRUE)]}.  
Two additional elements are added: \code{marginals.bym} for the marginal distributions of the spatial random
effects, and \code{marginals.fitted.bym} for the marginals of the fitted values.  }
\item{data}{A \code{data.frame} or \code{SpatialPolygonsDataFrame} containing posterior means
and quantiles of the spatial random effect and fitted values.}
\item{parameters}{Prior and posterior distributions of the two covariance parameters, 
and a table summary with posterior quantiles of all model parameters. 
}
}

\author{
Patrick Brown
}


\seealso{
\url{http://r-inla.org}, 
\Sexpr[results=rd]{c(
  '\\\\code{inla}', 
  '\\\\code{\\\\link[INLA]{inla}}'
  )[1+requireNamespace('INLA', quietly=TRUE)]
  }
\code{\link[geostatsp]{glgm}}, \code{\link[diseasemapping]{getSMR}} 
}

\examples{

data('kentucky')

# must have an internet connection to do the following
\dontrun{
	larynxRates= cancerRates("USA", year=1998:2002,site="Larynx")
	dput(larynxRates)
}

larynxRates = structure(c(0, 0, 0, 0, 0, 0, 1e-06, 6e-06, 2.3e-05, 4.5e-05, 
9.9e-05, 0.000163, 0.000243, 0.000299, 0.000343, 0.000308, 0.000291, 
0.000217, 0, 0, 0, 0, 0, 1e-06, 1e-06, 3e-06, 8e-06, 1.3e-05, 
2.3e-05, 3.5e-05, 5.8e-05, 6.8e-05, 7.5e-05, 5.5e-05, 4.1e-05, 
3e-05), .Names = c("M_0", "M_5", "M_10", "M_15", "M_20", "M_25", 
"M_30", "M_35", "M_40", "M_45", "M_50", "M_55", "M_60", "M_65", 
"M_70", "M_75", "M_80", "M_85", "F_0", "F_5", "F_10", "F_15", 
"F_20", "F_25", "F_30", "F_35", "F_40", "F_45", "F_50", "F_55", 
"F_60", "F_65", "F_70", "F_75", "F_80", "F_85"), 
site = "Larynx", area = "USA, SEER", year = "1998-2002")

# get rid of under 10's
larynxRates = larynxRates[-grep("_(0|5)$",names(larynxRates))]

kentucky = getSMR(kentucky, larynxRates, larynx, regionCode="County")

if( require("spdep", quietly=TRUE)) {

kBYM = bym(observed ~ offset(logExpected) + poverty, kentucky, 
	 priorCI = list(sdSpatial=c(0.1, 5), sdIndep=c(0.1, 5)),
	 control.mode=list(theta=c(3.52, 3.35),restart=TRUE))

	kBYM$par$summary
	
	if(requireNamespace('geostatsp', quietly=TRUE))
		kBYM$data$exc1 = geostatsp::excProb(
			kBYM$inla$marginals.fitted.bym, log(1.2)
			)
}  else {
	kBYM = list()
}




if(require('mapmisc', quietly=TRUE) & length(kBYM$data$fitted.exp)){

thecol = colourScale(kBYM$data$fitted.exp, 
	breaks=5, dec=1, opacity = 0.7)

map.new(kBYM$data)

\dontrun{
kmap = openmap(kBYM$data)
plot(kmap,add=TRUE)
}

plot(kBYM$data, col=thecol$plot,add=TRUE)
legendBreaks("topleft", thecol)

}

}

\keyword{ Bayesian inference }
\keyword{ disease mapping }

