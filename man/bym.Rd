\name{bym-methods}
\docType{methods}
\alias{bym}
\alias{bym-methods}
\alias{bym,formula,ANY,ANY,missing-method}
\alias{bym,formula,ANY,missing,missing-method}
\alias{bym,formula,SpatVector,NULL,character-method}
\alias{bym,formula,SpatVector,missing,character-method}
\alias{bym,formula,SpatVector,matrix,character-method}
\alias{bym,formula,data.frame,matrix,character-method}


\title{
Fit the BYM model
}
\description{
Uses inla to fit a Besag, York and Mollie disease mapping model
}
\usage{
\S4method{bym}{formula,ANY,ANY,missing}(formula,data,adjMat,region.id,...)
\S4method{bym}{formula,ANY,missing,missing}(formula,data,adjMat,region.id,...)
\S4method{bym}{formula,SpatVector,NULL,character}(formula, data, adjMat, region.id, ...)
\S4method{bym}{formula,SpatVector,missing,character}(formula, data, adjMat, region.id, ...)
\S4method{bym}{formula,SpatVector,matrix,character}(formula,data,adjMat,region.id,...)
\S4method{bym}{formula,data.frame,matrix,character}(
formula,data,adjMat,region.id,
prior=list(sd=c(0.1,0.5),propSpatial=c(0.5,0.5)),
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
\item{prior}{named list of vectors specifying priors, see Details}
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

The \code{prior} 
has elements named \code{sd} and \code{propSpatial}, which 
specifies \code{model="bym2"} should be used with penalized complexity priors.  
The \code{sd} element gives a prior for the marginal standard deviation
\eqn{\sigma_0 =\sqrt{\sigma_1^2+\sigma_2^2}}{\sigma_0 = sqrt(\sigma_1^2+\sigma_2^2)}.  
This prior is approximately exponential, and \code{prior$sd = c(1, 0.01)} specifies a 
 prior probability \eqn{pr(\sigma_0 > 1) = 0.01}.  
 The \code{propSpatial} element gives the prior for the ratio
 \eqn{\phi = \sigma_1/\sigma_0}.  Having \code{prior$propSpatial = c(0.5, 0.9)} implies 
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
\code{\link[diseasemapping]{getSMR}} 
}

\examples{

data('kentucky')
kentucky = terra::unwrap(kentucky)

# get rid of under 10s
larynxRates = larynxRates[-grep("_(0|5)$",names(larynxRates))]

kentucky = getSMR(kentucky, larynxRates, larynx, regionCode="County")

kBYM <- try(
  bym(
    observed ~ offset(logExpected) + poverty, kentucky,
	  prior= list(sd=c(0.1, 0.5), propSpatial=c(0.5, 0.5))
   ), silent=TRUE)

if(length(grep("parameters", names(kBYM)))) {
  kBYM$parameters$summary
}



if( require("mapmisc", quietly=TRUE) && length(grep("parameters", names(kBYM))) ) {
  thecol = colourScale(kBYM$data$fitted.exp, breaks=5, dec=1, style="equal")
  terra::plot(kBYM$data, col=thecol$plot)
  legendBreaks("topleft", thecol)
}


}


\keyword{ Bayesian inference }
\keyword{ disease mapping }

