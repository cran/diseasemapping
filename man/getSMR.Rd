\name{getSMR-methods}
\docType{methods}
\alias{getSMR}
\alias{getSMR-methods}
\alias{getSMR,SpatialPolygonsDataFrame,ANY,ANY,ANY,ANY-method}
\alias{getSMR,list,ANY,ANY,ANY,ANY-method}
\alias{getSMR,data.frame,ANY,missing,missing,missing-method}

\alias{getSMR,data.frame,ANY,data.frame,missing,missing-method}
\alias{getSMR,data.frame,ANY,data.frame,character,missing-method}
\alias{getSMR,data.frame,ANY,missing,character,missing-method}

\alias{getSMR,data.frame,ANY,data.frame,character,character-method}


\title{Calculate the standardized mortality/morbidity ratios}
\description{
Calculates the rate of observe value over expected value. 
It will also merge back the observed value, expected value and the ratio back to the population data set. 
}

\usage{
\S4method{getSMR}{SpatialPolygonsDataFrame,ANY,ANY,ANY,ANY}(
popdata, model, casedata, regionCode , regionCodeCases , 
area.scale = 1,  sex=c('m','f'),...
)
\S4method{getSMR}{list,ANY,ANY,ANY,ANY}(
popdata, model, casedata, regionCode , regionCodeCases , 
area.scale=1, sex=c('m','f'), ...
)
\S4method{getSMR}{data.frame,ANY,missing,missing,missing}(
popdata, model, casedata, regionCode , regionCodeCases , 
area.scale = 1,  sex=c('m','f'),...
)

\S4method{getSMR}{data.frame,ANY,data.frame,missing,missing}(
popdata, model, casedata, regionCode , regionCodeCases , 
area.scale = 1,  sex=c('m','f'),...
)
\S4method{getSMR}{data.frame,ANY,data.frame,character,missing}(
popdata, model, casedata, regionCode , regionCodeCases , 
area.scale = 1,  sex=c('m','f'),...
)
\S4method{getSMR}{data.frame,ANY,missing,character,missing}(
popdata, model, casedata, regionCode , regionCodeCases , 
area.scale = 1,  sex=c('m','f'),...
)

\S4method{getSMR}{data.frame,ANY,data.frame,character,character}(
popdata, model, casedata, regionCode , regionCodeCases , 
area.scale = 1,  sex=c('m','f'),...
)

}
\arguments{
  \item{popdata}{ the name of population data set  }
  \item{model}{rates, either fitted model (usually a \code{glm} object), or a vector of rates. }
  \item{casedata}{ the name of case data set }
  \item{regionCode}{ the name of district area column in population data set}
  \item{regionCodeCases}{the name of district area column in case data set }
  \item{area.scale}{ scale the unit of area. e.g $10^6$: if your spatial coordinates are metres 
  and you want intensity in cases per km2 }
  \item{sex}{possible subsetting of cases and population, set \code{sex='f'} for females only.}
 \item{\dots}{ additional arguments.   When \code{popdata} is a \code{list}, arguments 
 can be 
  \code{personYears} (logical, convert rates to person years), 
  \code{years} (a vector with the year of each dataset), or
  \code{year.range} (two dimensional vector with first and last year)
 }
}
\details{
If \code{model} is numeric, it's assumed to be a vector of rates, with the names of the elements corresponding to columns of the population data set.  Names do not need to match exactly (can have M in one set of names, male in another for instance).

Otherwise, \code{model} is passed to the \code{predict} function.  
}
\value{
  Returns a new population data set contains expected number of cases, observed number of cases and SMR. 
  It has the same format as the population data set which put into the function.
}
 
\examples{
data(kentucky)

kentucky2 = getSMR(kentucky, larynxRates, larynx, 
		regionCode="County")

data.frame(kentucky2)[1:10,grep("^F|^M", names(kentucky2), invert=TRUE)]

}
 