\name{getSMR}
\alias{getSMR}
\alias{getSMR.data.frame}
\alias{getSMR.list}
%\alias{getSMR.PostgreSQLConnection}
\alias{getSMR.SpatialPolygonsDataFrame}
\title{Calculate the standardized mortality/morbility ratios   }
\description{
The getSMR function calculates the rate of observe value over expected value. 
It will also merge back the observed value, expected value and the ratio back to the population data set. 
}
\usage{
getSMR(popdata, model, casedata, regionCode = "CSDUID", regionCodeCases = "CSD2006", area = FALSE, area.scale = 1,  ...)
\method{getSMR}{data.frame}(popdata, model, casedata, regionCode = "CSDUID", regionCodeCases = "CSD2006", area = FALSE, area.scale = 1, formatPop=TRUE,...)
\method{getSMR}{list}(popdata, model, casedata = NULL, regionCode = "CSDUID", regionCodeCases = "CSD2006", area = FALSE, area.scale = 1, years = NULL, year.range = NULL, ...)
\method{getSMR}{SpatialPolygonsDataFrame}(popdata, ...)
}
\arguments{
  \item{popdata}{ the name of population data set  }
  \item{model}{ the estimated coefficients we calculated from getRates function. }
  \item{\dots}{additional arguments}
  \item{casedata}{ the name of case data set }
  \item{regionCode}{ the name of district area column in population data set}
  \item{regionCodeCases}{the name of district area column in case data set }
  \item{area}{ if \code{TRUE}, calculate the expected number of cases per area in each region }
  \item{area.scale}{  control the unit of area. e.g $10^6$: per square kilometers }
  \item{formatPop}{if \code{TRUE}, then format the population }
  \item{years}{the vector of years}
  \item{year.range}{the range of the year vector} 
  
}
\details{
  
 If there is no rate calculates from the getRates function, 
 it will leave NA in the merged population data set. 
}
\value{
  Returns a new population data set contains expected number of cases, observed number of cases and SMR. 
  It has the same format as the population data set which put into the function.
}
\references{  }
\author{ }
\note{ 
}
\seealso{ }
\examples{
data(casedata)
data(popdata)
therates = getRates(casedata, popdata, ~age*sex,breaks=c(seq(0, 80, by=10), Inf) )
thesmr = getSMR(popdata, therates, casedata)
head(thesmr@data)
}
\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
