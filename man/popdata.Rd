\name{popdata}
\alias{popdata}
\docType{data}
\title{Ontario 2006 population by census subdivision}
\description{
 Data set contains the information of population, by age, sex, and census subdivision.
}
\usage{data(popdata)}
\format{
A SpatialPolygonsDataFrame object, which needs the \code{sp} package for full functionality.  
}
\details{
This data is from the 2006 Census of canada offering by Statistics Canada web site, \code{www12.statcan.gc.ca/english/census06/data/highlights/agesex/Index_PR.cfm?Lang=E&Geo=CSD&Table=1} 
}


\examples{
data(popdata)
head(popdata@data)
## maybe str(popdata) ; plot(popdata) ...
}
\keyword{datasets}
