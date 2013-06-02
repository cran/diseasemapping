\name{kentucky}
\alias{kentucky}
\alias{larynx}
\alias{larynxRates}
\docType{data}
\title{Larynx cancer cases and population in Kentucky}
\description{
 Data set contains the information of population, by age, sex, and census subdivision.
}
\usage{data(kentucky)}
\format{A SpatialPolygonsDataFrame of Kentucky county boundaries and populations, a data frame 
of case numbers for each county, and a vector of cancer rates by age and sex group.

}
\details{
The cancer case data are obtained from www.cancer-rates.info and are for a single  
deliberately unspecified year. 
}



\examples{
data(kentucky)

\dontrun{
larynxRates= cancerRates("USA", year=1998:2002,site="Larynx")
}

kentucky2 = getSMR(kentucky, larynxRates, larynx, 
		regionCode="County")

breaks = c(0,1,seq(2, ceiling(max(kentucky2$SMR)),by=2))
thecol = terrain.colors(length(breaks)-1)

plot(kentucky2, col = thecol[cut(kentucky2$SMR, breaks,include.lowest=TRUE)] )
legend("topleft", pch=15, pt.cex=2.5, adj=c(0,15), 
  legend=rev(breaks), col=c(NA, rev(thecol)))
}
\keyword{datasets}
