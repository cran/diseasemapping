\name{kentucky}
\alias{kentucky}
\alias{kentuckyTract}
\alias{larynx}
\alias{larynxRates}
\docType{data}
\title{Larynx cancer cases and population in Kentucky}
\description{
 Data set contains the information of population, by age, sex, and census subdivision.
}
\usage{data('kentucky')}
\format{A \code{SpatialPolygonsDataFrame} of Kentucky boundaries and populations, 
 case numbers for each county, and a vector of cancer rates by age and sex group.

}
\details{

\code{larynx} is a \code{data.frame} of cancer case counts by county, 
obtained from \url{http://www.cancer-rates.info} and are for a single  
deliberately unspecified year. 

\code{kentucky} contains country boundaries and populations.

\code{kentuckyTract} contains census tract boundaries and populations.

}



\examples{
library('terra')
data('kentucky')
kentucky = terra::unwrap(kentucky)

head(larynx)
10^5*larynxRates[paste(c("M","F"), 50, sep="_")]



kentucky2 = getSMR(kentucky, larynxRates, larynx, 
		regionCode="County")

names(kentucky2)
length(kentucky2)

data('kentuckyTract')
kentuckyTract = unwrap(kentuckyTract)

length(kentuckyTract)


if(require('mapmisc', quietly=TRUE)) {
	mycol = colourScale(kentucky2$SMR, 
		breaks=10, dec=-log10(0.5), style='quantile')
	plot(kentucky2, col=mycol$plot, border='#00000040')
	legendBreaks('topright', mycol)
} else {
	terra::plot(kentucky2)
}

breaks = c(0,1,seq(2, ceiling(max(kentucky2$SMR,na.rm=TRUE)),by=2))
thecol = terrain.colors(length(breaks)-1)

plot(kentucky2, col = thecol[cut(kentucky2$SMR, 
	breaks,include.lowest=TRUE)] )

legend("topleft", pch=15, pt.cex=2.5, adj=c(0,15), 
  legend=rev(breaks), col=c(NA, rev(thecol)))

}

\keyword{datasets}
