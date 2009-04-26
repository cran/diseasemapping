\name{getBreaks}
\alias{getBreaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Age Breaks}
\description{
getBreaks will return a list contains age breaks, ages in the population data set, 
sex in the population data set, and age sex groups will be used in the formatPopulation function.
}
\usage{
getBreaks(colNames, breaks = NULL)
}

\arguments{
  \item{colNames}{names from the population data set}
  \item{breaks}{the age breaks, i.e breaks =seq(0, 80, by= 10) }
}
\details{

}

\references{  }
\author{}
\note{
}
\seealso{  }
\examples{
data(popdata)
ageBreaks = getBreaks(names(popdata), breaks=c(seq(0, 80, by=10), Inf))
ageBreaks
}
\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
