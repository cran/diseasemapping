\name{getBreaks}
\alias{getBreaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Age Breaks}
\description{
An internal function to return a list contains age breaks, ages in the population data set, 
sex in the population data set, and age sex groups will be used in the \code{formatPopulation} function.
}
\usage{
getBreaks(colNames, breaks = NULL)
}

\arguments{
  \item{colNames}{names from the population data set}
  \item{breaks}{the age breaks, i.e breaks =seq(0, 80, by= 10) }
}
 
\examples{
data(kentucky)
ageBreaks = getBreaks(names(kentucky), breaks=c(seq(0, 80, by=10), Inf))
ageBreaks
}
 
