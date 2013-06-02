\name{referencepop}
\alias{referencepop}
\docType{data}
\title{ Standard Canadian population data set from year 1991.}
\description{
  A data set contains population and age sex group from year 1991.
}
\usage{data(referencepop)}
\format{
Data frame with columns POPULATION, sex, and age  for the Canada 1991 population.
}
\details{
  data frame with rows representing age-sex groups, first column giving proportion of 
  Canada 1991 population in that group, and subsequent columns giving sex  and start of
  age range for each group
}


\examples{
data(referencepop)
head(referencepop)
sum(referencepop$POPULATION)
}
\keyword{datasets}
