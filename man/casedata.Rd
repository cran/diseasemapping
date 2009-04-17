\name{casedata}
\alias{casedata}
\docType{data}
\title{Data set contains the number of cases information}
\description{
 Cases of Hepatitis Z in Ontario.
}
\usage{data(casedata)}
\format{
data frame 
}
\details{
This dataset refers to cases of Hepatitis Z in Ontario for the years 1916 to 1918, giving the number of cases in each census subdivision by age, sex and year.  For reasons of privacy, any counts between 1 and 5 have been changed to 1.
  
}

\examples{
data(casedata)
head(casedata)
table(casedata$cases)
tapply(casedata$cases, casedata$age, sum)

## maybe str(casedata) ; plot(casedata) ...
}
\keyword{datasets}
