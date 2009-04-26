\name{formatCases}
\alias{formatCases}
\title{ Format the disease case data set}
\description{
 The formatCases funtion formats the case data set. 
 Changes other formats of age and sex group to three columns: age, ageNumeric and sex. 
}
\usage{
formatCases(casedata, ageBreaks = NULL, years = NULL, aggregate.by = NULL)
}
\arguments{
  \item{casedata}{disease cases data set, usually a data.frame contains age and sex and number of cases.}
  \item{ageBreaks}{results from getBreaks function. }
  \item{years}{if it contains multiple years, define which years will be included in. }
  \item{aggregate.by}{if want to view the data set from a macro way, could aggregate the data set by age or sex or other variables. }
}
\details{
After using formatCases function, the age columns will be a "character" column contains ages in a cut format, i.e [50,55), denotes age 50. 
The cut breaks will get from the breaks from population data set or define by user.  
The original age column will changed to ageNumeric columns as factors. 
The sex column will have two levels "M" and "F" as factors. 
If "aggregate.by" is not NULL, the number of cases will be sum up by the groups define in aggregate.by. 
}
\value{
formatCases function will return a data frame. 
}
\references{ }
\author{Patrick Brown  }
\note{ 

}
\seealso{  }
\examples{
data(casedata)
data(popdata)
head(casedata)
caseformat <- formatCases(casedata, ageBreaks = getBreaks(names(popdata@data)))
head(caseformat)
caseformatagg <- formatCases(casedata, ageBreaks = getBreaks(names(popdata@data)), aggregate.by=c("age", "sex"))
head(caseformatagg)

}

\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
