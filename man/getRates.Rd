\name{getRates}
\alias{getRates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the estimated coefficients of age and sex group from the glm model }
\description{
  The getRates function calculates the estimated coefficient of the age and sex group from the case and population data set. It fits 
a glm model with poisson distribution by default.
}
\usage{
getRates(casedata, popdata, formula, family = poisson, minimumAge = 0, 
    maximumAge = 100, S = c("M", "F"), years = NULL, year.range = NULL, 
    case.years = grep("^year$", names(casedata), ignore.case = TRUE, 
        value = TRUE)[1], breaks = NULL) 

}

\arguments{
  \item{casedata}{ case data set name }
  \item{popdata}{population data set name }
  \item{formula}{ the glm model you want to fit. ie. $~age*sex$~ }
  \item{family}{the distribution to fit the model}
  \item{minimumAge}{ the lower boundry of the age, default is 0 }
  \item{maximumAge}{the higher boundry of the age, default is 100  }
  \item{S}{the group of sex you want to choose, default is both male and female. 
            "M" denotes male; "F" denotes female }
  \item{years}{  a vector of census years  }
  \item{year.range}{study period: a vector of two elements, starting dates and ending dates  }
  \item{case.years}{ variable name in the case data which contains time  }
  \item{breaks}{ the age breaks }
}
\details{
 It fits a glm model with poisson or binomial distribution over case and population data sets. 
 If there is no data set in some age and sex group, an NA will show there.
}
\value{
 A summary of the glm model contains set of estimated coefficients for different age and sex groups.
}
\references{  }
\author{ Patrick Brown  }
\note{ 
}
\seealso{ }
\examples{
data(casedata)
data(popdata)
therates = getRates(casedata, popdata, ~age*sex,breaks=c(seq(0, 80, by=10), Inf))
therates
}

\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
