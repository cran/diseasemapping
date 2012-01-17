\name{getRates}
\alias{getRates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the estimated coefficients of age and sex group from the glm model }
\description{
  The getRates function calculates the estimated coefficient of the age and sex group from the case and population data set. It fits 
a glm model with Poisson distribution by default.
}
\usage{
getRates(casedata, popdata, formula, family = poisson, minimumAge = 0, 
    maximumAge = 100, S = c("M", "F"), years = NULL, year.range = NULL, 
    case.years = grep("^year$", names(casedata), ignore.case = TRUE, 
        value = TRUE), fit.numeric=NULL, breaks = NULL) 
}
\arguments{
  \item{casedata}{A data frame of case data, with columns corresponding to variables in \code{formula}.  
  Assumed to be one row per case, unless a column called \code{y} or \code{cases} or \code{count} is
  included, in which case this column gives the number of cases per row. }
  \item{popdata}{population data set}
  \item{formula}{ the glm model you want to fit. ie. \code{~age*sex} }
  \item{family}{the distribution to fit the model}
  \item{minimumAge}{ the lower boundary of the age, default is 0 }
  \item{maximumAge}{the higher boundary of the age, default is 100  }
  \item{S}{vector of sexes to include in the analysis.  Defaults to both "M" and "F" }
  \item{years}{  a vector of census years  }
  \item{year.range}{study period: a vector of two elements, starting dates and ending dates  }
  \item{case.years}{ variable name in the case data which contains time  }
  \item{fit.numeric}{ the variables which needed to be changed from factor to numeric}
  \item{breaks}{ the age breaks }
}
\details{
 It fits a glm model with Poisson or binomial distribution over case and population data sets. 
 If there is no data set in some age and sex group, an NA will show there.
}
\value{
 A summary of the glm model contains set of estimated coefficients for different age and sex groups.
}


\author{ Patrick Brown  }


\examples{
data(casedata)
data(popdata)
therates = getRates(casedata, popdata, ~age*sex,breaks=c(seq(0, 80, by=10), Inf))
therates
}

 
