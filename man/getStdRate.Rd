\name{getStdRate}
\alias{getStdRate}
\title{Calculate the standardized rate }
\description{
 A function to calculate the standard rate according to the Canadian standard population data set from year 1991. 
}
\usage{
getStdRate(relativeRate, model, referencePopulation, scale = 1e+05)
}
\arguments{
  \item{relativeRate}{ the relative cancer rate calculated by glmmBUGS of different sex and age group of people from ontario . }
  \item{model}{ the estimated cancer rate calculated by glm model of different sex and age group of people from ontario. }
  \item{referencePopulation}{the standard Canadian population from year 1991. }
  \item{scale}{ the unit of the population. }
}

\author{Lutong Zhou }

\examples{
data(referencepop)
data(casedata)
data(popdata)
# get the age sex group rate:
model = getRates(casedata, popdata, ~age*sex, S=c("M","F"))
# calculate the relative rate using glmmBUGS:
ontario = getSMR(popdata, model, casedata,regionCode="CSDUID", regionCodeCases="CSD2006")
newpop <- getStdRate(ontario$SMR, model, referencepop, scale=100000)

}
 
