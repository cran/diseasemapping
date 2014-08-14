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
  \item{model}{Model to standardize to, either \code{glm} model output or a vector of rates by age and sex group}
  \item{referencePopulation}{population to standardize to } 
  \item{scale}{compute the expected rate per `scale' number of people. }
}

\author{Lutong Zhou }

\examples{
data(kentucky)

kentucky2 = getSMR(kentucky, larynxRates, larynx, 
		regionCode="County")

data(referencepop)
newpop <- getStdRate(kentucky2$SMR, larynxRates, referencepop, scale=100000)

newpop[1:10]
}
 
