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
  \item{relativeRate}{ the relateve cancer rate calculated by glmmBUGS of different sex and age group of people from ontario . }
  \item{model}{ the estimated cancer rate calculated by glm model of different sex and age group of people from ontario. }
  \item{referencePopulation}{the standard Canadian population from year 1991. }
  \item{scale}{ the unit of the population. }
}
\details{
  
}
\value{
  
}
\references{  }
\author{Lutong Zhou }
\note{ 
}
\seealso{ }
\examples{
data(referencepop)
data(casedata)
data(popdata)
# get the age sex group rate:
model = getRates(casedata, popdata, ~age*sex, S=c("M","F"))
# calculate the relative rate using glmmBUGS:
ontario = getSMR(popdata, model, casedata)
\dontrun{
library(spdep)
popDataAdjMat <- poly2nb(popdata, popdata[["CSDUID"]])
library(glmmBUGS)
forBugs = glmmBUGS(formula = observed + logExpected ~ 1, effects= "CSDUID",
                   family="poisson", spatial =popDataAdjMat, data=ontario@data)
startingValues=forBugs$startingValues
source("getInits.R")
library(R2WinBUGS)
ontarioResult = bugs(forBugs$ragged, getInits, parameters.to.save=names(getInits()),
                     model.file="model.bug", n.chain=3, n.iter=1000, n.burnin=100,
                     n.thin=10, program="winbugs")
ontarioParams=restoreParams(ontarioResult, forBugs$ragged)
ontarioSummary = summaryChain(ontarioParams)
ontario = mergeBugsData(ontario, ontarioSummary)
ontario@data = ontario@data[!is.na(ontario$SMR),]
relativeRate = ontario@data[,grep("Fitted", names(ontario), value=TRUE)]
newpop <- getStdRate(relativeRate, model, referencepop, scale=100000)
}
}
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
