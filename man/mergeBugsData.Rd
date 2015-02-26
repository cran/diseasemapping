\name{mergeBugsData-methods}
\docType{methods}
\alias{mergeBugsData}
\alias{mergeBugsData-methods}
\alias{mergeBugsData,data.frame-method}
\alias{mergeBugsData,SpatialPolygonsDataFrame-method}

\title{Merge results from BUGS into a data.frame or SPDF}
\description{
 merge the result from bugs function
}
\usage{
\S4method{mergeBugsData}{data.frame}(
x, bugsSummary, by.x = NULL, newcol = "mean", ...
)
\S4method{mergeBugsData}{SpatialPolygonsDataFrame}(
x, bugsSummary, by.x = NULL, newcol = "mean", ...
)
}

\arguments{
  \item{x}{spatial polygon object i.e population data set (popdata)}
  \item{bugsSummary}{ posterior distribution result from summaryChain function }
  \item{by.x}{ the common term from the spatial polygon object and the bugs function result  }
  \item{newcol}{the summary statistic that to be merged back to the data frame  }
  \item{\dots}{ additional arguments }
}

\author{ Patrick Brown }


\examples{
#data(popdata)
#newdata = c("3560102"=2, "3560104"=3)
#popdatatry = mergeBugsData(popdata, newdata, by.x="CSDUID")

# if the data set is a spatial polygons data frame: 
#popdatatry = mergeBugsData.SpatialPolygonsDataFrame(popdata, newdata, by.x="CSDUID")

# if the data set is a data frame
#popdatatry = mergeBugsData.data.frame(popdata, newdata, by.x="CSDUID")

\dontrun{
library(glmmBUGS)
data(popdata)
data(casedata)

therates = getRates(casedata, popdata, ~age*sex)
ontario = getSMR(popdata, therates, casedata)
ontario@data = ontario@data[,c("CSDUID","observed","logExpected")]
library(spdep)
popDataAdjMat = poly2nb(ontario, ontario[["CDSUID"]])
library(glmmBUGS)
forBugs = glmmBUGS(formula=observed + logExpected ~ 1,
  effects="CSDUID", family="poisson", spatial=popDataAdjMat,
  data=ontario@data)

startingValues = forBugs$startingValues

source("getInits.R")
library(R2WinBUGS)
ontarioResult = bugs(forBugs$ragged, getInits, parameters.to.save = names(getInits()),
    model.file="model.bug", n.chain=3, n.iter=100, n.burnin=10, n.thin=2,
      program="winbugs", debug=TRUE)

data(ontarioResult)
ontarioParams = restoreParams(ontarioResult, forBugs$ragged)
ontarioSummary = summaryChain(ontarioParams)
ontario = mergeBugsData(ontario, ontarioSummary)
}

}
