\name{nbToInlaGraph}
\alias{nbToInlaGraph}


\title{
Write a graph file for INLA
}
\description{
Writes a graph file from an adjacency matrix suitable for use with INLA.
}

\usage{
nbToInlaGraph(adjMat, graphFile="graph.dat", region.id = attributes(adjMat)$region.id)
}

\arguments{
\item{adjMat}{An object of class \code{nb} containing the adjacency matrix.}
\item{graphFile}{name of file to save adjacency matrix to.}
\item{region.id}{names of regions}
}

\details{
	This function correctly handles regions which have zero neighbours.
}

\value{
	A vector of names and indices
}

\author{
Patrick Brown
}


\seealso{
\code{\link[terra]{adjacent}}
}

\examples{

data('kentucky')
kentucky = terra::unwrap(kentucky)
# remove all the neighbours Ballard county
kSub = kentucky[-c(2,20,79),]


adjMat = terra::adjacent(kSub)
attributes(adjMat)$region.id = kSub$County

nFile = tempfile()
nbRes = nbToInlaGraph(adjMat, nFile)

# Ballard is region 3
nbRes['Ballard']
# note that Ballard has no neighbours
table(adjMat[,'from']==3)

cat(readLines(nFile, n=5), sep='\n')

\donttest{
# there will be a warning about zero neighbours
junk = bym(poverty ~ 1, data=kSub, family='gaussian')
}


}

