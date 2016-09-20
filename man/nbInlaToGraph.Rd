\name{nbToInlaGraph}
\alias{nbToInlaGraph}


\title{
Write a graph file for INLA
}
\description{
Writes a graph file from an adjacency matrix suitable for use with INLA.
}

\usage{
nbToInlaGraph(adjMat, graphFile="graph.dat")
}

\arguments{
\item{adjMat}{An object of class \code{nb} containing the adjacency matrix.}
\item{graphFile}{name of file to save adjacency matrix to.}
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
\code{\link[spdep]{poly2nb}},
\code{\link[spdep]{nb2INLA}} 
}

\examples{

data('kentucky')
# remove all the neighbours Ballard county
kSub = kentucky[-c(2,20,79),]

if( require("spdep", quietly=TRUE)) {


adjMat = poly2nb(
	kSub, 
	row.names=kSub$County,
	queen=FALSE
)

nFile = tempfile()
nbRes = nbToInlaGraph(adjMat, nFile)

# Ballard is region 3
nbRes['Ballard']
# note that Ballard has no neighbours
adjMat[[3]]

cat(readLines(nFile, n=5), sep='\n')

\dontrun{
# there will be a warning about zero neighbours
junk = bym(poverty ~ 1, data=kSub, family='gaussian')
}

}

}

