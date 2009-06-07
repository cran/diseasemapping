\name{computeArea}
\alias{computeArea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute areas of polygons}
\description{
  Computes the areas of all the polygons of a spatialPolygonsDataFrame object
}
\usage{
computeArea(sp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp}{A spatialPolygonsDataFrame}
}

\value{
  A list, each element containing the area of a polygon.
}
\seealso{}
\examples{
data(popdata)
theAreas = computeArea(popdata)
# histogram of areas of ontario CDS's
hist(unlist(theAreas))



}