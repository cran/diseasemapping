\name{expandCensus}
\Rdversion{1.1}
\alias{expandCensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Given population data at each census year and range of the years that the data spans, this function
returns a new list of population data for each consective year assuming population is constant between two census years
}
\description{
Use if census are not done at each year
}
\usage{
expandCensus(popdata, year.range = range(as.integer(names(popdata))), years = as.integer(names(popdata)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{popdata}{
A list of population data at each census year
}
  \item{year.range}{
The range of the years which the census covers
}
  \item{years}{
The years that the census data are collected
}
}
\details{

}
\value{

}
\references{
%% ~put references to the literature/web site here ~
}
\author{  Patrick Brown

}
\note{

}



\seealso{
}

\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{ }% __ONLY ONE__ keyword per line
