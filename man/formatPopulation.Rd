\name{formatPopulation}
\alias{formatPopulation}
\alias{formatPopulation.data.frame}
\alias{formatPopulation.list}
%\alias{formatPopulation.PostgreSQLConnection}
\alias{formatPopulation.SpatialPolygonsDataFrame}
\title{Format the population data set }
\description{
 The formatCases funtion formats the population data set. 
Reshape the population data set to "long" format, 
add in 4 columns : GROUP, POPULATION, sex and age. 
}
\usage{
formatPopulation(popdata, aggregate.by = NULL, breaks = NULL, ...)
\method{formatPopulation}{data.frame}(popdata, aggregate.by=NULL, breaks=NULL,...)
\method{formatPopulation}{list}(popdata, aggregate.by=NULL, breaks = NULL, 
  years=as.integer(names(popdata)), year.range=NULL,  time="YEAR", personYears=TRUE, 
  S=c("M", "F"),...)
%formatPopulation.PostgreSQLConnection(dbname, user="postgres",driver="PostgreSQL",aggregate=FALSE)
\method{formatPopulation}{SpatialPolygonsDataFrame}(popdata, ...)
}

\arguments{
  \item{popdata}{ population data set. It can be a data frame, list, database connection, or spatial polygon data frame }
%  \item{dbname}{name of data base}
  \item{aggregate.by}{ if want to view the data set from a macro way, 
     could aggregate the data set by age or sex or other variables }
%  \item{user}{}
  \item{breaks}{ age breaks the user want to use. i.e breaks = c(10, 20, 30 ,40, 60, Inf). }
%  \item{driver}{}
  \item{years}{the years vector}
  \item{aggregate}{if \code{TRUE} then aggregate the population data}
  \item{year.range}{the range of the year vector}
  \item{time}{the time variable, i.e years}
  \item{personYears}{convert populations to person-years}
  \item{S}{ can choose both or either of each sex.}
  \item{\dots}{ additional arguments. }
}
\details{
  After using the \code{formatPopulation} function, it will return the population data set in the same class as the original data set. 
i.e if a spatial polygon data frame has been put into the \code{formatPopulation} function, it will return a spatial polygon data frame. 
If \code{aggregate.by} is not NULL, the number of cases will be sum up by the groups define in aggregate.by. 
The "Group" column contains information of sex and age groups,in the format of M.55, denotes male, year 55. 
The "POPULATION" column is a numeric column, denotes the size of population for the particular age and sex group.
The "age" column will be a "character" column contains ages in a cut format. i.e [50,55), denotes age 50. 
The cut breaks will get from the breaks of population data set or define by user.   
The sex column will have two levels "M" and "F" as factors.
}


\author{ Patrick Brown }
\note{
 If \code{breaks} is not specified, the function will aggregate by "age" as default.
}

\examples{
data(kentucky)
head(kentucky@data)
poptry <- formatPopulation(kentucky, breaks = c(seq(0, 80, by=10), Inf))
head(poptry)
poptryagg <- formatPopulation(kentucky, breaks = c(seq(0, 80, by=10), Inf), 
  aggregate.by=c("sex", "age"))
head(poptryagg)
}
 
