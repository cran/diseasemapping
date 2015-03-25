\name{popdata}
\alias{popdata}
\docType{data}
\title{Ontario 2006 population by census subdivision}
\description{
 Data set contains the information of population, by age, sex, and census subdivision.
}
\usage{data(popdata)}
\format{
A SpatialPolygonsDataFrame object, which needs the \code{sp} package for full functionality.  
}
\details{
This data is from the 2006 Census of canada offering by Statistics Canada web site, \code{www12.statcan.gc.ca/english/census06/data/highlights/agesex/Index_PR.cfm?Lang=E&Geo=CSD&Table=1} 
}


\examples{
data(popdata)
head(popdata@data)
\dontrun{
library(sp)
spplot(popdata, zcol='F.50_54', breaks=9, col=rainbow(8))
}
\dontrun{

library('raster')
library('sp')
bfile = tempfile(fileext='.zip')
download.file(
paste('http://www12.statcan.gc.ca/census-recensement/',
'2011/geo/bound-limit/files-fichiers/gcsd000a06a_e.zip',
sep=''),
bfile)
unzip(bfile, exdir=tempdir())
sfile = grep('shp$',unzip(bfile, list=TRUE)$Name, value=TRUE)
popdata = shapefile(file.path(tempdir(),sfile))
popdata$PRNAME = iconv(popdata$PRNAME, 'UTF-8', 'latin1')
popdata = popdata[grep("^Ont", popdata$PRNAME),]
popdataS= rgeos::gSimplify(popdata, 0.01, topologyPreserve=TRUE)
popdata = SpatialPolygonsDataFrame(popdataS, popdata@data)
projection(popdata) = CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')

pfile = tempfile(fileext='zip')

download.file(
paste('https://www12.statcan.gc.ca/census-recensement/',
'2011/dp-pd/prof/details/download-telecharger/comprehensive/',
'comp_download.cfm?CTLG=92-591-XE&FMT=CSV301&Lang=E&Tab=1&',
'Geo1=PR&Code1=01&Geo2=PR&Code2=01&Data=Count&SearchText=&',
'SearchType=Begins&SearchPR=01&B1=All&Custom=&TABID=1', sep=''),
pfile, method='curl')
unzip(pfile, exdir=tempdir())
ofile = grep('ONT', unzip(pfile, list=TRUE)$Name,value=TRUE)
opop = read.table(file.path(tempdir(),ofile),header=F,skip=3, 
sep=',', nrows=163210,stringsAsFactors=FALSE)
opop= opop[grep("^([[:digit:]]|to| )+ years( and over)?$", opop[,7]),]
opop = opop[,c(1,4,7,11,13)]
colnames(opop) = c('id','name','var','M','F')
opop[,'var'] = gsub(" to ", "_", opop[,'var'])
opop[,'var'] = gsub(" years( and over)?", "", opop[,'var'])
opop[,'var'] = gsub("[[:space:]]", "", opop[,'var'])
opop2 = reshape(opop, direction='wide', 
idvar=c('id','name'),
timevar='var', v.names=c('M','F'))


popdata = sp::merge(popdata, opop2, by.x='CSDUID', by.y='id')
popdata=popdata[,c('CSDUID', grep("^(M|F)", names(popdata), value=TRUE))]



save(popdata, file=
'/home/patrick/workspace/diseasemapping/pkg/diseasemapping/data/popdata.RData',
compress='xz')
}
}
\keyword{datasets}
