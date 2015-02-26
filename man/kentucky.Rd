\name{kentucky}
\alias{kentucky}
\alias{larynx}
\alias{larynxRates}
\docType{data}
\title{Larynx cancer cases and population in Kentucky}
\description{
 Data set contains the information of population, by age, sex, and census subdivision.
}
\usage{data(kentucky)}
\format{A SpatialPolygonsDataFrame of Kentucky county boundaries and populations, a data frame 
of case numbers for each county, and a vector of cancer rates by age and sex group.

}
\details{
The cancer case data are obtained from www.cancer-rates.info and are for a single  
deliberately unspecified year. 
}



\examples{
data('kentucky')

head(larynx)
10^5*larynxRates[paste(c("M","F"), 50, sep="_")]


kentucky2 = getSMR(kentucky, larynxRates, larynx, 
		regionCode="County")

names(kentucky2)

if(require('mapmisc', quietly=TRUE)) {
	mycol = colourScale(kentucky2$SMR, 
		breaks=10, dec=-log10(0.5), style='quantile')
	map.new(kentucky2)
	plot(kentucky2, col=mycol$plot, border='#00000040',add=TRUE)
	legendBreaks('topright', mycol)
	scaleBar(kentucky2,"topleft")
} else {
	plot(kentucky2)
}

breaks = c(0,1,seq(2, ceiling(max(kentucky2$SMR,na.rm=TRUE)),by=2))
thecol = terrain.colors(length(breaks)-1)

plot(kentucky2, col = thecol[cut(kentucky2$SMR, 
	breaks,include.lowest=TRUE)] )

legend("topleft", pch=15, pt.cex=2.5, adj=c(0,15), 
  legend=rev(breaks), col=c(NA, rev(thecol)))

\dontrun{
	# the data were created with
	larynxRates= cancerRates("USA", year=1998:2002,site="Larynx")
	
	load(url("http://biogeo.ucdavis.edu/data/gadm2/R/USA_adm2.RData"))
	kentucky = gadm[gadm$NAME_1 =="Kentucky",]

	# population data
	download.file(
	"http://www.census.gov/popest/data/counties/asrh/2011/files/CC-EST2011-ALLDATA-21.csv",
		destfile = 
	"/store/patrick/spatialData/C-EST2011-ALLDATA-21.csv")
	# file layout
	download.file(
	"http://www.census.gov/popest/data/counties/asrh/2011/files/CC-EST2011-ALLDATA.pdf",
		destfile = "/store/patrick/spatialData/kentuckyPopFormat.pdf")

	kpop = read.table(
	"/store/patrick/spatialData/C-EST2011-ALLDATA-21.csv",
		header=TRUE,as.is=TRUE,sep=",")	
	kpop = kpop[kpop$YEAR==1 & kpop$AGEGRP != 0, ]
	names(kpop) = gsub("^TOT_","", names(kpop))
	names(kpop) = gsub("(EM)?ALE$","", names(kpop))
	 
	kpop$age = (kpop$AGEGRP-1)*5 
	kpop$County =  gsub(" County$", "", kpop$CTYNAME)
	kpop = kpop[,c("County","age","M","F")]
	kpop2 = reshape(kpop,direction="wide", idvar="County",
		v.names=c("M","F"), timevar="age")
	rownames(kpop2) = kpop2$County

		
	# poverty
	download.file(
	paste(
	"http://www.ers.usda.gov/ReportExport.aspx?reportPath=/State_Fact_Sheets/",
	"PovertyReport&stat_year=2011&stat_type=0&fips_st=21&",
	"exportType=EXCEL&exportName=PovertyReport",
	sep=""),
	destfile="/store/patrick/spatialData/poverty.xls")
	library('gdata')
	kpov = read.xls("/store/patrick/spatialData/poverty.xls",
		header=TRUE,skip=3)
	kpov = kpov[!is.na(kpov$Percent),c("FIPS.", "Name","Percent")]
	rownames(kpov) = kpov$Name
	kpop2$poverty = kpov[rownames(kpop2), "Percent"]
	
	# merge population and spatial data
	kdata = kpop2[kentucky$NAME_2,]
	rownames(kdata) = NULL
	kentucky = SpatialPolygonsDataFrame(
		polygons(kentucky),
		data=kdata,match.ID=FALSE)

larynx <- structure(
list(County = c("Hickman", "Caldwell", "Anderson", 
"Marion", "Wayne", "Lincoln", "Livingston", "Montgomery", "Adair", 
"Henderson", "Knox", "Martin", "Monroe", "Wolfe", "Breathitt", 
"Fleming", "Woodford", "Garrard", "Bracken", "Barren", "Lawrence", 
"Logan", "Clark", "Scott", "Madison", "Oldham", "Clay", "Russell", 
"Shelby", "Letcher", "Campbell", "Graves", "Johnson", "Metcalfe", 
"Pulaski", "Bullitt", "Knott", "Boyd", "Ohio", "Bath", "Butler", 
"Todd", "Mercer", "Green", "Greenup", "Larue", "Calloway", "Webster", 
"Morgan", "Pendleton", "Mason", "Hardin", "Lewis", "McCreary", 
"Spencer", "Union", "Marshall", "Jessamine", "Henry", "Trigg", 
"Pike", "Nelson", "Jefferson", "Floyd", "Bourbon", "McCracken", 
"Boone", "Kenton", "Grayson", "Taylor", "Hopkins", "Boyle", "Meade", 
"Fayette", "Daviess", "Harlan", "Warren", "Christian", "Magoffin", 
"Carter", "Hart", "Lee", "Elliott", "Edmonson", "Crittenden", 
"Leslie", "Laurel", "Cumberland", "Menifee", "Fulton", "Carlisle", 
"McLean", "Owsley", "Carroll", "Estill", "Harrison", "Owen", 
"Breckinridge", "Nicholas", "Bell", "Trimble", "Allen", "Rowan", 
"Simpson", "Perry", "Powell", "Rockcastle", "Hancock", "Robertson", 
"Franklin", "Washington", "Casey", "Clinton", "Lyon", "Muhlenberg", 
"Ballard", "Gallatin", "Whitley", "Grant", "Jackson", "Breathitt", 
"Nicholas", "Bracken", "Todd", "Magoffin", "Pendleton", "Metcalfe", 
"Webster", "Leslie", "Henry", "Union", "Adair", "Casey", "Pike", 
"Jessamine", "Nelson", "Garrard", "Pulaski", "Meade", "Harlan", 
"Floyd", "Carter", "Shelby", "Barren", "Franklin", "Boyd", "Jefferson", 
"Fayette", "Hopkins", "Kenton", "Warren", "Bullitt", "Knox", 
"Butler", "Bourbon", "Elliott", "Johnson", "Estill", "Boone", 
"Boyle", "Breckinridge", "Bell", "Crittenden", "Cumberland", 
"Daviess", "Edmonson", "Calloway", "Caldwell", "Anderson", "Ballard", 
"Bath", "Allen", "Graves", "Clinton", "Fleming", "Fulton", "Gallatin", 
"Grayson", "Hardin", "Lincoln", "Green", "Greenup", "Hancock", 
"Grant", "Harrison", "Laurel", "Larue", "Henderson", "Hickman", 
"Lewis", "Letcher", "Hart", "Lawrence", "Lee", "Jackson", "Ohio", 
"Taylor", "Owen", "Monroe", "Madison", "Trigg", "Spencer", "Washington", 
"Trimble", "Scott", "Simpson", "Mason", "Marion", "Lyon", "Logan", 
"McCracken", "Mercer", "Menifee", "McCreary", "Marshall", "Martin", 
"McLean", "Woodford", "Wolfe", "Knott", "Montgomery", "Perry", 
"Powell", "Christian", "Clark", "Campbell", "Carlisle", "Carroll", 
"Robertson", "Rockcastle", "Rowan", "Russell", "Wayne", "Whitley", 
"Oldham", "Muhlenberg", "Owsley", "Livingston", "Morgan", "Clay"
), 
Cases = c(2, 3, 3, 3, 4, 3, 1, 3, 3, 5, 3, 2, 2, 1, 2, 2, 
2, 2, 1, 4, 2, 3, 4, 3, 7, 4, 2, 2, 3, 2, 6, 3, 1, 1, 5, 4, 1, 
3, 1, 1, 1, 1, 1, 1, 3, 1, 2, 1, 1, 1, 1, 4, 1, 1, 1, 1, 2, 2, 
1, 1, 3, 1, 28, 2, 1, 2, 3, 5, 1, 1, 2, 1, 1, 6, 2, 1, 2, 1, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 2, 2, 1, 3, 1, 1, 1, 
1, 1, 1, 1, 1, 11, 3, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
sex = c("M", 
"M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", 
"M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", 
"M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", 
"M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", 
"M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", 
"M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", 
"M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", 
"M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", 
"M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", 
"M", "M", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", 
"F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", 
"F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", 
"F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", 
"F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", 
"F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", 
"F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", 
"F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", 
"F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", 
"F", "F", "F", "F", "F")), 
.Names = c("County", "Cases", "sex"), 
row.names = 1:240, class = "data.frame")
larynx$age=NA

save(kentucky, larynx, larynxRates, 
	file="~/workspace/diseasemapping/pkg/diseasemapping/data/kentucky.RData",
	compress='xz')
}
}
\keyword{datasets}
