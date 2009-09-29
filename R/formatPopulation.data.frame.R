formatPopulation <- function(popdata, aggregate.by=NULL, breaks=NULL, ...) {
  UseMethod("formatPopulation")
}

`formatPopulation.data.frame` <-
function(popdata, aggregate.by=NULL, breaks=NULL,...) {

#popdata <- popdata@data
ageBreaks = getBreaks(names(popdata), breaks)

####reshape the popdata:
poplong = reshape(popdata,  varying=ageBreaks$oldNames, direction="long",
	v.names="POPULATION", timevar="GROUP", times = ageBreaks$newNames)
# create age and sex variables
agecol = grep("^age$", names(poplong), value=TRUE, ignore.case=TRUE)
sexcol = grep("^sex$", names(poplong), value=TRUE, ignore.case=TRUE)


if("GROUP" %in% names(poplong)) {
    if(!length(sexcol)){
 
  poplong$sex = factor(toupper(substr(poplong$GROUP, 1, 1)))}
     if(!length(agecol)){
   ageNumeric = as.numeric(substr(poplong$GROUP, 3, 4))
   poplong$age = cut(ageNumeric, ageBreaks$breaks, right=FALSE)
  }else {
  warning("no age and sex variables found or no group variable found in popdata")
  }
}

poplong$id<-NULL
row.names(poplong)<-NULL
#names(poplong)<-toupper(names(poplong))

if(!is.null(aggregate.by)) {

  popa <- aggregate(poplong$POPULATION, poplong[, aggregate.by, drop=FALSE], 
    sum, na.rm=TRUE)

  # change x column name to 'population'
  names(popa)[names(popa)=="x"] = "POPULATION"
  names(popa)[names(popa)=="poplong[, aggregate.by]"] = aggregate.by
  poplong<-popa

}

if(length(sexcol)) poplong[,sexcol] <- toupper(poplong[,sexcol])
attributes(poplong)$breaks = ageBreaks





#poplong <- poplong[!is.na(poplong$POPULATION),  ]
poplong



}
