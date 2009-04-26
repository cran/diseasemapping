#years: a vector of census years, default will use the names of popdata list
#year.range: A two-element numerical vector indicting the starting and ending year of study period, default set to range of casedata
#case.years: the variable name that stores the year variable of case file

`getRates` <-
function(casedata, popdata, formula, family=poisson, minimumAge=0,
   maximumAge=100, S=c("M", "F"), years=NULL, year.range=NULL,
   case.years=grep("^year$", names(casedata), ignore.case=TRUE, value=TRUE)[1], breaks=NULL){

# check the formula is one sided
  
if(attributes(terms(formula))$response)
  warning("formula should be one sided")

morethanoneyear = class(popdata)=="list"

#is SP or not
if(morethanoneyear){
  isSP = (class(popdata[[1]])== "SpatialPolygonsDataFrame")
}else{
  isSP = (class(popdata)== "SpatialPolygonsDataFrame")
}

#if years not supplied, use the names of list
if(is.null(years) & morethanoneyear ){
  years = as.integer(names(popdata))
}

#factors we need to aggregate by
theterms = (rownames(attributes(terms(formula))$factors))

pops <- formatPopulation(popdata, aggregate.by= theterms, breaks=breaks, getoff=FALSE) 
          

##format case data
#casedata = formatCases(casedata, ageBreaks=attributes(pops)$breaks, aggregate.by = theterms)

#if ranges not supplied, use the year ranges of case files
if(is.null(year.range) & morethanoneyear){
  year.range = range(casedata[[case.years]])
} 



# keep only the desired sex in the dataset
if(length(S)==1) {
  if(length(grep("^sex$", theterms, ignore.case=TRUE)))
    warning("sex is in the model but only one sex is being used")
  casedata=casedata[casedata[[
    grep("^sex$", names(casedata), value=TRUE, ignore.case=TRUE)
      ]]==S,]
}

#format case data

casedata = formatCases(casedata, ageBreaks=attributes(pops)$breaks, aggregate.by = theterms)


#find number of cases per group
casecol = grep("^cases$", names(casedata), value=TRUE, ignore.case=TRUE)
if(!length(casecol)) {
  #there is no case col
  casecol = "cases"
  casedata[,casecol] = 1
}


##### merge case data set and shape data set according to the same Year and DA2001
by.x =  paste("^", theterms, "$", sep="")
by.x = paste(by.x, collapse="|")
by.x = paste("(", by.x, ")", sep="")
by.pop = grep(by.x, names(pops), ignore.case=TRUE, value=TRUE)


newdata <- merge(casedata, pops, by.x = theterms, by.y = by.pop)
if (morethanoneyear){
####find Popoluation census year
#a vector of all years


times<-c(year.range[1],sort(years),year.range[2])

inter<-diff(times)/2 #mid points
#sum of consective mid points
nseq<-1:length(inter)-1
mseq<-2:length(inter)
interval<-inter[mseq] + inter[nseq]

names(interval)<-names(table(newdata$YEAR))         
newdata$yearsForCensus = interval[as.character(newdata$YEAR)]
newdata$POPULATION = newdata$POPULATION  * newdata$yearsForCensus 
newdata$YEAR= factor(newdata$YEAR, levels = unique(newdata$YEAR))
}


newdata = newdata[newdata$POPULATION>0,]
newdata$logpop = log(newdata$POPULATION)
	

# make the age group with the most cases as the base line
agevar =  grep("^age$", theterms, ignore.case=TRUE, value=TRUE)
if(length(agevar)==1) {
  agetable = tapply(newdata$CASES, newdata[[agevar]], sum)
  agetable = names(sort(agetable, decreasing=TRUE))
  newdata[[agevar]] = factor(as.character(newdata[[agevar]]),
    levels= agetable)
}


sexvar = grep("^sex", theterms, ignore.case=TRUE, value=TRUE)
if(length(sexvar) == 1){
newdata[[sexvar]] = factor(newdata[[sexvar]])
}

# add cases and logpop to formula
formula1 = update.formula(formula, CASES ~ offset(logpop) + .)
#return(newdata, formula1)

model = glm(formula1, family=family, data=newdata)

model$sexSubset = S

#attributes(model)$years = ageBreaks$breaks
attributes(model)$breaks = attributes(pops)$breaks
                                  

return(model)
}

