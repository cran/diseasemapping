#years: a vector of census years, default will use the names of popdata list
#year.range: A two-element numerical vector indicting the starting and ending year of study period, default set to range of casedata
#case.years: the variable name that stores the year variable of case file
#' @export
`getRates` <-
function(casedata, popdata, formula, family='poisson', minimumAge=0,
   maximumAge=100, S=c("M", "F"), years=NULL, year.range=NULL,
   case.years=grep("^year$", names(casedata), ignore.case=TRUE, value=TRUE),
   fit.numeric=NULL ,breaks=NULL){
 
# check the formula is one sided
  
if(attributes(stats::terms(formula))$response) {
	casecol = rownames(attributes(stats::terms(formula))$factors)[
			attributes(stats::terms(formula))$response
			]
} else {
	casecol = NULL
}
attributes(casedata)$casecol = casecol

morethanoneyear = class(popdata)=="list"



#is SP or not
if(morethanoneyear){
  isSP = (class(popdata[[1]])== "SpatVector")
}else{
  isSP = (class(popdata)== "SpatVector")
}

#if years not supplied, use the names of list
if(is.null(years) & morethanoneyear ){
  years = as.integer(names(popdata))
}

#factors we need to aggregate by
theterms=setdiff(gsub("^s\\(|,([[:alnum:]]|=|[[:space:]]|,|\\$|\\[|\\])+\\)$|\\)$", "", 
		rownames(attributes(stats::terms(formula))$factors)), casecol)

if(all(c('population',theterms) %in% colnames(popdata))) {
	# data is already long
	pops = stats::aggregate(
			popdata[,'population',drop=FALSE],
			popdata[,theterms],
			sum,
			na.rm=TRUE
			)
	theBreaks = sort(unique(pops$age))
	pops$age = cut(pops$age, theBreaks, right=FALSE)
	attributes(pops)$breaks = list(
			breaks = theBreaks,
			sex = unique(pops$sex)
			)
			
} else {
	pops <- formatPopulation(popdata, aggregate.by= theterms, 
			breaks=breaks, personYears=FALSE,S=S)
}


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
# add covariates to case data from the population data
# if they are any missing from the case data.

# check to see if any of the terms in the model aren't in the case data
termsToAdd = NULL
for(D in theterms) {
	if(!length(grep(D, names(casedata), ignore.case=TRUE)))
		termsToAdd = c(termsToAdd, D)
}
	

if(length(termsToAdd) ) {
  if(morethanoneyear)
    warning("All covariates must be added to case data if popdata is a list")
   colsToTry = popdata[,-grep("^(M|F)[[:digit:]]", names(popdata))]
     if(isSP) colsToTry = terra::values(colsToTry)

  # find column to merge on
  commonCol = names(casedata)[(names(casedata) %in% names(popdata))]
  if(!length(commonCol)) {
     # no common columns, look for factors
     # find columns with at least as many different entries as terms
    Npop = dim(colsToTry)[1]
    Nincases = apply(casedata, 2, function(qq) length(unique(qq)))   
    Ninpop = apply(colsToTry, 2, function(qq) length(unique(qq)))
    
    colCase = names(Nincases)[order(abs(Nincases - Npop),decreasing=FALSE)[1]]
    colPop = names(Ninpop)[order(abs(Ninpop - Npop),decreasing=FALSE)[1]]
    
     
  } else {
    colPop = colCase = commonCol[1]
  }
  message(paste0("adding variable ", toString(termsToAdd), " to case data using \n",
    "columns ", colPop, " and ", colCase, ".\n", 
   "add variable to casedata manually if this is not correct\n"))
    
  casedata = merge(casedata, colsToTry[,c(colPop,termsToAdd)], 
    by.x=colCase, by.y=colPop)

} # end if termsToAdd
 
casedata = formatCases(casedata, 
		ageBreaks=attributes(pops)$breaks, 
  aggregate.by = theterms)

casecol = attributes(casedata)$casecol

##### merge case data set and shape data set according to the same Year and DA2001
by.x =  paste("^", theterms, "$", sep="")
by.x = paste(by.x, collapse="|")
by.x = paste("(", by.x, ")", sep="")
by.pop = grep(by.x, names(pops), ignore.case=TRUE, value=TRUE)

##make them same order
by.pop<-by.pop[order(by.pop)]
theterms<-theterms[order(theterms)]

 
newdata <- merge(casedata, pops, by.x = theterms, by.y = by.pop,all.x=TRUE)

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
} # end more than one year

popCol = grep("^population$", colnames(newdata), ignore.case=TRUE, value=TRUE)
newdata = newdata[!is.na(newdata[,popCol]), ]
newdata = newdata[newdata[,popCol]>0,]
newdata$logpop = log(newdata[,popCol])

newdata[is.na(newdata[,casecol]),casecol] <-0
	

# make the age group with the most cases as the base line
agevar =  grep("^age$", theterms, ignore.case=TRUE, value=TRUE)
if(length(agevar)==1) {
  agetable = tapply(newdata[,casecol], newdata[[agevar]], sum,na.rm=TRUE)
  agetable = names(sort(agetable, decreasing=TRUE))
  newdata[[agevar]] = factor(as.character(newdata[[agevar]]),levels= agetable)
}


sexvar = grep("^sex", theterms, ignore.case=TRUE, value=TRUE)
if(length(sexvar) == 1){
newdata[[sexvar]] = factor(newdata[[sexvar]])
}

#change factor to numeric
if(!is.null(fit.numeric)){
    for (i in 1:length(fit.numeric)){
    toChange = grep(paste("^",fit.numeric[i],"$",sep=""),names(newdata),value=TRUE,
      ignore.case=TRUE)
    newdata[,toChange] = as.numeric(as.character(newdata[,toChange]))
  }
}


#todel <- as.formula(paste(".~.-",sexvar,"-",agevar,":",sexvar,sep=""))
#if(length(S)==1) formula=update.formula(formula, todel)

# add cases and logpop to formula
formula1 = stats::update.formula(
    formula, 
		stats::as.formula(paste(casecol," ~ offset(logpop) + ."))
)
#return(newdata, formula1)

#if s() is in formula, use GAM
useGam <- length(grep("s\\(",formula))>0
if(useGam & requireNamespace("mgcv", quietly = TRUE)) {
  model <-  mgcv::gam(formula = formula1, family = family, data=newdata)
} else {
  model <-  stats::glm(formula = formula1, family = family, data=newdata)
}

model$sexSubset = S
model$data<-newdata
#attributes(model)$years = ageBreaks$breaks
attributes(model)$breaks = attributes(pops)$breaks
                                  

return(model)
}

