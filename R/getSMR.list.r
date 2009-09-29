getSMR.list <- function(popdata, model, casedata = NULL, regionCode = "CSDUID", 
              regionCodeCases = "CSD2006", area = FALSE, area.scale = 1, 
              years = NULL, personYears=TRUE,year.range = NULL,...){
#  lennon's stuff
        #isSP = (class(popdata[[1]]) == "SpatialPolygonsDataFrame")
   
    if (is.null(years)) {
        years = as.integer(names(popdata))
    }        

    if(personYears){
    year.range = range(names(popdata))
    times <- c(year.range[1], sort(years), year.range[2])
    times <- as.numeric(times)
    inter <- diff(times)/2
    nseq <- 1:length(inter) - 1
    mseq <- 2:length(inter)
    interval <- inter[mseq] + inter[nseq]
    names(interval) <- names(popdata)
    }else{
     interval<-rep(1,length(popdata))
    }

    names(interval) <- names(popdata)
      
  yearVar = grep("year",names(attributes((terms(model)))$dataClasses) ,value=TRUE,ignore.case=TRUE)  # find year var in the model
  #If year is not in the model as factor
  if(length(yearVar)==0) yearVar="YEAR"
  
  caseYearVar = grep("year",names(casedata),value=TRUE,ignore.case=TRUE)
  caseSexVar =grep("^sex$",names(casedata),value=TRUE,ignore.case=TRUE)
  if (length(model$sexSubset) == 1) warning("only one sex is being used:",model$sexSubset)


  
  result=list()

  #for(Dyear in seq(1, length(popdata))) {
  for(Dyear in names(popdata)) {
  
    #Compute area
    if(area & class(popdata[[Dyear]]) == "SpatialPolygonsDataFrame") {
      areas <- unlist(computeArea(popdata[[Dyear]]))
      popdata[[Dyear]]$sqk <- areas * area.scale
    }

  
    # scaling factor to convert to person years
    attributes(popdata[[Dyear]]@data)$popScale = interval[names(interval)==Dyear]
    
    # add year column
    #popdata[[Dyear]][[yearVar]] = factor(rep(as.character(names(popdata[Dyear])), length(popdata[[Dyear]][1]),
    #levels=levels(model$data[,yearVar])))
    
    popdata[[Dyear]][[yearVar]] = as.integer(Dyear)


    caseThisYear = casedata[casedata[[caseYearVar]]==Dyear,]

    if (length(model$sexSubset) == 1) {      
        caseThisYear = caseThisYear[caseThisYear[[caseSexVar]] == model$sexSubset, ]
    }

    if(dim(caseThisYear)[1]==0) caseThisYear<-NULL   

    cat("computing",Dyear,"\n")

    result[[Dyear]] = getSMR(popdata[[Dyear]], model, caseThisYear,regionCode =regionCode,
                     regionCodeCases = regionCodeCases, years = years, year.range = year.range,
                     area = area, area.scale = area.scale)
    
  }
  names(result) = years

         
  return(result)       

}

