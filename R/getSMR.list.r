getSMR.list <- function(popdata, model, casedata = NULL, regionCode = "CSDUID", 
              regionCodeCases = "CSD2006", area = FALSE, area.scale = 1, 
              years = NULL, year.range = NULL, ...){ 
#  lennon's stuff
        #isSP = (class(popdata[[1]]) == "SpatialPolygonsDataFrame")
   
    if (is.null(years)) {
        years = as.integer(names(popdata))
    }        


    year.range = range(names(popdata))
    times <- c(year.range[1], sort(years), year.range[2])
    times <- as.numeric(times)
    inter <- diff(times)/2
    nseq <- 1:length(inter) - 1
    mseq <- 2:length(inter)
    interval <- inter[mseq] + inter[nseq]
    names(interval) <- names(popdata)

  yearVar = grep("year",names(model$xlevels),value=TRUE,ignore.case=TRUE)  # find year var in the model
  caseYearVar = grep("year",names(casedata),value=TRUE,ignore.case=TRUE)
  
  result=list()

  #for(Dyear in seq(1, length(popdata))) {
  for(Dyear in names(popdata)) {
  
    #Compute area
    if(area & class(popdata[[Dyear]]) == "SpatialPolygonsDataFrame") {
      areas <- unlist(computeArea(popdata[[Dyear]]))
      popdata[[Dyear]]$sqk <- areas * area.scale
    }

  
    # scaling factor to convert to person years
    attributes(popdata[[Dyear]])$popScale = interval[names(interval)==Dyear]
    
    # add year column
    popdata[[Dyear]][[yearVar]] = factor(rep(as.character(names(popdata[Dyear])), length(popdata[[Dyear]][1]),
    levels=levels(model$data[,yearVar])))
    


    caseThisYear = casedata[casedata[[caseYearVar]]==Dyear,]
    
    result[[Dyear]] = getSMR(popdata[[Dyear]], model, caseThisYear,regionCode =regionCode,
                     regionCodeCases = regionCodeCases, years = years, year.range = year.range,
                     area = area, area.scale = area.scale)
    
  }
  names(result) = years

         
  return(result)       

}

