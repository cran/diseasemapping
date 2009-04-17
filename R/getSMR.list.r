getSMR.list <- function(popdata, model, casedata = NULL, regionCode = "CSDUID", 
              regionCodeCases = "CSD2006", area = FALSE, area.scale = 1, years = NULL, year.range = NULL, ...){ 
#  lennon's stuff
        isSP = (class(popdata[[1]]) == "SpatialPolygonsDataFrame")
   
    if (is.null(years)) {
        years = as.integer(names(popdata))
    }        
    if (area & isSP) {
       
            areas <- lapply(lapply(popdata, area), as.numeric)
            for (i in 1:length(popdata)) {
                popdata[[i]]$sqk <- areas[[i]] * area.scale
            }
    }
    
    poplong <- formatPopulation(popdata, breaks = attributes(model)$breaks$breaks,
      years = model$xlevels$YEAR, mustAggregate = FALSE, year.range=year.range)

       
    ll<-split(poplong,poplong$YEAR)     
         
    ##list if df
    listpop<-lapply(ll, getSMR, casedata=casedata, model, regionCode =regionCode,
                     regionCodeCases = regionCodeCases, years = years, year.range = year.range,
                     area = area, area.scale = area.scale,formatPop=FALSE)
                                
         #if input is  list of sp, return list of sp
        if (isSP) {
            for (i in 1:length(years)) {
                popdata[[i]]@data = listpop[[i]]
            }
            listpop<-popdata
         }   
        #else return list of df
  
    listpop
}

