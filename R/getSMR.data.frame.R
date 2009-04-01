`getSMR` <- function(popdata, model, casedata, regionCode = "CSDUID",
    regionCodeCases = "CSD2006", area = FALSE, area.scale = 1, ...){
       UseMethod("getSMR")
 }
 
getSMR.data.frame <- function(popdata, model, casedata, regionCode = "CSDUID",
    regionCodeCases = "CSD2006", area = FALSE, area.scale = 1, formatPop=TRUE,...){
#  getSMR(popdata@data, ...)
    if(formatPop){
    poplong <- formatPopulation(popdata, breaks=attributes(model)$breaks$breaks, 
      mustAggregate = FALSE)
     }else{poplong<-popdata}
     
    popBreaks = attributes(poplong)$breaks
     
     # get rid of zero populations,because they dont lead to rates of exactly zero
     poplong=poplong[poplong[,
      grep("^population$", names(poplong), value=TRUE, ignore.case=TRUE)]>0, ]     
    #changes poplong names to be consistent with model
    agevar<-grep("^age$",names(model$xlevels),value=TRUE,ignore.case=TRUE)
    sexvar<-grep("^sex$",names(model$xlevels),value=TRUE,ignore.case=TRUE)
    yearvar<-grep("^year$",names(model$xlevels),value=TRUE,ignore.case=TRUE)

    agevar1<-grep("^age$",names(poplong),value=TRUE,ignore.case=TRUE)
    sexvar1<-grep("^sex$",names(poplong),value=TRUE,ignore.case=TRUE)
    yearvar1<-grep("^year$",names(poplong),value=TRUE,ignore.case=TRUE)
    
    if(length(agevar) & length(agevar1)){
     names(poplong[[agevar1]])=agevar
    }
    if(length(sexvar) & length(sexvar1)){
     names(poplong[[sexvar1]])=sexvar
    }
    if(length(yearvar) & length(yearvar1)){
     names(poplong[[yearvar1]])=yearvar
    }

    
    if (length(model$sexSubset) == 1) {
        poplong = poplong[poplong$sex == model$sexSubset, ]
    }
    
      
      offsetvar<- grep("logpop",names(model$data) ,value=TRUE, ignore.case=TRUE)
    poplong[[offsetvar]] = log(poplong$POPULATION)
#    poplong[[sexvar]]= factor(poplong[[sexvar]])
#    poplong[[agevar]] = factor(poplong[[agevar]])
    #names(poplong) <- tolower(names(poplong))
    for (Dlevel in names(model$xlevels)) {
        alllevels = levels(poplong[[Dlevel]])
        if (!all(alllevels %in% model$xlevels[[Dlevel]])) {
            tokeep = poplong[[Dlevel]] %in% model$xlevels[[Dlevel]]
            poplong = poplong[tokeep, ]
        }
    }
    interactNA <- names(model$coefficients)[is.na(model$coefficients)]
    if (length(interactNA) > 0) {
        interact <- grep(":", interactNA, value = TRUE)
        poplong$param <- paste(paste("age", poplong$age, sep = ""),
            paste("sex", poplong$sex, sep = ""), sep = ":")
        poplong = poplong[!poplong$param %in% interact, ]
        poplong$param <- NULL
    }
    if(length(yearvar)) {   
           agg<-c(yearvar,agevar, sexvar,offsetvar)
           }else{
                agg<-c(agevar, sexvar, offsetvar)
           }

    poplong$expected <- predict(model, poplong[, agg], type = "response")

    
    poplong <- aggregate(poplong$expected, list(poplong[[regionCode]]), sum)
    rownames(poplong) = poplong[,1]
    poplong=poplong[poplong[,2] > 0,]


    # merge results back in to the population data
    # the merge function changes the order, so can't use it.
    popdata$expected = NA
    rownames(popdata) = popdata[,regionCode]

    popdata[rownames(poplong), "expected"] = poplong[,2]
    
    if (area & ("sqk" %in% names(popdata) ) ) {
        popdata$expected_sqk <- 
          popdata$expected/popdata$sqk
    }
    popdata$logExpected = log(popdata$expected)

    if (!is.null(casedata)) {
       casedata = formatCases(casedata, ageBreaks=popBreaks)
       casedata = casedata[
          as.character(casedata[, regionCodeCases]) %in% 
             rownames(popdata), ]
       casecol = grep("^cases$", names(casedata), value = TRUE,
            ignore.case = TRUE)
       if (!length(casecol)) {
            casecol = "cases"
            cases[, casecol] = 1
       }

       casedata <- aggregate(casedata[[casecol]], 
         list(casedata[[regionCodeCases]]), sum)
       names(casedata) = c(regionCodeCases, "observed")

      popdata$observed = NA
      popdata[as.character(casedata[,1]),"observed"] = casedata[,2]
   
    #   popdata = merge(popdata, casedata, 
     #    by.x=regionCode, by.y =regionCodeCases, all.x = TRUE)
        
       popdata$SMR <- popdata$observed/popdata$expected
   }

   popdata
}


