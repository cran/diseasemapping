`getSMR` <- function(popdata, model, casedata, regionCode,
    regionCodeCases, area = FALSE, area.scale = 1, ...){
       UseMethod("getSMR")
 }
 
getSMR.data.frame <- function(popdata, model, casedata=NULL, 
		regionCode = intersect(names(popdata), names(casedata))[1],
    regionCodeCases=regionCode, area=FALSE, area.scale=1, ...){




    if(is.numeric(model)) {
    # model is a vector of rates
        # check breaks for groups, make sure they line up
        rateBreaks =getBreaks(names(model))
        popBreaks = getBreaks(names(popdata))

#return(list(r=rateBreaks, p=popBreaks))
        noPop = ! popBreaks$newNames %in% rateBreaks$newNames
        if(any(noPop))
          warning(paste("population group(s)", toString(popBreaks$oldNames[noPop]),
           "ignored") )
        noRate = ! rateBreaks$newNames %in% popBreaks$newNames  
        if(any(noRate))
          warning(paste("rate group(s)", toString(rateBreaks$oldNames[noRate]),
           "ignored\n") )
  
        popGroups = popBreaks$oldNames[!noPop]
        names(popGroups) = popBreaks$newNames[!noPop]

        rateGroups = rateBreaks$oldNames[!noRate]
        names(rateGroups) = rateBreaks$newNames[!noRate]

       
        popdata$expected = as.vector(
            as.matrix(popdata[,popGroups]) %*% model[rateGroups[names(popGroups)]]
          )
	if(!is.null(regionCode)) {
        rownames(popdata) = as.character(popdata[,regionCode])
    }
	
    } else {
    # use the predict method on the model

    poplong <- formatPopulation(popdata, breaks=attributes(model)$breaks$breaks)
      
    p<-grep("^population$", names(poplong), value=TRUE, ignore.case=TRUE)  
     

    poplong[is.na(poplong[,p]),p] <- 0
    
    popBreaks = attributes(poplong)$breaks
     
     # get rid of zero populations,because they dont lead to rates of exactly zero
    ## check the class of the POPULATION column and make sure it is numeric: 
	# converted to character first in case it's a factor
		# (convert factor names to integers, not factor levels)
 
    if(class(poplong[,p]) != "numeric"){
     poplong[,p] <- as.numeric(as.character(poplong[,p]))
    }

     poplong=poplong[poplong[,
      grep("^population$", names(poplong), value=TRUE, ignore.case=TRUE)]>0, ]     
    #changes poplong names to be consistent with model
    agevar<-grep("^age$",names(attributes((terms(model)))$dataClasses),value=TRUE,ignore.case=TRUE)
    sexvar<-grep("^sex$",names(attributes((terms(model)))$dataClasses),value=TRUE,ignore.case=TRUE)
    yearvar<-grep("^year$",names(attributes((terms(model)))$dataClasses),value=TRUE,ignore.case=TRUE)

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

      # get rid of a gender
    if (length(model$sexSubset) == 1) {
        poplong = poplong[poplong[[sexvar1]] == model$sexSubset, ]
    }
    
      
      offsetvar<- grep("^offset",names(attributes((terms(model)))$dataClasses)  ,value=TRUE, ignore.case=TRUE)
       offsetvar<- substr(offsetvar,8,nchar(offsetvar)-1)
      
      
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
    if(length(yearvar1)) {
           agg<-c(yearvar1,agevar, sexvar,offsetvar)
           }else{
                agg<-c(agevar, sexvar, offsetvar)
           }

    # multiply population by popScale, to make it in person years
    if(any(names(attributes(popdata))=="popScale")) {
      poplong[,offsetvar]=     poplong[,offsetvar] + log(attributes(popdata)$popScale)
    }
    

    poplong$expected <- predict(model, poplong,#[,agg],
     type = "response")
    
     poplong <- aggregate(poplong$expected, list(poplong[[regionCode]]), sum)
    rownames(poplong) = as.character(poplong[,1])
    poplong=poplong[poplong[,2] > 0,]


    # merge results back in to the population data
    # the merge function changes the order, so can't use it.
    popdata$expected = NA
    rownames(popdata) = as.character(popdata[,regionCode])

    popdata[rownames(poplong), "expected"] = poplong[,2]

    } # done predicting rates from model
    
    
    if (area & ("sqk" %in% names(popdata) ) ) {
        popdata$expected_sqk <- popdata$expected/popdata$sqk
        popdata$logExpected_sqk = log(popdata$expected_sqk)
    }

    
    popdata$logExpected = log(popdata$expected)
    popdata$observed = 0
    # change NA's and -Inf in logExpected to zeros, so that it can be used in models
    popdata$logExpected[is.na(popdata$logExpected)] = 0
    popdata$logExpected[popdata$logExpected==-Inf] = 0
    

    if (!is.null(casedata)) {

# find column with cases
		casecol = grep("^cases$|^count$|^y$", names(casedata), value=TRUE, ignore.case=TRUE)
		if(length(casecol)>1) {
			casecol=casecol[1]
			warning("more than one column which could be interpreted as case numbers, using ", casecol)
		}
		
		if(!length(casecol)) {
			#there is no case col
			casecol = "cases"
			casedata[,casecol] = 1
		}
		
       casedata = casedata[
          as.character(casedata[, regionCodeCases]) %in% 
             rownames(popdata), ]
	   
      casedata <- aggregate(casedata[[casecol]], 
         list(casedata[[regionCodeCases]]), sum)
       names(casedata) = c(regionCodeCases, "observed")


      popdata[as.character(casedata[,1]),"observed"] = casedata[,2]
   
      # change 0's in expected to NA, so SMR is NA
    theexpected = popdata$expected
    theexpected[theexpected==0] = NA
    
     popdata$SMR <- popdata$observed/theexpected
   }

   popdata
}
