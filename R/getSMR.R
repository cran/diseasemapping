#' @export
setGeneric('getSMR', 
    function(popdata, model, casedata, 
        regionCode,
        regionCodeCases, 
        area.scale=1, 
        sex=c('m','f'), ...
) {
      standardGeneric("getSMR")
    }
)

setMethod("getSMR", 
		signature("SpatialPolygonsDataFrame", 'ANY','ANY', "ANY", "ANY"),
		function(
				popdata, model,  casedata, 
				regionCode,
				regionCodeCases, area.scale=1, 
				sex=c('m','f'), ...
		) {
			
			area = !length(grep("longlat", popdata@proj4string@projargs))
			
			if (area & !("surfaceArea" %in% names(popdata) ) ) {
				popdata$surfaceArea = sapply(methods::slot(popdata, "polygons"), 
						methods::slot, "area")
			}	
			
			popdata$idForGetSmr = 1:length(popdata)
			theSP = popdata
			popdata = data.frame(popdata)
			
			forSp <- methods::callGeneric()
			
			theSP@data = forSp[pmatch(
							theSP$idForGetSmr,
							forSp$idForGetSmr
					),]
			theSP = theSP[,grep("^idForGetSmr$", names(theSP), invert=TRUE)]    
			theSP
		}
)


# set regionCode and when there's no case data
setMethod("getSMR", 
		signature("data.frame", "ANY", "missing", 'missing', 'missing'),
		function(popdata, model, casedata, 
				regionCode,
				regionCodeCases, area.scale=1, 
				sex=c('m','f'), ...
		){
			regionCode = grep("^id([[:digit:]]|[[:punct:]])", names(popdata), 
					value=TRUE, ignore.case=TRUE)
			
			if(!length(regionCode))
				regionCode = names(popdata)[1]
			
			regionCode = max(regionCode)
			
			methods::callGeneric(popdata, model, 
					regionCode=regionCode,
					area.scale=area.scale, 
					sex=sex, ...)
		
		}
)


# set regionCode and regionCodeCases with case data
setMethod("getSMR", 
		signature("data.frame", "ANY", "data.frame", 'missing', 'missing'),
		function(popdata, model, casedata, 
				regionCode,
				regionCodeCases, area.scale=1, 
				sex=c('m','f'), ...
		){
			
			regionCodeCases=regionCode = intersect(names(popdata), 
					names(casedata)
			)[1]
			
			
			if(!length(regionCodeCases))
				warning('cant find region code in case data')
			
			methods::callGeneric(popdata, model, casedata, 
					regionCode,
					regionCodeCases, area.scale=1, 
					sex=c('m','f'), ...)
			
		}
)


# set regionCode present and regionCodeCases missing
setMethod("getSMR", 
		signature("data.frame", "ANY", "data.frame", 'character', 'missing'),
		function(popdata, model, casedata, 
				regionCode,
				regionCodeCases, area.scale=1, 
				sex=c('m','f'), ...
		){

			regionCodeCases=grep(
					regionCode, names(casedata), 
					value=TRUE,
					ignore.case=TRUE)[1]
			
			if(!length(regionCodeCases))
				warning('cant find region code in case data')
			
			methods::callGeneric(popdata, model, casedata, 
					regionCode,
					regionCodeCases, area.scale=1, 
					sex=c('m','f'), ...)
			
		}
)





# no case data, just get expected
setMethod("getSMR", 
		signature("data.frame", "ANY", "missing", "character", "missing"),
		function(popdata, model, casedata, 
				regionCode,
				regionCodeCases, area.scale=1, 
				sex=c('m','f'), ...
		){
			getExpected(
					popdata, model,
					area.scale,	sex
			)
		}
)


setMethod("getSMR", 
		signature("data.frame", "ANY", "data.frame", "character", "character"),
		function(popdata, model, casedata, 
				regionCode,
				regionCodeCases, area.scale=1, 
				sex=c('m','f'), ...
		){
# have case data, compute SMR

			popdata = methods::callGeneric(
					popdata, model,
					area.scale=area.scale,	sex=sex,
					...
			)
			
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
							as.character(popdata[[regionCode]]), ]
			
			casedata <- stats::aggregate(casedata[,casecol], 
					list(casedata[,regionCodeCases]), sum)
			names(casedata) = c(regionCodeCases, "observed")
			
			popdata[['observed']] = casedata[
					pmatch(
							as.character(popdata[[regionCode]]),
							as.character(casedata[[regionCodeCases]])
					),
					'observed']
			
			# change 0's in expected to NA, so SMR is NA
			theexpected = popdata$expected
			theexpected[theexpected==0] = NA
			
			popdata$SMR <- popdata$observed/theexpected
			
			popdata
		}
)

setMethod("getSMR", 
		signature("list", 'ANY', 'ANY', "ANY", "ANY"),
		function(popdata, model, casedata, regionCode,
				regionCodeCases, area.scale=1,
				sex=c('m','f'), ...){  
			
			# years = NULL, personYears=TRUE,year.range = NULL,...){
			
			dots = list(...)
			
			if (is.null(dots$years)) {
				dots$years = as.integer(names(popdata))
			}
			years = dots$years
			
			yearVar="YEAR"
			if(!is.vector(model) & class(model)[1]!='list' ) {
				yearVar = grep("year", names(attributes((stats::terms(model)))$dataClasses) ,
						value=TRUE,ignore.case=TRUE)  # find year var in the model
				if (length(model$sexSubset) == 1) {
					message("only one sex is being used:",model$sexSubset)
					sex = model$sexSubset
				}
			}
			
			if(!missing(casedata)) {
				caseYearVar = grep("year",names(casedata),value=TRUE,ignore.case=TRUE)
				caseSexVar = grep("^sex$",names(casedata),value=TRUE,ignore.case=TRUE)
			} else {
				caseSexVar = caseYearVar = NULL
			}
			
			year.range = dots$year.range
			personYears = dots$personYears
			if(!length(personYears)) {
				personYears = TRUE
			}
			
			if(personYears){
				# if year.range is missing, use year range of the cases
				if(is.null(year.range) & !missing(casedata)) {
					if(!length(caseYearVar))
						warning("year.range unspecified and no year column in case data")
					year.range = range(as.numeric(as.character(casedata[,caseYearVar])), na.rm=TRUE)
				}
				
				times <- c(year.range[1], sort(years), year.range[2])
				times <- as.numeric(times)
				inter <- diff(times)/2
				interval <- c(0, inter) + c(inter, 0)
			} else{ # not person years
				interval<-rep(1,length(popdata))
			}
			
			names(interval) <- names(popdata)
			
			result=list()
			
			popdataOrig = popdata
			if(!missing(casedata)) {
				casedataOrig = casedata
			}
			
			for(Dyear in names(popdataOrig)) {
				
				popdata = popdataOrig[[Dyear]]
				
				if(!missing(casedata)) {
					casedata = casedataOrig[casedataOrig[[caseYearVar]]==Dyear,]
					
					if (length(sex) == 1) {      
						casedata = casedata[casedata[[caseSexVar]] %in% sex, ]
					}
					if(dim(casedata)[1]==0) casedata<-NULL   
				} else { # casedata is missing
					caseThisYear = NULL
				}
				
				# Compute area if spdf not in long-lat
				if(length(grep("Spatial*DataFrame", class(popdata)))) {
					if(!length(grep("longlat", popdata@proj4string@projargs))) {
						areas <- unlist(computeArea(popdata))
						popdata$sqk <- areas * area.scale
					}
				}
				
				# scaling factor to convert to person years
				popScale = interval[names(interval)==Dyear]
					
				popdata[[yearVar]] = as.integer(Dyear)
				
				
				result[[Dyear]] = methods::callGeneric()
				
			}
			names(result) = years
			
			return(result)       
			
		}
)

