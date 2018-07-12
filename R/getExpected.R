aggPopLong = function(popdata, poplong, sex, area.scale) {
	
# called by getSMR df numeric and getSMR df glm  
	
	if('sex' %in% names(poplong)){
		sex = paste("^", paste(sex, collapse='$|^'), "$", sep='')
		poplong= poplong[grep(sex, poplong$sex, ignore.case=TRUE),]
	}	
	
	poplong <- stats::aggregate(
			poplong$expected, 
			list(poplong[['idForAgg']]), 
			sum, na.rm=TRUE)
	rownames(poplong) = as.character(poplong[,1])
	poplong=poplong[poplong[,2] > 0,]
	
	
# merge results back in to the population data
# the merge function changes the order, so can't use it.
	popdata$expected = NA
	rownames(popdata) = as.character(popdata[,'idForAgg'])
	
	popdata[rownames(poplong), "expected"] = poplong[,2]
	
	surfaceArea = grep("^surfacearea$", names(popdata), ignore.case=TRUE, value=TRUE)
	if (length(surfaceArea) ) {
		popdata$expected_surfaceArea <- area.scale *popdata$expected/popdata[[surfaceArea[1] ]]
		popdata$logExpected_surfaceArea = log(popdata$expected_surfaceArea)
	}
	
	popdata$logExpected = log(popdata$expected)
	
	popdata = popdata[,
			grep('^idForAgg$', colnames(popdata), invert=TRUE)]
	
	
	popdata
	
}

setGeneric('getExpected', 
		function(popdata, model,
				area.scale=1, 
				sex=c('m','f')) {
			standardGeneric("getExpected")
		}
)

setMethod("getExpected", 
		signature("data.frame", "numeric"),
		function(popdata, model, 
				area.scale=1, 
				sex=c('m','f')){
			# model is a vector of rates
			popdata$idForAgg = 1:nrow(popdata)
			
			# check breaks for groups, make sure they line up
			rateBreaks = getBreaks(names(model))
			popBreaks = getBreaks(names(popdata))
			
			newBreaks = getBreaks(intersect(rateBreaks$newNames, popBreaks$newNames))
			
			
			newModel = data.frame(
					age=rateBreaks$age, 
					sex=rateBreaks$sex, 
					rate=model)
			newModel = formatCases(newModel, newBreaks)
			newModel = stats::aggregate(newModel$rate, 
					newModel[,c('age','sex')],mean,na.rm=T)
			rownames(newModel) = 
					paste(newModel$sex, 
							newModel$age, 
							sep='.')
			
			poplong = formatPopulation(popdata, breaks=newBreaks$breaks)
			
			poplong$expected = poplong$POPULATION *  
					newModel[paste(poplong$sex, poplong$age, sep='.')
							,'x']
			
			aggPopLong(popdata, poplong, sex, area.scale) 
		}
)

setMethod("getExpected", 
		signature("data.frame", "glm"),
		function(popdata, model,
				area.scale=1, 
				sex=c('m','f')){
			
			
			allVars= gsub("offset[(]|[)]$", "", rownames(attributes(model$terms)$factors)[-1])
			if(all(allVars %in% colnames(popdata))) {
				
				# data doesn't need formatting
				theFactors = attributes(model$terms)$dataClasses
				theFactors = names(theFactors[theFactors == 'factor'])
				# get rid of levels 
				for(D in theFactors) {
					popdata = popdata[as.character(popdata[,D]) %in% 
									levels(model$data[,D]), ]
				}
				popdata$expected = predict(model, popdata, type='response')
				if(length(model$aggregate))
					popdata = stats::aggregate(
							popdata[,'expected', drop=FALSE],
							popdata[,model$aggregate, drop=FALSE],
							sum, na.rm=TRUE)
				
				return(popdata)
			}
			
			popdata$idForAgg = 1:nrow(popdata)
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
			agevar<-grep("^age$",names(attributes((stats::terms(model)))$dataClasses),value=TRUE,ignore.case=TRUE)
			sexvar<-grep("^sex$",names(attributes((stats::terms(model)))$dataClasses),value=TRUE,ignore.case=TRUE)
			yearvar<-grep("^year$",names(attributes((stats::terms(model)))$dataClasses),value=TRUE,ignore.case=TRUE)
			
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
			
			
			offsetvar<- grep("^offset",names(attributes((stats::terms(model)))$dataClasses)  ,value=TRUE, ignore.case=TRUE)
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
				poplong[,offsetvar] = poplong[,offsetvar] + 
						log(attributes(popdata)$popScale)
			}
			
			
			poplong$expected <- stats::predict(model, poplong,#[,agg],
					type = "response")
			
			aggPopLong(popdata, poplong, sex, area.scale) 
			
		}
)




setMethod("getExpected", 
		signature("data.frame", 'list'),
		function(popdata, model, 
				area.scale=1, 
				sex=c('m','f')){
			
			keepCols = c(
					'expected','expected_surfaceArea',
					'logExpected_surfaceArea', 'logExpected'
			)
			
			
			if(!length(names(model))){
				names(model) = as.character(1:length(model))
			}
			
			result =  as.data.frame(matrix(NA, nrow(popdata), 0))
			
			
			modelFull = model
			for(Drate in names(modelFull)){
				model=modelFull[[Drate]]
				toBind = methods::callGeneric()
				toBind = toBind[,intersect(colnames(toBind), keepCols)]
				names(toBind) = paste(
						names(toBind), Drate, sep="_"
				)
				result = cbind(result,toBind)
			}
			
			popdata = cbind(popdata, result)
			popdata
		}
)

setMethod("getExpected", 
		signature("list", 'list'),
		function(popdata, model, 
				area.scale=1, 
				sex=c('m','f')){
			
			Nmax = max(c(length(popdata),length(model)))
			
			if(!length(names(model))){
				names(model) = as.character(round(seq(from=1,to=Nmax, len=length(model))))
			}
			if(!length(names(popdata))){
				names(popdata) = as.character(round(seq(from=1,to=Nmax, len=length(popdata))))
			}
			
			Scensus = as.numeric(names(popdata))
			Srates = as.numeric(names(model))
			if(any(is.na(Scensus)) | any(is.na(Srates)))
				warning("can't turn names of popdata or model into numbers")
			
			yearBreaks = c(
					-Inf,
					Scensus[-1] - diff(Scensus)/2,
					Inf
			)
			
			rateCut = Scensus[
					as.integer(cut(Srates, breaks=yearBreaks))
			]
			
			popDataFull = popdata
			modelFull = model
			
			
			result = list()
			
			
			for(Dcensus in Scensus){
				popdata = popDataFull[[as.character(Dcensus)]]
				model = modelFull[
						rateCut==Dcensus
				]
				result[[
						as.character(Dcensus)
				]] = methods::callGeneric()
			}
			result
		}
)
