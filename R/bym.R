lcOneRow = function(thisrow, idxCol=NULL) {
	thisrow = thisrow[!is.na(thisrow)]
	if(length(thisrow)) {
		thisrow = sapply(thisrow, function(qq) list(list(weight=qq)))
		for(D  in idxCol)
			thisrow[[D]] = list(
					weight=1, 
					idx=thisrow[[D]]$weight
			)
		for(D in names(thisrow))
			thisrow[[D]] = thisrow[D]
		names(thisrow) = paste("v", 1:length(thisrow), sep="")
	}
	thisrow
}

setClass('nb', 
		representation(
		)
)

#' @export
`nbToInlaGraph` = function(adjMat, graphFile="graph.dat", region.id = attributes(adjMat)$region.id)
{

	if(!length(region.id)) {
		region.id = seq(1, max(adjMat[,1]))
	}
	region.index = 1:length(region.id)

	missingRegions = setdiff(region.index, adjMat[,'from'])

	adjMat = rbind(adjMat, cbind(from=missingRegions, to=rep(0, length(missingRegions))))
	adjMat = adjMat[order(adjMat[,'from']),]
	nbList = split(adjMat[,'to'], adjMat[,'from'])

	# get rid of zeros, they mean a region has no neighbours
	nbList = lapply(nbList, function(x) x[x>0])
	
	nbLength = unlist(lapply(nbList, length))
	inlaVec = unlist(lapply(nbList, paste, collapse=' '))
	inlaGraph = c(length(nbList), paste(1:length(nbList), nbLength, inlaVec))
	inlaGraph = paste(inlaGraph, collapse='\n')

	cat(paste(inlaGraph, '\n',sep=''), file=graphFile)

	names(region.index) = as.character(region.id)
	
	attributes(region.index)$Nneighbours = nbLength

	invisible(region.index)
}

#' @importClassesFrom sp SpatialPolygonsDataFrame
#' @export
setGeneric('bym', 
		function(
				formula, data, adjMat=NULL, region.id,
				...) {
			standardGeneric("bym")
		}
)

bym.spartan = function(
		formula, data, adjMat=NULL, region.id,
		...) {	
	
	region.idX="region.id"
	data[[region.idX]] = 1:length(data)
	
	methods::callGeneric(
			formula , data ,
			adjMat , region.idX,
			...
	)
}


setMethod("bym", 
		signature("formula","ANY","ANY", "missing"),
	bym.spartan
		)

	setMethod("bym", 
				signature("formula","ANY","missing", "missing"),
				bym.spartan
	)
		

bym.needAdjmat = function(
			formula, data, adjMat=NULL, region.id,
		...) {	
 	
		if(missing(region.id)) {
			region.id = 'region.id'
			terra::values(data)[,region.id]=1:length(data)
		}
	# intersects works better than queen, queen misses some neighbours
 	adjMatNB = terra::adjacent(data, type='intersects')
 	attributes(adjMatNB)$region.id = unname(unlist(data[[region.id]]))

	methods::callGeneric(
			formula=formula, data=data,
			adjMat=adjMatNB, region.id=region.id,
			...
	)
	
	}	
setMethod("bym", 
			signature("formula", "SpatVector", "missing","character"),
			bym.needAdjmat		
)
	
setMethod("bym", 
			signature("formula", "SpatVector", "NULL","character"),
			bym.needAdjmat		
)

setMethod("bym", 
		signature("formula", "SpatVector", "matrix","character"),
		function(
				formula, data, adjMat,region.id,
				...) {	
			
			
			result = methods::callGeneric(
					formula=formula, data=terra::values(data),
					adjMat=adjMat,region.id=region.id,
					...
			)


			if(any(names(result)=="data")) {		
				# merge data back into SPDF
				resultSV <- terra::merge(data[,region.id], result$data)
				result$data = resultSV
			}
			
			result
			
		}	
)
			

bym.data.frame = function(formula, data, adjMat,		region.id, 
		prior=list(sd = c(0.1, 0.5), propSpatial = c(0.5, 0.5)), 
		family="poisson",
		formula.fitted=formula,
		...
	) {
		#neighbourhood structure
		graphfile=tempfile()
		
		# if using windows, replace back slashes with forward slashes...
		graphfile = gsub("\\\\", "/", graphfile)
		

		if(missing(region.id)  & is.null(attributes(adjMat)$region.id) ) {
			if(nrow(data) == max(adjMat[,1])) {
			attributes(adjMat)$region.id = data[[region.id]]
				} else {
				warning("adjMat should have a region.id attribute")
			}
		}

		region.index = diseasemapping::nbToInlaGraph(adjMat, graphfile)
		
		# check for regions without neighbours
		badNeighbours = which(
				attributes(region.index)$Nneighbours < 1
		)
		if(length(badNeighbours)){
			if(length(badNeighbours) == length(data))
				stop('No spatial regions are neighbours of each other.')
			warning('There are ', length(badNeighbours), ' regions without neighbours, consider removing these.')
		}
		
		
		# check for data regions missing from adj mat
		data[,region.id] = as.character(data[,region.id])
		if(!all(data[[region.id]] %in% names(region.index))  )
			warning("regions in data missing from adjacency matrix")
		data$region.indexS = data$region.indexI = region.index[data[[region.id]]]
		
		
	
		precPrior = list()
		
		# priors
		if(all(c("sd","propSpatial")%in% names(prior))) {
			
			# pc priors
			# if length zero, set to default
			if(!length(prior$sd)){
				prior$sd = 1
			}	
			if(!length(prior$propSpatial)){
				prior$propSpatial = c(0.5)
			}	
			
		
			# if length 1, assume u provided and alpha = 0.95
			for(Dpar in c('sd','propSpatial')){
				if(length(prior[[Dpar]])==1) 
					prior[[Dpar]] = c(prior[[Dpar]], 0.05)
				if(is.null(names(prior[[Dpar]]))) 
					names(prior[[Dpar]])[1:2] = c('u','alpha')
			}	
	

			
		
			bymTerm = paste(
					".~.+f(region.indexS, model='bym2', graph='",
					graphfile,
					"', hyper = list(theta1=list(prior='pc.prec', param=c(",
					paste(prior[["sd"]], collapse=","), 
					")), theta2=list(prior = 'pc', param=c(", 		
					paste(prior[["propSpatial"]], collapse=","),
					"))))",
					sep="")
			
		} else if(all(c("sdSpatial","sdIndep")%in% names(prior))) {
		
		for(D in c("sdSpatial","sdIndep")) {
			obj1 = sort(prior[[D]]^-2)
		
			startMean = mean(obj1)
			startSD = diff(obj1)/4
			startShape = startSD^2/startMean^2
			
	cifun = function(pars, obj) {
			theci = 	stats::pgamma(obj, shape=pars[1], rate=pars[2],log.p=T)
			
			(log(0.025) - theci[1])^2 +
					(2*(log(0.975) - theci[2]))^2		
			
		}


			precPrior2=stats::optim(c(startShape,startShape/startMean), cifun, 
					lower=c(0.000001,0.0000001),method="L-BFGS-B", obj=obj1,
					control=list(parscale=c(startShape,startShape/startMean)))
			precPrior[[D]] = precPrior2$par
			names(precPrior[[D]] ) = c("shape","rate")
			
			stats::pgamma(obj1, shape= precPrior[[D]]["shape"], rate=precPrior[[D]]["rate"],log.p=F)
			stats::pgamma(obj1, shape= precPrior[[D]]["shape"], rate=precPrior[[D]]["rate"],log.p=T)
			log(c(0.025, 0.975))
			precPrior2
			stats::pgamma(obj1, shape=precPrior[[D]]["shape"], rate=precPrior[[D]]["rate"],log.p=T)
			log(c(0.025, 0.975)) 
			1/sqrt(stats::qgamma(c(0.975,0.025), shape=precPrior[[D]]["shape"], rate=precPrior[[D]]["rate"]))
			prior[[D]]
			
		} 
		
		#f(CSDUID, model = "bym", graph.file = "nb.graph", 
		#		+     param = c(prior.iid, prior.besag), values = CSDUID)
		
		
	bymTerm = paste(
			".~.+f(region.indexS, model='bym', graph='",
			graphfile,
			"', hyper = list(theta2=list(param=c(",
			paste(precPrior[["sdSpatial"]], collapse=","), 
			")), theta1=list(param=c(", 		
			paste(precPrior[["sdIndep"]], collapse=","),
			"))))",
			sep="")
	} else {
		warning("prior needs elements sdSpatial and sdIndep, or sd and propSpatial")
	}	
	
	allVars = all.vars(formula)
	formula = stats::update(formula, stats::as.formula(bymTerm))

	if(!all(allVars %in% names(data))) {
		warning(paste("missing variables", toString(setdiff(allVars, names(data)))))
	}

	# INLA doesn't like missing values in categorical covariates
	# remove rows with NA's 
	if(is.matrix(data[[allVars[1]]])) {
		# response is a matrix, don't look for NA's
		anyNA = which(apply(data[,allVars[-1], drop=FALSE], 1, function(qq) any(is.na(qq))))
	} else {
		anyNA = which(apply(data[,allVars, drop=FALSE], 1, function(qq) any(is.na(qq))))
	}
	if(length(anyNA)) {
		data = data[-anyNA, ]
	}

	##################
	# linear combinations
	###################
		# fitted values

# get rid of left side of formula
#  formulaForLincombs =  formulaRhs(formula.fitted,char=TRUE)
formulaForLincombs = base::format(formula.fitted)
# if there is a line break in the formula, 
# format(formula) will create a vector
formulaForLincombs = paste(formulaForLincombs, collapse="")
formulaForLincombs = toString(formulaForLincombs)
#formulaForLincombs = gsub("^.*~", "", toString(formulaForLincombs))
theTilde = gregexpr("[~]",  formulaForLincombs)[[1]][1]
formulaForLincombs = substr(formulaForLincombs, theTilde+1, nchar(formulaForLincombs))

# get rid of f(stuff) in formula
formulaForLincombs =
		gsub("f\\([[:print:]]*\\)", "", formulaForLincombs)
formulaForLincombs = gsub(
		"\\+[[:space:]]?\\+([[:space:]]?\\+)?", "+",
		formulaForLincombs)
# get rid of offset(stuff)
formulaForLincombs =
		gsub("offset\\([[:print:]]*\\)[[:space:]]?($|\\+)", "", formulaForLincombs)

# convert multiple + to a single +
formulaForLincombs = gsub(
		"\\+[[:space:]]?\\+([[:space:]]?\\+)?", "+",
		formulaForLincombs)
# strip out trailing or leading +
formulaForLincombs = gsub("\\+[[:space:]]+?$|^[[:space:]]?\\+[[:space:]]+", "", formulaForLincombs)


	startIndex = length(region.index)
 
	# if there are covariates
	if(nchar(trimws(formulaForLincombs)) & formulaForLincombs != "1" &
			!length(grep("^[[:space:]]+$", formulaForLincombs))
		) { #make linear combinations
 
		formulaForLincombs=stats::as.formula(
			paste("~", paste(c("1",formulaForLincombs),collapse="+"))
		)
 
		# remove regions in the data set twice
		theDuplicated = duplicated(data$region.indexS)
		notDuplicated = which(!theDuplicated)
		
		# reorder the matrix by region ID
		dataOrder = data[notDuplicated,]
		dataOrder = dataOrder[order(dataOrder$region.indexS),]

		
		lincombFrame = stats::model.frame(formulaForLincombs, dataOrder,
				na.action=stats::na.omit)

		
		SregionFitted = dataOrder[rownames(lincombFrame),"region.indexS"]
		names(SregionFitted) = dataOrder[rownames(lincombFrame),region.id]
		
		
		lincombMat = stats::model.matrix(formulaForLincombs, lincombFrame)
		
		lincombMat[lincombMat==0]= NA
		lincombMat = cbind(lincombMat, region.indexS = dataOrder[rownames(lincombMat),"region.indexS"])
		
		if(!dim(lincombMat)[1])
			warning("the dataset appears to have no rows")

		lcFitted <- apply(lincombMat, 1, lcOneRow, idxCol="region.indexS")
		names(lcFitted) = paste("fitted_", dataOrder[rownames(lincombMat), region.id],sep="")

		inlaLincombs = lcFitted
		

	} else { # add only intercept to predictions because no inla or no covariates
		formulaForLincombs = ~1
		lincombMat = data.frame(x=rep(1,length(region.index)))
		SregionFitted = region.index
		inlaLincombs = list()
		for(D in 1:length(region.index)) {	
			inlaLincombs[[D]] = 
					 list(
							list("(Intercept)" = list(weight=1)),
							list(region.indexS=
											list(idx=region.index[D], weight=1))
					)
		}
		names(inlaLincombs) =
				paste("fitted",names(region.index),sep="_")
	}
 


	##########################
	# run inla!		
	########################
	if(requireNamespace("INLA", quietly=TRUE)) { # not enough to have requireNamespace
 
			
		inlaRes = INLA::inla(formula, data=data, 
				family=family,
			lincomb=inlaLincombs, ...)
	} else{
		inlaRes = 
				list(logfile="INLA is not installed.  \n see www.r-inla.org")
	}
 	
	if(all(names(inlaRes)=="logfile"))
		return(c(list(formula=formula, data=data,
						family=family, graphfile=graphfile,
						lincomb=inlaLincombs, 
						ldots = list(...)),
						inlaRes)
	)
	
	inlaRes$graphfile = graphfile

	# posterior distributions of random effect (spatial + independent)
	Sbym = seq(1, length(region.index))
	thebym = inlaRes$summary.random$region.indexS[Sbym,]
	
	inlaRes$marginals.bym = inlaRes$marginals.random$region.indexS[
			Sbym
			]
	

		# E(exp(U)  | Y)
		meanExp = unlist(
				lapply(inlaRes$marginals.bym, 
						function(qq) {
							sum(
									exp(qq[,"x"])*c(0,diff(qq[,"x"]))*qq[,"y"]	
							)
						})
		) # end unlist
		meanExp[meanExp==Inf]=NA
		thebym = cbind(thebym, exp=meanExp)
			
	thebym = thebym[,!names(thebym) %in% c("ID","kld")]
	colnames(thebym) = paste("random.",colnames(thebym),sep="")
	names(inlaRes$marginals.bym) = rownames(thebym) = 
			names(region.index)
	

	# make sure they're in the correct order
	thebym = thebym[names(region.index),]
	inlaRes$marginals.bym = inlaRes$marginals.bym[names(region.index)] 

	# fitted values, some regions dont have them if covariates are missing
	theFitted = inlaRes$summary.lincomb.derived[
			grep("^fitted_", rownames(inlaRes$summary.lincomb.derived)),]
	inlaRes$marginals.fitted.bym = inlaRes$marginals.lincomb.derived[
			grep("^fitted_", names(inlaRes$marginals.lincomb.derived))
			]
			
	names(inlaRes$marginals.fitted.bym)= 
			gsub("^fitted_", "", names(inlaRes$marginals.fitted.bym))
	rownames(theFitted) = gsub("^fitted_", "", rownames(theFitted))
	
	# E(exp(lambda)  | Y)
	meanExp = unlist(
					lapply(inlaRes$marginals.fitted.bym, 
							function(qq) {
								sum(
										exp(qq[,"x"])*c(0,diff(qq[,"x"]))*qq[,"y"]	
								)
							})
				) # end unlist
	meanExp[meanExp==Inf]=NA
	theFitted = cbind(theFitted, exp = meanExp[rownames(theFitted)])
	
	# E inv logit(lincombs)

	if(length(grep("^binomial$",inlaRes$.args$family))) {
		invlogit=unlist(
				lapply(inlaRes$marginals.fitted.bym, 
						function(qq) {
							eqqx = exp(qq[,"x"])
							sum(
									eqqx/(1+eqqx)*c(0,diff(qq[,"x"]))*qq[,"y"]	
							)
						})
		)
		theFitted = cbind(theFitted, invlogit = invlogit[rownames(theFitted)])
	}

	
	theFitted = theFitted[,!names(theFitted) %in% c("ID","kld")]
	colnames(theFitted) = paste("fitted.",colnames(theFitted),sep="")

	# merge fitted falue summary into BYM
	thebym = cbind(thebym, matrix(NA, dim(thebym)[1], dim(theFitted)[2],
					dimnames=list(NULL, names(theFitted))))
	thebym[rownames(theFitted), colnames(theFitted)] = theFitted
	thebym[,region.id] = names(region.index)
	
	
	# make fitted marginals list have same names and order as fitted radom
	notInFitted = names(inlaRes$marginals.bym) [
		!names(inlaRes$marginals.bym)%in%names(inlaRes$marginals.fitted.bym)]

	toAdd = replicate(length(notInFitted),NULL,simplify=FALSE)
	names(toAdd) = notInFitted
	inlaRes$marginals.fitted.bym = c(inlaRes$marginals.fitted.bym, 
			 toAdd)
	 inlaRes$marginals.fitted.bym = inlaRes$marginals.fitted.bym[
			 names(inlaRes$marginals.bym) 
			 ]
	
	# the parameters
	params=list()
	params$summary = inlaRes$summary.fixed
	
	params$summary = cbind(params$summary, 
			meanExp = unlist(
					lapply(inlaRes$marginals.fixed,
							function(qq) {
								sum(
										exp(qq[,"x"])*c(0,diff(qq[,"x"]))*qq[,"y"]	
								)
							}
					)
	))
	
	
	
	quantNames = grep("quant$", colnames(params$summary), value=TRUE)
	revQuant = rev(quantNames)	
	
	sdNames = paste(c(sd="Precision",propSpatial="Phi"), "for region.indexS")
	names(sdNames)= c('sd','propSpatial')
	
	if(all( names(sdNames) %in% names(prior))){
		# model='bym2'
	 # sd
		params$sd = list(
				params.intern = prior$sd,
				prior = 1/sqrt(
					INLA::inla.pc.qprec(c(0.975,0.025),  
						u = prior$sd['u'], alpha = prior$sd['alpha'])
				)
		)
			
		imname = grep(
				"^Precision.*region.indexS",
				rownames(inlaRes$summary.hyperpar),
				value=TRUE
		)
		Dname = 'sd'
		params[[Dname]]$posterior=
				inlaRes$marginals.hyperpar[[
						imname	
				]]
		params[[Dname]]$posterior[,"y"] = params[[Dname]]$posterior[,"y"] * 2*  
				params[[Dname]]$posterior[,"x"]^(3/2) 
		params[[Dname]]$posterior[,"x"] = 1/sqrt(params[[Dname]]$posterior[,"x"])  
		params[[Dname]]$posterior = 
				params[[Dname]]$posterior[seq(dim(params[[Dname]]$posterior)[1],1),]		
		
		
		precLim = range( inlaRes$marginals.hyperpar[[
						imname
				]][,1] ) 
		sdSeq = seq(0, 2*max(params[[Dname]]$prior), len=1000)
		precSeq = sdSeq^(-2)
				
		params[[Dname]]$prior=cbind(
				x=sdSeq,
				y=INLA::inla.pc.dprec(precSeq, 
						u = prior$sd['u'], alpha = prior$sd['alpha']
				) * 2 * (precSeq)^(3/2) 
		)

		thesummary = inlaRes$summary.hyperpar[imname, ,drop=FALSE]
		
		thesummary[,quantNames] = 1/sqrt(thesummary[,revQuant])
		
		thesummary[,"mean"] =sum(
				1/sqrt(inlaRes$marginals.hyperpar[[imname]][,"x"])*
						c(0,diff(inlaRes$marginals.hyperpar[[imname]][,"x"]))*
						inlaRes$marginals.hyperpar[[imname]][,"y"]
		)
		thesummary[,"sd"] = NA
    thesummary[,'mode'] =  1/sqrt(thesummary[,'mode'])
		rownames(thesummary) = Dname
		
		donthave = colnames(params$summary)[
				!colnames(params$summary) %in% colnames(thesummary)]
		
		thesummary = cbind(thesummary, matrix(NA, nrow(thesummary), length(donthave),
						dimnames=list(NULL, donthave)))
  	
		params$summary = rbind(params$summary, 
				thesummary[,colnames(params$summary),drop=FALSE]
		)		
		
		
		# propSpatial
		
		Deffect = which(unlist(lapply(inlaRes$all.hyper$random, function(qq) qq$hyperid))==
						'region.indexS')

		priorProp = matrix(scan(text=gsub("[[:alpha:]]+:", "", 
								inlaRes$all.hyper$random[[Deffect]]$hyper$theta2$prior),
						quiet=TRUE), 
				ncol=2, dimnames = list(NULL, c( "xTrans","logDensTrans")))
		priorProp = cbind(priorProp, orX = exp(priorProp[,'xTrans']))
		priorProp = cbind(priorProp, 
				x = priorProp[,'orX']/ (1+priorProp[,'orX'])
		)
		
		# what the prior integrates to	
		const = sum(exp(priorProp[,'logDensTrans']+
								c(NA,log(abs(diff(priorProp[,'xTrans']))))),
				na.rm=TRUE)		
	# add to it integrates to 1		
		
	diffTrans = c(NA, diff(priorProp[,'xTrans']))
	diffX = c(NA, diff(priorProp[,'x']))
	
 # jacobian is p*(1-p)
		priorProp = cbind(priorProp, 
				logDens = priorProp[,'logDensTrans'] - log(priorProp[,'x']) - log(1-priorProp[,'x'])
#						log(abs(diffTrans)) - 
#						log(abs(diffX)) - log(const)
		)
		priorProp[1, 'logDens'] = priorProp[2, 'logDens']

		priorProp = cbind(priorProp, 
				dens = exp(priorProp[,'logDens'] )
				)
		priorProp = cbind(priorProp, 
						cDens = cumsum(priorProp[,'dens'] *
  									c(0, abs(diff(priorProp[,'x']))))
		)
		
		imname = grep(
				"^Phi.*region.indexS",
				rownames(inlaRes$summary.hyperpar),
				value=TRUE
		)
		phiSeq = seq(0,1,len=1000)
		
		params$propSpatial = list(
				params.intern = prior$propSpatial,
				priorCI = stats::approx(priorProp[,'cDens'], 
						priorProp[,'x'], c(0.025, 0.975))$y,
				posterior = inlaRes$marginals.hyperpar[[imname]],
				prior = as.data.frame(stats::approx(
						priorProp[,'x'], priorProp[,'dens'], phiSeq
						))
		)
		
		thesummary = inlaRes$summary.hyperpar[imname, ,drop=FALSE]
		donthave = colnames(params$summary)[
				!colnames(params$summary) %in% colnames(thesummary)]
		
		thesummary = cbind(thesummary, matrix(NA, nrow(thesummary), length(donthave),
						dimnames=list(NULL, donthave)))
		params$summary = rbind(params$summary, 
				thesummary[,colnames(params$summary),drop=FALSE]
		)		

		
		
	} else {
		# model = 'bym'
	sdNames = c(S='spatial', I='iid')
	for(D in c("S","I")) {
		
		Dname = grep( paste("^sd",D,sep=""),names(prior),value=TRUE)
		
		params[[Dname]] = list(
				userprior=prior[[Dname]], 
				prior = 
						1/sqrt(
								stats::qgamma(c(0.975,0.025), 
										shape=precPrior[[Dname]]["shape"], 
										rate=precPrior[[Dname]]["rate"])),
				params.intern=precPrior[[Dname]])
	
		
		imname = grep(
				paste("Precision.*region.indexS.*", sdNames[D], ' component',sep=""),
				rownames(inlaRes$summary.hyperpar),
				value=TRUE
		)
		
		params[[Dname]]$posterior=
				inlaRes$marginals.hyperpar[[
					imname	
				]]
		params[[Dname]]$posterior[,"y"] = params[[Dname]]$posterior[,"y"] * 2*  
				params[[Dname]]$posterior[,"x"]^(3/2) 
		params[[Dname]]$posterior[,"x"] = 1/sqrt(params[[Dname]]$posterior[,"x"])  
		params[[Dname]]$posterior = 
				params[[Dname]]$posterior[seq(dim(params[[Dname]]$posterior)[1],1),]		
		
		
		precLim = range( inlaRes$marginals.hyperpar[[
								imname
						]][,1] ) 
		precLim = precLim * c(0.8, 1.2)		
		sdLim = 1/sqrt(precLim)
		sdSeq = seq(min(sdLim), max(sdLim), len=1000)
		precSeq = sdSeq^(-2)
		params[[Dname]]$prior=cbind(
				x=sdSeq,
				y=stats::dgamma(precSeq, shape=precPrior[[Dname]]["shape"], 
						rate=precPrior[[Dname]]["rate"]) *2* (precSeq)^(3/2) 
		)
		
		thesummary = inlaRes$summary.hyperpar[imname, ,drop=FALSE]
		thesummary[,quantNames] = 1/sqrt(thesummary[,revQuant])
		
		thesummary[,"mean"] =sum(
				1/sqrt(inlaRes$marginals.hyperpar[[imname]][,"x"])*
						c(0,diff(inlaRes$marginals.hyperpar[[imname]][,"x"]))*
						inlaRes$marginals.hyperpar[[imname]][,"y"]
		)
		thesummary[,"sd"] = NA
    thesummary[,'mode'] =  1/sqrt(thesummary[,'mode'])
		rownames(thesummary) = Dname
 
		
		donthave = colnames(params$summary)[
				!colnames(params$summary) %in% colnames(thesummary)]

		thesummary = cbind(thesummary, matrix(NA, nrow(thesummary), length(donthave),
							dimnames=list(NULL, donthave)))
 
			
		params$summary = rbind(params$summary, 
				thesummary[,colnames(params$summary),drop=FALSE]
				)		
	}
	}
	
	rownames(params$summary) = gsub(
			"^Phi[[:space:]]+for[[:space:]]+region.*", "propSpatial", 
			rownames(params$summary)
	)
	


	return(list( inla=inlaRes, data=thebym, parameters=params))
}

setMethod("bym", 
		signature("formula", "data.frame", "matrix","character"),
		bym.data.frame
)
