precToSd = function(densmat) {
  densmat[,"y"] = densmat[,"y"] * 2*  densmat[,"x"]^(3/2) 
  densmat[,"x"] = 1/sqrt(densmat[,"x"])  
  densmat
}



getStRes = function(x, Squant = c(0.5, 0.025, 0.975)) {

qCols = paste0(Squant, 'quant')

  result = list(marginals = list())

  result$summary = x$summary.hyperpar[
    grep("Precision", rownames(x$summary.hyperpar), invert=TRUE), ]

for(Dpar in names(x$summary.random)) {
  result$marginals[[paste('SD for', Dpar)]] = 
    precToSd(x$marginals.hyperpar[[paste('Precision for', Dpar)]])

  result$summary[paste('SD for', Dpar), ] = NA
  result$summary[paste('SD for', Dpar), qCols] = 
  INLA::inla.qmarginal(
    Squant, 
    result$marginals[[paste('SD for', Dpar)]] )

  result$summary[paste('SD for', Dpar), 'mean'] = 
  INLA::inla.emarginal(
    c,
    result$marginals[[paste('SD for', Dpar)]] )

  result$summary[paste('SD for', Dpar), 'mode'] = 
    result$marginals[[paste('SD for', Dpar)]][
    which.max(result$marginals[[paste('SD for', Dpar)]][,2]),1]
}
  rownames(result$summary) = gsub("GroupRho", 'AR coef',rownames(result$summary) )
  rownames(result$summary) = gsub("Phi ", 'Spatial frac ',rownames(result$summary) )



stEffect = grep("Group", rownames(x$summary.hyperpar), value=TRUE)
stEffect = gsub("GroupRho for ", "", stEffect)

spatialEffect = grep("Phi for", rownames(x$summary.hyperpar), value=TRUE)
spatialEffect = gsub("Phi for ", "", spatialEffect)
spatialEffect = setdiff(spatialEffect, stEffect)


theCol = '0.5quant'
spatialMedian = x$summary.random[[spatialEffect]][,theCol]
spatialMedian = spatialMedian[seq(1, length(spatialMedian)/2)]

stMedian = x$summary.random[[stEffect]][,Squant]
stMedian = array(stMedian, c(
  length(spatialMedian),length(Squant),
  length(stMedian) / (length(Squant)*length(spatialMedian)) ))
result$spatialMedian = spatialMedian
result$stMedian = stMedian

    result
}

stbym = function(
    formula, data, adjMat, region.id, time.id, prior, ...
    ) {

  if(length(region.id) == 1) region.id = c(data=region.id, spatial=region.id)
  if(! all(c('data','spatial') %in% names(region.id))) {
    names(region.id)[1:2] = c('data','spatial')
  }
  region.id = region.id[c('data','spatial')]

  if(class(adjMat) == 'SpatialPolygons') {
    theRegionId = data.frame(xx = 1:length(adjMat))
    names(theRegionId) = region.id[2]
    rownames(theRegionId) = names(adjMat)
    adjMat = sp::SpatialPolygonsDataFrame(adjMat, theRegionId)
  }
  if(class(adjMat) == 'SpatialPolygonsDataFrame') {
      adjMat = spdep::poly2nb(adjMat, row.names=adjMat@data[[ region.id[2] ]])
  }

  if(missing(prior)) prior = list()

  if(! 'spaceTime' %in% names(prior)) {
    prior$spaceTime = list(
      ar = list(prior='pccor0', param=c(0.1, 0.5)),
      sd = c(u=1, alpha=0.5),
      propSpatial = c(u=0.5, alpha=0.5)
      )
  }

  graphFileST =tempfile()
  regionST = diseasemapping::nbToInlaGraph(adjMat, graphFileST)

  data$regionST = regionST[as.character(data[[region.id['data']]])]

  timeStFac = factor(data[[ time.id ]])
  data$timeST = as.integer(timeStFac)


  stFormula = paste0(
      '.~.+f(regionST, ',
      'model=\'bym2\', graph=\'',
      graphFileST, '\', ',
      'hyper = list(',
        'theta1 = list(prior=\'pc.prec\', param=c(',
        prior$spaceTime$sd[1],',', prior$spaceTime$sd[2],')),',
        'theta2 = list(prior = \'pc\', param = c(',
        prior$spaceTime$propSpatial[1],',', prior$spaceTime$propSpatial[2],
        '))),',
      'group=timeST, control.group = list(model=\'ar1\',',
      'hyper =  list(theta=list(prior=\'', prior$spaceTime$ar$prior, '\',',
      'param=c(',
      prior$spaceTime$ar$param[1],',', prior$spaceTime$ar$param[2],
      ')))) )')
  formulaOrig = formula
  formulaUpdated = stats::update.formula(formula, stFormula)

  fromBym = bym(
    formula.fitted=formulaOrig,
    formula = formulaUpdated,
    data = data,
    adjMat = adjMat,
    region.id = region.id['data'],
    prior = prior[setdiff(names(prior), 'spaceTime')],
    ...
    )

  fitStFullRes = getStRes(fromBym$inla)
  fromBym$parameters$sdSpaceTime = list(
    posterior = fitStFullRes$marginals$'SD for regionST',
    summarySpaceTime = fitStFullRes$summary)

stMedianMat = fromBym$inla$summary.random$regionST
stMedian = array(unlist(stMedianMat), 
  c(
  length(adjMat), 
  2, 
  nlevels(timeStFac),
  ncol(stMedianMat)) )

dimnames(stMedian) = list(
  region = attributes(adjMat)$region.id,
  effect = c('bym','indep'),
  time = levels(timeStFac),
  quantile = colnames(stMedianMat)
  )
stMedian = stMedian[,'bym',,setdiff(dimnames(stMedian)[[4]], 'ID')]
  
fromBym$inla$summary.random$regionST = stMedian

fromBym
}
