havePackages = c(
  'INLA' = requireNamespace('INLA', quietly=TRUE)
)

print(havePackages)

if(havePackages) {
	INLA::inla.setOption(num.threads=2)
	# not all versions of INLA support blas.num.threads
	try(INLA::inla.setOption(blas.num.threads=2), silent=TRUE)
	print(INLA::inla.getOption('num.threads'))
}


library('diseasemapping')
data('kentucky')
kentucky = terra::unwrap(kentucky)





kentucky = getSMR(kentucky, larynxRates, larynx,
		regionCode="County")

kentuckyAdjMat = terra::adjacent(kentucky, type='intersects')
attributes(kentuckyAdjMat)$region.id = kentucky$County

if(all(havePackages)){

  kBYM = bym(
			formula = observed ~ offset(logExpected) + poverty,
      data=terra::values(kentucky),
			adjMat = kentuckyAdjMat,
      prior = list(sd=c(0.1, 0.5), propSpatial=c(0.5, 0.5)),
			region.id='County',
			control.predictor=list(compute=TRUE)
  )
	kBYM$parameters$summary

	if(!interactive()) pdf("priorPostKentucky.pdf")
	plot(kBYM$parameters$sd$posterior, type='l', xlim = c(0, 1), lwd=2)
	lines(kBYM$parameters$sd$prior, col='blue')
	legend('topright', lty=1, col=c('black','blue'), legend=c('posterior','prior'))
	if(!interactive()) dev.off()

	
	


# also try no covariate or adj mat

kBYM = bym(
		formula = observed ~ offset(logExpected),
		prior = list(sd = 0.5, propSpatial = 0.5),
		data=kentucky)



 	kBYM$data$exc1 =  unlist(lapply(kBYM$inla$marginals.fitted.bym, INLA::inla.pmarginal, q=log(1.2)))

kBYM$par$summary

if(require('mapmisc', quietly=TRUE)) {

colFit = colourScale(kBYM$data$fitted.exp,
		breaks=6, dec=3)
	
plot(kBYM$data, col=colFit$plot)
legendBreaks('topleft', colFit)

colExc = colourScale(kBYM$data$exc1 ,
		style='fixed',
		breaks=c(0, 0.2, 0.8,0.9, 1), 
		col=rainbow, rev=TRUE
	)

	plot(kBYM$data, col=colExc$plot)
	legendBreaks('topleft', colExc)
 		
}

	
# subtract a few regions

kBYM = bym(
    formula=observed ~ offset(logExpected) + poverty,
    data=terra::values(kentucky)[-(1:4),],  
 	  adjMat = kentuckyAdjMat, region.id="County",
		prior = list(sd=0.1, propSpatial = 0.1))
 

kBYM$par$summary

# intercept only, no offset


kBYM = bym(data=kentucky,  formula=observed ~ 1,
		prior = list(sd=1, propSpatial=0.9))

kBYM$par$summary


if(require('mapmisc', quietly=TRUE)) {
	
	colFit = colourScale(kBYM$data$fitted.exp,
			breaks=6, dec=1)
	
	plot(kBYM$data, col=colFit$plot)
	legendBreaks('topleft', colFit)
	
}

 

# give spdf but some regions have no data
# but keep the 'county' column as is
kentucky[1:2,-grep("County", names(kentucky))] = NA 

kBYM = bym(observed ~ offset(logExpected) + poverty,
		kentucky, 
		region.id="County",
		prior = list(sd=0.1, propSpatial = 0.5))

 
kBYM$par$summary


# missing value in a categorical variable

pCuts = quantile(kentucky$poverty, na.rm=TRUE)
kentucky$povertyFac = cut(kentucky$poverty, 
		breaks = pCuts,
		labels = letters[seq(1,length(pCuts)-1)])
kentucky$povertyFac[c(2,34,100)] = NA

kBYM = bym(
		formula = observed ~ offset(logExpected) + povertyFac,
		data = kentucky, 
		region.id="County",
		prior = list(sd=c(0.1, 0.5), propSpatial=c(0.1, 0.5))
)


kBYM$par$summary
}


