## ----knitr, include=FALSE-----------------------------------------------------
require('knitr')
opts_chunk$set(out.width='0.48\\textwidth', fig.align='default', fig.height=3, fig.width=5, tidy=FALSE)
knit_hooks$set(source  =function(x, options) {
  paste0(c('\\begin{verbatim}', x, '\\end{verbatim}', ''),
      collapse = '\n')
})
knit_hooks$set(chunk = function(x, options) {x})
hook_output <- function(x, options) {
  if (knitr:::output_asis(x, options)) return(x)
  paste0('\\begin{verbatim}\n', x, '\\end{verbatim}\n')
}
knit_hooks$set(output = hook_output)
knit_hooks$set(message = hook_output)
knit_hooks$set(warning = hook_output)
knit_hooks$set(plot = function (x, options) 
    {
      paste0(knitr::hook_plot_tex(x, options), "\n")
   })

## ----inla, include=FALSE------------------------------------------------------
library('terra')
library('mapmisc')
if(requireNamespace('INLA')) {
	INLA::inla.setOption(num.threads=2)
	# not all versions of INLA support blas.num.threads
	try(INLA::inla.setOption(blas.num.threads=2), silent=TRUE)
}

## ----packages-----------------------------------------------------------------
require('diseasemapping')
data('kentucky')
kentucky = terra::unwrap(kentucky)

## ----rates--------------------------------------------------------------------
if(FALSE) {
	# must have an internet connection to do the following
	larynxRates= cancerRates("USA", year=1998:2002,site="Larynx")
	dput(larynxRates)
} else {
	larynxRates = structure(c(0, 0, 0, 0, 1e-06, 6e-06, 2.3e-05, 4.5e-05, 9.9e-05, 
					0.000163, 0.000243, 0.000299, 0.000343, 0.000308, 0.000291, 0.000217, 
					0, 0, 0, 1e-06, 1e-06, 3e-06, 8e-06, 1.3e-05, 2.3e-05, 3.5e-05, 
					5.8e-05, 6.8e-05, 7.5e-05, 5.5e-05, 4.1e-05, 3e-05), .Names = c("M_10", 
					"M_15", "M_20", "M_25", "M_30", "M_35", "M_40", "M_45", "M_50", 
					"M_55", "M_60", "M_65", "M_70", "M_75", "M_80", "M_85", "F_10", 
					"F_15", "F_20", "F_25", "F_30", "F_35", "F_40", "F_45", "F_50", 
					"F_55", "F_60", "F_65", "F_70", "F_75", "F_80", "F_85"))
	
}

## ----getRid-------------------------------------------------------------------
larynxRates = larynxRates[grep("_(0|5)$",names(larynxRates), invert=TRUE)]

## ----smr----------------------------------------------------------------------
kentucky = diseasemapping::getSMR(
    popdata=kentucky, 
    model = larynxRates, 
    casedata=larynx, 
    regionCode="County")

## ----bymPc--------------------------------------------------------------------
kBYMpc = try(
bym(
	formula = observed ~ offset(logExpected) + poverty,
  kentucky,
	prior = list(
      sd=c(u=1, alpha=0.05), 
      propSpatial = c(u=0.5, alpha=0.8)),
  verbose=TRUE), silent=TRUE)

## ----bymPcTry, include=FALSE--------------------------------------------------
if(class(kBYMpc) == 'try-error') 
	kBYMpc = list()

## ----summaryPc----------------------------------------------------------------
if(!is.null(kBYMpc$parameters))
	knitr::kable(kBYMpc$parameters$summary[,c(1,3,5)], digits=3)

## ----priorPostPc, fig.cap="PC priors variance parameters", fig.height=4, fig.width=3, fig.subcap=c("sd","prop spatial"), echo=FALSE----

if(!is.null(kBYMpc$parameters)) {

par(mar = c(3,3,0,0))	
plot(kBYMpc$parameters$sd$posterior, type='l', 
		xlim=c(0,1))
lines(kBYMpc$parameters$sd$prior, col='blue')
abline(v=kBYMpc$parameters$sd$params.intern[1], lty=3)
legend('topright', lty=1, col=c('black','blue'), legend=c('posterior','prior'), bg='white')

plot(kBYMpc$parameters$propSpatial$posterior, 
  type='l', 
		xlim=c(0, 1))
lines(kBYMpc$parameters$propSpatial$prior, col='blue')
abline(v=kBYMpc$parameters$propSpatial$params.intern[1], lty=3)
legend('topright', lty=1, col=c('black','blue'), legend=c('posterior','prior'))
} else {
	plot(1:10, type='n')
	text(5,5,'inla is not installed')
	plot(1:10, type='n')
	text(5,5,'inla is not installed')
}

## ----maps, fig.cap='Random effects and fitted values', fig.subcap=c('fitted','random'), echo=FALSE, out.width='0.9\\textwidth'----

if(require('mapmisc', quietly=TRUE) & !is.null(kBYMpc$parameters)) {
	
thecex=0.6	
	

colFitPc = colourScale(kBYMpc$data$fitted.exp, breaks=6, dec=1, style='equal')

plot(kBYMpc$data, col=colFitPc$plot)
legendBreaks('left', colFitPc, cex=thecex)


colRpc = colourScale(kBYMpc$data$random.mean,breaks=12, dec=-log10(0.05), style='equal')


plot(kBYMpc$data, col=colRpc$plot)
legendBreaks('left', colRpc, cex=thecex)


} else {
	plot(1:10, type='n')
	text(5,5,'inla is not installed')
	plot(1:10, type='n')
	text(5,5,'inla is not installed')
}

