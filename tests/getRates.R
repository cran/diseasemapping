library('diseasemapping')
data('casedata')
data('popdata')
popdata = terra::unwrap(popdata)

therates = getRates(casedata, popdata, ~age*sex,breaks=c(seq(0, 80, by=10), Inf))
therates
