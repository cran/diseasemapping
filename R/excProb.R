excProb = function(...){
	
	if(requireNamespace('geostatsp', quietly=TRUE)){
		result = geostatsp::excProb(...)
	}  else {
		result = NA
	}
	result
	
}