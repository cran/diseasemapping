#return a list of areas
`area`<-function(sp){
  allareas=lapply(sp@polygons,function(x) x@area )
  lapply(allareas,sum)
}