#return a list of areas
`computeArea`<-function(sp){
  allareas=lapply(sp@polygons,function(x) x@area )
  lapply(allareas,sum)
}