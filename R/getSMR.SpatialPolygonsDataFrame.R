`getSMR.SpatialPolygonsDataFrame` <- function(popdata, ...) {

popdata@data <- getSMR(popdata@data, ...)
popdata
}
