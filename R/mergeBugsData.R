mergeBugsData = function(x, bugsSummary, by.x = NULL, newcol="mean", ...) {
  UseMethod("mergeBugsData")
}

mergeBugsData.SpatialPolygonsDataFrame =
  function(x, bugsSummary, by.x=NULL, newcol="mean", ...) {
    x@data = mergeBugsData(x@data, bugsSummary, by.x, newcol, ...)
    x
}

mergeBugsData.data.frame = function(x, bugsSummary,
   by.x=NULL, newcol="mean", ...) {

  if(!is.list(bugsSummary)) {
    if(is.vector(bugsSummary)  ) {
      bugsSummary = data.frame(bugsSummary)
      names(bugsSummary) = newcol

    }
    bugsSummary = list( R= bugsSummary)

    if(is.null(by.x))
          warning("by.x should be provided if bugsSummary isn't a list")
  }
  effects = grep("(^R|^Fitted)", names(bugsSummary), value=TRUE)

  if(is.null(by.x)) {
    effectString = gsub("^FittedRate", "", effects)
    effectString = gsub("Spatial$", "", effectString)
    effectString = gsub("^R", "", effectString)
    effectString = unique(effectString)
    if(length(effectString) > 1) {
      warning("more than one effect level found")
    }
    if(!any(names(x) == effectString))
        stop("can't find the effect in the data provided")
       by.x = effectString
  }

  newNames = c(outer(effects, newcol, FUN=paste, sep="."))

  newcols = matrix(NA, dim(x)[1], length(newNames),
    dimnames = list(as.character(x[[by.x]]), newNames)  )
  for(D in effects) {
    colsFromSummary = newcol[newcol %in% colnames(bugsSummary[[D]])]
    newColsHere = paste(D, colsFromSummary, sep=".")
    whichrows = rownames(bugsSummary[[D]])[
        rownames(bugsSummary[[D]]) %in% rownames(newcols)
      ]
    newcols[whichrows, newColsHere] =
      bugsSummary[[D]][whichrows,colsFromSummary]
  }
  if(ncol(newcols)==1)
    colnames(newcols) = gsub("^R\\.", "", colnames(newcols))

  x[,colnames(newcols)] = newcols
  if(dim(newcols)[2]==1)
    x[,colnames(newcols)] = as.vector(x[,colnames(newcols)])  
  
  x
}
