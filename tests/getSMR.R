library('diseasemapping')
data('kentucky')
kentucky = terra::unwrap(kentucky)

head(larynx)
10^5*larynxRates[paste(c("M","F"), 50, sep="_")]

kentucky1 = getSMR(terra::values(kentucky), larynxRates)
kentucky1[1:4,c(1,2,grep("expected", names(kentucky1),ignore.case=TRUE))]

kentucky1 = getSMR(kentucky, larynxRates)
kentucky1[1:4,c(1,2,grep("expected", names(kentucky1),ignore.case=TRUE))]

if(require('mapmisc', quietly=TRUE)) {
#  kmap = openmap(kentucky)
  col = colourScale(
      kentucky1$expected,
      style='fixed',  
      breaks=c(0:5,max(kentucky1$expected)), 
    dec=0,opacity=c(0.6,1)
  )

  plot(kentucky1, col=col$plot)
  legendBreaks('topleft', col)
}

junk = getSMR(kentucky, larynxRates, regionCode='junk')


kentucky2 = getSMR(terra::values(kentucky), larynxRates, 
    larynx, 
    regionCode="County")
kentucky2[1:4,c(1,2,grep("expected|observed", names(kentucky2),ignore.case=TRUE))]


kentucky2 = getSMR(kentucky, 
    larynxRates, 
    casedata=larynx, 
    regionCode="County")
terra::values(kentucky2)[1:4,c(1,2,grep("expected|observed", names(kentucky2),ignore.case=TRUE))]

if(require('mapmisc', quietly=TRUE)) {


  col = colourScale(
        kentucky2$observed,
        col='RdYlBu',
        style='quantile', 
        breaks=12, dec=0,opacity=c(0.6,1),
        rev=TRUE
    )
    
    plot(kentucky1, col=col$plot)
    legendBreaks('topleft', col)

    
  col = colourScale(
      kentucky2$expected,
  style='fixed', 
  col=col$col,
  breaks=col$breaks,opacity=c(0.6,1)
  )
  
  plot(kentucky2, col=col$plot)
  legendBreaks('topleft', col)
  
}


kentucky3 = getSMR(terra::values(kentucky), 
    model=list(larynxRates, larynxRates*2)
)
kentucky3[1:4,c(1,2,grep("expected|observed", names(kentucky3),ignore.case=TRUE))]

kentucky3 = getSMR(kentucky, 
    model=list('1990'=larynxRates, '1991'=larynxRates*2)
)
terra::values(kentucky3)[1:4,c(1,2,grep("expected|observed", names(kentucky3),ignore.case=TRUE))]

modelList = list()
for (D in 3:12) {
  modelList[[
      as.character(D)
      ]] = larynxRates*D/5
}


kentucky4 = getSMR(
		popdata=list(
        '5'=kentucky, '10'=kentucky
    ),
    model=modelList
)
terra::values(kentucky4[[1]])[1:4,c(1,2,grep("expected|observed", 
            names(kentucky4[[1]]),ignore.case=TRUE))
]
terra::values(kentucky4[[2]])[1:4,c(1,2,grep("expected|observed", 
            names(kentucky4[[2]]),ignore.case=TRUE))
]
 
if(require('mapmisc', quietly=TRUE)) {
  
  
  col = colourScale(
      kentucky4[[1]]$expected_6,
      col='RdYlBu',
      style='quantile', 
      breaks=15, dec=0,opacity=c(0.6,1),
      rev=TRUE
  )
  
  plot(kentucky4[[1]], col=col$plot)
  legendBreaks('topleft', col)
  
  
  col = colourScale(
      kentucky4[[2]]$expected_11,
      col='RdYlBu',
      style='quantile', 
      breaks=15, dec=0,opacity=c(0.6,1),
      rev=TRUE
  )
  
  plot(kentucky4[[2]], col=col$plot)
  legendBreaks('topleft', col)
  
}
