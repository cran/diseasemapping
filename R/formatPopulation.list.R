`formatPopulation.list`<-function(popdata, aggregate.by=NULL, breaks = NULL, 
  years=as.integer(names(popdata)), year.range=NULL,  time="YEAR", 
  personYears=TRUE, S=c("M", "F"),...) {
    
   time<-toupper(time)
    
   #If aggregate, then see if YEAR is there or not, if so, remove it
   if(!is.null(aggregate.by)){ 
#   agg<-aggregate.by<-toupper(aggregate.by)
   agg<-aggregate.by
   byYear<- time %in% aggregate.by
   if(byYear){aggregate.by<-aggregate.by[-which(aggregate.by==time)]}
   }


  
   listpop<-lapply(popdata,formatPopulation,aggregate.by, breaks= breaks)
  
   breaks = attributes(listpop[[1]])$breaks
   
   
   listdataframe<-lapply(listpop,as.data.frame)
   #if did not aggregate, then the data frames will have differnt columns
   pop<-NULL
   for (i in 1:length(listdataframe)){
    temp<-listdataframe[[i]]
    temp[,time]<-years[i]
    pop<-rbind(pop,temp)
   }
   
   attributes(pop)$breaks = breaks
   
   if(personYears){
    if (is.null(year.range)) {
      year.range = range(pop[,time])
    }
        times <- c(year.range[1], sort(years), year.range[2])
        times <- as.numeric(times)
        inter <- diff(times)/2
        nseq <- 1:length(inter) - 1
        mseq <- 2:length(inter)
        interval <- inter[mseq] + inter[nseq]
        names(interval) <- names(table(pop[,time]))
        pop$yearsForCensus = interval[as.character(pop[,time])]
        pop$POPULATION = pop$POPULATION * pop$yearsForCensus
        pop[,time] = factor(pop[,time], levels = unique(pop[,time]))
        pop[,time] = factor(pop[,time])

        pop <- pop[!is.na(pop$POPULATION),  ]
        pop <- pop[pop$POPULATION > 0,  ]
    }
   
   
   #keep desired sex in the data
sexcol = grep("^sex$", names(pop), value=TRUE, ignore.case=TRUE)
if(length(S)==1)  pop=pop[as.character(pop[[sexcol]]) %in% S,]
  
   
   pop
   #if(byYear) {p<- aggregate(pop$POPULATION,pop[,agg,drop=FALSE],sum)}
   #names(p)[names(p)=="x"] = "POPULATION"
   #p
}