`formatCases` <- function(casedata, ageBreaks=NULL, years=NULL, aggregate.by=NULL) {

# are there age and sex columns?
agecol = grep("^age$", names(casedata), value=TRUE, ignore.case=TRUE)
sexcol = grep("^sex$", names(casedata), value=TRUE, ignore.case=TRUE)

if(length(agecol) & length(sexcol)){
casedata$age = casedata[[agecol]]
casedata$sex = casedata[[sexcol]]
}else{
# if not, is there an age_sex_group column, in rif format?  use it to create age and sex

groupvar = grep("^AGE_SEX_GROUP$", names(casedata), value=TRUE, ignore.case=TRUE)
#agecol = grep("^age$", names(casedata), value=TRUE, ignore.case=TRUE)
#sexcol = grep("^sex$", names(casedata), value=TRUE, ignore.case=TRUE)

if(length(groupvar)){
    if(length(grep("_", casedata[[groupvar]], ignore.case=TRUE ))){
# if sex column is missing, creat it
        if(!length(sexcol)) {
            casedata$sex = substr(casedata[[groupvar]], 1, 1)
            sexcol = "sex"
            }
# if age column is missing, creat it
        if(!length(agecol)){
            n <- regexpr("_", casedata[[groupvar]], fixed = TRUE)
            casedata$age = substr(casedata[[groupvar]], n-2, n-1)
            agecol = "age"
            }
            }else{
# if groupvar is in rif format:
       if(!length(sexcol)) {
            casedata$sex = factor(substring(casedata[[groupvar]], 1, 1), levels=c("1", "2"), labels=c("M", "F"))
            sexcol = "sex"
            }
       if(!length(agecol)){
            casedata$age = 5*as.numeric(substring(casedata[[groupvar]], 2, 3))
            agecol = "age"
            }
            }
}
}


if(!is.null(ageBreaks)){
  casedata$ageNumeric = casedata[[agecol]]
 casedata$age = as.character(cut(as.numeric(as.character(casedata$ageNumeric)),
    ageBreaks$breaks, right=FALSE))
  attributes(casedata)$breaks = ageBreaks
}else{
    casedata$ageNumeric = casedata[[agecol]]
    casedata$age =  as.factor(as.character(cut(as.numeric(as.character(casedata$ageNumeric)), 
                                 sort(as.numeric(unique(casedata[[agecol]]))), right=FALSE)))
    attributes(casedata)$breaks = ageBreaks
}

# aggregate, if necessary
if(!is.null(aggregate.by) & length(aggregate.by)) {
   casedata$cases <- 1
   popa = aggregate(casedata$cases, casedata[, aggregate.by, drop=FALSE], sum)
   names(popa)[names(popa)=="x"] = "CASES"
   casedata <- popa
}


attributes(casedata)$breaks = ageBreaks
casedata
}

