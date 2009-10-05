
CanadaCancerRates = function(area = "ontario",
   years=2000:2005,   site="colon", sex=c("M", "F")) {

  sexes = c("M"=1, "F"=2)[sex]
  
  


  areaCodes = 
  c("Canada"="00",
  "Newfoundland and Labrador"=10,
"Prince Edward Island"=11,
"Nova Scotia"=12,
"New Brunswick"=13,
"Quebec"=24,
"Ontario"=35,
"Manitoba"=46,
"Saskatchewan"=47,
"Alberta"=48,
"British Columbia"=59,
"Yukon"=60,
"Northwest Territories"=61,
"Nunavut"=62)


area = areaCodes[grep(paste("^", area[1], sep=""), names(areaCodes), ignore.case=TRUE)]


cancerCodes =   c(
      "799" ="All Sites",
               "782"="Acute Lymphocytic Leukemia",
               "784"="Acute Myeloid Leukemia",
               "718"="Anus",
               "798"="Bladder",

               "740"="Bone",
               "772"="Brain",
               "745"="Breast",
               "710"="Buccal cavity and Pharynx, Other",
               "747"="Cervix Uteri",
               "783"="Chronic Lymphocytic Leukemia",

               "785"="Chronic Myeloid Leukemia",
               "716"="Colon excluding Rectum",
               "748"="Corpus Uteri",
               "722"="Digestive System, Other",
               "776"="Endocrine Glands, Other",
               "712"="Esophagus",

               "770"="Eye",
               "705"="Floor of Mouth",
               "720"="Gallbladder",
               "751"="Genital Organs, Other, Female",
               "764"="Genital Organs, Other, Male",
               "706"="Gum and Other Mouth",

               "778"="Hodgkin Lymphoma",
               "709"="Hypopharynx",
               "456"="Kaposi Sarcoma",
               "767"="Kidney",
               "731"="Larynx",
               "786"="Leukemia, Other",

               "702"="Lip",
               "719"="Liver",
               "732"="Lung, Bronchus",
               "743"="Melanoma of Skin",
               "455"="Mesothelioma",
               "780"="Multiple Myeloma",

               "707"="Nasopharynx",
               "773"="Nervous System, Other",
               "779"="Non-Hodgkin Lymphoma",
               "708"="Oropharynx",
               "750"="Ovary",
               "721"="Pancreas",

               "763"="Penis",
               "761"="Prostate",
               "717"="Rectum",
               "733"="Respiratory System, Other",
               "704"="Salivary Glands",
               "744"="Skin, Other",

               "714"="Small Intestine",
               "741"="Soft tissue (incl. Heart)",
               "713"="Stomach",
               "762"="Testis",
               "775"="Thyroid",
               "703"="Tongue",

               "768"="Ureter",
               "769"="Urinary System, Other",
               "749"="Uterus, not Otherwise Specified",
               "787"="Other, Ill-Defined/Unknown")

site = cancerCodes[grep(paste("^", site, sep=""), cancerCodes, ignore.case=TRUE)]

codes = names(site)

  
  result = NULL
  
  for(Dyear in years) {
  
  ratesThisYear=NULL
  
  for(Dsex in names(sexes)) {
    ratesD = 0
    for(Dsite in codes) {
  
  fromScan = read.table(paste("http://dsol-smed.phac-aspc.gc.ca/dsol-smed/cancer/cgi-bin/cancerchart2/chartdata.tsv?DATA_TYPE=R&AGE_GROUPS=A%3BB%3BC%3BD%3BE%3BF%3BG%3BH%3BI%3BJ%3BK%3BL%3BM%3BN%3BO%3BP%3BQ%3BR%3BS&CAUSE2=",
    Dsite, "&YEAR2=", 
    substr(Dyear, 3, 4), "&AREA2=", area, "&SEX2=", sexes[Dsex], 
    "&CAGE2=View+Chart&SCALE=LINEAR&OUTPUT=DATA", sep=""), 
    sep="\t", quote="\"", as.is=TRUE, na.string="r.d.")

    ageGroups = paste(Dsex, fromScan[1,-1], sep="" )
    ratesD = ratesD + as.numeric(fromScan[2,-1])
    names(ratesD) = ageGroups
      
    }
    
    ratesThisYear = c(ratesThisYear, ratesD)
  }  
     result = rbind(result, ratesThisYear)
  
  }
  result[is.na(result)] = 0
  result = apply(result, 2, mean)/100000
  attributes(result)$area = area
  attributes(result)$years = years
  attributes(result)$site=site
  
  class(result)="incidenceRates"
  
  result
}