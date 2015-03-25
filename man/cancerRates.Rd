\name{cancerRates}
\alias{cancerRates}
\title{
Download cancer incidence rates from the International Agency for Research on Cancer (IARC)
}
\description{
Rates by age and sex group are retreived from http://ci5.iarc.fr/CI5plus/ci5plus.htm
}
\usage{
cancerRates(area = "canada", year=2000, sex=c("M", "F"), site="Lung") 
}
\arguments{
  \item{area}{
Region to retrieve rates from,  
}
  \item{year}{
year or sequence of years to retrieve data from, within the period 1978 to 2002
}
  \item{site}{
a vector of cancer sites, see details}
  \item{sex}{
\code{"M"} or \code{"F"} for male or female rates only, \code{c("M","F")} (the default) for both sexes.}
}
\details{
\code{area} must be one of Canada, Norway,
Latvia,
Lithuania,
Iceland,
Finland,
Estonia,
Denmark,
"Czech Republic",
"Costa Rica",
USA,
Iowa,
"New Mexico"
or the Canadian provinces of
Newfoundland, Prince Edward Island, 
Nova Scotia,
Ontario,
Manitoba,
Saskatchewan,
Alberta, and
British Columbia.  Alternately an integer specifying a registry code from http://ci5.iarc.fr.

\code{site} must be one or more of
All Sites,
Oral Cavity and Pharynx,
Oesophagus.
Stomach,
Colon,
Rectum and Anus,
Liver,
Gallbladder,
Pancreas,
Larynx,
Lung,
Bone,
Melanoma of skin,
Prostate \bold{(Males only)},
Testis \bold{(Males only)},
Breast \bold{(Females only)},
Cervix uteri \bold{(Females only)},
Corpus uteri \bold{(Females only)},
Ovary and other uterine adnexa \bold{(Females only)},
Kidney,
Bladder,
Eye,
Brain and Central Nervous System,
Thyroid,
Non-Hodgkin Lymphoma,
Hodgkin Lymphoma,
Multiple myeloma,
Leukaemia.
}

\value{
vector of cancer rates, by age and sex group
}

\examples{

# won't work if offline or if the iarc web site is down

if(interactive() | Sys.info()['user'] =='patrick') {	 
  qcLungF=cancerRates(area="canada", year=2001:2002, site="lung", sex="F")
} else {
	qcLungF = structure(c(0, 5e-06, 0, 0, 5e-06, 1e-05, 0, 3.4e-05, 9.6e-05, 
0.000211, 0.000559, 0.001289, 0.002003, 0.002508, 0.002728, 0.003189, 
0.002792, 0.001905), .Names = c("F_0", "F_5", "F_10", "F_15", 
"F_20", "F_25", "F_30", "F_35", "F_40", "F_45", "F_50", "F_55", 
"F_60", "F_65", "F_70", "F_75", "F_80", "F_85"), site = "Lung", 
area = "Canada", year = "2001-2002")
}
qcLungF

data('popdata')

qcLungExp = getSMR(popdata, qcLungF)

names(qcLungExp)

if(require('mapmisc', quietly=TRUE)) {

	mycol = colourScale(qcLungExp$expected, 
		breaks=9, dec=0, transform='log', style='equal')
	plot(popdata[1:400,])
	plot(qcLungExp, col=mycol$plot, border='#00000040',add=TRUE)
	legendBreaks('topright', mycol)
}

}
                                                       