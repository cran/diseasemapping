\name{cancerRates}
\Rdversion{1.1}
\alias{cancerRates}
%- Also NEED an '\alias' for EACH other topic documented here.
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
\code{area} must be one of Camada, Norway,
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
vector of cancer rates, by age and sex group}
\examples{
\dontrun{
# won't work if offline or if the iarc web site is down
qcLungF=cancerRates(area="canada", year=2001:2002, site="lung", sex="F")
data(popdata)
qcLungExp = getSMR(popdata, qcLungF)
spplot(qcLungExp, "logExpected")
}
}
                                                       