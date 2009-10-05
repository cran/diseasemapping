\name{CanadaCancerRates}
\Rdversion{1.1}
\alias{CanadaCancerRates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retreive cancer incidence rates from the Public Health Agency of Canada
}
\description{
Rates by age and sex group are retreived from http://dsol-smed.phac-aspc.gc.ca/dsol-smed/cancer/c_age_e.html
}
\usage{
CanadaCancerRates(area = "ontario", years = 2000:2005, site = "colon", sex=c("M", "F"))
}
\arguments{
  \item{area}{
Canadian province to retreive rates from 
}
  \item{years}{
years to retrieve data from, must be after 1992
}
  \item{site}{
a vector of cancer sites, see details}
  \item{sex}{
\code{"M"} or \code{"F"} for male or female rates only, \code{c("M","F")} (the default) for both sexes.}
}
\details{
\code{area} must be one of  Canada, Newfoundland and Labrador, Prince Edward Island, Nova Scotia,
New Brunswick,
Quebec,
Ontario,
Manitoba,
Saskatchewan,
Alberta,
British Columbia,
Yukon,
Northwest Territories,
Nunavut.

\code{site} must be one or more of
All Sites,
Acute Lymphocytic Leukemia,
Acute Myeloid Leukemia,
Anus,
Bladder,
Bone,
Brain,
Breast,
Buccal cavity and Pharynx, Other,
Cervix Uteri,
Chronic Lymphocytic Leukemia,
Chronic Myeloid Leukemia,
Colon excluding Rectum,
Corpus Uteri,
Digestive System, Other,
Endocrine Glands, Other,
Esophagus,
Eye,
Floor of Mouth,
Gallbladder,
Genital Organs, Other, Female,
Genital Organs, Other, Male,
Gum and Other Mouth,
Hodgkin Lymphoma,
Hypopharynx,
Kaposi Sarcoma,
Kidney,
Larynx,
Leukemia, Other,
Lip,
Liver,
Lung, Bronchus,
Melanoma of Skin,
Mesothelioma,
Multiple Myeloma,
Nasopharynx,
Nervous System, Other,
Non-Hodgkin Lymphoma,
Oropharynx,
Ovary,
Pancreas,
Penis,
Prostate,
Rectum,
Respiratory System, Other,
Salivary Glands,
Skin, Other,
Small Intestine,
Soft tissue (incl. Heart),
Stomach,
Testis,
Thyroid,
Tongue,
Ureter,
Urinary System, Other,
Uterus, not Otherwise Specified,
Other, Ill-Defined/Unknown.

}
\value{
vector of cancer rates, by age and sex group}
\examples{
qcLungF=CanadaCancerRates(area="quebec", years=2001:2003, site="lung", sex="F")
data(popdata)
qcLungExp = getSMR(popdata, qcLungF)
\dontrun{
spplot(qcLungExp, "logExpected")
}
}
