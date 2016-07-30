\name{get_cancer_studies}
\alias{get_cancer_studies}
\title{
Retrieve a List of Cancer Studies Available at TCGA
}
\description{
Retrieve information about the different TCGA studies that are available at cBioPortal. Information include a cancer_study_id, a name of the study and a description for each study. }
\usage{
get_cancer_studies()
}
\value{
Data Frame including one study per row and three columns.
}
\references{
http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/
}
\author{
Damiano Fantini
}
\examples{
all_studies <- get_cancer_studies()
message(paste("There are", nrow(all_studies), "studies currently available..."))
if(ncol(all_studies) >= 2) {
  head(all_studies[,1:2])
}
}
