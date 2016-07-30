\name{get_clinical_data}
\alias{get_clinical_data}
\title{
Retrieve Clinical Information from a TCGA Study
}
\description{
Retrieve Information about the Patients included in a TCGA Study of Interest. Each patient is associates with a case_id. Each case_id is accompained by a set of clinical information that may include sex, age, therapeutic regimen, Tumor Staging, vital status and others. NA are allowed. 
}
\usage{
get_clinical_data(case_id = NULL)
}
\arguments{
  \item{case_id}{String corresponding to the case_list identifier of a specific list of cases of interest
}
}
\value{
Data Frame including one row per patient/case/sample
}
\references{
http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/
}
\author{
Damiano Fantini
}
\examples{
clinic_data <- get_clinical_data("blca_tcga_all")
if (nrow(clinic_data) >= 6 & ncol(clinic_data) >= 5) {
clinic_data[1:6,1:5]
hist(as.numeric(clinic_data$AGE), 
col = "darkorange", 
xlab = "Age", 
main = "Bladder Cancer, age of diagnosis")
}
}