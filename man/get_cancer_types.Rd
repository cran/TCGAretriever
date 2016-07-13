\name{get_cancer_types}
\alias{get_cancer_types}
\title{
Retrieve a List of Cancer Types as Defined by the TCGA Guidelines
}
\description{
Retrieve information about the different types of cancer that may be included in TCGA Studies. Information include Identifier and Cancer Name.
}
\usage{
get_cancer_types()
}
\value{
Data Frame with one row per cancer type and two columns
}
\references{
http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/
}
\author{
Damiano Fantini
}
\examples{
all_canc <- get_cancer_types()
message(paste("There are", nrow(all_canc), "types on cancer defined at TCGA..."))
head(all_canc)
}
