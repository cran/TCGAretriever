\name{get_protein_data}
\alias{get_protein_data}
\title{Retrieve Protein Expression Data from a TCGA Study
}
\description{
TCGA includes Information about Protein Expression measured by reverse-phase protein arrays. Antibody Information can be exported together with Expression Data. All expression data will be retrieved for all available protein targets.
}
\usage{
get_protein_data(case_id = NULL, array_info = TRUE)
}
\arguments{
  \item{case_id}{ String corresponding to the Identifier of the Case List of Interest}
  \item{array_info}{Logical. If TRUE, Antibody Information will also be exported}
}
\value{
Data Frame with one gene (protein target) per row
}
\references{
http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/
}
\author{
Damiano Fantini
}
\examples{
# Protein Expression Only
get_protein_data("blca_tcga_sequenced", FALSE)[1:20,1:7]
#
# Example including Antibody Information
get_protein_data("blca_tcga_sequenced", TRUE)[1:20,1:5]
}
