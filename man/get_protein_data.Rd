% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TCGAretr_scr.R
\name{get_protein_data}
\alias{get_protein_data}
\title{Retrieve Protein Expression Data from a TCGA Study}
\usage{
get_protein_data(case_id = NULL, array_info = TRUE)
}
\arguments{
\item{case_id}{String corresponding to the Identifier of the Case List of Interest}

\item{array_info}{Logical. If TRUE, Antibody Information will also be exported}
}
\value{
Data Frame with one gene (protein target) per row
}
\description{
TCGA includes Information about Protein Expression measured by reverse-phase protein arrays. Antibody Information can be exported together with Expression Data. All expression data will be retrieved for all available protein targets.
}
\examples{
# Protein Expression Only
blca_protein <- get_protein_data("blca_tcga_sequenced", FALSE)
if (nrow(blca_protein) > 10 & ncol(blca_protein) > 8) {
  blca_protein[1:8,1:8]
} else {
  message("Server may be down, please try again later...")
}
#
# Example including Antibody Information
blca_protein <- get_protein_data("blca_tcga_sequenced", TRUE)
if (nrow(blca_protein) > 10 & ncol(blca_protein) > 8) {
  blca_protein[1:8,1:8]
} else {
  message("Server may be down, please try again later...")
}

}
\references{
\itemize{
  \item \url{http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/}
  \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
}
}
\author{
Damiano Fantini, \email{damiano.fantini@gmail.com}
}
