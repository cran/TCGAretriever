% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tcgaretriever_v_1_9e.R
\docType{package}
\name{TCGAretriever-package}
\alias{TCGAretriever}
\alias{TCGAretriever-package}
\title{Retrieve Genomic and Clinical Data from CBioPortal}
\description{
The Cancer Genome Atlas (TCGA) is a scientific and medical program 
aimed at improving our understanding of Cancer Biology. 
Part of the TCGA Datasets (free-access tier) are hosted
on cBioPortal, which is an open-access, open-source resource for 
interactive exploration of multidimensional cancer genomics data sets. 
TCGAretriever helps accessing and downloading TCGA data via the 
cBioPortal API. Features of TCGAretriever are: 
\itemize{
  \item it is simple and reliable
  \item it is tailored for downloading large volumes of data
}
}
\examples{
# List available Adenoid Cystic Carcinoma (acyc) Studies.
# Set `dryrun = FALSE` (default option) in production!
all_studies <- get_cancer_studies(dryrun = TRUE)
acyc_ids <- grep('acyc', all_studies$studyId, value = TRUE)
print(acyc_ids)


# List blca_tcga profiles.
# Set `dryrun = FALSE` (default option) in production!
get_genetic_profiles(csid = 'blca_tcga', dryrun = TRUE)


# List blca_tcga case lists.
# Set `dryrun = FALSE` (default option) in production!
get_case_lists(csid = 'blca_tcga', dryrun = TRUE)


# Retrieve expression data.
# Set `dryrun = FALSE` (default option) in production!
my_genes <- c('PTEN', 'TP53')
get_molecular_data(case_list_id = 'blca_tcga_3way_complete', 
                   gprofile_id = 'blca_tcga_rna_seq_v2_mrna', 
                   glist = my_genes, dryrun = TRUE)


}
\references{
\itemize{
  \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
  \item \url{http://www.cbioportal.org/} 
  \item \url{https://www.cancer.gov/ccg/research/genome-sequencing/tcga}
}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
}

}
\author{
Damiano Fantini, \email{damiano.fantini@gmail.com}
}
\keyword{internal}
