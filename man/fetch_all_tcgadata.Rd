\name{fetch_all_tcgadata}
\alias{fetch_all_tcgadata}
\title{
Recursively Fetch All Data Included in a TCGA Study Subset
}
\description{
Recursively query TCGA to retrieve large volumes of data corresponding to a high number of genes (up to the entire genome). Data are returned as a data frame that can be easily manipulated for further analyses. 
}
\usage{
fetch_all_tcgadata(case_id = NULL, gprofile_id = NULL, glist = NULL, mutations = FALSE)
}
\arguments{
  \item{case_id}{String corresponding to the identifier of the TCGA Case List of interest }
  \item{gprofile_id}{String corresponding to the identifier of the TCGA Profile of interest }
  \item{glist}{Character vector containing one or more gene identifiers (ENTREZID or the OFFICIAL SYMBOL can be used) }
  \item{mutations}{Logical. If TRUE, extended mutation data are fetched instead of the standard TCGA data }
}
\value{
Data Frame containing the desired TCGA data. Typically, rows are genes and columns are cases. If Extended Mutation data are retrieved (mutations = TRUE), rows correspond to individual mutations while columns are populated with details about the mutations
}
\references{
http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/
}
\author{
Damiano Fantini
}
\examples{
# Mutations occurring on TP53 and PTEN genes in the bladder cancer study
# Returns 1 data frame: rows = genes; columns = cases
fetch_all_tcgadata("blca_tcga_all", "blca_tcga_mutations", c("PTEN", "TP53"), mutation = FALSE)
# Extended mutations occurring on TP53 and PTEN genes in the bladder cancer study
# Returns 1 data frame: rows = mutations; columns = extended information
fetch_all_tcgadata("blca_tcga_all", "blca_tcga_mutations", c("PTEN", "TP53"), mutation = TRUE)
}
