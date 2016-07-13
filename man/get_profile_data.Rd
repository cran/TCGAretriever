\name{get_profile_data}
\alias{get_profile_data}
\title{
Retrieve TCGA Data corresponding to a Specific Genetic Profile of Interest
}
\description{
Retrieve Data corresponding to a Genetic Profile of interest from a given TCGA Study. This function is the workhorse of the TCGAretriever package and can be used to fetch data concerning several genes at once. For larger queries, the use of the fetch_all_tcgadata() function is mandatory
}
\usage{
get_profile_data(case_id = NULL, gprofile_id = NULL, glist = NULL)
}
\arguments{
  \item{case_id}{ String corresponding to the Identifier of a list of cases }
  \item{gprofile_id}{ String corresponding to the Identifier of a genetic Profile of interest }
  \item{glist}{ Character vector including one or more gene identifiers (ENTREZID or OFFICIAL_SYMOL)}
}
\value{
Data Frame with one row per gene and one column per case/sample
}
\references{
http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/
}
\author{
Damiano Fantini
}
\examples{
get_profile_data("blca_tcga_all", "blca_tcga_mutations", c("TP53", "E2F1"))
}