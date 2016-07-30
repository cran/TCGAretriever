\name{get_ext_mutation}
\alias{get_ext_mutation}
\title{
Retrieve Extended Information About DNA Mutations from TCGA
}
\description{
Query TCGA for Data about DNA Sequence Variations (Mutations) identified by exome sequencing projects. The function will retrieve an extensive set of information for each mutation that was identified in the set of cases of interest. The function can only handle a limited number of query genes. For larger queries, use the fetch_all_tcgadata() function.
}
\usage{
get_ext_mutation(case_id = NULL, gprofile_id = NULL, glist = NULL)
}
\arguments{
  \item{case_id}{String corresponding to the Identifier of the case_list of interest }
  \item{gprofile_id}{String corresponding to the Identifier of the Genetic Profile of Interest }
  \item{glist}{Character vector including Gene Identifiers (ENTREZID or OFFICIAL_SYMBOL) }
}
\value{
Data Frame inluding one row per mutation
}
\references{
http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/
}
\author{
Damiano Fantini
}
\examples{
tp53_mutats <- get_ext_mutation("blca_tcga_all", "blca_tcga_mutations", "TP53")
if(ncol(tp53_mutats) >= 6 & nrow(tp53_mutats) >= 10){
tp53_mutats[1:10,1:6]
}
}