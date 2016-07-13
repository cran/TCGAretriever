\name{get_case_lists}
\alias{get_case_lists}
\title{
Retrieve All Case List Available for a Specific TCGA Study
}
\description{
TCGA keeps track of which samples were analyzed by which technique within a given Study. Sample identifiers are organized in lists of cases (samples/patients) and are associated with a case_list identifier. The function retrieves information about the case lists available for a given TCGA Study.   
}
\usage{
get_case_lists(csid = NULL)
}
\arguments{
  \item{csid}{String corresponding to the Identifier of the TCGA Study of Interest}
}
\value{
Data Frame including one row per case_list and five columns 
}
\references{
http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/
}
\author{
Damiano Fantini
}
\examples{
all_case_lists <- get_case_lists("blca_tcga")
all_case_lists[,1:3]
}