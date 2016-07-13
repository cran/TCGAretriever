\name{get_genetic_profiles}
\alias{get_genetic_profiles}
\title{
Retrieve Genetic Profiles for a TCGA Study of Interest
}
\description{
Retrieve Information about all genetic profiles associated with a TCGA Study of interest. Each TCGA Study includes one or more kind of molecular analyses whose results are referred to as genetic profiles. 
}
\usage{
get_genetic_profiles(csid = NULL)
}
\arguments{
  \item{csid}{String corresponding to the cancer study id of interest}
}
\value{
Data Frame including one row per genetic profile and six columns
}
\references{
http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/
}
\author{
Damiano Fantini
}
\examples{
get_genetic_profiles("blca_tcga")
}
