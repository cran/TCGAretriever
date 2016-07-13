\name{expand_cases}
\alias{expand_cases}
\title{
Explode TCGA Case Identifiers from a TCGA Study
}
\description{
Each TCGA Study includes one or more Case Lists, that are lists of sample/patient identifiers. All Case Lists for a Study of Interest are retrieved and the individual case identifiers belonging to each of those lists are exploded and returned.  
}
\usage{
expand_cases(csid = NULL)
}
\arguments{
  \item{csid}{String corresponding to a TCGA Cancer Study identifier }
}
\value{
List containing as many elements as TCGA Case Lists available for a given TCGA Study. Each element is a list containing two elements: a) String corresponding to the Id of the case list as defined by TCGA. b)Character vector including all the case Ids belonging to the corresponding Case List. }
\references{
http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/
}
\author{
Damiano Fantini
}
\examples{
expand_cases("blca_tcga")
}

