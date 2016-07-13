\name{basic_tcga_query}
\alias{basic_tcga_query}
\title{
TCGA Core Query Engine
}
\description{
Core Function that queries the URL provided as argument (typically a cbioportal.org URL). The function halts until the content has been completely downloaded and returns a data frame.  
}
\usage{
basic_tcga_query(my_url)
}
\arguments{
  \item{my_url}{
String. Typically, a URL pointing to the cBioPortal API. 
}
}
\details{
This is a core function invoked by other functions in the package.
}