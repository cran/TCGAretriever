\name{get_protein_info}
\alias{get_protein_info}
\title{
Retrieve Information on Antibodies Used for Protein Levels Determination
}
\description{
Retrieve information on antibodies used by reverse-phase protein arrays (RPPA) to measure protein/phosphoprotein levels.
}
\usage{
get_protein_info(csid = NULL, array_type = "protein_level", glist = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{csid}{String corresponding to the Cancer Study Identifier }
  \item{array_type}{String, c("protein_level", "phosphorylation"). Retrieve information about antibodies used for detecting total protein levels or phosphorilated levels of the protein product of the gene of interest }
  \item{glist}{ Character vector including one or more gene identifiers (ENTREZID or OFFICIAL_SYMBOL)}
}
\value{ Data Frame having one antibody per row and four columns }
\references{
http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/
}
\author{
Damiano Fantini
}
\examples{
get_protein_info("blca_tcga", glist = c("TP53", "PTEN", "E2F1", "AKT1"))
get_protein_info("blca_tcga", "phosphorilation", c("TP53", "PTEN", "E2F1", "AKT1"))
}
