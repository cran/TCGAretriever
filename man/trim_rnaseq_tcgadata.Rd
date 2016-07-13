\name{trim_rnaseq_tcgadata}
\alias{trim_rnaseq_tcgadata}
\title{
RNAseq Expression Data Trimmer
}
\description{
Groom RNAseq gene expression data frames by removing cases and/or genes having a large fraction of missing values (NaN). The ratio of missing to non-missing values is adjustable by the user. Genes having very low read counts will also be filtered out. }
\usage{
trim_rnaseq_tcgadata(rnaseq_df, smplcov = 0.95, gncov = 0.45, low_val_thr = 0.45)
}
\arguments{
  \item{rnaseq_df}{Data Frame including RNAseq expression values. Typically, the output of a fetch_all_tcgadata() call}
  \item{smplcov}{numeric, indicates the minimum non-NA ratio for sample inclusion}
  \item{gncov}{numeric, indicates the minimum non-NA ratio for gene inclusion }
  \item{low_val_thr}{numeric, indicates the minimum average expression count value for gene inclusion}
}
\value{
Data Frame resulting from trimming rnaseq_df.
}
\references{
http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/
}
\author{
Damiano Fantini
}
\examples{
my_exprs <- t(sapply(1:10, (function(i){c(i, paste("gene", i, sep = "_"))})))
my_exprs <- cbind (my_exprs, matrix(1:80, nrow = 10)/80)
my_exprs[7, 4:10] <- "NaN"
my_exprs[1:10, 7] <- "NaN"
print("Low count treshold")
trim_rnaseq_tcgadata(data.frame(my_exprs), low_val_thr = 0.3)
print("High count treshold")
trim_rnaseq_tcgadata(data.frame(my_exprs), low_val_thr = 0.5)
}
