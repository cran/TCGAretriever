## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(graphics)
library(utils)

## ----fig.align='center', fig.width=6, fig.height=5-----------------------
library(TCGAretriever)

# Obtain a list of cancer studies from cBio
all_studies <- get_cancer_studies()

# Find published TCGA datasets
keep <- grepl("tcga_pub$", all_studies[,1])
tcga_studies <- all_studies[keep, ]

# Show results
head(tcga_studies[, 1:2])

## ----fig.align='center', fig.width=6, fig.height=5-----------------------
# Define the cancer study id: brca_tcga_pub
my_csid <- "brca_tcga_pub"

## ----fig.align='center', fig.width=6, fig.height=5-----------------------
# Obtain genetic profiles
blca_pro <- get_genetic_profiles(csid = my_csid)
head(blca_pro[, 1:2], n = 8)

## ----fig.align='center', fig.width=6, fig.height=5-----------------------
# Obtain cases 
blca_cas <- get_case_lists(csid = my_csid)
head(blca_cas[, 1:2])

## ----fig.align='center', fig.width=6, fig.height=5-----------------------
# Define a set of genes of interest
q_genes <- c("TP53", "MDM2", "E2F1", "EZH2")
q_cases <- "brca_tcga_pub_complete"
rna_prf <- "brca_tcga_pub_mrna"
mut_prf <- "brca_tcga_pub_mutations"

# Download RNA
brca_RNA <- TCGAretriever::get_profile_data(case_id = q_cases, gprofile_id = rna_prf, glist = q_genes)

## ----fig.align='center', fig.width=6, fig.height=5-----------------------
head(brca_RNA[, 1:5])

## ----fig.align='center', fig.width=6, fig.height=5-----------------------
# Set SYMBOLs as rownames
# Note that you may prefer to use the tibble package for this
rownames(brca_RNA) <- brca_RNA$COMMON
brca_RNA <- brca_RNA[, -c(1,2)]

# Download mutations (simple)
brca_MUT <- TCGAretriever::get_profile_data(case_id = q_cases, gprofile_id = mut_prf, glist = q_genes)
rownames(brca_MUT) <- brca_MUT$COMMON
brca_MUT <- brca_MUT[, -c(1,2)]

# Show results
brca_RNA[,1:6]
brca_MUT[,1:6]

## ----fig.align='center', fig.width=6, fig.height=5-----------------------
# Note that the columns (cases) are identical 
# and have the same order in both data.frames
sum(colnames(brca_MUT) != colnames(brca_RNA))

## ----fig.align='center', fig.width=5.5, fig.height=5.8-------------------
# Coerce to data.frame with numeric features 
df <- data.frame(t(brca_RNA), stringsAsFactors = FALSE)
for(i in 1:ncol(df)) { df[, i] <- as.numeric(df[, i])}

# Visualize the correlation between EZH2 and E2F1
with(df, 
     plot(E2F1, EZH2, 
          pch = 19, cex = 0.5, main = "E2F1-EZH2 correlation in BRCA"))

## ----fig.align='center', fig.width=10, fig.height=5----------------------
# Coerce to data.frame with numeric features 
df <- data.frame(t(brca_RNA), stringsAsFactors = FALSE)
for(i in 1:ncol(df)) { df[, i] <- as.numeric(df[, i])}
df$TP53.status <- as.factor(ifelse(brca_MUT["TP53",] == "NaN", "WT", "MUT"))

# Split data based on TP53.status 
lst <- split(df, f = df$TP53.status)

# Visualize the correlation between MDM2 and TP53 by P53 mutation status
par(mfrow = c(1, 2))
for(x in names(lst)) {
  with(lst[[x]], 
       plot(TP53, MDM2, 
            pch = 19, cex = 0.5, 
            xlim = c(-2.6, 2.2), ylim = c(-1.6, 3.6),
            main = paste0("MDM2-vs-P53 in ", x , " P53 tumors")))
}

