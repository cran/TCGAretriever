## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, results = "markup", fig.align = "center", fig.width = 6, fig.height = 5)
library(TCGAretriever)
library(reshape2)
library(ggplot2)

## ----results='markup'---------------------------------------------------------
library(TCGAretriever)
library(reshape2)
library(ggplot2)

# Obtain a list of cancer studies from cBio
all_studies <- get_cancer_studies()

# Find published TCGA datasets
keep <- grepl("tcga_pub$", all_studies[,1])
tcga_studies <- all_studies[keep, ]

# Show results
show_head(tcga_studies, 6, 2)

## ----results='markup'---------------------------------------------------------
# Define the cancer study id: brca_tcga_pub
my_csid <- "brca_tcga_pub"

# Obtain genetic profiles
blca_pro <- get_genetic_profiles(csid = my_csid)
show_head(blca_pro, 8, 2)

## ----results='markup'---------------------------------------------------------
# Obtain cases 
blca_cas <- get_case_lists(csid = my_csid)
show_head(blca_cas, 8, 2)

## -----------------------------------------------------------------------------
# Define a set of genes of interest
q_genes <- c("TP53", "CLDN7", "E2F1", "EZH2")
q_cases <- "brca_tcga_pub_complete"
rna_prf <- "brca_tcga_pub_mrna"
mut_prf <- "brca_tcga_pub_mutations"

# Download Clinical Data
brca_cli <- get_clinical_data(case_id = q_cases)

# Download RNA
brca_RNA <- get_profile_data(case_id = q_cases, 
                             gprofile_id = rna_prf, 
                             glist = q_genes, 
                             force_numeric = TRUE)

## ----results='markup'---------------------------------------------------------
# Show results
show_head(brca_RNA, 4, 4)

## -----------------------------------------------------------------------------
# Set SYMBOLs as rownames
# Note that you may prefer to use the tibble package for this
rownames(brca_RNA) <- brca_RNA$COMMON
brca_RNA <- brca_RNA[, -c(1,2)]

# Round numeric vals to 3 decimals
for (i in 1:ncol(brca_RNA)) {
  brca_RNA[, i] <- round(brca_RNA[, i], digits = 3)
}

# Download mutations (simple)
brca_MUT <- get_profile_data(case_id = q_cases, 
                             gprofile_id = mut_prf, 
                             glist = q_genes)

rownames(brca_MUT) <- brca_MUT$COMMON
brca_MUT <- brca_MUT[, -c(1,2)]

# Show results
show_head(brca_RNA, 4, 4)

## -----------------------------------------------------------------------------
show_head(brca_MUT, 4, 4)

## -----------------------------------------------------------------------------
# Note that the columns (cases) are identical 
# and have the same order in both data.frames
sum(colnames(brca_MUT) != colnames(brca_RNA))

## ----fig.width=5, fig.height=5------------------------------------------------
# Visualize the correlation between EZH2 and E2F1
df <- data.frame(sample_id = colnames(brca_RNA), 
                 EZH2 = as.numeric(brca_RNA['EZH2', ]), 
                 E2F1 = as.numeric(brca_RNA['E2F1', ]), 
                 stringsAsFactors = FALSE)

ggplot(df, aes(x = EZH2, y = E2F1)) +
  geom_point(color = 'gray60', size = 0.75) +
  theme_bw() +
  geom_smooth(method = 'lm', color = 'red2', 
              size=0.3, fill = 'gray85') +
  ggtitle('E2F1-EZH2 correlation in BRCA') + 
  theme(plot.title = element_text(hjust = 0.5))

## ----fig.width=9, fig.height=5------------------------------------------------
# Coerce to data.frame with numeric features 
xpr_df <- data.frame(sample_id = colnames(brca_RNA), 
                     CLDN7 = as.numeric(brca_RNA['CLDN7', ]),
                     TP53 = as.numeric(brca_RNA['TP53', ]),
                     stringsAsFactors = FALSE)

mut_df <- data.frame(
  sample_id = colnames(brca_RNA), 
  TP53.status = as.factor(ifelse(brca_MUT["TP53",] == "NaN", "WT", "MUT")),
  stringsAsFactors = FALSE)

df <- dplyr::inner_join(xpr_df, mut_df, by='sample_id')

# Visualize the correlation between EZH2 and E2F1
ggplot(df, aes(x = TP53, y = CLDN7)) +
  geom_point(color = 'gray60', size = 0.75) +
  facet_grid(cols = vars(TP53.status)) +
  theme_bw() +
  geom_smooth(mapping = aes(color = TP53.status), 
              method = 'lm', size=0.3, fill = 'gray85') +
  ggtitle('E2F1-EZH2 correlation in BRCA') + 
  theme(plot.title = element_text(hjust = 0.5))

## ----message = FALSE, warning = FALSE, eval=TRUE------------------------------
sessionInfo()

