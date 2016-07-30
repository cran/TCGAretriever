trim_rnaseq_tcgadata <-
function(rnaseq_df, smplcov = 0.95, gncov = 0.45, low_val_thr = 0.45) {
  full_matrix <- rnaseq_df
  cols_to_keep <- !sapply(1:ncol(full_matrix), (function(j){
    sum(full_matrix[,j] == "NaN")/nrow(full_matrix) > smplcov
  }))
  full_matrix <- full_matrix[,cols_to_keep]
  #
  # remove non-mapped genes
  rows_to_keep <- !sapply(1:nrow(full_matrix), (function(j){
    my_ratio <- sum(full_matrix[j,3:ncol(full_matrix)] == "NaN")/(ncol(full_matrix)-2) > gncov
    if(j%%2000 == 0) { message(paste(j, " genes were processed so far..."))}
    my_ratio
  }))
  full_matrix <- full_matrix[rows_to_keep,]
  #
  full_matrix <- full_matrix[apply(full_matrix[,-c(1,2)], 1, (function(r){sum(as.numeric(r), na.rm = T)})) > (ncol(full_matrix) - 2) * low_val_thr ,]
  return(full_matrix)
}
