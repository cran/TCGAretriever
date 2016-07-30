fetch_all_tcgadata <-
function (case_id = NULL, gprofile_id = NULL, glist = NULL, mutations = FALSE) {
  if (length(glist) < 500) {
    chunk_gene <- list()
    chunk_gene[["1"]] <- glist
  } else {
    chunk_gene <- split(glist, ceiling(seq_along(glist)/200))
  }
  tmp_output <- lapply(1:length(chunk_gene), (function(i) {
    if (mutations == TRUE) {
      tmp_res <- get_ext_mutation(case_id, gprofile_id, 
                                  paste(chunk_gene[[i]], collapse = "+"))
    } else {
      tmp_res <- get_profile_data(case_id, gprofile_id, 
                                  paste(chunk_gene[[i]], collapse = "+"))
    }
    if (i%%10 == 0) {
      message(paste("Genes processed: ", i * 200, "...", 
                    sep = ""))
    }
    tmp_res
  }))
  message("combining everything together...")
  final_out <- do.call(rbind, tmp_output)
  final_out <- final_out[!sapply(1:nrow(final_out), (function(i) {
    sum(is.na(final_out[i, ])) == ncol(final_out)
  })), ]
  rownames(final_out) <- NULL
  return(final_out)
}
