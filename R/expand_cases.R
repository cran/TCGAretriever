expand_cases <-
function(csid = NULL) {
  if (!is.null(csid)){
    all_cases <- get_case_lists(csid = csid)
    lapply(1:nrow(all_cases), (function(i){
      #
      result <- list()
      result[['case_list_id']] <- as.character(all_cases[i,1])
      #
      my_cases <- as.character(all_cases[i,5])
      my_cases <- unlist(strsplit(my_cases, "TCGA"))
      my_cases <- my_cases[-1]
      case_id <- paste("TCGA", sub("[[:space:]]", "", my_cases), sep = "")
      result[['case_id']] <- case_id
      #
      result
      #
    }))
  }
}
