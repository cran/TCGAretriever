basic_tcga_query <-
function(my_url) {
  my_get <- "error"
  while (my_get[1] == "error" | my_get[2] != 200) {
    my_get <- tryCatch(httr::GET(my_url, timeout(5)), error = function(e) { "error" })
  }
  my_lines <- unlist(strsplit(httr::content(my_get, "text"), "\\n"))
  result <- do.call(rbind,sapply(1:length(my_lines), (function(x){
    strsplit(my_lines[x], "\\t")
  })))
  rownames(result) <- NULL
  # get rid of comment rows
  rows_to_keep <- !regexpr("\\#(.)+",result[,1]) == 1
  result <- result[rows_to_keep,]
  if(class(result) != "matrix") {
    result <- matrix(result, ncol = length(result), nrow = 1)
    result <- rbind(result, NA)
  }
  colnames(result) <- gsub("\\.", "-", result[1,])
  #
  if(nrow(result) == 2) {
    return(t(data.frame(result[-1,], stringsAsFactors = FALSE) ))
  } else {
    return(data.frame(result[-1,], stringsAsFactors = FALSE))
  }
}
