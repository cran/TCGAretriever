make_groups <-
function(num_vector, groups, group_labels = NULL, desc = FALSE) {
  if(is.numeric(num_vector) & is.numeric(groups) & length(num_vector) > 5 & groups[1] > 1){
    #
    my_order <- order(num_vector, decreasing = desc)
    groups <- as.integer(groups)
    batch_size <- as.integer(length(num_vector)/groups)
    my_ranks <- do.call(c,lapply(1:groups, (function(x){
      if(x != groups){
        rep(x,batch_size)  
      } else {
        rep(x, length(num_vector) - (length(1:(groups - 1)) * batch_size))
      }
    })))
    tmp_ranks <- rep(0, length(num_vector))
    tmp_ranks[my_order] <- my_ranks
    result <- cbind(num_vector, tmp_ranks)
    colnames(result) <- c("value", "rank")
    if(!is.null(names(num_vector))){
      rownames(result) <- names(num_vector)
    }
    result <- data.frame(result, stringsAsFactors = FALSE)
    if (!is.null(group_labels) & length(group_labels) == groups){
      lab_ranks <- sapply(tmp_ranks, (function(i){
        group_labels[i]
      }))
      result$labels <- lab_ranks
    }
    return(result)
  }
}
