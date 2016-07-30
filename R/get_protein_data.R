get_protein_data <-
function(case_id = NULL, array_info = TRUE) {
  if(!is.null(case_id) & !is.null(array_info)){
    my_url <- "http://www.cbioportal.org/webservice.do?cmd=getProteinArrayData"
    my_url <- paste(my_url, "&case_set_id=", case_id, sep = "")
    if(array_info == TRUE){
      array_ctrl <- 1  
    } else {
      array_ctrl <- 0
    }
    my_url <- paste(my_url, "&array_info=", array_ctrl, sep = "")
    return(basic_tcga_query(my_url))
  } else {
    message("Missing case study ID")
  }
}
