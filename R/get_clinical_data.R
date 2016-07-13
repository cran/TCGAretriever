get_clinical_data <-
function(case_id = NULL) {
  if(!is.null(case_id)){
    my_url <- "http://www.cbioportal.org/webservice.do?cmd=getClinicalData"
    my_url <- paste(my_url, "&case_set_id=", case_id, sep = "")
    return(basic_tcga_query(my_url))
  } else {
    message("Missing case set ID")
  }
}
