get_profile_data <-
function(case_id = NULL, gprofile_id = NULL, glist = NULL){
  if(!is.null(case_id) & !is.null(gprofile_id) & !is.null(glist) ){
    my_url <- "http://www.cbioportal.org/webservice.do?cmd=getProfileData"
    my_url <- paste(my_url, "&case_set_id=", case_id, sep = "")
    my_url <- paste(my_url, "&genetic_profile_id=", paste(gprofile_id, collapse = "+"), sep = "")
    my_url <- paste(my_url, "&gene_list=", paste(glist, collapse = "+"), sep = "")
    result <- basic_tcga_query(my_url)
    colnames(result) <- gsub("\\.", "-", colnames(result))
    return(result)
  } else {
    message("Missing cancer study ID / genetic profile ID / gene list")
  }
}
