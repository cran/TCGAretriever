get_protein_info <-
function(csid = NULL, array_type = "protein_level", glist = NULL) {
  if(!is.null(csid) & !is.null(array_type) & !is.null(glist) ){
    my_url <- "http://www.cbioportal.org/webservice.do?cmd=getProteinArrayInfo"
    my_url <- paste(my_url, "&cancer_study_id=", csid, sep = "")
    if(array_type != "protein_level"){
      array_type <- "phosphorylation"  
    }
    my_url <- paste(my_url, "&protein_array_type=", array_type, sep = "")
    my_url <- paste(my_url, "&gene_list=", paste(glist, collapse = "+"), sep = "")
    return(basic_tcga_query(my_url))
  } else {
    message("Missing cancer study ID / protein array type / gene list")
  }
}
