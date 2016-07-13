get_cancer_types <-
function() {
  my_url <- "http://www.cbioportal.org/webservice.do?cmd=getTypesOfCancer"
  return(basic_tcga_query(my_url))
}
