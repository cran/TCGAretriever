get_cancer_studies <-
function() {
  my_url <- "http://www.cbioportal.org/webservice.do?cmd=getCancerStudies"
  return(basic_tcga_query(my_url))
}
