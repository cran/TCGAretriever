# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ TCGAretriever ver 1.5 ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#' TCGA Core Query Engine
#' 
#' Core Function that queries the URL provided as argument (typically a cbioportal.org URL). 
#' The function halts until the content has been completely downloaded and returns a data frame.
#'   
#' @param my_url string. Typically, a URL pointing to the cBioPortal API. 
#' 
#' @details This is a core function invoked by other functions in the package.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \itemize{
#'   \item \url{http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/}
#'   \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' }
#'  
#' @importFrom httr GET content timeout
#' 
#' @keywords intenal
basic_tcga_query <- function(my_url) {
  
  curWarn <- options()$warn
  options(warn = -1)
  # options(warn = curWarn)
  
  my_get <- "error"
  while (my_get[1] == "error" | my_get[2] != 200) {
    my_get <- tryCatch(httr::GET(my_url, httr::timeout(5)), error = function(e) { "error" })
  }
  my_lines <- unlist(strsplit(httr::content(my_get, "text"), "\\n"))
  result <- do.call(rbind,sapply(1:length(my_lines), (function(x){
    strsplit(my_lines[x], "\\t")
  })))
  rownames(result) <- NULL
  # get rid of comment rows
  rows_to_keep <- !regexpr("\\#(.)+",result[,1]) == 1
  result <- result[rows_to_keep,]
  
  #if(class(result) != "matrix") {
  if(!is.matrix(result)) {  
    result <- matrix(result, ncol = length(result), nrow = 1)
    result <- rbind(result, NA)
  }
  colnames(result) <- gsub("\\.", "-", result[1,])
  #
  options(warn = curWarn)
  if(nrow(result) == 2) {
    return(t(data.frame(result[-1,], stringsAsFactors = FALSE) ))
  } else {
    return(data.frame(result[-1,], stringsAsFactors = FALSE))
  }
}


#' Explode TCGA Case Identifiers from a TCGA Study
#' 
#' Each TCGA Study includes one or more "case lists". These are lists of 
#' sample/patient identifiers. All case lists of a study of interest are 
#' retrieved and the individual case identifiers are expanded and returned  
#' 
#' @param csid string corresponding to a TCGA Cancer Study identifier
#' 
#' @return list containing as many elements as TCGA case lists 
#' available for a given TCGA Study. Each element is a list containing two elements: 
#' \itemize{
#'   \item a string corresponding to the Id of the case list as defined by TCGA
#'   \item character vector including all case IDs corresponding to the case list
#' }
#' 
#' @examples 
#' expand_cases("blca_tcga")
#' 
#' @export
expand_cases <-function(csid = NULL) {
  
  curWarn <- options()$warn
  options(warn = -1)
  # options(warn = curWarn)
  
  OUT <- NULL
  
  if (!is.null(csid)){
    all_cases <- get_case_lists(csid = csid)
    OUT <- lapply(1:nrow(all_cases), (function(i){
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
  
  options(warn = curWarn)
  return(OUT)
}


#' Recursively Fetch All Data Included in a TCGA Study Subset
#' 
#' Recursively query TCGA to retrieve large volumes of data corresponding to a 
#' high number of genes (up to the entire genome). Data are returned as a data frame 
#' that can be easily manipulated for further analyses. 
#' 
#' @param case_id string corresponding to the identifier of the TCGA Case List of interest 
#' @param gprofile_id string corresponding to the identifier of the TCGA Profile of interest 
#' @param glist character vector including one or more gene identifiers (ENTREZID or the OFFICIAL SYMBOL can be used) 
#' @param mutations logical. If TRUE, extended mutation data are fetched instead of the standard TCGA data 
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \itemize{
#'   \item \url{http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/}
#'   \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' }
#'  
#' 
#' @return 
#' A data.frame is returned, including the desired TCGA data. 
#' Typically, rows are genes and columns are cases. 
#' If "extended mutation" data are retrieved (mutations = TRUE), 
#' rows correspond to individual mutations 
#' while columns are populated with mutation features
#' 
#' 
#' @examples 
#' # Mutations occurring on TP53 and PTEN genes in the bladder cancer study
#' # Returns 1 data frame: rows = genes; columns = cases
#' fetch_all_tcgadata("blca_tcga_all", "blca_tcga_mutations", c("PTEN", "TP53"), mutation = FALSE)
#' # Extended mutations occurring on TP53 and PTEN genes in the bladder cancer study
#' # Returns 1 data frame: rows = mutations; columns = extended information
#' fetch_all_tcgadata("blca_tcga_all", "blca_tcga_mutations", c("PTEN", "TP53"), mutation = TRUE)
#' 
#' 
#' @export
fetch_all_tcgadata <- function (case_id = NULL, gprofile_id = NULL, glist = NULL, mutations = FALSE) {
  
  curWarn <- options()$warn
  options(warn = -1)
  # options(warn = curWarn)
  
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
  
  options(warn = curWarn)
  return(final_out)
}


#' Retrieve a List of Cancer Studies Available at TCGA
#' 
#' Retrieve information about the different TCGA studies that are available at cBioPortal. 
#' Information include a cancer_study_id, a name of the study and a description for each study. 
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \itemize{
#'   \item \url{http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/}
#'   \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' }
#'  
#' @return 
#' Data Frame including one study per row and three columns.
#' 
#' @examples 
#' all_studies <- get_cancer_studies()
#' message(paste("There are", nrow(all_studies), "studies currently available..."))
#' if(ncol(all_studies) >= 2) {
#'   head(all_studies[,1:2])
#' }
#' 
#' @export
get_cancer_studies <- function() {
  
  curWarn <- options()$warn
  options(warn = -1)
  # options(warn = curWarn)
  
  OUT <- NULL
  my_url <- "http://www.cbioportal.org/webservice.do?cmd=getCancerStudies"
  OUT <- try({basic_tcga_query(my_url)}, silent = TRUE)
  
  options(warn = curWarn)
  return(OUT)
}


#' Retrieve a List of Cancer Types as Defined by the TCGA Guidelines
#' 
#' Retrieve information about the different types of cancer that may be 
#' included in TCGA Studies. Information include Identifier and Cancer Name.
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \itemize{
#'   \item \url{http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/}
#'   \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' }
#'  
#' 
#' @return 
#' A data.frame with one row per cancer type and two columns
#' 
#' @examples
#' all_canc <- get_cancer_types()
#' message(paste("There are", nrow(all_canc), "types on cancer defined at TCGA..."))
#' head(all_canc)
#' 
#' @export
get_cancer_types <- function() {
  
  curWarn <- options()$warn
  options(warn = -1)
  # options(warn = curWarn)
  
  OUT <- NULL
  my_url <- "http://www.cbioportal.org/webservice.do?cmd=getTypesOfCancer"
  OUT <- try({basic_tcga_query(my_url)}, silent = TRUE)
  
  options(warn = curWarn)
  return(OUT)
}


#' Retrieve All Case List Available for a Specific TCGA Study
#' 
#' TCGA keeps track of which samples were analyzed by which technique within a 
#' given Study. Sample identifiers are organized in lists of cases (samples/patients) 
#' and are associated with a case_list identifier. The function retrieves information 
#' about the case lists available for a given TCGA Study.   
#' 
#' @param csid String corresponding to the Identifier of the TCGA Study of Interest
#' 
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \itemize{
#'   \item \url{http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/}
#'   \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' }
#'  
#' @return Data Frame including one row per case_list and five columns 
#' 
#' @examples 
#' all_case_lists <- get_case_lists("blca_tcga")
#' if(ncol(all_case_lists) >= 3) {
#'   all_case_lists[,1:3]
#' }
#' 
#' @export
get_case_lists <- function(csid = NULL) {
  
  curWarn <- options()$warn
  options(warn = -1)
  # options(warn = curWarn)
  
  if(!is.null(csid)){
    my_url <- "http://www.cbioportal.org/webservice.do?cmd=getCaseLists"
    my_url <- paste(my_url, "&cancer_study_id=", tolower(as.character(csid)), sep = "")
    OUT <- basic_tcga_query(my_url)
  } else {
    message("Missing cancer study ID")
    OUT <- NULL
  }
  
  options(warn = curWarn)
  return(OUT)
}


#' Retrieve Clinical Information from a TCGA Study
#' 
#' Retrieve Information about the Patients included in a TCGA Study of Interest. 
#' Each patient is associates with a case_id. Each case_id is accompained by a set 
#' of clinical information that may include sex, age, therapeutic regimen, Tumor Staging, 
#' vital status and others. NA are allowed. 
#' 
#' @param case_id string corresponding to the case_list identifier of a specific list of cases of interest
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \itemize{
#'   \item \url{http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/}
#'   \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' }
#'
#' 
#' @return data.frame including one row per patient/case/sample
#' 
#' @examples 
#' clinic_data <- get_clinical_data("blca_tcga_all")
#' if (nrow(clinic_data) >= 6 & ncol(clinic_data) >= 5) {
#'   clinic_data[1:6,1:5]
#'   hist(as.numeric(clinic_data$AGE), 
#'   col = "darkorange", 
#'   xlab = "Age", 
#'   main = "Bladder Cancer, age of diagnosis")
#' }
#' 
#' @export
get_clinical_data <- function(case_id = NULL) {
  
  curWarn <- options()$warn
  options(warn = -1)
  # options(warn = curWarn)
  
  if(!is.null(case_id)){
    my_url <- "http://www.cbioportal.org/webservice.do?cmd=getClinicalData"
    my_url <- paste(my_url, "&case_set_id=", case_id, sep = "")
    OUT <- basic_tcga_query(my_url)
  } else {
    message("Missing case set ID")
    OUT <- NULL
  }

  options(warn = curWarn)
  return(OUT)
}



#' Retrieve Extended Information About DNA Mutations from TCGA
#' 
#' Query TCGA for Data about DNA Sequence Variations (Mutations) identified by exome sequencing projects. 
#' The function will retrieve an extensive set of information for each mutation that was identified in 
#' the set of cases of interest. The function can only handle a limited number of query genes. 
#' For larger queries, use the fetch_all_tcgadata() function.
#'
#' @param case_id string corresponding to the Identifier of the case_list of interest 
#' @param gprofile_id string corresponding to the Identifier of the Genetic Profile of Interest 
#' @param glist character vector including Gene Identifiers (ENTREZID or OFFICIAL_SYMBOL) 
#' 
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \itemize{
#'   \item \url{http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/}
#'   \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' }
#'  
#' 
#' @return data Frame inluding one row per mutation
#' 
#' @examples 
#' tp53_mutats <- get_ext_mutation("blca_tcga_all", "blca_tcga_mutations", "TP53")
#' if(ncol(tp53_mutats) >= 6 & nrow(tp53_mutats) >= 10){
#'   tp53_mutats[1:10,1:6]
#' }
#' 
#' @export
get_ext_mutation <- function(case_id = NULL, gprofile_id = NULL, glist = NULL){
  
  curWarn <- options()$warn
  options(warn = -1)
  # options(warn = curWarn)
  
  if(!is.null(gprofile_id) & !is.null(glist) ){
    my_url <- "http://www.cbioportal.org/webservice.do?cmd=getMutationData"
    my_url <- paste(my_url, "&case_set_id=", case_id, sep = "")
    my_url <- paste(my_url, "&genetic_profile_id=", paste(gprofile_id, collapse = "+"), sep = "")
    my_url <- paste(my_url, "&gene_list=", paste(glist, collapse = "+"), sep = "")

    OUT <- basic_tcga_query(my_url)
  } else {
    message("Missing cancer study ID / genetic profile ID / gene list")
    OUT <- NULL
  }
  
  options(warn = curWarn)
  return(OUT)
}



#' Retrieve Genetic Profiles for a TCGA Study of Interest
#' 
#' Retrieve Information about all genetic profiles associated with a TCGA Study of interest.
#' Each TCGA Study includes one or more kind of molecular analyses whose results are 
#' referred to as genetic profiles. 
#' 
#' @param csid string corresponding to the cancer study id of interest
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \itemize{
#'   \item \url{http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/}
#'   \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' }
#' 
#' @return 
#' data.frame including one row per genetic profile and six columns
#' 
#' @examples 
#' get_genetic_profiles("blca_tcga")
#' 
#' @export
get_genetic_profiles <- function(csid = NULL){
  
  curWarn <- options()$warn
  options(warn = -1)
  # options(warn = curWarn)
  
  if(!is.null(csid)){
    my_url <- "http://www.cbioportal.org/webservice.do?cmd=getGeneticProfiles"
    my_url <- paste(my_url, "&cancer_study_id=", tolower(as.character(csid)), sep = "")
    OUT <- basic_tcga_query(my_url)
  } else {
    message("Missing cancer study ID")
    OUT <- NULL
  }
  
  options(warn = curWarn)
  return(OUT)
}


#' Retrieve TCGA Data corresponding to a Specific Genetic Profile of Interest
#' 
#' Retrieve Data corresponding to a Genetic Profile of interest from a given TCGA Study. 
#' This function is the workhorse of the TCGAretriever package and can be used to fetch data 
#' concerning several genes at once. For larger queries, the use of the fetch_all_tcgadata() 
#' function is mandatory
#' 
#' @param case_id String corresponding to the Identifier of a list of cases 
#' @param gprofile_id String corresponding to the Identifier of a genetic Profile of interest 
#' @param glist Character vector including one or more gene identifiers (ENTREZID or OFFICIAL_SYMOL)
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \itemize{
#'   \item \url{http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/}
#'   \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' }
#'  
#' @return 
#' data.frame with one row per gene and one column per case/sample
#' 
#' @examples 
#' get_profile_data("blca_tcga_all", "blca_tcga_mutations", c("TP53", "E2F1"))
#' 
#' @export
get_profile_data <- function(case_id = NULL, gprofile_id = NULL, glist = NULL){
  
  curWarn <- options()$warn
  options(warn = -1)
  # options(warn = curWarn)
  
  if(!is.null(case_id) & !is.null(gprofile_id) & !is.null(glist) ){
    my_url <- "http://www.cbioportal.org/webservice.do?cmd=getProfileData"
    my_url <- paste(my_url, "&case_set_id=", case_id, sep = "")
    my_url <- paste(my_url, "&genetic_profile_id=", paste(gprofile_id, collapse = "+"), sep = "")
    my_url <- paste(my_url, "&gene_list=", paste(glist, collapse = "+"), sep = "")
    result <- basic_tcga_query(my_url)
    colnames(result) <- gsub("\\.", "-", colnames(result))
  } else {
    message("Missing cancer study ID / genetic profile ID / gene list")
    result <- NULL
  }
  
  options(warn = curWarn)
  return(result)
}


#' Retrieve Protein Expression Data from a TCGA Study
#' 
#' TCGA includes Information about Protein Expression measured by reverse-phase protein arrays. Antibody Information can be exported together with Expression Data. All expression data will be retrieved for all available protein targets.
#' 
#' @param case_id String corresponding to the Identifier of the Case List of Interest
#' @param array_info Logical. If TRUE, Antibody Information will also be exported
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \itemize{
#'   \item \url{http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/}
#'   \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' }
#'  
#' 
#' @return 
#' Data Frame with one gene (protein target) per row
#' 
#' @examples 
#' # Protein Expression Only
#' blca_protein <- get_protein_data("blca_tcga_sequenced", FALSE)
#' if (nrow(blca_protein) > 10 & ncol(blca_protein) > 8) {
#'   blca_protein[1:8,1:8]
#' } else {
#'   message("Server may be down, please try again later...")
#' }
#' #
#' # Example including Antibody Information
#' blca_protein <- get_protein_data("blca_tcga_sequenced", TRUE)
#' if (nrow(blca_protein) > 10 & ncol(blca_protein) > 8) {
#'   blca_protein[1:8,1:8]
#' } else {
#'   message("Server may be down, please try again later...")
#' }
#' 
#' @export
get_protein_data <- function(case_id = NULL, array_info = TRUE) {
  
  curWarn <- options()$warn
  options(warn = -1)
  # options(warn = curWarn)
  
  if(!is.null(case_id) & !is.null(array_info)){
    my_url <- "http://www.cbioportal.org/webservice.do?cmd=getProteinArrayData"
    my_url <- paste(my_url, "&case_set_id=", case_id, sep = "")
    if(array_info == TRUE){
      array_ctrl <- 1  
    } else {
      array_ctrl <- 0
    }
    my_url <- paste(my_url, "&array_info=", array_ctrl, sep = "")
    OUT <- basic_tcga_query(my_url)
  } else {
    message("Missing case study ID")
    OUT <- NULL
  }
  
  options(warn = curWarn)
  return(OUT)
}


#' Retrieve Information on Antibodies Used for Protein Levels Determination
#' 
#' Retrieve information on antibodies used by reverse-phase protein arrays (RPPA) 
#' to measure protein/phosphoprotein levels.
#' 
#' @param csid String corresponding to the Cancer Study Identifier 
#' @param array_type String, c("protein_level", "phosphorylation"). Retrieve information 
#' about antibodies used for detecting total protein levels or phosphorilated levels of 
#' the protein product of the gene of interest 
#' @param glist Character vector including one or more gene identifiers (ENTREZID or OFFICIAL_SYMBOL)
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \itemize{
#'   \item \url{http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/}
#'   \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' }
#' 
#' @return data frame having one antibody per row and four columns 
#' 
#' @examples 
#' info1 <- get_protein_info("blca_tcga", glist = c("TP53", "PTEN", "E2F1", "AKT1"))
#' if (nrow(info1) > 0) {
#'   message("Total protein levels information")
#'   info1
#' } else {
#'   message("Server may be down, please try again later...")
#' }
#' #
#' info2 <- get_protein_info("blca_tcga", "phosphorilation", c("TP53", "PTEN", "E2F1", "AKT1"))
#' if (nrow(info2) > 0) {
#'   message("Phospho-protein levels information")
#'   info2
#' } else {
#'   message("Server may be down, please try again later...")
#' }
#' 
#' @export
get_protein_info <- function(csid = NULL, array_type = "protein_level", glist = NULL) {
  
  curWarn <- options()$warn
  options(warn = -1)
  # options(warn = curWarn)
  
  if(!is.null(csid) & !is.null(array_type) & !is.null(glist) ){
    my_url <- "http://www.cbioportal.org/webservice.do?cmd=getProteinArrayInfo"
    my_url <- paste(my_url, "&cancer_study_id=", csid, sep = "")
    if(array_type != "protein_level"){
      array_type <- "phosphorylation"  
    }
    my_url <- paste(my_url, "&protein_array_type=", array_type, sep = "")
    my_url <- paste(my_url, "&gene_list=", paste(glist, collapse = "+"), sep = "")
    OUT <- basic_tcga_query(my_url)
  } else {
    message("Missing cancer study ID / protein array type / gene list")
    OUT <- NULL
  }
  
  options(warn = curWarn)
  return(OUT)
}


#' Split Numeric Vectors in Groups
#' 
#' Assign each element of a numeric vector to a group. Grouping is based on ranks: 
#' numeric values are sorted and then split in 2 or more groups. Values may be sorted 
#' in an increasing or decreasing fashion. The vector is returned in the original order. 
#' Labels may be assigned to each groug.
#' 
#' @param num_vector numeric vector. It includes the values to be assigned to the different groups
#' @param groups integer. The number of groups that will be generated
#' @param group_labels character vector. Labels for each group. Note that 
#' the length of group_labels has to be equal to the number of groups
#' @param desc logical. If TRUE, the sorting is applied in a decreasing fashion
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \itemize{
#'   \item \url{http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/}
#'   \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' }
#'  
#' @return 
#' data.frame including the vector provided as argument in the original order ("value") 
#' and the grouping vector ("rank"). If labels are provided as an argument, group labels 
#' are also included in the data.frame ("labels"). 
#' 
#' @examples
#' exprs_geneX <- c(19.1,18.4,22.4,15.5,20.2,17.4,9.4,12.4,31.2,33.2,18.4,22.1)
#' groups_num <- 3
#' groups_labels <- c("high", "med", "low")
#' make_groups(exprs_geneX, groups_num, groups_labels, desc = TRUE)
#' 
#' @export
make_groups <- function(num_vector, groups, group_labels = NULL, desc = FALSE) {
  
  curWarn <- options()$warn
  options(warn = -1)
  # options(warn = curWarn)
  
  result <- NULL
  
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
  } 
  
  options(warn = curWarn)
  return(result)
}



#' Retrieve Genomic and Clinical Data from cBioPortal
#'
#' @description 
#' The Cancer Genome Atlas (TCGA) is a program aimed at improving our understanding 
#' of Cancer Biology. Several TCGA Datasets are available online. TCGAretriever 
#' helps accessing and downloading TCGA data via the cBioPortal API. 
#' Features of TCGAretriever are: 
#' \itemize{
#'   \item it is simple and reliable
#'   \item it is tailored for downloading large volumes of data
#' }
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \itemize{
#'   \item \url{http://www.biotechworld.it/bioinf/2016/07/11/tcga-data-via-tcgaretriever/}
#'   \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#'   \item \url{http://www.cbioportal.org/} 
#'   \item \url{http://cancergenome.nih.gov/abouttcga/}
#' }
#'  
#' 
#' @keywords internal
"_PACKAGE"


# Internal scr // Roxygenize
#
# 
# setwd("~/Documents/r_pack_dev/TCGAretriever/ver_1_5/build_pack/")
# package.skeleton(name = "TCGAretriever", code_files = "TCGAretr_scr.R")
# dir.create("TCGAretriever/vignettes/")
# file.copy(from = "DESCRIPTION", to = "TCGAretriever/DESCRIPTION", overwrite = TRUE)
# file.copy(from = "getting_started_with_TCGAretriever.Rmd", to = "TCGAretriever/vignettes/getting_started_with_TCGAretriever.Rmd", overwrite = TRUE)
# file.remove("TCGAretriever/NAMESPACE")
# file.remove("TCGAretriever/Read-and-delete-me")
# file.remove(paste0("TCGAretriever/man/", dir("TCGAretriever/man/")))
# roxygen2::roxygenize("TCGAretriever/")



  