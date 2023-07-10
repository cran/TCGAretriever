# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~ TCGAretriever ver 1.7 ~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Return a Subset of an Object 
#'
#'
#' Show the first part of an object (vector, matrix, data frame, or list). For 
#' two-dimensional objects, the operation returns up to a certain number of 
#' rows and columns as indicated by the user.
#' 
#'  
#' @param x an object.
#' @param i an integer vector of length 1. Maximum number of rows (or elements)
#' to be shown. 
#' @param j an integer vector of length 1. Maximum number of columns 
#' (or inner elements) to be shown. 
#' 
#' @return an object (subset) of the same class as `x`. NULL is returned if the
#' header cannot be extracted.
#' 
#' @details This function is a simple error-compliant version 
#' of a header function. NULL is returned if a header cannot be computed (e.g., 
#' if a funciton is passed as object `x`.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' 
#'  
#' @examples 
#' my_x <- data.frame(A=1:5, B=2:6, C=3:7)
#' show_head(my_x, 2, 2)
#' 
#' 
#' @export
show_head <- function(x, i=6, j=6) {
  
  inner_fx_1d <- function(x, k) {
    
    stopifnot(is.list(x) | is.data.frame(x) | is.matrix(x) | is.vector(x), 
              is.numeric(k), length(k) == 1, !is.na(k), k > 0)
    
    if (is.list(x) || is.vector(x)) {
      y <- x[seq(1, k, by = 1)]
    } else {
      y <- x[seq(1, k, by = 1), ]
    }
    return(y)
  }

  inner_fx_2d <- function(x, i, j) {
    
    stopifnot(is.list(x) | is.data.frame(x) | is.matrix(x) | is.vector(x), 
              is.numeric(i), length(i) == 1, !is.na(i), i > 0, 
              is.numeric(j), length(j) == 1, !is.na(j), j > 0)
    
    i <- round(i, digits = 0)
    j <- round(j, digits = 0)
    
    if (is.data.frame(x) || is.matrix(x)) {
      # PART 0
      exp_c <- NULL
      exp_r <- NULL
      # PART 1
      if (ncol(x) > 0) {
        exp_c <- seq(1, j, by = 1)
        exp_c <- exp_c[exp_c <= ncol(x)]
      }
      # PART 2
      if (nrow(x) > 0) {
        exp_r <- seq(1, i, by = 1)
        exp_r <- exp_r[exp_r <= nrow(x)]
      }
      
      if (!is.null(exp_r) && !is.null(exp_c)) {
        y <- x[exp_r, exp_c]
      } else if (!is.null(exp_r)) {
        y <- x[exp_r, ]
      } else if (!is.null(exp_c)) {
        y <- x[, exp_c]
      } else {
        y <- x
      }
      
    } else if (is.list(x)) {
      
      if (length(x) > 0) {
        exp_c <- seq(1, i, by = 1)
        exp_c <- exp_c[exp_c <= length(x)]
        y <- x[exp_c]
  
        # Nested
        for (k in seq(1, length(y), by = 1)) {
          y[[k]] <- inner_fx_1d(y[[k]], k=j)
          message('.')
        }
        
      } else {
        y <- x
      }
     
    } else if (is.vector(x)) {
      y <- inner_fx_1d(x, k=i)
      
      
    } else {
      y <- NULL
    }
    
    return(y)
  }
  
  # Run
  y <- tryCatch({inner_fx_2d(x, i, j)}, error = function(e) { NULL })
  return(y)
}





#' TCGA Core Query Engine
#' 
#' Core Function that queries the URL provided as argument (typically a cbioportal.org URL). 
#' The function halts until the content has been completely downloaded and returns a data frame.
#'   
#' @param my_url string. Typically, a URL pointing to the cBioPortal API. 
#' 
#' @return data.frame including data retrieved from cBioPortal.
#' 
#' @details This is a core function invoked by other functions in the package.
#'
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' 
#'  
#' @examples 
#' my_url <- "http://www.cbioportal.org/webservice.do?cmd=getCancerStudies"
#' x <- TCGAretriever:::basic_tcga_query(my_url)
#' show_head(x, 5, 2) 
#' 
#'  
#' @importFrom httr GET content timeout
#' 
#' @keywords intenal
basic_tcga_query <- function(my_url) {
  
  curWarn <- options()$warn
  options(warn = -1)
  my_get <- c("error", 100)
  
  y <- tryCatch({
    
    while (my_get[1] == "error" | my_get[2] != 200) {
      my_get <- tryCatch(httr::GET(my_url, httr::timeout(5)), 
                         error = function(e) { "error" }) }
    
    my_lines <- strsplit(httr::content(my_get, "text"), split = '\n')[[1]]
    my_lines <- my_lines[!grepl('^#', my_lines)]
    my_lines <- strsplit(my_lines, split = "\t")
    result <- do.call(rbind, my_lines[-1])
    rownames(result) <- NULL
    exp_colnames <- as.character(my_lines[[1]])
    
    # Sometimes, exp_colnames may NOT match the expected dimension
    problm_cl <- c('xvar_link', 'xvar_link_pdb', 
                   'functional_impact_score', 'xvar_link_msa')
    if (length(exp_colnames) == ncol(result)) {
      colnames(result) <- exp_colnames
    } else if (sum(problm_cl %in% exp_colnames) > 0) {
      new_colnames <- exp_colnames[!exp_colnames %in% problm_cl]
      miss_clms <- ncol(result) - length(new_colnames) 
      if (miss_clms > 0) {
        new_colnames <- c(new_colnames, paste0('extra_col_', 1:miss_clms))  
      }
      colnames(result) <- new_colnames
    } else if (length(exp_colnames) > ncol(result)) {
      colnames(result) <- exp_colnames[1:ncol(result)]
    } 
    as.data.frame(result)
  }, error = function(e) { data.frame() })
    
  return(y)
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
#' expand_cases("mel_tsam_liang_2017")
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' 
#' 
#' @export
expand_cases <-function(csid = "mel_tsam_liang_2017") {

  stopifnot(!is.null(csid),
            length(csid) == 1,
            !is.na(csid), 
            is.character(csid))
  
  curWarn <- options()$warn
  options(warn = -1)

  OUT <- tryCatch({
    all_cases <- get_case_lists(csid = csid)
    
    TMP <- lapply(1:nrow(all_cases), (function(i){
      #
      result <- list()
      result[['case_list_id']] <- as.character(all_cases[i,1])
      #
      my_cases <- as.character(all_cases[i,5])
      my_cases <- strsplit(my_cases, "[[:space:]]")[[1]]
      my_cases <- gsub('(^[[:space:]]+)|([[:space:]]+$)', '', my_cases)
      my_cases <- my_cases[nchar(my_cases) > 0]
      my_cases <- my_cases[-1]
      result[['case_id']] <- my_cases
      result
    }))
    TMP
  }, error = function(e) { list() })
  
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
#' @param force_numeric logical. Shall columns including numeric values be coerced to numeric.
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' 
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
#' x <- fetch_all_tcgadata("blca_tcga_mutations", "blca_tcga_mutations", 
#'                         c("PTEN", "TP53"), mutation = FALSE)
#' # Extended mutations occurring on TP53 and PTEN genes in the bladder cancer study
#' # Returns 1 data frame: rows = mutations; columns = extended information
#' fetch_all_tcgadata("blca_tcga_all", "blca_tcga_mutations", c("PTEN", "TP53"), mutation = TRUE)
#' 
#' 
#' @export
fetch_all_tcgadata <- function (case_id = 'blca_tcga_sequenced', 
                                gprofile_id = 'blca_tcga_mutations', 
                                glist = c('PTEN', 'TP53'), 
                                mutations = FALSE, force_numeric = FALSE) {
  
  # Check 1  
  stopifnot(!is.null(case_id), is.character(case_id), length(case_id) == 1, 
            !is.na(case_id), nchar(case_id) > 0,
            !is.null(gprofile_id), is.character(gprofile_id), 
            length(gprofile_id) == 1, !is.na(gprofile_id), 
            nchar(gprofile_id) > 0, !is.null(glist), 
            is.vector(glist), length(glist) > 0, is.logical(mutations))

  # Check 2
  glist <- glist[!is.na(glist)]
  glist <- glist[nchar(glist) > 0]
  glist <- unique(glist)
  stopifnot(length(glist) > 0)

  # proceed  
  curWarn <- options()$warn
  options(warn = -1)

  # Handle batches
  if (length(glist) < 500) {
    chunk_gene <- list()
    chunk_gene[["1"]] <- glist
  } else {
    chunk_gene <- split(glist, ceiling(seq_along(glist)/200))
  }
  
  tmp_output <- lapply(1:length(chunk_gene), (function(i) {
    
    if (mutations == TRUE) {
      tmp_res <- get_ext_mutation(
        case_id = case_id, gprofile_id = gprofile_id, 
        glist = paste(chunk_gene[[i]], collapse = "+"))
    } else {
      tmp_res <- get_profile_data(
        case_id = case_id, gprofile_id = gprofile_id, 
        glist = paste(chunk_gene[[i]], collapse = "+"), 
        force_numeric = FALSE)
    }
    
    if (i %% 10 == 0) {
      message(paste("Genes processed: ", i * 200, "...", 
                    sep = ""))
    }
    
    tmp_res
  }))
  message("combining everything together...")
  final_out <- do.call(rbind, tmp_output)
  r2rem <- do.call(c, lapply(1:nrow(final_out), (function(i) {
    sum(is.na(final_out[i, ])) == ncol(final_out)
  })))
  final_out <- final_out[!r2rem, ]
  rownames(final_out) <- NULL
  
  if (force_numeric){
    if (ncol(final_out) > 0 && nrow(final_out) > 0) {
      for (i in seq(1, ncol(final_out), by = 1)) {
        v_i <- (!is.na(suppressWarnings(as.numeric(final_out[, i])))) | 
          final_out[, i] == 'NaN'
        rat_i <- sum(v_i, na.rm = TRUE) / nrow(final_out)
        if (rat_i >= 0.7) {
          final_out[, i] <- suppressWarnings(as.numeric(final_out[, i]))
        }
      }
    }
  }
  
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
#' \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' 
#'  
#' @return 
#' Data Frame including one study per row and three columns.
#' 
#' @examples 
#' all_studies <- get_cancer_studies()
#' message(paste("There are", nrow(all_studies), "studies currently available..."))
#' 
#' @export
get_cancer_studies <- function() {
  
  curWarn <- options()$warn
  options(warn = -1)
  my_url <- "http://www.cbioportal.org/webservice.do?cmd=getCancerStudies"
  OUT <- tryCatch({as.data.frame(basic_tcga_query(my_url))}, 
                  error = function(e) { data.frame()})

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
#' \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' 
#'  
#' 
#' @return 
#' A data.frame with one row per cancer type and two columns
#' 
#' @examples
#' all_canc <- get_cancer_types()
#' message(paste("There are", nrow(all_canc), "cancer indications..."))
#' 
#' @export
get_cancer_types <- function() {
  
  curWarn <- options()$warn
  options(warn = -1)
  my_url <- "http://www.cbioportal.org/webservice.do?cmd=getTypesOfCancer"
  OUT <- tryCatch({as.data.frame(basic_tcga_query(my_url))},
                  error = function(e) { data.frame()})
  
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
#' \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' 
#'  
#' @return Data Frame including one row per case_list and five columns 
#' 
#' @examples 
#' blca_case_lists <- get_case_lists("blca_tcga")
#' blca_case_lists
#' 
#' @export
get_case_lists <- function(csid = "blca_tcga") {

  stopifnot(!is.null(csid), length(csid) == 1, 
            !is.na(csid), is.character(csid))
  
  curWarn <- options()$warn
  options(warn = -1)

  my_url <- "http://www.cbioportal.org/webservice.do?cmd=getCaseLists"
  my_url <- paste(my_url, "&cancer_study_id=", 
                  tolower(as.character(csid)), sep = "")
  
  OUT <- tryCatch({as.data.frame(basic_tcga_query(my_url))}, 
                  error = function(e) { data.frame() })

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
#' \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' 
#'
#' 
#' @return data.frame including one row per patient/case/sample
#' 
#' @examples 
#' clinic_data <- get_clinical_data("blca_tcga_all")
#' clinic_data
#' 
#' @export
get_clinical_data <- function(case_id = 'blca_tcga_all') {

  stopifnot(!is.null(case_id), length(case_id) == 1, 
            !is.na(case_id), is.character(case_id), 
            nchar(case_id) > 0)
  
  curWarn <- options()$warn
  options(warn = -1)
  
  my_url <- "http://www.cbioportal.org/webservice.do?cmd=getClinicalData"
  my_url <- paste(my_url, "&case_set_id=", case_id, sep = "")
  OUT <- tryCatch({as.data.frame(basic_tcga_query(my_url))}, 
                  error = function(e) { data.frame() })

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
#' \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' 
#'  
#' 
#' @return data Frame inluding one row per mutation
#' 
#' @examples 
#' my_case_id <- 'blca_tcga_sequence'
#' my_gprofile_id <- 'blca_tcga_mutations'
#' my_gene <- 'TP53'
#' tp53_mutats <- get_ext_mutation(my_case_id, my_gprofile_id, my_gene)
#' if(ncol(tp53_mutats) >= 6 & nrow(tp53_mutats) >= 10){
#'   tp53_mutats[1:10,1:6]
#' }
#' 
#' @export
get_ext_mutation <- function(case_id = 'blca_tcga_sequence', 
                             gprofile_id = 'blca_tcga_mutations', 
                             glist = c('TP53', 'PTEN')){

  stopifnot(!is.null(case_id), length(case_id) == 1, 
            !is.na(case_id), is.character(case_id), nchar(case_id) > 0, 
            !is.null(gprofile_id), length(gprofile_id) == 1, 
            !is.na(gprofile_id), is.character(gprofile_id), 
            nchar(gprofile_id) > 0, !is.null(glist), length(glist) > 0) 
  
  curWarn <- options()$warn
  options(warn = -1)

  my_url <- "http://www.cbioportal.org/webservice.do?cmd=getMutationData"
  my_url <- paste(my_url, "&case_set_id=", case_id, sep = "")
  my_url <- paste(my_url, "&genetic_profile_id=", 
                  paste(gprofile_id, collapse = "+"), sep = "")
  my_url <- paste(my_url, "&gene_list=", 
                  paste(glist, collapse = "+"), sep = "")
  
  OUT <- tryCatch({basic_tcga_query(my_url)}, 
                  error = function(e) { data.frame() })

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
#' \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' 
#' 
#' @return 
#' data.frame including one row per genetic profile and six columns
#' 
#' @examples 
#' get_genetic_profiles("blca_tcga")
#' 
#' @export
get_genetic_profiles <- function(csid = NULL){
  
  stopifnot(!is.null(csid), length(csid) == 1, 
            !is.na(csid), is.character(csid)) 
  
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
#' Retrieve Data corresponding to a Genetic Profile of interest from a study of interest.
#' This function is the workhorse of the TCGAretriever package and can be used to fetch data 
#' concerning several genes at once. For larger queries, the use of the fetch_all_tcgadata() 
#' function is mandatory.
#' 
#' @param case_id String corresponding to the Identifier of a list of cases 
#' @param gprofile_id String corresponding to the Identifier of a genetic Profile of interest 
#' @param glist Character vector including one or more gene identifiers (ENTREZID or OFFICIAL_SYMOL)
#' @param force_numeric logical. Shall numeric data be coerced to numeric?
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' 
#'  
#' @return 
#' data.frame with one row per gene and one column per case/sample
#' 
#' @examples 
#' get_profile_data("blca_tcga_all", "blca_tcga_mutations", c("TP53", "E2F1"))
#' 
#' @export
get_profile_data <- function(case_id = 'blca_tcga_all', 
                             gprofile_id = 'blca_tcga_mutations', 
                             glist = c("TP53", "E2F1"), 
                             force_numeric = FALSE ){

  # Check Args, part 1
  stopifnot(!is.null(case_id), length(case_id) == 1, !is.na(case_id), 
            !is.null(gprofile_id), length(gprofile_id) == 1, 
            !is.na(gprofile_id), !is.null(glist))

  # Check Args, part 1
  glist <- unique(glist)
  glist <- glist[!is.na(glist)]
  stopifnot(length(glist) > 0, length(glist) <= 500)  
  
  curWarn <- options()$warn
  options(warn = -1)

  # Build Query
  my_url <- "http://www.cbioportal.org/webservice.do?cmd=getProfileData"
  my_url <- paste(my_url, "&case_set_id=", case_id, sep = "")
  my_url <- paste(my_url, "&genetic_profile_id=", paste(gprofile_id, collapse = "+"), sep = "")
  my_url <- paste(my_url, "&gene_list=", paste(glist, collapse = "+"), sep = "")
  result <- tryCatch({as.data.frame(basic_tcga_query(my_url))}, 
                     error = function(e) { data.frame() })
  
  if (force_numeric){
    if (ncol(result) > 0 && nrow(result) > 0) {
      for (i in seq(1, ncol(result), by = 1)) {
        v_i <- (!is.na(suppressWarnings(as.numeric(result[, i])))) | 
          result[, i] == 'NaN'
        rat_i <- sum(v_i, na.rm = TRUE) / nrow(result)
        if (rat_i >= 0.7) {
          result[, i] <- suppressWarnings(as.numeric(result[, i]))
        }
      }
    }
  }

  options(warn = curWarn)
  return(result)
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
#' \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#' 
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
  
  OUT <- tryCatch({
    if(is.numeric(num_vector) & is.numeric(groups) & 
       length(num_vector) > 5 & groups[1] > 1){
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
    result
  }, error = function(e) { NULL })
  
  options(warn = curWarn)
  return(OUT)
}



#' Retrieve Genomic and Clinical Data from CBioPortal
#'
#' @description 
#' The Cancer Genome Atlas (TCGA) is a scientific and medical program 
#' aimed at improving our understanding of Cancer Biology. 
#' Several TCGA Datasets are available online, and some are hosted
#' on cBioPortal, which is an open-access, open-source resource for 
#' interactive exploration of multidimensional cancer genomics data sets. 
#' TCGAretriever helps accessing and downloading TCGA data via the 
#' cBioPortal API. Features of TCGAretriever are: 
#' \itemize{
#'   \item it is simple and reliable
#'   \item it is tailored for downloading large volumes of data
#' }
#' 
#' @author Damiano Fantini, \email{damiano.fantini@@gmail.com}
#' @references 
#' \itemize{
#'   \item \url{https://www.data-pulse.com/dev_site/TCGAretriever/}
#'   \item \url{http://www.cbioportal.org/} 
#'   \item \url{http://cancergenome.nih.gov/abouttcga/}
#' }
#'  
#' 
#' @keywords internal
"_PACKAGE"


