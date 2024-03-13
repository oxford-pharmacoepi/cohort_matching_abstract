name <- "matched"
targetCohortName <- NULL




generateMatchedCohortSet <- function(cdm,
                                     name,
                                     targetCohortName,
                                     targetCohortId = NULL,
                                     matchSex = TRUE,
                                     matchYearOfBirth = TRUE,
                                     ratio = 1){
  # validate initial input
  validateInput(
    cdm = cdm, name = name, targetCohortName = targetCohortName,
    targetCohortId = targetCohortId, matchSex = matchSex,
    matchYearOfBirth = matchYearOfBirth, ratio = ratio
  )
  
  # table prefix
  #tablePrefix <- randomPrefix()
  
  # get the number of cohorts
  n <- getNumberOfCohorts(cdm, targetCohortName)
  
  if (n == 0) {
    cdm[[name]] <- cdm[[targetCohortName]] %>%
      dplyr::compute(name = name, temporary = FALSE) |>
      omopgenerics::newCohortTable()
    
  } else {
    # get target cohort id
    targetCohortId <- getTargetCohortId(cdm, targetCohortId, targetCohortName)
    
    # Create the cohort name with cases and controls of the targetCohortId
    cdm <- getNewCohort(cdm, name, targetCohortName, targetCohortId, n)
    
    # Exclude cases from controls
    cdm <- excludeCases(cdm, name, targetCohortId, n)
    
    # get matched tables
    matchCols <- getMatchCols(matchSex, matchYearOfBirth)
    
    if(!is.null(matchCols)){
      # Exclude individuals without any match
      cdm <- excludeNoMatchedIndividuals(cdm, name, matchCols, n)
      
      # Match as ratio was infinite
      cdm <- infiniteMatching(cdm, name, targetCohortId)
      
      # Delete controls that are not in observation
      cdm <- checkObservationPeriod(cdm, name, targetCohortId, n)
      
      # Check ratio
      cdm <- checkRatio(cdm, name, ratio, targetCohortId, n)
      
      # Check cohort set ref
      cdm <- checkCohortSetRef(cdm, name, targetCohortName, matchSex, matchYearOfBirth, targetCohortId, n)
      
      # Rename cohort definition ids
      cdm <- renameCohortDefinitionIds(cdm, name)
      
    } else {
      # TO DO
    }
  }
  # Return
  return(cdm)
}

#' @noRd
validateInput <- function(cdm,
                          name,
                          targetCohortName,
                          targetCohortId,
                          matchSex,
                          matchYearOfBirth,
                          ratio) {
  errorMessage <- checkmate::makeAssertCollection()
  # Check cdm class
  data_check   <- any("cdm_reference" == class(cdm))
  checkmate::assertTRUE(data_check, add = errorMessage)
  if(!isTRUE(data_check)){
    errorMessage$push(glue::glue("- cdm input must be a cdm object"))
  }
  # Check if targetCohortName is a character
  targetCohortName_format_check <- any(class(targetCohortName) %in% c("character"))
  checkmate::assertTRUE(targetCohortName_format_check, add = errorMessage)
  if(!isTRUE(targetCohortName_format_check)){
    errorMessage$push(glue::glue("- targetCohortName input must be a string"))
  }
  # Check if targetCohortName length
  targetCohortName_length_check <- length(targetCohortName) == 1
  checkmate::assertTRUE(  targetCohortName_length_check, add = errorMessage)
  if(!isTRUE(  targetCohortName_length_check)){
    errorMessage$push(glue::glue("- targetCohortName input must have length equal to 1"))
  }
  # Check if targetCohortName is within the cdm object
  targetCohortName_check <- targetCohortName %in% names(cdm)
  checkmate::assertTRUE(targetCohortName_check, add = errorMessage)
  if(!isTRUE(targetCohortName_check)){
    errorMessage$push(glue::glue("- cdm input has not table named {targetCohortName}"))
  }
  # Check if observation period is within the cdm object
  observation_period_check <- "observation_period" %in% names(cdm)
  checkmate::assertTRUE(observation_period_check , add = errorMessage)
  if(!isTRUE(observation_period_check)){
    errorMessage$push(glue::glue("- cdm input has not table named 'observation_period'"))
  }
  # Check if targetCohortId is a numeric value
  if(!is.null(targetCohortId)){
    targetCohortId_format_check <- any(class(targetCohortId) %in% c("numeric","double","integer"))
    checkmate::assertTRUE(targetCohortId_format_check, add = errorMessage)
    if(!isTRUE(targetCohortId_format_check)){
      errorMessage$push(glue::glue("- targetCohortId input must be numeric"))
    }
  }
  # Check if targetCohortId is in the cohort_definition_id
  if(!is.null(targetCohortId)){
    rows <- cdm[[targetCohortName]] %>% dplyr::filter(.data$cohort_definition_id %in% targetCohortId) %>% dplyr::tally() %>% dplyr::pull()
    targetCohortId_check <- rows != 0
    checkmate::assertTRUE(targetCohortId_check, add = errorMessage)
    if(!isTRUE(targetCohortId_check)){
      errorMessage$push(glue::glue("- {name} table does not containg '{targetCohortId}' as a cohort_definition_id"))
    }
  }
  # Check if ratio is > 0
  ratio_check <- ratio > 0
  checkmate::assertTRUE(ratio_check, add = errorMessage)
  if(!isTRUE(ratio_check)){
    errorMessage$push(glue::glue("- ratio parameter must be > 0 "))
  }
  
  checkmate::reportAssertions(collection = errorMessage)
  return(invisible(TRUE))
}









renameCohortDefinitionIds <- function(cdm, name){
  new_cohort_set <- cdm[[name]] %>%
    omopgenerics::settings() %>%
    dplyr::mutate(cohort_definition_id_new = .data$target_cohort_id) %>%
    dplyr::arrange(.data$cohort_definition_id_new) %>%
    dplyr::mutate(cohort_definition_id_new = dplyr::row_number())
  
  new_cohort_attrition <- cdm[[name]] %>%
    omopgenerics::attrition() %>%
    dplyr::inner_join(
      new_cohort_set %>% dplyr::select("cohort_definition_id","cohort_definition_id_new"),
      by = "cohort_definition_id"
    ) %>%
    dplyr::select(-"cohort_definition_id") %>%
    dplyr::rename("cohort_definition_id" = "cohort_definition_id_new") %>%
    dplyr::relocate("cohort_definition_id")
  
  new_cohort_count <- cdm[[name]] %>%
    omopgenerics::cohortCount() %>%
    dplyr::inner_join(
      new_cohort_set %>% dplyr::select("cohort_definition_id","cohort_definition_id_new"),
      by = "cohort_definition_id"
    ) %>%
    dplyr::select(-"cohort_definition_id") %>%
    dplyr::rename("cohort_definition_id" = "cohort_definition_id_new") %>%
    dplyr::relocate("cohort_definition_id")
  
  new_cohort <- cdm[[name]] %>%
    dplyr::inner_join(
      new_cohort_set %>% dplyr::select("cohort_definition_id","cohort_definition_id_new"),
      by = "cohort_definition_id",
      copy = TRUE
    ) %>%
    dplyr::select(-"cohort_definition_id") %>%
    dplyr::rename("cohort_definition_id" = "cohort_definition_id_new") %>%
    dplyr::relocate("cohort_definition_id") %>%
    dplyr::compute(name = name, temporary = FALSE)
  
  new_cohort_set <- new_cohort_set %>%
    dplyr::select(-"cohort_definition_id") %>%
    dplyr::rename("cohort_definition_id" = "cohort_definition_id_new") %>%
    dplyr::relocate("cohort_definition_id")
  
  cdm[[name]] <- omopgenerics::newCohortTable(
    table = new_cohort,
    cohortAttritionRef =  new_cohort_attrition,
    cohortSetRef = new_cohort_set
  )
  
  return(cdm)
}