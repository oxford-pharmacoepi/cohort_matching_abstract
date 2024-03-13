generateMatchedCohortSet_mah <- function(cdm,
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
    }
  }
  # Return
  return(cdm)
}


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

getNumberOfCohorts <- function(cdm, targetCohortName){
  tic(msg = "getNumberOfCohorts")
  # Read number of cohorts
  n <- cdm[[targetCohortName]] %>%
    dplyr::summarise(v = max(.data$cohort_definition_id, na.rm = TRUE)) %>%
    dplyr::pull("v") # number of different cohorts
  
  if(is.na(n)){# Empty table, number of cohorts is 0
    n <- 0
  }
  toc(log = TRUE)
  return(n)
}

getTargetCohortId <- function(cdm, targetCohortId, targetCohortName){
  tic(msg = "getTargetCohortId")
  if(is.null(targetCohortId)){
    targetCohortId <-cdm[[targetCohortName]] %>%
      omopgenerics::settings() %>%
      dplyr::arrange(.data$cohort_definition_id) %>%
      dplyr::pull("cohort_definition_id")
  }
  toc(log = TRUE) 
  return(targetCohortId)

}

setInitialControlAttriton <- function(cdm, ids) {
  tic(msg = "setInitialControlAttrition")
  num_records <- cdm[["person"]] %>%
    dplyr::tally() %>%
    dplyr::pull(.data$n)
  num_subjects <- cdm[["person"]] %>%
    dplyr::distinct(.data$person_id) %>%
    dplyr::tally() %>%
    dplyr::pull(.data$n)
  
  toc(log = TRUE)
  return(
    dplyr::tibble(
      cohort_definition_id = ids,
      number_records = num_records,
      number_subjects = num_subjects,
      reason_id = 1,
      reason = "Subjects in the database",
      excluded_records = 0,
      excluded_subjects = 0
    )
  )

}

getNewCohort <- function(cdm, name, targetCohortName, targetCohortId, n){
  tic(msg = "getNewCohort, part 1")
  # Create controls cohort
  temp_name <- "temp_ctr_ids"
  cdm <- omopgenerics::insertTable(
    cdm = cdm,
    name = temp_name,
    table = dplyr::tibble(cohort_definition_id = targetCohortId+n)
  )
  
  controls <- cdm[[temp_name]] %>%
    dplyr::cross_join(cdm[["person"]] %>%
                        dplyr::select("subject_id" = "person_id")) %>%
    dplyr::compute()
  cdm <- omopgenerics::dropTable(cdm, temp_name)
  toc(log = TRUE)
  
  # Create table with controls + cases (all cases existing in the cohort, without considering the targetCohortId)
  tic(msg = "getNewCohort, part 2")
  all <- controls %>%
    dplyr::inner_join(
      cdm$observation_period %>%
        dplyr::group_by(.data$person_id) %>%
        dplyr::filter(.data$observation_period_start_date == min(.data$observation_period_start_date, na.rm = TRUE)) %>%
        dplyr::filter(.data$observation_period_end_date == max(.data$observation_period_end_date, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::select(
          "subject_id" = "person_id",
          "cohort_start_date" = "observation_period_start_date",
          "cohort_end_date" = "observation_period_end_date"
        ),
      by = "subject_id"
    ) %>%
    dplyr::union_all(
      cdm[[targetCohortName]] %>%
        dplyr::filter(.data$cohort_definition_id %in% .env$targetCohortId)
    ) %>%
    dplyr::compute(name = name, temporary = FALSE)
  toc(log = TRUE)
  
  # settings
  tic(msg = "getNewCohort, part 3")
  cohort_set_ref <- cdm[[targetCohortName]] %>%
    omopgenerics::settings() %>%
    dplyr::filter(.data$cohort_definition_id %in% .env$targetCohortId) %>%
    dplyr::slice(rep(1:dplyr::n(), times = 2)) %>%
    dplyr::group_by(.data$cohort_definition_id) %>%
    dplyr::mutate(
      cohort_name = dplyr::if_else(dplyr::row_number() == 2, paste0(.data$cohort_name,"_matched"), .data$cohort_name),
      cohort_definition_id = dplyr::if_else(dplyr::row_number() == 2, .data$cohort_definition_id+.env$n, .data$cohort_definition_id)
    ) %>%
    dplyr::ungroup()
  toc(log = TRUE)
  
  # attrition
  tic(msg = "getNewCohort, part 4")
  cohort_attrition <- cdm[[targetCohortName]] %>%
    omopgenerics::attrition() %>%
    dplyr::filter(.data$cohort_definition_id %in% .env$targetCohortId) %>%
    dplyr::union_all(setInitialControlAttriton(cdm, targetCohortId+n))
  toc(log = TRUE)
  tic(msg = "getNewCohort, part 5")
  
  cdm[[name]] <- omopgenerics::newCohortTable(
    table = all,
    cohortAttritionRef = cohort_attrition %>% dplyr::as_tibble(),
    cohortSetRef = cohort_set_ref
  )
  toc(log = TRUE)
  return(cdm)

}

excludeCases <- function(cdm, name, targetCohortId, n){
  tic(msg = "excludeCases, part 1")
  # For each target cohort id
  for(targetCohortId_i in targetCohortId){
    # Controls
    controls <- cdm[[name]] %>%
      dplyr::filter(.data$cohort_definition_id == targetCohortId_i+.env$n) %>%
      dplyr::anti_join(
        # Cases
        cdm[[name]] %>%
          dplyr::select("subject_id","cohort_definition_id") %>%
          dplyr::filter(.data$cohort_definition_id == targetCohortId_i) %>%
          dplyr::mutate(cohort_definition_id = targetCohortId_i + .env$n),
        by = c("subject_id", "cohort_definition_id")
      ) %>%
      dplyr::compute()
    
    cdm[[name]] <- cdm[[name]] %>%
      # Delete the controls
      dplyr::filter(.data$cohort_definition_id != targetCohortId_i + .env$n) %>%
      # Add the new controls set
      dplyr::union_all(controls) %>%
      dplyr::compute(name = name, temporary = FALSE)
  }
  toc(log = TRUE)
  
  # Record attrition
  tic(msg = "excludeCases, part 2")
  cdm[[name]] <- cdm[[name]] %>%
    CDMConnector::record_cohort_attrition("Exclude cases",
                                          cohortId = c(targetCohortId+n))%>%
    dplyr::compute(name = name, temporary = FALSE)
  toc(log = TRUE) 
  return(cdm)

}

getMatchCols <- function(matchSex, matchYearOfBirth){
  tic(msg = "getMatchCols")
  # Obtain matched columns
  matchCols <- c()
  if(matchSex){
    matchCols <- append(matchCols, "gender_concept_id")
  }
  if(matchYearOfBirth){
    matchCols <- append(matchCols, "year_of_birth")
  }  
  toc(log = TRUE)
  return(matchCols)

}

excludeNoMatchedIndividuals <- function(cdm, name, matchCols, n){
  tic(msg = "excludeNotMatchedIndividuals, part 1")
  cdm[[name]] <- cdm[[name]] %>%
    # Append matchcols
    dplyr::left_join(
      cdm[["person"]] %>%
        dplyr::select("subject_id" = "person_id", dplyr::all_of(matchCols)),
      by = c("subject_id")
    ) %>%
    dplyr::compute(name = name, temporary = FALSE)
  toc(log = TRUE)
  
  # Create column group id
  tic(msg = "excludeNotMatchedIndividuals, part 2")
  cdm[[name]] <- cdm[[name]] %>%
    dplyr::inner_join(
      cdm[[name]] %>%
        dplyr::select(dplyr::all_of(matchCols)) %>%
        dplyr::distinct() %>%
        dplyr::mutate(group_id = dplyr::row_number()),
      by = c(matchCols)
    ) %>%
    dplyr::select(-dplyr::all_of(matchCols)) %>%
    # Create target definition id column
    dplyr::mutate(target_definition_id =
                    dplyr::if_else(
                      .data$cohort_definition_id <= .env$n,
                      .data$cohort_definition_id,
                      .data$cohort_definition_id - .env$n
                    )) %>%
    dplyr::compute(name = name, temporary = FALSE)
  toc(log = TRUE)
  
  # Exclude individuals that do not have any match
  tic(msg = "excludeNotMatchedIndividuals, part 3")
  cdm[[name]] <- cdm[[name]] %>%
    dplyr::inner_join(
      cdm[[name]] %>%
        dplyr::mutate(
          "cohort_definition_id" = dplyr::if_else(
            .data$target_definition_id == .data$cohort_definition_id,
            .data$cohort_definition_id + .env$n,
            .data$cohort_definition_id - .env$n
          )
        ) %>%
        dplyr::select("cohort_definition_id", "target_definition_id", "group_id") %>%
        dplyr::distinct(),
      by = c("target_definition_id", "group_id", "cohort_definition_id")
    ) %>%
    dplyr::compute(name = name, temporary = FALSE) %>%
    CDMConnector::record_cohort_attrition("Exclude individuals that do not have any match")
  toc(log = TRUE)
  return(cdm)
}

infiniteMatching <- function(cdm, name, targetCohortId){
  tic(msg = "infiniteMatching, part 1")
  # Create pair id to perform a random match
  cdm[[name]] <- cdm[[name]] %>%
    dplyr::mutate(id = dbplyr::sql("random()")) %>%
    dplyr::group_by(.data$cohort_definition_id, .data$group_id) %>%
    dbplyr::window_order(.data$id) %>%
    dplyr::mutate(pair_id = dplyr::row_number()) %>%
    dplyr::select(-"id") %>%
    dplyr::ungroup() %>%
    dplyr::compute(name = name, temporary = FALSE)
  toc(log = TRUE)
  
  tic(msg = "infiniteMatching, part 2")
  cdm[[name]] <- cdm[[name]] %>%
    dplyr::inner_join(
      # Calculate the maximum number of cases per group
      cdm[[name]] %>%
        dplyr::filter(.data$cohort_definition_id %in% .env$targetCohortId) %>%
        dplyr::group_by(.data$cohort_definition_id, .data$group_id) %>%
        dplyr::mutate(max_cases = max(.data$pair_id, na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::select("group_id", "target_definition_id", "max_cases") %>%
        dplyr::distinct(),
      by = c("group_id", "target_definition_id")
    ) %>%
    # Calculate the maximum ratio per group
    dplyr::mutate(id = (.data$pair_id-1) %% .data$max_cases + 1) %>%
    dplyr::mutate(pair_id = .data$id) %>%
    dplyr::select(-"max_cases", -"id") %>%
    dplyr::compute(name = name, temporary = FALSE)
  toc(log = TRUE)
  
  tic(msg = "infiniteMatching, part 3")
  # Perform random matches with ratio 1:Inf
  cdm[[name]] <- cdm[[name]] %>%
    dplyr::select(-"cohort_start_date", -"cohort_end_date") %>%
    dplyr::inner_join(
      # Cohort start date and end date of cases
      cdm[[name]] %>%
        dplyr::filter(.data$cohort_definition_id %in% targetCohortId) %>%
        dplyr::select("pair_id", "group_id", "target_definition_id", "cohort_start_date", "cohort_end_date"),
      by = c("pair_id", "group_id", "target_definition_id")
    ) %>%
    dplyr::distinct() %>%
    dplyr::compute(name = name, temporary = FALSE)
  toc(log = TRUE)
  return(cdm)
}

checkObservationPeriod <- function(cdm, name, targetCohortId, n){
  tic(msg = "checkObservationPeriod")
  cdm[[name]] <- cdm[[name]] %>%
    # Add future observation
    PatientProfiles::addFutureObservation() %>%
    # Remove those with no future observation
    dplyr::filter(!is.na(.data$future_observation)) %>%
    # Cohort end date is, if the cohort is the one targeted, cohort end date, if not future observation - cohort start date
    dplyr::mutate(cohort_end_date = dplyr::if_else(
      .data$cohort_definition_id %in% .env$targetCohortId,
      .data$cohort_end_date,
      !!CDMConnector::dateadd("cohort_start_date", "future_observation")
    )) %>%
    # remove future observation column
    dplyr::select(-"future_observation") |>
    dplyr::compute(name = name, temporary = FALSE)
  
  cdm[[name]] <- cdm[[name]] |>
    dplyr::inner_join(
      cdm[[name]] |>
        dplyr::mutate(cohort_definition_id = dplyr::if_else(cohort_definition_id == target_definition_id,
                                              cohort_definition_id + n,
                                              cohort_definition_id - n)) |>
        dplyr::select("cohort_definition_id", "target_definition_id", "group_id", "pair_id"),
      by = c("cohort_definition_id", "target_definition_id", "group_id", "pair_id")
    ) |>
    # group by target definition id, group_id, and pair_id
    # dplyr::group_by(.data$target_definition_id, .data$group_id, .data$pair_id) %>%
    # dplyr::filter(dplyr::n() > 1) %>%
    # dplyr::ungroup() %>%
    dplyr::compute(name = name, temporary = FALSE) %>%
    CDMConnector::record_cohort_attrition("Exclude individuals that are not in observation", cohortId = targetCohortId + n) %>%
    CDMConnector::record_cohort_attrition("Exclude individuals that their only pair is not in observation", cohortId = targetCohortId)
  toc(log = TRUE)
   return(cdm)
}

checkRatio <- function(cdm, name, ratio, targetCohortId, n){
  tic(msg = "checkRatio")
  if (ratio == Inf) {
    cdm[[name]] <- cdm[[name]] %>%
      dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") %>%
      dplyr::compute(name = name, temporary = FALSE)
  } else {
    cdm[[name]] <- cdm[[name]] %>%
      dplyr::group_by(.data$pair_id, .data$group_id, .data$target_definition_id) %>%
      dbplyr::window_order(.data$cohort_definition_id) %>%
      dplyr::filter(dplyr::row_number() <= .env$ratio+1) %>%
      dplyr::ungroup() %>%
      dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") %>%
      dplyr::compute(name = name, temporary = FALSE) %>%
      CDMConnector::record_cohort_attrition("Exclude individuals that do not fulfil the ratio", cohortId = targetCohortId+n)
  }
  toc(log = TRUE)
  
  return(cdm)
}


checkCohortSetRef <- function(cdm, name, targetCohortName, matchSex, matchYearOfBirth, targetCohortId, n){
  tic(msg = "checkCohortSetRef")
  cohort_set_ref <- cdm[[name]] %>%
    omopgenerics::settings() %>%
    dplyr::mutate(target_cohort_name  = .env$targetCohortName) %>%
    dplyr::mutate(match_sex           = .env$matchSex) %>%
    dplyr::mutate(match_year_of_birth = .env$matchYearOfBirth) %>%
    dplyr::mutate(match_status        = dplyr::if_else(.data$cohort_definition_id %in% .env$targetCohortId, "target", "matched")) %>%
    dplyr::mutate(target_cohort_id    = dplyr::if_else(.data$cohort_definition_id %in% .env$targetCohortId, .data$cohort_definition_id, .data$cohort_definition_id-n))
  
  cdm[[name]] <- omopgenerics::newCohortTable(
    table = cdm[[name]],
    cohortSetRef = cohort_set_ref
  )
  toc(log = TRUE)
  
  return(cdm)
}

renameCohortDefinitionIds <- function(cdm, name){
  tic(msg = "renameCohortDefinitionIds, part 1")
  new_cohort_set <- cdm[[name]] %>%
    omopgenerics::settings() %>%
    dplyr::mutate(cohort_definition_id_new = .data$target_cohort_id) %>%
    dplyr::arrange(.data$cohort_definition_id_new) %>%
    dplyr::mutate(cohort_definition_id_new = dplyr::row_number())
  toc(log = TRUE)
  
  tic(msg = "renameCohortDefinitionIds, part 2")
  new_cohort_attrition <- cdm[[name]] %>%
    omopgenerics::attrition() %>%
    dplyr::inner_join(
      new_cohort_set %>% dplyr::select("cohort_definition_id","cohort_definition_id_new"),
      by = "cohort_definition_id"
    ) %>%
    dplyr::select(-"cohort_definition_id") %>%
    dplyr::rename("cohort_definition_id" = "cohort_definition_id_new") %>%
    dplyr::relocate("cohort_definition_id")
  toc(log = TRUE)
  
  tic(msg = "renameCohortDefinitionIds, part 3")
  new_cohort_count <- cdm[[name]] %>%
    omopgenerics::cohortCount() %>%
    dplyr::inner_join(
      new_cohort_set %>% dplyr::select("cohort_definition_id","cohort_definition_id_new"),
      by = "cohort_definition_id"
    ) %>%
    dplyr::select(-"cohort_definition_id") %>%
    dplyr::rename("cohort_definition_id" = "cohort_definition_id_new") %>%
    dplyr::relocate("cohort_definition_id")
  toc(log = TRUE)
  
  tic(msg = "renameCohortDefinitionIds, part 4")
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
  toc(log = TRUE)
  
  tic(msg = "renameCohortDefinitionIds, part 5")
  new_cohort_set <- new_cohort_set %>%
    dplyr::select(-"cohort_definition_id") %>%
    dplyr::rename("cohort_definition_id" = "cohort_definition_id_new") %>%
    dplyr::relocate("cohort_definition_id")
  toc(log = TRUE)
  
  tic(msg = "renameCohortDefinitionIds, part 5")
  cdm[[name]] <- omopgenerics::newCohortTable(
    table = new_cohort,
    cohortAttritionRef =  new_cohort_attrition,
    cohortSetRef = new_cohort_set
  )
  toc(log = TRUE)
  
  return(cdm)
}