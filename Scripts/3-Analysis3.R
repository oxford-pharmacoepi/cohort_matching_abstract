randomPrefix <- function(n = 5) {
  paste0(
    "temp_", paste0(sample(letters, 5, TRUE), collapse = ""), "_", collapse = ""
  )
}

getNumberOfCohorts <- function(cdm, targetCohortName){
  # Read number of cohorts
  n <- cdm[[targetCohortName]] %>%
    dplyr::summarise(v = max(.data$cohort_definition_id, na.rm = TRUE)) %>%
    dplyr::pull("v") # number of different cohorts
  
  if(is.na(n)){# Empty table, number of cohorts is 0
    n <- 0
  }
  return(n)
}

getTargetCohortId <- function(cdm, targetCohortId, targetCohortName){
  if(is.null(targetCohortId)){
    targetCohortId <-cdm[[targetCohortName]] %>%
      omopgenerics::settings() %>%
      dplyr::arrange(.data$cohort_definition_id) %>%
      dplyr::pull("cohort_definition_id")
  }
  
  return(targetCohortId)
}

setInitialControlAttriton <- function(cdm, ids) {
  num_records <- cdm[["person"]] %>%
    dplyr::tally() %>%
    dplyr::pull(.data$n)
  num_subjects <- cdm[["person"]] %>%
    dplyr::distinct(.data$person_id) %>%
    dplyr::tally() %>%
    dplyr::pull(.data$n)
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
  # Create controls cohort
  temp_name <- "temp_ctr_ids"
  cdm <- omopgenerics::insertTable(
    cdm   = cdm,
    name  = temp_name,
    table = dplyr::tibble(cohort_definition_id = targetCohortId+n)
  )
  controls <- cdm[[temp_name]] %>%
    dplyr::cross_join(cdm[["person"]] %>%
                        dplyr::select("subject_id" = "person_id")) %>%
    dplyr::compute()
  cdm <- omopgenerics::dropTable(cdm, temp_name)
  
  # Create table with controls + cases (all cases existing in the cohort, without considering the targetCohortId)
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
  
  # settings
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
  
  # attrition
  cohort_attrition <- cdm[[targetCohortName]] %>%
    omopgenerics::attrition() %>%
    dplyr::filter(.data$cohort_definition_id %in% .env$targetCohortId) %>%
    dplyr::union_all(setInitialControlAttriton(cdm, targetCohortId+n))
  
  cdm[[name]] <- omopgenerics::newCohortTable(
    table = all,
    cohortAttritionRef = cohort_attrition %>% dplyr::as_tibble(),
    cohortSetRef = cohort_set_ref
  )
  
  return(cdm)
}

excludeCases <- function(cdm, name, targetCohortId, n){
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
  
  # Record attrition
  cdm[[name]] <- cdm[[name]] %>%
    CDMConnector::record_cohort_attrition("Exclude cases",
                                          cohortId = c(targetCohortId+n))%>%
    dplyr::compute(name = name, temporary = FALSE)
  
  return(cdm)
}



getMatchCols <- function(matchSex, matchYearOfBirth){
  # Obtain matched columns
  matchCols <- c()
  if(matchSex){
    matchCols <- append(matchCols, "gender_concept_id")
  }
  if(matchYearOfBirth){
    matchCols <- append(matchCols, "year_of_birth")
  }
  return(matchCols)
}





excludeNoMatchedIndividuals <- function(cdm, name, matchCols, n){
  cdm[[name]] <- cdm[[name]] %>%
    # Append matchcols
    dplyr::left_join(
      cdm[["person"]] %>%
        dplyr::select("subject_id" = "person_id", dplyr::all_of(matchCols)),
      by = c("subject_id")
    ) %>%
    dplyr::compute(name = name, temporary = FALSE)
  
  # Create column group id
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
  
  # Exclude individuals that do not have any match
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
  
  return(cdm)
}

infiniteMatching <- function(cdm, name, targetCohortId){
  # Create pair id to perform a random match
  cdm[[name]] <- cdm[[name]] %>%
    dplyr::mutate(id = dbplyr::sql("random()")) %>%
    dplyr::group_by(.data$cohort_definition_id, .data$group_id) %>%
    dbplyr::window_order(.data$id) %>%
    dplyr::mutate(pair_id = dplyr::row_number()) %>%
    dplyr::select(-"id") %>%
    dplyr::ungroup() %>%
    dplyr::compute(name = name, temporary = FALSE)
  
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
  
  return(cdm)
}

checkObservationPeriod <- function(cdm, name, targetCohortId, n){
  cdm[[name]] <- cdm[[name]] %>%
    PatientProfiles::addFutureObservation() %>%
    dplyr::filter(!is.na(.data$future_observation)) %>%
    dplyr::mutate(cohort_end_date = dplyr::if_else(
      .data$cohort_definition_id %in% .env$targetCohortId,
      .data$cohort_end_date,
      !!CDMConnector::dateadd("cohort_start_date", "future_observation")
    )) %>%
    dplyr::select(-"future_observation") %>%
    dplyr::group_by(.data$target_definition_id, .data$group_id, .data$pair_id) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup() %>%
    dplyr::compute(name = name, temporary = FALSE) %>%
    CDMConnector::record_cohort_attrition("Exclude individuals that are not in observation", cohortId = targetCohortId + n) %>%
    CDMConnector::record_cohort_attrition("Exclude individuals that their only pair is not in observation", cohortId = targetCohortId)
  return(cdm)
}



checkRatio <- function(cdm, name, ratio, targetCohortId, n){
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
  
  return(cdm)
}



checkCohortSetRef <- function(cdm, name, targetCohortName, matchSex, matchYearOfBirth, targetCohortId, n){
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
  
  return(cdm)
}