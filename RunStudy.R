library("dplyr")
library("DrugUtilisation")
library("duckdb")
library("PatientProfiles")
library("IncidencePrevalence")
library("devtools")
library("tictoc")
library("tidyr")

denominator_concept_id <- list("headache" = 378253)
denominator_cohort  <- "denominator"
factor <- 100000

# Create the denominator cohort
cdm <- generateConceptCohortSet(cdm,
                                conceptSet = denominator_concept_id,
                                name       = denominator_cohort,
                                overwrite  = TRUE) 
cdm[[denominator_cohort]] <- cdm[[denominator_cohort]] |> newCohortTable()

# Large scale characterization -------------------------------------------------
info(logger, "START LARGE SCALE CHARACTERISATION - ORIGINAL COHORT")
tic(msg = "- Large scale characterisation of the cases")

# Large scale characterization of the cases (conditions, drug_exposures) 
info(logger, "- Large scale characteristics of the cases")
largeScaleCharacteristics_analysis1 <- cdm[[denominator_cohort]] |>
  summariseLargeScaleCharacteristics(
    strata = list(),
    window = list(c(-Inf, -366), c(-365, -31), c(-30,-1)),
    eventInWindow   = "condition_occurrence",
    episodeInWindow = NULL,
    indexDate  = "cohort_start_date",
    censorDate = NULL,
    includeSource    = TRUE,
    minimumFrequency = 0.005,
    excludedCodes = NULL
  )
x <- toc(log = TRUE)
info(logger, x$callback_msg)

# Calculate counts using achilles ----------------------------------------------
info(logger, "START LARGE SCALE CHARACTERISATION USING ACHILLES")
tic(msg = "- Large scale characterisation using achilles")

# Select concepts
info(logger, "- Select conditions")
concepts <- largeScaleCharacteristics_analysis1 |>
  # Concepts in largeScaleCharacterisation
  select("concept_name" = "variable_name", "additional_level") |>
  distinct() |>
  separate(additional_level, c("table_name", "type", "analysis", "concept_id"), sep = " and ") |>
  mutate(concept_id = as.numeric(concept_id)) |>
  distinct() |>
  inner_join(
    # Counts of the conditions in Achilles tables
    cdm$achilles_results |>
      filter(analysis_id == 401) |>
      select("concept_id" = "stratum_1","counts" = "count_value") |>
      mutate(concept_id = as.numeric(concept_id)),
    by = "concept_id",
    copy = TRUE
  ) |>
  select("concept_name", "table_name", "concept_id", "counts")

# Calculate estimates
concepts |>
  mutate(n_people      = 17216081,
         n_denominator = cdm[[denominator_cohort]] |> tally() |> pull() |> as.numeric(),
         factor        = factor) |>
  mutate(estimate = counts/n_people*n_denominator/factor*100)




# Number of person-years in all the database
info(logger, "- Calculate person-years")
person_years <- cdm$observation_period %>%
  mutate(days_in_observation = !!datediff(start = "observation_period_start_date",
                                          end   = "observation_period_end_date",
                                          interval = "day")) |>
  mutate(years_in_observation = days_in_observation/365) |>
  group_by(person_id) |>
  mutate(person_number = row_number()) |>
  filter(person_number == 1) |>
  mutate("person_days"  = days_in_observation) |>
  mutate("person_years" = years_in_observation) |>
  ungroup() |>
  summarise(person_years = sum(person_years),
            person_days  = sum(person_days),
            observation_period_start_date = min(observation_period_start_date),
            observation_period_end_date   = max(observation_period_end_date),
            n_persons    = sum(person_number)) |>
  collect()

# Calculate the uniform frequencies of the concepts
info(logger, "- Calculate the uniform frequencies for each concept")
counts |>
  mutate(person_years = person_years$person_years) |>
  mutate(cdm_name = "CPRD GOLD",
         variable_level = "N/A",
         estimate_name  = "incidence",
         estimate_type  = "numeric",
         estimate_value = count_value/person_years)
  
  


# Calculate incidence using cohort-matching ------------------------------------
#devtools::install_github("oxford-pharmacoepi/CohortConstructor")
library("CohortConstructor")
cdm <- cdm |>
  generateMatchedCohortSet(
  name = "matched",
  targetCohortName = denominator_cohort,
  targetCohortId   = NULL,
  matchSex = TRUE,
  matchYearOfBirth = TRUE,
  ratio = 1
)

cdm$matched <- cdm$matched |>
  mutate(subject_id = as.numeric(subject_id),
         cohort_end_date = as.Date(cohort_end_date))

analysis_table_3 <- estimateIncidence(
  cdm,
  denominatorTable = "matched",
  outcomeTable     = outcome_cohort,
  denominatorCohortId = NULL,
  outcomeCohortId     = NULL,
  interval = "overall",
  completeDatabaseIntervals = FALSE,
  outcomeWashout = Inf,
  repeatedEvents = FALSE,
  minCellCount   = 5,
  strata = list(),
  includeOverallStrata = TRUE,
  returnParticipants = FALSE
) |>
  select("analysis_id","n_persons","person_days","n_events","incidence_start_date",
         "incidence_end_date","person_years","incidence_100000_pys","incidence_100000_pys_95CI_lower",
         "incidence_100000_pys_95CI_upper") |>
  mutate(analysis_id = "3")

# Overall result ---------------------------------------------------------------
analysis_table_1 |>
  union_all(
    analysis_table_2
  ) |>
  union_all(
    analysis_table_3
  )

