# Analysis 1: Calculate incidence of concept_id_2 in cohort_name_1 -------------
# Flag cohorts:
cdm[[cohort_name_1]] <- cdm[[cohort_name_1]] |>
  addCohortIntersect(
    targetCohortTable = cohort_name_2,
    indexDate = "cohort_start_date",
    flag  = TRUE,
    date  = FALSE,
    count = FALSE,
    days  = FALSE,
    targetStartDate = "cohort_start_date",
    order     = "first",
    nameStyle = paste0(cohort_name_2),
    window = list(c(-Inf,0))
  ) |>
  rename("outcome" = cohort_name_2)

# Calculate proportion of people that have condition_id_2 in condition_id_1 cohort
table_analysis <- cdm[[cohort_name_1]] |>
  inner_join(
    cdm[["observation_period"]] |> 
      select("subject_id" = "person_id",
             "observation_period_start_date",
             "observation_period_end_date"),
    by = "subject_id"
  ) |>
  summarise("analysis_id" = 1,
            "observation_period_start_date" = min(observation_period_start_date, na.rm = TRUE),
            "observation_period_end_date"   = max(observation_period_end_date,   na.rm = TRUE),
            "number_of_cases"       = sum(outcome, na.rm = TRUE),
            "number_of_individuals" = n()) |>
  mutate("proportion_analysis" = number_of_cases/number_of_individuals*100)

# Calculate incidence 
incidence_analysis <- estimateIncidence(
  cdm,
  denominatorTable = cohort_name_1,
  outcomeTable     = cohort_name_2,
  denominatorCohortId = NULL,
  outcomeCohortId     = NULL,
  interval = "overall",
  completeDatabaseIntervals = TRUE, # Si faig incidence rate, això és important?
  outcomeWashout = Inf,
  repeatedEvents = FALSE,
  minCellCount = 5,
  strata = list(),
  includeOverallStrata = TRUE,
  temporary = FALSE,
  returnParticipants = FALSE
)
