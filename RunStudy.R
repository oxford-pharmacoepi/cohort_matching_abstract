library("dplyr")
library("DrugUtilisation")
library("duckdb")
library("PatientProfiles")
library("IncidencePrevalence")
library("devtools")
library("tictoc")
library("tidyr")
library("CDMConnector")
library("flextable")

source(here("Scripts/Functions.R"))

# denominator_concept_id <- list("headache" = 378253)
denominator_cohort  <- "hpv_cin23"
factor <- 1 #100000/17267137

# Create the denominator cohort
# cdm <- generateConceptCohortSet(cdm,
#                                 conceptSet = denominator_concept_id,
#                                 name       = denominator_cohort,
#                                 overwrite  = TRUE) 
# cdm[[denominator_cohort]] <- cdm[[denominator_cohort]] |> newCohortTable()

# Large scale characterization -------------------------------------------------
info(logger, "START LARGE SCALE CHARACTERISATION - ORIGINAL COHORT")
tic(msg = "- Large scale characterisation of the cases")

# Large scale characterization of the cases (conditions only) 
info(logger, "- Large scale characteristics of the cases")

lsc_analysis1 <- cdm[[denominator_cohort]] |>
  summariseLargeScaleCharacteristics(
    strata = list(),
    window = list( c(-365, -31), c(-30,-1)),
    eventInWindow   = "condition_occurrence",
    episodeInWindow = NULL,
    indexDate  = "cohort_start_date",
    censorDate = NULL,
    includeSource    = FALSE,
    minimumFrequency = 0.005,
    excludedCodes = NULL
  )

x1 <- toc(log = TRUE)
info(logger, x1$callback_msg)

# Calculate counts using achilles ----------------------------------------------
info(logger, "START LARGE SCALE CHARACTERISATION USING ACHILLES")
tic(msg = "- Large scale characterisation using achilles")

# Select concepts
info(logger, "- Select conditions")
concepts <- lsc_analysis1 |>
  # Concepts in largeScaleCharacterisation
  select("concept_name" = "variable_name", "additional_level") |>
  distinct() |>
  separate(additional_level, c("table_name", "type", "analysis", "concept_id"), sep = " and ") |>
  mutate(concept_id = as.numeric(concept_id)) |>
  left_join(
    # Counts of the conditions in Achilles tables
    cdm$achilles_results |>
      filter(analysis_id == 401) |>
      select("concept_id" = "stratum_1","stratum_2", "counts" = "count_value") |>
      mutate(concept_id = as.numeric(concept_id)),
    by = "concept_id",
    copy = TRUE
  ) |>
  select("concept_name", "table_name", "concept_id", "counts") |>
  mutate(counts = counts*factor)


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
            n_people     = sum(person_number)) |>
  collect()


# Calculate estimates and format table
lsc_analysis2 <- concepts |>
  mutate(n_people    = person_years$n_people,
         n_phenotype = cdm[[denominator_cohort]] |> tally() |> pull() |> as.numeric(),
         person_days = person_years$person_days,
         n_condition = counts) |>
  slice(rep(1:n(), each = 2)) |>
  group_by(concept_name) |>
  mutate(n = row_number()) |>
  mutate(variable_level = if_else(n == 1, "-365 to -31", "-30 to -1"))  |>
  mutate(estimate_name  = "percentage",
         estimate_value = if_else(n == 1, 335*n_condition/person_days*100, 30*n_condition/person_days*100))

x2 <- toc(log = TRUE)
info(logger, x2$callback_msg) 

# Calculate percentage using cohort-matching ------------------------------------
info(logger, "START MATCHING")
tic(msg = "- Large scale characterisation using matched cohorts")

#devtools::install_github("oxford-pharmacoepi/CohortConstructor", force = TRUE)
library("CohortConstructor")
cdm <- cdm |>
  generateMatchedCohortSet_mah(
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

lsc_analysis3 <- cdm[["matched"]] |>
  summariseLargeScaleCharacteristics(
    strata = list(c("cohort_definition_id")),
    window = list(c(-365, -31), c(-30,-1)),
    eventInWindow   = "condition_occurrence",
    episodeInWindow = NULL,
    indexDate  = "cohort_start_date",
    censorDate = NULL,
    includeSource    = FALSE,
    minimumFrequency = 0.005,
    excludedCodes = NULL
  )

x3 <- toc(log = TRUE)
info(logger, x3$callback_msg)


# Compute tables ---------------------------------------------------------------
a1 <- lsc_analysis1 |>
  filter(estimate_name == "percentage") |>
  select("variable_name", "variable_level","Original cohort (%)" = "estimate_value")
  
a2 <- lsc_analysis2 |>
  select("variable_name" = "concept_name", "variable_level","Achilles approximation (%)" = "estimate_value")

a3 <-  lsc_analysis3 |>
  filter(estimate_name == "percentage",
         strata_level  != "overall") |>
  select("variable_name", "variable_level","strata_level","estimate_value") |>
   mutate(strata_level = if_else(strata_level == 1, "Cohort (%)", "Cohort matched (%)")) |>
   pivot_wider(id_cols = c("variable_name","variable_level"),
               names_from  = "strata_level",
               values_from = "estimate_value")
    
tab <- a1 |>
  full_join(a2) |>
  full_join(a3) |>
  arrange(desc(`Original cohort (%)`)) |>
  group_by(variable_name) |>
  mutate(n = row_number()) |>
  mutate(n = if_else(n > 1, 0,n)) |>
  ungroup() |>
  mutate(n = cumsum(n)) |>
  group_by(variable_name) |>
  mutate(n = min(n)) |>
  arrange(n,variable_name,variable_level) |>
  select(-n) |>
  mutate(`Original cohort (%)` = round(as.numeric(`Original cohort (%)`),2),
         `Achilles approximation (%)` = round(as.numeric(`Achilles approximation (%)`),2),
         `Cohort (%)` = round(as.numeric(`Cohort (%)`),2),
         `Cohort matched (%)` = round(as.numeric(`Cohort matched (%)`),2))
  
tab <- rbind(tibble("variable_name"  = "Computational time",
                    "variable_level" = " ",
                    `Original cohort (%)` = gsub(" elap.*","",gsub(".*: ","",x1$callback_msg)),
                    `Achilles approximation (%)` = gsub(" elap.*","",gsub(".*: ","",x2$callback_msg)),
                    `Cohort (%)`          = gsub(" elap.*","",gsub(".*: ","",x3$callback_msg)),
                    `Cohort matched (%)`  = " "), 
             tab) |>
  rename("Conditions" = "variable_name") |>
  rename("Windows"    = "variable_level")

vec <- tab |>
  group_by(Conditions) |>
  mutate(n = row_number()) |> 
  mutate(n_max = max(n)) |>
  mutate(n = if_else(n == max(n), n, 0)) |> 
  ungroup() |>
  mutate(n = cumsum(n)) |>
  select(n) |>
  distinct() |>
  filter(n != 0) |>
  pull(n)


tab |>
  flextable() |>
  align(j = c(3,4,5,6), align = "center", part = "all") |>
  bold(i = 1, bold = TRUE, part = "header") |>
  hline(i = c(vec), part = "body") |>
  merge_at(i = 1, j = c(1,2), part = "body") |>
  merge_at(i = 1, j = c(5,6), part = "body") |>
  width(j = 1, width = 8, unit = "cm") |>
  width(j = 2, width = 3, unit = "cm") |>
  merge_v(j = ~Conditions) 
  
  
  
                






