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
factor <- 100000/17267137

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
  select("variable_name", "additional_level") |>
  distinct() |>
  separate(additional_level, c("table_name", "type", "analysis", "concept_id"), sep = " &&& ") |>
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
  select("variable_name", "table_name", "concept_id", "counts") |>
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
  group_by(variable_name) |>
  mutate(n = row_number()) |>
  mutate(variable_level = if_else(n == 1, "-365 to -31", "-30 to -1"))  |>
  mutate(estimate_name  = "percentage",
         estimate_value = if_else(n == 1, 335*n_condition/person_days*100, 30*n_condition/person_days*100),
         group_level = paste0(denominator_cohort, "_achilles"))

x2 <- toc(log = TRUE)
info(logger, x2$callback_msg) 

# Calculate percentage using cohort-matching ------------------------------------
info(logger, "START MATCHING")
tic(msg = "- Large scale characterisation using matched cohorts")

#devtools::install_github("oxford-pharmacoepi/CohortConstructor", force = TRUE)
library("CohortConstructor")
# Alternative matching
cdm[["matched_alternative"]] <- cdm[[denominator_cohort]] |>
  # Add sex and year of birth
  inner_join(
    cdm[["person"]] |>
      select("subject_id" = "person_id", "gender_concept_id", "year_of_birth"),
    by = "subject_id"
  ) |>
  rename("subject_id_reference" = subject_id) |>
  # Matching
  inner_join(
    # Remove people that is in the cohort
    cdm[["person"]] |>
      anti_join(
        cdm[[denominator_cohort]] |>
          select("person_id" = "subject_id"),
        by = "person_id"
      ) |>
      select("subject_id_comparator" = "person_id", "gender_concept_id", "year_of_birth"),
    by = c("gender_concept_id", "year_of_birth")
  ) |>
  group_by(subject_id_reference) |> 
  filter(row_number() == 1) |> 
  ungroup() |>
  group_by(subject_id_comparator) |>
  filter(row_number() == 1) |> 
  ungroup() |>
  group_by(subject_id_comparator, subject_id_reference) |>
  filter(row_number() == 1) |> 
  ungroup() |>
  compute(name = "matched_alternative", temporary = FALSE)

cdm[["matched_alternative"]] <- cdm[["matched_alternative"]] |>
  rename("subject_id" = "subject_id_comparator") |>
  addFutureObservation() |> 
  filter(!is.na(future_observation)) |>
  rename("subject_id_comparator" = "subject_id") |>
  compute(name = "matched_alternative", temporary = FALSE)
  
cdm[["matched_alternative"]] <- cdm[["matched_alternative"]] |>
  select("cohort_definition_id", "subject_id" = "subject_id_reference", "cohort_start_date", "cohort_end_date") |>
  full_join(
    cdm[["matched_alternative"]] |>
      select("cohort_definition_id", "subject_id" = "subject_id_comparator", "cohort_start_date", "cohort_end_date") |>
      mutate("cohort_definition_id" = 2)
  ) |>
  compute(name = "matched_alternative", temporary = FALSE)

# Original matching:
# cdm <- cdm |>
#   generateMatchedCohortSet(
#   name = "matched",
#   targetCohortName = denominator_cohort,
#   targetCohortId   = NULL,
#   matchSex = TRUE,
#   matchYearOfBirth = TRUE,
#   ratio = 1
# )

# cdm$matched <- cdm$matched |>
#   mutate(subject_id = as.numeric(subject_id),
#          cohort_end_date = as.Date(cohort_end_date))

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

# Export tables ----------------------------------------------------------------
write.csv(lsc_analysis1, here("Results",paste0(writePrefix,"lsc_analysis1.csv")))
write.csv(lsc_analysis2, here("Results",paste0(writePrefix,"lsc_analysis2.csv")))
write.csv(lsc_analysis3, here("Results",paste0(writePrefix,"lsc_analysis3.csv")))

# Calculate standarised mean differences ---------------------------------------
# Calculate asmd between cohort and achilles
c1 <- lsc_analysis1 |>
  filter(estimate_name == "percentage") |>
  select(
    "group_level_reference" = "group_level",  
    "variable_name","variable_level",
    "estimate_value_reference" = "estimate_value"
  ) |>
  mutate("pi" = estimate_value_reference) |>
  inner_join(
    lsc_analysis2 |> 
      select(
        "group_level_comparator" = "group_level",
        "variable_name", "variable_level",
        "estimate_value_comparator" = "estimate_value",
      ) |>
      mutate("pk" = estimate_value_comparator),
    by = c("variable_name", "variable_level")
  ) |>
  mutate(pk = if_else(is.na(pk), 0, as.numeric(pk)/100)) |>
  mutate(pi = if_else(is.na(pi), 0, as.numeric(pi)/100)) |>
  mutate(
    "smd" = as.character((pk - pi) / sqrt((pk*(1-pk) + pi*(1-pi))/2)),
    "strata_name_reference" = cohorts$group_level[1],
    "strata_name_comparator" = cohorts$group_level[2]
  ) |>
  select(-c("pi", "pk"))


# Calculate asmd between matched cohorts
x <- lsc_analysis3 |>
  filter(estimate_name == "percentage") 
cohorts <- x |>
  select("group_level", "strata_name", "strata_level") |>
  distinct() |>
  filter(strata_name != "overall")

xi <- x |> 
  filter(group_level == cohorts$group_level[1], strata_name != "overall")
xk <- x |>
  filter(group_level == cohorts$group_level[2], strata_name != "overall")

c2 <- xi |>
  filter(strata_level == 1, estimate_type == "percentage") |>
  select(
    "group_level_reference" = "group_level", 
    "strata_level_reference" = "strata_level", "variable_name", 
    "estimate_value_reference" = "estimate_value", "variable_level", "additional_name", 
    "additional_level"
  ) |>
  mutate("pi" = estimate_value_reference) |>
  inner_join(
    xk |> 
      select(
        "group_level_comparator" = "group_level",
        "strata_level_comparator" = "strata_level", "variable_name", 
        "estimate_value_comparator" = "estimate_value", "variable_level", "additional_name", 
        "additional_level"
      ) |>
      mutate("pk" = estimate_value_comparator),
    by = c("variable_name", "variable_level", "additional_name", "additional_level")
  ) |>
  mutate(pk = if_else(is.na(pk), 0, as.numeric(pk)/100)) |>
  mutate(pi = if_else(is.na(pi), 0, as.numeric(pi)/100)) |>
  mutate(
    "smd" = as.character((pk - pi) / sqrt((pk*(1-pk) + pi*(1-pi))/2)),
    "strata_name_reference" = cohorts$group_level[1],
    "strata_name_comparator" = cohorts$group_level[2]
  ) |>
  select(-c("pi", "pk"))

# Compute tables ---------------------------------------------------------------
# Top 5 lookig at the % in the original cohort
table1 <- lsc_analysis1 |>
  filter(estimate_type == "percentage") |>
  mutate(estimate_value = as.numeric(estimate_value)) |>
  select("variable_name", "variable_level", "estimate_value") |>
  arrange(desc(estimate_value)) |>
  filter(row_number() < 11)

table2 <- c1 |>
  select("variable_name", "variable_level", "estimate_value_reference", 
         "estimate_value_comparator", "smd") |>
  mutate(asmd = abs(as.numeric(smd))) |>
  arrange(desc(asmd)) |>
  filter(row_number() < 11) |>
  select(-"smd")

table3 <- c2 |>
  select("variable_name", "variable_level", "estimate_value_reference", 
         "estimate_value_comparator", "smd") |>
  mutate(asmd = abs(as.numeric(smd))) |>
  arrange(desc(asmd)) |>
  filter(row_number() < 11) |>
  select(-"smd")

# Compare tables ---------------------------------------------------------------
t1 <- table1 |>
  rename("Condition" = "variable_name", "Window" = "variable_level", "Percentage (%)" = "estimate_value") |>
  rename_all(.funs = ~paste0("Original cohort_",.))
  
t2 <- table2 |>
  rename("Condition" = "variable_name", "Window" = "variable_level", 
         "Percentage (%)_Original cohort" = "estimate_value_reference",
         "Percentage (%)_Achilles"        = "estimate_value_comparator",
         "ASMD" = "asmd") |>
  rename_all(.funs = ~paste0("Achilles table comparison_",.))

t3 <- table3 |>
  rename("Condition" = "variable_name", "Window" = "variable_level", 
         "Percentage (%)_Original cohort" = "estimate_value_reference",
         "Percentage (%)_Matched cohort"  = "estimate_value_comparator",
         "ASMD" = "asmd") |>
  rename_all(.funs = ~paste0("Matched cohorts comparison_",.))

# Export comparison tables -----------------------------------------------------
write.csv(t1, here("Results",paste0(writePrefix,"table1.csv")))
write.csv(t2, here("Results",paste0(writePrefix,"table2.csv")))
write.csv(t3, here("Results",paste0(writePrefix,"table3.csv")))


