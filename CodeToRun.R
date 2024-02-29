rm(list = ls())
library("dplyr")
library("DrugUtilisation")
library("duckdb")
library("PatientProfiles")
library("IncidencePrevalence")

cdm <- mockDrugUtilisation(n = 10000)
denominator_concept_id <- list("concept" = 378253)
denominator_cohort  <- "random_1"
outcome_concept_id  <- list("concept" = 317009)
outcome_cohort <- "random_2"

# Create the general cohort
cdm <- generateConceptCohortSet(cdm,
                                conceptSet = denominator_concept_id,
                                name       = denominator_cohort,
                                overwrite  = TRUE)

# Create a cohort with the condition under study
cdm <- generateConceptCohortSet(cdm,
                                conceptSet = outcome_concept_id,
                                name       = outcome_cohort,
                                overwrite  = TRUE
)

# Calculate overall incidence --------------------------------------------------
incidence_table <- estimateIncidence(
  cdm,
  denominatorTable = denominator_cohort,
  outcomeTable     = outcome_cohort,
  denominatorCohortId = NULL,
  outcomeCohortId     = NULL,
  interval = "overall",
  completeDatabaseIntervals = TRUE,
  outcomeWashout = Inf,
  repeatedEvents = FALSE,
  minCellCount   = 5,
  strata = list(),
  includeOverallStrata = TRUE,
  temporary = TRUE, 
  returnParticipants = FALSE
)


# Calculate incidence using achilles -------------------------------------------




