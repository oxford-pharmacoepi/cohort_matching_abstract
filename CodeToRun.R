rm(list = ls())
library("dplyr")
library("DrugUtilisation")
library("duckdb")
library("PatientProfiles")
library("IncidencePrevalence")

cdm <- mockDrugUtilisation(n = 10000)
concept_id_1  <- list("concept" = 378253)
cohort_name_1 <- "random_1"
concept_id_2  <- list("concept" = 317009)
cohort_name_2 <- "random_2"

# Create the general cohort
cdm <- generateConceptCohortSet(cdm,
                                conceptSet = concept_id_1,
                                name = cohort_name_1,
                                overwrite = TRUE)

# Create a cohort with the condition under study
cdm <- generateConceptCohortSet(cdm,
                                conceptSet = concept_id_2,
                                name = cohort_name_2,
                                overwrite = TRUE
)
