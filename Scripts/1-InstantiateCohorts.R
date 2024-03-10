# instantiate covariates
covariates <- readCohortSet(here("Cohorts"))
cdm <- generateCohortSet(
  cdm = cdm,
  cohortSet = covariates,
  overwrite = TRUE,
  name = "hpv_cin23"
)
