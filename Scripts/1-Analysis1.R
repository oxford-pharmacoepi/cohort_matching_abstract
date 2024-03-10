info(logger, "START LARGE SCALE CHARACTERISATION - ORIGINAL COHORT")
tic(msg = "- Large scale characterisation of the cases")

# Large scale characterization of the cases (conditions, drug_exposures) 
info(logger, "- Large scale characteristics of the cases")

lsc_analysis1 <- cdm[[denominator_cohort]] |>
  summariseLargeScaleCharacteristics(
    strata = list(),
    window = list( c(-365, -31), c(-30,-1)),
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