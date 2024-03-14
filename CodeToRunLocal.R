# Load required packages -------------------------------------------------------
library(DBI)
library(CDMConnector)
library(log4r)
library(here)
library(readr)
library(zip)
library(RPostgres)

# Connect to database ----------------------------------------------------------
# please see examples how to connect to the database here:
# https://darwin-eu.github.io/CDMConnector/articles/a04_DBI_connection_examples.html
server_dbi <- Sys.getenv("server_dbi")
user       <- Sys.getenv("user")
port       <- Sys.getenv("port")
host       <- Sys.getenv("host")

#Connect to the database
db <- dbConnect(RPostgres::Postgres(),
                dbname = server_dbi,
                port = port,
                host = host,
                user = user,
                password = Sys.getenv("password"))


# name of the schema where cdm tables are located
cdmSchema <- "public_100k"

# name of a schema in the database where you have writing permission
writeSchema <- "results"

# combination of at least 5 letters + _ (eg. "abcde_") that will lead any table
# written in the write schema
writePrefix <- "mah_cohort_matching_100k"

# name of the database, use acronym in capital letters (eg. "CPRD GOLD")
dbName <- "CPRD GOLD"

# minimum number of counts to be reported
minCellCount <- 5

achilles_schema <- "results"

# Create cdm object ------------------------------------------------------------
cdm <- cdm_from_con(db, 
                    cdm_schema = cdmSchema, 
                    write_schema = c(schema = writeSchema, 
                                     prefix = writePrefix),
                    cdm_name = dbName,
                    achilles_schema = achilles_schema, cohort_tables = "hpv_cin23"
)

# Create log file
resultsFolder <- here("Results")
log_file <- paste0(resultsFolder, paste0("/",writePrefix,"log.txt"))
logger <- create.logger(logfile = log_file, level = "INFO")
info(logger = logger, "START RUN STUDY")


# Instantiate cohort -----------------------------------------------------------
covariates <- readCohortSet(here("Cohorts"))
cdm <- generateCohortSet(
  cdm = cdm,
  cohortSet = covariates,
  overwrite = TRUE,
  name = "hpv_cin23"
)

