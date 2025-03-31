# Make sure to set the working directory to simulationStudy folder
source("R/src/AFM.R")
source("R/src/montecarlo_M1.R")
source("R/src/crossValidation.R")
library(purrr)
library(furrr)
library(tidyr)
library(dplyr)
options(future.rng.onMisuse = "ignore")


mcSize <- 100
ncores <- 30 # Number of cores for parallel computing
plan(multisession, workers = ncores )
fit_args <- expand_grid(mc = 1:mcSize, n = c(50, 100, 500, 1000 ) )

future_pwalk(fit_args, fit_M1 )




