# Make sure to set the working directory to simulationStudy folder
source("R/src/AFM.R")
source("R/src/montecarlo_M2.R")
source("R/src/crossValidation.R")
library(purrr)
library(furrr)
library(tidyr)
library(dplyr)
options(future.rng.onMisuse = "ignore")


mcSize <- 100
p_vec <- c(20, 50, 100, 200, 500)
q_vec <- c(1, 3)


ncores <- 30  # Number of cores for parallel computing
plan(multisession, workers = ncores )
fit_args <- expand_grid(mc = 1:mcSize, p = p_vec, q = q_vec )

future_pwalk(fit_args, fit_M2, save = TRUE, a = 2, b = 5, K = 8, lgrid = 10 )






