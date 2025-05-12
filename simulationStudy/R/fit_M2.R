# Simulation Study: Model 2 Fitting Script
# ----------------------------------------
# This script runs Monte Carlo simulations to fit SQL for Model 1 across multiple sample sizes and repetitions.
# Results are saved in the `simulationStudy/simulations/` directory for later analysis.

library(purrr)
library(furrr)
library(tidyr)
library(dplyr)
library(here)
source( here( "simulationStudy", "R", "src", "AFM.R") )
source( here( "simulationStudy", "R", "src", "montecarlo_M2.R") )
source( here( "simulationStudy", "R", "src", "crossValidation.R") )


# Configure Parallel Processing
options(future.rng.onMisuse = "ignore")
mcSize <- 100
p_vec <- c(20, 50, 100, 200, 500)
q_vec <- c(1, 3)


# Set the number of cores for parallel computing (a good rule of thumb is to use half of the available cores)
ncores <- round( availableCores() / 2 )
plan(multisession, workers = ncores )
fit_args <- expand_grid(mc = 1:mcSize, p = p_vec, q = q_vec )

# Path to directory containing simulated datasets
sim_folder <- here("simulationStudy", "simulations/")

# Montecarlo runs (results will be saved in sim_folder directly)
future_pwalk(fit_args, fit_M2, sim_folder = sim_folder, save = TRUE, a = 2, b = 5, K = 8, lgrid = 10 )






