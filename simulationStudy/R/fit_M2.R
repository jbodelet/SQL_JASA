library(purrr)
library(furrr)
library(tidyr)
library(dplyr)
library(here)
source( here( "simulationStudy", "R", "src", "AFM.R") )
source( here( "simulationStudy", "R", "src", "montecarlo_M2.R") )
source( here( "simulationStudy", "R", "src", "crossValidation.R") )
options(future.rng.onMisuse = "ignore")



mcSize <- 100
p_vec <- c(20, 50, 100, 200, 500)
q_vec <- c(1, 3)


# Set the number of cores for parallel computing (a good rule of thumbe is to use half of the available cores)
ncores <- round( availableCores() / 2 )
plan(multisession, workers = ncores )
fit_args <- expand_grid(mc = 1:mcSize, p = p_vec, q = q_vec )

sim_folder <- here("simulationStudy", "simulations/")
future_pwalk(fit_args, fit_M2, sim_folder = sim_folder, save = TRUE, a = 2, b = 5, K = 8, lgrid = 10 )






