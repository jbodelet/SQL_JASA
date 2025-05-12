library(purrr)
library(furrr)
library(tidyr)
library(dplyr)
library(here)
source( here( "simulationStudy", "R", "src", "AFM.R") )
source( here( "simulationStudy", "R", "src", "montecarlo_M1.R") )
source( here( "simulationStudy", "R", "src", "crossValidation.R") )
options(future.rng.onMisuse = "ignore")



mcSize <- 100
# Set the number of cores for parallel computing (a good rule of thumbe is to use half of the available cores)
ncores <- round( availableCores() / 2 )
plan(multisession, workers = ncores )
fit_args <- expand_grid(mc = 1:mcSize, n = c(50, 100, 500, 1000 ) )

sim_folder <- here("simulationStudy", "simulations/")
future_pwalk(fit_args, fit_M1, sim_folder = sim_folder )




