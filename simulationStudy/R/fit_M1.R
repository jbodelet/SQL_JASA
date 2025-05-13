#------------------------------------------------------------
# Script: fit_M1.R
# Purpose: This script runs Monte Carlo simulations to fit SQL for Model 1 across multiple sample sizes and repetitions.
# Usage:
#   Results are saved in the `simulationStudy/simulations/` directory for later analysis.
#   It is best practice to open the repository as an Rstudio project (open SQL_JASA.Rproj in top-level directory).
#   Alternatively, set the working directory as the top-level directory and make sure that here package works automatically.
# Running time without parallelization: 2 hours
#------------------------------------------------------------

# Load required packages (see Reproducibility section in README).
library(purrr)
library(purrr)
library(furrr)
library(tidyr)
library(dplyr)
library(here)

# Load functions:
source( here( "simulationStudy", "R", "src", "AFM.R") )
source( here( "simulationStudy", "R", "src", "montecarlo_M1.R") )
source( here( "simulationStudy", "R", "src", "crossValidation.R") )


# Configure Parallel Processing
options(future.rng.onMisuse = "ignore")
ncores <- round( availableCores() / 2 )  # Set the number of cores for parallel computing (a good rule of thumb is to use half of the available cores)
plan(multisession, workers = ncores )  

# Set montecarlo parameters:
mcSize <- 100 # Montecarlo size
fit_args <- expand_grid(mc = 1:mcSize, n = c(50, 100, 500, 1000 ) )

# Path to directory containing simulated datasets
sim_folder <- here("simulationStudy", "simulations/")

# Montecarlo runs (results will be saved in sim_folder directly)
future_pwalk(fit_args, fit_M1, sim_folder = sim_folder )




