#------------------------------------------------------------
# Script: simulation_table_M1.R
# Purpose: This script generates Table 1 (Model 1) in the manuscript
# Usage:
#   This script aggregates Monte Carlo simulation results for Model 1,
#   computes median and standard deviation of MSE metrics, and
#   outputs a summary table (both printed and LaTeX via stargazer)
#   It is best practice to open the repository as an Rstudio project (open SQL_JASA.Rproj in top-level directory).
#   Alternatively, set the working directory as the top-level directory and make sure that here package works automatically.
# Running time: < 1 min
#------------------------------------------------------------

# Load required packages (see Reproducibility section in README).
library(purrr)
library(tidyr)
library(dplyr)
library(here)

# Load functions:
source( here( "simulationStudy", "R", "src", "montecarlo_M1.R") )


# Define Simulation Parameters
mcSize <- 100
n_vec <- c( 50, 100, 500, 1000)
fit_args <- expand_grid(mc = 1:mcSize, n = n_vec )

# Load and Aggregate Simulation Results
mse_data <- pmap_dfr(fit_args, possibly( get_mc_results_M1, otherwise = NULL ), folder = here("simulationStudy", "simulations", "M1/") )

# Compute Summary Statistics
mse_table <- mse_data %>% group_by(n) %>% 
  summarise(across( c(sql_z, sql_g, vae_z, vae_g), med_and_sd ) ) %>% ungroup()

# Print the summary table
print( mse_table )

#  Export LaTeX Table via stargazer
stargazer::stargazer(mse_table, summary = FALSE, rownames = FALSE)










