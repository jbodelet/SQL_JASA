#------------------------------------------------------------
# Script: simulation_table_M2.R
# Purpose: This script generates Table 2 (Model 2) in the manuscript
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
library(gtools)
library(here)

# Load functions:
source( here( "simulationStudy", "R", "src", "montecarlo_M2.R") )


# Define Simulation Parameters
mcSize <- 100
p_vec <- c(20, 50, 100, 200, 500)
q_vec <- c(1, 3)
fit_args <- expand_grid(mc = 1:mcSize, p = p_vec, q = q_vec )

# Load and Aggregate Simulation Results
mse_data <- pmap_dfr(fit_args, possibly( get_mc_results_M2, otherwise = NULL ), comp_vae = TRUE, sim_folder = here( "simulationSudy", "simulations/" ) )


# Compute Summary Statistics
mse_table <- mse_data %>% group_by(q,p) %>% 
  summarise(across( c(sql_z, sql_g, vae_z, vae_g), med_and_sd ) )


# Print the summary table
print( mse_table )


#  Export LaTeX Table via stargazer
stargazer::stargazer(mse_table, summary = FALSE, rownames = FALSE)



