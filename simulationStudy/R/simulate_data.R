#------------------------------------------------------------
# Script: simulate_data.R
# Purpose: This script generates simulated datasets for Model 1 and Model 2 as described in Section 4 of the manuscript.
# Usage:
#   The code below generates datasets for simulations, create subfolders in simulationStudy/simulations and save the simulated.
#   It is best practice to open the repository as an Rstudio project (open SQL_JASA.Rproj in top-level directory).
#   Alternatively, set the working directory as the top-level directory and make sure that here package works automatically.
# Running time: 10 minutes (most of the time is taken my M2_q3)
#------------------------------------------------------------


# Load required packages (see Reproducibility section in README).
library(dplyr)
library(here)

# Load functions:
source( here( "simulationStudy", "R", "src", "data_simulations.R") )
source( here( "simulationStudy", "R", "src", "generate_data.R") )


# file to store simulated data
sim_path <- here("simulationStudy", "simulations/" ) 

# Check:
dir.exists(sim_path) # should be true

# Model 1 
generate_M1(sim_path)  

# MODEL 2 with 1 latent (q=1)
generate_M2_q1(sim_path) 

# MODEL 2 with 3 latents (q=3)
generate_M2_q3(sim_path ) 









