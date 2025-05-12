# Simulation Study Data Generation Script
# --------------------------------------
# This script generates simulated datasets for Model 1 and Model 2
# as described in Section 4 of the manuscript. All outputs are
# saved in the `simulationStudy/simulations/` directory.
#
# Make sure to pen the project root as an RStudio project to ensure here() works.
library(here)
source( here( "simulationStudy", "R", "src", "data_simulations.R") )
source( here( "simulationStudy", "R", "src", "generate_data.R") )


sim_path <- here("simulationStudy", "simulations/" ) # file to store simulated data

dir.exists(sim_path) # should be true

# The code below generates datasets for simulations, create subfolders in simulationStudy/simulations and save the simulated

# Model 1 
generate_M1(sim_path)  

# MODEL 2 with 1 latent (q=1)
generate_M2_q1(sim_path) 

# MODEL 2 with 3 latents (q=3)
generate_M2_q3(sim_path ) 









