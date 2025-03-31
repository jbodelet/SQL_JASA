# Make sure to set the working directory to simulationStudy folder
source("R/src/data_simulations.R")
source("R/src/generate_data.R")


sim_path <- "./simulations/" # file to store simulated data

dir.exists(sim_path) # should be true

# Generate datasets for simulations
# The data is stored in the simulations folder.

generate_M1(sim_path)  # Model 1 

generate_M2_q1(sim_path) # MODEL 2 with 1 latent (q=1)

generate_M2_q3(sim_path ) # MODEL 2 with 3 latents (q=3)









