library(purrr)
library(tidyr)
library(dplyr)
library(gtools)
library(here)
source( here( "simulationStudy", "R", "src", "montecarlo_M2.R") )

sim_folder = here( "simulationSudy", "simulations/" )

mcSize <- 100

p_vec <- c(20, 50, 100, 200, 500)

q_vec <- c(1, 3)

fit_args <- expand_grid(mc = 1:mcSize, p = p_vec, q = q_vec )

mse_data <- pmap_dfr(fit_args, possibly( get_mc_results_M2, otherwise = NULL ), comp_vae = TRUE, sim_folder = folder )

mse_table <- mse_data %>% group_by(q,p) %>% 
  summarise(across( c(sql_z, sql_g, vae_z, vae_g), med_and_sd ) )

print( mse_table )

# latex output:
stargazer::stargazer(mse_table, summary = FALSE, rownames = FALSE)



