# Make sure to set the working directory to simulationStudy folder
source("R/src/montecarlo_M1.R")
library(purrr)
library(tidyr)
library(dplyr)


folder = "./simulations/"

mcSize <- 100

n_vec <- c( 50, 100, 500, 1000)

fit_args <- expand_grid(mc = 1:mcSize, n = n_vec )

mse_data <- pmap_dfr(fit_args, possibly( get_mc_results_M1, otherwise = NULL ) )

mse_table <- mse_data %>% group_by(n) %>% 
  summarise(across( c(sql_z, sql_g, vae_z, vae_g), med_and_sd ) ) %>% ungroup()

print( mse_table )

# latex output:
stargazer::stargazer(mse_table, summary = FALSE, rownames = FALSE)










