devtools::install_github("jbodelet/SQL/sql")
library(sql)
library(RGAN)
library(torch)
library(transport)
library(MonteCarlo)
library(parallel)
library(future)
library(reshape)
library(ggplot2)
library(Rfast)
library(patchwork)
source("src/architectures.R")
source("src/gan_training.R")
source("src/simulations.R")
source("src/mcSim.R")


mc_run <- function(n, p){
  # Simulations:
  epochs <- 800
  Niter <- 4
  std <- 0.5
  sim <- simulate_model(n, p, std = std, FALSE )
  # GAN:
  #=======
  out_gan <- gan_wrapper(sim, gen_std = std, lr = 5e-4, epochs = 400, Niter = 2 )
  out_gan <- gan_wrapper(sim, out = out_gan, gen_std = std, lr = 1e-4, epochs = epochs, Niter = Niter )
  
  # SQL:
  #=======
  if(p>5){ d <- 5 }else{ d <- 4 }
  out_sql <- get_sql(sim$data, sim$trueGenerator, sim$grid, d = d)
  
  # Output:
  #=========
  out <- c(out_gan[1:3], out_sql[1:3])
  names(out) <- c( paste0("gan_", names(out[1:3])), paste0("sql_", names(out[4:6])) )
  return(out)
}

# Parameters:
param_list <- list(n = c(200, 2000 ), p = c(5, 10, 20, 40, 80) )

set.seed(42)
out <- MonteCarlo(mc_run, nrep = 50, param_list = param_list, ncpus = 35 )
saveRDS(out, "montecarlo_output.RDS")

res <- out$results




#=========================
# Plot results: n = 200 
#=========================


n_index <- 1

long_data_w <- get_combined_data(res$gan_W[n_index,,], res$sql_W[n_index,,], param_list$p)
long_data_k <- get_combined_data(res$gan_ks[n_index,,], res$sql_ks[n_index,,], param_list$p)

plot_w <- ggplot(long_data_w, aes(x = p, y = value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(x = "Dimensionality", y = "Wasserstein" ) +
  theme_minimal() +
  theme(
    # legend.text = element_blank(),      # Adjust legend text size
    legend.position = "none",
    legend.title = element_blank(),
    axis.text = element_text(size = 18),        # Adjust axis text size
    axis.title = element_text(size = 18)        # Adjust axis title size
  )

plot_k <- ggplot(long_data_k, aes(x = p, y = value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(x = "Dimensionality", y = "Kolmogorov" ) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 22),      # Adjust legend text size
    legend.title = element_blank(),
    axis.text = element_text(size = 14),        # Adjust axis text size
    axis.title = element_text(size = 18)        # Adjust axis title size
  )

final_plot <- plot_w|plot_k
final_plot
ggsave("boxplot_sql_gan_n200.png", plot = final_plot, width = 10, height = 4, dpi = 300, units = "in")






#=========================
# Plot results: n = 2000 
#=========================


n_index <- 2

long_data_w <- get_combined_data(res$gan_W[n_index,,], res$sql_W[n_index,,], param_list$p)
long_data_k <- get_combined_data(res$gan_ks[n_index,,], res$sql_ks[n_index,,], param_list$p)

plot_w <- ggplot(long_data_w, aes(x = p, y = value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(x = "Dimensionality", y = "Wasserstein" ) +
  theme_minimal() +
  theme(
    # legend.text = element_blank(),      # Adjust legend text size
    legend.position = "none",
    legend.title = element_blank(),
    axis.text = element_text(size = 18),        # Adjust axis text size
    axis.title = element_text(size = 18)        # Adjust axis title size
  )

plot_k <- ggplot(long_data_k, aes(x = p, y = value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(x = "Dimensionality", y = "Kolmogorov" ) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 22),      # Adjust legend text size
    legend.title = element_blank(),
    axis.text = element_text(size = 14),        # Adjust axis text size
    axis.title = element_text(size = 18)        # Adjust axis title size
  )

final_plot <- plot_w|plot_k
final_plot
ggsave("boxplot_sql_gan_n2000.png", plot = final_plot, width = 10, height = 4, dpi = 300, units = "in")






