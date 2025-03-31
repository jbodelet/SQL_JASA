mc_run <- function(n, p){
  # GAN doesn't know the true standard deviation
  # Simulations:
  sim <- simulate_model(n, p, std = 1.5, FALSE )
  # GAN:
  out_gan <- gan_superfit(sim$data, sim$trueGenerator, sim$grid, gen_std = NULL, epochs = 400, Nepochs = 10)
  # mse_gam <- get_mse_g(sim$trueGenerator, out_gan$Ghat)
  # SQL:
  out_sql <- get_sql(sim$data, sim$trueGenerator, sim$grid, d = 6)
  # mse_sql <- get_mse_g(sim$trueGenerator, out_sql$Ghat)
  # Output:
  return(list(W_sql = out_sql$wasserstein_dist, W_gan = out_gan$wasserstein_dist,
              ks_sql = mean( out_sql$ks ), ks_gan = mean( out_gan$ks )
              # , mse_sql = mse_sql, mse_gam = mse_gam
              ) 
         )
}


  

mc_run_std <- function(n, p){
  # GAN knows the true standard deviation
  # Simulations:
  sim <- simulate_model(n, p, std = 0.1, FALSE )
  grid <- sim$grid
  gen_std = sim$std
  # GAN:
  out_gan <- gan_superfit(sim$data, sim$trueGenerator, sim$grid, gen_std = gen_std, epochs = 400, Nepochs = 20)
  out_sql <- get_sql(sim$data, sim$trueGenerator, sim$grid)
  return(list(W_sql = out_sql$wasserstein_dist, W_gan = out_gan$wasserstein_dist,
              ks_sql = mean( out_sql$ks ), ks_gan = mean( out_gan$ks ) ) )
}


mc_run_test <- function(n, p){
  # GAN doesn't know the true standard deviation
  # Simulations:
  sim <- simulate_model(n, p, std = 1.5, FALSE, FALSE )
  grid <- sim$grid
  gen_std = NULL
  # GAN:
  out_gan <- gan_superfit(sim$data, sim$trueGenerator, sim$grid, gen_std = gen_std, epochs = 4, Nepochs = 1)
  out_sql <- get_sql(sim$data, sim$trueGenerator, sim$grid)
  return(list(W_sql = out_sql$wasserstein_dist, W_gan = out_gan$wasserstein_dist,
              ks_sql = mean( out_sql$ks ), ks_gan = mean( out_gan$ks ) ) )
}








