# batch = batch size
# lr = learning rate
gan_wrapper <- function( sim, out = NULL, max_std = 2, lr = 1e-4, batch = 200, epochs = 200, Niter = 20, gen_std = NULL ){
  data <- sim$data
  if(is.null(out)){
    distances <- list()
    time <- system.time(
      res <- gan_fit(data, lr = lr, gen_std = gen_std, max_std = max_std, epochs = 10, batch = batch)          
    )/10
  }else{
    distances <- out$distances
    res <- out$res
    time <- out$time
  }
  cat("Estimated time:", time[3] * epochs * Niter )
  iter_index <- 1:Niter + length(distances)
  for(it in iter_index){
    res <- gan_fit(data, res$generator, disc = res$discriminator, lr = lr, gen_std = gen_std, epochs = epochs, batch = batch)
    distances[[it]] <- unlist( get_distances(res$generator, sim$trueGenerator) )
    # print( paste("iter:", it) )
    print( round( distances[[it]], 3 ) )
  }
  #===========
  # Output:
  #===========
  distances_mat <- do.call(rbind, distances)
  out <- as.list( apply( distances_mat, 2, min ) )
  out$res <- res
  out$Ghat <- get_expected_generator(res$generator, sim$grid, ncol(data))
  out$distances <- distances
  out$time <- time # time per epoch
  return(out)
}



# print( cat( "wasserstein: ", distances[it, 1], "Kolmogorov:", distances[it, 2], "mse:", distances[it, 3] ) )
get_distances <- function(generator, trueG){
  grid <- get_grid(nrow(trueG))
  Ghat <- get_expected_generator(generator, grid, ncol(trueG))
  get_output(Ghat, trueG )
}






get_grid <- function(m, scale = "normal"){
  as.matrix( qnorm( 1:m / (m+1) ) )
}



gan_fit <- function( data, gen = NULL, disc = NULL, epochs = 2000, gen_std = NULL, max_std = 2, plot_dim = NULL, lr = 1e-4, batch = 200){
  n <- nrow(data)
  p <- ncol(data)
  noise_dim <- p + 1
  if( n < batch){
    batch <- n
  }
  if(is.null(plot_dim)){
    plot_dim <- sample(2:p, 1)
    # plot_dim <- min(p,7)
  }
  #==================
  # Architecture:
  #==================
  if(is.null(gen)){
    # if null, initiate the generator
    if(is.null(gen_std)){# with unknown standard deviation
      gen <- get_generator(p, max_std)    
    }else{# with known standard deviation
      gen <- get_generator_std(p, gen_std)
    }
  }
  if(is.null(disc)){
    disc <- get_discriminator(p)    
  }
  WGAN_weight_clipper(disc, clip_values = c(-0.01, 0.01))
  
  #==================
  # Hyperparameters:
  #==================
  optimizer_gen <- optim_adam(gen$parameters, lr = lr )
  optimizer_disc <- optim_adam(disc$parameters, lr = lr )
  
  #==================
  # Training:
  #==================
  res <- gan_trainer(
    data,
    generator = gen,
    discriminator = disc,
    generator_optimizer = optimizer_gen,
    discriminator_optimizer = optimizer_disc,
    noise_dim = noise_dim,
    epochs = epochs,
    batch = batch,
    value_function = "wasserstein",
    plot_progress = TRUE,
    plot_interval = 400,
    plot_dimensions = c(1,plot_dim)
  )
  return(res)
}




get_sql <- function(data, trueG, grid, d= 4){
  res <- sql::SQL(data, q = 1, d = d, greedy = TRUE, use_Rfast = TRUE)
  Ghat <- sql:::predict.sql(res, grid)[[1]]
  out <- get_output(Ghat, trueG )
  out$res <- res
  return(out)
}

get_expected_generator <- function(gen, grid, p){
  m <- nrow(grid)
  input <- cbind( grid, matrix(0, ncol = p, nrow = m) )
  return( as.matrix( gen(torch_tensor(input)) ) )
}



get_output <- function(Ghat, trueG){
  wasserstein_dist <- my_wasserstein(Ghat, trueG)
  ks <- get_kolmogorov_distance(Ghat, trueG)
  mse <- get_mse_g(trueG, Ghat)
  return( list( W = wasserstein_dist, ks = mean(ks), mse = mse ) )
}


get_kolmogorov_distance <- function(x,y){
  p <- ncol(x)
  return( sapply( 1:p, function(j) ks.test(x[,j], y[, j])$statistic ) )
}






get_mse_g <- function(trueG, Ghat){
  m <- nrow(Ghat)
  mse_g1 <- mean( ( trueG - Ghat )^2 )
  mse_g2 <- mean( ( trueG - Ghat[m:1, ] )^2 )
  return( min(mse_g1, mse_g2) )
}


plot_generator <- function( gen, plot_dim = NULL, m = 1000){
  if(is.null(plot_dim)){
    plot_dim <- sample(1:p, 1)
  }
  grid <- as.matrix(qnorm(1:m / (m+1)) )
  # Compute true:
  Gtrue <- as.matrix( gfunc_toy(grid) )
  p <- ncol(Gtrue)
  # Compute estimated generator:
  Ghat <- get_expected_generator(grid, p)  
  # Plot:
  plot( Gtrue[,plot_dim] ~ grid, type = "l", lwd = 4, col = "gray")
  lines( Ghat[,plot_dim] ~ grid, type = "l", lwd = 2, col = "blue")  
}


my_wasserstein <- function(x, y){
  n <- nrow(x)
  a <- rep(1/n, n)
  # cost_matrix <- sapply( 1:n, function(j) sapply(1:n, function(i) sqrt(mean( (x[i,]- y[j,] )^2 ) ) ) )
  cost_matrix <- sqrt( ( t( t( - 2 * x %*% t(y) + diag(x %*% t(x)) ) + diag(y %*% t(y)) ) ) / ncol(x) )
  out <- wasserstein(a, a, costm = cost_matrix)
  return(out)
}
