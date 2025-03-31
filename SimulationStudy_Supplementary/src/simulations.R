simulate_model <- function(n, p, std = 0.05, model1 = FALSE, by_max = FALSE, latent = NULL){
  # Get generator:
  #================
  if(model1){
    gfunc <- gfunc_model1 
  }else{
    gfunc <- gfunc_model2(p)
  }
  gfunc <- rescale_generator(gfunc, by_max)
  # true generator (compute for checking results):
  m <- 1000
  grid <- qnorm(1:m/(m+1))
  trueGenerator <- gfunc(grid)[,1:p]  # keep only the first p variables
  
  # Simulate data:
  #================
  if(is.null(latent)){
    latent <- rnorm(n) 
  }
  X <- gfunc(latent)[,1:p]
  eps <- matrix( rnorm(p*n, sd = std), ncol = p )
  data <- X + eps
  sim <- list( data = data, common = X, latent = latent, std = std, trueGenerator = trueGenerator, 
               grid = as.matrix(grid), model1 = model1)
  return(sim)
}


# get mean:
rescale_generator <- function(gfunc, by_max = TRUE){
  m <- 1000
  rand_grid <- rnorm(m)
  Gmean <- colMeans(gfunc(rand_grid))
  Gmean_func <- function(m) t(Gmean) %x% rep(1,m)
  if(by_max){
    # get max
    Gmax <- colMeans( t(replicate( 10, apply( abs( gfunc(rand_grid)- Gmean_func(m) ), 2, max) )) ) 
    Gmax[Gmax == 0] <- 1
    Gmax_func <- function(m) t(Gmax) %x% rep(1,m)
    # scaled generator:
    gfunc_scaled <- function(z) 0.9 * ( gfunc(z) - Gmean_func(length(z)) ) / Gmax_func(length(z))    
  }else{
    # get sd
    Gsd <- sqrt( colMeans( ( gfunc(rand_grid)- Gmean_func(m) )^2 ) )
    Gsd[Gsd == 0] <- 1
    Gsd_func <- function(m) t(Gsd) %x% rep(1,m)
    # scaled generator:
    gfunc_scaled <- function(z) ( gfunc(z) - Gmean_func(length(z)) ) / Gsd_func(length(z))
  }
  return(gfunc_scaled)
}


gfunc_model1 <- function(zz ){
  J <- function(z, b) z^2 * (abs(z)<b) + (abs(z) - b^2) * (abs(z) >= b)
  t( sapply( zz, function(z){
    out <- rep(0, 12 )
    out[1] <- z
    out[2] <- z^2
    out[3] <- cos(pi * z /1.5 )
    out[4] <- sin(pi * z / 2)
    out[5] <- tanh(1.5* z)
    out[6] <- sin(pi * z / 2) / (2- sin(pi *z/2) )
    out[7] <- sin(pi*z*0.4)^3
    out[8] <- J( z, 0.5)
    out[9] <- exp(-z^2)
    out[10] <- z * cos(3.5 * z)
    out[11] <- exp(z) / (1 + exp(4* z))
    out[12] <- z * 0
    return( out )
  }) )
}

gfunc_model2 <- function( p, M = 4 ){
  a <- matrix( rnorm(p * M ), ncol = p ) / (1:M)
  b <- matrix( rnorm(p * M ), ncol = p ) / (1:M)
  f <- function(zz) t(sapply( zz, function(z){ 2 * colSums( a * cos( pi * 1:M * z / 4 ) + b * sin( pi * 1:M * z / 4 ) ) } ))
  # f <- function(zz) t(sapply( zz, function(z){ 
  #   constant <- sqrt( colSums( a^2 + b^2 ) )
  #   colSums( a * cos( pi * 1:M * z / 4 ) + b * sin( pi * 1:M * z / 4 ) ) /  constant
  #   } ) )
  return(f)
}


a <- matrix(rnorm(100), nrow = 5)
b <- matrix(rnorm(100), nrow = 5)

get_boxplot <- function(a,b, p_index, method = c("GAN", "SQL"), statistics = "MSE"){
  # a and b should be matrices, where the columns are the repetitions
  a_df <- a %>% as.data.frame()
  b_df <- b %>% as.data.frame()
  
  # Add a column to indicate the sample size
  a_df$p <- factor(p_index)
  b_df$p <- factor(p_index)
  
  # Melt the data frames from wide to long format
  a_long <- melt(a_df, id.vars = 'p')
  b_long <- melt(b_df, id.vars = 'p')
  
  # Add a column to indicate the method
  a_long$Method <- method[1]
  b_long$Method <- method[2]
  
  # Combine both datasets
  combined_data <- rbind(a_long, b_long)
  
  # Create the boxplot with ggplot2
  ggplot(combined_data, aes(x = p, y = value, fill = Method)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    labs(x = "Dimensionality", y = statistics ) +
    theme_minimal() +
    theme(
      legend.text = element_text(size = 22),      # Adjust legend text size
      legend.title = element_blank(),
      axis.text = element_text(size = 18),        # Adjust axis text size
      axis.title = element_text(size = 24)        # Adjust axis title size
    )
}



get_combined_data <- function(a,b, p_index, method = c("GAN", "SQL") ){
  # a and b should be matrices, where the columns are the repetitions
  a_df <- a %>% as.data.frame()
  b_df <- b %>% as.data.frame()
  
  # Add a column to indicate the sample size
  a_df$p <- factor(p_index)
  b_df$p <- factor(p_index)
  
  # Melt the data frames from wide to long format
  a_long <- melt(a_df, id.vars = 'p')
  b_long <- melt(b_df, id.vars = 'p')
  
  # Add a column to indicate the method
  a_long$Method <- method[1]
  b_long$Method <- method[2]
  
  # Combine both datasets
  return( rbind(a_long, b_long) )
}









