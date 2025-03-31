generate_M1 <- function( sim_path = "../simulations/", n = c(50, 100, 500, 1000), mcSize = 100 ){  
  set.seed(42)  # set seed for reproducibility
  sd_x <- 1
  sd_z <- 1
  theta <- 0
  main_path <- paste0(sim_path, "M1/")
  if( !dir.exists(main_path) ){
    dir.create(main_path)    
  }

  
  gen_Z <- function(n, q, sd_z ) {
    Z <- matrix(rnorm(n*q, 0, sd_z), n, q)
    colnames(Z) <- paste0("z", 1:ncol(Z))
    Z
  }
  g <- list(
    g1 = function(z) z,
    g2 = function(z) z**2/2,
    g3 = function(z) 2*cos(pi*z/1.5),
    g4 = function(z) 2*sin(.5*pi*z),
    g5 = function(z) 4*tanh(1.5*z),
    g6 = function(z) (3*sin(.5*pi*z))/(2 - sin(.5*pi*z)),
    g7 = function(z) 4*sin(.4*pi*z)**3,
    g8 = function(z) {
      delta <- 0.5
      z <- ifelse(abs(z) <= delta, 1/2 * z**2, delta*(abs(z)-1/2*delta))
      4*z - 2
    },
    g9 = function(z) exp(-z**2)*4,
    g10 = function(z) z * cos(3.5*z),
    g11 = function(z) 10*exp(z)/(1+exp(4*z))-2,
    g12 = function(z) rep(0, length(z))
  )
  
  gen_data <- function(n, g, theta, sd_z, sd_x) {
    p <- length(g)
    Z <- gen_Z(n, 1, sd_z )
    X <- sapply(g, function(gj){
      gj(Z) 
    })
    X <- X + gen_Z(n, p, sd_x )
    colnames(X) <- paste0("x", 1:ncol(X))
    list(x=X, z=Z)
  }
  # checking the functions
  check_g <- function(g, which=1:length(g)){
    par(ask=T)
    zgrid <- seq(-3, 3, l=101)
    for (j in 1:length(which)) {
      plot(zgrid, g[[which[j]]](zgrid), main=paste0("g",which[j]))
    }
    par(ask=F)
  }
  # Create path and simulate:
  for (ni in n) {
    subpath <- paste0(main_path, "data_n",ni,"_var", round(sd_x**2,1),"/")
    dir.create(subpath)
    dir.create(paste0(subpath, "x"))
    dir.create(paste0(subpath, "z"))
    dir.create(paste0(subpath, "vae_fit"))
    dir.create(paste0(subpath, "vae_fit/g"))
    dir.create(paste0(subpath, "vae_fit/z"))
    dir.create(paste0(subpath, "AFM_fit"))
    dir.create(paste0(subpath, "AFM_fit/g"))
    dir.create(paste0(subpath, "AFM_fit/z"))
    
    # 2. save g
    z <- qnorm((1:200)/201)
    g_eval <- sapply(g, function(gj){
      gj(z)
    })
    g_eval <- cbind(z, g_eval)
    
    file <- paste0(subpath, "g.csv")
    write.csv(x = g_eval, file=file, row.names=F)
    
    for (sim_num in 1:mcSize) {
      data <- gen_data(ni, g, theta, sd_z, sd_x)
      
      file <- paste0(subpath, "/x/M1_n", ni ,"_var", round(sd_x**2,1) , "_", sim_num, ".csv")
      write.csv(x = data$x, file=file, row.names=F)
      
      file <- paste0(subpath, "/z/M1_n", ni ,"_var", round(sd_x**2,1) , "_", sim_num, ".csv")
      write.csv(x = data$z, file=file, row.names=F)
    }
  }
}

initiate_dir_M2 <- function(p, q, main_path){
  subpath <- paste0(main_path, "data_p", p, "/" )
  dir.create(subpath)
  dir.create(paste0(subpath, "x"))
  dir.create(paste0(subpath, "z"))
  dir.create(paste0(subpath, "vae_fit"))
  dir.create(paste0(subpath, "vae_fit/z"))
  dir.create(paste0(subpath, "AFM_fit"))
  dir.create(paste0(subpath, "AFM_fit/z"))
  
  if( q == 1 ){
    dir.create(paste0(subpath, "g"))
    dir.create(paste0(subpath, "vae_fit/g"))
    dir.create(paste0(subpath, "AFM_fit/g"))
  }else{
    for(l in 1:q){
      dir.create(paste0(subpath, "g", l) )
      dir.create(paste0(subpath, "vae_fit/g", l) )    
      dir.create(paste0(subpath, "AFM_fit/g", l) )    
    }
  }
}



generate_M2_q1 <- function( sim_path = "../simulations/", p_vec = c(20, 50, 100, 200, 500 ), mcSize = 100 ){
  # create M2_q1 folder:
  main_path <- paste0(sim_path, "M2_q1/")
  if( !dir.exists(main_path) ){
    dir.create(main_path)    
  }
  
  simulate_M2_q1 <- function( p, sim_num, sde = 1.5 ){
    n <- 200
    grid <- qnorm(1:n / (n+1) )
    index <- paste0("_p", p )
    subpath <- paste0(main_path, "data", index, "/" )
    sim <- simulate_npfac(n, p, sde )
    x <- sim$data
    colnames(x) <- paste0("x", 1:p )
    z <- data.frame( z1 = sim$factor_z )
    g <- data.frame( z = grid, sim$gfunc_z(grid) )
    colnames(g)[-1] <- paste0("g", 1:p )
    # write csv:
    output_name <- paste0( "M2_q1", index, "_", sim_num )
    write.csv( x, file = paste0( subpath, "x/", output_name, ".csv" ), row.names = FALSE )
    write.csv( z, file = paste0( subpath, "z/", output_name, ".csv" ), row.names = FALSE )
    write.csv( g, file = paste0( subpath, "g/", output_name, ".csv" ), row.names = FALSE )
    invisible()
  }
  for(p in p_vec){
    initiate_dir_M2(p, q = 1, main_path )
    for(mc in 1:mcSize){
      simulate_M2_q1( p, mc )
    }
  }
}



generate_M2_q3 <- function( sim_path = "../simulations/", p_vec = c(20, 50, 100, 200, 500 ), mcSize = 100 ){
  # create M2_q1 folder:
  main_path <- paste0(sim_path, "M2_q3/")
  if( !dir.exists(main_path) ){
    dir.create(main_path)    
  }
  q <- 3
  simulate_M2_q3 <- function(p, sim_num, funcs, main_path ){
    n <- 200
    grid <- qnorm(1:200 / 201 )
    index <- paste0("_p", p )
    subpath <- paste0( main_path, "data", index, "/" )
    sim <- simulate_afm(n, p, q = q, sde = 1.5, funcs = funcs )
    x <- sim$data
    colnames(x) <- paste0("x", 1:p )
    z <- data.frame( sim$factor_z )
    colnames(z) <- paste0("z", 1:q)
    g_list <- lapply( sim$gfunc_z, function(g){
      g <- data.frame( z = grid, g(grid) ) 
      colnames(g)[-1] <- paste0("g", 1:p )
      return(g)
    })
    # write csv:
    output_name <- paste0( "M2_q3", index, "_", sim_num )
    write.csv( x, file = paste0( subpath, "x/", output_name, ".csv" ), row.names = FALSE )
    write.csv( z, file = paste0( subpath, "z/", output_name, ".csv" ), row.names = FALSE )
    lapply(1:q, function(l){
      write.csv( g_list[[l]], file = paste0( subpath, "g", l, "/", output_name, ".csv" ), row.names = FALSE )    
    })   
    invisible()
  }
  funcs <- replicate( q, replicate( max(p_vec) , gen_func(rnorm(8)), simplify = FALSE ), simplify = FALSE )
  
  for(p in p_vec ){
    initiate_dir_M2(p, q = 3, main_path )
    funcs_p <- lapply( funcs, function(f) f[1:p] )
    for(mc in 1:mcSize){
      simulate_M2_q3( p, mc, funcs = funcs_p, main_path = main_path)     
    }
  }
}





