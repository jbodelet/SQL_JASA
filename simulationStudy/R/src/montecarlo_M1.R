fit_M1 <- function(mc, n, sim_folder = "./simulationStudy/simulations/", lgrid = 50 ){
  set.seed(42)
  grid <- qnorm(1:200/ 201)
  label_n <- paste0("_n", n, "_var1" )
  filePath  <- paste0( sim_folder, "M1/")
  filePath_2data <- paste0( filePath, "data",label_n )
  file_name <- paste0( "M1", label_n, "_", mc,".csv" )
  data <- as.matrix( read.csv( paste0( filePath_2data, "/x/", file_name ) ) )
  
  # FIT:
  K <- 12
  b <- 4
  a <- exp(-5)
  gcv <- get_GCV( data, q = 1, K = K, Niter = 30, a = a, b = b, lgrid = lgrid, smooth_df = 8 )
  lambda <- gcv$lambda_min2
  abline(v = lambda)
  afm <- AFM( data, q = 1, K = K, Niter = 100, method = "Greedy", lambda =  lambda )
  # save:
  zdata <- data.frame( "z1" = afm$factor_z)
  write.csv( zdata, file = paste0( filePath_2data, "/AFM_fit/z/", file_name ), row.names = F )
  gdata <- as.data.frame( afm$g_eval[[1]] )
  names(gdata) <- paste0("g", 1:12)
  write.csv( cbind(z = grid, gdata), file = paste0( filePath_2data, "/AFM_fit/g/", file_name ), row.names = F )
  
}




get_mc_results_M1 <- function(mc, n, folder = "./simulationStudy/simulations/M1/" ){
  get_mse <- function(mc, n, method ){
    # get data:
    label_n <- paste0("_n", n, "_var1" )
    filePath_2data <- paste0( folder, "data", label_n, "/" )
    g0 <- as.matrix( read.csv( paste0( filePath_2data, "g.csv") ) )
    z <- as.matrix( read.csv( paste0( filePath_2data, "z/M1", label_n, "_", mc, ".csv") ) )
    zhat <- as.matrix( read.csv( paste0(filePath_2data, method, "/z/", "M1",label_n, "_", mc,".csv") ) )
    ghat <- as.matrix( read.csv( paste0(filePath_2data, method, "/g/", "M1",label_n, "_", mc,".csv") ) )
    # scale:
    g0 <- scale(g0[, -1], scale = F)
    ghat <- scale(ghat[, -1], scale = F)
    # mse up to a sign:
    mse_z <- min( mean( ( z - zhat )^2 ), mean( ( z + zhat )^2 ) )
    mse_g <- min( mean( ( g0 - ghat )^2 ), mean( ( g0 - ghat[200:1, ] )^2 ) )
    out <- list( mse_z, mse_g)
    if(method == "AFM_fit") names(out) <- paste0("sql", c("_z", "_g") )
    if(method == "vae_fit") names(out) <- paste0("vae", c("_z", "_g") )
    return(out)
  }
  out <- data.frame(n = n, mc = mc, get_mse(mc, n, "AFM_fit") )
  out <- cbind( out, get_mse(mc, n, "vae_fit") )
  # output:
  return(out)
}


med_and_sd <- function(x) {
  out <- round( c( median(x), mad(x) ), digits = 3)
  paste( out[1], sprintf("(%.3f)", out[2] ) )
}

