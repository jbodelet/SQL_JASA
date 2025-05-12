fit_M2 <- function(mc, p, q, K = 12, save = FALSE, a = -2, b = 4, lgrid = 50, sim_folder = "./simulationStudy/simulations/" ){
  set.seed(42)
  get_path <- path_M2_function(mc, p, q, sim_folder )
  data <- as.matrix( read.csv( get_path("x") ) )
  gcv <- get_GCV(data, q = q, K = K, Niter = 30, lgrid = lgrid, a = a, b = b )
  lambda <- gcv$lambda_min2
  afm <- AFM( data, q = q, K = K, Niter = 50, method = "Greedy", lambda = lambda )
  # rotation:
  trueZ <- as.matrix( read.csv( get_path("z") ) )
  if(save){
    write_output(afm, get_path )
  }
}



path_M2_function <- function(mc, p, q, sim_folder){
  prefix <- paste0(sim_folder, "M2_q", q, "/data_p", p, "/" )
  postfix <- paste0( "/M2_q", q, "_p", p, "_", mc,".csv" )
  function(...) paste0( prefix, ..., postfix )
}

write_output <- function(afm, get_path){
  # z:
  q <- length( afm$g_eval )
  p <- ncol( afm$g_eval[[1]] )
  zdata <- data.frame( afm$factor_z )
  names(zdata) <- paste0("z", 1:q)
  write.csv( zdata, file = get_path("AFM_fit/z"), row.names = F )
  # g:
  grid <- qnorm(1:200 / 201)
  lapply(1:q, function(l){
    gdata <- as.data.frame( as.matrix( afm$g_eval[[l]] ) )  
    names(gdata) <- paste0("g", 1:p )
    if(q == 1 ){
      path <- get_path( "AFM_fit/g" ) 
    }else{
      path <- get_path( "AFM_fit/g", l)      
    }
    write.csv( cbind(z = grid, gdata), file = path, row.names = F )
  })
}





#=====================
# get results:
#=====================


get_mc_results_M2 <- function(mc, p, q, comp_vae = TRUE, sim_folder = "./simulationStudy/simulations/" ){
  
  folder <- paste0(sim_folder, "M2_q", q, "/data_p", p, "/" )
  file <- paste0( "/M2_q", q, "_p", p, "_", mc,".csv" )
  
  # True parameters:
  true_z <- load_z(folder, file, "")
  true_g <- load_g(folder, file, "", q)
  
  out <- data.frame(q = q, p = p, mc = mc )
  out <- cbind( out, get_mc_mse(folder, file, "AFM_fit/", true_z, true_g ) )
  out <- cbind( out, get_mc_mse(folder, file, "vae_fit/", true_z, true_g ) )
  
  # output:
  return(out)
}


load_z <- function(folder, file, method){
  path <- paste0(folder, method, "z", file )
  as.matrix( read.csv( path ) )
}


load_g <- function(folder, file, method, q = 1 ){
  if( q == 1 ){
    path <- paste0(folder, method, "g", file )
    g <- list( as.matrix( read.csv( path ) )[, -1] )
  }else{
    g <- lapply( 1:q, function(l){
      path <- paste0(folder, method, "g", l, file )
      as.matrix( read.csv( path ) )[, -1]
    })      
  }
  return(g)
}


# Permutes and change sign of fit to best fit to the true
get_mc_mse <- function(folder, file, method, true_z, true_g ){
  est_z <- load_z(folder, file, method)
  q <- ncol(true_z)
  n <- nrow(true_z)
  est_g <- load_g(folder, file, method, q)
  
  # Find the best permutation
  perm <- optimal_permutation(est_z, true_z )
  
  # Find the corresponding signs
  signs <- sign( diag( cor( est_z[, perm], true_z ) ) )
  
  # permute z and find signs:
  if(q== 1){
    est_z <- est_z * signs
  }else{
    est_z <- est_z[, perm, drop = FALSE] %*% diag(signs)
  }

  # permute functions:
  est_g <- est_g[perm]
  
  # inverse functions:
  for(l in 1:q ){
    if(signs[l] < 0 ){
      est_g[[l]] <- est_g[[l]][ n:1, ]
    }
  }
  
  # compute mse:
  mse_z <- mean( (est_z - true_z )^2 )
  mse_g <- mean( sapply( 1:q, function(l) mean( ( est_g[[l]] - true_g[[l]])^2 ) ) )
  
  # output:
  out <- list( mse_z, mse_g)
  if(method == "AFM_fit/") names(out) <- paste0("sql", c("_z", "_g") )
  if(method == "vae_fit/") names(out) <- paste0("vae", c("_z", "_g") )
  return(out)
}

optimal_permutation <- function( est_z, true_z){
  q <- ncol(est_z)
  perms <- permutations(n = q, r = q, v = 1:q )
  out <- apply(perms, 1, function(perm) mean( diag( abs( cor( est_z[, perm], true_z ) ) ) ) )
  return( perms[which.max(out), ] )
}



med_and_sd <- function(x) {
  out <- round( c( median(x), mad(x) ), digits = 3)
  paste( out[1], sprintf("(%.3f)", out[2] ) )
}








