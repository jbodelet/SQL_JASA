library(cvTools)
library(purrr)
library(furrr)

#' @export
get_GCV <- function(x, q, K, Niter, lambdaGrid = NULL, a = -4, b = 3, lgrid = 50, ncores = 1, smooth_df = 6 ){
  if(is.null(lambdaGrid)){
    lambdaGrid <- exp(seq(a, b, l = lgrid ) )    
  }
  if(ncores == 1 ){
    out <- lapply( lambdaGrid, function(lam){ 
      AFM(x, lambda = lam, q = q, K = K, Niter = Niter, method = "Greedy" ) 
    }
    )    
  }else{
    plan(multisession, workers = ncores)
    out <- furrr::future_map(lambdaGrid, function(lam){ 
      AFM(x, lambda = lam, q = q, K = K, Niter = Niter, method = "Greedy" ) 
    })
  }
  mse <- sapply(out, function(x) tail(x$mse, 1) )
  df <- sapply(out, function(x) x$df)
  gcv <- sapply(out, function(x) x$gcv)
  
  # smooth version:
  smoothed <- smooth.spline(mse ~ lambdaGrid, df = smooth_df)
  mse2 <- smoothed$y
  gcv2 <- mse2 / (1- df / nrow(x) )^2
  
  # Plot:
  plot(gcv ~ lambdaGrid)
  lines( gcv2  ~ lambdaGrid)
  lambda_min <- lambdaGrid[ which.min(gcv) ]
  lambda_min2 <- lambdaGrid[ which.min(gcv2) ]
  return(list(errors = gcv, lambdaGrid = lambdaGrid, lambda_min = lambda_min, lambda_min2 = lambda_min2 ) )
}






#' @export
cross_validate <- function (data, lambdaGrid = 1:20, ncores = 1, nfold = 5, q = 1, K = 12, Niter = 20, random_p = NULL){
  folds <- cvFolds(NROW(data), K = nfold)
  compute_cv_error <- function(lambda) {
    errors <- sapply(1:nfold, function(fold) {
      trainSet <- data[folds$which != fold, ]
      testSet <- data[folds$which == fold, ]
      fit <- AFM(trainSet, q, lambda = lambda, K, Niter, fastProd = TRUE, method = "Greedy", random_p = random_p)
      factors <- estimate_factors(data, fit$gfunc_z, method = "Greedy", Niter = 2)$factor_z
      factors <- as.matrix(factors[folds$which == fold, ])
      pred <- as.matrix(Reduce("+", lapply(1:q, function(l) fit$gfunc_z[[l]](factors[, l]))))
      return(mean((testSet - pred)^2))
    })
    out <- data.frame( cv_m = mean(errors), cv_sd = sd(errors) / sqrt( nfold ) )
    return( out )
  }
  if (ncores > 1) {
    plan(multisession, workers = ncores)
    mse <- furrr::future_map_dfr(lambdaGrid, compute_cv_error)
  }
  else {
    mse <- map_dfr(lambdaGrid, compute_cv_error)
  }
  return(list(mse = mse, lambdaGrid = lambdaGrid))
}


#' @export
estimate_factors <- function( x, gfunc_z, Niter = 10, method = "Hungarian", fastProd = FALSE, random_p = NULL ){
  # Initialization:
  q <- length(gfunc_z)
  n <- nrow(x)
  grid <- 1:n / (n+1)
  factor <- get_initial_value(x, q )
  ord <- apply( factor, 2, order )
  P <- apply( ord, 2, get_permuationMatrix )
  G <- lapply(gfunc_z, function(g) g(qnorm(grid)) )
  mse <- list()
  # Algorithm:
  for( it in 1:Niter ){
    P <- est_P(x, P, G, method, fastProd, random_p )
  }
  # Output:
  factor <- sapply( P, function(p) as.numeric( p %*% grid ) )
  return( list( factor_u = factor, factor_z = apply( factor, 2, qnorm ), P = P ) )
}