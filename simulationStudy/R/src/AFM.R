library(Matrix)
library(irlba)

#' @export
AFM <- function( x, q = 1, K = 4, Niter = 10, lambda = 0.1, method = "Hungarian", fastProd = FALSE, random_p = NULL, Kvec = NULL ){
  # Initialization:
  if(is.null(Kvec)){
    Kvec <- round( seq( 4, K, length.out = Niter ) ) # vector of the number of splines 
  }
  n <- nrow(x)
  grid <- 1:n / (n+1)
  factor <- get_initial_value(x, q )
  ord <- apply( factor, 2, order )
  P <- apply( ord, 2, get_permuationMatrix )
  mse <- list()
  # Algorithm: it <- 1
  for( it in 1:Niter ){
    splines <- get_splines_basis( K = Kvec[it] + 1 )
    psi <- get_psi(splines$basis, grid, center = TRUE )
    d2psi <- splines$d2psi[1:Kvec[it], 1:Kvec[it]]
    beta <- est_beta( x, psi, P, K = Kvec[it], d2psi = d2psi, lambda = lambda )
    G <- lapply( beta, function(b) psi %*% b )
    P <- est_P(x, P, G, method, fastProd, random_p )
    mse_l <- sapply(1:q, function(l) get_mse(x, P[l], G[l] ) )
    ev <- 1- mse_l / Matrix::mean( ( x - mean(x) )^2 )
    factor_ordering <- q + 1 - rank( ev )
    P <- P[factor_ordering]
    G <- G[factor_ordering]
    beta <- beta[factor_ordering]
    # get mse:
    mse[[it]] <- get_mse(x, P, G)
  }
  
  # Output:
  factor <- sapply( P, function(p) as.numeric( p %*% grid ) )
  df <- get_degree_of_freedom(psi, P, K = Kvec[it], d2psi = d2psi, lambda = lambda )
  gcv <- mse[[it]] / (1 - df/n)^2
  gfunc_z <- lapply( beta, function(b){
    function(z){
      get_psi(splines$basis, pnorm(z) ) %*% b 
    }
  })
  g_eval <- lapply( gfunc_z, function(g) as.matrix( g(qnorm(1:200/201 ) ) ) )
  return( list( factor_u = factor, factor_z = apply( factor, 2, qnorm ), P = P, G = G, beta = beta,
                gfunc_z = gfunc_z, splines = splines, mse = unlist( mse ), K = K, method = method, q = q,
                g_eval = g_eval, lambda = lambda, Niter = Niter, df = df, gcv = gcv ) )
}





#' @export
get_mse <- function(x, P, G ){
  q <- length(P)
  Matrix::mean( ( x - Reduce( '+', lapply( 1:q, function(l) P[[l]] %*% G[[l]] ) ) )^2 )
}



#' @export
get_explainedVar <- function(x, P, G){
  q <- length(G)
  mse_l <- sapply(1:q, function(l) get_mse(x, P[1:l], G[1:l] ) )
  totalVar <- Matrix::mean( ( x - mean(x) )^2 )
  ev <- 1- mse_l / totalVar
  return(ev)
}

#' @export
set_as_ref <- function( est, ref ){
  q <- ncol( ref )
  n <- nrow(ref)
  # 1) permutations:
  if(q > 1 ){
    per <- gtools::permutations(q, q)
    per_star <- per[ which.max( apply( per, 1, function(u){
      sum( abs( diag( cor( ref, est$factor_z[, u ] ) ) ) )
    })), ]
  }else{
    per_star <- 1  
  }
  factor_z <- est$factor_z[, per_star]
  gfunc_z <- est$gfunc_z[per_star]
  g_eval <- est$g_eval[per_star]    
  # 2) switch signs:
  signs <- sign( diag( t( factor_z ) %*% ref ) )
  if(q>1){
    factor_z <- factor_z %*% diag(signs)     
  }else{
    factor_z <- factor_z * signs
  }
  for(l in 1:q){
    if( signs[l] < 0 ){
      g_eval[[l]] <- g_eval[[l]][200:1, ]
    }
  }
  return( list(factor_z = factor_z, g_eval = g_eval,
               mse = est$mse, K = est$K ) )
}







#==============================
# HELPER FUNCTIONS:
#==============================


get_splines_basis <- function( K, range = c(0,1), degree = 3 ){
  nknots <- K - 2
  knots <- orthogonalsplinebasis::expand.knots( seq( range[1], range[2], length.out = nknots ) )
  basis  <-  orthogonalsplinebasis::SplineBasis(knots = knots, order = degree + 1 )
  d_basis <- orthogonalsplinebasis::deriv( basis)
  d2psi <- orthogonalsplinebasis::OuterProdSecondDerivative(basis)
  return( list( basis = basis, d_basis = d_basis, d2psi = d2psi ) )
}

get_psi <- function(basis, grid, center = TRUE ){
  # when centered splines are used, we remove the last basis to ensure linear independence
  psi <- orthogonalsplinebasis::evaluate( basis, grid )
  if(center){
    psi <- scale( psi, scale = FALSE )
    psi <- psi[, 1:(ncol(psi) - 1 ) ]
  }
  return(psi)
}



est_beta <- function( x, psi, P, K, d2psi, lambda = 0.1 ){
  q <- length(P)
  Omega <- do.call(Matrix::bdiag, replicate(q, d2psi, simplify = F) )
  psimat <- do.call( cbind, lapply( P, function(A) A %*% psi ) )
  psi_inverse <- Matrix::t( psimat ) %*% psimat + lambda * Omega
  # check if it is singular:
  if( matrixcalc::is.singular.matrix( as.matrix(psi_inverse) ) ){
    psi_inverse <- psi_inverse + (lambda / 2) * diag(q*K)
  }
  Proj <- Matrix::solve( psi_inverse, Matrix::t( psimat ) )
  B <- Proj %*% x
  return( lapply( 1:q, function(l) B[ (l-1) * K + 1:K, ] ) )
}


get_degree_of_freedom <- function( psi, P, K, d2psi, lambda = 0.1 ){
  q <- length(P)
  Omega <- do.call(Matrix::bdiag, replicate(q, d2psi, simplify = F) )
  psimat <- do.call( cbind, lapply( P, function(A) A %*% psi ) )
  psi_inverse <- Matrix::t( psimat ) %*% psimat + lambda * Omega
  # check if it is singular:
  if( matrixcalc::is.singular.matrix( as.matrix(psi_inverse) ) ){
    psi_inverse <- psi_inverse + (lambda / 2) * diag(q*K)
  }
  Proj <- Matrix::solve( psi_inverse, Matrix::t( psimat ) )
  S <- psimat %*% Proj
  return( sum(diag(as.matrix(S) ) ) )
}





solve_AssignmentProblem <- function( A, B, method = "Hungarian", fastProd = FALSE ){
  get_cost <- function( A, B ){
    Da <- rowSums( A^2 )
    Db <- rowSums( B^2 )
    if( fastProd ){
      t( t( - 2 * Rfast::mat.mult( A, t(B) ) + Da  ) + Db )
    } else{
      t( t( - 2 * A %*% t( B ) + Da  ) + Db )
    }
  }
  A <- as.matrix(A) # check that it is indeed numeric
  B <- as.matrix(B)
  cost <- get_cost( A, B )
  if( method == "Hungarian" ){
    return( clue::solve_LSAP( as.matrix( cost ) ) )
  }
  if( method == "Greedy" ){
    return( rank( apply( t( cost ), 2, which.min ), ties.method = "random" ) )
  }
}




get_initial_value <- function( x, q ){
  uniformize <- function( x ){ rank( x, ties.method = "random" ) / ( length(x) + 1 ) }
  pca <- irlba::prcomp_irlba( x, n = q, scale = F, center = F )
  if( q == 1 ){ matrix( uniformize( pca$x[,1] / pca$sdev[1] )  )
  } else{ apply( pca$x[, 1:q] %*% solve( diag( pca$sdev[1:q] ) ), 2, uniformize ) }
}



get_permuationMatrix <- function( ord ) Matrix::t( Matrix::sparseMatrix( seq_along( ord ), ord, x = 1 ) )

get_R <- function( l, x, P, G ){
  q <- length(P)
  if( q > 1 ){
    R <- x - Reduce( '+', lapply( (1:q)[-l], function(k) P[[k]] %*% G[[k]] ) )
  }else{ R <- x } 
}

est_P <- function(x, P, G, method, fastProd, random_p){
  q <- length(P)
  for( l in 1:q ){ P[[l]] <- est_P_l( l, x, P, G, method, fastProd, random_p ) }
  return(P)
}

est_P_l <- function( l, x, P, G, method, fastProd, random_p ){
  q <- length(P)
  R <- get_R( l, x, P, G )
  if( is.null( random_p ) ){
    ord <- solve_AssignmentProblem( R, G[[l]], method, fastProd )      
  }else{
    ind <- sort( sample(1:ncol(R), random_p ) )
    ord <- solve_AssignmentProblem( R[, ind ], G[[l]][, ind ], method, fastProd )
  }
  Matrix::t( get_permuationMatrix( ord ) )
}