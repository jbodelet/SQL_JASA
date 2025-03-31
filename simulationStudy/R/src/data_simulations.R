simulate_npfac <- function(n, p, sde = 1){
  factor <- matrix( rnorm(n) )
  funcs <- replicate( p , gen_func(rnorm(8)), simplify = FALSE )
  gfunc_z <- function(z){
    matrix( sapply( funcs, function(g) g(z) ), ncol = p )
  }
  eps <- matrix( rnorm( n * p, sd = sde ), ncol = p )
  x_expected <- gfunc_z(factor)
  x <- x_expected + eps
  SNR <- mean( apply( x_expected, 2, var ) / sde^2 )
  return( list( data = x, factor_z = factor, gfunc_z = gfunc_z, SNR = SNR ) )
}


simulate_afm <- function(n, p, q = 2, sde = 1, funcs = NULL ){
  # generate factors:
  factor <- matrix( rnorm(n*q ), ncol = q )
  factor <- apply(factor, 2, scale)
  factor <- factor %*% solve( expm::sqrtm( var(factor) ) )
  # functions:
  if(is.null(funcs)){
    funcs <- replicate( q, replicate( p , gen_func(rnorm(8)), simplify = FALSE ), simplify = FALSE )    
  }
  gfunc_z <- lapply( 1:q, function(l){
    function(z){
      matrix( sapply( funcs[[l]], function(g) g(z) ), ncol = p )
      }
  })
  # model:
  eps <- matrix( rnorm( n * p, sd = sde ), ncol = p )
  x_expected <- Reduce( '+', lapply( 1:q, function(l) gfunc_z[[l]](factor[,l]) ) )
  x <- x_expected + eps
  SNR <- mean( apply( x_expected, 2, var ) / sde^2 )
  return( list( data = x, factor_z = factor, gfunc_z = gfunc_z, SNR = SNR ) )
}



gen_func <- function(theta){
  theta <- theta / mean( theta^2 )
  m <- length(theta) / 2
  ind <- 1:m
  a <- theta[ind] / ind 
  b <- theta[ m + ind] / ind 
  f <- Vectorize( function(z){ sum( a * cos( 2 * pi * ind * z / 8 ) + b * sin( 2 * pi * ind * z / 8 ) ) / 2  } )
  grid <- qnorm(1:200/201)
  mean_f <- mean( f(grid) )
  return( function(z) f(z) - mean_f )
}

plotfac <- function(factor1, factor2){
  par( mfrow = c(1,1) )
  plot(factor1, type = "l")
  lines(factor2, col = "red")  
}

getSign <- function(fac1, fac2){
  sign( sum( fac1 * fac2 ) )
}

switchSign <- function(gfunc_z, s ){
  function(z) gfunc_z(s * z)
}



plot_generators <- function(..., ylim = c(-3.5, 3.5) ){
  gen_list <- list(...)
  nb_gen <- length(gen_list)
  p <- ncol( gen_list[[1]](1) )
  rand_index <- sample( 1:p, min(9,p) )
  zz <- qnorm(1:200 /201)
  colors <- c("gray", "black", "red", "blue", "purple", "orange")
  par( mfrow = c( 3, 3) )
  for( j in rand_index ){
    plot( gen_list[[1]](zz)[, j] ~ zz, lwd = 4, col = "gray", ylim = ylim )
    if( nb_gen > 1 ){
      lapply( 2:nb_gen, function(k){
        lines( gen_list[[k]]( zz )[, j] ~ zz, lwd = 2, col = colors[k] ) 
      })
    }
  }
}

rotate <- function(ghat, zhat, z){
  q <- length(ghat)
  # permute:
  per <- gtools::permutations(q, q)
  per_star <- per[ which.max( apply( per, 1, function(u) sum( abs( diag( cor( z, zhat[, u ] ) ) ) ) ) ), ]
  zhat <- zhat[, per_star]
  ghat <- ghat[ per_star ]
  nrow_ghat <- nrow(ghat[[1]])
  # change sign:
  signs <- sign( diag( t( zhat ) %*% z ) )
  zhat <- zhat %*% diag(signs)
  for(l in 1:q){
    if(signs[l] < 0){
      ghat[[l]][, -1] <- ghat[[l]][nrow_ghat:1, -1]
    }
  }
  return(list(ghat = ghat, zhat = zhat))
}
