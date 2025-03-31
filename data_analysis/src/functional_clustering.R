get_groups <- function( gfunc, ng ){
  D <- dist(t(gfunc) )
  clust <- hclust( D )
  return( list( groups = cutree( clust, ng ), clust = clust, D = D, ng = ng ) )
}


plot_func_clusters <- function(gfunc2, groups2, grid, ...){
  ng <- max(groups2)
  lapply( 1:ng, function(g){
    matplot( grid, gfunc2[ , groups2 == g ], type = "l", col = "gray", lty = 1, 
             ylab = paste0( "cluster ", g ), yaxis = c(0,0.5, 1), 
             xaxt = "n", yaxt = "n" )
    lines( rowMeans( gfunc2[ , groups2 == g ] ) ~ grid, col = "black", lwd = 2 )
  })
}

get_most_important_genes <- function(G, nb ){
  G[, names( sort(colMeans(G^2), decreasing = T)[1:nb] ) ]
}



