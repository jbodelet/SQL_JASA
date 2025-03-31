str2vec <- function(str){
  numbers_str <- strsplit(gsub("\\[|\\]", "", str), split = ",")[[1]] # Remove square brackets and split by comma
  numbers <- sapply(numbers_str, function(x) as.numeric(gsub(",", "", x))) # Replace commas and convert to numeric
  names(numbers) <- paste0("M", 1:length(numbers))
  return(numbers)
}

get_explainedVar <- function(x, P, G){
  q <- length(G)
  mse <- rep(0, q)
  for(l in 1:q ){
    pred <- Reduce( '+', lapply( 1:l, function(m) P[[m]] %*% G[[m]] ) )
    mse[l] <- Matrix::mean( ( x - pred )^2 )
  }
  totalVar <- Matrix::mean( ( x - mean(x) )^2 )
  ev <- 1- mse / totalVar
  return(ev)
}

get_accuracy <- function(x, y, testindex){
  dat <- data.frame(y = y, x )
  testset <- dat[testindex,]
  trainset <- dat[-testindex,]
  svm.model <- svm(y ~ ., data = trainset, kernel = "radial", cost = 2, gamma = 1)
  svm.pred <- predict(svm.model, testset[,-1])
  accuracy <- 1 - hamming.distance( as.numeric( svm.pred), as.numeric(testset$y) )/N
  return(accuracy)
}

cross_validate_K <- function(data, Kgrid = 4:10, ncores = 2, nfold = 5, q = 1, 
                             lambda = 0, Niter = 20, random_p = NULL){
  library(cvTools)
  library(purrr)
  library(furrr)
  folds <- cvFolds(NROW(data), K = nfold)
  compute_cv_error <- function(K) {
    mean(sapply(1:nfold, function(fold) {
      trainSet <- data[folds$which != fold, ]
      testSet <- data[folds$which == fold, ]
      fit <- AFM(trainSet, q, lambda = lambda, K = K, Niter, 
                 fastProd = TRUE, method = "Greedy", random_p = random_p)
      factors <- estimate_factors(data, fit$gfunc_z, method = "Greedy", 
                                  Niter = 2)$factor_z
      factors <- as.matrix(factors[folds$which == fold, 
      ])
      pred <- as.matrix(Reduce("+", lapply(1:q, function(l) fit$gfunc_z[[l]](factors[, 
                                                                                     l]))))
      return(mean((testSet - pred)^2))
    }))
  }
  plan(multisession, workers = ncores)
  mse <- furrr::future_map_dbl(Kgrid, compute_cv_error)
  return(list(mse = mse, Kgrid = Kgrid))
}



create_composite_plot <- function(data, x_var, y_var, top_title, y_lab) {
  
  # 1. Top label plot: displays the top title
  top_plot <- ggplot() +
    theme_void() +
    labs(title = top_title) +
    theme(
      plot.title = element_text(size = 11, hjust = 0.5),
      plot.margin  = margin(b = 0)
    )
  
  # 2. Main scatter plot
  # We use aes_string so that x and y can be passed as strings.
  main_plot <- data %>%
    ggplot(aes_string(x = x_var, y = y_var, fill = "pheno")) +
    geom_point(
      alpha = 0.7,
      shape  = 21,       # allows separate fill and outline
      size   = 2.5,
      stroke = 0.5,
      color  = "gray40"  # fixed outline color (so no legend is created from 'color')
    ) +
    labs(x = NULL, y = y_lab) +
    scale_fill_manual(
      values = custom_colors,
      name   = "Pheno",
      guide  = guide_legend(override.aes = list(shape = 21))
    ) +
    theme_minimal() +
    theme(
      panel.grid   = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.position = "right"
    )
  
  # 3. Bottom density plot
  bottom_plot <- data %>%
    ggplot(aes_string(x = x_var, fill = "pheno", color = "pheno")) +
    geom_density(alpha = 0.4, show.legend = FALSE) +  # turn off legend here
    scale_y_reverse() +
    scale_fill_manual(
      values = custom_colors,
      name   = "Pheno"
    ) +
    scale_color_manual(
      values = custom_colors,
      guide  = "none"  # hide color legend
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x       = element_blank(),
      axis.ticks.x      = element_blank(),
      panel.grid        = element_blank(),
      panel.background  = element_blank(),
      plot.background   = element_blank(),
      axis.text.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      axis.title.x      = element_blank(),
      legend.position   = "none",
      plot.margin       = margin(t = -10)
    )
  
  # 4. Stack the three components vertically.
  composite <- top_plot / main_plot / bottom_plot +
    plot_layout(heights = c(0.1, 3, 1))
  
  return(composite)
}




