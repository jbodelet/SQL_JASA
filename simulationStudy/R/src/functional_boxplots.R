load_g_est <- function(path, g_true= NULL, estimator_name=NULL) {
  files <- list.files(paste0(path, "g/"))
  g_est <- lapply(seq_along(files), function(sim_i) {
    cat("\nProcessing sim ", sim_i, "/", length(files))
    g <- read_csv(paste0(path, "g/", files[[sim_i]]),
                  show_col_types=FALSE,
                  name_repair="minimal")
    if (!is_null(g_true)) {
      g_true <- g_true %>% select(starts_with("g"))
      g_est <- g %>% select(starts_with("g")) %>% scale(scale = FALSE)
      change_column_sign <- mean( colMeans((g_est - g_true)**2) ) > mean( colMeans((g_est[nrow(g_est):1,] - g_true)**2) )
      if(change_column_sign) g_est <- g_est[nrow(g_est):1, ]
      g <- as_tibble(cbind(g %>% select(z), g_est))
    }
    g %>% mutate(sim_number=sim_i)
    
  })
  g_est <- do.call(rbind, g_est)
  if(!is.null(estimator_name))
    g_est <- g_est %>% mutate(estimator=estimator_name)
  g_est
}


change_sign <- function(g_est, g_true){
  g_est <- g_est %>% mutate(z=round(z,2))
  g_true <- g_true %>% mutate(z=round(z,2))
  
  index <- g_est %>% select(z, sim_number, estimator)
  g_true <- full_join(index, g_true%>%select(-sim_number, -estimator), by="z")
  
  diff_squared <- as_tibble(cbind(index, (g_est %>% select(starts_with("g")) - g_true %>% select(starts_with("g")))**2))
  diff_squared <- diff_squared %>% group_by(sim_number) %>% summarize(across(c(-z, -estimator), function(g)rep(mean(g), length(g))))%>% ungroup() %>% select(-sim_number)
  diff_neg_squared <- as_tibble(cbind(index, (g_est %>% select(starts_with("g")) + g_true %>% select(starts_with("g")))**2))
  diff_neg_squared <- diff_neg_squared %>% group_by(sim_number) %>% summarize(across(c(-z, -estimator), function(g)rep(mean(g), length(g))))%>% ungroup() %>% select(-sim_number)
  change_sign <- (diff_neg_squared < diff_squared) * (-1) + (diff_neg_squared >= diff_squared) * (1)
  
  as_tibble(cbind(index, g_est %>%select(starts_with("g")) * change_sign))
}

load_g_true <- function(path, estimator_name="true"){
  read_csv(paste0(PATH_TO_SIM, "g.csv")) %>% mutate(across(everything(), ~ . - mean(.))) %>% mutate(sim_number=0, estimator=estimator_name)
}



get_boxplot <- function(x, q) {
  stats <- boxplot(x)$stats
  tibble(x=quantile(x, q))
}


functional_boxplot_transform <- function(data) {
  data <- data %>% group_by(z) %>% 
    summarize(across(
      c(sim_number, estimator),
      function(x) x[1]
    ), across(
      -c(sim_number, estimator),
      function(x) boxplot(x, plot=FALSE)$stats
    ))
  data <- data %>% mutate(stats = c("limit_lo", "q1", "median", "q3", "limit_up")) %>% ungroup()
  data <- data %>% pivot_longer(starts_with("g"), names_to = "g", values_to = "values")
  data <- data %>% pivot_wider(names_from="stats", values_from="values")
  data
}
