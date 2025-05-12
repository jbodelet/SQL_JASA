library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(here)
setwd(here("simulationStudy", "R"))

PATH_TO_SIM = here( "simulationStudy", "simulations", "M1", "data_n100_var1/")


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


g_true <- load_g_true(PATH_TO_SIM, estimator_name="true")
g_VAE <- load_g_est( paste0(PATH_TO_SIM, "vae_fit/"), g_true=g_true, estimator_name="VAE")
g_SQL <- load_g_est( paste0(PATH_TO_SIM, "AFM_fit/"), g_true=g_true, estimator_name="SQL")

data_true <- g_true %>% pivot_longer(-c(z, sim_number, estimator), names_to = "g", values_to="values")
data_VAE <- g_VAE %>% functional_boxplot_transform()
data_SQL <- g_SQL %>% functional_boxplot_transform()


# define the colors for VAE, SQL, and true
colors = c(brewer.pal(8, "Set1")[c(1,2)], "black")
ggplot(mapping=aes(x=z)) +
  # Plot VAE
  geom_ribbon(aes(ymin=limit_lo, ymax=limit_up, fill="VAE"), alpha=.1, data=data_VAE) +
  geom_ribbon(aes(ymin=q1, ymax=q3, fill="VAE"), alpha=.3, data=data_VAE) +
  
  
  # Plot SQL
  geom_ribbon(aes(ymin=limit_lo, ymax=limit_up, fill="SQL"), alpha=.1, data=data_SQL) +
  geom_ribbon(aes(ymin=q1, ymax=q3, fill="SQL"), alpha=.3, data=data_SQL) +
  
  # Emphasize some lines
  geom_line(aes(y=median, colour="VAE"), data=data_VAE, size=.5)+
  geom_line(aes(y=q1, colour="VAE"), data=data_VAE, size=.2, alpha=.5)+
  geom_line(aes(y=q3, colour="VAE"), data=data_VAE, size=.2, alpha=.5)+
  
  geom_line(aes(y=median, colour="SQL"), data=data_SQL, size=.5)+
  geom_line(aes(y=q1, colour="SQL"), data=data_SQL, size=.2, alpha=.5)+
  geom_line(aes(y=q3, colour="SQL"), data=data_SQL, size=.2, alpha=.5)+
  
  # Plot true line
  geom_line(aes(y=values, colour="Truth"), size=.5,  data=data_true)+
  
  # Graph parameters
  facet_wrap(~factor(g, levels=paste0("g", 1:12))) +
  scale_colour_manual(aesthetics = "colour", 
                      values = c("VAE" = colors[1], "SQL" = colors[2], "Truth" = colors[3]),
                      breaks=c("Truth"), 
                      name="")+
  scale_fill_manual(aesthetics = "fill",
                    values = c("VAE" = colors[1], "SQL" = colors[2], "Truth" = colors[3]),
                    breaks=c("VAE", "SQL"),
                    name= "") +
  xlab("z") +
  ylab("g(z)")+
  theme_bw() +
  theme(text = element_text(size=16),
        legend.spacing.y = unit(-.2, 'cm'),
        legend.margin = margin(0,0,0,0),
        legend.background = element_blank())




