# Reproduce Figure 1 in the manuscript (Functional Boxplots for Model M1)
# -----------------------------------------------------
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(here)
source( here( "simulationStudy", "R", "src", "functional_boxplots.R") )

PATH_TO_SIM = here( "simulationStudy", "simulations", "M1", "data_n100_var1/")


# Load True and Estimated Functions
g_true <- load_g_true(PATH_TO_SIM, estimator_name="true")
g_VAE <- load_g_est( paste0(PATH_TO_SIM, "vae_fit/"), g_true=g_true, estimator_name="VAE")
g_SQL <- load_g_est( paste0(PATH_TO_SIM, "AFM_fit/"), g_true=g_true, estimator_name="SQL")

# Prepare Data for Plotting
data_true <- g_true %>% pivot_longer(-c(z, sim_number, estimator), names_to = "g", values_to="values")
data_VAE <- g_VAE %>% functional_boxplot_transform()
data_SQL <- g_SQL %>% functional_boxplot_transform()


# define the colors for VAE, SQL, and true
colors = c(brewer.pal(8, "Set1")[c(1,2)], "black")



# Plot Functional Clustering Results
#====================================
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




