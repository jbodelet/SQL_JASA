# Make sure your working directory is the data_analysis/
library(dplyr)
library(tidyverse)
library(purrr)
library(furrr)
library(ggplot2)
library(patchwork)
library(irlba)  # PCA for large data
library(e1071)
library(cowplot)
library(gridGraphics)
library(cvTools)
library(gtools)
library(here)
source(here("data_analysis", "src", "estimation.R"))
source(here("data_analysis", "src", "dataAnalysis.R"))



#=================
# 1) Load DATA:
#=================

# download data:
dataset_url <- "https://archive.ics.uci.edu/static/public/401/gene+expression+cancer+rna+seq.zip"
zip_file <- "gene_expression_cancer_rna_seq.zip"
download.file(dataset_url, destfile = zip_file, mode = "wb")
unzip(zip_file, exdir = "gene_expression_data")
tar_gz_file <- list.files("gene_expression_data", pattern = "\\.tar\\.gz$", full.names = TRUE)
untar(tar_gz_file[1], exdir = ".")  # Extract inside the working directory


# extract gene data:
folder <- "TCGA-PANCAN-HiSeq-801x20531/"
raw <- read.csv2(file = paste0(folder, "data.csv" ), sep = ",")
raw$gene_0 <- NULL
gdata <- apply(raw[, -1], 2, as.numeric )
gdata <- scale( gdata )  # rescale data
gdata <- gdata[, complete.cases( t(gdata) ) ]  # remove genes with missing values

# phenotype:
pheno <- read.csv2(file = paste0(folder, "labels.csv" ), sep = ",")
pheno <- as.factor( pheno[, - 1] )



#================================================
# 2) Estimate the Additive Latent Variable Model: 
#================================================

if(1){
  lambda <- 0.49 # The obtained tuning parameter from our cross-validations
}else{ # Running this cross-validation is time-consuming
  library(cvTools)
  source(here("data_analysis", "src", "crossValidation.R"))
  ncores <- 4 # number of cores to use
  lambdaGrid <- seq(0.001, 8, l = 16 )
  gcv <- get_GCV( gdata, q = 2, K = 12, Niter = 20, ncores = 35 )
  plot(gcv$errors ~ gcv$lambdaGrid )
  lambda <- gcv$lambda_min2
}


set.seed(11)
fit <- AFM( gdata, q = 2, K = 12, lambda = lambda, method = "Greedy", Niter = 20, fastProd = TRUE )



#==================================
# 3) PLOT OF THE LATENT SPACE
#==================================


pca <- irlba::prcomp_irlba(gdata, n= 2) # PCA

data2plot <- data.frame( pheno, sql = fit$factor_z, pca$x )


# Plot results:
custom_colors <- c( "#3498DB", "#E67E22", "#2ECC71", "#E74C3C", "#9B59B6" )

sql_plot <- create_composite_plot(data2plot, "sql.1", "sql.2", "First latent", "Second latent")

pca_plot <- create_composite_plot(data2plot, "PC1", "PC2", "First component", "Second component")

plot_latent_space <- (sql_plot | pca_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


legend_sql <- ggplot() +
  geom_text(aes(x = 0.5, y = 0.5, label = "(a) SQL"), size = 5) +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0))

legend_pca <- ggplot() +
  geom_text(aes(x = 0.5, y = 0.5, label = "(b) PCA"), size = 5) +
  theme_void() +
  theme(plot.margin = margin(0,0,0,0))

legend_row <- legend_sql + legend_pca + plot_layout(ncol = 2)

final_plot <- plot_latent_space / legend_row + plot_layout(heights = c(10, 1))

final_plot

ggsave(here("data_analysis", "latent_space.png"), plot = final_plot, width = 8, height = 4, dpi = 300, units = "in")



#========================
# 4) Explained variance: 
#========================

qvec <- 1:20
for( q in qvec ){ # This may be time-consuming
  temp <- AFM( gdata, q = q, K = 12, lambda = 1, method = "Greedy", Niter = 30, fastProd = TRUE )
  filename <- paste0("fit_q", q, ".RDS")
  saveRDS(temp, file = here( "data_analysis", "output", filename ) )
}
totalVar <- mean( (gdata - mean(gdata))^2 )
mse <- sapply( 1:20, function(q) min( readRDS( here( "data_analysis", "output", paste0("fit_q", q, ".RDS") ) )$mse ) )
ev <- 1 - mse / totalVar

# PCA:
pca <- prcomp(gdata)
eigs <- pca$sd^2
ev_pca <- cumsum(eigs[1:20])/ sum(eigs)


# EXPLAINED VARIANCE PLOT
#========================================================================================
data.frame(Factor = 1:20, sql = ev, pca = ev_pca )  %>% 
  ggplot() +
  geom_line(aes(x = Factor, y = sql, color = "SQL"), size = 1 ) +
  geom_point(aes(x = Factor, y = sql, color = "SQL"), size = 3) +
  geom_line(aes(x = Factor, y = pca, color = "PCA"), size = 1 ) +
  geom_point(aes(x = Factor, y = pca, color = "PCA"), size = 3) +
  labs(x = "Number of Factors",
       y = "Cumulative Explained Variance") +
  theme_minimal() +
  theme(panel.grid.major = element_line(linetype = "dotted")) +
  scale_color_manual(name = "Model", values = c("SQL" = "blue", "PCA" = "red")) +
  scale_y_continuous(limits = c(0, 0.9))
#========================================================================================





#=====================
# 5) Classification:
#=====================




# use fit obtained in step 2
factor_z <- fit$factor_z

pca <- irlba::prcomp_irlba(gdata, n= 2) # PCA
pcs <- pca$x

# Training set:
index <- 1:nrow(gdata)
N <- trunc(length(index)/5)
set.seed(42)
accuracy_sql <- replicate(100, get_accuracy(factor_z, pheno, sample(index, N) ) )
accuracy_pca <- replicate(100, get_accuracy(pcs, pheno, sample(index, N) ) )

summary(accuracy_sql); sd(accuracy_sql)
summary(accuracy_pca); sd(accuracy_pca)





#=============================
# 6) FUNCTIONAL CLUSTERING:
#=============================

source(here("data_analysis", "src", "functional_clustering.R"))


# Get important genes from the fit object obtained at step 2:
gene_functions <- get_most_important_genes(fit$g_eval[[1]], nb = 3000 )

# functional clustering:
set.seed(42)
clusters <- get_groups( gene_functions, ng = 10 )
grid <- qnorm(1:200 / 201 )
clusteredFunctions <- list(G = gene_functions, clusters = clusters$groups, grid = grid)

# Save the JSON string to a file:
clusteredFunctionsJSON <- jsonlite::toJSON(clusteredFunctions)
write(clusteredFunctionsJSON, file = "./clusteredFunctions_10.json" )

# save latent space:
write.csv(data.frame(pheno, sql = fit$factor_z ), file = "cancerDataResults.csv")

# Figure 1 in the Supplementary Material can be obtained by running figures_functional_clusters.py
# conda activate tensorflow
# python figures_functional_clusters.py
# the code will produce a pdf called "plot_functionalClusters_10.pdf"


