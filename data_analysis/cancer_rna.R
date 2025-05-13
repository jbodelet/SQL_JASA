#------------------------------------------------------------
# Script: cancer_rna.R
# Purpose: Master script  for reproducing all the data analysis in Section 5 of the manuscript and Section E in the Supplementary Materials
# Usage:
#   This file contains instructions for reproducing all the results for the data analysis, including numbers and figures contained in the paper. 
#   It handles data download, preprocessing, model estimation, plotting the latent space, computing explained variance,
#   classification accuracy, and functional clustering of genes.
#   It is best practice to open the repository as an Rstudio project (open SQL_JASA.Rproj in top-level directory).
#   Alternatively, set the working directory as the top-level directory and make sure that here package works automatically.
# Running time:
#   Fitting the SQL method with pre-selected lambda takes approximatively 4 minutes runtime.
#   Running the cross-validation is computationally intensive and takes about an hour.
#   Obtaining the explained variances takes about 3 hours.
#------------------------------------------------------------


# Load required packages (see Reproducibility section in README).
library(dplyr)
library(tidyverse)
library(purrr)
library(furrr)
library(ggplot2)
library(patchwork)
library(irlba)
library(e1071)
library(cowplot)
library(gridGraphics)
library(cvTools)
library(gtools)
library(here)

# Load functions necessary for the data analysis:
source(here("data_analysis", "src", "estimation.R"))
source(here("data_analysis", "src", "dataAnalysis.R"))



#====================================================
# 1) Download, unzip, and prepare the RNA-Seq Data
#====================================================

# download data the PanCanAtlas RNA-Seq dataset (801 samples x ~20k genes)
dataset_url <- "https://archive.ics.uci.edu/static/public/401/gene+expression+cancer+rna+seq.zip"
zip_file <- here( "data_analysis", "gene_expression_cancer_rna_seq.zip" )
download.file(dataset_url, destfile = zip_file, mode = "wb")
unzip(zip_file, exdir = here( "data_analysis", "gene_expression_data") )

# Extract the tar.gz archive contained within the unzipped folder
tar_gz_file <- list.files(here( "data_analysis", "gene_expression_data" ), pattern = "\\.tar\\.gz$", full.names = TRUE)
untar(tar_gz_file[1], exdir = here( "data_analysis") )  # Extract inside the working directory

# Read gene expression CSV, remove index column, and convert to numeric matrix
folder <- here( "data_analysis", "TCGA-PANCAN-HiSeq-801x20531/")
raw <- read.csv2(file = paste0(folder, "data.csv" ), sep = ",")
raw$gene_0 <- NULL

# Convert expression values to numeric
gdata <- apply(raw[, -1], 2, as.numeric )

# Rescale data
gdata <- scale( gdata )  

# Remove genes (columns) with missing values
gdata <- gdata[, complete.cases( t(gdata) ) ]

# Load phenotype (tumor type) labels and convert to factor
pheno <- read.csv2(file = paste0(folder, "labels.csv" ), sep = ",")
pheno <- as.factor( pheno[, - 1] )



#================================================
# 2) Estimate the Additive Latent Variable Model: 
#================================================

# The cross-validation takes about an hour on most machines (without parallelization)
# To run the cross-validation, set notrun <- FALSE
# Otherwise, to save time, we provide the tuning parameter (lambda) used in our data analysis
notrun <- TRUE

if(notrun){
  # The obtained tuning parameter from our cross-validations:
  lambda <- 0.49 
}else{ 
  # Run the cross-validations:
  library(cvTools)
  source(here("data_analysis", "src", "crossValidation.R"))
  ncores <- 4 # Note: adjust ncores based on your machine
  lambdaGrid <- seq(0.001, 8, l = 16 )
  gcv <- get_GCV( gdata, q = 2, K = 12, Niter = 20, ncores = 35 )
  plot(gcv$errors ~ gcv$lambdaGrid )
  lambda <- gcv$lambda_min2  # obtained tuning parameter
}

# Fit the SQL method with selected lambda (approx. 4 minutes runtime)
set.seed(11)
fit <- AFM( gdata, q = 2, K = 12, lambda = lambda, method = "Greedy", Niter = 20, fastProd = TRUE )



#============================================
# 3) Reproduce Figure 3 (Latent space)
#============================================
# Approximate run time: < 1min

# Compute PCA with 2 principal component vectors:
pca <- irlba::prcomp_irlba(gdata, n= 2)

# Prepare data frame for plotting
data2plot <- data.frame( pheno, sql = fit$factor_z, pca$x )

# Generate composite scatter plots 
sql_plot <- create_composite_plot(data2plot, "sql.1", "sql.2", "First latent", "Second latent")

pca_plot <- create_composite_plot(data2plot, "PC1", "PC2", "First component", "Second component")

plot_latent_space <- (sql_plot | pca_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


# Arrange plots side by side with shared legend
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

print( final_plot )

# Save final latent space figure
ggsave(here("data_analysis", "latent_space.png"), plot = final_plot, width = 8, height = 4, dpi = 300, units = "in")



#=========================================================
# 4) Compute and Plot Explained Variance for SQL and PCA
#=========================================================
# Approximate run time (without parallelization): 3 hours
qvec <- 1:20

# For each number of factors, estimate SQL and save results
# This may be time-consuming
for( q in qvec ){
  temp <- AFM( gdata, q = q, K = 12, lambda = 1, method = "Greedy", Niter = 30, fastProd = TRUE )
  filename <- paste0("fit_q", q, ".RDS")
  saveRDS(temp, file = here( "data_analysis", "output", filename ) )
}

# For each number of factors, load saved AFM fit and compute MSE
mse <- sapply( 1:20, function(q) min( readRDS( here( "data_analysis", "output", paste0("fit_q", q, ".RDS") ) )$mse ) )

# Compute explained variance for SQL
totalVar <- mean( (gdata - mean(gdata))^2 )
ev <- 1 - mse / totalVar

# Compute explained variance for PCA:
pca <- prcomp(gdata)
eigs <- pca$sd^2
ev_pca <- cumsum(eigs[1:20])/ sum(eigs)



# Reproduce Figure 2 (explained Variance for SQL vs. PCA)
#=========================================================


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






#================================================
# 5) Classification Accuracy for SQL vs. PCA:
#================================================
# Approximate run time: < 1min

# use fit obtained in step 2
factor_z <- fit$factor_z

# Get PCA scores:
pca <- irlba::prcomp_irlba(gdata, n= 2)
pcs <- pca$x

# Compare 2-dimensional representations via repeated random train/test splits
index <- 1:nrow(gdata)
N <- trunc(length(index)/5)
set.seed(42)
accuracy_sql <- replicate(100, get_accuracy(factor_z, pheno, sample(index, N) ) )
accuracy_pca <- replicate(100, get_accuracy(pcs, pheno, sample(index, N) ) )

# Summarize and print accuracy results
summary(accuracy_sql); sd(accuracy_sql)
summary(accuracy_pca); sd(accuracy_pca)





#=====================================================================================================
# 6) Perform Functional Clustering of Top Genes and Reproduce Figure 1 in the Supplementary Materials
#=====================================================================================================
# Reproducing Figure 1 in Supplementary Materials requires performing the four step below (three with R and one with Python)
# It requires as input the clustered functions (clusteredFunctions_10.json) 
# and the obtained latent space (cancerDataResults.csv)
# Approximate run time: < 1min

source(here("data_analysis", "src", "functional_clustering.R"))

# 1. Identify top 3000 genes by importance from SQL fit
gene_functions <- get_most_important_genes(fit$g_eval[[1]], nb = 3000 )

# 2. Cluster gene functions into 10 groups
set.seed(42)
clusters <- get_groups( gene_functions, ng = 10 )
grid <- qnorm(1:200 / 201 )
clusteredFunctions <- list(G = gene_functions, clusters = clusters$groups, grid = grid)

# 3. Save Results:
clusteredFunctionsJSON <- jsonlite::toJSON(clusteredFunctions)
write(clusteredFunctionsJSON, file = here("data_analysis","clusteredFunctions_10.json" ) ) # clusters
write.csv(data.frame(pheno, sql = fit$factor_z ), file = here("data_analysis","cancerDataResults.csv") ) # latent space

# 4. Run Python file:
# # Figure 1 in the Supplementary Material can be obtained by running figures_functional_clusters.py
# # Once the files clusteredFunctions_10.json and cancerDataResults.csv are saved, 
# # run the following commands in a terminal to produce Figure 1 in the Supplementary Materials
# conda activate tensorflow
# python figures_functional_clusters.py
# # the code will produce a pdf called "plot_functionalClusters_10.pdf"


