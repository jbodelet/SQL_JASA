# Statistical Quantile Learning for Large Additive Latent Variable Models

This repository contains all the code to reproduce the simulation results and real data analysis from the paper "Statistical Quantile Learning for Large Additive Latent Variable Models".
Our main goal in maintaining this repository is to provide a transparent account of the methods used to generate our empirical results. The software in its current form is not intended for practical analytical use. Instead, we recommend using the latest recent version of SQL, which is maintained as a package at [https://github.com/jbodelet/SQL](https://github.com/jbodelet/SQL).

## Organization

The repository is organized as follows:

### data_analysis

The file cancer_rna.R contains the code to reproduce all the results for the data analysis in Section 5.
The data is automatically downloaded and unzipped from [https://archive.ics.uci.edu](https://archive.ics.uci.edu).


### simulationStudy
This folder contains the code to reproduce the simulation study in Section 4 of the manuscript.
Running this script on a desktop or laptop is very time-consuming, even with parallel processing.
Please ensure that you keep the folder structure exactly as provided in the repository—including all folder names—since the code relies on these specific names to locate and save data. 
Changing or removing any folder may cause errors.

Simulations are performed with the following steps.

1. Navigate to the R folder and run the file `simulate_data.R`. This will generate the data for the simulations and save it in the `simulations` folder.

2. Still in the R folder, run the files `fit_M1.R`, `fit_M2.R`. This will fit the different methods on the simulated data and save the results in the `simulations` folder.

3. Navigate to the python folder and train the VAE. The specific steps are detailed in `python/README.md`.

4. Return to the R folder. Run the file `simulation_table_M1.R` to compute Table 1 and `simulation_table_M2.R` to compute Table 2.

5. Create Figure 1 using the file `Figure1.R`.



### simulationStudy_Supplementary

The main.R file contains the code to conduct the additional simulations presented in Section F of the Supplementary Materials.

## Reproducibility

### Supporting software requirements

Version of primary software used: R version 4.1.3 or superior. Python version 3.12.7 or superior.

Libraries and dependencies used by the code R packages: dplyr (1.1.4) furrr (0.3.1) purrr (1.0.2)
ggplot2 (3.5.1) gtools (3.9.5) cvTools (0.3.3) RGAN (0.1.1) torch (0.13.0) transport (0.15.4) MonteCarlo
(1.0.6) parallel (4.1.2) reshape (0.8.9) Rfast (2.1.0) patchwork (1.3.0) tidyverse (2.0.0) irlba (2.3.5.1) e1071
(1.7.16) cowplot (1.1.3) gridGraphics (0.5.1) Matrix (1.4.0)

Python libraries: streamlit tensorflow scipy numpy matplotlib pandas Jupyter seaborn requests scikit-learn

For the additional simulation study (simulationStudy_Supplementary), you need to install the SQL package:
```
devtools::install_github("jbodelet/SQL/sql")
``








