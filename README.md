# Statistical Quantile Learning for Large Additive Latent Variable Models

This repository contains all the code to reproduce the simulation results and real data analysis from the paper "Statistical Quantile Learning for Large Additive Latent Variable Models".
Our main goal in maintaining this repository is to provide a transparent account of the methods used to generate our empirical results. The software in its current form is not intended for practical analytical use. In contrast, we recommend to use the most recent version of SQL which is available and maintained as a package at 👉[https://github.com/jbodelet/SQL](https://github.com/jbodelet/SQL).

## Organisation

The repository is organised as follows:

### data_analysis

The file cancer_rna.R contains the code to reproduce all the results for the data analysis in Section 5.
The data is automatically downloaded and unzipped from [https://archive.ics.uci.edu](https://archive.ics.uci.edu).
You will need the following R packages:
```
dplyr
tidyverse
ggplot2
Matrix
irlba
cvTools
e1071
cowplot
gridGraphics
```


### simulationStudy
This folder contains the code to reproduce the simulation study in Section 4 of the manuscript.
Running this script on a desktop or laptop is very time-consuming, even with parallel processing.
To run the simulations it is important to keep the folder structure.

Simulations are performed with the following steps.

1. Go to the R folder and run the file `simulate_data.R`. This will generate the data for the simulations and save it in the `simulations` folder.

2. Still in the R folder, run the files `fit_M1.R`, `fit_M2_q1.R`, `fit_M2_q3.R`. This will fit the different methods on the simulated data and save the results in the `simulations` folder.

3. Go to the python folder and train the VAE. The specific steps are detailed in `python/README.md`.

4. Go back to the R folder. Run the file `simulation_table_M1.R` to compute Table 1 and `simulation_table_M2.R` to compute Table 2.

5. Create Figure 1 using the file `figure1.R`.

You will need the following R packages to run the simulations:
```
purrr
furrr
tidyr
dplyr (1.1.4)
ggplot2
gtools
cvTools
```


### simulationStudy_Supplementary

The main.R file contains the code to conduct the additional simulations in Section F of the Supplementary Materials.
You will need the following R packages to run the simulations:
```
RGAN
torch
transport
MonteCarlo
parallel
future
reshape
ggplot2
sql
Rfast
patchwork
```










