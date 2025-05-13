# Statistical Quantile Learning for Large Additive Latent Variable Models

This repository contains all the code to reproduce the simulation results and real data analysis from the paper "Statistical Quantile Learning for Large Additive Latent Variable Models".
Our main goal in maintaining this repository is to provide a transparent account of the methods used to generate our empirical results. The software in its current form is not intended for practical analytical use. Instead, we recommend using the latest recent version of SQL, which is maintained as a package at [https://github.com/jbodelet/SQL](https://github.com/jbodelet/SQL).


## Data Access

### Description of Data

We use an RNA-Seq (HiSeq) dataset derived from the Pan-Cancer Atlas (PanCanAtlas) Initiative. 
This dataset contains RNA sequencing measurements used to study gene expression across multiple cancer types. 
It includes expression levels for 20,263 genes (with non-null expression) from a total of 801 patients. 
The dataset comprises five tumor types:
* Breast cancer (n = 300)
* Kidney cancer (n = 146)
* Colon cancer (n = 78)
* Lung cancer (n = 141)
* Prostate cancer (n = 136)
    


### Accessing Data

The dataset can be manually downloaded. Go to 
[https://archive.ics.uci.edu/dataset/401/gene+expression+cancer+rna+seq](https://archive.ics.uci.edu/dataset/401/gene+expression+cancer+rna+seq) and click on Download.

For convenience, you can instead run the `cancer_rna.R` script (in the data_analysis folder), which automatically downloads and unzips the dataset.



## Organization

The repository is organized as follows:

### data_analysis

The script `cancer_rna.R` reproduces all the analyses and figures for the Cancer RNA-Seq data as described in Section 5 of the manuscript 
and Section E in the Supplementary Materials.
It handles data download, preprocessing, model estimation, plotting the latent space, computing explained variance,
classification accuracy, and functional clustering of genes.



### simulationStudy
This folder contains the code to reproduce the simulation study in Section 4 of the manuscript.
Running this script on a desktop or laptop is very time-consuming, even with parallel processing.
Please ensure that you keep the folder structure exactly as provided in the repository—including all folder names—since the code relies on these specific names to locate and save data. 
Changing or removing any folder may cause errors.


Simulations are performed with the following steps.

1. Navigate to the R folder and run the file `simulate_data.R`. This will generate the data for the simulations and save it in the `simulations` folder.
Specifically, it creates folders M1, M2_q1, M3_q3 for Model M1, Model M2 with one latent variable, and Model M3 for 3 latent variables respectively.
The simulated data and true parameters (the latent variables and the generators) are saved for each scenario (sample size, for M1, and feature size, for M2).

2. Still in the R folder, run the files `fit_M1.R`, `fit_M2.R`. 
This will fit SQL to all the simulated data and save the results in subfolders (called AFM_fit) contained in `simulations`.

3. Navigate to the python folder and train the VAE. The specific steps are detailed in `python/README.md`.
This will fit VAE to all the simulated data and save the results in subfolders (called vae_fit) contained in `simulations`.

4. Return to the R folder. Run the file `simulation_table_M1.R`.
This will extract the simulated data, the SQL fits, and VAE fits contained in `simulations` and compute Table 1 in the manuscript. 

5. Run the file `simulation_table_M2.R`.
This will extract the simulated data, the SQL fits, and VAE fits contained in `simulations` and compute Table 2 in the manuscript. 

6. Create Figure 1 using the script called `Figure1.R`.



### simulationStudy_Supplementary

The `main.R` file contains the code to conduct the additional simulations presented in Section F of the Supplementary Materials.



## Reproducibility

### Supporting software requirements

All the R analyses were run using the following environment:
- Platform: x86_64-pc-linux-gnu
- Operating System: Ubuntu 22.04.5 LTS
- R version: 4.4.3 (2025-02-28)
- Python version: 3.12.7





### Running time

We indicate the estimated running time (without parallelization) for the following machine:
```
R version 4.3.3 (2024-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 24.04.2 LTS
```

#### data_analysis
- Data analysis: 4 minutes
- Cross-validation: 1 hour
- Explained variance: 2 hours
- Rest of the code: < 1min

#### SimulationStudy

- Simulating data: 10 minutes.
- Running SQL (R code): 5 hours
- Running VAE (Python code): 50 hours
- Rest of the code: < 1min

#### SimulationStudy_Supplementary

- All simulations: 50 hours


### Dependencies

We list below the packages and dependencies used for each folder.

#### data_analysis
```
here (1.0.1)
dplyr (1.1.4) 
tidyverse (2.0.0) 
purrr (1.0.2)
furrr (0.3.1) 
ggplot2 (3.5.1) 
patchwork (1.3.0) 
irlba (2.3.5.1) 
e1071 (1.7.16) 
cowplot (1.1.3) 
gridGraphics (0.5.1) 
gtools (3.9.5) 
cvTools (0.3.3) 
Rfast (2.1.0) 
```

#### SimulationStudy

R packages:

```
here (1.0.1)
dplyr (1.1.4) 
tidyverse (2.0.0) 
purrr (1.0.2)
furrr (0.3.1) 
gtools (3.9.5) 
cvTools (0.3.3) 
irlba (2.3.5.1) 
Matrix (1.4.0)
Rfast (2.1.0)
stargazer (5.2.3)
gridExtra (2.3)
RColorBrewer (1.1.3)
```

Python libraries:

```
python=3.10
numpy=1.26.4
pandas=2.2.3
matplotlib=3.9.1
seaborn=0.13.2
scipy=1.15.2
scikit-learn=1.6.1
jupyter
requests=2.32.3
tensorflow=2.12.1
keras=2.12.0
tensorflow-probability=0.20.0
pip
```


#### SimulationStudy_Supplementary

```
sql (1.0)
here (1.0.1)
future (1.40.0)
ggplot2 (3.5.1) 
Rfast (2.1.0) 
patchwork (1.3.0)
RGAN (0.1.1) 
torch (0.13.0) 
transport (0.15.4) 
MonteCarlo (1.0.6) 
parallel (4.1.2) 
reshape (0.8.9) 
```

To install the sql package, please run:
```
devtools::install_github("jbodelet/SQL/sql@v1.0")
```







