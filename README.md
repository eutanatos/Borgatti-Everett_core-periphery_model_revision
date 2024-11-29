# Revising the Borgatti-Everett Core-Periphery Model

[![R language](https://img.shields.io/badge/language-R-blue)](https://www.r-project.org/)
[![NordForsk funding](https://img.shields.io/badge/funding-NordForsk-green)](https://www.nordforsk.org/projects/network-dynamics-ethnic-integration)
[![Research Council of Finland funding](https://img.shields.io/badge/funding-Research_Council_of_Finland-green)](https://research.fi/en/results/funding/81442)

This repository contains the data and code to replicate the results of the article:
- Estévez, J. L., & Nordlund, C. (2024). 'Revising the Borgatti-Everett Core-Periphery Model: Inter-Categorical Density Blocks and Partially Connected Cores'. _Social Networks_.

## Purpose of the repository

This repository aims to provide transparency and reproducibility for the results reported in the article.
The main results reported in the paper were conducted in [socnet.se](https://socnet.se/), the software developed to implement the modifications presented in the paper.
Here is the [code](https://github.com/joseluisesna/Borgatti-Everett_core-periphery_model_revision/blob/main/socnet%20script.txt) used to produce the same results. The results, in csv format can be found [here](https://github.com/joseluisesna/Borgatti-Everett_core-periphery_model_revision/blob/main/socnet%20output.csv). The analyses of such results, as well as the simulations and further tests were run in R. For that, we wrote similar functions using the core_periphery() function in [NetUtils](https://schochastics.github.io/netUtils/) (0.8.2) as a template.

## Software requirements

The results were obtained using R version 4.4.1 in RStudio 2024.04.02-764. 
The analysis relies on several dependencies, listed below with their respective versions:
- [data.table](https://rdatatable.gitlab.io/data.table/) (1.15.4)
- [tidyverse](https://www.tidyverse.org/) (2.0.0)
- [stringr](https://stringr.tidyverse.org/) (1.5.1)
- [ggplot2](https://ggplot2.tidyverse.org/) (3.5.1)
- [purrr](https://purrr.tidyverse.org/) (1.0.2)
- [ggpubr](https://rpkgs.datanovia.com/ggpubr/) (0.6.0)
- [patchwork](https://patchwork.data-imaginist.com/) (1.2.0)
- [scales](https://scales.r-lib.org/) (1.3.0)
- [igraph](https://r.igraph.org/) (2.0.3)
- [sna](https://cran.r-project.org/web/packages/sna/index.html) (2.7-2)
- [NetUtils](https://schochastics.github.io/netUtils/) (0.8.2)
- [wCorr](https://cran.r-project.org/web/packages/wCorr/index.html) (1.9.8)
- [glmmTMB](https://cran.r-project.org/web/packages/glmmTMB/index.html) (1.1.10)
- [microbenchmark](https://cran.r-project.org/web/packages/microbenchmark/index.html) (1.4.10)
- [ggeffects](https://strengejacke.github.io/ggeffects/) (1.7.1)

## File list

The repository includes four R scripts and a data folder:
- [01_Illustration figures.R](https://github.com/joseluisesna/Borgatti-Everett_core-periphery_model_revision/blob/main/01_Illustration%20figures.R)
- [02_Analyses of SocNet results.R](https://github.com/joseluisesna/Borgatti-Everett_core-periphery_model_revision/blob/main/02_Analyses%20of%20SocNet%20results.R)
- [03_Functions.R](https://github.com/joseluisesna/Borgatti-Everett_core-periphery_model_revision/blob/main/03_Functions.R)
- [04_Simulations.R](https://github.com/joseluisesna/Borgatti-Everett_core-periphery_model_revision/blob/main/04_Simulations.R) 

## Instructions for use

- Ensure you have the required version of R and RStudio installed.
- Install the necessary packages using the specified versions.
- Notice that script [02_Analyses of SocNet results.R](https://github.com/joseluisesna/Borgatti-Everett_core-periphery_model_revision/blob/main/02_Analyses%20of%20SocNet%20results.R) requires loading the data that is in the [networks folder](https://github.com/joseluisesna/Borgatti-Everett_core-periphery_model_revision/tree/main/networks) and the [output from the socnet.se software](https://github.com/joseluisesna/Borgatti-Everett_core-periphery_model_revision/blob/main/socnet%20output.csv)

## Funding

This research was supported by NordForsk through the funding to "The Network Dynamics of Ethnic Integration", project number [105147](https://www.nordforsk.org/projects/network-dynamics-ethnic-integration). 
José Luis Estévez was also supported by the Research Council of Finland, grant numbers 
[360022](https://research.fi/en/results/funding/81442),
[364382](https://research.fi/en/results/funding/81092), 
[364386](https://research.fi/en/results/funding/81095), and 
[364371](https://research.fi/en/results/funding/81099).

## Citation

- Estévez, J. L., & Nordlund, C. (2024). 'Revising the Borgatti-Everett Core-Periphery Model: Inter-Categorical Density Blocks and Partially Connected Cores'. _Social Networks_.

## Contact information

For any questions, please contact:
- José Luis Estévez (jose.estevez@helsinki.fi)
- Carl Nordlund (carl.nordlund@liu.se)
