# MvGGE Codes (Archived)

The codes from this repository have been converted into an R package **[MvGGE](https://github.com/SauMStats/MvGGE)**, which can now be installed and used directly in R.

[![R-CMD-check](https://github.com/SauMStats/MvGGE/workflows/R-CMD-check/badge.svg)](https://github.com/SauMStats/MvGGE/actions)  
[![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)  
[![Code Size](https://img.shields.io/github/languages/code-size/SauMStats/MvGGE)](https://github.com/SauMStats/MvGGE)  
[![Last Commit](https://img.shields.io/github/last-commit/SauMStats/MvGGE)](https://github.com/SauMStats/MvGGE/commits/main)

## Overview

**MvGGE** is an R package implementing a multi-phenotype approach for joint testing of main genetic effects and gene-environment (GxE) interactions, as described in the paper:

> Mishra, S., & Majumdar, A. (2025).  
> *A Multi-Phenotype Approach to Joint Testing of Main Genetic and Gene-Environment Interaction Effects*.  
> Statistics in Medicine, 44, e70253. [https://doi.org/10.1002/sim.70253](https://doi.org/10.1002/sim.70253)

- **Continuous phenotypes**: Multivariate multiple linear regression (MMLR) with Wilks' Lambda test.
- **Binary or mixed phenotypes**: Generalized estimating equations (GEE) under seemingly unrelated regressions (SUR) with Wald tests.

Key features:
- Unified wrapper function `test_MvGGE()` for automatic phenotype type detection (bivariate support for mixed/binary; multi-phenotype for continuous).
- Simulation tools to generate data under various scenarios (e.g., varying MAF, effect sizes, pleiotropy).
- **_Note_**: _The current version supports bivariate phenotypes only._

This package is ideal for researchers in genetic epidemiology, biostatistics, and related fields analyzing GWAS or similar data.

## Installation

Install the development version from GitHub:

```r
# Install remotes if not already installed
if (!require("remotes")) install.packages("remotes")

remotes::install_github("SauMStats/MvGGE")
```

Load the package:

```r
library(MvGGE)
```

### Dependencies
- R (>= 3.6.0)
- MASS (for multivariate normal simulations)
- stats (base R, for GLM and linear models)

No additional installations are required beyond these.

## Usage

### Quick Start
The core function is `test_MvGGE(Y1, Y2, G, E)`, which automatically detects phenotype types and performs a joint test. 

```r
# Simulate bivariate data (continuous Y1/Y2, binary Yb1/Yb2)
set.seed(123)
sim_data <- simulate_data(
  n       = 1000,            # sample size
  maf     = 0.25,            # minor allele frequency
  f       = 0.2,             # environmental exposure pervalance
  rho     = 0.2,             # correlation between continuous phenotypes (Y1, Y2)
  beta_g  = c(0.1, 0.1),     # genetic effect sizes for Y1 and Y2
  beta_e  = c(0.2, 0.2),     # environmental effect sizes for Y1 and Y2
  beta_ge = c(0.2, 0.2),     # G×E interaction effect sizes for Y1 and Y2
  tau1    = 0.8,             # threshold for binary phenotype Yb1
  tau2    = 0.8              # threshold for binary phenotype Yb2
)


# Test 1: Bivariate Continuous
result_cont <- test_MvGGE(sim_data$Y1, sim_data$Y2, sim_data$G, sim_data$E)
print(result_cont)

# Test 2: Bivariate Binary
result_bin <- test_MvGGE(sim_data$Yb1, sim_data$Yb2, sim_data$G, sim_data$E)
print(result_bin)

# Test 3: Mixed (Continuous Y1, Binary Yb2)
result_mixed <- test_MvGGE(sim_data$Y1, sim_data$Yb2, sim_data$G, sim_data$E)
print(result_mixed)
```

Output example (for mixed case):
```
> summary(result_mixed)
MvGGE Analysis Summary
======================

Analysis: Mixed (Y1 Continuous, Y2 Binary) 
Sample Size: 1000 
Method: GEE-Mixed 
Test Statistic: 20.0868 
P-value: 0.0004801
Significant: Yes  # (at 0.05)

```


## License
This project is licensed under the GPL-3 License. See the [LICENSE](https://github.com/SauMStats/MvGGE/blob/main/LICENSE) file in the main package repository.
