

# Multivariate Joint Test (GGE Test)

## Overview
This repository contains R code used for real data analysis and simulations as part of the paper titled **"A Multivariate Approach to Joint Testing of Main Genetic and Gene-Environment Interaction Effects"**. The main focus of the repository is to demonstrate the power and effectiveness of the joint testing of genetic (G) and gene-environment (GxE) interaction effects on multivariate phenotypes using the **GGE test**. The method is applied to lipid phenotypes using data from the UK Biobank and simulated phenotypes in different scenarios.

### Abstract
Gene-environment (GxE) interactions play a critical role in shaping complex phenotypes. However, due to the often subtle nature of GxE effect sizes, statistical power in GxE studies can be limited. To address this, we propose a multivariate joint testing procedure that enhances power by combining pleiotropic effects of the main genetic and GxE effects on multivariate phenotypes. We apply the method to lipid phenotypes such as HDL, LDL, and triglycerides (TG), with sleep duration as the environmental factor.

### Simulation Studies
The repository also includes simulations for the paper, where three bivariate phenotype cases are demonstrated:
1. **Bivariate Continuous Case** (File: `ContPhenSim.R`): Both phenotypes are continuous.
2. **Mixed Phenotype Case** (File: `MixedPhenSim.R`): One phenotype is continuous, and the other is binary.
3. **Bivariate Binary Case** (File: `BinPhenSim.R`): Both phenotypes are binary.

These simulations compare the GGE test results with univariate and multivariate marginal genetic tests (G test), marginal GxE tests, and univariate joint tests.

## Files
- **`mv_joint_test.R`**: The main R script for the real data analysis, focusing on the joint genetic and GxE effects for lipid phenotypes using UK Biobank data.
  
- **`ContPhenSim.R`**: Simulation for the continuous phenotype case. This script demonstrates the GGE test for two continuous phenotypes, comparing it with alternative tests (univariate and multivariate).

- **`MixedPhenSim.R`**: Simulation for the mixed phenotype case (one continuous, one binary). The script implements the GGE test for this mixed phenotype scenario, comparing it with alternative tests (univariate and multivariate).

- **`BinPhenSim.R`**: Simulation for the binary phenotype case, where both phenotypes are binary, applying the GGE test and comparing it with other competing tests.

- **`README.md`**: The explanation of the repository contents, code functionality, and simulations conducted for both real and simulated data.



## Methodology
The **GGE test** evaluates the joint effect of the main genetic and the GxE interaction across multiple phenotypes. In each simulation, we compare the power of the GGE test against:
- **Univariate Marginal Genetic Tests**
- **Multivariate Marginal Genetic Test (G test)**
- **Univariate GxE Tests (univariate GxE tests)**
- **Multivariate GxE Tests (GE test)**
- **Univariate Joint GxE Tests (gge tests)**

The simulation results help validate the robustness and increased power of the proposed multivariate joint testing approach in comparison to traditional methods.

<!-- ## Details -->
<!-- 1. **Real Data Analysis**:   -->
<!--    - To perform joint testing on the lipid phenotypes from the UK Biobank, run `mv_joint_test.R`. -->
<!--    - The script will compute and compare G, GxE, and GGE effects on lipid traits (HDL, LDL, TG). -->

<!-- 2. **Simulations**:   -->
<!--    - To simulate datasets and assess the performance of the GGE test in different bivariate cases, run the relevant script: -->
<!--      - `ContPhenSim.R` for continuous phenotypes. -->
<!--      - `MixedPhenSim.R` for mixed phenotypes (one continuous, one binary). -->
<!--      - `BinPhenSim.R` for binary phenotypes. -->

<!-- 3. **Dependencies**:   -->
<!--    Ensure the following R packages are installed: -->
<!--    - `data.table` -->
<!--    - `MASS` -->
<!--    - `broom` -->
<!--    - `dplyr` -->

<!-- 4. **Customization**:   -->
<!--    Simulation parameters such as sample size, minor allele frequency (MAF), and the coefficient matrix for G and GxE effects can be customized in the respective scripts. -->


<!-- ## Citation -->
<!-- If you use this code or methodology in your work, please cite the paper:   -->
<!-- **"A Multivariate Approach to Joint Testing of Main Genetic and Gene-Environment Interaction Effects"**. -->

<!-- ## Contact -->
<!-- For any questions or further clarification, feel free to reach out to the author. -->

<!-- --- -->

<!-- This repository is continually updated to include further analysis and code improvements. We welcome feedback and contributions. -->
