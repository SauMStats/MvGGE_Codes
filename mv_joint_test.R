# Multivariate Joint GGE test script for continuous phenotypes (lipids).
# Title: A multivariate approach to joint testing of main genetic and gene-environment interaction effects
# Author: Saurabh Mishra, Arunabha Majumdar
#--------------------------------------------------------------------------------------------------------------------#

# Description:
# This script performs real data analysis for evaluating joint genetic (G) and gene-environment interaction (GxE) 
# effects on lipid phenotypes (HDL, LDL, triglycerides) using UK Biobank data. The environmental factor is sleep duration.

# The analysis includes:
# 1. Univariate, bivariate, and trivariate tests for:
#    - Main genetic effects (G effects) on the phenotypes.
#    - Gene-environment interaction effects (GxE effects), incorporating sleep duration.
#    - Joint tests of genetic and GxE effects.
# 2. P-values are calculated for each test scenario, including marginal genetic tests, marginal GxE tests, 
#    and joint genetic and GxE interaction tests.

# The script handles missing genotype data with imputation and standardization, 
# and parallelizes the computations across all SNPs for efficiency.


## Requirements:
# Load required libraries
# Load genotype and phenotype data ( After all the quality control)

# Main analysis steps:
# 1. Extract genotype and phenotype data
# 2. Perform the univariate and multivariate regression on models for G, GE, and GGE tests
# 3. Compile p-values for each hypothesis under comparison for all these models.
# 4. Compare the results of the GGE test with those of other competing tests.

## -----------------------------------------------------------------------------------------##

# Load required packages
required_packages <- c("readr", "data.table", "parallel", "snpStats", "dplyr", "broom")
install_if_missing <- function(p) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p, repos = "http://cran.us.r-project.org")
    library(p, character.only = TRUE)
  }
}
lapply(required_packages, install_if_missing)

# Load the packages
library(snpStats)
library(dplyr)
library(broom)
library(parallel)
library(data.table)
library(reader)


# Load phenotype data
phenotype_data <- fread(file.path(pheno_dir, "Pheno_data.txt"), header = TRUE)

# Function to impute missing genotypes and standardize the matrix
mat_na <- function(G) {
  if (any(is.na(G))) {
    for (i in which(colSums(is.na(G)) != 0)) {
      G[is.na(G[, i]), i] <- mean(G[, i], na.rm = TRUE)
    }
  }
  scale(G)  # Standardizing the matrix
}


# Define file paths for genotype data
ukb_bed_file <- file.path(geno_dir, paste0("geno_file_name.bed"))
ukb_bim_file <- file.path(geno_dir, paste0("geno_file_name.bim"))
ukb_fam_file <- file.path(geno_dir, paste0("geno_file_name.fam"))

# Load genotype data
geno <- read.plink(bed = ukb_bed_file, bim = ukb_bim_file, fam = ukb_fam_file)


# To run the following code, we used the phenotype file, which contains data on lipid phenotypes 'HDL', 'LDL', and 'triglycerides'
# and environmental variable 'sleep duration.'


# GWAS GxE Analysis Function: Generalized and Refactored
gwas_gxe_analysis <- function(snp_index, geno, phenotype_data) { # SNP index is the index for running the analysis for all the SNPs in our study.
  
  # Extract genotypes for the given SNP index
  Geno <- as.data.frame(as(geno$genotypes[, snp_index], "numeric"), row.names = FALSE)
  G <- mat_na(Geno)  # Handle missing data
  E <- phenotype_data$SleepDur
  Y1 <- phenotype_data$HDL
  Y2 <- phenotype_data$LDL
  Y3 <- phenotype_data$TG
  
  # Interaction term between genotype and environment
  GE <- G * E
  Data <- data.frame(Y1, Y2, Y3, G, E, GE)
  DataG <- data.frame(Y1, Y2, Y3, G)
  
  # Helper function: Extract p-values from linear models
  extract_p_value <- function(model, coefficient) {
    summary(model)$coefficients[coefficient, 4]
  }
  
  # Perform G tests (Univariate and Multivariate)
  p_values_G <- sapply(1:3, function(i) extract_p_value(lm(DataG[[i]] ~ G, data = DataG), 2))
  names(p_values_G) <- c("pH_G", "pL_G", "pT_G")

  # Multivariate G model
  model_G <- lm(cbind(Y1, Y2, Y3) ~ G, data = DataG)
  pMV_G <- anova(model_G, update(model_G, . ~ . - G), test = "Wilks")$`Pr(>F)`[2]
  
  # Bivariate G models
  pHL_G <- anova(lm(cbind(Y1, Y2) ~ G, data = DataG), update(lm(cbind(Y1, Y2) ~ G, data = DataG), . ~ . - G), test = "Wilks")$`Pr(>F)`[2]
  pLT_G <- anova(lm(cbind(Y2, Y3) ~ G, data = DataG), update(lm(cbind(Y2, Y3) ~ G, data = DataG), . ~ . - G), test = "Wilks")$`Pr(>F)`[2]
  pHT_G <- anova(lm(cbind(Y1, Y3) ~ G, data = DataG), update(lm(cbind(Y1, Y3) ~ G, data = DataG), . ~ . - G), test = "Wilks")$`Pr(>F)`[2]

  ####################################################################################
  # Perform GE tests
  p_values_GE <- sapply(1:3, function(i) extract_p_value(lm(Data[[i]] ~ G + E + GE, data = Data), 4))
  names(p_values_GE) <- c("pH_GE", "pL_GE", "pT_GE")
  
  model_GE <- lm(cbind(Y1, Y2, Y3) ~ G + E + GE, data = Data)
  pMV_GE <- anova(model_GE, update(model_GE, . ~ . - GE), test = "Wilks")$`Pr(>F)`[2]
  
  # Bivariate GE models
  pHL_GE <- anova(lm(cbind(Y1, Y2) ~ G + E + GE, data = Data), update(lm(cbind(Y1, Y2) ~ G + E + GE, data = Data), . ~ . - GE), test = "Wilks")$`Pr(>F)`[2]
  pLT_GE <- anova(lm(cbind(Y2, Y3) ~ G + E + GE, data = Data), update(lm(cbind(Y2, Y3) ~ G + E + GE, data = Data), . ~ . - GE), test = "Wilks")$`Pr(>F)`[2]
  pHT_GE <- anova(lm(cbind(Y1, Y3) ~ G + E + GE, data = Data), update(lm(cbind(Y1, Y3) ~ G + E + GE, data = Data), . ~ . - GE), test = "Wilks")$`Pr(>F)`[2]
  
  # Perform GGE tests
  p_values_GGE <- sapply(1:3, function(i) {
    model <- lm(Data[[i]] ~ G + E + GE, data = Data)
    anova(model, update(model, . ~ . - G - GE))$`Pr(>F)`[2]
  })
  names(p_values_GGE) <- c("pH_GGE", "pL_GGE", "pT_GGE")
  
  model_updated_GGE <- update(model_GE, . ~ . - G - GE)
  pMV_GGE <- anova(model_GE, model_updated_GGE, test = "Wilks")$`Pr(>F)`[2]
  
  # Bivariate GGE models
  pHL_GGE <- anova(lm(cbind(Y1, Y2) ~ G + E + GE, data = Data), update(lm(cbind(Y1, Y2) ~ G + E + GE, data = Data), . ~ . - G - GE), test = "Wilks")$`Pr(>F)`[2]
  pLT_GGE <- anova(lm(cbind(Y2, Y3) ~ G + E + GE, data = Data), update(lm(cbind(Y2, Y3) ~ G + E + GE, data = Data), . ~ . - G - GE), test = "Wilks")$`Pr(>F)`[2]
  pHT_GGE <- anova(lm(cbind(Y1, Y3) ~ G + E + GE, data = Data), update(lm(cbind(Y1, Y3) ~ G + E + GE, data = Data), . ~ . - G - GE), test = "Wilks")$`Pr(>F)`[2]
  
  # Compile results
  results <- data.frame(
    Chr = geno$map$chromosome[snp_index],
    BP = geno$map$position[snp_index],
    SNP = geno$map$snp.name[snp_index],
    pH_G = p_values_G["pH_G"], pL_G = p_values_G["pL_G"], pT_G = p_values_G["pT_G"], pMV_G = pMV_G,
    pHL_G = pHL_G, pHT_G = pHT_G, pLT_G = pLT_G,
    pH_GE = p_values_GE["pH_GE"], pL_GE = p_values_GE["pL_GE"], pT_GE = p_values_GE["pT_GE"], pMV_GE = pMV_GE,
    pHL_GE = pHL_GE, pHT_GE = pHT_GE, pLT_GE = pLT_GE,
    pH_GGE = p_values_GGE["pH_GGE"], pL_GGE = p_values_GGE["pL_GGE"], pT_GGE = p_values_GGE["pT_GGE"], pMV_GGE = pMV_GGE,
    pHL_GGE = pHL_GGE, pHT_GGE = pHT_GGE, pLT_GGE = pLT_GGE
  )
  
  return(results)
}


## -------------------------------------------------------------------------------------------------------------- ##

  ## Result's columns contain:
                        
# Chr: Chromosome number for the SNP.
# BP: Base pair position of the SNP on the chromosome.
# SNP: Name or identifier of the SNP.

# pH_G: P-value for genetic effect on HDL.
# pL_G: P-value for genetic effect on LDL.
# pT_G: P-value for genetic effect on triglycerides.

# pMV_G: P-value for multivariate genetic test across all phenotypes.

# pHL_G: P-value for genetic effect on HDL and LDL.
# pHT_G: P-value for genetic effect on HDL and triglycerides.
# pLT_G: P-value for genetic effect on LDL and triglycerides.

# pH_GE: P-value for GxE effect on HDL.
# pL_GE: P-value for GxE effect on LDL.
# pT_GE: P-value for GxE effect on triglycerides.

# pMV_GE: P-value for multivariate GxE test across all phenotypes.

# pHL_GE: P-value for GxE effect on HDL and LDL.
# pHT_GE: P-value for GxE effect on HDL and triglycerides.
# pLT_GE: P-value for GxE effect on LDL and triglycerides.

# pH_GGE: P-value for joint G and GxE effect on HDL.
# pL_GGE: P-value for joint G and GxE effect on LDL.
# pT_GGE: P-value for joint G and GxE effect on triglycerides.

( Our proposed ones)

# pMV_GGE: P-value for multivariate joint G and GxE test across all phenotypes. 

# pHL_GGE: P-value for joint G and GxE effect on HDL and LDL.
# pHT_GGE: P-value for joint G and GxE effect on HDL and triglycerides.
# pLT_GGE: P-value for joint G and GxE effect on LDL and triglycerides.
         
