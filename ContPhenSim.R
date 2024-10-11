#' Simulate Bivariate Continuous Phenotypes and Run Power Analysis
#'
#' This script simulates a dataset for two continuous phenotypes (Y1 and Y2)
#' using genetic (G), environmental (E), and interaction (GE) terms. It performs
#' linear regression models to assess the effects of G, GE, and GGE interactions,
#' and estimates statistical power for these tests across multiple simulations.
#'
#' @param n.sims Number of simulation iterations.
#' @param n Sample size for each simulation.
#' @param B Coefficient matrix for G, E, and GE effects.
#' @param mu Mean vector for multivariate normal distribution of residuals.
#' @param sigma Covariance matrix for the multivariate normal distribution.
#' @param p Minor allele frequency (MAF) for the G variable (genotype).
#' @param f Probability of exposure for the E variable (environment).
#' @param alpha Significance threshold (e.g., 0.05).
#' @return A data table with power estimates for G, GE, and GGE tests.



## Packages
library(data.table)
library(MASS)
library(broom)
library(dplyr)

############# Functions ################
# Simulate Dataset Function
sim.dataset <- function(n, B, mu, sigma, p, f) {
  # Create a dataset with genotype (G) and environment (E)
  simulated_data <- data.table(
    G = rbinom(n, 2, p),   # Genotype G with binomial distribution (0, 1, 2)
    E = rbinom(n, 1, f)    # Environment E (binary: 0 or 1)
  )
  
  # Interaction term (G * E)
  simulated_data[, GE := G * E]
  
  # Combine effects of G, E, and GE using the coefficient matrix (B)
  XX <- as.matrix(simulated_data)
  
  # Generate phenotypes (Y1 and Y2) as a linear combination of predictors with random noise
  YY <- XX %*% B + mvrnorm(n, mu = mu, Sigma = sigma)
  colnames(YY) <- c("Y1", "Y2")
  
  # Add Y1 and Y2 to the dataset
  simulated_data[, `:=`(Y1 = YY[, 1], Y2 = YY[, 2])]
  
  return(simulated_data)
}


# Power Function
power.fn <- function(n.sims, n, B, mu, sigma, p, f, alpha) {
  results <- vector("list", n.sims)
  
  for (dataset.i in 1:n.sims) {
    simulated_df <- sim.dataset(n, B, mu, sigma, p, f)
    
    simulated_df_G <- simulated_df[, .(Resid_Y1 = residuals(lm(Y1 ~ E, data = simulated_df)),
                                       Resid_Y2 = residuals(lm(Y2 ~ E, data = simulated_df)),
                                       G, E, GE)]
    
    model1 <- tidy(lm(Resid_Y1 ~ G, data = simulated_df_G))
    sig.bg1 <- model1$p.value[2] < alpha
    
    model1 <- tidy(lm(Y1 ~ G + E + GE, data = simulated_df))
    sig.bge1 <- model1$p.value[4] < alpha
    
    model2 <- tidy(lm(Resid_Y2 ~ G, data = simulated_df_G))
    sig.bg2 <- model2$p.value[2] < alpha
    
    model2 <- tidy(lm(Y2 ~ G + E + GE, data = simulated_df))
    sig.bge2 <- model2$p.value[4] < alpha
    
    model1 <- lm(Y1 ~ G + E + GE, data = simulated_df)
    model.updated1 <- update(model1, . ~ . - G - GE)
    anova.fit.Y1 <- anova(model1, model.updated1)
    sig.bgge1 <- anova.fit.Y1$`Pr(>F)`[2] < alpha
    
    model2 <- lm(Y2 ~ G + E + GE, data = simulated_df)
    model.updated2 <- update(model2, . ~ . - G - GE)
    anova.fit.Y2 <- anova(model2, model.updated2)
    sig.bgge2 <- anova.fit.Y2$`Pr(>F)`[2] < alpha
    
    model <- lm(cbind(Resid_Y1, Resid_Y2) ~ G, data = simulated_df_G)
    model.updated.G <- update(model, . ~ . - G)
    anova.fit.G <- anova(model, model.updated.G, test = "Wilks")
    sig.wilks.G <- anova.fit.G$`Pr(>F)`[2] < alpha
    
    model <- lm(cbind(Y1, Y2) ~ G + E + GE, data = simulated_df)
    model.updated.GE <- update(model, . ~ . - GE)
    anova.fit.GE <- anova(model, model.updated.GE, test = "Wilks")
    sig.wilks.GE <- anova.fit.GE$`Pr(>F)`[2] < alpha
    
    model.updated.GGE <- update(model, . ~ . - G - GE)
    anova.fit.GGE <- anova(model, model.updated.GGE, test = "Wilks")
    sig.wilks.GGE <- anova.fit.GGE$`Pr(>F)`[2] < alpha
    
    results[[dataset.i]] <- list(sig.bg1 = sig.bg1, sig.bge1 = sig.bge1, sig.bg2 = sig.bg2, 
                                 sig.bge2 = sig.bge2, sig.bgge1 = sig.bgge1, sig.bgge2 = sig.bgge2, 
                                 sig.wilks.G = sig.wilks.G, sig.wilks.GE = sig.wilks.GE, sig.wilks.GGE = sig.wilks.GGE)
  }
  
  results_dt <- rbindlist(results)
  
  coef.power <- results_dt[, .(
    bg1.power = mean(sig.bg1),
    bg2.power = mean(sig.bg2),
    bge1.power = mean(sig.bge1),
    bge2.power = mean(sig.bge2),
    bg_ge1.power = mean(sig.bgge1),
    bg_ge2.power = mean(sig.bgge2),
    Wilks.power.G = mean(sig.wilks.G),
    Wilks.power.GE = mean(sig.wilks.GE),
    Wilks.power.GGE = mean(sig.wilks.GGE)
  )]
  
  return(coef.power)
}


# Set Seed for Reproducibility
set.seed(12345)
alpha <- 0.05

# Define Simulation Parameters
r <- 0.3
mu1 <- 0; s1 <- sqrt(0.5)
mu2 <- 0; s2 <- sqrt(0.5)

# Parameters for Bivariate Normal Distribution
mu <- c(mu1, mu2) # Mean
S <- matrix(c(s1^2, s1 * s2 * r, s1 * s2 * r, s2^2), ncol = 2) # Covariance Matrix

n.sims <- 1000
n <- 1000
p <- 0.2
f <- 0.1
B <- matrix(c(0.1, 0.15, 0.3, 0.1, 0.15, 0.3), ncol = 2) 

# Calculate Power
power.sim <- power.fn(n.sims, n, B, mu, S, p, f, alpha)
print(power.sim)

# ----------------------------------------------------------------------------------
## Details of results columns:

# bg1.power: Power for detecting the genetic effect (G) on the first phenotype (Y1).
# bg2.power: Power for detecting the genetic effect (G) on the second phenotype (Y2).

# bge1.power: Power for detecting the gene-environment interaction effect (GE) on the first phenotype (Y1).
# bge2.power: Power for detecting the gene-environment interaction effect (GE) on the second phenotype (Y2).
# This measures the ability to detect the significant interaction between G and E in phenotype  

# bg_ge1.power: Power for detecting the combined genetic and interaction effects (GGE) on the first phenotype (Y1).
# bg_ge2.power: Power for detecting the combined genetic and interaction effects (GGE) on the second phenotype (Y2).
# This assesses the joint effect of both G and GE terms on Y1 (or Y2).

# Wilks.power.G: Multivariate power for detecting the joint effect of G on both phenotypes (Y1 and Y2).
# This is based on Wilks' test, which simultaneously assesses the effect of G across both phenotypes.

# Wilks.power.GE: 
# Multivariate power for detecting the gene-environment interaction effect (GE) on both phenotypes (Y1 and Y2).
# This measures the ability to detect GE effects across both Y1 and Y2 using Wilks' test.

# Wilks.power.GGE: ( Our proposed hypothesis)
# Multivariate power for detecting the joint effect of both G and GE (GGE) on the two phenotypes (Y1 and Y2).
# This is the power of the GGE test for detecting the combined genetic and interaction effects across both phenotypes using Wilks' test.

