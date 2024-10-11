
# Multivariate Generalized Estimating Equations (GEE) for Binary Phenotypes
# This code simulates mixed phenotype data and performs GEE analysis for multivariate tests
# in marginal genetic (G), gene-environment interaction (GE), and Joint gene and gene-environment interaction (GGE) models.

# Required libraries
library(data.table)
library(broom)
library(MASS)
library(dplyr)

# Function to simulate the full dataset with binary phenotypes
# Inputs: 
#   n - number of individuals
#   p - minor allele frequency for the genetic variable (G)
#   f - frequency for the environmental variable (E)
#   B - matrix of regression coefficients
#   mu - mean vector for multivariate normal distribution
#   S - covariance matrix for multivariate normal distribution


sim.dataset.full <- function(n, p, f, B, mu, S) {
  
  # Simulate genotype (G) and environment (E), and compute their interaction (GE)
  simulated_data2 <- data.frame(
    G = rbinom(n, 2, p),  # Simulate genotype as binomial distribution (0, 1, or 2 copies of the minor allele)
    E = rbinom(n, 1, f)   # Simulate environmental exposure as binary (0 or 1)
  ) %>%
    mutate(GE = G * E)    # Interaction term
  
  # Convert to matrix format for further analysis
  X <- as.matrix(simulated_data2)
  
  # Generate phenotype values using linear model plus noise
  YY = X %*% B + mvrnorm(n, mu = mu, Sigma = S)
  colnames(YY) <- c("Y1", "Y2")
  
  # Binarize the continuous outcomes Y1 and Y2 based on their 80th percentiles
  simulated_data.df <- cbind(YY, simulated_data2)
  Y1.cutoff = quantile(simulated_data.df$Y1, 0.8)
  Y2.cutoff = quantile(simulated_data.df$Y2, 0.8)
  simulated_data.df$Y1 <- as.numeric(simulated_data.df$Y1 > Y1.cutoff)
  simulated_data.df$Y2 <- as.numeric(simulated_data.df$Y2 > Y2.cutoff)
  
  return(simulated_data.df)
}

# Note
#mutate(G = ifelse(G >= 1, 1, 0), GE = G * E)  # Dominant model: G is 1 if genotype is 1 or 2.
#mutate(G = ifelse(G == 2, 1, 0),GE = G * E)  # Recessive model: G is 1 if genotype is 2.


# GEE function for the full model, including GE and GGE tests
# Inputs:
#   n - number of individuals
#   response_mat - matrix of binary phenotypes (Y1 and Y2)
#   design_mat - design matrix including genetic and environmental variables
#   beta - initial coefficients for the model
#   r - correlation between phenotypes
#   phi - dispersion parameter for the model


GEE <- function(n, response_mat, design_mat, beta, r, phi) {
  # Initialize response variables and design matrix
  Y1 <- response_mat[, 1]
  Y2 <- response_mat[, 2]
  X <- as.matrix(design_mat)
  
  # Small regularization value to ensure numerical stability
  epsi <- 1e-5
  
  # Initialize summation matrices for score functions, information matrix, etc.
  B1 <- matrix(0, 8, 1)
  B2 <- matrix(0, 8, 8)
  B3 <- matrix(0, 8, 8)
  B4 <- matrix(0, 8, 8)
  
  for (i in 1:n) {
    # Compute probabilities and variances for Y1 and Y2 using logistic regression
    mu1 <- as.numeric(exp(X[i, ] %*% beta[1:4])) / (1 + as.numeric(exp(X[i, ] %*% beta[1:4]))) # Added as.numeric() to ensure scalar
    var1<-  mu1*(1- mu1)
    mu2 <- as.numeric(exp(X[i, ] %*% beta[5:8])) / (1 + as.numeric(exp(X[i, ] %*% beta[5:8])))
    var2 <- mu2 * (1 - mu2)
    
    # Construct design and variance matrices for the GEE
    Di <- rbind(
      c(var1*X[i, ], rep(0, 4)),
      c(rep(0, 4), var2 * X[i, ])
    )
    
    Vi<-matrix(c(var1,r*sqrt(var1*var2),r*sqrt(var1*var2),var2)
               ,2,2,byrow = TRUE)
    V_inv <- solve(Vi + epsi * diag(2))
    
    Yi_mui <- c(Y1[i] - mu1, Y2[i] - mu2)

    # Update summation matrices    
    Ui <- t(Di) %*% V_inv %*% Yi_mui
    B1 <- B1 + Ui
    B2 <- B2 + Ui %*% t(Ui)
    B3 <- B3 + t(Di) %*% V_inv %*% Di # Info matrix for Fisher Scoring
    B4 <- B4 + t(Di) %*% V_inv %*% Di # Ustar Calculation for Sandwich estimator 
    # Note: Info Mat and Ustar are same in GEE. In EGEE they are different.
  }
  
  # Compute covariance matrix for the parameter estimates
  VCovMat <- solve(B4 + epsi * diag(8)) %*% B2 %*% solve(B4 + epsi * diag(8))
  
  # Return the results as a list
  return(list(Usum = B1, InfoMat = B3, UUtSum = B2, Ustar = B4, VCovMat = VCovMat))  
  
}


# GEE function for marginal G tests (testing genetic effects without interactions)
GEE.G <- function(n, response_mat, design_mat, beta.G, r, phi) {

    # Initialize response variables and design matrix
  Y1 <- response_mat[, 1]
  Y2 <- response_mat[, 2]
  X <- as.matrix(design_mat)
  
  # Regularization value to prevent singular matrices
  epsi <- 1e-5
  
  # Initialize matrices for score functions and information
  B1 <- matrix(0, 4, 1)
  B2 <- matrix(0, 4, 4)
  B3 <- matrix(0, 4, 4)
  B4 <- matrix(0, 4, 4)
  
  for (i in 1:n) {  
    # Calculate probabilities and variances for Y1 and Y2
    mu1 <- as.numeric(exp(X[i, ] %*% beta.G[1:2])) / (1 + as.numeric(exp(X[i, ] %*% beta.G[1:2]))) # Added as.numeric() to ensure scalar
    var1<-  mu1*(1- mu1)
    mu2 <- as.numeric(exp(X[i, ] %*% beta.G[3:4])) / (1 + as.numeric(exp(X[i, ] %*% beta.G[3:4])))
    var2 <- mu2 * (1 - mu2)

    # Build design and variance matrices    
    Di <- rbind(
      c(var1*X[i, ], rep(0, 2)),
      c(rep(0, 2), var2 * X[i, ])
    )
    
    Vi<-matrix(c(var1,r*sqrt(var1*var2),r*sqrt(var1*var2),var2)
               ,2,2,byrow = TRUE)
    V_inv <- solve(Vi + epsi * diag(2))
    
    Yi_mui <- c(Y1[i] - mu1, Y2[i] - mu2)
    
    # Update summation matrices
    Ui <- t(Di) %*% V_inv %*% Yi_mui
    B1 <- B1 + Ui
    B2 <- B2 + Ui %*% t(Ui)
    B3 <- B3 + t(Di) %*% V_inv %*% Di
    B4 <- B4 + t(Di) %*% V_inv %*% Di
  }

  # Compute covariance matrix for genetic effects
  VCovMat <- solve(B4 + epsi * diag(4)) %*% B2 %*% solve(B4 + epsi * diag(4))
  
  # Return the results as a list
  return(list(Usum = B1, InfoMat = B3, UUtSum = B2, Ustar = B4, VCovMat = VCovMat))  

}



# Power calculation function for simulating multiple datasets
# and assessing the statistical power of G, GE, and GGE tests
# Inputs:
#   n.sims - number of simulations
#   n - number of individuals
#   p - minor allele frequency
#   f - environmental exposure frequency
#   B - regression coefficients
#   mu - mean vector
#   S - covariance matrix
#   alpha - significance level
power.fn <- function(n.sims, n, p, f, B, mu, S, alpha) {
  
  # Initialize a table to store results from each simulation
  power <- data.table(
    sig.Y1.G = rep(NA, n.sims),
    sig.Y2.G = rep(NA, n.sims),
    sig.Y1_Y2.G = rep(NA, n.sims),
    sig.Y1.GE = rep(NA, n.sims),
    sig.Y2.GE = rep(NA, n.sims),
    sig.Y1_Y2.GE = rep(NA, n.sims),
    sig.Y1.GGE = rep(NA, n.sims),
    sig.Y2.GGE = rep(NA, n.sims),
    sig.Y1_Y2.GGE = rep(NA, n.sims)
  )
  
  # Regularization by adding small value epsilon at the diagonal
  epsi <- 1e-5
  max_iter<-100
  
  for (dataset.i in 1:n.sims) {
    Data <- sim.dataset.full(n, p, f, B, mu, S)
 
    ###################### POWER FOR G-TEST ######################################
    
    pheno<-  Data[,c(1,2)]
    X <- cbind(1, Data$G)
    
    model1 <- glm(Y1 ~ G, data = Data, family = "binomial")
    power$sig.Y1.G[dataset.i] <- tidy(model1)$p.value[2] < alpha
    model2 <- glm(Y2 ~ G, data = Data, family = "binomial") # When using the residual in Y2, it is not 0 1 type and hence was throwing an error.
    power$sig.Y2.G[dataset.i] <- tidy(model2)$p.value[2] < alpha
    
    pearson_resid1 <- residuals(model1, type = "pearson")
    pearson_resid2 <- residuals(model2, type = "pearson")
    
    phi <- (sum(pearson_resid1^2) + sum(pearson_resid2^2)) / (2 * n - 4)
    r <- sum(pearson_resid1 * pearson_resid2) / (phi * (2 * n - 4))
    
    beta.oldGEE.G <- rep(0.1, 4)
    iter <- 0
    
    repeat {
      GEE.G_result <- GEE.G(n, pheno, X, beta.oldGEE.G, r, phi)
      beta.newGEE.G <- beta.oldGEE.G + solve(GEE.G_result$InfoMat + epsi * diag(4)) %*% GEE.G_result$Usum
      
      if (sum((beta.newGEE.G - beta.oldGEE.G)^2) / 4 <= 1e-7 || iter >= max_iter) break
      
      # calculate Pearson residuals for Binary Phenotype1
      pred_prob1<-exp(X%*%beta.newGEE.G[1:2])/ (1+exp(X%*%beta.newGEE.G[1:2]))
      pearson_resid1<-(pheno[1] - pred_prob1) / sqrt(pred_prob1 * (1 - pred_prob1))
      
      # calculate Pearson residuals for Binary Phenotype2
      pred_prob<-exp(X%*%beta.newGEE.G[3:4])/ (1+exp(X%*%beta.newGEE.G[3:4]))
      pearson_resid2 <- (pheno[2] - pred_prob) / sqrt(pred_prob * (1 - pred_prob))
      
      phi<-(sum(pearson_resid1^2)+sum(pearson_resid2^2))/(2 * n - 4)
      r<-(sum(pearson_resid1*pearson_resid2))/(phi * (2 * n - 4))
      
      beta.oldGEE.G <- beta.newGEE.G
      iter <- iter + 1
    }
    
    VarCovMat <- GEE.G_result$VCovMat
    VarCov.gen <- VarCovMat[c(2, 4), c(2, 4)]
    
    beta.G <- beta.newGEE.G[c(2, 4)]
    Wald.G <- t(beta.G) %*% solve(VarCov.gen + epsi * diag(2)) %*% beta.G
    power$sig.Y1_Y2.G[dataset.i] <- pchisq(Wald.G, df = 2, lower.tail = FALSE) < alpha
    
    ###################### POWER FOR GE and GGE-TESTS ######################################
    #pheno<-  Data[,c(1,2)]
    X <- as.matrix(cbind(1, Data[, -c(1, 2)]))
    
    #Y1
    model1 <- glm(Y1 ~ G + E + GE, data = Data, family = "binomial")
    model.updated <- update(model1, . ~ . - GE)
    anova.fit.Y1 <- anova(model1, model.updated, test = "Chisq")
    power$sig.Y1.GE[dataset.i] <- anova.fit.Y1$`Pr(>Chi)`[2] < alpha
    
    # For GGE
    model.updated <- update(model1, . ~ . - G - GE)
    anova.fit.Y1 <- anova(model1, model.updated, test = "Chisq")
    power$sig.Y1.GGE[dataset.i] <- anova.fit.Y1$`Pr(>Chi)`[2] < alpha
    
    # Y2
    model2 <- glm(Y2 ~ G + E + GE, data = Data, family = "binomial")
    model.updated <- update(model2, . ~ . - GE)
    anova.fit.Y2 <- anova(model2, model.updated, test = "Chisq")
    power$sig.Y2.GE[dataset.i] <- anova.fit.Y2$`Pr(>Chi)`[2] < alpha
    # GGE
    model.updated <- update(model2, . ~ . - G - GE)
    anova.fit.Y2 <- anova(model2, model.updated, test = "Chisq")
    power$sig.Y2.GGE[dataset.i] <- anova.fit.Y2$`Pr(>Chi)`[2] < alpha
    
    pearson_resid1 <- residuals(model1, type = "pearson")
    pearson_resid2 <- residuals(model2, type = "pearson")
    
    phi <- (sum(pearson_resid1^2) + sum(pearson_resid2^2)) / (2 * n - 8)
    r <- sum(pearson_resid1 * pearson_resid2) / (phi * (2 * n - 8))
    
    beta.oldGEE <- rep(0.1, 8)
    iter <- 0
    
    repeat {
      GEE_result <- GEE(n, pheno, X, beta.oldGEE, r, phi)
      beta.newGEE <- beta.oldGEE + solve(GEE_result$InfoMat + epsi * diag(8)) %*% GEE_result$Usum
      
      if (sum((beta.newGEE - beta.oldGEE)^2) / 8 <= 1e-7 || iter >= max_iter) break
      
      # calculate Pearson residuals for Binary Phenotype1
      pred_prob1<-exp(X%*%beta.newGEE[1:4])/ (1+exp(X%*%beta.newGEE[1:4]))
      pearson_resid1<-(pheno[1] - pred_prob1) / sqrt(pred_prob1 * (1 - pred_prob1))
      
      # calculate Pearson residuals for Binary Phenotype
      pred_prob<-exp(X%*%beta.newGEE[5:8])/ (1+exp(X%*%beta.newGEE[5:8]))
      pearson_resid2 <- (pheno[2] - pred_prob) / sqrt(pred_prob * (1 - pred_prob))
      
      phi<-(sum(pearson_resid1^2)+sum(pearson_resid2^2))/(2 * n - 8) #(N-p), N=N1+N2, n=8
      r<-(sum(pearson_resid1*pearson_resid2))/(phi * (2 * n - 8))
      
      beta.oldGEE <- beta.newGEE
      iter <- iter + 1
    }
    
    VarCovMat <- GEE_result$VCovMat
    VarCov.ge <- VarCovMat[c(4, 8), c(4, 8)]
    
    beta.GE <- beta.newGEE[c(4, 8)]
    Wald.GE <- t(beta.GE) %*% solve(VarCov.ge + epsi * diag(2)) %*% beta.GE
    power$sig.Y1_Y2.GE[dataset.i] <- pchisq(Wald.GE, df = 2, lower.tail = FALSE) < alpha
    
    ###################### POWER FOR GGE-TEST ######################################
 
    VarCov.gge <- VarCovMat[c(2, 4, 6, 8), c(2, 4, 6, 8)]
    
    beta.GGE <- beta.newGEE[c(2, 4, 6, 8)]
    Wald.GGE <- t(beta.GGE) %*% solve(VarCov.gge + epsi * diag(4)) %*% beta.GGE
    power$sig.Y1_Y2.GGE[dataset.i] <- pchisq(Wald.GGE, df = 4, lower.tail = FALSE) < alpha
  }
  
  power.sim <- power[, .(
    power.Y1_G = mean(sig.Y1.G, na.rm = TRUE),
    power.Y2_G = mean(sig.Y2.G, na.rm = TRUE),
    power.Y1_Y2_G = mean(sig.Y1_Y2.G, na.rm = TRUE),
    power.Y1_GE = mean(sig.Y1.GE, na.rm = TRUE),
    power.Y2_GE = mean(sig.Y2.GE, na.rm = TRUE),
    power.Y1_Y2_GE = mean(sig.Y1_Y2.GE, na.rm = TRUE),
    power.Y1_GGE = mean(sig.Y1.GGE, na.rm = TRUE),
    power.Y2_GGE = mean(sig.Y2.GGE, na.rm = TRUE),
    power.Y1_Y2_GGE = mean(sig.Y1_Y2.GGE, na.rm = TRUE)
  )]
  
  return(power.sim)#power)#power)#
  
  
  
}

###############################################################################

t<-Sys.time()
# Example usage
set.seed(12345)
r<-0.3
mu1 <- 0; s1 <- sqrt(0.5)
mu2 <- 0; s2 <- sqrt(0.5)

# Parameters for bivariate normal distribution
mu <- c(mu1,mu2) # Mean
S <- matrix(c(s1^2, s1*s2*r, s1*s2*r, s2^2),
            2) # Covariance matrix


n.sims <- 10#0
n <- 100#00
p <- 0.2
f <- 0.1
B <- matrix(c(0.15, 0.2, 0.4, 0.15, 0.2, 0.4), ncol = 2)
alpha <- 0.05

power.sim <- power.fn(n.sims, n, p, f, B, mu, S, alpha)
print(power.sim)
Sys.time()-t


####### Comparsion ###########
# For the above simulation settings, the MixPhenOldRegularized took 2.963578 mins
# wheareas MixPhenNewOptimized took just 2.076434 mins.

# One problem is still there that sometimes the resulats are giving NA values of in the result data frame.
# This is happening for different parameter setting.

