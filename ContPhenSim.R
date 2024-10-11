##############
library(data.table)
library(MASS)
library(broom)
library(dplyr)


############# Functions ################

# Simulating Dataset Function

sim.dataset <- function(n, B, mu, sigma, p, f) {
  simulated_data <- data.table(
    G = rbinom(n, 2, p),
    E = rbinom(n, 1, f)
  )
  
  simulated_data[, GE := G * E]
  XX <- as.matrix(simulated_data)
  YY <- XX %*% B + mvrnorm(n, mu = mu, Sigma = sigma)
  colnames(YY) <- c("Y1", "Y2")
  
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

# Define Parameters
r <- 0.3
mu1 <- 0; s1 <- sqrt(0.5)
mu2 <- 0; s2 <- sqrt(0.5)

# Parameters for Bivariate Normal Distribution
mu <- c(mu1, mu2) # Mean
S <- matrix(c(s1^2, s1 * s2 * r, s1 * s2 * r, s2^2), ncol = 2) # Covariance Matrix

# Simulation Parameters
n.sims <- 500
n <- 1000
p <- 0.2
f <- 0.1

# Adjust B to have the correct dimensions
B <- matrix(c(0.15, 0.2, 0.4, 0.15, 0.2, 0.4), ncol = 2) # B should be 3 x 2

# Run the Simulation
simulated_data <- sim.dataset(n, B, mu, S, p, f)

# Calculate Power
power.sim <- power.fn(n.sims, n, B, mu, S, p, f, alpha)
print(power.sim)
