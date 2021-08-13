### Purpose: simulation 1 for interactions vs stratification
### Author: S Bauldry
### Date: Aug 13, 2021


### Load packages and set working directory
setwd("~/desktop")
library(tidyverse)
library(MASS)
library(ggplot2)
library(ggpubr)

### Simulation 1: function for data generation and parameter estimation
f1_sim1 <- function(nobs = 1000,  # sample size
                    r1   = 0.3,   # correlation between x1 and x2
                    r2   = 0.2,   # correlation between x1 and pr group member
                    r3   = 0.2,   # correlation between x2 and pr group member
                    pr   = 0.4,   # group membership threshold
                    a1   = 2,     # group 1 intercept
                    b11  = 2,     # group 1 effect of x1
                    b12  = 2,     # group 1 effect of x2
                    e1   = 1,     # group 1 error term
                    a0   = 0,     # group 0 intercept
                    b01  = 1,     # group 0 effect of x1
                    b02  = 1,     # group 0 effect of x2
                    e0   = 2) {   # group 0 error term
                    
  # simulate data
  x_mu    <- rbind(0, 0, 0)
  x_sigma <- matrix( c(1, r1, r2, r1, 1, r3, r2, r3, 1), nrow = 3, ncol = 3)
  x       <- mvrnorm(n = nobs, mu = x_mu, Sigma = x_sigma)
  g       <- ifelse(x[,3] >= quantile(x[,3], probs = pr), 1, 0)
  y       <- ifelse(g == 1, 
                    a1 + b11*x[,1] + b12*x[,2] + rnorm(nobs, 0, e1), 
                    a0 + b01*x[,1] + b02*x[,2] + rnorm(nobs, 0, e0))
  x       <- cbind(x, g, y)
  
  colnames(x) <- c("x1", "x2", "x3", "gr", "y")
  d           <- tibble( data.frame(x) )
  
  # fit models with interaction term and stratified
  m1 <- lm(y ~ x1 + x2 + gr + x1:gr, data = d)
  m2 <- lm(y ~ x1 + x2, data = subset(d, gr == 0))
  m3 <- lm(y ~ x1 + x2, data = subset(d, gr == 1))
  
  # extract and return desired estimates
  ib1 <- m1$coefficients[[2]]
  ib2 <- m1$coefficients[[5]]
  sb0 <- m2$coefficients[[2]]
  sb1 <- m3$coefficients[[2]]
  
  c(ib1, ib2, sb0, sb1)
}


### Simulation 1: function for simulation parameters
### Note: need to add function parameters if want to vary parameters in f1
f2_sim1 <- function(nsim = 1000, # number of simulations
                    nobs = 1000, # sample size
                    b11  = 2,    # population group 1 effect of x1
                    b01  = 1) {  # population group 0 effect of x1
  
  # run simulation
  sim <- replicate(n = nsim, f1_sim1(nobs = nobs, b11 = b11, b01 = b01))
  
  # calculate percent bias, RMSEs, and return desired estimates
  ib_0_bias <- 100*abs( mean(sim[1,]) - b01 )/b01
  ib_0_rmse <- sqrt( mean( (sim[1,] - b01)^2 ) )
  
  ib_1_bias <- 100*abs( mean(sim[1,] + sim[2,]) - b11 )/b11
  ib_1_rmse <- sqrt( mean( (sim[1,] + sim[2,] - b11)^2 ) )
  
  sb_0_bias <- 100*abs( mean(sim[3,]) - b01 )/b01
  sb_0_rmse <- sqrt( mean( (sim[3,] - b01)^2 ) )
  
  sb_1_bias <- 100*abs( mean(sim[4,]) - b11 )/b11
  sb_1_rmse <- sqrt( mean( (sim[4,] - b11)^2 ) )
  
  c(nobs, ib_0_bias, ib_0_rmse, ib_1_bias, ib_1_rmse, sb_0_bias, sb_0_rmse, 
    sb_1_bias, sb_1_rmse)
}


### Run simulation 1
sim1_res <- matrix(0, nrow = 10, ncol = 9)
i        <- 0
for(n in seq(100, 1000, 100)) {
  i <- i + 1
  sim1_res[i,] <- f2_sim1(nsim = 1000, nobs = n)
}


### Graph results
colnames(sim1_res) <- c("n", "ib_0_bias", "ib_0_rmse", "ib_1_bias", "ib_1_rmse", 
                        "sb_0_bias", "sb_0_rmse", "sb_1_bias", "sb_1_rmse")
res <- as_tibble(sim1_res)

f1 <- ggplot(data = res, mapping = aes(x = n)) +
  geom_line(mapping = aes(y = ib_0_bias)) +
  geom_line(mapping = aes(y = ib_1_bias), linetype = "dashed") +
  scale_y_continuous(name = "percent bias", limits = c(0, 20)) +
  scale_x_continuous(name = "sample size", breaks = seq(100, 1000, 100)) +
  ggtitle("Percent bias using interaction") +
  theme_light()

f2 <- ggplot(data = res, mapping = aes(x = n)) +
  geom_line(mapping = aes(y = ib_0_rmse)) +
  geom_line(mapping = aes(y = ib_1_rmse), linetype = "dashed") +
  scale_y_continuous(name = "RMSE", limits = c(0, 0.6)) +
  scale_x_continuous(name = "sample size", breaks = seq(100, 1000, 100)) +
  ggtitle("RMSE using interaction") +
  theme_light()

f3 <- ggplot(data = res, mapping = aes(x = n)) +
  geom_line(mapping = aes(y = sb_0_bias)) +
  geom_line(mapping = aes(y = sb_1_bias), linetype = "dashed") +
  scale_y_continuous(name = "percent bias", limits = c(0, 20)) +
  scale_x_continuous(name = "sample size", breaks = seq(100, 1000, 100)) +
  ggtitle("Percent bias using stratification") +
  theme_light()

f4 <- ggplot(data = res, mapping = aes(x = n)) +
  geom_line(mapping = aes(y = sb_0_rmse)) +
  geom_line(mapping = aes(y = sb_1_rmse), linetype = "dashed") +
  scale_y_continuous(name = "RMSE", limits = c(0, 0.6)) +
  scale_x_continuous(name = "sample size", breaks = seq(100, 1000, 100)) +
  ggtitle("RMSE using stratification") +
  theme_light()

fig1 <- ggarrange(f1, f3, f2, f4)
fig1

