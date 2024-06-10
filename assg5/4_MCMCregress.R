# 1122 Statistical Computation and Simulation
# Assignment 5

# 4. Use “MCMCregress” in the module MCMCpack
# to obtain MCMC estimation of regression analysis. 
# Duplicate the analysis in the lecture notes
# and apply the MCMC on the “bikes.csv” data. 
# Compare your results with the regular simple linear regression.

library(MCMCpack)
library(readxl)

bikes <- read_excel("assg5_motor_bikes.xlsx", sheet = "Bikes")

plot(bikes$temp_feel, bikes$riders_total)

md1 <- lm(riders_total ~ temp_feel, data = bikes)
summary(md1)

md2 <- MCMCregress(riders_total ~ temp_feel, 
            b0 = 0, b1 = 0.1, sigma.mu = 5, sigma.var = 25, 
            verbose = 1000, data = bikes)

plot(md2)
raftery.diag(md2)
summary(md2)
