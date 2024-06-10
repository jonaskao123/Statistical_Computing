# 1122 Statistical Computation and Simulation
# Assignment 5

# 6. You can use “MCMClogit” in the module MCMCpack
# to obtain MCMC estimation of logistic regression analysis. 
# Conduct the logistic regression via the “glm” and MCMC, 
# using the data “birthwt”, and comment on the results you found.
# (Note: data("birthwt", package = "MASS"))

library(MCMCpack)
library(magrittr)
data("birthwt", package = "MASS")

corrplot::corrplot(cor(birthwt))


md1 <- glm(low ~ ., family = binomial(link = "logit"), data = birthwt)
summary(md1)
plot(md1)

glm(low ~ smoke + ptl + ht, family = binomial(link = "logit"), data = birthwt) %>% 
  summary


logpriorfun <- function(beta, location, scale){
  sum(dcauchy(beta, location, scale, log=TRUE))
}

## default improper uniform prior
posterior <- MCMClogit(low~age+as.factor(race)+smoke, data=birthwt)
plot(posterior)
summary(posterior)

## multivariate normal prior
posterior <- MCMClogit(low~age+as.factor(race)+smoke, b0=0, B0=.001,
                       data=birthwt)
plot(posterior)
summary(posterior)

posterior <- MCMClogit(low ~ smoke + ptl + ht, b0=0, B0=.001,
                       data=birthwt)
summary(posterior)


## user-defined independent Cauchy prior
logpriorfun <- function(beta){
  sum(dcauchy(beta, log=TRUE))
}

posterior <- MCMClogit(low~age+as.factor(race)+smoke,
                       data=birthwt, user.prior.density=logpriorfun,
                       logfun=TRUE)
plot(posterior)
summary(posterior)


## user-defined independent Cauchy prior with additional args
logpriorfun <- function(beta, location, scale){
  sum(dcauchy(beta, location, scale, log=TRUE))
}

posterior <- MCMClogit(low~age+as.factor(race)+smoke,
                       data=birthwt, user.prior.density=logpriorfun,
                       logfun=TRUE, location=0, scale=10)
plot(posterior)
summary(posterior)
