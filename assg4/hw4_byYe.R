# Statistical Computation and Simulation ####

library(readODS)
library(reshape2)
library(easynls)
library(matlib)
library(pracma)
library(marqLevAlg)

# 2. Gompertz Model ####

## mortgage rate data ####
y3s2 <- read_ods("y3s2-00000.ods", skip = 3)

colnames(y3s2)[1:4] <- c("Y", "year", "gender", "total")
male <- unique(y3s2$gender)[3]

mr <- y3s2[y3s2$gender %in% c(male) & y3s2$year %in% c(2019:2021), -c(1,3,4)]
colnames(mr) <- c("year", "0", "2", "7", 
                  "12", "17" , "22", "27", "32", "37",
                  "42", "47" , "52", "57", "62", "67",
                  "72", "77" , "82", "87", "92", "97", "100")
# rownames(mr) <- c(2019:2021)

mr <- melt(mr, id = "year", variable.name = "age", value.name = "rate")
mr$age <- as.numeric(as.character(mr$age))
mr$rate <- as.numeric(mr$rate)

mr <- mr[mr$age != 0, ]

plot(mr$age, mr$rate, 
     main = "Mortgage Rate in 2019~2021",
     xlab = "age", ylab = "rate")

## Fit nonlinear models {easynls} ####
md1 <- nlsfit(mr[c("age", "rate")], model = 10, 
              start = c(a = 2000, b = 40, c = 0.03))

nlsplot(mr[c("age", "rate")], model = 10, 
        start = c(a = 2000, b = 40, c = 0.03),
        xlab = "age", ylab = "rate")

## Nonlinear Least Squares ####
md2 <- nls(rate ~ a * exp(-b * exp(-c * age)), data = mr, 
           start = list(a = 2000, b = 40, c = 0.03))
summary(md2)

plot(mr$age, mr$rate,
     main = "Mortgage Rate and Gompertz Model",
     xlab = "age", ylab = "rate")
lines(mr$age, predict(md2))

mr.d <- data.frame(x = mr$age, y = mr$rate/1000)
md3 <- nls(y ~ exp(-b * exp(-c * x)), data = mr.d,
           start = list(b = 50, c = 0.05))

summary(md3)
plot(mr.d$x, mr.d$y,
     main = "Mortgage Rate and Gompertz Model",
     xlab = "age", ylab = "rate/1000")
lines(mr.d$x, predict(md3))

md4 <- nls(y ~ b * c^x, data = mr.d,
           start = list(b = 0.1, c = 1))

nls(y ~ b * c^x, data = mr.d,
    start = list(b = 0.1, c = 1))

summary(md4)
plot(mr.d$x, log(mr.d$y),
     main = "Mortgage Rate and Gompertz Model",
     xlab = "age", ylab = "log(rate)")
lines(mr.d$x, log(predict(md4)))

## MyNewton for Gompertz Model ####
mr.d0 <- mr.d[mr.d$x > 0.02, ]

x <- mr.d0$x
y <- mr.d0$y
C <- 1
B <- 0.001

ssr <- function(y, x, B, C){
  return(sum((y - (B * C^x))^2))
}

ssr(y, x, B, C)

bc <- c(B, C)

gd <- function(y, x, B, C){
  dB <- (-2)*sum((y - (B * C^x)) * C^x)
  dC <- (-2)*sum((y - (B * C^x)) * B * x * C^(x-1))
  return(c(dB, dC))
}

g <- gd(y, x, B, C)
g

Hes <- function(y, x, B, C){
  dB2 <- 2 * sum(C^(2*x))
  dC2 <- 2 * sum((B * x * (C^(x-2))) * ((B*((2*x)-1)*(C^x)) - y*(x-1)))
  # dC2 <- 2 * sum((y - B * C^x) * B * x^2 * C^(x-1))
  dBC <- 2 * sum(x * (C^(x-1)) * ((2*B*C^(x)) - y))
  # dBC <- 2 * sum(B * x * C^(2*x-1))
  H <- matrix(c(dB2, dBC, dBC, dC2), nrow = 2)
  return(H)
}

H <- Hes(y, x, B, C)

bc - solve(H) %*% g

bcs <- c(rep = NULL, B = NULL, C = NULL, res = NULL)
bc <- c(B, C)
x <- mr.d0$x
y <- mr.d0$y

for (i in 1:10000){

  g <- gd(y, x, bc[1], bc[2])
  H_1 <- Hes(y, x, bc[1], bc[2])
  
  bc <- bc - solve(H) %*% g
  res <- ssr(y, x, bc[1], bc[2])
  bcs <- rbind(bcs, c(i, bc, res))
}

bcs

## MyLM (Levenberg-Marquardt algorithm) ----

library(pracma)

# tolerance = t, λ = l 
MyLM <- function(f, y, x, bc0, t, l=10, r=10) { 
  
  bc <- bc0
  iter <- 0
  fbc <- function(bc){f(y = y, x = x, bc = bc)}
  print("here")
  while (TRUE) {
    H <- hessian(fbc, bc)
    G <- grad(fbc, bc)
    dk <- inv(H + l * diag(nrow(H))) %*% G   # dk <- solve(H + l * diag(nrow(H)), G)
    bc1 <- bc - dk   # update rule
    print(iter)  # iteration
    # print(l) # λ
    print(bc1) # x1, x2
    print(G)  # ∇f(x)
    print(dk) # d1, d2
    if (Norm(G) < t) break
    l <- ifelse(fbc(bc1) < fbc(bc), l / r, l * r)
    iter <- iter + 1
    bc <- bc1 # update the old point 
  }
  return(bc)
}

f <- function(y, x, bc) {
  return(sum((y - (bc[1] * bc[2]^x))^2))
}

bc0 <- c(0.1, 1)
bc <- MyLM(f, y, x, bc0, t=1e-3, l=0.01, r=1.5)

f(y, x, bc)

## marqLevAlg----
bc0 <- c(0.1, 1)
f <- function(y, x, bc) {
  return(sum((y - (bc[1] * bc[2]^x))^2))
}

estim <- marqLevAlg(bc0, fn = f, x = as.matrix(mr.d$x), y = mr.d$y)

estim$b

pred.y <- log(estim$b[1] * estim$b[2]^mr.d$x)

plot(mr.d$x, log(mr.d$y),
     main = "Mortgage Rate and Gompertz Model with LM method",
     xlab = "age", ylab = "rate/1000")
lines(mr.d$x, pred.y)

## Conjugate Gradient Method----
fn <- function(y = y, x = x, bc) {
  return(sum((y - (bc[1] * bc[2]^x))^2))
}
bc0 <- c(0.1, 1)
## Steepest Descent ####

fn <- function(y, x, bc) {
  return(sum((y - (bc[1] * bc[2]^x))^2))
}
bc0 <- c(0.1, 1)

MySteepestDescent <- function(fn, bc0, x, y, max_iter = 1000, tol = 1e-6, step_size = 0.01) {
  params <- bc0
  
  f <- function(params){fn(y = y, x = x, bc = params)}
  for (iter in 1:max_iter) {
    # Calculate gradient
    # grad <- numerical_grad(fn, y, x, params, n)
    G <- gradient(f, params)
    
    # Update parameters
    params <- params - step_size * G
    
    # Check convergence
    if (max(abs(grad)) < tol) {
      break
    }
  }
  
  params
}

MySteepestDescent(fn, bc0, mr.d$x, mr.d$y)

## Steepest Descent: Steepest Descent ----
fn <- function(bc) {
  y <- mr.d$y
  x <- mr.d$x
  return(sum((y - (bc[1] * bc[2]^x))^2))
}
bc0 <- c(0.0001, 1.1)

fn(bc0)

stpd <- steep_descent(bc0, fn, maxiter = 1000)
stpd

pred.y <- log(stpd$xmin[1] * stpd$xmin[2]^mr.d$x)

plot(mr.d$x, log(mr.d$y),
     main = "Mortgage Rate and Gompertz Model with LM method",
     xlab = "age", ylab = "rate/1000")
lines(mr.d$x, pred.y)

######################################################################

# 4. Evaluate CDF of Normal ----

## Important Sampling ----
x <- c(-6, -5, -4, -3.5, -3, -2.5, -2)
n <- c(100, 1000, 5000)

round(pnorm(x), 9)
pnorm(-2)

n <- 100
x <- -1.5

u <- runif(n)
hist(u)

y <- -(sqrt(3)/pi) * log((1 + exp(-pi*x/sqrt(3))) * u / (1 - u))

hist(y, breaks = 100)

gy <- (pi * exp(-pi*y/sqrt(3))) / (sqrt(3) * (1 + exp(-pi*y/sqrt(3)))^2)

# phiy <- exp(-(y^2)/2) / (sqrt(2) * pi)
phiy <- dnorm(y)

est <- ((1/n) / (1 + exp(-pi*x/sqrt(3)))) * sum(phiy / gy)
est

ImportantSample <- function(x, n){
  u <- runif(n)
  # y <- -(sqrt(3)/pi) * log((1 + exp(-pi*x/sqrt(3))) * u / (1 - u))
  y <- - (sqrt(3)/pi) * log((1 + exp(-pi*x/sqrt(3))) / u - 1)
  gy <- (pi * exp(-pi*y/sqrt(3))) / (sqrt(3) * (1 + exp(-pi*y/sqrt(3)))^2)
  phiy <- dnorm(y)
  est <- ((1/n) / (1 + exp(-pi*x/sqrt(3)))) * sum(phiy / gy)
  return(est)
}

ImportantSample(-2, 100)

x <- c(-6, -5, -4, -3.5, -3, -2.5, -2)
n <- c(100, 1000, 5000)
res_is <- outer(x, n, Vectorize(function(x, n){ImportantSample(x = x, n = n)}))

res_is <- cbind(res_is, pnorm(x))
dimnames(res_is) <- list(paste0("x=", x), c(paste0("n=", n), "pnorm"))

round(res_is, 10)

## Monte Carlo Integration----

x <- 5
n <- 1000
u <- runif(n, 0, 1/x)
1 - pnorm(x)

y <- u^2 / (sqrt(2*pi)) * exp(-1/(2*u^2))

mean(y)

MCInteg <- function(x, n){
  u <- runif(n, 1/x, 0)
  y <- (1 / (sqrt(2*pi))) * exp(-1/(2*u^2)) / (abs(x)*u^2)
  est <- sum(y) / n
  return(est)
}

x <- c(-6, -5, -4, -3.5, -3, -2.5, -2)
n <- c(100, 1000, 5000)
res_mci <- outer(x, n, Vectorize(function(x, n){MCInteg(x = x, n = n)}))

res_mci <- cbind(res_mci, pnorm(x))
dimnames(res_mci) <- list(paste0("x=", x), c(paste0("n=", n), "pnorm"))

round(res_mci, 10)

#################################################################


# 6. Estimate Exponential Distribution ----

## Simple Sample ----
1 - pexp(21.6/15)
1 - pgamma(21.6, 15, 1)
w <- c(1, 2, 3, 4, 5)
a <- 21.6
times <- 1000

x <- rexp(5)
sum(x * w) >= a

n <- 100
ys <- NULL
for (i in 1:n){
  x <- rexp(5)
  y <- as.numeric(sum(x * w) >= a)
  ys <- c(ys, y)
}
sum(ys) / n

ExpSample <- function(n){
  n <- 100
  ys <- NULL
  for (i in 1:n){
    x <- rexp(5)
    y <- as.numeric(sum(x * w) >= a)
    ys <- c(ys, y)
  }
  est <- sum(ys) / n
  return(est)
}

ExpSample(100)
ests <- sapply(1:1000, function(x){ExpSample(100)})

mean(ests)
sd(ests)

hist(ests)

## Sample Gamma distribution----

g <- rgamma(100, shape = 15, rate = 1)

sum(g > a) / 100
ests <- sapply(1:1000, function(x){sum(rgamma(100, 15, 1) > a) / 100})

mean(ests)
sd(ests)
hist(ests)

## New Understand ----

lambda <- 1
num <- 5
critical <- 21.6

mean_s <- sum(seq(1, num))
var_s <- sum((seq(1, num))^2)
sd_s <- sqrt(var_s)

std_c <- (critical - mean_s) / sd_s

# theoretical value
# P(Z > 0.89)
1 - pnorm(std_c)

## 1. ----

p <- 1 - pnorm(std_c)
n <- 100

sum(rnorm(n) > std_c) / n

Rtheta1 <- function(n, std_c){
  sum(rnorm(n) > std_c) / n
}

Rtheta1(n, std_c)

## 2. ----

sum(abs(rnorm(n)) > std_c) / (2 * n)

Rtheta2 <- function(n, std_c){
  sum(abs(rnorm(n)) > std_c) / (2*n)
}

Rtheta2(n, std_c)

## 3. ----

1 - 2*p
n <- 1000
u <- runif(n, min = 0, max = std_c)

0.5*(1 - sum(2 * std_c * (1 / sqrt(2*pi)) * exp(-(u^2) / 2)) / n)

Rtheta3 <- function(n, std_c){
  u <- runif(n, min = 0, max = std_c)
  0.5*(1 - sum(2 * std_c * (1 / sqrt(2*pi)) * exp(-(u^2) / 2)) / n)
}

Rtheta3(n, std_c)

## 4. ----

y <- runif(n, min = 0, max = 1/std_c)

sum((1 / std_c) * (1 / sqrt(2*pi)) * exp(-1 / (2*y^2)) / (y^2)) / n

Rtheta4 <- function(n, std_c){
  y <- runif(n, min = 0, max = 1/std_c)
  sum((1 / std_c) * (1 / sqrt(2*pi)) * exp(-1 / (2*y^2)) / (y^2)) / n
}

Rtheta4(n, std_c)

Rtheta_list <- c(rtheta1=Rtheta1, rtheta2=Rtheta2, rtheta3=Rtheta3, rtheta4=Rtheta4)
names(Rtheta_list)

n <- c(100, 1000, 5000)

EstTheta <- function(Rtheta, n, std_c, time = 1000, var = TRUE){
  ests <- replicate(time, Rtheta(n, std_c))
  if (var == TRUE){
    return(sd(ests))
  }else{
    return(mean(ests))
  }
}

replicate(1000, Rtheta4(n, std_c))

mean_theta <- outer(n, Rtheta_list, 
      Vectorize(function(n, Rtheta){EstTheta(Rtheta, n = n, std_c = std_c, var = FALSE)}))
dimnames(mean_theta) <- list(paste0("N=", n), paste("method", 1:4))
print(mean_theta)

sd_theta <- outer(n, Rtheta_list, 
                    Vectorize(function(n, Rtheta){EstTheta(Rtheta, n = n, std_c = std_c)}))
dimnames(sd_theta) <- list(paste0("N=", n), paste("method", 1:4))
sd_theta <- round(sd_theta, 4)
print(sd_theta)


