# Statistical Computing and Simulation ####
# Assignment 2 ####

library(magrittr)
library(markdown)


## 2. Box-Muller, Polar, Ratio-of-uniform ####

### Box-Muller ####
BoxMuller <- function(n){
  u1 <- runif(n)
  u2 <- runif(n)
  theta <- 2 * pi * u1
  e <- -log(u2)
  r <- sqrt(2 * e)
  
  x <- r * cos(theta)
  y <- r * sin(theta)
  res <- cbind(x, y)
  return(res)
}

# set.seed(1)
box_sample <- BoxMuller(10000)

par(mfrow = c(1,2))
hist(box_sample[,"x"], main = "Box-Muller (10,000 runs)", xlab = "x")
hist(box_sample[,"y"], main = "Box-Muller (10,000 runs)", xlab = "y")

par(mfrow = c(1,1))
plot(box_sample[,"x"], box_sample[,"y"],
     main = "Scatter Plot of Box-Muller", xlab = "x", ylab = "y")

ks.test(box_sample[,"x"], "pnorm")
ks.test(box_sample[,"y"], "pnorm")

BoxMuller1 <- function(n){
  res <- BoxMuller(n)
  return(res[,1])
}


KSTestRep <- function(sim, n = 10000, freq = 1000, dist = "pnorm", alpha = 0.05){
  x <- sapply(1:freq, function(n){sim(n)})
  is_norm <- sapply(x, function(x){ks.test(x, dist)[["p.value"]] >= alpha})
  res <- (table(is_norm, dnn = "reject freq"))
  
  return(res)
}

KSTestRep(BoxMuller1, n = 10000, freq = 1000)

### Polar ####

Polar <- function(n){
  v1 <- runif(2*n, min = -1, max = 1)
  v2 <- runif(2*n, min = -1, max = 1)
  w <- v1^2 + v2^2
  v1 <- v1[w < 1][1:n]
  v2 <- v2[w < 1][1:n]
  w <- w[w < 1][1:n]
  
  c <- sqrt((-2)*log(w) / w)
  x <- c*v1
  y <- c*v2
  return(cbind(x, y))
}

Polar1 <- function(n){
  x <- Polar(n)[, "x"]
  return(x)
}

polar_sample <- Polar(10000)

par(mfrow = c(1,2))
hist(polar_sample[,"x"], main = "Polar Method (10,000 runs)", xlab = "x")
hist(polar_sample[,"y"], main = "Polar Method (10,000 runs)", xlab = "y")

par(mfrow = c(1,1))
plot(polar_sample[,"x"], polar_sample[,"y"],
     main = "Scatter Plot of Polar Method", xlab = "x", ylab = "y")

ks.test(polar_sample[,"x"], "pnorm")

KSTestRep(Polar1, n = 10000, freq = 1000)

rbind(BoxMuller = KSTestRep(BoxMuller1, n = 10000, freq = 1000),
      Polar = KSTestRep(Polar1, n = 10000, freq = 1000))


### Ratio-of-uniforms ####

RatioUnif <- function(n){
  u1 <- runif(10*n)
  u2 <- runif(10*n)
  v <- sqrt(2*exp(1))*(2*u2 - 1)
  x <- v / u1
  z <- (x^2) / 4
  x <- x[(z <= ((0.259 / u1) + 0.35)) & (z <= ((-1)*log(u1)))][1:n]
  return(x)
  
}

ratiounif_sample <- RatioUnif(10000)

par(mfrow = c(1,1))
hist(ratiounif_sample, main = "Ratio of Uniform (10,000 runs)", xlab = "x")

ks.test(ratiounif_sample, "pnorm")

KSTestRep(RatioUnif, n = 10000, freq = 1000)

reject_tb <- rbind(BoxMuller = KSTestRep(BoxMuller1),
                   Polar = KSTestRep(Polar1),
                   RatioUnif = KSTestRep(RatioUnif))

colnames(reject_tb) <- c("Reject", "NonReject")

print(reject_tb)

### Congruential Generators ####


CongGen <- function(n, a, c, m){
  u <- runif(1, max = m)
  res <- c(NULL)
  for(i in 1:n){
    u <- (a*u + c) %% m
    res[i] <- u
  }
  return(res/m)
}

BoxMullerCong <- function(n, a, c, m){
  u1 <- CongGen(n, a, c, m)
  u2 <- CongGen(n, a, c, m)
  theta <- 2 * pi * u1
  e <- -log(u2)
  r <- sqrt(2 * e)
  
  x <- r * cos(theta)
  y <- r * sin(theta)
  res <- cbind(x, y)
  return(res)
}

boxcong_sample <- BoxMullerCong(n = 10000, a = 131, c = 0, m = 2^35)
print(range(boxcong_sample[, "x"]))

## 4. Rejection Algorithm (tv-variate) ####

tvVariate <- function(n){
  tv <- c()
  while(length(tv) < n){
    u <- runif(1)
    x <- (u < 0.5)*(1/(4*u - 1)) + (u >= 0.5)*(4*u - 3)
    v <- (u < 0.5)*(u/(x^2)) + (u >= 0.5)*(u) 
    if ((v < (1 - abs(x)/2)) | (v > (1 + (x^2)/v)^(-(v+1)/2))){
      tv <- c(tv, x)
    }
  }
  return(tv)
}

tvariate_sample <- tvVariate(10000)
hist(tvariate_sample, main = "t-Variate Method (10,000 runs)")

range(tvariate_sample)
dtv <- density(tvariate_sample)

hist(tvariate_sample, freq = FALSE, xlim = c(-3, 3), ylim = c(0, 0.4), 
     xlab = "x", main = "Distribution of t-Variate Method (10,000 runs)")
curve(dnorm(x), add = TRUE)

ks.test(tvariate_sample, "pnorm")

reject_tb4 <- KSTestRep(tvVariate)
dimnames(reject_tb4) <- list(kVariateRejectTest = c("Reject", "NonReject"))
print(reject_tb4)



## 6. AR model ####

data(lynx)

### AR(1) ####
ar1_md <- arima(lynx, order = c(1, 0, 0))
ar1_md

### AR(2) ####
ar2_md <- arima(lynx, order = c(2, 0, 0))
ar2_md

### My Program　####

MyAR <- function(x, order){
  data <- x
  data_colname <- paste0("x", c("",1:order))
  data <- cbind(x = data, x_1 = lag(x, -1))
  order <- order - 1
  
  while (order > 0) {
    data <- cbind(data, lag(data[,ncol(data)], -1))
    order <- order - 1
  }
  colnames(data) <- data_colname
  data <- na.omit(data)
  
  md <- lm(x ~ ., data = data)
  
  return(md)
}
MyAR(lynx, 1)


x <- lynx
data <- x
data <- cbind(data, lag(x))
data
ncol(data)
data %>% class

data <- array(x, dim = c(length(x), 1))
cbind(x = data, lag(data[,ncol(data)]))
