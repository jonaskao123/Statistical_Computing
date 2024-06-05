library(plyr)
library(dplyr)
library(readxl)
library(tidyverse)
library(ggplot2)
library(KernSmooth)
library(broom)
library(magrittr)
library(igraph)

SEED <- 123

########## 1. ######

set.seed(SEED)
random_unif = runif(100)
random_norm1 = rnorm(100, -2, 1)
random_norm2 = rnorm(100, 2 ,1)

mixed_data = (random_unif < 0.5) * random_norm1 + (random_unif > 0.5) * random_norm2

###. histogram density estimator
histest = function(x, h) {
  w = function(x, a, b) {
    if (x <= b & x >= a) {return(1)}
    else {return(0)}}
  n = length(x)
  sx = seq(min(x), max(x), by = h)
  a = sx[-length(sx)]
  b = sx[-1]
  ni = NULL
  for (j in 1:length(a)) {
    ni[j] = sum(x <= b[j] & x >= a[j])}
  t1 = NULL
  for (i in sort(x)){
    t0 = NULL
    for (j in 1:length(a)){
      wei = w(i, a[j], b[j])
      t0 = c(t0,wei)}
    y = 1/n * sum(ni/h * t0)
    t1 = c(t1,y)}
  return(t1)}

###. naive density estimator 
naiveest = function(x, h) {
  w = function(y) {
    if (abs(y) < 1) {return(1/2)}
    else {return(0)}
  }
  
  n = length(x)
  sx = seq(min(x), max(x), length = 500)
  t1 = NULL
  
  for (i in sx) {
    t0 = NULL
    for (j in x) {
      wei = w((i - j)/h)
      t0 = c(t0, wei)
    }
    y = 1/n * sum(1/h * t0)
    t1 = c(t1, y)
  }
  return(t1)
}

###. kernel density estimator
kernelest = function(x, h) {
  w <- function(y) {dnorm(y)}
  n <- length(x)
  sx <- seq(min(x), max(x), length = 500)
  
  t1 <- NULL
  for (i in sx) {
    t0 <- NULL
    for (j in x) {
      wei <- w((i - j)/h)
      t0 <- c(t0, wei)
    }
    y <- 1/n * sum(1/h * t0)
    t1 <- c(t1, y)
  }
  return(t1)
}

y1 <- histest(mixed_data, 0.8)
xa <- seq(min(mixed_data), max(mixed_data), length = 500)
y2 <- naiveest(mixed_data, 0.8) 
y3 <- kernelest(mixed_data,0.8) 

par(mfrow = c(1,1))
plot(sort(mixed_data), y1, typ = "l", xlab = "x", ylab = "f(x)", xlim = c(-5, 5),
     ylim = c(0, 0.3), lwd = 2)
matplot(cbind(xa,xa),cbind(y2,y3),typ=c("l", "l"),
        col = c(2, 4),lty = c(2:3), lwd = 2, main = " Three Density Estimates ", add = TRUE)
legend(1,0.3,c("histogram(h = 0.8)","naive(h = 0.8)","kernel(h = 0.8)")
       ,col=c(1,2,4),lty = c(1:3),cex = 0.8)

########## 3. ######

###. kernel smooth
set.seed(SEED)
a <- seq(0, 2*pi, length = 100)
b <- sin(a) + rnorm(100, 0, 0.09)
c <- sin(a)
data <- ksmooth(a, b, kernel = "normal", bandwidth = 0.1) %>% as.data.frame()
data <- cbind(b, c, data) %>% as.data.frame()
ggplot(data,aes(x = x)) + labs(title = "kernel smooth method",x = "x",y = "sin(x)")+
  geom_point(aes(y = b),col = "#333344")+
  geom_vline(xintercept = 0, size=1)+
  geom_hline(yintercept = 0, size=1)+
  geom_line(aes(y = y),col = "#808080", lwd=1)+
  scale_x_continuous(breaks = c(0:2*pi))+
  theme(panel.grid.major = element_line(NA), panel.grid.minor
        = element_line(NA))+
  theme(panel.background = element_rect(color = "black",size = 2))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  theme(legend.title=element_text(size = 24))+
  theme(legend.text=element_text(size = 20))

MSE <- c()
for (i in 1:1000) {
  b <- NULL
  b <- sin(seq(0, 2*pi, length = 100)) + rnorm(100, 0, 0.09)
  b1 <- ksmooth(a, b, kernel = "normal", bandwidth = 0.1)$y
  MSE[i] <- mean((b1 - c)^2)
}
MSE <- mean(MSE)
cat("MSE of kernel smooth:", MSE)

###. LOWESS
set.seed(SEED)
a <- seq(0, 2*pi, length = 100)
b <- sin(a) + rnorm(100, 0, 0.09)
c <- sin(a)
# 計算 f 值為 0.23
lowess(x = a, y = b, f = 0.23)

data <- lowess(x = a, y = b, f = 0.23) %>% as.data.frame()
data <- cbind(b,c,data) %>% as.data.frame()
ggplot(data,aes(x = x))+ labs(title = "LOWESS smooth method", x = "x", y = "sin(x)") +
  geom_point(aes(y = b)) +
  geom_vline(xintercept = 0,size = 1) +
  geom_hline(yintercept = 0,size = 1) +
  geom_line(aes(y = y),col = "#808080",lwd = 1) +
  scale_x_continuous(breaks = c(0:2*pi)) +
  theme(panel.background = element_rect(colour = "black",size = 2)) +
  theme(panel.grid.major = element_line(NA),panel.grid.minor
        =element_line(NA))+
  theme(plot.title = element_text(size = 30, face = "bold")) +
  theme(legend.title = element_text(size = 24))+
  theme(legend.text = element_text(size = 20))

MSE <- c()
for (i in 1:1000) {
  b <- NULL
  b <- sin(seq(0, 2*pi, length = 100)) + rnorm(100, 0, 0.09)
  b1 <- lowess(x = a, y = b, f = 0.23)$y
  MSE[i] <- mean((b1 - c)^2)
}
MSE <- mean(MSE)
cat("MSE of LOWESS:", MSE)

###. Running means
set.seed(SEED)
a <- seq(0, 2*pi, length = 100)
b <- sin(seq(0, 2*pi, length = 100)) + rnorm(100, 0, 0.09)
c <- sin(seq(0, 2*pi, length = 100))
mse = c()
for (k in 1:20){
  r <- running_mean(b, binwidth = k)
  x = NULL
  for(i in 1:(100 - k + 1)){
    x[i] <- mean(a[i:(i + k - 1)])
  }
  mse[k] <- mean((sin(x) - r)^2)
  num = which(mse == min(mse))
}

b <- running_mean(b, binwidth = which(mse == min(mse)))
model <- lm(b ~ poly(seq(0,2*pi, length = 100 - num + 1), 3))
y <- fitted.values(model)
data <- cbind(a[1:length(b)],b,c[1:length(b)]) %>% as.data.frame()
colnames(data) <- c("x","running mean","sin(x)")
data2 <- gather(data,key = "type",value = "value",2:3)
data2
data2$type %<>% as.factor()
ggplot(data2) + labs(title = "Running mean")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  xlim(0,2*pi) + ylim(-1,1) +
  geom_vline(xintercept = 0,size = 1)+
  geom_hline(yintercept = 0,size = 1)+
  geom_line(mapping = aes(x = x,y = value, color = type),lwd = 1.5)+
  geom_point(mapping = aes(x = x, y = value), color = "black",size = 1)+
  theme(legend.text = element_text(size = 16))+
  theme(legend.position = c(0.8,0.8))+
  theme(legend.background = element_rect(size = 0.5, linetype = "solid",fill
                                         = "#FFFFF0",colour = "black"))+
  theme(panel.background = element_rect(color = '#000000',size = 2))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  theme(legend.title = element_text(size = 24))+
  theme(legend.text = element_text(size = 20))
# 1,000 simulation runs (k=2 )
a <- seq(0,2*pi, length = 100)
b <- sin(seq(0,2*pi, length = 100)) + rnorm(100, 0, 0.09)
c <- sin(seq(0,2*pi, length = 100))

# MSE
c <- sin((a[-1] + a[-100])/2)
for (i in 1:1000) {
  b <- NULL
  b <- sin(seq(0,2*pi, length=100)) + rnorm(100,0,0.09)
  b1 <- running_mean(b, binwidth=2)
  MSE[i] <- mean((b1-c)^2)
}
MSE <- mean(MSE)
cat("MSE of Running mean:", MSE)

########## 5. ######
life_table_df <- read_excel('2020生命表.xlsx', col_names = FALSE)
mortality_df <- read_excel('2020死亡率.xlsx')



prior_means <- life_table_df[[2]]
death_counts <- life_table_df[[3]]
total_population <- life_table_df[[4]]
data_means <- as.numeric(mortality_df[1, ]) / 1000

NormalNormal <- function(prior_variance, data_variance) {
  posterior_variance = 1 / (1/prior_variance + 1/data_variance)
  posterior_means = posterior_variance * (prior_means/prior_variance + data_means/data_variance)
  
  posterior_df <- data.frame(
    Age = 0:(length(posterior_means) - 1),
    Prior_Mean = prior_means,
    Data_Mean = data_means,
    Posterior_Mean = posterior_means,
    Posterior_Variance = posterior_variance
  )
  
  error <- posterior_df$Posterior_Mean - posterior_df$Prior_Mean
  return(sum(abs(error)))
}

BetaBinomial <- function(alpha_prior, beta_prior) {
  alpha_posterior = alpha_prior + death_counts
  beta_posterior = beta_prior + total_population - death_counts
  posterior_means = alpha_posterior / (alpha_posterior + beta_posterior)
  
  result_df = data.frame(
    Age = 0:(length(posterior_means) - 1),
    Official_Mortality_Rate = prior_means,
    Death_counts = round(death_counts, 0),
    Total_population = round(total_population, 0),
    Posterior_Mean = posterior_means)
  
  error <- result_df$Posterior_Mean - result_df$Official_Mortality_Rate
  return(sum(abs(error)))
}

cat("Suppose prior_variance = 0.01 and data_variance = 0.01, error:",
    NormalNormal(0.01, 0.01))
cat("Suppose prior_variance = 0.02 and data_variance = 0.01, error:",
    NormalNormal(0.02, 0.01))
cat("Suppose prior_variance = 0.01 and data_variance = 0.02, error:",
    NormalNormal(0.01, 0.02))

cat("Suppose alpha_prior = 2 and beta_prior = 20, error:",
    BetaBinomial(2, 20))
cat("Suppose alpha_prior = 2 and beta_prior = 10, error:",
    BetaBinomial(2, 10))
cat("Suppose alpha_prior = 2 and beta_prior = 5, error:",
    BetaBinomial(2, 5))
