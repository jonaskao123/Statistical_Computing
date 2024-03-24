SEED <- 123

###### 1.(a)

library(stats)

theta_values <- c(0, 0.05, 0.1, 0.15, 0.2)

set.seed(SEED)
for (theta in theta_values) {
  correlation1_values <- c()
  correlation2_values <- c()
  
  for (i in 1:1000) {
    arima_sim_data <- arima.sim(model = list(order = c(2,0,0), ar = c(theta, theta)), n = 100, sd = 1)
    correlation1 <- cor(arima_sim_data[1:99], arima_sim_data[2:100])
    correlation2 <- cor(arima_sim_data[1:98], arima_sim_data[3:100])
    correlation1_values <- c(correlation1_values, correlation1)
    correlation2_values <- c(correlation2_values, correlation2)
  }
  
  avg_correlation1 <- mean(correlation1_values)
  avg_correlation2 <- mean(correlation2_values)
  
  print(paste("When theta is", theta, ", the average correlation between x_i and x_i+1 is :", round(avg_correlation1, 5)), quote = FALSE)
  print(paste("When theta is", theta, ", the average correlation between x_i and x_i+1 is :", round(avg_correlation2, 5)), quote = FALSE)
}


###### 1.(b)

library(randtests)


gapTest <- function(data, a, b) {
  gap <- function(data, a, b) {
    n <- length(data)
    x <- c(1:n) * (a < data & data < b)
    x1 <- x[x > 0]
    y <- x1[-1] - x1[-length(x1)] - 1
    return(table(y))
  }
  
  gap_counts <- gap(data, a, b)
  
  vec <- gap_counts
  
  while (vec[length(vec)] < 5) {
    vec <- c(vec[1:(length(vec) - 2)], vec[(length(vec) - 1)] + vec[length(vec)])
  }
  
  result_vector <- numeric()
  for (i in c(0:(length(vec) - 1))) {
    a1 <- (b-a) * (1-b+a)^i
    result_vector <- c(result_vector, a1)
  }
  
  result_vector[length(result_vector)] <- 1 - sum(result_vector[1:(length(result_vector) - 1)])
  
  chisq.test(vec, p = result_vector)
}

numRej1 <- function(n, theta) {
  t1err=0
  for (i in 1:n){
    set.seed(1+i)
    data <- arima.sim(model = list(order = c(2,0,0), ar = c(theta, theta)), n = 100, sd = 1)
    data <- data+abs(min(data))
    data <- data/max(data)
    x = gapTest(data,47/100,97/100)
    if ((x$p.value)<=0.05) (t1err=t1err+1)
  }
  results <- (t1err/n)*100
  cat("number of rejection is",(t1err/n)*100,"%")
}

set.seed(SEED)
for (theta in theta_values) {
  numRej1(1000,theta)
}

permuteTest = function(data, k) {
  permute <- function(data, k) {
    y=rep(10,k)^c((k-1):0)
    x=matrix(data,ncol = k,byrow = T)
    x1=apply(x, 1, rank)
    yy=apply(x1*y, 2, sum)
    return(table(yy))
  }
  
  permute_count <- permute(data, k)
  prob <- rep(1/factorial(k), factorial(k))
  chisq.test(permute_count, p = prob)
}

numRej2 <- function(n,theta) {
  t1err = 0
  for (i in 1:n) {
    set.seed(1 + i)
    data <- arima.sim(model = list(order = c(2,0,0), ar = c(theta, theta)), n = 300, sd = 1)
    if (permuteTest(data, 3)$p.value <= 0.05) {
      t1err = t1err + 1
    }
  }
  results <- (t1err/n)*100
  cat("number of rejection is",(t1err/n)*100,"%")
}

set.seed(SEED)
for (theta in theta_values) {
  numRej2(1000 , theta)
}


numRej3 <- function(n,theta) {
  t1err=0
  for (i in 1:n) {
    data<-arima.sim(model = list(order = c(2,0,0), ar = c(theta, theta)), n = 100, sd = 1)
    x=runs.test(data)
    if ((x$p.value) <= 0.05) (t1err=t1err+1)
  }
  cat("number of rejection is",(t1err/n)*100,"%")
}

set.seed(SEED)
for (theta in theta_values) {
  numRej3(1000, theta)
}

###### 3.

generatePoisson <- function(lambda, starting_point) {
  
  if (starting_point == "mean") {
    i <- lambda
  } else if (starting_point == "0") {
    i <- 0
  } else {
    stop("Invalid starting point. Choose 'mean' or '0'.")
  }
  
  while (TRUE) {
    U <- runif(1)  
    cdf <- ppois(i, lambda)
    if (U >= ppois(i, lambda)) {
      i <- i + 1 
      cdf <- ppois(i, lambda)
    } else {
      return(i) 
    }
  }
}

set.seed(SEED)
results_from_0 <- replicate(100, generatePoisson(10, "0"))

set.seed(SEED)
results_from_mean <- replicate(100, generatePoisson(10, "mean"))

mean_steps_from_0<- mean(results_from_0)
mean_steps_from_mean <- mean(results_from_mean) - 10

cat("Average number of steps starting from 0:", mean_steps_from_0, "\n")
cat("Average number of steps starting from mean:", mean_steps_from_mean, "\n")


###### 5.(a)
A <- matrix(c(1, 0.5, 0.25, 0.125,
            0.5, 1, 0.5, 0.25,  
            0.25, 0.5, 1, 0.5,
            0.125, 0.25, 0.5, 1), nrow = 4, byrow=TRUE)

myChol=function(A) {
  m = nrow(A)
  for (i in 1:(m - 1)) {
    A[i, i]=sqrt(A[i, i] - sum(A[0:(i - 1),i]^2))
    for (j in (i + 1):m) {
      A[i, j] = (A[i, j]-sum(A[0:(i-1),i]*A[0:(i - 1),j]))/A[i, i]
    }
  }
  A[m, m] = sqrt(A[m, m] - sum(A[0:(m - 1), m]^2))
  for (j in 1:m - 1) {
    for (i in (j + 1):m) {
      A[i,j] = 0
    }
  }
  return(A)
}

chol(A)
myChol(A)


###### 5.(b)

eigen_values <- eigen(A)$values
eigen_vectors <- eigen(A)$vectors
eigen_vectors %*% diag(eigen_values) %*% t(eigen_vectors) 

Q <- qr.Q(qr(A))
R <- qr.R(qr(A))
Q %*% R

U <- svd(A)$u
V <- svd(A)$v
D <- svd(A)$d

U %*% diag(D) %*% t(V)

