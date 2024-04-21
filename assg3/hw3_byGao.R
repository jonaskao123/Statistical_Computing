library(ggplot2)
library(tidyverse)
library(gridExtra)

###### 2. ##########

## data preprocess. ######
data <- round(read.csv('1970~2005台灣男性死亡率.csv'), 2)
row.names(data) <- data[, 1]
data <- data[, -1]

train <- data[1:31, ]
test <- data[32:36, ]

train_matrix <- as.matrix(log(train))
test_matrix <- as.matrix(log(test))

alphax <- NULL
for (i in 1:17) {
  alphax[i] <- mean(train_matrix[,i])
}

for (i in 1:17) {
  train_matrix[, i] <- train_matrix[, i] - alphax[i]
}

## predictions and calculate MAPE. ######
PreMape <- function(time_component, age_component, singular_value) {
  train_x <- 1970:2000
  lm_model <- lm(time_component*singular_value ~ train_x)
  beta <- as.matrix(unlist(lm_model$coefficients))
  
  pre_y <- beta[1] + beta[2]*matrix(2001:2005, ncol = 1)
  lnmxt <- pre_y %*% matrix(age_component, ncol = 17) +
    matrix(rep(alphax, time=5), ncol=17, byrow = TRUE)
  
  pre_mxt <- exp(lnmxt)
  test_real <- exp(test_matrix)
  mape <- mean(as.vector(abs((pre_mxt - test_real) / test_real)))
  print(mape)
}

## SVD. ######
train_svd <- svd(train_matrix)
singular_value_1 <- train_svd$d[1]
singular_value_1**2 / sum((train_svd$d)**2)

timecompo_svd <- train_svd$u[, 1]
agecompo_svd <- train_svd$v[, 1]

PreMape(timecompo_svd, agecompo_svd, singular_value_1)

## PCA. ######
train_pca <- prcomp(train_matrix) 
summary(train_pca)

timecompo_pca <- train_pca$x[,1]
agecompo_pca <- train_pca$rotation[,1]

PreMape(timecompo_pca, agecompo_pca, train_pca$sdev[1])

###### 4. ##########
MWW <- function(x, y) {
  n_1 <- length(x)
  n_2 <- length(y)
  r_1 <- rank(c(x, y))[1:n_1]
  r_2 <- rank(c(x, y))[-(1:n_1)]
  w_1 <- sum(r_1)
  w_2 <- sum(r_2)
  u_1 <- n_1*n_2 + n_1*(n_1 + 1) / 2 - w_1
  u_2 <- n_1*n_2 + n_2*(n_2 + 1) / 2 - w_2
  u <- min(u_1, u_2)
  e_u <- n_1*n_2 / 2
  var_u <- n_1*n_2*(n_1 + n_2 + 1) / 12
  z = (u - e_u) / sqrt(var_u)
  return(z)
}

CriticalPlot <- function (n_1, n_2) {
  critical <- NULL
  for (i in 1:10000) {
    critical <- c(critical, MWW(rnorm(n_1), rnorm(n_2)))
  }
  
  plot <- critical %>% data.frame(critical) %>%
    ggplot(aes(x = critical))+
    geom_histogram(binwidth = 0.3)+
    theme(axis.title = element_text(size = 10))+
    labs(title = paste("n1 =", n_1, ", n2 =", n_2, "in 10,000 times"), x = 'Critical value', y = 'Frequency')
  return(plot)
}

plot_1 <- CriticalPlot(2, 2)
plot_2 <- CriticalPlot(5, 5)
plot_3 <- CriticalPlot(8, 8)
plot_4 <- CriticalPlot(10, 10)

grid.arrange(plot_1, plot_2, plot_3, plot_4, ncol = 2)

###### 6. ##########

conventional <- c(65, 79, 90, 75, 61, 85, 98, 80, 97, 75)
new <- c(90, 98, 73, 79, 84, 81, 98, 90, 83, 88)

## permutation test. ######
PermuteTest <- function(conventional, new, n_permutations) {
  permuted_diffs <- numeric(n_permutations)
  all_data <- c(conventional, new)
  observed_diff <- mean(new) - mean(conventional)
  
  for (i in 1:n_permutations) {
    permuted_labels <- sample(c(rep("conventional", 10), rep("new", 10)))
    permuted_diffs[i] <- mean(all_data[permuted_labels == "new"]) - mean(all_data[permuted_labels == "conventional"])
  }
  
  p_value <- sum(permuted_diffs >= observed_diff) / n_permutations
  cat("Observed mean difference:", observed_diff, "\n")
  cat("p-value:", p_value, "\n")
}

PermuteTest(conventional, new, 10000)

## non-parametric bootstrap test. ######
t_1 <- NULL
for (i in 1:10000) {
  sample_con <- sample(conventional, size = 10, replace = T)
  sample_new <- sample(new, size = 10, replace = T)
  m_1 <- mean(sample_con) - mean(sample_new)
  t_1 <- c(t_1, m_1)
}

hist(t_1)
c((mean(conventional) - mean(new)) - qt(0.975, 18)*sd(t_1),
  (mean(conventional) - mean(new)) + qt(0.975, 18)*sd(t_1))

## parametric bootstrap test. ######
t_2 <- NULL
for (i in 1:10000) {
  rnorm_con <- rnorm(10, mean = mean(conventional), sd = sd(conventional))
  rnorm_new <- rnorm(10, mean = mean(new), sd = sd(new))
  m_2 <- mean(rnorm_con) - mean(rnorm_new)
  t_2 <- c(t_2, m_2)
}

hist(t_2)
c((mean(rnorm_con) - mean(rnorm_new)) - qt(0.975, 18)*sd(t_2),
  (mean(rnorm_con) - mean(rnorm_new)) + qt(0.975, 18)*sd(t_2))

## parametric test. ######
t.test(conventional, new, paired = T, alternative = "two.sided")
