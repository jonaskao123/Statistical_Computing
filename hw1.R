library(magrittr)
library(extraDistr)


midSquare <- function(size = 1){
  initSample <- function(){sample(0:999999, size = 1)}

  num <- initSample()
  nums <- c(num)
  
  for (i in seq(size)){
    num <- num^2 %>% 
      # format(scientific = FALSE, width=12) %>%
      format(scientific = FALSE) %>%
      substr(start = 4, stop = 9) %>%
      # ifelse(. == "      ", initSample(), .) %>%
      ifelse(. == "", initSample(), .) %>%
      as.integer()
    
    num <- ifelse(num == nums[length(nums)], initSample(), num)
    nums <- c(nums, num)
  }
  
  return(nums)
}

midSquare(1000)

num <- 89000
num^2 %>% 
  format(scientific = FALSE) %>%
  # format(scientific = FALSE) %>% 
  substr(start = 4, stop = 9)
  

randSeq <- midSquare(10000)
randSeq <- rnorm(10000, 0, 999999)
randSeq <- rdunif(10000, 0, 999999)
randSeq <- runif(10000, 0, 999999)

ks.test(randSeq, "pdunif", 0, 999999)
ks.test(randSeq, "punif", 0, 999999)


chisq.test(table(randSeq))

############################################################

initSample <- function(){sample(0:999999, size = 1)}


mod <- c(30269, 30307, 30323)
mul <- c(172, 171, 170)

xyz <- sapply(mod, function(x){sample(0:x, size = 1)})

u <- sum(xyz/mod)

xyz <- (xyz * mul) %% mod
xyz


LCG <- function(size, mod, mul, seed = NULL){
  nums <- c()
  if (!is.null(seed)){set.seed(seed)}
  
  xyz <- sapply(mod, function(x){sample(0:x, size = 1)})
  
  for (i in seq(size)){
    u <- sum(xyz / mod)
    nums <- c(nums, u)
    xyz <- (xyz * mul) %% mod
  }
  
  return(nums)
}

randomSeq <- LCG(size = 10000, mod = mod, mul = mul, seed = 354)
randomSeq <- runif(n = 10000, 0, 3)
randomSeq <- rdunif(n = 10000, 0, 3)

ks.test(randomSeq, "punif", 0, 3)
ks.test(randomSeq, "pdunif", 0, 3)


####################################################################

library(nortest)

ad.test(rnorm(10))
ad.test(rt(10, df = 10))
ad.test(rt(10, df = 20))["method"]


n <- c(10, 50, 100)
rnorm(n)

c(
  sapply(n, function(n){ad.test(rnorm(n))[["p.value"]]}),
  sapply(n, function(n){ad.test(rt(n, df = 10))[["p.value"]]}),
  sapply(n, function(n){ad.test(rt(n, df = 20))[["p.value"]]})
)


cbind(
  sapply(n, function(n){ad.test(rnorm(n))}),
  sapply(n, function(n){ad.test(rt(n, df = 10))}),
  sapply(n, function(n){ad.test(rt(n, df = 20))})
) %>% t


normSeq <- rnorm(20)
t10Seq <- rt(10, df = 10)

ad.test(normSeq)
ad.test(t10Seq)

testList <- c(ad.test, cvm.test, lillie.test, pearson.test, sf.test)
testList[[1]](rnorm(10))

class(testList)


res <- rbind(
  sapply(testList, function(test){test(rnorm(10))[["p.value"]]}),
  sapply(testList, function(test){test(rt(10, df = 10))[["p.value"]]}),
  sapply(testList, function(test){test(rt(10, df = 20))[["p.value"]]})
)





colnames(res) <- c(c("ad", "cvm", "lillie", "pearson", "sf"))
rownames(res) <- c("N(0,1)", "t(df=10)", "t(df=20)")
res

resList <- list()
c(resList, n10 = list(res))




norTest <- function(testList, n){
  res <- rbind(
    sapply(testList, function(test){test(rnorm(n))[["p.value"]]}),
    sapply(testList, function(test){test(rt(n, df = 10))[["p.value"]]}),
    sapply(testList, function(test){test(rt(n, df = 20))[["p.value"]]})
  )
  
  colnames(res) <- c(c("ad", "cvm", "lillie", "pearson", "sf"))
  rownames(res) <- c("N(0,1)", "t(df=10)", "t(df=20)")
  
  return(res)
}

res <- norTest(testList, n = 10)

replicate(n = 10000, norTest(testList, n = 10) <= 0.05) %>% apply(MARGIN = c(1,2), sum)


nTrial <- replicate(n=10, sapply(n, function(n){norTest(testList, n)}))

countReject <- function(n_trial, n, alpha){
  replicate(n = n_trial, 
            norTest(testList, n = n) <= alpha) %>% 
    apply(MARGIN = c(1, 2), sum)
}

countReject(n_trial = 1000, n = 10, alpha = 0.05)

lapply(n, countReject, n_trial = 1000, n = n, alpha = 0.05)
rejList <- lapply(n, function(n){countReject(n_trial = 1000, n = n, alpha = 0.05)})
names(rejList) <- paste0("n=", n)
View(rejList)


resList <- lapply(n, function(n){norTest(testList, n)})
names(resList) <- paste0("n=", n)
resList




################################################################
setwd("C:/Users/zuoch/OneDrive/Desktop/Statistical Computing and Simulation")
# pi1M <- read.delim("Pi1MDP.txt", sep="")

# library(data.table)
library(magrittr)
library(extraDistr)
library(randtoolbox)

# pi1M <- fread("Pi1MDP.txt", sep="")
# pi1K <- fread("Pi1KDP.txt", sep="", data.table = FALSE, header = FALSE)
# pi1K <- read.delim("Pi1KDP.txt", sep="", header = FALSE)

pi1K <- readLines("Pi1KDP.txt", warn = FALSE) %>% 
  strsplit(split = "") %>% 
  unlist() %>% 
  as.integer()

pi1M <- readLines("Pi1MDP.txt", warn = FALSE) %>% 
  strsplit(split = "") %>% 
  unlist() %>% 
  as.integer()

table(pi1M) %>% 
  chisq.test()


ks.test(pi1M, "pdunif", 0, 9)

dunifSeq <- rdunif(1000, 0, 9)
ks.test(dunifSeq, "pdunif", 0, 9)

ks.test(pi1K, dunifSeq)

gap.test(dunifSeq/9)
gap.test(runif(100))[["expected"]]

gap.test(pi1M/9)
testGap <- gap.test(pi1M/9, echo = TRUE)
print(testGap)



tb <- table(pi1M)
testChi <- tb %>% chisq.test()
ifelse(testChi[["p.value"]] < 0.05, FALSE, TRUE)

print(testChi)

isRandom <- function(seq, alpha = 0.05){
  tb <- table(pi1M)
  testChi <- tb %>% chisq.test()
  isUnif <- ifelse(testChi[["p.value"]] < alpha, FALSE, TRUE)
  
  r <- diff(range(seq))
  isIndept <- ifelse(gap.test(seq/r)[["p.value"]] < alpha, FALSE, TRUE)
  
  return(isUnif & isIndept)
}


diff(range(dunifSeq))

isRandom(dunifSeq)
isRandom(pi1M)

ks.test()
