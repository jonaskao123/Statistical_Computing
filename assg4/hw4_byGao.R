SEED <- 123

########## 1. ######

###. Hyperbolic
Hyperbolic <- function(q_x) {
  (-q_x^2/((1 - q_x)*log(1 - q_x)) - 0.04)^2}

nlminb(start = 0.04, objective = Hyperbolic)

###. UDD
UniDenth <- function(q_x) {
  (q_x/(1 - (1/2)*q_x) - 0.04)^2}

nlminb(start = 0.04, objective = UniDenth)

###. Constant Force
ConForce <- function(q_x) {
  (-log(1 - q_x) - 0.04)^2}

nlminb(start = 0.04, objective = ConForce)

#.資料來源:https://users.stat.ufl.edu/~rrandles/sta4930/4930lectures/chapter3/chapter3R.pdf
#.在每種假設下，定義一組函數，用來計算m_x和q_x之間的平方差。
#.然後使用nlminb函數找固定m_x下，最小化平方差的q_x。

########## 3. ######

Likelihood <- function(theta) {
  estimator = (1997/(2 + theta)) - (1810/(1 - theta)) + 32/theta
  return(estimator)}

DerivLikeli <- function(theta) {
  estimator = (-1997/(2 + theta)^2) + (1810/(1 - theta)^2) - 32/theta^2
  return(estimator)}

###. Newton
NewtonMethod <- function(estimator, threshold) {
  e <- 1
  iter <- 0
  while(e > threshold) {
    y = estimator - Likelihood(estimator)/DerivLikeli(estimator)
    e = abs(y - estimator)
    estimator = y
    iter = iter + 1}
  return(data.frame(ans = y, Likelihood = Likelihood(y), iteration = iter))}

###. Secant
SecantMethod <- function(a, b, threshold) {
  if(Likelihood(a)*Likelihood(b) > 0) {
    return('No solution by SecantMethod.')}
  
  else{
    iter <- 0
    e <- 1
    while(e > threshold) {
      new = b - Likelihood(b)*(b - a)/(Likelihood(b) - Likelihood(a))
      a = b
      b = new
      e = abs(Likelihood(new))
      iter = iter + 1}}
  return(data.frame(ans = new, Likelihood = Likelihood(new), iteration = iter))}

###. Ridder
RidderMethod = function(a, b, threshold) {
  if(Likelihood(a)*Likelihood(b) > 0) {
    return('No solution by RidderMethod.')}
  
  else{
    iter <- 0
    e <- 1
    while(e > threshold) {
      c = (a + b)/2
      new = c + (c - a)*(sign(Likelihood(a) - Likelihood(b))*Likelihood(c)) / sqrt(Likelihood(c)^2 - Likelihood(a)*Likelihood(b))
      
      if(Likelihood(new) < 0) {
        a = a
        b = new}
      
      else{
        a = new
        b = b}
      
      e = abs(Likelihood(new))
      iter = iter + 1
    }}
  return(data.frame(ans = new, Likelihood = Likelihood(new), iteration = iter))}

NewtonMethod(0.05, 10^-6)
SecantMethod(0.01, 0.09, 10^-6)
RidderMethod(0.01, 0.09, 10^-6)

#. 定義概似函數及其導數，再使用三種方法分別進行估計。
#. 牛頓法:透過概似函數值及導數值迭代，收斂速度快，對初始值選擇敏感。
#. 割線法:透過概似函數值及兩個初始值迭代，計算方便，收斂速度慢。
#. 里德爾法:類似於割線法，但會透過概似函數的符號變化來修正逼近值。
#. 每個方法可以制定初始值和誤差門檻，最後結果為估計值、概似函數值以及迭代次數。

########## 5. ######

###. Numerical Integration
funct = function(x) exp(x^2)
value = round(c(integrate(funct, 0, 1))$value, 6)

###. Monte Carlo Integration
MonteCar = function(n) {
  random = runif(n)
  random_value = funct(random)
  estimate = mean(random_value)
  standerr = sd(random_value)/sqrt(n)
  err = abs(estimate - value)
  r1 = round(c(estimate, standerr, err), 6)
  return(r1)}

###. Antithetic Method
AntiMeth = function(n) {
  random = runif(n)
  random_anti = c(random, 1 - random)
  random_value = funct(random_anti)
  estimate = round(mean(random_value), 6)
  standerr = sd(funct(random) + funct(1 - random))/sqrt(2*n)
  err = abs(estimate - value)
  r2 = round(c(estimate, standerr, err), 6)
  return(r2)}

###. Stratified Sampling Method
StratSamp = function(n, m, N) {
  L = seq(0, 1, length = m + 1)
  Vtheta_S <- matrix(0, N, 2)
  for(i in 1:N) {
    random_value <- funct(runif(n))
    Vtheta_S[i,1] <- mean(random_value)
    theta_MCJ <- c()
    for (j in 1:m) {
      theta_MCJ[j] <- mean(funct(runif(n/m, L[j], L[j+1])))}
    Vtheta_S[i,2] <- mean(theta_MCJ)}
  estimate = round(apply(Vtheta_S, 2, mean)[2], 6)
  standerr = apply(Vtheta_S, 2, sd)[2]
  err = abs(estimate - value)
  r3 = round(c(estimate, standerr, err), 6)
  return(r3)}

Com = function(n) {
  result = cbind(c(value, "", ""), MonteCar(n), 
                 AntiMeth(n/2), StratSamp(n, 4, 1000))
  colnames(result) = c("Numerical ", "Monte Carlo ","Antithetic ", 
                       "Stratified Sampling ")
  rownames(result) = c("Estimate", "SE", "Error")
  return(result)}

set.seed(SEED)
list('n = 50' = Com(50), 'n = 100'= Com(100),
     'n = 500'= Com(500),'n = 1000' = Com(1000))


#. 數值積分:使用integrate函數計算在範圍[0, 1]內的積分。
#. 蒙特卡羅積分:通過服從0到1均勻分配的隨機樣本，帶入函數計算平均值，得到估計積分。
#. 對立數據方法:透過補數，來減少蒙特卡羅積分的估計誤差與標準誤差。
#. 分層抽樣方法:將積分區間分成多個小區間，並在每個子區間內均勻抽樣。
