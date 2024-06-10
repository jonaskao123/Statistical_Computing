# 1122 Statistical Computation and Simulation
# Assignment 5

# 2. Crime data are available in many countries and we can use them the explore
# whether there are hot spots and/or peak seasons. Explore the reported cases of
# stolen motorcycles provided by Taipei City and evaluate which month(s) has the
# highest reported cases of stolen motorcycles (via density estimation methods).
# (Bonus: Explore if there are hot spots in stolen motorcycles.)


# packages
library(readxl)
library(dplyr)
library(lubridate)
require(graphics)
library(smoothr)
library(gplots)
# data

motor <- read_excel("stolen_motorcycles.xlsx",
                    sheet = "motorcycles")

colnames(motor) <- c("pno", "case", "date", "time", "location")
motor[motor$date == "110420", "date"] <- 1110420

motor <- mutate(motor,
                year = as.numeric(substr(date, 1, 3)) + 1911,
                month = as.numeric(substr(date, 4, 5)),
                day = as.numeric(substr(date, 6, 7))) %>% 
  mutate(date = ymd(sprintf("%04d%02d%02d", year, month, day))) %>% 
  mutate(mon_num = interval(min(date), date) %/% months(1) + 1)

st_num <- table(motor$mon_num) %>% 
  data.frame() %>% 
  mutate(Var1 = as.numeric(Var1))
colnames(st_num) <- c("mon", "num")

num_ts <- ts(table(motor$mon_num))

plot(st_num$mon, st_num$num)
ts.plot(num_ts, xlab = "month", ylab = "frequence", 
        main = "Number of Stolen Motorcycles each Month")

ksmooth(st_num$mon, st_num$num)

num_m <- as.matrix(st_num)

num_smooth <- smooth_ksmooth(num_m)
plot(st_num$mon, st_num$num)
lines(num_smooth)

num_smooth_h10 <- smooth_ksmooth(num_m, bandwidth = 10)
num_smooth_h20 <- smooth_ksmooth(num_m, bandwidth = 20)
num_smooth_h60 <- smooth_ksmooth(num_m, bandwidth = 60)
num_smooth_h40 <- smooth_ksmooth(num_m, bandwidth = 40)

plot(st_num$mon, st_num$num,
     main = "Smoothed Data with Kernel Smoother",
     xlab = "month", ylab = "frequence")
lines(num_smooth_h10, col = "red")
lines(num_smooth_h20, col = "blue")
lines(num_smooth_h40, col = "purple")
lines(num_smooth_h60, col = "green")

legend("topright", legend = c("h=10", "h=20", "h=40", "h=60"),
       col = c("red", "blue", "purple", "green"), lwd = 2)

max_mon <- num_smooth_h60[which.max(num_smooth_h60[,2]), 1]
min(motor$date) + months(as.integer(max_mon))

num_spline <- smooth.spline(num_m, df = 1)
plot(num_spline, type = "l")

num_lowess <- lowess(num_m)
plot(num_lowess)

plot(num ~ mon, data = num_m, 
     main = "Smoothed Data with Spline and LOWESS",
     xlab = "month", ylab = "frequence")
lines(num_spline, col = "red", lwd = 2)
lines(num_lowess, col = "blue", lwd = 2)
legend("topright", legend = c("spline", "LOWESS"),
       col = c("red", "blue"), lwd = 2)
