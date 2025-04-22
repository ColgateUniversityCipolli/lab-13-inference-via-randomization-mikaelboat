############################################
# LAB 13
############################################
# load libraries
############################################
library(tidyverse)
library(e1071)
############################################
# Q1
finches.dat <- read_csv("zebrafinches.csv")
further.dat <- finches.dat$further

mu.0 <- 0
t.further <- t.test(x = further.dat,
                    alternative = "less",
                    mu = mu.0)
# t-value
t <- t.further$statistic
n <- length(further.dat)

gauss.pdf <- dnorm(x = t)
gauss.cdf <- pnorm(q = t)

prob <- ((skewness(further.dat) / sqrt(n)) * ((2*t^2 + 1)/6)) * gauss.pdf

##############################################
tVals <- seq(-10, 10, 0.1)
t.errors = tibble(t = numeric(length(tVals)), error = numeric(length(tVals)))
i = 1

for(t in tVals){
  t.errors$t[i] <- t
  gauss.pdf <- dnorm(x = t)
  gauss.cdf <- pnorm(q = t)
  curr.error <- ((skewness(further.dat) / sqrt(n)) * ((2*t^2 + 1)/6)) * gauss.pdf
  t.errors$error[i] <- curr.error
  i <- i+1
}

first.p <- ggplot(x = t.errors) +
  geom_line(aes(x = t.errors$t, y = t.errors$error)) +
  theme_bw() +
  labs(y = "Error rate",
       x = "t value")

################################################
# t-value
t <- qnorm(0.05)
n <- length(further.dat)
gauss.pdf <- dnorm(x = t)
gauss.cdf <- pnorm(q = t)
alpha <- 0.05
skew <- skewness(further.dat)

n <- (((skew/(6*(0.1*alpha))) * ((2*t^2) + 1) * gauss.pdf))^2

################################################
# Q2
# re sampling for further data
n <- length(further.dat)
R <- 10000
resamples.further <- tibble(t.stat = numeric(R))
s <- sd(further.dat)

for(i in 1:R){
  
  curr.sample <- sample(x = further.dat,
                        size = n,
                        replace = T)
  t.stat <- mean(curr.sample) / (s/sqrt(n))
  
  resamples.further$t.stat[i] <- t.stat
}
###################################################
# further bootstrapping
delta <- mean(resamples.further$t.stat) - 0

resamples.null.further <- resamples.further |>
  mutate(t.shifted = resamples.further$t.stat - delta)

# re sampling for closer data
closer.dat <- finches.dat$closer
n <- length(closer.dat)
R <- 10000
s < sd(closer.dat)
resamples.closer <- tibble(t.stat = numeric(R))

for (i in 1:R){
  
  curr.sample <- sample(x = closer.dat,
                        size = n,
                        replace = T)
  
  t.stat <- mean(curr.sample) / (s/sqrt(n))
                      
  resamples.closer$t.stat[i] <- curr.tVal$statistic
}
# closer bootstrapping
delta <- mean(resamples.closer$t.stat) - 0

resamples.null.closer <- resamples.closer |>
  mutate(t.shifted = resamples.closer$t.stat - delta)


# re sampling for difference data
difference.dat <- finches.dat$diff
n <- length(closer.dat)
R <- 10000
s <- sd(difference.dat)
resamples.difference <- tibble(t.stat = numeric(R))

for (i in 1:R){
  
  curr.sample <- sample(x = difference.dat,
                        size = n,
                        replace = T)
  
  t.stat <- mean(curr.sample) / (s/sqrt(n))
  
  resamples.difference$t.stat[i] <- curr.tVal$statistic
}
# diff bootstrapping
delta <- mean(resamples.difference$t.stat) - 0

resamples.null.diff <- resamples.difference |>
  mutate(t.shifted = resamples.difference$t.stat - delta)

#####################################
