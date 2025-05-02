################################################################################
# LAB 13
################################################################################
# load libraries
################################################################################
library(tidyverse)
library(e1071)
library(boot.pval)
library(boot)
################################################################################
# question 1
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

################################################################################
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

################################################################################
# t-value
t <- qnorm(0.05)
n <- length(further.dat)
gauss.pdf <- dnorm(x = t)
gauss.cdf <- pnorm(q = t)
alpha <- 0.05
skew <- skewness(further.dat)

size <- (((skew/(6*(0.1*alpha))) * ((2*t^2) + 1) * gauss.pdf))^2

################################################################################
# Question 2
################################################################################
# re-sampling for further data
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

################################################################################
# further bootstrapping
f.delta <- mean(resamples.further$t.stat) - 0

resamples.null.further <- resamples.further |>
  mutate(t.shifted = resamples.further$t.stat - f.delta)

# re sampling for closer data
closer.dat <- finches.dat$closer
n <- length(closer.dat)
R <- 10000
s <- sd(closer.dat)
resamples.closer <- tibble(t.stat = numeric(R))

for (i in 1:R){
  
  curr.sample <- sample(x = closer.dat,
                        size = n,
                        replace = T)
  
  t.stat <- mean(curr.sample) / (s/sqrt(n))
                      
  resamples.closer$t.stat[i] <- t.stat
}

################################################################################
# closer bootstrapping
c.delta <- mean(resamples.closer$t.stat) - 0

resamples.null.closer <- resamples.closer |>
  mutate(t.shifted = resamples.closer$t.stat - c.delta)


################################################################################
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
  
  resamples.difference$t.stat[i] <- t.stat
}
#################################################################################
# diff bootstrapping
d.delta <- mean(resamples.difference$t.stat) - 0

resamples.null.diff <- resamples.difference |>
  mutate(t.shifted = resamples.difference$t.stat - d.delta)

################################################################################
# 2(b) compute bootstrap p-value
################################################################################
# further
further.boot.p <- mean(resamples.null.further$t.shifted <= f.delta )
further.t.p <- (t.test(x = further.dat, mu = 0, alternative = "less"))$p.value

# closer
closer.boot.p <- mean(resamples.null.closer$t.shifted <= c.delta)
closer.t.p <- (t.test(x = closer.dat, mu = 0, alternative = "greater"))$p.value

# difference
low <- mean(resamples.null.diff$t.shifted <= d.delta)
high <- mean(resamples.null.diff$t.shifted >= d.delta)
diff.boot.p <- low + high
diff.t.p <- (t.test(x = difference.dat, mu = 0, alternative = "two.sided"))$p.value

################################################################################
# 2(c) 5th percentile
################################################################################
further.boot.ptl <- quantile(resamples.null.further$t.shifted, 0.05)
further.t.ptl <- qt(0.05, df = n - 1)

closer.boot.ptl <- quantile(resamples.null.closer$t.shifted, 0.95)
closer.t.ptl <- qt(0.95, df = n - 1)

diff.boot.ptl <- quantile(resamples.null.diff$t.shifted, 0.05)
diff.t.ptl <- qt(0.05, df = n - 1)

################################################################################
# 2(d) Bootstrap CI
################################################################################
# boot function
boot.t.stat <- function(d, i){
  n <- length(d[i])
  s <- sd(d[i])
  mean(d[i]) / (s/sqrt(n))
}

################################################################################
# further bootstrapping
further.boots <- boot(data = further.dat,
              statistic = boot.t.stat,
              R = R)

further.boot.ci <- boot.ci(further.boots, type = "bca")
further.boot.p <- boot.pval(further.boots)


# closer bootstrapping
closer.boots <- boot(data = closer.dat,
                      statistic = boot.t.stat,
                      R = R)

closer.boot.ci <- boot.ci(closer.boots, type = "bca")
closer.boot.p <- boot.pval(closer.boots)

# difference bootstrapping
difference.boots <- boot(data = difference.dat,
                      statistic = boot.t.stat,
                      R = R)

difference.boot.ci <- boot.ci(difference.boots, type = "bca")
difference.boot.p <- boot.pval(difference.boots)


################################################################################
# Question 3
################################################################################
# randomization for further data
shifted.further <- further.dat - mu.0
# perform (randomization) shuffling
further.rand <- tibble(t.stat = numeric(R))
s <- sd(further.dat)
n <- length(further.dat)

for(i in 1:R){
  curr.rand <- shifted.further *
    sample(x = c(-1, 1),
           size = length(shifted.further),
           replace = T)
  
  further.rand$t.stat[i] <- mean(curr.rand) / (s/sqrt(n))
}

further.rand <- further.rand |>
  mutate(t.stat = t.stat + mu.0) 


################################################################################
# randomization for closer data
shifted.closer <- closer.dat - mu.0
# perform (randomization) shuffling
closer.rand <- tibble(t.stat = numeric(R))
s <- sd(closer.dat)
n <- length(closer.dat)

for(i in 1:R){
  curr.rand <- shifted.closer *
    sample(x = c(-1, 1),
           size = length(shifted.closer),
           replace = T)
  
  closer.rand$t.stat[i] <- mean(curr.rand) / (s/sqrt(n))
}
# shift back
closer.rand <- closer.rand |>
  mutate(t.stat = t.stat + mu.0) 

################################################################################
# randomization for difference data
shifted.difference <- difference.dat - mu.0
# perform (randomization) shuffling
diff.rand <- tibble(t.stat = numeric(R))
s <- sd(difference.dat)
n <- length(difference.dat)

for(i in 1:R){
  curr.rand <- shifted.difference *
    sample(x = c(-1, 1),
           size = length(shifted.difference),
           replace = T)
  
  diff.rand$t.stat[i] <- mean(curr.rand) / (s/sqrt(n))
}
# shift back
diff.rand <- diff.rand |>
  mutate(t.stat = t.stat + mu.0) 

################################################################################
# randomization p-value
################################################################################
# further data
delta <- abs(mean(further.dat) - mu.0)
low <- mu.0 - delta # mirror
high<- mu.0 + delta   # xbar

further.rand.p <- mean(further.rand$t.stat <= low) +
  mean(further.rand$t.stat >= high)

################################################################################
# closer data
delta <- abs(mean(closer.dat) - mu.0)
low <- mu.0 - delta # mirror
high<- mu.0 + delta  # xbar

closer.rand.p <- mean(closer.rand$t.stat <= low) +
  mean(closer.rand$t.stat >= high)

################################################################################
# difference data
delta <- abs(mean(difference.dat) - mu.0)
low <- mu.0 - delta # mirror
high<- mu.0 + delta   # xbar

difference.rand.p <- mean(diff.rand$t.stat <= low) +
  mean(diff.rand$t.stat >= high)

################################################################################
# c) confidence interval
################################################################################
# for further data
R <- 1000
mu0.iterate <- 0.001
starting.point <- mean(further.dat)

mu.lower <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- further.dat - mu.lower
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower) # shifting back
  
  # p-value 
  delta <- abs(mean(further.dat) - mu.lower)
  low <- mu.lower - delta # mirror
  high<- mu.lower + delta   # xbar
  p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high)
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower <- mu.lower - mu0.iterate
  }
}

mu.upper <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- further.dat - mu.upper
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  
  # thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper) # shifting back
  
  # p-value 
  delta <- abs(mean(further.dat) - mu.upper)
  low <- mu.upper - delta # mirror
  high<- mu.upper + delta   # xbar
  p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high)
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper <- mu.upper + mu0.iterate
  }
}

further.rand.ci <- c(mu.lower, mu.upper)

################################################################################
# for closer data
R <- 1000
mu0.iterate <- 0.001
starting.point <- mean(closer.dat)

mu.lower <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- closer.dat - mu.lower
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower) # shifting back
  
  # p-value 
  delta <- abs(mean(closer.dat) - mu.lower)
  low <- mu.lower - delta # mirror
  high<- mu.lower + delta   # xbar
  p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high)
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower <- mu.lower - mu0.iterate
  }
}

mu.upper <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- closer.dat - mu.upper
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper) # shifting back
  
  # p-value 
  delta <- abs(mean(closer.dat) - mu.upper)
  low <- mu.upper - delta # mirror
  high<- mu.upper + delta   # xbar
  p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high)
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper <- mu.upper + mu0.iterate
  }
}

closer.rand.ci <- c(mu.lower, mu.upper)

################################################################################
# for difference data
R <- 1000
mu0.iterate <- 0.001
starting.point <- mean(difference.dat)

mu.lower <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- difference.dat - mu.lower
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower) # shifting back
  
  # p-value 
  delta <- abs(mean(difference.dat) - mu.lower)
  low <- mu.lower - delta # mirror
  high<- mu.lower + delta # xbar
  p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high)
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower <- mu.lower - mu0.iterate
  }
}

mu.upper <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- difference.dat - mu.upper
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper) # shifting back
  
  # p-value 
  delta <- abs(mean(difference.dat) - mu.upper)
  low <- mu.upper - delta # mirror
  high<- mu.upper + delta   
  p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high)
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper <- mu.upper + mu0.iterate
  }
}

difference.rand.ci <- c(mu.lower, mu.upper)

################################################################################
# question 4
################################################################################
# using the closer data for the hypothesis test
R <- 10000
n <- 25
alpha <- 0.05
mu <- 0
type1_fixed <- numeric(R)
type1_resample <- numeric(R)

for (i in 1:R) {
  #x <- closer.dat
  x <- rnorm(n, mean = mu)
  s_fixed <- sd(x)
  
  # using original sd
  t_fixed <- mean(x) / (s_fixed / sqrt(n))
  
  # using resample sd
  resample <- sample(x, 
                     size = n, 
                     replace = TRUE)
  
  s_resample <- sd(resample)
  t_resample <- mean(resample) / (s_resample / sqrt(n))

  type1_fixed[i] <- (t_fixed < qt(0.05, n - 1))
  type1_resample[i] <- (t_resample < qt(0.05, n - 1))
}

mean(type1_fixed)    
mean(type1_resample)


################################################################################
# confidence interval coverage
R <- 1000
n <- 25
mu <- 0
coverage_fixed <- numeric(R)
coverage_resample <- numeric(R)

# function for CI using original sd
boot.stat.fixed <- function(d, i){ 
  s_fixed <- sd(x)
  mean(d[i]) / (s_fixed / sqrt(n))
}

# function for CI using resample sd
boot.stat.resample <- function(d, i) {
  x_i <- d[i]
  mean(x_i) / (sd(x_i) / sqrt(n))
}

for (i in 1:R) {
  x <- rnorm(n, mean = mu)
  
  boots_fixed <- boot(data = x, 
                      statistic = boot.stat.fixed, 
                      R = 1000)
  ci_fixed <- boot.ci(boots_fixed, 
                      type = "bca")$bca[4:5]
  coverage_fixed[i] <- (mu >= ci_fixed[1] & mu <= ci_fixed[2])
  
  # boot ci using resample sd
  boots_resample <- boot(data = x, statistic = boot.stat.resample, R = 1000)
  ci_resample <- boot.ci(boots_resample, type = "bca")$bca[4:5]
  coverage_resample[i] <- (mu >= ci_resample[1] & mu <= ci_resample[2])
}

mean(coverage_fixed)    
mean(coverage_resample)  
