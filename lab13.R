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

n <- (((skew/(6*(0.1*alpha))) * ((2*t^2) + 1) * gauss.pdf))^2

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
s < sd(closer.dat)
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
further.boot.p <- mean(resamples.null.further$t.shifted <= f.delta)
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

closer.boot.ptl <- quantile(resamples.null.closer$t.shifted, 0.05)
closer.t.ptl <- qt(0.05, df = n - 1)

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













