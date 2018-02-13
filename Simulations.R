library(pendensity)
library(tidyverse)
library(RColorBrewer)

set.seed(1)

GKDE <- function(x, data.points, bandwidth = 1.06*sd(data.points)*length(data.points)^(-1/5)){
  f.hat <- NULL
  for(i in 1:length(x)){
    f.hat[i] <- mean(dnorm(x[i], mean = data.points, sd = bandwidth))
  }
  return(f.hat)
}

rGKDE <- function(n, data.points, bandwidth = 1.06*sd(data.points)*length(data.points)^(-1/5)){
  rgkde <- rnorm(n, sample(data.points, size = n, replace = TRUE), bandwidth)
  return(rgkde)
}

PDE <- function(x, data.points, pde){
  c.k <- pde$results$ck[1,]
  mu <- pde$splines$MeanW
  sd <- pde$splines$Stand.abw
  f.hat <- NULL
  for(i in 1:length(x)){
    f.hat[i] <- sum(c.k*dnorm(x[i], mean = mu, sd = sd))
  }
  return(f.hat)
}

rPDE <- function(n, data.points, pde){
  K <- pde$spline$K
  c.k <- pde$results$ck[1,]
  mu <- pde$splines$MeanW
  sd <- pde$splines$Stand.abw
  mix.prop <- sample(1:K, prob = c.k, size = n, replace = TRUE)
  rpde <- rnorm(n=n, mean=mu[mix.prop], sd=sd[mix.prop])
  return(rpde)
}

dGausMix <- function(x, prob, mus, vars){
  return(0.3*dnorm(x, mean = mus[1],sd = sqrt(vars[1])) +
           0.5*dnorm(x, mean = mus[2],sd = sqrt(vars[2])) +
           0.2*dnorm(x, mean = mus[3], sd = sqrt(vars[3])))
}

rGausMix <- function(n, prob, mus, vars){
  mix.prop <- sample(1:3, prob = prob, size = n, replace=TRUE)
  rnorm(n=n, mean=mus[mix.prop], sd=sqrt(vars[mix.prop]))
}

SE.GKDE.N <- function(x, data.points, bandwidth = 1.06*sd(data.points)*length(data.points)^(-1/5)){
  f.hat <- NULL
  for(i in 1:length(x)){
    f.hat[i] <- mean(dnorm(x[i], mean = data.points, sd = bandwidth))
  }
  return((f.hat - dnorm(x))^2)
}

SE.PDE.N <- function(x, data.points, pde){
  c.k <- pde$results$ck[1,]
  mu <- pde$splines$MeanW
  sd <- pde$splines$Stand.abw
  f.hat <- NULL
  for(i in 1:length(x)){
    f.hat[i] <- sum(c.k*dnorm(x[i], mean = mu, sd = sd))
  }
  return((f.hat - dnorm(x))^2)
}

SE.GKDE.GMix <- function(x, data.points, bandwidth = 1.06*sd(data.points)*length(data.points)^(-1/5),
                         prob, mus, vars){
  f.hat <- NULL
  for(i in 1:length(x)){
    f.hat[i] <- mean(dnorm(x[i], mean = data.points, sd = bandwidth))
  }
  return((f.hat - dGausMix(x, prob = prob, mus = mus, vars = vars))^2)
}

SE.PDE.GMix <- function(x, data.points, prob, mus, vars, pde){
  c.k <- pde$results$ck[1,]
  mu <- pde$splines$MeanW
  sd <- pde$splines$Stand.abw
  f.hat <- NULL
  for(i in 1:length(x)){
    f.hat[i] <- sum(c.k*dnorm(x[i], mean = mu, sd = sd))
  }
  return((f.hat - dGausMix(x, prob = prob, mus = mus, vars = vars))^2)
}

SE.GKDE.G <- function(x, data.points, bandwidth = 1.06*sd(data.points)*length(data.points)^(-1/5),
                      alpha, beta){
  f.hat <- NULL
  for(i in 1:length(x)){
    f.hat[i] <- mean(dnorm(x[i], mean = data.points, sd = bandwidth))
  }
  return((f.hat - dgamma(x, alpha, beta))^2)
}

SE.PDE.G <- function(x, data.points, alpha, beta, pde){
  c.k <- pde$results$ck[1,]
  mu <- pde$splines$MeanW
  sd <- pde$splines$Stand.abw
  f.hat <- NULL
  for(i in 1:length(x)){
    f.hat[i] <- sum(c.k*dnorm(x[i], mean = mu, sd = sd))
  }
  return((f.hat - dgamma(x, alpha, beta))^2)
}


Tests.f <- function(Sample.Size, N.Sim, prob, mus, vars, alpha, beta){
  Results.PDE <- list(x1 = list(ks = numeric(N.Sim), mise = numeric(N.Sim)),
                      x2 = list(ks = numeric(N.Sim), mise = numeric(N.Sim)),
                      x3 = list(ks = numeric(N.Sim), mise = numeric(N.Sim)))
  Results.GKDE <- list(x1 = list(ks = numeric(N.Sim), mise = numeric(N.Sim)),
                       x2 = list(ks = numeric(N.Sim), mise = numeric(N.Sim)),
                       x3 = list(ks = numeric(N.Sim), mise = numeric(N.Sim)))
  for(i in 1:N.Sim){
    x1.sim <- rnorm(Sample.Size)
    x2.sim <- rGausMix(Sample.Size, prob, mus, vars)
    x3.sim <- rgamma(Sample.Size, alpha, beta)
    pde.x1 <- pendensity(x1.sim~1, base = "gaussian")
    pde.x2 <- pendensity(x2.sim~1, base = "gaussian")
    pde.x3 <- pendensity(x3.sim~1, base = "gaussian")
    Results.PDE$x1$ks[i] <- ks.test(x1.sim,
                                    rPDE(n = Sample.Size,
                                         data.points = x1.sim,
                                         pde = pde.x1))$statistic
    Results.PDE$x1$mise[i] <- integrate(SE.PDE.N,
                                        data.points = x1.sim,
                                        pde = pde.x1,
                                        -Inf, Inf)$value
    Results.PDE$x2$ks[i] <- ks.test(x2.sim,
                                    rPDE(n = Sample.Size,
                                         data.points = x2.sim,
                                         pde = pde.x2))$statistic
    Results.PDE$x2$mise[i] <- integrate(SE.PDE.GMix,
                                        data.points = x2.sim,
                                        prob = prob,
                                        mus = mus,
                                        vars = vars,
                                        pde = pde.x2,
                                        -Inf, Inf)$value
    Results.PDE$x3$ks[i] <- ks.test(x3.sim,
                                    rPDE(n = Sample.Size,
                                         data.points = x3.sim,
                                         pde = pde.x3))$statistic
    Results.PDE$x3$mise[i] <- integrate(SE.PDE.G,
                                        data.points = x3.sim,
                                        alpha = alpha, beta = beta,
                                        pde = pde.x3,
                                        -Inf, Inf)$value
    Results.GKDE$x1$ks[i] <- ks.test(x1.sim,
                                     rGKDE(n = Sample.Size,
                                           data.points = x1.sim))$statistic
    Results.GKDE$x1$mise[i] <- integrate(SE.GKDE.N,
                                         data.points = x1.sim,
                                         -Inf, Inf)$value
    Results.GKDE$x2$ks[i] <- ks.test(x2.sim,
                                     rGKDE(n = Sample.Size,
                                           data.points = x2.sim))$statistic
    Results.GKDE$x2$mise[i] <- integrate(SE.GKDE.GMix,
                                         data.points = x2.sim,
                                         prob = prob, mus = mus, vars = vars,
                                         -Inf, Inf)$value
    Results.GKDE$x3$ks[i] <- ks.test(x3.sim,
                                     rGKDE(n = Sample.Size,
                                           data.points = x3.sim))$statistic
    Results.GKDE$x3$mise[i] <- integrate(SE.GKDE.G,
                                         data.points = x3.sim,
                                         alpha = alpha, beta = beta,
                                         -Inf, Inf)$value
  }
  
  ks <- c(lapply(Results.PDE$x1, mean)$ks,
          lapply(Results.PDE$x2, mean)$ks,
          lapply(Results.PDE$x3, mean)$ks,
          lapply(Results.GKDE$x1, mean)$ks,
          lapply(Results.GKDE$x2, mean)$ks,
          lapply(Results.GKDE$x3, mean)$ks)
  mise <- c(lapply(Results.PDE$x1, mean)$mise,
            lapply(Results.PDE$x2, mean)$mise,
            lapply(Results.PDE$x3, mean)$mise,
            lapply(Results.GKDE$x1, mean)$mise,
            lapply(Results.GKDE$x2, mean)$mise,
            lapply(Results.GKDE$x3, mean)$mise)
  return(list(KS = ks, MISE = mise))
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n.cols = 3
cols = gg_color_hue(n.cols)


# Vector with sample sizes for different tests

N <- c(250, 250, 500, 1000)

# Vector with number of iterations for each test

R <- c(1, rep(1000, 3))

# Generation of data considering one shot experiment with n = 250

n <- N[4]

# First scenario, X ~ f1(x)

x1 <- rnorm(n)
pde1 <- pendensity(x1~1, base = "gaussian")

# Second scenario, X ~ Gaussian Mixture Model with 
# (mu1 = 0, sigma1 = 1), (mu2 = 5, sigma2 = 2) and
# (mu3 = 10, sigma3 = 3). Mixing proportions p1, p2 and p3
# are 0.3, 0.5 and 0.2, respectively

prob = c(0.3,0.5,0.2)
mus = c(0,5,10)
vars = c(1,2,3)

x2 <- rGausMix(n=n, prob = prob, mus = mus, vars = vars)

pde2 <- pendensity(x2~1, base = "gaussian")

# Third scenario, X ~ Gamma(2,10)

a <- 2
b <- 10

x3 <- rgamma(n, a, b)

pde3 <- pendensity(x3~1, base = "gaussian")

# Simulation of n.KS random variables from each density estimator and
# true density for cases X ~ f1, f2 and f3. A dataframe is created to
# store the random variables.

n.KS <- 10000

f1.data <- data.frame(X = c(rPDE(n.KS, x1, pde1),
                            rGKDE(n.KS, x1),
                            rnorm(n.KS)),
                      Density = c(rep("PDE1", n.KS),
                                  rep("GKDE1", n.KS),
                                  rep("f1", n.KS)))

f2.data <- data.frame(X = c(rPDE(n.KS, x2, pde2),
                            rGKDE(n.KS, x2),
                            rGausMix(n.KS,
                                     prob = prob,
                                     mus = mus,
                                     vars = vars)),
                      Density = c(rep("PDE2", n.KS),
                                  rep("GKDE2", n.KS),
                                  rep("f2", n.KS)))

f3.data <- data.frame(X = c(rPDE(n.KS, x3, pde3),
                            rGKDE(n.KS, x3),
                            rgamma(n.KS, a, b)),
                      Density = c(rep("PDE3", n.KS),
                                  rep("GKDE3", n.KS),
                                  rep("f3", n.KS)))

# Experiments

Results.df <- data.frame(f.x = rep(c("PDE1", "PDE2", "PDE3",
                                     "GKDE1", "GKDE2", "GKDE3"), 4),
                         KS.test = numeric(24),
                         MISE = numeric (24),
                         R = c(rep(R[1], 6), rep(R[2], 6), rep(R[3], 6),
                               rep(R[4], 6)),
                         n = c(rep(N[1], 6), rep(N[2], 6), rep(N[3], 6),
                               rep(N[4], 6)))


for(i in 1: length(N)){
  tmp <- Tests.f(Sample.Size = N[i],
                 N.Sim = R[i],
                 prob = prob, mus = mus, vars = vars,
                 alpha = a, beta = b)
  Results.df$KS.test[Results.df$R == R[i] & Results.df$n == N[i]] <- tmp$KS
  Results.df$MISE[Results.df$R == R[i] & Results.df$n == N[i]] <- tmp$MISE
}

Results.df$R <- as.factor(Results.df$R)
Results.df$n <- as.factor(Results.df$n)
Results.df$KS.test <- round(Results.df$KS.test, 4)
Results.df$MISE <- round(Results.df$MISE, 4)

# Plots of the empirical density function (ECDF) for the
# n = 250 random variables generated from the density estimators
# and the true distributions f1, f2 and f3

TextSize <- 20

ggplot(data = f1.data, aes(X, colour = Density)) +
  stat_ecdf(size = 0.5) +
  scale_colour_manual(values = cols,
                      labels = expression(f[1](x), GKDE[1], PDE[1])) +
  labs(x = paste("Cumulative density plots under Scenario I with n =", n)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = TextSize, face = "bold"),
        legend.position = c(0.1, 0.8),
        title = element_text(face = "bold"),
        axis.title.y = element_blank(),
        axis.title = element_text(face = "bold", size = TextSize))

ggsave(paste("ECDF1_n",n,".pdf", sep = ""))

ggplot(data = f2.data, aes(X, colour = Density)) +
  stat_ecdf(size = 0.5) +
  scale_colour_manual(values = cols,
                      labels = expression(f[2](x), GKDE[2], PDE[2])) +
  labs(x = paste("Cumulative density plots under Scenario II with n =", n)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = TextSize, face = "bold"),
        legend.position = c(0.1, 0.8),
        title = element_text(face = "bold"),
        axis.title.y = element_blank(),
        axis.title = element_text(face = "bold", size = TextSize))

ggsave(paste("ECDF2_n",n,".pdf", sep = ""))

ggplot(data = f3.data, aes(X, colour = Density)) +
  stat_ecdf(size = 0.5) +
  scale_colour_manual(values = cols,
                      labels = expression(f[3](x), GKDE[3], PDE[3])) +
  labs(x = paste("Cumulative density plots under Scenario III with n =", n)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = TextSize, face = "bold"),
        legend.position = c(0.1, 0.8),
        title = element_text(face = "bold"),
        axis.title.y = element_blank(),
        axis.title = element_text(face = "bold", size = TextSize))

ggsave(paste("ECDF3_n",n,".pdf", sep = ""))

# Plots of the estimated density functions and the true
# density functions f1, f2 and f3

ggplot(data.frame(x = c(-4,4)), aes(x)) +
  stat_function(fun = function(x) PDE(x, data.points = x1, pde = pde1),
                aes(colour = "PDE"), size = 0.8) +
  stat_function(fun = function(x) GKDE(x, data.points = x1),
                aes(colour = "GKDE"), size = 0.8) + 
  stat_function(fun = function(x) dnorm(x),
                aes(colour = "f1"), size = 0.8) +
  scale_colour_manual(values = cols,
                      labels = expression(f[1](x), GKDE[1], PDE[1])) +
  labs(x = paste("Density plots under Scenario I with n =", n)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = TextSize, face = "bold"),
        legend.position = c(0.1, 0.8),
        title = element_text(face = "bold"),
        axis.title.y = element_blank(),
        axis.title = element_text(face = "bold", size = TextSize))

ggsave(paste("PDF1_n",n,".pdf", sep = ""))

ggplot(data.frame(x = c(-5,15)), aes(x)) +
  stat_function(fun = function(x) PDE(x, data.points = x2, pde = pde2),
                aes(colour = "PDE"), size = 0.8) +
  stat_function(fun = function(x) GKDE(x, data.points = x2),
                aes(colour = "GKDE"), size = 0.8) + 
  stat_function(fun = function(x) dGausMix(x, prob = prob, mus = mus,
                                           vars = vars),
                aes(colour = "f2"), size = 0.8) +
  scale_colour_manual(values = cols,
                      labels = expression(f[2](x), GKDE[2], PDE[2])) +
  labs(x = paste("Density plots under Scenario II with n =", n)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = TextSize, face = "bold"),
        legend.position = c(0.1, 0.8),
        title = element_text(face = "bold"),
        axis.title.y = element_blank(),
        axis.title = element_text(face = "bold", size = TextSize))

ggsave(paste("PDF2_n",n,".pdf", sep = ""))

ggplot(data.frame(x = c(0,1)), aes(x)) +
  stat_function(fun = function(x) PDE(x, data.points = x3, pde = pde3),
                aes(colour = "PDE"), size = 0.8) +
  stat_function(fun = function(x) GKDE(x, data.points = x3),
                aes(colour = "GKDE"), size = 0.8) + 
  stat_function(fun = function(x) dgamma(x, a, b),
                aes(colour = "f2"), size = 0.8) +
  scale_colour_manual(values = cols,
                      labels = expression(f[3](x), GKDE[3], PDE[3])) +
  labs(x = paste("Density plots under Scenario III with n =", n)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = TextSize, face = "bold"),
        legend.position = c(0.8, 0.8),
        title = element_text(face = "bold"),
        axis.title.y = element_blank(),
        axis.title = element_text(face = "bold", size = TextSize))

ggsave(paste("PDF3_n",n,".pdf", sep = ""))