library(MASS)
library(LaplacesDemon)


N <- 500
omega.true <- c(1.0, -0.4)
lamda.true <- 1.0
beta.true <- c(-1.0, 1.0)
ran.sigma.true <- 0.5


eta.true <- 0.75
delta.true <- 0.98
# eta.true <- 0.9
# delta.true <- 0.9

# error.sigma <- sqrt(6) / pi  # SD = 1
error.sigma <- 1.5 * sqrt(6) / pi  # SD = 1.5
# error.sigma <- 2.0 * sqrt(6) / pi  # SD = 2.0
e.gamma <- 0.5772156

tmb.sigma <- 1.5 * sqrt(6) / pi
tmb.gamma <- 0.5772156

# TMB.true <- rlaplace(N, 0, 1.5 / sqrt(2)) + 1
# TMB.true <- tmb.sigma * (tmb.gamma + log(-log(1 - runif(N, 0, 1)))) + 1
TMB.true <- rnorm(N, 3, 1.0)
# TMB.true <- rbinom(N, 1, 0.5)
# TMB.true <- ifelse(TMB.true == 0, rnorm(N, 0, 0.8), rnorm(N, 2, 1.0))


X <- runif(N, 0, 1)
cov <- cbind(X, TMB.true)

# error <- error.sigma * (e.gamma + log(-log(1 - runif(N, 0, 1))))
# error <- rlaplace(N, 0, 1.0 / sqrt(2))
error <- rnorm(N, 0, 1.0)
# error <- rbinom(N, 1, 0.5)
# error <- ifelse(error == 0, rnorm(N, -1.5, 1.0), rnorm(N, 1.5, 1.0))
TMB.obs <- TMB.true + error

random.u <- rnorm(N, 0, ran.sigma.true)

eta1 <- as.matrix(cov %*% omega.true + random.u)
eta2 <- as.matrix(cov %*% beta.true + random.u)

time.C <- runif(N, 0, 3)
time.T <- rep(0, N)
status <- rep(1 ,N)
for (i in 1:N) {
  time.T[i] <- ((-log(1 - runif(1, 0, 1))) / exp(eta1[i]))^(1 / lamda.true)
  if(time.T[i] > time.C[i]){
    status[i] <- 0
    time.T[i] <- time.C[i]
  }
}

P <- exp(eta2) / (1 + exp(eta2))
RR <- rep(0, N)
RR.true <- rep(0, N)
for(i in 1:N){
  RR.true[i] <- rbinom(1, 1, P[i])
  if(RR.true[i] == 1){
    flag <- sample(x = c(0, 1), size = 1, prob = c(1 - eta.true, eta.true))
    RR[i] <- flag
  }else{
    flag <- sample(x = c(0, 1), size = 1, prob = c(delta.true, 1 - delta.true))
    RR[i] <- flag
  }
}


dat <- cbind(time.T, status, RR.true, RR, X, TMB.true, TMB.obs)
dat <- as.data.frame(dat)
write.table(dat,
            file = "**.csv",
            row.names = F,
            col.names = T,
            sep = ',')