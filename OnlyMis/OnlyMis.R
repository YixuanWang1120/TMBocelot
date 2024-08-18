library(MASS)


dat <- read.csv("***.csv", header = TRUE)
time.T <- dat[, 1]
status <- dat[, 2]
RR <- dat[, 4]
X <- dat[, 5]
TMB.obs <- dat[, 7]
TMB.true <- dat[, 6]
RR.true <- dat[, 3]

N <- nrow(dat)



simu.times <- 1

lamda.set <- rep(0, simu.times)
omega1.set <- rep(0, simu.times)
omega2.set <- rep(0, simu.times)
beta1.set <- rep(0, simu.times)
beta2.set <- rep(0, simu.times)
ran.sigma.set <- rep(0, simu.times)
eta.set <- rep(0, simu.times)
delta.set <- rep(0, simu.times)

SE.lamda.set <- rep(0, simu.times)
SE.omega1.set <- rep(0, simu.times)
SE.omega2.set <- rep(0, simu.times)
SE.beta1.set <- rep(0, simu.times)
SE.beta2.set <- rep(0, simu.times)
SE.ran.sigma.set <- rep(0, simu.times)
SE.eta.set <- rep(0, simu.times)
SE.delta.set <- rep(0, simu.times)

simu.iter <- 1


while(simu.iter <= simu.times){
  cat("simu:", simu.iter, "\n")
  MCMC.iter.max <- 2000
  
  lamda.mcmc <- rep(0, MCMC.iter.max + 1)
  omega1.mcmc <- rep(0, MCMC.iter.max + 1)
  omega2.mcmc <- rep(0, MCMC.iter.max + 1)
  beta1.mcmc <- rep(0, MCMC.iter.max + 1)
  beta2.mcmc <- rep(0, MCMC.iter.max + 1)
  
  random.u.mcmc <- matrix(0, nrow = MCMC.iter.max + 1, ncol = N)
  ran.sigma.mcmc <- rep(0, MCMC.iter.max + 1)
  
  
  eta.mcmc <- rep(0, MCMC.iter.max + 1)
  delta.mcmc <- rep(0, MCMC.iter.max + 1)
  RR.mcmc <- matrix(0, nrow = MCMC.iter.max + 1, ncol = N)
  
  hype.a <- 0.001
  hype.b <- 0.001
  
  hype.sigma <- 5
  
  hype.misclass.eta <- c(55.8, 19.2)
  hype.misclass.delta <- c(72.5, 2.5)
  
  # hype.misclass.eta <- c(72.5, 2.5)
  # hype.misclass.delta <- c(72.5, 2.5)
  
  lamda.mcmc[1] <- init.lamda <- 1.0
  omega1.mcmc[1] <- init.omega1 <- 1.0
  omega2.mcmc[1] <- init.omega2 <- -1.0
  beta1.mcmc[1] <- init.beta1 <- -0.8
  beta2.mcmc[1] <- init.beta2 <- 0.4
  ran.sigma.mcmc[1] <- init.ran.sigma <- 0.5
  eta.mcmc[1] <- init.eta <- 0.75
  delta.mcmc[1] <- init.delta <- 0.98
  init.TMB <- TMB.obs
  
  
  # C.hazard <- matrix(c(0.1, 0, 0, 0, 0.5, 0, 0, 0, 0.5), nrow = 3)
  C.lamda <- 1
  C.omega1 <- 0.1
  C.omega2 <- 0.1
  C.beta1 <- 0.1
  C.beta2 <- 0.1
  C.mis <- matrix(c(0.1, 0, 0, 0.1), nrow = 2)
  
  random.u.mcmc[1, ] <- init.random.u <- rnorm(N, 0, 0.1)
  C.u <- rep(0.5, N)
  
  
  
  i <- 1
  threshold <- 20
  
  repeat{
    
    temp.lamda <- rnorm(1, init.lamda, sqrt(C.lamda))
    line.expression <- X * init.omega1 + init.TMB * init.omega2 + init.random.u
    accept <- min(1, (prod(((temp.lamda * time.T^(temp.lamda - 1))^status * exp(-time.T^temp.lamda * exp(line.expression)))
                           / ((init.lamda * time.T^(init.lamda - 1))^status * exp(-time.T^init.lamda * exp(line.expression))))
                      * (dgamma(temp.lamda, hype.a, hype.b) / dgamma(init.lamda, hype.a, hype.b)))
    )
    
    u <- runif(1, 0, 1)
    if(!is.na(accept) && u < accept){
      init.lamda <- temp.lamda
    }
    lamda.mcmc[i + 1] <- init.lamda
    # cat("init.lamda","\n")
    if(i > threshold){
      mean.old <- mean(lamda.mcmc[1 : i])
      mean.new <- (mean.old * i + init.lamda) / (i + 1)
      
      C.lamda <- (i - 1) / i * C.lamda + 2.4 * 2.4 / i * (
        i * mean.old^2 - (i + 1) * mean.new^2 + init.lamda^2 + .Machine$double.eps)
    }
    
    
    temp.omega1 <- rnorm(1, init.omega1, sqrt(C.omega1))
    accept <- min(1, (prod((exp(temp.omega1 * X * status) * exp(-time.T^init.lamda * exp(X * temp.omega1 + init.TMB * init.omega2 + init.random.u)))
                           / (exp(init.omega1 * X * status) * exp(-time.T^init.lamda * exp(X * init.omega1 + init.TMB * init.omega2 + init.random.u))))
                      * (dnorm(temp.omega1, 0, hype.sigma) / dnorm(init.omega1, 0, hype.sigma)))
    )
    u <- runif(1, 0, 1)
    if(!is.na(accept) && u < accept){
      init.omega1 <- temp.omega1
    }
    omega1.mcmc[i + 1] <- init.omega1
    
    if(i > threshold){
      mean.old <- mean(omega1.mcmc[1 : i])
      mean.new <- (mean.old * i + init.omega1) / (i + 1)
      
      C.omega1 <- (i - 1) / i * C.omega1 + 2.4 * 2.4 / i * (
        i * mean.old^2 - (i + 1) * mean.new^2 + init.omega1^2 + .Machine$double.eps)
    }
    
    temp.omega2 <- rnorm(1, init.omega2, sqrt(C.omega2))
    accept <- min(1, (prod((exp(temp.omega2 * init.TMB * status) * exp(-time.T^init.lamda * exp(X * init.omega1 + init.TMB * temp.omega2 + init.random.u)))
                           / (exp(init.omega2 * init.TMB * status) * exp(-time.T^init.lamda * exp(X * init.omega1 + init.TMB * init.omega2 + init.random.u))))
                      * (dnorm(temp.omega2, 0, hype.sigma) / dnorm(init.omega2, 0, hype.sigma)))
    )
    u <- runif(1, 0, 1)
    if(!is.na(accept) && u < accept){
      init.omega2 <- temp.omega2
    }
    omega2.mcmc[i + 1] <- init.omega2
    
    if(i > threshold){
      mean.old <- mean(omega2.mcmc[1 : i])
      mean.new <- (mean.old * i + init.omega2) / (i + 1)
      
      C.omega2 <- (i - 1) / i * C.omega2 + 2.4 * 2.4 / i * (
        i * mean.old^2 - (i + 1) * mean.new^2 + init.omega2^2 + .Machine$double.eps)
    }
    
    temp.beta1 <- rnorm(1, init.beta1, sqrt(C.beta1))
    
    accept <- min(1,
                  (((prod(((1 - init.delta + (init.eta + init.delta - 1) * exp(X * temp.beta1 + init.TMB * init.beta2 + init.random.u) / (1 + exp(X * temp.beta1 + init.TMB * init.beta2 + init.random.u)))
                           / (1 - init.delta + (init.eta + init.delta - 1) * exp(X * init.beta1 + init.TMB * init.beta2 + init.random.u) / (1 + exp(X * init.beta1 + init.TMB * init.beta2 + init.random.u))))^RR))
                    * (prod(((init.delta + (1 - init.eta - init.delta) * exp(X * temp.beta1 + init.TMB * init.beta2 + init.random.u) / (1 + exp(X * temp.beta1 + init.TMB * init.beta2 + init.random.u)))
                             / (init.delta + (1 - init.eta - init.delta) * exp(X * init.beta1 + init.TMB * init.beta2 + init.random.u) / (1 + exp(X * init.beta1 + init.TMB * init.beta2 + init.random.u))))^(1 - RR))))
                   * (dnorm(temp.beta1, 0, hype.sigma) / dnorm(init.beta1, 0, hype.sigma)))
    )
    
    u <- runif(1, 0, 1)
    
    if(!is.na(accept) && u < accept){
      init.beta1 <- temp.beta1
    }
    beta1.mcmc[i + 1] <- init.beta1
    # cat("beta1: ", init.beta1, "\n")
    
    if(i > threshold){
      mean.old <- mean(beta1.mcmc[1 : i])
      mean.new <- (mean.old * i + init.beta1) / (i + 1)
      
      C.beta1 <- (i - 1) / i * C.beta1 + 2.4 * 2.4 / i * (
        i * mean.old^2 - (i + 1) * mean.new^2 + init.beta1^2 + .Machine$double.eps)
      
    }
    
    
    temp.beta2 <- rnorm(1, init.beta2, sqrt(C.beta2))
    
    accept <- min(1,
                  ((prod(((1 - init.delta + (init.eta + init.delta - 1) * exp(X * init.beta1 + init.TMB * temp.beta2 + init.random.u) / (1 + exp(X * init.beta1 + init.TMB * temp.beta2 + init.random.u)))
                          / (1 - init.delta + (init.eta + init.delta - 1) * exp(X * init.beta1 + init.TMB * init.beta2 + init.random.u) / (1 + exp(X * init.beta1 + init.TMB * init.beta2 + init.random.u))))^RR))
                   * (prod(((init.delta + (1 - init.eta - init.delta) * exp(X * init.beta1 + init.TMB * temp.beta2 + init.random.u) / (1 + exp(X * init.beta1 + init.TMB * temp.beta2 + init.random.u)))
                            / (init.delta + (1 - init.eta - init.delta) * exp(X * init.beta1 + init.TMB * init.beta2 + init.random.u) / (1 + exp(X * init.beta1 + init.TMB * init.beta2 + init.random.u))))^(1 - RR)))
                   * (dnorm(temp.beta2, 0, hype.sigma) / dnorm(init.beta2, 0, hype.sigma)))
    )
    u <- runif(1, 0, 1)
    if(!is.na(accept) && u < accept){
      init.beta2 <- temp.beta2
    }
    beta2.mcmc[i + 1] <- init.beta2
    # cat("beta2: ", init.beta2, "\n")
    
    if(i > threshold){
      mean.old <- mean(beta2.mcmc[1 : i])
      mean.new <- (mean.old * i + init.beta2) / (i + 1)
      
      C.beta2 <- (i - 1) / i * C.beta2 + 2.4 * 2.4 / i * (
        i * mean.old^2 - (i + 1) * mean.new^2 + init.beta2^2 + .Machine$double.eps)
      
    }
    
    # temp.random.u <- mvrnorm(1, init.random.u, diag(C.u))
    temp.random.u <- rnorm(N, init.random.u, sqrt(C.u))
    temp1 <- init.omega1 * X + init.omega2 * init.TMB
    temp2 <- X * init.beta1 + init.TMB * init.beta2
    temp.accept <- ( ((exp(temp.random.u * status) * exp(-time.T^init.lamda * exp(temp1 + temp.random.u))) / (exp(init.random.u * status) * exp(-time.T^init.lamda * exp(temp1 + init.random.u))))
                     * (((1 - init.delta + (init.eta + init.delta - 1) * exp(temp2 + temp.random.u) / (1 + exp(temp2 + temp.random.u)))) / ((1 - init.delta + (init.eta + init.delta - 1) * exp(temp2 + init.random.u) / (1 + exp(temp2 + init.random.u)))))^RR
                     * (((init.delta + (1 - init.eta - init.delta) * exp(temp2 + temp.random.u) / (1 + exp(temp2 + temp.random.u)))) / ((init.delta + (1 - init.eta - init.delta) * exp(temp2 + init.random.u) / (1 + exp(temp2 + init.random.u)))))^(1 - RR)
                     * (dnorm(temp.random.u, 0, init.ran.sigma) / dnorm(init.random.u, 0, init.ran.sigma))
    )
    accept <- ifelse(temp.accept > 1, 1, temp.accept)
    u <- runif(N, 0, 1)
    
    init.random.u <- ifelse(!is.na(accept) & u < accept, temp.random.u, init.random.u)
    random.u.mcmc[i + 1, ] <- init.random.u
    
    if(i > threshold){
      mean.old <- colMeans(random.u.mcmc[1:i, ])
      mean.new <- (mean.old * i + init.random.u) / (i + 1)
      
      C.u <- (i - 1) / i * C.u + 2.4 * 2.4 / i * (
        i * mean.old * mean.old - (i + 1) * mean.new * mean.new + init.random.u * init.random.u + .Machine$double.eps)
    }
    
    
    
    ran.sigma.mcmc[i + 1] <- init.ran.sigma <- 1 / sqrt(rgamma(1, 0.001 + N / 2, 0.001 + sum(init.random.u^2) / 2))
    # cat("sigma: ", init.ran.sigma, "\n")
    
    temp.mis <- mvrnorm(1, c(init.eta, init.delta), C.mis)
    temp.eta <- temp.mis[1]
    temp.delta <- temp.mis[2]
    temp.expression <- X * init.beta1 + init.TMB * init.beta2 + init.random.u
    
    accept <- min(1,  (((prod(((1 - temp.delta + (temp.eta + temp.delta - 1) * exp(temp.expression) / (1 + exp(temp.expression)))
                               / (1 - init.delta + (init.eta + init.delta - 1) * exp(temp.expression) / (1 + exp(temp.expression))))^RR))
                        * (prod(((temp.delta + (1 - temp.eta - temp.delta) * exp(temp.expression) / (1 + exp(temp.expression)))
                                 / (init.delta + (1 - init.eta - init.delta) * exp(temp.expression) / (1 + exp(temp.expression))))^(1 - RR))))
                       * (dbeta(temp.eta, hype.misclass.eta[1], hype.misclass.eta[2]) / dbeta(init.eta, hype.misclass.eta[1], hype.misclass.eta[2]))
                       * (dbeta(temp.delta, hype.misclass.delta[1], hype.misclass.delta[2]) / dbeta(init.delta, hype.misclass.delta[1], hype.misclass.delta[2])))
    )
    
    u <- runif(1, 0, 1)
    if(!is.na(accept) && u < accept){
      if(temp.delta > 0.95 && abs(temp.eta - 0.75) < 0.05){
        init.eta <- temp.eta
        init.delta <- temp.delta
      }
      
    }
    
    eta.mcmc[i + 1] <- init.eta
    delta.mcmc[i + 1] <- init.delta
    
    if(i > threshold){
      mean.old <- c(mean(eta.mcmc[1 : i]), mean(delta.mcmc[1 : i]))
      mean.new <- (mean.old * i + c(init.eta, init.delta)) / (i + 1)
      
      C.mis <- (i - 1) / i * C.mis + (2.4 * 2.4) / (2 * i) * (
        i * mean.old %*% t(mean.old) - (i + 1) * mean.new %*% t(mean.new) + c(init.eta, init.delta) %*% t(c(init.eta, init.delta)) + .Machine$double.eps * diag(2))
      
    }
    
    init.p <- exp(X * init.beta1 + init.TMB * init.beta2) / (1 + exp(X * init.beta1 + init.TMB * init.beta2))
    init.p01 <- (1 - init.p) * (1 - init.delta) / ((1 - init.p) * (1 - init.delta) + init.p * init.eta)
    init.p10 <- init.p * (1 - init.eta) / (init.p * (1 - init.eta) + (1 - init.p) * init.delta)
    
    for(k in 1:N){
      if(RR[k] == 1){
        init.m01 <- rbinom(1, size = 1, init.p01[k])
        init.m11 <- 1 - init.m01
        RR.mcmc[i + 1, k] <- init.m11
      }else{
        init.m10 <- rbinom(1, size = 1, init.p10[k])
        init.m00 <- 1 - init.m10
        RR.mcmc[i + 1, k] <- init.m10
      }
      
    }
    
    
    
    # cat("i: ", i, '\n')
    i <- i + 1
    
    if(i > MCMC.iter.max){
      break
    }
    
  }
  j <- seq(from = 1000, to = MCMC.iter.max, by = 1)
  
  lamda.set[simu.iter] <- mean(lamda.mcmc[j])
  omega1.set[simu.iter] <- mean(omega1.mcmc[j])
  omega2.set[simu.iter] <- mean(omega2.mcmc[j])
  beta1.set[simu.iter] <- mean(beta1.mcmc[j])
  beta2.set[simu.iter] <- mean(beta2.mcmc[j])
  ran.sigma.set[simu.iter] <- mean(ran.sigma.mcmc[j])
  eta.set[simu.iter] <- mean(eta.mcmc[j])
  delta.set[simu.iter] <- mean(delta.mcmc[j])
  
  SE.lamda.set[simu.iter] <- sqrt(var(lamda.mcmc[j]))
  SE.omega1.set[simu.iter] <- sqrt(var(omega1.mcmc[j]))
  SE.omega2.set[simu.iter] <- sqrt(var(omega2.mcmc[j]))
  SE.beta1.set[simu.iter] <- sqrt(var(beta1.mcmc[j]))
  SE.beta2.set[simu.iter] <- sqrt(var(beta2.mcmc[j]))
  SE.ran.sigma.set[simu.iter] <- sqrt(var(ran.sigma.mcmc[j]))
  SE.eta.set[simu.iter] <- sqrt(var(eta.mcmc[j]))
  SE.delta.set[simu.iter] <- sqrt(var(delta.mcmc[j]))
  
  simu.iter <- simu.iter + 1

  
}



plot(x = 1:MCMC.iter.max, y = lamda.mcmc[1:MCMC.iter.max], type = "l")
plot(x = 1:MCMC.iter.max, y = omega1.mcmc[1:MCMC.iter.max], type = "l")
plot(x = 1:MCMC.iter.max, y = omega2.mcmc[1:MCMC.iter.max], type = "l")
plot(x = 1:MCMC.iter.max, y = beta1.mcmc[1:MCMC.iter.max], type = "l")
plot(x = 1:MCMC.iter.max, y = beta2.mcmc[1:MCMC.iter.max], type = "l")
plot(x = 1:MCMC.iter.max, y = ran.sigma.mcmc[1:MCMC.iter.max], type = "l")
plot(x = 1:MCMC.iter.max, y = eta.mcmc[1:MCMC.iter.max], type = "l")
plot(x = 1:MCMC.iter.max, y = delta.mcmc[1:MCMC.iter.max], type = "l")

index <- 1

result.lamda <- mean(lamda.set[index])
result.omega1 <- mean(omega1.set[index])
result.omega2 <- mean(omega2.set[index])
result.beta1 <- mean(beta1.set[index])
result.beta2 <- mean(beta2.set[index])
result.sigma <- mean(ran.sigma.set[index])
result.eta <- mean(eta.set[index])
result.delta <- mean(delta.set[index])
result.u <- colMeans(random.u.mcmc[1000:2000, ])
cat("---------------------\n",
    "lamda: ", mean(lamda.set[index]), "SE: ", mean(SE.lamda.set[index]), "\n",
    
    "omega1: ", mean(omega1.set[index]), "SE: ", mean(SE.omega1.set[index]), "\n", 
    
    "omega2: ", mean(omega2.set[index]), "SE: ", mean(SE.omega2.set[index]),"\n",
    
    "beta1: ", mean(beta1.set[index]), "SE: ", mean(SE.beta1.set[index]),"\n", 
    
    "beta2: ", mean(beta2.set[index]), "SE: ", mean(SE.beta2.set[index]),"\n",
    
    "random effect parameter: ", mean(ran.sigma.set[index]),"SE: ", mean(SE.ran.sigma.set[index]),"\n",
    
    "eta: ", mean(eta.set[index]),  "SE: ", mean(SE.eta.set[index]),"\n",
    
    "detla: ", mean(delta.set[index]), "SE: ", mean(SE.delta.set[index]),"\n")

cat("-----------------------------------------------------------------------------------------\n")

TMB <- TMB.obs
temp.RR <- colMeans(RR.mcmc[1000:2000, ])

P_R <- function() {  
  p <- ((1 + exp(-X * result.beta1 - result.beta2 * TMB - result.u)) ^ (-1))  
  return(p)  
}  

# 定义P_T函数  
T0 <- median(time.T)
P_T <- function() {  
  p <- (exp(-T0 ^ result.lamda * exp(X * result.omega1 + result.omega2 * TMB + result.u)))  
  return(p)  
}  

# 定义p_random_effect函数  
p_random_effect <- function() {  
  p <- ((2 * pi * result.sigma^2) ^ 0.5) ^ (-1) * exp(-(1 / 2) * (result.u^2) / (result.sigma^2))  
  return(p)  
}  

P_J <- function() {
  f <- P_R() * P_T() * p_random_effect()
  return(f)
}
pro <- P_J()

# TMB-cat计算阈值
source("TMB-cat.R")
source("ROC.R")
source("Median_thres.R")
cat("----------------------------------------------TMB-cat----------------------------------------------\n")
result <- TMB.cat.computeThres(pro, TMB, time.T, status, temp.RR)
# cat("----------------------------------------------Median----------------------------------------------\n")
# Median.computeThres(TMB.obs, time.T, status, RR)

thres <- result$thres
group <- rep(0, N)
group <- ifelse(TMB >= thres, 1, 0)
df.sur <- data.frame(
  time = time.T * 12,
  event = status,
  group = group
)
cat(">=0: ", sum(df.sur[ which(df.sur[, 1] >= 0), 3]), length(df.sur[ which(df.sur[, 1] >= 0), 3]) - sum(df.sur[ which(df.sur[, 1] >= 0), 3]))
cat(">=10: ", sum(df.sur[ which(df.sur[, 1] >= 10), 3]), length(df.sur[ which(df.sur[, 1] >= 10), 3]) - sum(df.sur[ which(df.sur[, 1] >= 10), 3]))
cat(">=20: ", sum(df.sur[ which(df.sur[, 1] >= 20), 3]), length(df.sur[ which(df.sur[, 1] >= 20), 3]) - sum(df.sur[ which(df.sur[, 1] >= 20), 3]))
cat(">=30: ", sum(df.sur[ which(df.sur[, 1] >= 30), 3]), length(df.sur[ which(df.sur[, 1] >= 30), 3]) - sum(df.sur[ which(df.sur[, 1] >= 30), 3]))
cat(">=0: ", sum(df.sur[ which(df.sur[, 1] >= 40), 3]), length(df.sur[ which(df.sur[, 1] >= 40), 3]) - sum(df.sur[ which(df.sur[, 1] >= 40), 3]))

