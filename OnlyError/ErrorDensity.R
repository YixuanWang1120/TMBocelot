density.error <- function(x, c, theta, sigma, alpha, var.obs){
  
  n <- length(x)
  cof <- 0.5
  
  sigma.a <- 9 / 4
  sigma.b <- 5 / 8 * var.obs
  
  # sigma.a <- 0.001
  # sigma.b <- 0.001
  
  
  c.class <- table(c)
  c.name <- as.numeric(dimnames(c.class)$c)
  c.number <- as.vector(c.class)
  
  c.theta <- rep(0, length(c.number))
  c.sigma <- rep(0, length(c.number))
  
  for(i in 1:length(c.name)){
    c.theta[i] <- theta[which(c == i)[1]]
    c.sigma[i] <- sigma[which(c == i)[1]]
  }
  
  
  for(i in 1:N){
    pre.c <- c[i]
    c.number[pre.c] <- c.number[pre.c] - 1
    
    p <- c(c.number, alpha) / (n - 1 + alpha)
    K <- length(c.name)
    curr.c <- sample(x = 1:(K+1), prob = p, size = 1)
    
    if(curr.c == K + 1){
      
      temp.sigma <- 1 / sqrt(rgamma(1, sigma.a, sigma.b))
      temp.theta <- rnorm(1, 0, sqrt(cof) * temp.sigma)
      accept.prop <- min(1, dnorm(x[i], temp.theta, temp.sigma) / dnorm(x[i], theta[i], sigma[i]))
      u <- runif(1, 0, 1)
      if(!is.na(accept.prop) && u < accept.prop){
        c[i] <- curr.c
        theta[i] <- temp.theta
        sigma[i] <- temp.sigma
        
        c.name <- c(c.name, K + 1)
        c.number <- c(c.number, 1)
        c.theta <- c(c.theta, temp.theta)
        c.sigma <- c(c.sigma, temp.sigma)
      }else{
        c.number[pre.c] <- c.number[pre.c] + 1
      }
      
    }else{
      
      if(curr.c == pre.c){
        c.number[pre.c] <- c.number[pre.c] + 1
      }else{
        
        temp.theta <- c.theta[curr.c]
        temp.sigma <- c.sigma[curr.c]
        
        accept.prop <- min(1, dnorm(x[i], temp.theta, temp.sigma) / dnorm(x[i], theta[i], sigma[i]))
        u <- runif(1, 0, 1)
        if(!is.nan(accept.prop) && u < accept.prop){
          c[i] <- curr.c
          theta[i] <- temp.theta
          sigma[i] <- temp.sigma
          
          c.number[curr.c] <- c.number[curr.c] + 1
        }else{
          c.number[pre.c] <- c.number[pre.c] + 1
        }
        
      }
    }
    
  }
  
  c.update <- c
  for(i in 1:length(c.number)){
    if(c.number[i] == 0){
      for(j in (i + 1):length(c.number)){
        index <- which(c == j)
        if(length(index) != 0){
          c.update[index] <- c.update[index] - 1
        }
        
      }
      
    }
  }
  
  c <- c.update
  
  unique.class <- unique(c)
  flag <- 0
  update.theta <- rep(0, n)
  update.sigma <- rep(0, n)
  
  for (i in 1:length(unique.class)){
    
    
    index <- which(c == unique.class[i])
    K <- length(index)
    
    # theta
    update.theta[index] <- temp.theta <- 
      rnorm(1, cof * sum(x[index]) / (K * cof + 1), sqrt(cof / (K * cof + 1)) * sigma[index[1]])
    
    #sigma
    temp <- 0
    for (j in index) {
      temp <- temp + (temp.theta - x[j])^2
    }
    update.sigma[index] <- 1 / sqrt(rgamma(1, sigma.a + length(index) / 2 + 1 / 2, sigma.b + temp / 2 + temp.theta^2 / (2 * cof)))
  }
  theta <- update.theta
  sigma <- update.sigma
  
  k <- length(unique.class)
  a <- 0.001
  b <- 0.001
  med <- rbeta(1, alpha + 1, n)
  sel <- sample(x = c(0, 1), prob = c(a + k - 1, n * (b - log(med))), size = 1)
  if(sel == 0){
    alpha <- rgamma(1, a + k, b - log(med))
  }else{
    alpha <- rgamma(1, a + k - 1, b - log(med))
  }
  
  list(c = c, theta = theta, sigma = sigma, alpha = alpha)
}