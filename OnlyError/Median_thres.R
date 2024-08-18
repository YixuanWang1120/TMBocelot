library(stats)
library(survival)
Median.computeThres <- function(TMB, time.T, status, R){
  
  final_thres <- median(TMB)
  cat("Median thres: ", final_thres, "\n")


  group <- rep(0, N)
  group <- ifelse(TMB >= final_thres, 1, 0)
  df.sur <- data.frame(
    time = time.T,
    event = status,
    group = group
  )
  
  log.rank <- survdiff(Surv(df.sur$time, df.sur$event) ~ df.sur$group)
  cat("Survival P-value: ", log.rank$pvalue, "\n")
  # 加载stats包
  
  fit <- survfit(Surv(df.sur$time, df.sur$event) ~ df.sur$group)
  print(fit)
  plot(fit)
  group1 <- group2 <- c(-1)
  for(l in 1:N){
    if(TMB[l] >= final_thres){
      group1 <- c(group1, R[l])
    }else{
      group2 <- c(group2, R[l])
    }
  }
  group1 <- group1[-1]
  group2 <- group2[-1]
  
  result.mannwhitneyu <- wilcox.test(group1, group2)
  
  cat("ORR P-value: ", result.mannwhitneyu$p.value, "\n")
  cat("Low: ", sum(group2) / length(group2), "\n")
  cat("Hign: ", sum(group1) / length(group1), "\n")
  
}
