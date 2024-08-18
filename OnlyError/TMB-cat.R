library(stats)
library(survival)
TMB.cat.computeThres <- function(pro, TMB, time.T, status, R){
  
  pro_tmb <- cbind(pro, TMB)
  pro_tmb <- pro_tmb[order(pro_tmb[, 2]), ]
  
  F.list <- rep(-1, N)
  P.list <- rep(10, N)
  
  for(t in 2:N){
    
    lab <- c(rep('A', t - 1), rep('B', N - t + 1))
    pro_lab <- as.data.frame(cbind(pro_tmb[, 1], lab))
    names(pro_lab) <- c("prob", "level")
    aov.manu <- aov(prob ~ level, data=pro_lab)
    F.list[t] <- summary(aov.manu)[[1]][["F value"]][1]
    P.list[t] <- summary(aov.manu)[[1]][["Pr(>F)"]][1]
    
    # pro1 <- pro_tmb[1:(t - 1), 1]
    # pro2 <- pro_tmb[t : N, 1]
    # boxplot(pro1, pro2, main="Box Plot of Group 1 and Group 2", xlab="Group", ylab="Value")
    # t_test_result <- t.test(pro1, pro2)  
    # print(t_test_result)
  }
  index <- which(P.list == min(P.list))
  thres <- pro_tmb[index, 2]
  cat("TMB-cat thres: ", thres, "\n")
  final_thres <- thres
  
  
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
  list(thres = final_thres)
}