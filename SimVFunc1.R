# =================================================================================================
# Purpose: 
# 1) generate simulated true cohort, misclassified cohort and validation cohort 
# comprised of true exposure, observed exposure and true outcome values
# 20% true exposed prevalence (E) and 80% true unexposed prevalence
# 2% cumulative incidence of disease in an exposed group
# 1% cumulative incidence of disease in an unexposed group (R)
# assume 60% sensitivity (Se) and 98% specificity (Sp) of exposure classification
# 2) compare 5 bias-adjusted scenarios when the beta distributions are 
#    a. informed using the observed validation data directly 
#    b. centered on the mean (i.e., with a prior of beta (α = 1, β =1)) or on the mode (i.e., without a prior) 
#    c. with or without continuity correction solutions to the zero-cell problem in a validation sub-study with sparse data.
# 
# Echo 10/25/2021 
# history of update:
# updated 12/13/2021 1) add the random error to adjusted RRs
#                    2) add the coverage probability: for adjusted RRs, use MLE as 
#                       point estimate and use Marshall's variance, resulting in 
#                       60% coverage probability, which is not used
# updated 12/23/2021 add a sub-loop that uses each iteration with 100 iterations to 
#                    sample from the beta distributions for the PPV/NPV (cases and non-cases), 
#                    and calculate bias adjusted result with resampled conventional standard error
# updated 12/28/2021 change the sub-loop from mean to median for both result and variance
# updated 2/9/2022 continuity correction for the zero cells in validation cohort; 
#                  drop the finite and infinite values
# updated 3/30/2022 vectorize the code instead of looping
# updated 5/1/2022 add more scenarios: 
# a.population size 100k, exposure prevalence 20%, risk 0.1, RR=2, SE=0.6, SP=0.98, 
# 100K iterations each with 10K subiterations, and for each iterations, inform PPV/NPV as 
# (1) validation data as is, (2) validation data +1 (always), (3) validation data +0.5 (always), 
# (4) validation data +0.5 (only if 0 cell in validation data), (5) validation data +1 (only if  0 cell in validation data).
# b.	same five scenarios as above, except set population size to 10K instead of 100K
# this script is the (1) validation data as is
# updated 5/24/2022 corrected an error in the function adjest.fun:
# when calculating D0E0_adj1 and D1E0_adj1: the dimensions of variables are inconsistent
# updated 7/11/2022 1) drop the Step 4b (bias-adjusted estimate with prior) as this has been included in other scenarios;
#                   2) add the bias (difference between median adjusted RR and truth on the log scale) for each of the 100,000 iterations;
#                   3) add the MSE (mean squared error for 100,000 iterations) using the median point estimate 
#                      from each of the 100,000 iterations and the truth on the log scale;
#                   4) add the variables for original alphas and betas for npv and ppv of cases and non-cases
#                   5) add more parameters (mean, min, max, 95% interval) to the list saved to results of iterations
# =================================================================================================
# This is scenario 1: the conventional method based on validation data without a prior or a continuity correction
set.seed(789)
#install.packages("data.table")
library(data.table) # setDT function
#install.packages("grr")
library(grr) # matches function
library(tidyverse)
#library(tidyr)
#install.packages("openxlsx")
library(openxlsx)
#install.packages("pryr") # for memory check mem_used() mem_change()
library(pryr)

SimVFunc <- function(E, R, RR, Se, Sp, n, num.rep, num.subiter){
  # create a matrix for the final results
  # the number of columns and column names might need to be changed if revised
  final <- data.frame(matrix(NA, length(num.rep)*length(num.subiter), 39))
  names(final) <- c("E", "n", "iter", "sub_iter", 
                    as.vector(t(outer(c("convRR", "convV", "misRR", "misV", 
                                        "adj1RR_med", "adj1V_med", "adj1bias", "SE", "SP", "min"),
                                      c("_median", "_LCL95", "_UCL95"), paste0))),
                    "count", "conv_Coverage", "mis_Coverage", "adj1_Coverage", "adj1_MSE")

  unexpD = R
  expD = R*RR
  
  SeAlpha <- 100*Se
  SeBeta <- 100*(1-Se)
  
  SpAlpha <- 100*Sp
  SpBeta <- 100*(1-Sp)
  
  for (k in seq_along(num.rep)) {
    for (h in seq_along(num.subiter)){
      
        # Step 1: simulate true cohort
        ## 2 * 2 table for true cohort
        D1E1 <- rbinom(num.rep[k], n*E, expD) # or: qbinom(runif(1), n*E, expD)
        D0E1 <- n*E - D1E1
        D1E0 <- rbinom(num.rep[k], n*(1-E), unexpD)
        D0E0 <- n*(1-E) - D1E0
        
        ## conventional true RR
        convRR <- (D1E1/(n*E))/(D1E0/(n*(1-E)))
        convV <- (D0E1/D1E1)/(n*E) + (D0E0/D1E0)/(n*(1-E))
        ## check if the Wald-based 95% CI cover the truth
        convLower <- exp(log(convRR) - 1.96*sqrt(convV))
        convUpper <- exp(log(convRR) + 1.96*sqrt(convV))
        convCoverage <- ifelse((expD/unexpD >= convLower) & (expD/unexpD <= convUpper), 1, 0)
        
        # Step 2: simulate misclassified cohort
        ## randomly draw a value from the sensitivity and specificity distribution
        ## assume non-differential misclassification
        SeR <- rbeta(num.rep[k], SeAlpha, SeBeta) # rbeta can deal with vectors
        SpR <- rbeta(num.rep[k], SpAlpha, SpBeta) 
        
        ## true/classified numbers for each group 
        D1_E1E1 <- rbinom(num.rep[k], D1E1, SeR)
        D1_E1E0 <- D1E1 - D1_E1E1
        D1_E0E0 <- rbinom(num.rep[k], D1E0, SpR)
        D1_E0E1 <- D1E0 - D1_E0E0
        
        D0_E1E1 <- rbinom(num.rep[k], D0E1, SeR)
        D0_E1E0 <- D0E1 - D0_E1E1
        D0_E0E0 <- rbinom(num.rep[k], D0E0, SpR)
        D0_E0E1 <- D0E0 - D0_E0E0
        
        misD1E1 <- D1_E1E1 + D1_E0E1 
        misD1E0 <- D1_E1E0 + D1_E0E0
        misD0E1 <- D0_E1E1 + D0_E0E1
        misD0E0 <- D0_E1E0 + D0_E0E0 
        
        misRR <- (misD1E1/(misD1E1+misD0E1))/(misD1E0/(misD1E0+misD0E0))
        misV <- (misD0E1/misD1E1)/(misD0E1+misD1E1) + (misD0E0/misD1E0)/(misD0E0+misD1E0)
        ## check if the Wald-based 95% CI cover the truth
        misLower <- exp(log(misRR) - 1.96*sqrt(misV))
        misUpper <- exp(log(misRR) + 1.96*sqrt(misV))
        misCoverage <- ifelse((expD/unexpD >= misLower) & (expD/unexpD <= misUpper), 1, 0)
        
        # Step 3: validation cohort, balanced design
        ## minimum interior cell frequency
        sample <- mapply(min, misD1E1, misD1E0, misD0E1, misD0E0)
        ## for cases
        # validation data as is (leave alone 0 cells)
        ppvD1Alpha <- ifelse(misD1E1==sample, D1_E1E1, mapply(rbinom,1,sample, D1_E1E1/(D1_E1E1+D1_E0E1)))
        ppvD1Beta <- sample - ppvD1Alpha
        ppvD1Alpha_new <- ppvD1Alpha
        ppvD1Beta_new <- ppvD1Beta
        
        npvD1Beta <- ifelse(misD1E0==sample, D1_E1E0, mapply(rbinom,1,sample, D1_E1E0/(D1_E1E0+D1_E0E0)))
        npvD1Alpha <- sample - npvD1Beta
        npvD1Alpha_new <- npvD1Alpha
        npvD1Beta_new <- npvD1Beta
        
        ## for non-cases
        # validation data as is (leave alone 0 cells)
        ppvD0Alpha <- ifelse(misD0E1==sample, D0_E1E1, mapply(rbinom,1,sample, D0_E1E1/(D0_E1E1+D0_E0E1)))
        ppvD0Beta <- sample - ppvD0Alpha
        ppvD0Alpha_new <- ppvD0Alpha
        ppvD0Beta_new <- ppvD0Beta
        
        npvD0Beta <- ifelse(misD0E0==sample, D0_E1E0, mapply(rbinom,1,sample, D0_E1E0/(D0_E1E0+D0_E0E0)))
        npvD0Alpha <- sample - npvD0Beta
        npvD0Alpha_new <- npvD0Alpha
        npvD0Beta_new <- npvD0Beta
        
        rm(D1E1, D0E1, D1E0, D0E0)
        #print(memory.size(max = T))
        print(mem_used())
        
        # Step 4a: bias-adjusted estimate, no prior 
        
        adjest1 <- adjest.fun(num.rep[k], num.subiter[h], convRR, expD, unexpD, ppvD1Alpha_new, ppvD1Beta_new, npvD1Alpha_new, npvD1Beta_new, ppvD0Alpha_new, ppvD0Beta_new, npvD0Alpha_new, npvD0Beta_new,
                                 misD1E1, misD1E0, misD0E1, misD0E0)
        adj1RR_med <- adjest1[[1]]
        adj1RR_min <- adjest1[[2]]
        adj1RR_max <- adjest1[[3]]
        adj1Lower <- adjest1[[4]]
        adj1Upper <- adjest1[[5]]
        adj1RR_avg <- adjest1[[6]]
        adj1V_med <- adjest1[[7]]
        adj1V_avg <- adjest1[[8]]
        adj1Coverage2 <- adjest1[[9]]
        
        rm(adjest1)
        #print(memory.size(max = T))
        print(mem_used())
        
        # # Step 4b: bias-adjusted estimate, with prior (deleted, not needed)
      
      result <- data.frame(convRR, convV, misRR, misV, 
                           adj1RR_med, adj1V_med, adj1bias = log(adj1RR_med) - log(expD/unexpD),
                           SE = SeR, SP = SpR, min = sample, 
                           ppvD1Alpha_new, ppvD1Beta_new, npvD1Alpha_new, npvD1Beta_new,ppvD0Alpha_new, ppvD0Beta_new, npvD0Alpha_new, npvD0Beta_new) #these npv and ppv are old values
      # complete.cases only exclude NA and NaN, but using is.finite also exclude Inf and -Inf
      summary <- sapply(result[is.finite(rowSums(result)), ], summary.fun)
      summary <- setDT(as.data.frame(summary), keep.rownames = "parameter")[]
      summary$parameter <- replace(summary$parameter, summary$parameter == "LCL95.2.5%", "LCL95")
      summary$parameter <- replace(summary$parameter, summary$parameter == "UCL95.97.5%", "UCL95")
      summary$index <-ave(seq_along(summary$parameter), summary$parameter, FUN=seq_along)
      summary <- dcast(summary, index~factor(parameter, levels = unique(parameter)), value.var = c("convRR", "convV", "misRR", "misV", "adj1RR_med", "adj1V_med",
                                                                                                  "adj1bias", "SE", "SP", "min"))
      # get the min of all count variables
      summary <- summary %>% mutate(count = apply(select(summary, matches("(_count)$")), 1, min, na.rm = TRUE))
      summary[ , grep("_count", names(summary))] <- list(NULL)
      summary <- subset(summary,select = -index)
      coverage <- data.frame(conv_Coverage = mean(convCoverage, na.rm = TRUE), mis_Coverage = mean(misCoverage, na.rm = TRUE), adj1_Coverage = mean(adj1Coverage2, na.rm = TRUE))
      adj1mse <- mse(log(adj1RR_med), log(expD/unexpD))
      final[(k-1)*length(num.subiter) + h,] <- cbind(E, n, num.rep[k], num.subiter[h], summary, coverage, adj1mse)
    }
      
  }
  result <- data.frame(convRR, convLower, convUpper, convV, misRR, misLower, misUpper, misV, adj1RR_med, adj1RR_min, adj1RR_max, adj1Lower, adj1Upper, adj1RR_avg,
                       adj1V_med, adj1V_avg, adj1bias = log(adj1RR_med) - log(expD/unexpD),
                       SE = SeR, SP = SpR, min = sample, 
                       ppvD1Alpha_new, ppvD1Beta_new, npvD1Alpha_new, npvD1Beta_new,ppvD0Alpha_new, ppvD0Beta_new, npvD0Alpha_new, npvD0Beta_new)
  final <- list(final,result) # to see the values for each iteration
  print(mem_used())
  return(final)
}

summary.fun <- function(x){
  c(median = median(x, na.rm = TRUE), LCL95 = quantile(x, 0.025, na.rm = TRUE), UCL95 = quantile(x, 0.975, na.rm = TRUE), count = sum(!is.na(x)))
}

adjest.fun <- function(num.rep, num.subiter, convRR, expD, unexpD, ppvD1Alpha, ppvD1Beta, npvD1Alpha, npvD1Beta, ppvD0Alpha, ppvD0Beta, npvD0Alpha, npvD0Beta,
                       misD1E1, misD1E0, misD0E1, misD0E0){
  ## run num.subiter iterations to sample from the ppv and npv distribution
  ppvD1_adj1R <- npvD1_adj1R <- ppvD0_adj1R <- npvD0_adj1R <- D1E1_adj1 <- 
    D1E0_adj1 <- D0E1_adj1 <- D0E0_adj1 <- adj1RR <- adj1V <- adj1RR_nr <- c()
  # mapply takes either lists or vectors, but not matrix
  ppvD1_adj1R <- mapply(rbeta, num.subiter, ppvD1Alpha, ppvD1Beta) # a matrix, where nrow = #sub-iteration, ncol = #iteration
  npvD1_adj1R <- mapply(rbeta, num.subiter, npvD1Alpha, npvD1Beta)
  
  ppvD0_adj1R <- mapply(rbeta, num.subiter, ppvD0Alpha, ppvD0Beta)
  npvD0_adj1R <- mapply(rbeta, num.subiter, npvD0Alpha, npvD0Beta)
  
  
  D1E1_adj1 <- adjcell(num.rep, num.subiter,misD1E1, ppvD1_adj1R) + adjcell(num.rep, num.subiter,misD1E0, 1-npvD1_adj1R)
  rm(ppvD1_adj1R, npvD1_adj1R)
  misD1E1 <- matrix(misD1E1, nrow=num.subiter, ncol = num.rep, byrow=TRUE)
  misD1E0 <- matrix(misD1E0, nrow=num.subiter, ncol = num.rep, byrow=TRUE)
  D1E0_adj1 <- misD1E1 + misD1E0 - D1E1_adj1
  rm(misD1E1, misD1E0)
  
  
  D0E1_adj1 <- adjcell(num.rep, num.subiter,misD0E1, ppvD0_adj1R) + adjcell(num.rep, num.subiter,misD0E0, 1-npvD0_adj1R)
  rm(ppvD0_adj1R, npvD0_adj1R)
  misD0E1 <- matrix(misD0E1, nrow=num.subiter, ncol = num.rep, byrow=TRUE)
  misD0E0 <- matrix(misD0E0, nrow=num.subiter, ncol = num.rep, byrow=TRUE)
  D0E0_adj1 <- misD0E1 + misD0E0 - D0E1_adj1
  rm(misD0E1, misD0E0)
  
  ## calculate adjusted RR without including random error
  adj1RR_nr <- (D1E1_adj1/(D1E1_adj1+D0E1_adj1)) / (D1E0_adj1/(D1E0_adj1+D0E0_adj1))
  adj1V <- (D0E1_adj1/D1E1_adj1)/(D0E1_adj1+D1E1_adj1) + (D0E0_adj1/D1E0_adj1)/(D0E0_adj1+D1E0_adj1)
  rm(D1E1_adj1, D1E0_adj1, D0E1_adj1, D0E0_adj1)
  ## add random error (which is retrieved by sampling a z value from standard normal distribution) to the adjusted RR
  rand_err <- mapply(rnorm, rep(num.subiter,num.rep), rep(0,num.rep), rep(1,num.rep))
  adj1RR <- exp(log(adj1RR_nr) - rand_err*sqrt(adj1V))
  rm(adj1RR_nr, rand_err)
  adj1RR[!is.finite(adj1RR)] <- NA
  adj1RR_med <- apply(adj1RR, 2, median, na.rm = TRUE)
  adj1RR_min <- apply(adj1RR, 2, min, na.rm = TRUE)
  adj1RR_max <- apply(adj1RR, 2, max, na.rm = TRUE)
  adj1RR_avg <- apply(adj1RR, 2, mean, na.rm = TRUE)
  adj1V[!is.finite(adj1V)] <- NA
  adj1V_med <- apply(adj1V, 2, median, na.rm = TRUE)
  adj1V_avg <- apply(adj1V, 2, mean, na.rm = TRUE)
  rm(adj1V)
  adj1Lower <- apply(adj1RR, 2, quantile, 0.025, na.rm = TRUE)
  adj1Upper <- apply(adj1RR, 2, quantile, 0.975, na.rm = TRUE)
  # calculate the MSE for the 10,000 sub-iterations
  
  # adj1mse <- mapply(mse, data.frame(adj1RR), convRR) # a vector of mse values
  rm(adj1RR)
  adj1Coverage <- ifelse((expD/unexpD >= adj1Lower) & (expD/unexpD <= adj1Upper), 1, 0) 
  return(list(adj1RR_med, adj1RR_min, adj1RR_max, adj1Lower, adj1Upper, adj1RR_avg, adj1V_med, adj1V_avg, adj1Coverage))
}

adjcell <- function(num.rep, num.subiter, misD1E1, ppvD1_adj1R){
  misD1E1 <- matrix(misD1E1, nrow=num.subiter, ncol = num.rep, byrow=TRUE)
  cell_vector <- mapply(rbinom, 1, misD1E1, ppvD1_adj1R) # result is simplified into a vector by column
  return(matrix((cell_vector), nrow=num.subiter, ncol = num.rep)) # convert the vector back to matrix by column
}

mse <- function(log_adjRR, log_trueRR){
  return(mean((log_adjRR-log_trueRR)^2,na.rm = TRUE))
}
#---------------------------------------------------------------------------------------------------------------#
# run function
# population size 100k
final <- SimVFunc(E = 0.2, R = 0.01, RR = 2, Se = 0.6,  Sp = 0.98, n = 100000, num.rep = 100000, num.subiter =10000)
write.xlsx(final[[1]], file = paste0("S1_result100k_", Sys.Date(), ".xlsx"))
write.xlsx(final[[2]], file = paste0("S1_Data100k_", Sys.Date(), ".xlsx"))
# population size 10k
final2 <- SimVFunc(E = 0.2, R = 0.01, RR = 2, Se = 0.6,  Sp = 0.98, n = 10000, num.rep = 100000, num.subiter =10000)
write.xlsx(final2[[1]], file = paste0("S1_result10k_", Sys.Date(), ".xlsx"))
write.xlsx(final2[[2]], file = paste0("S1_Data10k_", Sys.Date(), ".xlsx"))

