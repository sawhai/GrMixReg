rm(list=ls())
library(data.table)
library(dplyr)
library(MASS)
library(flexmix)
library(ggplot2)
library(latex2exp)


source('modules/GMR_data_gen.R')
source('modules/fit_GMR.R')
source('modules/network_commons.R')
source('modules/MMCL_fns.R')

d <- 2            # The dimension of data, i.e., number of covariates.
K <- 2            # The number of components. This code needs K <= d+1, since we generate equidistant components beta's
nobs <- 200       # Total number of observations
test_perc <- 0.2  # Percentage of data used for testing

total_num_runs <- 20  # number of replications. Reduce to speed up the simulation.

# combinations used in simulation
runs <- expand.grid(run_id=1:total_num_runs, 
                    bet_dist=8, 
                    noise_lev=seq(2,10,2),  # Noise level
                    G=5,                    # Number of groups in each component
                    nobs=nobs,
                    K=K,
                    d=d)

rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))
#####################
keep_cols <- c(paste0('x',1:d),'Y')

for (r in 1:nrow(runs)) {
  cat(sprintf('Run %4d out of %d\r', r, nrow(runs))) 
  run <- runs[r,]
  
  K <- run$K
  G <- run$G
  nobs <- run$nobs
  R <- K*G          # Total number of groups
  nr <- nobs / R    # Number of observ. per group
  
  out <- data_gen(K, nobs, G, run$bet_dist, run$d, run$noise_lev, VERB=F)
  dat <- out$data
  tru_bets <- out$bets
  
  # Split into training and test
  # Randomly pick %test_perc of observations from each group (idx) for testing
  nr_tst <- round(nr * test_perc)
  if (nr_tst == 0) {
    warning(sprintf('%%%2.1f of observations is zero. Choosing at least 1 observations for test.',test_perc*100))
    nr_tst <- 1
  }
  tst_idx <- sample(1:nr,nr_tst)
  #nr_tr <- nr - nr_tst
  dat_tr = dat[,.SD[!tst_idx], by = idx] #Training data 
  dat_tst = dat[,.SD[tst_idx], by = idx] #Test data
  setcolorder(dat_tr,names(dat)) #Change the order of the columns
  setcolorder(dat_tst,names(dat))
  tru_obs_labels <- dat_tr[,tru.label] #  true obs. labels :note the change to label"s" from label
  tru_grp_labels <- out$grp_labels # true group labels
  
  # Fit GMR
  fit <- fit_grp_mix_reg(dat_tr, K=K, d=d, VERB=F)
  
  est_grp_labels <- label_mat2vec(fit$tau) # cluster assignment for groups
  est_obs_labels <- est_grp_labels[dat_tr$idx] # cluster assignment for individual obs.
  est_bets <- do.call(cbind,fit$beta) # estimated beta's
  
  tau <- fit$tau
  predict_gmr(dat_tst, tau, est_bets)
  runs[r,"rmse"] <- rmse( dat_tst[, Yh],  dat_tst[, Mu] )
  
  runs[r,"nmi"] <-  compute_mutual_info(tru_obs_labels, est_obs_labels)
  
  #MMCL------
  mmcl_fit <- fit_mmcl(dat_tr,K,R,d,keep_cols,10)
  #######################################################################
  mmcl_est_labels <- mmcl_fit$result$Best.rmse.Model 
  runs[r,"nmi_mmcl"] <-  compute_mutual_info(tru_grp_labels, mmcl_est_labels)
  tau_mmcl <- matrix(rep(0,K*R),ncol=2)
  for(i in 1:K){
    tau_mmcl[which(mmcl_est_labels==i),i] <- 1
  }
  dat_tst[,Yh:=NULL]
  mmcl_est_bets <- matrix(unlist(mmcl_fit$betah),nrow = K)
  predict_gmr(dat_tst,tau_mmcl,mmcl_est_bets)
  runs[r,"rmse_mmcl"] <- rmse( dat_tst[, Yh],  dat_tst[, Mu] )
}

# Calculate the averages
runs2 <- runs %>%
  mutate(bet_dist = factor(bet_dist)) %>%
  group_by(noise_lev) %>% 
  summarize(avg_nmi = mean(nmi, na.rm = T), 
            avg_mmcl_nmi = mean(nmi_mmcl, na.rm = T),
            avg_rmse = mean(rmse, na.rm = T),
            avg_mmcl_rmse = mean(rmse_mmcl, na.rm = T)) %>% 
  gather(nmi, nmi_val, -c(1,4,5)) %>%
  gather(rmse, rmse_val, 2:3)


xlab = '$\\sigma_k$'
mmcl_ggplot(runs2, aes(noise_lev, rmse_val, color=rmse), title='Average RMSE', xlab=xlab,c("MMCL++", "GMR"))
mmcl_ggplot(runs2, aes(noise_lev, nmi_val, color=nmi), title='Average NMI', xlab=xlab,c("MMCL++", "GMR"))
