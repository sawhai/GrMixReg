library(dplyr)
# library(clue)
library(data.table)
library(ggplot2)
library(latex2exp)
source('modules/GMR_data_gen4.R')
source('modules/fit_GMR.R')
source('modules/network_commons.R')

d <- 2            # The dimension of data, i.e., number of covariates.
K <- 2            # The number of components. This code needs K <= d+1, since we generate equidistant components beta's
nobs <- 600       # Total number of observations
G <- 10           # Number of groups in each component
Rtot <- K*G       # Total number of groups
nr <- nobs / Rtot # Number of observ. per group
test_perc <- 0.2 # Percentage of data used for testing

meanColSSQ <- function(X,Y) mean(colSums((X-Y)^2))
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))

# combinations used in simulation
total_num_runs <- 5
runs <- expand.grid(run_id=1:total_num_runs, bet_dist=c(4,8,12), noise_lev=seq(2,10,2))

# Run the simulations
dt <- system.time(
for (r in 1:nrow(runs)) {
      cat(sprintf('Run %4d out of %d\n', r, nrow(runs))) 
      run <- runs[r,]
      #out <- data_gen(K, N, R, run$bet_dist, d, run$noise_lev, VERB=F)
      out <- data_gen(K, nobs, G, run$bet_dist, d, run$noise_lev, VERB=F)
      dat <- out$data
      tru_bets <- out$bets
      # Split into training and test
      # Randomly pick 6 observations from each group (idx) for testing
      nr_tst <- round(nr * test_perc)
      tst_idx <- sample(1:nr,nr_tst)
      nr_tr <- nr - nr_tst
      dat_tr = dat[,.SD[!tst_idx], by = idx] #Training data 
      dat_tst = dat[,.SD[tst_idx], by = idx] #Test data
      setcolorder(dat_tr,names(dat)) #Change the order of the columns
      setcolorder(dat_tst,names(dat))
      tru_obs_labels <- dat_tr[,tru.label] #  true obs. labels :note the change to label"s" from label
      tru_grp_labels <- out$grp_labels # true group labels
      dat_tr[, tru.label := NULL]
      dat_tst[, tru.label := NULL]

      # Fit gmR
      fit <- fit_grp_mix_reg(dat_tr, K=K, d=d, n.gr=Rtot, nr=nr_tr, VERB=F)
      # fit <- fit_grp_mix_reg(dat, K=K, d=d, n.gr=Rtot, nr=nr_tr, VERB=F)
      
      est_grp_labels <- label_mat2vec(fit$tau) # cluster assignment for groups
      est_obs_labels <- est_grp_labels[dat_tr$idx] # cluster assignment for individual obs.
      est_bets <- do.call(cbind,fit$beta) # estimated beta's
      
      tau <- fit$tau
      predict_gmr(dat_tst, tau, est_bets)
      runs[r,"rmse"] <- rmse( dat_tst[, Yh],  dat_tst[, Y] )
      
      # Fit simple lm model as a baseline
      lm_fit <- lm(as.formula(paste("Y~ 0 +", paste(paste0("x", 1:d), collapse="+"))), dat_tr)
      runs[r,"rmse_lm"] <- rmse( predict(lm_fit, dat_tst),  dat_tst[, Y] )
      
      runs[r,"nmi"] <-  compute_mutual_info(tru_obs_labels, est_obs_labels)
      # runs[r,"nmi"] <-  compute_mutual_info(tru_grp_labels, est_grp_labels) 
      runs[r,"n_iter"] <- fit$n.itr
      runs[r,"beta_err"] <- meanColSSQ(tru_bets[ , tru_grp_labels], est_bets[ , est_grp_labels])
}
)["elapsed"]

cat("Total runtime = ", dt, "(s)\n")

# Calculate the averages
runs2 <- runs %>%
  mutate(bet_dist = factor(bet_dist)) %>%
  group_by(bet_dist, noise_lev) %>% 
  summarize(avg_nmi = mean(nmi, na.rm = T), 
            avg_n_iter=mean(n_iter, na.rm = T),
            avg_beta_err = mean(beta_err, na.rm = T),
            avg_rmse = mean(rmse, na.rm = T),
            avg_rmse_lm = mean(rmse_lm, na.rm = T))

# Plot
custom_ggplot <- function( data, aesth_map, title) {
  ggplot(data, aesth_map) + 
    geom_point(size=4.5,shape=1) +
    geom_line(lty=2) +
    theme_bw() + 
    labs(x=TeX('$\\sigma_k$'), y=TeX(title))+
    theme(text = element_text(size=25),
          panel.grid.major=element_line(colour='gray75'),
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15))+
    scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))#+scale_y_continuous(limits = c(0,1.5))
}

custom_ggplot(runs2, aes(noise_lev, avg_nmi, color=bet_dist), 'Average NMI' )
#ggsave('nmi.pdf')
custom_ggplot(runs2, aes(noise_lev, avg_rmse, color=bet_dist), 'Average RMSE' )
#ggsave('rmse.pdf')
custom_ggplot(runs2, aes(noise_lev, avg_n_iter, color=bet_dist), 'Average Number of Iterations' ) + scale_y_continuous(trans='log10')
#ggsave('niter.pdf')
custom_ggplot(runs2, aes(noise_lev, avg_beta_err, color=bet_dist), 'Average error $\\beta$' ) # + scale_y_continuous(trans='log10')
#ggsave('beta_err')
custom_ggplot(runs2, aes(noise_lev, avg_rmse_lm, color=bet_dist), 'Average RMSE (LM)' )
