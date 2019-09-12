source('modules/GMR_data_gen4.R')
source('modules/fit_GMR2.R')
source('modules/network_commons.R')

meanColSSQ <- function(X,Y) mean(colSums((X-Y)^2))
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))

# Execute main simulation loop
run_simulations <- function(runs, test_perc=0.2) {
  cat('Running simullations...\n')
  dt <- system.time(
    for (r in 1:nrow(runs)) {
      cat(sprintf('Run %4d out of %d\r', r, nrow(runs))) 
      run <- runs[r,]
      
      K <- run$K
      G <- run$G
      nobs <- run$nobs
      R <- K*G          # Total number of groups
      nr <- nobs / R    # Number of observ. per group
      
      #out <- data_gen(K, N, R, run$bet_dist, d, run$noise_lev, VERB=F)
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
      dat_tr[, tru.label := NULL]
      dat_tst[, tru.label := NULL]
      
      # Fit gmR
      fit <- fit_grp_mix_reg(dat_tr, K=K, d=d, VERB=F)
      # fit <- fit_grp_mix_reg(dat, K=K, d=d, n.gr=R, nr=nr_tr, VERB=F)
      
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
  
  runs
}


# Plot
custom_ggplot <- function( data, aesth_map, title, xlab) {
  ggplot(data, aesth_map) + 
    geom_point(size=4.5,shape=1) +
    geom_line(lty=2) +
    theme_bw() + 
    labs(x=TeX(xlab), y=TeX(title))+
    theme(text = element_text(size=25),
          panel.grid.major=element_line(colour='gray75'),
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15))+
    scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3)) # +scale_y_continuous(limits = c(0,1.5))
}

