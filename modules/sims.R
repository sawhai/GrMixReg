library(data.table)
library(ggplot2)
library(latex2exp)
library(flexmix)
source('modules/GMR_data_gen.R')
source('modules/fit_GMR.R')
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
      #dat_tr[, tru.label := NULL]
      #dat_tst[, tru.label := NULL]
      
      # Fit GMR
      fit <- fit_grp_mix_reg(dat_tr, K=K, d=d, VERB=F)
      
      est_grp_labels <- label_mat2vec(fit$tau) # cluster assignment for groups
      est_obs_labels <- est_grp_labels[dat_tr$idx] # cluster assignment for individual obs.
      est_bets <- do.call(cbind,fit$beta) # estimated beta's
      
      tau <- fit$tau
      predict_gmr(dat_tst, tau, est_bets)
      runs[r,"rmse"] <- rmse( dat_tst[, Yh],  dat_tst[, Mu] )
      
      runs[r,"nmi"] <-  compute_mutual_info(tru_obs_labels, est_obs_labels)
      # runs[r,"nmi"] <-  compute_mutual_info(tru_grp_labels, est_grp_labels) 
      runs[r,"n_iter"] <- fit$n.itr
      runs[r,"beta_err"] <- meanColSSQ(tru_bets[ , tru_grp_labels], est_bets[ , est_grp_labels])
      
      # Fit simple lm model as a baseline
      lm_fit <- lm(as.formula(paste("Y~ 0 +", paste(paste0("x", 1:d), collapse="+"))), dat_tr)
      runs[r,"rmse_lm"] <- rmse( predict(lm_fit, dat_tst),  dat_tst$Mu )
      
      fmr_fit <- flexmix(as.formula(paste("Y~ 0 +", paste(paste0("x", 1:d), collapse="+"))), dat_tr, k = run$K)
      bet_fmr <- parameters(fmr_fit)[1:K,]
      
      runs[r,"nmi_fmr"] <- compute_mutual_info(tru_obs_labels, clusters(fmr_fit))
      
      xcol_names <- paste0("x", 1:d)
      dat_tst_df <- as.data.frame(dat_tst[, c(xcol_names,'Y'),with=F])
      dat_tst$fmr_clust <- clusters(fmr_fit, newdata = dat_tst_df)
      
      
      for(u in 1:nrow(dat_tst)){
       dat_tst[u,'yh_fmr'] <- as.matrix(dat_tst[u, xcol_names, with=F]) %*% bet_fmr[, dat_tst$fmr_clust[u]] 
      }
      
      runs[r,"rmse_fmr"] <- rmse(dat_tst$yh_fmr,  dat_tst$Mu )
    }
  )["elapsed"]
  cat("Total runtime = ", dt, "(s)\n")
  
  runs
}

# use either mean, LM or GMR for prediction
predict_for_k <- function(K, dat_tr, dat_tst, d) {
  if (K ==0) {
    yh <- mean(dat_tr[, Y]) # predict_by_mean(dat_tr[, Y], dat_tst[, Y])
  } else if (K ==1) {
    xcols <- paste0('x',1:d)
    cols <- c(xcols,'Y')
    lm.model <- lm(Y ~ 0 + ., dat_tr[, ..cols])
    yh <- predict(lm.model, dat_tst[, ..xcols])
  } else {
    # Fit GMR
    fit <- fit_grp_mix_reg(dat_tr, K=K, d=d, VERB=F)
    est_bets <- do.call(cbind,fit$beta) # estimated beta's
    tau <- fit$tau
    predict_gmr(dat_tst, tau, est_bets)
    yh <- dat_tst[, Yh]
  }
  
  rmse( yh,  dat_tst[, Y] )
}

# Execute main loop for testing "finding optimal K"
find_optimal_k_cv <- function(runs, test_perc=0.2) {
  cat('Running simullations...\n')
  dt <- system.time(
    for (r in 1:nrow(runs)) {
      cat(sprintf('Run %4d out of %d\r', r, nrow(runs)))
      run <- runs[r,]
      
      G <- run$G
      nobs <- run$nobs
      Ktru <- run$Ktru
      R <- Ktru*G       # Total number of groups
      nr <- nobs/R      # Number of observ. per group
      
      out <- data_gen(Ktru, nobs, G, run$bet_dist, run$d, run$noise_lev, VERB=F)
      #out <- data_gen(Ktru, nobs, G, bet_dist, d, noise_lev, VERB=F)
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
      dat_tr = dat[,.SD[!tst_idx], by = idx] #Training data 
      dat_tst = dat[,.SD[tst_idx], by = idx] #Test data
      setcolorder(dat_tr,names(dat)) #Change the order of the columns
      setcolorder(dat_tst,names(dat))
      
      runs[r,"rmse"] <- tryCatch(predict_for_k( run$Ktst, dat_tr, dat_tst, d=run$d ),  
                                 error=function(e){cat("ERROR : GMR failed to converge \n"); NA }) #conditionMessage(e), "\n")
      
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

