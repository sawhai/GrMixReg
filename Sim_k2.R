source('modules/GMR_data_gen4.R')
source('modules/fit_GMR2.R')
#source('modules/network_commons.R')
#library(flexmix)

meanColSSQ <- function(X,Y) mean(colSums((X-Y)^2))
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))

# #Function to predict by mean -------
# #Given true y (training), the prediction is mean(y_training)
# predict_by_mean <- function(true_ytrain,ytest){
#   rep(mean(true_ytrain),length(ytest))
# }
# 
# #Function to predict test data using a linear regression model ---------
# single_lm <- function(tr.data,tst.data){
#   cols <- c(paste0('x',1:(ncol(tst.data)-4)),'Y')
#   tr.data <- tr.data[, ..cols]
#   #tr.dat <- tr.data[,c('x1','x2','x3','x4','Y'),with=F]
#   lm.model <- lm(Y ~ .,tr.data)
#   
#   predict(lm.model,tst.data)
# }

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

# Execute main simulation loop
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



