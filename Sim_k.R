source('modules/GMR_data_gen4.R')
source('modules/fit_GMR2.R')
#source('modules/network_commons.R')
#library(flexmix)

meanColSSQ <- function(X,Y) mean(colSums((X-Y)^2))
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))

#Function to predict by mean -------
#Given true y (training), the prediction is mean(y_training)
predict_by_mean <- function(true_ytrain,ytest){
  yh <- rep(mean(true_ytrain),length(ytest))
  yh
}

#Function to predict test data using a linear regression model ---------
single_lm <- function(tr.data,tst.data){
  cols <- c(paste0('x',1:(ncol(tst.data)-4)),'Y')
  tr.data <- tr.data[, ..cols]
  #tr.dat <- tr.data[,c('x1','x2','x3','x4','Y'),with=F]
  lm.model <- lm(Y ~ .,tr.data)
  yh <- predict(lm.model,tst.data)
  yh
}

# Execute main simulation loop
optimal_k_cv <- function(runs,tru_k,min_k,max_k,test_perc=0.2) {
  cat('Running simullations...\n')
  r <- 1
  dt <- system.time(
    while (r <= nrow(runs)){
      run <- runs[r,]
      
      G <- run$G
      nobs <- run$nobs
      R <- K*G          # Total number of groups
      R_tru <- tru_k*G
      nr <- nobs / R    # Number of observ. per group
      nr_tru <- nobs/R_tru
      
      out <- data_gen(tru_k, nobs, G, run$bet_dist, run$d, run$noise_lev, VERB=F)
      #out <- data_gen(tru_k, nobs, G, bet_dist, d, noise_lev, VERB=F)
      dat <- out$data
      tru_bets <- out$bets
      
      # Split into training and test
      # Randomly pick %test_perc of observations from each group (idx) for testing
      nr_tst <- round(nr_tru * test_perc)
      if (nr_tst == 0) {
        warning(sprintf('%%%2.1f of observations is zero. Choosing at least 1 observations for test.',test_perc*100))
        nr_tst <- 1
      }
      tst_idx <- sample(1:nr_tru,nr_tst)
      #tst_idx <- c(nr-1,nr-2)
      #nr_tr <- nr - nr_tst
      dat_tr = dat[,.SD[!tst_idx], by = idx] #Training data 
      dat_tst = dat[,.SD[tst_idx], by = idx] #Test data
      setcolorder(dat_tr,names(dat)) #Change the order of the columns
      setcolorder(dat_tst,names(dat))

      for(n in min_k:max_k){
        #cat('Run ', r,'\n')
        tryCatch({
          cat(sprintf('Run %4d out of %d\r', r, nrow(runs)))
          
          K <- runs$K[r]
          if(K ==0 ){
            yh <- predict_by_mean(dat_tr[, Y], dat_tst[, Y])
            #runs[r,'K'] <- n
            runs[r,"rmse"] <- rmse( yh,  dat_tst[, Y] )
            r <- r +1
          }else if(K ==1){
            yh <- single_lm(dat_tr,dat_tst)
            #runs[t.r,'K'] <- n
            runs[r,"rmse"] <- rmse( yh,  dat_tst[, Y] )
            r <- r +1
          }else{ #For K in 2:8
            # Fit GMR
            fit <- fit_grp_mix_reg(dat_tr, K=K, d=d, VERB=F)

            est_bets <- do.call(cbind,fit$beta) # estimated beta's
            
            tau <- fit$tau
            predict_gmr(dat_tst, tau, est_bets)
            runs[r,"rmse"] <- rmse( dat_tst[, Yh],  dat_tst[, Y] ) 
            r <- r +1
          }
        }, error=function(e){cat("ERROR : GMR failed to converge \n")})
        
      }
      
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



