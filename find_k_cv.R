####################################################################################
# To find optimal K by Cross Validation
# This simulation generates synthetic data for K = 4 (and d = 4) and uses
# 80% of it to fit GMR for K = 2,3, ...,8. 20% is used to predict with the 
# trained model.
# It also predicts usin mean (K = 0), and using single linear regression (K=1)
#####################################################################################

rm(list=ls())
library(ggplot2)
library(data.table)
library(MASS)
require(gtools)
library(caret)
library(miscTools)
library(plyr)
library(matrixStats)
library(Metrics)
library(plotrix)
library(psych)
#Run the necessary scripts
source('modules/GMR_data_gen.R')
source('modules/fit_GMR.R')
source('modules/NMI.AMI.calc.fn.R')

K = 4 #Real K
N_all <- rep(3e2,K)  #Number of observations for each component
R <- rep(10 ,K) #Number of groups in each component
Rtot <- sum(R)  #Total number of groups
# 24 out of 30 (80%) observations per group for training, 6 for testing:
N <- rep(240,K) #Number of observations for each component for training
d <- 4
noise <- 6 #Noise level is set to 6
bet.dist <- c(4.6,6.9) #delta beta of 8 and 12
#K from 0 to 8 where 0 is prediction by mean and 1 is lm. 
#note: K is set to go from 0 to 6 since for K=7 and 8 the code is slow. It is because of the fact that to
# correct label switching, all the permutations of K should be checked and this number becomes large when K > 6
alter_k <- 0:6 
l1 <- length(alter_k); l2 <- length(bet.dist)
n.runs <- 20 #Number of runs

#Function to predict by mean -------
#Given true y (training), the prediction is mean(y_training)
predict_by_mean <- function(true_ytrain,ytest){
  yh <- rep(mean(true_ytrain),length(ytest))
  yh
}

#Function to predict test data using a linear regression model ---------
single_lm <- function(tr.data,tst.data){
  tr.dat <- tr.data[,c('X1','X2','X3','X4','Y'),with=F]
  lm.model <- lm(Y ~ .,tr.dat)
  yh <- predict(lm.model,tst.data)
  yh
}


#fin.res to hold the final result
fin.res <- data.table(N1 =rep(N[1],(n.runs*l1*l2)),N2 =N[2],R1 = R[1], R2 =R[2],Noise = noise, K = NaN,
                      Beta_dist=0,RMSE = 0)


t.r <- 1
res1 <- data.table(Groups= 1:Rtot, Clust = 0)
Conf <- list() #Holds the confusion matrices 
total_num_iter <- n.runs*l1*l2 #total number of iterations

#Main loop -----
for(n.rn in 1: n.runs){
  for(r in 1:l2){
    # Generate train and test data
    d1 <- data_gen(K,N_all,R,bet.dist[r],d,noise)
    X <- d1$X
    dat <- X[,names(X) %in% c('x1','x2','x3','x4','Y','idx'),with=F]
    setnames(dat, c("X1","X2","X3","X4","Y","idx"))
    dat[,X1:=scale(X1)]
    dat[,X2:=scale(X2)]
    dat[,X3:=scale(X3)]
    dat[,X4:=scale(X4)]
    
    # Randomly pick 6 observations from each group (idx) for testing
    tst_idx <- sample(1:(N_all[1]/R[1]),6)
    
    dat_tr = dat[,.SD[!tst_idx], by = idx] #Training
    tst = dat[,.SD[tst_idx], by = idx] #Test
    setcolorder(dat_tr,names(dat)) #Change the order of the columns
    setcolorder(tst,names(dat))
    
    
    X.tst <- tst[,c('X1','X2','X3','X4'),with=F]  #Test covariates
    Y.tst <- tst[,Y]  #Test responses
    idx.tst <- tst[,'idx',with=F]  #Test indecis
    
    for(n1 in alter_k){ #K from 0 to 8
      try({
        fin.res[t.r,K:=n1]  
        fin.res[t.r,Beta_dist := d1$dist]  #Delta beta
        cat(sprintf('Run %4d out of %d\n', t.r, total_num_iter))
        if(n1==0){  #Predict by mean
          yh <- predict_by_mean(dat_tr[,Y],Y.tst)
          y.err <- rmse(Y.tst,yh)
          fin.res[t.r,RMSE:=y.err] #RMSE error
          t.r <- t.r +1
        }else if(n1==1){
          yh <- single_lm(dat_tr,tst)
          y.err <- rmse(Y.tst,yh)
          fin.res[t.r,RMSE:=y.err]  #RMSE error
          t.r <- t.r +1
        }else{ #For K in 2:8
          # Applying the algorithm --------
          nr <- as.vector(table(dat_tr$idx))
          #Fit GMR using K components
          f1 <- fit_grp_mix_reg(dat_tr, K=n1, d=4, n.gr=Rtot, nr=nr, VERB = F)
          res1[,Clust:=apply(f1$tau,1,which.max)]  #Hard clustering
          
          #Save the true labels for each group
          th.cl <- X[,.(tru.label,idx),by=.(tru.label,idx)]
          # The columns are repeated (delete a pair)
          th.cl[,`:=`(tru.label=NULL,idx=NULL)]
          corr.cl.rate = vector()
          perm <- gtools::permutations(n1,n1,seq(1:n1))
          #Try all the permutations for cluster assignment and get the max of it
          for(ic in 1:(nrow(perm))){
            #swap the values for all the permutations
            a <- data.table(Clust = res1$Clust)
            a[,Clust := as.list(mapvalues(a$Clust,perm[(ic),], to = seq(1:n1)))]
            corr.cl.rate[ic] <- (sum(a$Clust==th.cl$tru.label))/Rtot
          }
          
          #which permutation achieved the max? 
          rm = which.max(corr.cl.rate)
          #Correct the labels 
          res1$Clust = mapvalues(res1$Clust,perm[(rm),], to = perm[1,])
          tau <- f1$tau[,perm[(rm),]]  #Correct the position of obtained taus
          #Correct the betas due to lable switching
          for(i in 1:n1){
            assign(paste0('beta',i),f1$beta[[perm[rm,i]]])
          }
          
          
          X.tst <- cbind(X.tst,idx.tst)
          yh <- list() #List that holds the yhat for each group
          for(i in 1:length(unique(X.tst$idx))){
            for(j in 1:n1){
              #y1 ,... , yk is Xnew * beta_rk *tau_rk
              assign(paste0('y',j),as.matrix(X.tst[idx==i,1:4]) %*% 
                       as.matrix(get(paste0('beta',j)))*tau[i,j])
            }
            #ytemp = sum from 1:K for y1 to yk obtained in previous loop
            ytemp <- rep(0,nrow(X.tst[idx==i]))
            for(k in 1:n1){
              ytemp <- ytemp+get(paste0('y',k))
            }
            # yh[[i]] is the MAP for Xnew that comes from group i
            yh[[i]] <- ytemp
          }
          
          yh <- unlist(yh)
          y.err <- rmse(Y.tst,yh) #RMSE error
          fin.res[t.r,RMSE:=y.err]
          t.r <- t.r +1
        }
        
      }, silent = T)
      
    }
  }
  
}

print(fin.res[,list(mean=mean(RMSE)),by=list(K,Beta_dist)])


# Plotting the result -------
library(latex2exp)
fin.res = na.omit(fin.res)
levels(fin.res$Beta_dist)[levels(fin.res$Beta_dist) == '0'] <- NA
fin.res$Beta_dist <- as.factor(fin.res$Beta_dist)
fin.res$K <- as.factor(fin.res$K)
nmi_dt2 <- fin.res[,lapply(.SD, mean, na.rm=TRUE), by=.(Beta_dist,K),.SDcols=c('RMSE')]
p1 <- ggplot(nmi_dt2,aes(K,RMSE,group=Beta_dist))+geom_point(size=4.5,shape=1,aes(color=Beta_dist))+geom_line(lty=2,aes(color=Beta_dist))+
  theme_bw()+labs(x='K', y='Average RMSE')+
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18))+
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))

p1
