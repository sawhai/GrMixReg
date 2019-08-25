# Fit GMR for K = 2, d = 2

rm(list=ls())
cat("\014")  #Clear the screen
library('MASS')
library('ggplot2')
library('caret')
library('miscTools')
library('plyr')
setwd("~/Documents/Gr.EM/Delivarable/")
source('GMR_data_gen.R')
source('fit_GMR.R')
source('NMI.AMI.calc.fn.R')
K = 2
N <- c(3e2, 3e2)  #Number of observations for each component
R <- c(10 ,10) #Number of groups in each component
Rtot <- sum(R)  #Total number of groups
cumsumR <- cumsum(c(0,R)) 
d <- 2
noise <- seq(2,10,2)
bet.dist <- c(2.3,4.6,6.9)
l1 <- length(noise); l2 <- length(bet.dist)
n.runs= 250
fin.res <- data.table(N1 =rep(N[1],(n.runs*l1*l2)),N2 =N[2],R1 = R[1], R2 =R[2],noise_level = 0,
                      r.bt.cl1.1 =0, r.bt.cl1.2 = 0, e.bt.cl1.1 = 0,e.bt.cl1.2=0,r.bt.cl2.1 =0,r.bt.cl2.2=0, 
                      e.bt.cl2.1 = 0,e.bt.cl2.2=0, corr.cl.rate = 0,bt.dist=0,num.itr = 0,nmi=0,b.err=0,RMSE = 0)

require('gtools')
perm <- permutations(K,K,seq(1:K))
t.r <- 1
cols <- c('x1','x2','Y')
res1 <- data.table(Groups= 1:Rtot, Clust = 0)
Conf <- list() #Holds the confusion matrices
for(n.rn in 1: n.runs){
  for(n in 1:l1){
    for(r in 1:l2){
      cat('Run    ', t.r,'\n', '---------------','\n')
      d1 <- data_gen(K,N,R,bet.dist[r],d,noise[n])
      fin.res[t.r,noise_level:=noise[n]]
      X <- d1$X
      dat <- X[,names(X) %in% c('x1','x2','Y','idx'),with=F]
      setnames(dat, c("X1","X2","Y","idx"))
      dat[,X1:=scale(X1)]
      dat[,X2:=scale(X2)]
      
      bets <- d1$bets
      fin.res[t.r,c('r.bt.cl1.1','r.bt.cl1.2'):=as.list(bets[,1])]
      fin.res[t.r,c('r.bt.cl2.1','r.bt.cl2.2'):=as.list(bets[,2])]
      fin.res[t.r,bt.dist := d1$dist]
      # Applying the algorithm --------
      nr <- as.vector(table(dat$idx))
      f1 <- fit_grp_mix_reg(dat, K=2, d=2, n.gr=Rtot, nr=nr)
      res1[,Clust:=apply(f1$tau,1,which.max)]
      
      cat('\n')
      #r1 <- data.table(Group = )
      #setnames(r1,c('Group','Clust','min.AIC'))
      fin.res[t.r,num.itr := f1$n.itr]
      #Save the true labels for each group
      th.cl <- X[,.(tru.label,idx),by=.(tru.label,idx)]
      # The columns are repeated (delete a pair)
      th.cl[,`:=`(tru.label=NULL,idx=NULL)]
      corr.cl.rate = vector()
      #Try all the permutations for cluster assignment and get the max of it
      for(ic in 1:(nrow(perm))){
        #swap the values for all the permutations
        a <- data.table(Clust = res1$Clust)
        a[,Clust := as.list(mapvalues(a$Clust,perm[(ic),], to = seq(1:K)))]
        corr.cl.rate[ic] <- (sum(a$Clust==th.cl$tru.label))/Rtot
      }
      
      #which permutation achieved the max? 
      rm = which.max(corr.cl.rate)
      res1$Clust = mapvalues(res1$Clust,perm[(rm),], to = perm[1,])
      fin.res$corr.cl.rate[t.r] <- corr.cl.rate[rm] #Maximum recovered cluster rate
      cl_agr <- f_rez(res1$Clust,th.cl$tru.label)
      fin.res$nmi[t.r] <- cl_agr[1]
      
      cm <- confusionMatrix(as.factor(res1$Clust),as.factor(th.cl$tru.label))
      cm <- cm$table
      Conf[[t.r]] <- cm
      #Correct the betas due to lable switching
      beta1 <- f1$beta[[rm]]
      beta2 <- f1$beta[[setdiff(c(1:2),rm)]]
      
      # u <- as.data.table(cbind(b2[,perm[rm,1],with=F],(b2[,perm[rm,2],with=F])))
      # setnames(u,c('beta1','beta2'))
      fin.res[t.r,c('e.bt.cl1.1','e.bt.cl1.2'):= as.list(beta1)]
      fin.res[t.r,c('e.bt.cl2.1','e.bt.cl2.2'):=as.list(beta2)]
      
      #Calculate the error for beta
      rr <- cm/Rtot
      br1 <- bets[,1]
      br2 <- bets[,2]
      d1 <- matrix(c(norm(br1-beta1,type='2'),norm(br1-beta2,type='2'),norm(br2-beta1,type='2'),norm(br2-beta2,type='2')),nrow=2)
      dr <- t(d1) %*% rr
      err <- sum(diag(dr))
      fin.res[t.r,b.err:= err]
      fin.res[t.r,b.err:=err]
      #Calculate the component weights
      for(i in 1:K){
        do.call('<-',list(paste0('N',i),sum(res1$Clust==i)))
      }
      for(i in 1:K){
        do.call('<-',list(paste0('pi',i),(get(paste0('N',i)))/Rtot))
      }
      tst <- data_gen(K,c(0.25*nrow(X),0.25*nrow(X)),R,bet.dist[r],d,noise[n])
      tst <- tst$X
      X.tst <- tst[,c('x1','x2'),with=F]
      Y.tst <- tst[,Y]
      for(i in 1:K){
        do.call('<-',list(paste0('y',i),get(paste0('pi',i)) * get(paste0('beta',i)) %*% t(X.tst)))
      }
      
      yh <- rep(0,length(Y.tst))
      for(i in 1:K){
        yh <- yh+get(paste0('y',i))
      }
      
      y.err <- rmse(Y.tst,yh)
      fin.res[t.r,RMSE:=y.err]
      t.r <- t.r +1
    }
  }
  
}

