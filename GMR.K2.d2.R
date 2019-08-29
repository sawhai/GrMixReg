# Fit GMR for K = 2, d = 2 and plot the result

rm(list=ls())
library(ggplot2)
library(data.table)
library(MASS) # mvnnorm
library(plyr) # mapvalues
library(caret) # confusionMatrix 
require(gtools)

rowMaxs <- function(X) apply(X, 1, function(row) max(row))
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))
# library(Metrics) # rmse

source('modules/GMR_data_gen.R')
source('modules/fit_GMR.R')
source('modules/NMI.AMI.calc.fn.R')


K <- 2
d <- 2

N <- c(3e2, 3e2)  #Number of observations for each component
R <- c(10 ,10) #Number of groups in each component
Rtot <- sum(R)  #Total number of groups
cumsumR <- cumsum(c(0,R)) 

noise <- seq(2,10,2)
bet.dist <- c(2.3,4.6,6.9)
l1 <- length(noise); l2 <- length(bet.dist)
n.runs <- 10
fin.res <- data.table(N1 =rep(N[1],(n.runs*l1*l2)),N2 =N[2],R1 = R[1], R2 =R[2],Noise = 0,
                      r.bt.cl1.1 =0, r.bt.cl1.2 = 0, e.bt.cl1.1 = 0,e.bt.cl1.2=0,r.bt.cl2.1 =0,r.bt.cl2.2=0, 
                      e.bt.cl2.1 = 0,e.bt.cl2.2=0, corr.cl.rate = 0,Beta_dist=0,num.itr = 0,NMI=0,Beta_err=0,RMSE = 0)


perm <- gtools::permutations(K,K,seq(1:K))
t.r <- 1
res1 <- data.table(Groups= 1:Rtot, Clust = 0)
Conf <- list() #Holds the confusion matrices
total_num_iter <- n.runs*l1*l2
for(n.rn in 1: n.runs){
  for(n in 1:l1){
    for(r in 1:l2){
      cat(sprintf('Run %4d out of %d\n', t.r, total_num_iter)) 
      
      d1 <- data_gen(K,N,R,bet.dist[r],d,noise[n])
      fin.res[t.r,Noise:=noise[n]]
      X <- d1$X
      dat <- X[,names(X) %in% c('x1','x2','Y','idx'),with=F]
      setnames(dat, c("X1","X2","Y","idx"))
      dat[,X1:=scale(X1)]
      dat[,X2:=scale(X2)]
      
      bets <- d1$bets
      fin.res[t.r,c('r.bt.cl1.1','r.bt.cl1.2'):=as.list(bets[,1])]
      fin.res[t.r,c('r.bt.cl2.1','r.bt.cl2.2'):=as.list(bets[,2])]
      fin.res[t.r,Beta_dist := d1$dist]
      
      
      # Applying the algorithm --------
      nr <- as.vector(table(dat$idx))
      f1 <- fit_grp_mix_reg(dat, K=2, d=2, n.gr=Rtot, nr=nr, VERB=F)
      res1[,Clust:=apply(f1$tau,1,which.max)]
      fin.res[t.r, num.itr := f1$n.itr]
      #Save the true labels for each group
      th.cl <- X[ , .(tru.label,idx), by=.(tru.label,idx)]
      # The columns are repeated (delete a pair)
      th.cl[,`:=`(tru.label=NULL, idx=NULL)]
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
      fin.res$NMI[t.r] <- cl_agr[1]
      
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
      fin.res[t.r,Beta_err:= err]
      fin.res[t.r,Beta_err:=err]
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

##############################################################################################################

# Plotting the result
library(latex2exp)
fin.res <- fin.res[complete.cases(fin.res), ]
#levels(fin.res$Beta_dist)[levels(fin.res$Beta_dist) == '0'] <- NA
fin.res$Noise <- as.factor(fin.res$Noise)
fin.res$Beta_dist <- as.factor(fin.res$Beta_dist)
# fin.res <- fin.res[Noise %in% c(2,4,6,8,10),]
nmi_dt2 <- fin.res[,lapply(.SD, mean, na.rm=TRUE), by=.(Noise,Beta_dist),.SDcols=c("NMI", "Beta_err",'RMSE','num.itr')]
# levels(nmi_dt2$Beta_dist) <- c('4','8','12')

p1 <- ggplot(nmi_dt2,aes(Noise,NMI,group=Beta_dist))+geom_point(size=4.5,shape=1,aes(color=Beta_dist))+geom_line(lty=2,aes(color=Beta_dist))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average NMI')+
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18))+
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))#+scale_y_continuous(limits = c(0.2,1))

p1

p2 <- ggplot(nmi_dt2,aes(Noise,Beta_err,group=Beta_dist))+geom_point(size=4.5,shape=1,aes(color=Beta_dist))+geom_line(lty=2,aes(color=Beta_dist))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y=TeX('Average error $\\beta$'))+
  theme(text = element_text(size=25),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))+
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))#+scale_y_continuous(limits = c(0,1.5))
p2
p3 <- ggplot(nmi_dt2,aes(Noise,num.itr,group=Beta_dist))+geom_point(size=4.5,shape=1,aes(color=Beta_dist))+geom_line(lty=2,aes(color=Beta_dist))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average Number of Iterations')+
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))+
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))+scale_y_continuous(trans='log2',limits = c(2,120))
p3
p4 <- ggplot(nmi_dt2,aes(Noise,RMSE,group=Beta_dist))+geom_point(size=4.5,shape=1,aes(color=Beta_dist))+geom_line(lty=2,aes(color=Beta_dist))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average RMSE')+
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18))+
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))

p4


