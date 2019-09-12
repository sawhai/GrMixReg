#Functions that are called for generating data
library(MASS)
#To generate response (called by data_gen)
grp_mix_reg_gen_response <- function(dat, d, bets, clust, noise_sig) {
  # Generate response. More general. Changes "dat" in place (pass by reference)
  K = clust
  n_k <- nrow(dat)
  # Mu_k <- as.matrix(dat[, 1:d, with=F]) %*% bets[K,]
  Mu_k <- as.matrix(dat[, 1:d, with=F]) %*% bets[ , K]  
  noise_k <- noise_sig*rnorm(n_k)
  y <- Mu_k + noise_k
  # dat[, Mu := Mu_k]
  # dat[, noise := noise_k]
  
  dat[, Y := y]
  
}

# #Data preparation ------
# # K   := number of components
# # # R := vector containing Number of groups in each component
# # # r := value to control the distance between betas
# # # noise_level := noise signal when generating y
# # # d := dimension (number) of covariates
# data_gen <- function(K, N, R, r, d, noise_level, normalize=T, VERB=T){
#   
#   G <- R[1]
#   nr <- N[1]/R[1]
#   # Rtot <- sum(R)  #Total number of groups
#   #cumsumR <- cumsum(c(0,R))

# --- Simulated data generation ---
# K   := number of components
# # G := number of groups in each component
# # nobs := total number of observations 
# # bet_dist := pairwise distance between betas (scalar)
# # noise_level := noise signal when generating y
# # d := dimension (number) of covariates
data_gen <- function(K, nobs, G, bet_dist, d, noise_level, normalize=T, VERB=T){
  
  R <- K*G # total number of groups
  if ( nobs %% R != 0) warning('nobs is not divisible by K*G ... rounding.')
  nobs_per_grp <- round(nobs / R)   # is nobs_per_grp = nr? if yes, do we need both loines?
  nr <- nobs / R # number of observations per group
  nobs_per_clust <- nobs / K
  
  mu <- matrix(0,ncol = d, nrow = K)  #Matrix to hold the means
  sigma <- lapply(1:K, function(k) gen_scaled_Wishart(50, d))
  B <- generate_equidistant_pts(d, K, bet_dist, VERB=VERB) #Get the points
  
  X <- data.table()
  for (k in 1:K) {
    #t <- data.table(mvrnorm(N[k],mu[k,],sigma[[k]])) # the covariates
    t <- data.table(mvrnorm(nobs_per_clust, mu[k,], sigma[[k]])) # the covariates
    if (normalize) t = t[, lapply(.SD, scale)]
    
    #idx <- rep(k,nrow(t))
    grp_mix_reg_gen_response(t,d,B,k,noise_level)
    t[,tru.label := k]
    grp_idx  <- (k-1)*G + rep(1:G, nr)
    # grp_idx2 <- rep((cumsumR[k]+1):cumsumR[k+1],N[1]/R[1])
    # print(sum(grp_idx - grp_idx2))
    t[,idx := grp_idx]
    
    X <- rbind(X,t)
    
  }
  grp_labels <- rep(1:K, each=G)
 
  tru_labels <- X[, tru.label]
  X[, tru.label := NULL]
  names(X)[1:K] <- sapply(1:K, function(i) paste0('x',i))
  return(list(bets = B, data = X, tru_labels = tru_labels, grp_labels=grp_labels))
  
}

#Function to generate covariance matrix using a normalized Wishart distribution
gen_scaled_Wishart <- function(N, p, Sigma =  diag(p)) {
  X <- mvrnorm(N, rep(0,p), Sigma)
  t(X) %*% X / N
}


rotation = function(x,y){
  u=x/sqrt(sum(x^2))
  
  v=y-sum(u*y)*u
  v=v/sqrt(sum(v^2))
  
  cost=sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2))
  
  sint=sqrt(1-cost^2);
  
  diag(length(x)) - u %*% t(u) - v %*% t(v) + 
    cbind(u,v) %*% matrix(c(cost,-sint,sint,cost), 2) %*% t(cbind(u,v))
}

scalar1 <- function(x) {x / sqrt(sum(x^2))}
##################################################################
# Funtion to generate equidistance points used to generate betas
# dim      := the dimension of the betas
# num.pts  := number of betas that we wish to generate
# rd       := the disired distance between the points
# VERB     := (True/Fale) Print diagnostics
generate_equidistant_pts <- function(dim, num.pts, rd, VERB=T){
  
  dns = 2   #initial distance between pair of points
  pns <- matrix(rep(0,dim*(dim+1)),ncol=(dim+1))   #matrix to hold the points on its columns
  pns[,1] <- c(1,rep(0,(dim-1))); pns[,2] <- c(-1,rep(0,(dim-1))) #initialize the first two points
  pn <- matrix(rep(0,dim*(dim+1)),ncol=(dim+1))  
  for(i in 2:dim){
    dis <- dns^2-(norm(pns[-i,1],type='2')^2) #dist: dns^2-two points except for the ith entry
    pns[i,i+1] <- sqrt(dis)-pns[i,i] #calculate the ith entry of the i+1th point
    p.pt <- apply(pns[,1:(i+1)],1,function(x) sum(x)/(i+1)) #pivot point of the calculated points
    pn[,1:(i+1)] <- apply(pns[,1:(i+1)],2,function(x) x-p.pt) #Shift to the origin
    rnz <- norm(pn[,1],type='2')  #unscaled radius of the hypersphere
    K <- rd/rnz  #Ratio to scale
    pns <- K*pn  #Scale to the desired radius
    rn <- norm(pns[,1],type='2')  #Scaled readius of the hypersphere
    dns <- norm(pns[,1]-pns[,2],type='2') #Distant between the scaled points
  }
  #Rotate
  #create two normalized random points to calculate the rotation matrix
  a1 <- scalar1(runif(dim,-2,2)); a2 <- scalar1(runif(dim,-2,2))  
  
  #Calculate the rotation matrix
  Rot.mat <- rotation(a1,a2)
  p <- Rot.mat %*% pns #Rotate the points
  
  #randomly pick whatever number of points you desire
  u <- sample(1:(dim+1),num.pts)
  p <- p[,u]
  
  # center the points and rescale to have the desired distance
  curr_dist <- sqrt(sum((p[,1] - p[,2])^2))
  p <- (rd/curr_dist) * sweep(p,1, rowMeans(p),"-")
  
  # print diagnostic: see if the points are really equi-distant
  if (VERB) {
    pwdist <- dist(t(p)) # pairwise distances of columns of B
    pwd_dist_cv <- sd(pwdist)/mean(pwdist)
    cat("Coeff of Var. of pairwaise distnces of points =", ifelse(is.na(pwd_dist_cv), 0, pwd_dist_cv),"\n")
    cat("dist(p[,1], p[,2]) =", pwdist[1],"\n")
  }  
    
  p
  
}

# Function to predict a hold out set using GMR
# k: number of clusters
# GMR_fit: a GMR object resulted from fitting GMR 
# test_set: a hold out set in the form generated from data_gen routine

predict_GMR <- function(k, GMR_fit, test_set){
  #Exctract betas
  for(i in 1:k){
    assign(paste0('beta',i),GMR_fit$beta[[i]])
  }
  yh <- list() #List that holds the yhat for each group
  for(i in 1:length(unique(test_set$idx))){
    for(j in 1:k){
      #y1 ,... , yk is Xnew * beta_rk *tau_rk
      assign(paste0('y',j),as.matrix(test_set[idx==i,!c('Y','idx'),with=F]) %*% 
               as.matrix(get(paste0('beta',j)))*GMR_fit$tau[i,j])
    }
    #ytemp = sum from 1:K for y1 to yk obtained in previous loop
    ytemp <- rep(0,nrow(test_set[idx==i]))
    for(h in 1:K){
      ytemp <- ytemp+get(paste0('y',h))
    }
    # yh[[i]] is the MAP for Xnew that comes from group i
    yh[[i]] <- ytemp
  }
  yh <- unlist(yh)
  yh
}



