#Functions that are called for generating data
#To generate response (called by data_gen)
grp_mix_reg_gen_response <- function(dat, d, bets, clust, noise_sig) {
  # Generate response. More general. Changes "dat" in place (pass by reference)
  K = clust
  n_k <- nrow(dat)
  Mu_k <- as.matrix(dat[, 1:d, with=F]) %*% bets[K,]  
  noise_k <- noise_sig*rnorm(n_k)
  
  dat[, Mu := Mu_k]
  dat[, noise := noise_k]
  
  dat[, Y := Mu + noise]
  
}

#Data preparation ------
# K:= #of components
# # R := vector containing Number of groups in each component
# # r := value to control the distance between betas
# # noise_level := noise signal when generating y
# # d := dimension (number) of covariates
data_gen <- function(K,N,R,r,d,noise_level){
  
  Rtot <- sum(R)  #Total number of groups
  cumsumR <- cumsum(c(0,R))
  mu <- matrix(0,ncol = d, nrow = K)  #Matrix to hold the means
  sig2 <- abs(rnorm(K,1,0))
  for(i in 1:K){
    do.call('<-',list(paste0('rn',i),abs(rnorm(1,0,sig2[i]))))
  }
  
  for(i in 1:K){
    do.call('<-',list(paste0('m',i),abs(rnorm(K,0,sig2[i]))))
  }
  
  sigma <- list()
  for(i in 1:K){
    sigma[[i]] <- gen_scaled_Wishart(50, d)
  }
  B <- generate_equidistant_pts(d,K,r) #Get the points
  #Normalize them
  dist <- round(norm(B[,1]-B[,2],type='2'))
  #Generate the grouped data
  X <- list()
  #index <- list()
  X <- data.table()
  for (k in 1:K) {
    t <- data.table(mvrnorm(N[k],mu[k,],sigma[[k]]))
    
    idx <- rep(k,nrow(t))
    grp_mix_reg_gen_response(t,d,B,k,noise_level)
    t[,tru.label :=k]
    t[,idx := rep((cumsumR[k]+1):cumsumR[k+1],N[1]/R[1])]
    X <- rbind(X,t)
    
  }
  names <- rep('x',K)
  for(i in 1:K){
    names[i] <- paste0('x',i)
  }
  names(X)[1:K] <- names
  return(list(bets = B, X = X,dist = dist))
  
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
#Funtion to generate equidistance points used to generate betas
#dim := the dimension of the betas
#num.pts := number of betas that we wish to generate
#rd := the disired distance between the points
generate_equidistant_pts <- function(dim,num.pts,rd){
  #desired distance 
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
  p
  
}

