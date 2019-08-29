# Function to fit GMR
fit_grp_mix_reg <- function(dat, K, d, n.gr, nr, max.itr=1000, pi=NA, tau_init=NA, tol=1e-6, VERB=T) {
  # Assuming columns 1:d of dat are X and column d+1 is Y
  # dat should be a data.table with columns labeld X1 X2 ... Y idx
  
  # n.gr is the same as Rtot
  # nr: vector of group sizes or a scalar if all are the same
  # d: dimension of the covariates
  
  if (is.na(pi)){
    pi = rep(1.,K)/K
  }
  
  # let us make "nr" a vector if it is not already one
  if (length(nr) != n.gr) {
    nr <- rep(nr, n.gr)
  }
  
  if (is.na(tau_init)){
    temp <- diag(K)
    tau_init <- temp[sample(K, size=n.gr, replace=T),]
    
    #tau_init <- matrix(1/K, nrow=n.gr, ncol=K)  # does not work! This is a fixed point of iterations.
  }
  
  if (VERB) cat('Computing per group covariances ... ')  # This is the time-consuming loop
  Sigh <- list()
  rhoh <- list()
  for (r in 1:n.gr) {
    X <- as.matrix( dat[idx==r, 1:d, with=F] )
    #y <- as.vector( dat[idx==r, d+1, with=F] )  # this does not give a vector
    y <- dat[idx == r,][[d+1]]
    
    Sigh[[r]] <- t(X) %*% X / nr[r]
    rhoh[[r]] <- as.vector(t(X) %*% y) /nr[r]
  }
  if (VERB) cat('Done. \n')
  
  # http://stats.stackexchange.com/questions/8605/column-wise-matrix-normalization-in-r
  row_normalize <- function(x) sweep(x, 1, rowSums(x), FUN="/")  # Can probably use matrix multiplication
  col_normalize <- function(x) sweep(x, 2, colSums(x), FUN="/")  # alternative t(t(x) / colSums(x))
  
  # http://stats.stackexchange.com/questions/126602/adding-a-value-to-each-element-of-a-column-in-r
  # http://stackoverflow.com/questions/24520720/subtract-a-constant-vector-from-each-row-in-a-matrix-in-r
  add_vec_to_each_row <- function(vec,mat) t(vec + t(mat))
  
  safe_exp <- function(X) exp(-rowMaxs(X)+X)
  
  #Initialize tau for each group
  tau <- tau_init
  n.itr <- 1
  if (VERB) cat(sprintf('%3s | %-7s \n', 'itr', 'tau_err'))
  for(itr in 1:max.itr){
    # for (r in 1:n.gr) {
    #   w[r,] <- nr[r] * tau[r,]
    # }
    #w <- diag(nr) %*% tau
    w <- nr * tau  # nr is a vector, this multiplies tau[r,] by nr[r], hopefully row-wise multiplication
    
    #notmalized w (w.check)
    #w.check <-sweep(w, 2, colSums(w), FUN="/")
    w.check <- col_normalize(w)
    
    Erk <- matrix(0, nrow=n.gr, ncol=K)
    Sigt <- list()
    rhot <- list()
    betah <- list()
    for (k in 1:K) {
      # http://stackoverflow.com/questions/34820827/r-how-to-multiply-list-elements-by-a-vector
      # http://stackoverflow.com/questions/11641701/sum-a-list-of-matrices
      Sigt[[k]] <- Reduce('+',Map('*', Sigh, w.check[,k]))  
      rhot[[k]] <- Reduce('+',Map('*', rhoh, w.check[,k]))
      
      betah[[k]] <- solve(Sigt[[k]], rhot[[k]]) 
      pi[k] <- mean(tau[,k])
      
      #Calculate Erk
      X <- as.matrix(dat[,1:d, with=F])
      pred <- X %*% betah[[k]]
      dat[, rs := (Y - pred)^2] # dat should have a column labeled "Y"
      Erk[,k] <- dat[, mean(rs), by=idx]$V1
      n.itr <- n.itr+1
    }
    
    sig2 <- diag(t(w.check) %*% Erk)
    
    #What is this fore? 
    gamma.r.k <- safe_exp(-0.5 * diag(nr) %*% add_vec_to_each_row( log(sig2), Erk %*% diag(1/sig2) ) ) %*% diag(pi)
    
    tau_old <- tau
    tau <- row_normalize(gamma.r.k)
    
    #delta <- norm(tau-tau_old) / max(norm(tau),1)
    delta <- norm(tau-tau_old, type="i") # maximum absolute row sum
    if (VERB) cat(sprintf('%3d | %5.2e \n', itr, delta))
    if (delta < tol) break 
    
  }
  
  return( list(tau=tau, beta=betah, sig2=sig2, pi=pi, Erk=Erk, tau_init=tau_init, n.itr = n.itr) )
}



ones <- function(m,n) matrix(rep(1,m*n),nrow=m)
zeros <- function(m,n) matrix(rep(0,m*n),nrow=m)
