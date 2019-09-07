# Author: A. Amini 

label_vec2mat <- function(z, K=NULL) {
  # Converts a vector of labels from {1,...,K} to a (hard) cluster matrix
  #     each row of the cluster matrix is the one-hot encoding of the corresponding label
  #     if K is not given, the maximum of z is used as K
  if (is.null(K)) K <- max(z) 
  diag(K)[z,]
}

label_mat2vec <- function(Z) {
  apply(Z, 1, which.max)
}

compute_confusion_matrix <- function (z, y, K=NULL) {
  # Computes the confusion matrix between labels y and z
  # z and y are sets of labels:
  #     should be either vectors or n x k cluster matrices matrices
  #     the i-th row of a cluster matrix specifies the assignement of data point i to each of the k clusters 
  #     cluster matrices can be soft (each row has nonnegative values that sum to 1)
  # K   number of labels in both z and y. Automatically computer if unspecifed.
  
  if (is.null(K)) K = max(c(z,y))
  t(label_vec2mat(z,K)) %*% label_vec2mat(y,K)
}

compute_mutual_info  <- function(z,y) {
  # Computes the normalized mutual information between two clusterings (i.e., two sets of labels or two partitions)
  # z and y are sets of labels:
  #   should be either vectors or n x k cluster matrices matrices
  #   the i-th row of a cluster matrix specifies the assignement of data point i to each of the k clusters 
  #   cluster matrices can be soft (each row has nonnegative values that sum to 1)

  CM = compute_confusion_matrix(z,y)
  normCM = CM / sum(CM); # normalized confusion matrix
  IDX = CM == 0 # index of zero elements of CM so that we can avoid them

  jointEnt = - sum( (normCM[!IDX])*log(normCM[!IDX]) )
  indpt = matrix(rowSums(normCM),ncol=1) %*% matrix(colSums(normCM),nrow=1)
  muInfo = sum(normCM[!IDX] * log(normCM[!IDX] / indpt[!IDX]) )

  muInfo / jointEnt
}