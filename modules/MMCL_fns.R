# Functions to apply MMCL++ to data


lm.fn <- function(df,d,cols){ #given data and dimension(d), fit a linear regression
  # cols is: "x columns and Y" in df
  fit <- lm(Y ~ .,data=df[, ..cols])
  fit
}

# To smartly pick the initial groups (for MMCL++)
pick_k_init_groups <- function(data,k,R,d,cols){
  # cols is: "x columns and Y" in data
  a <- sample(R,1) #Pick the first group at random
  init.grps<- rep(0,k) #to hold the final groups
  init.grps[1] <- a # The first one
  grs <- setdiff(1:R,a)
  c1 <- data[idx %in% a, ..cols]
  models <- list()
  models[[1]] <- lm(Y~.,c1)
  h = 2
  # to find the group that is furthest from the first group
  while(h <= k){
    j <- length(models)
    rmses <- data.frame(group = grs)
    
    col_in <- paste0('rmse.model',1:j)

    rmses[,c(col_in)] <- 0
    for(i in 1:length(grs)){
      df <- data[idx %in% i]
      for(v in 1:j){
        yhat <- predict(models[[v]],df)
        rmses[i,col_in[v]] <- rmse(yhat,df[,Y])
      }
    }
    
    if(h!=2){
      rmses$sum.rmse <- apply(rmses[,-1],1,sum)
      max.rmse <- which.max(rmses$sum.rmse)
      init.grps[h] <- rmses[max.rmse,'group']
      grs <- setdiff(grs,init.grps[h])
      
    }
    else{
      max.rmse <- which.max(rmses$rmse.model1)
      init.grps[h] <- rmses[max.rmse,'group']
      grs <- setdiff(grs,init.grps[h])
    }
    c2 <- data[idx %in% init.grps[h], ..cols]
    models[[h]] <- lm(Y~.,c2)
    h <- h+1
  }
  
  return(init.grps)
}


best.model <- function(df,models,cols){ # To find the model that best predicts the group in step 1
  min.rmse <- 1e20
  n <- length(models)
  rmses <- rep(0,n)
  for(i in 1:n){
    yhat <- predict(models[[i]],df[, ..cols])
    rmses[i] <- rmse(yhat,df$Y)
  }
  min.rmse <- min(rmses)
  Best.rmse.Model <- which.min(rmses)
  return(c(Best.rmse.Model,min.rmse))
}

best.model.2 <- function(df,models,res,cols){ # To find the model that best predicts the group in step 2
  min.rmse <- 1e20
  n <- length(models)
  mis.clus <- setdiff(seq(1:n),res$Best.rmse.Model)
  #Exclude the empty cluster (if any)
  clusts <- setdiff(seq(1:n),mis.clus)
  rmses <- rep(0,n)
  for(i in clusts){
    yhat <- predict(models[[i]],df[, ..cols])
    rmses[i] <- rmse(yhat,df$Y)
    
  }
  min.rmse <- min(rmses)
  Best.rmse.Model <- which.min(rmses)
  return(c(Best.rmse.Model,min.rmse))
}

#Generate K initial models
k_initial_lm_models <- function(data,grps,k,d,cols){
  models <- list()
  for(i in 1:k){
    models[[i]] <-lm.fn(data[data$idx == grps[i]],d,cols)
  }
  models
}

# Main function to perform MMCL
fit_mmcl <- function(data,k,R,d,cols,num.itr=50,tol=0.0001,verb=F){
  #MMCL------
  #######################################################################
  res1 <- data.frame(Groups= 1:R,Best.rmse.Model = NA, min.rmse=0)
  current.rmse <- 1e20
  new.rmse <- current.rmse - 1e19
  step1 <- T
  num.itr <- 10
  betah <- list()
  t <- 1
  init.grps <- pick_k_init_groups(data,k,R,d,cols)
  while(abs((new.rmse- current.rmse)/current.rmse) > tol && t < num.itr){
    #inside while loop
    current.rmse <- new.rmse
    if(step1){
      models <- k_initial_lm_models(data,init.grps,k,d,cols)
      for(j in 1:R){
        res1[j,2:3] <- best.model(data[idx == j],models,cols)
      }
      step1 <- F
    } else {
      grp_members <- lapply(1:k, function(x) res1$Groups[res1$Best.rmse.Model==x]) #Get members of each group
      empty_groups <- which(sapply(grp_members,length)==0) #Which groups have no members?
      if(length(empty_groups)!=0){ #If there are empty groups, assign one group to each component
        for(i in 1:length(empty_groups)){
          res1$Best.rmse.Model[i] <- empty_groups[i]
        }
        grp_members <- lapply(1:k, function(x) res1$Groups[res1$Best.rmse.Model==x]) #Get members of each group
      }
      
      # Get the data for each cluster (consisted of observations of groups assigned to each cluster)
      dat_grps <- lapply(grp_members, function(x) data[idx %in% x])
      
      models2 <- list() #To hold the models 
      for(i in 1:k){
        models2[[i]] <- lm.fn(dat_grps[[i]],d,cols)
      }
      
      for(j in 1:R){
        res1[j,2:3] <- best.model.2(data[idx == j],models,res1,cols)
      }
      
    }
    new.rmse <- sum(res1$min.rmse)
    if(verb){
      cat('itr    ',t,'\n')
      cat('New rmse =    ', new.rmse, '\n')
    }

    t <- t+1
  }
  # End MMCL -------
  ####################################################
  betah <- lapply(models2, function(x) coefficients(x)[-1]) #Extract estimated betas (MMCL)
  return(list(result=res1,betah=betah))
  
}

# Plot

mmcl_ggplot <- function( data, aesth_map, title, xlab,labs){
  ggplot(data, aesth_map) + 
    geom_point(size=4.5,shape=1) +
    geom_line(lty=2) +
    theme_bw() + 
    labs(x=TeX(xlab), y=TeX(title))+
    theme(text = element_text(size=25),
          legend.title=element_blank(), legend.position="bottom",
          panel.grid.major=element_line(colour='gray75'),
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15))+
    scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3),labels = labs)
}
