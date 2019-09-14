library(dplyr)
# library(clue)
library(data.table)
library(ggplot2)
library(latex2exp)
rm(list=ls())
d <- 4            # The dimension of data, i.e., number of covariates.
K <- 4
# The number of components. This code needs K <= d+1, since we generate equidistant components beta's
nobs <- 600       # Total number of observations
test_perc <- 0.2  # Percentage of data used for testing

total_num_runs <- 50   # number of replications. Reduce to speed up the simulation.

# combinations used in simulation
runs <- expand.grid(K=rep(0:8,total_num_runs),
                    bet_dist=c(8,12), 
                    noise_lev= 6,  # Noise level
                    G=10,                   # Number of groups in each component
                    nobs=nobs,
                    d=d)

# Run the simulations
source("Sim_k.R")
runs <- optimal_k_cv(runs, tru_k = K,test_perc = test_perc)

# Calculate the averages
runs2 <- runs %>%
   mutate(bet_dist = factor(bet_dist)) %>%
   group_by(bet_dist, K) %>% 
   summarize(avg_rmse = mean(rmse, na.rm = T)) 

custom_ggplot(runs2, aes(K, avg_rmse, color=bet_dist), title='Average RMSE', xlab='K')
ggsave('optimal_K.pdf')
