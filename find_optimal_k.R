library(data.table)
library(dplyr)
library(MASS)
library(flexmix)
library(ggplot2)
library(latex2exp)

rm(list=ls())     # clear the workspace

d <- 4            # The dimension of data, i.e., number of covariates.
K <- 4            # The number of (true) components. This code needs K <= d+1, since we generate equidistant components beta's
nobs <- 400       # Total number of observations
test_perc <- 0.2  # Percentage of data used for testing

total_num_runs <- 5   # number of replications. Reduce to speed up the simulation.

# combinations used in simulation
runs <- expand.grid(run_id=1:total_num_runs,
                    Ktru=K,   # true number of clusters
                    Ktst=0:8, # K's to test: Ktst=0 and Ktst=1 refer to prediction by mean and by linear regression, resp.
                    bet_dist=c(8,12), 
                    noise_lev= 6,  # Noise level
                    G=10,          # Number of groups in each component
                    nobs=nobs,
                    d=d)

# Run the simulations
source("modules/sims.R")
set.seed(1234)
runs <- find_optimal_k_cv(runs, test_perc = test_perc)

# Calculate the averages
runs2 <- runs %>%
   mutate(bet_dist = factor(bet_dist)) %>%
   group_by(bet_dist, Ktst) %>% 
   summarize(avg_rmse = mean(rmse, na.rm = T)) 

custom_ggplot(runs2, aes(Ktst, avg_rmse, color=bet_dist), title='Average RMSE', xlab='K')
ggsave('optimal_K.pdf')
