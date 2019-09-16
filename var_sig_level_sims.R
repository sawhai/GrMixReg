library(data.table)
library(dplyr)
library(MASS)
library(flexmix)
library(ggplot2)
library(latex2exp)

d <- 2            # The dimension of data, i.e., number of covariates.
K <- 2            # The number of components. This code needs K <= d+1, since we generate equidistant components beta's
nobs <- 600       # Total number of observations
test_perc <- 0.2  # Percentage of data used for testing

total_num_runs <- 10  # number of replications. Reduce to speed up the simulation.

# combinations used in simulation
runs <- expand.grid(run_id=1:total_num_runs, 
                    bet_dist=c(4,8,12), 
                    noise_lev=seq(2,10,2),  # Noise level
                    G=10,                   # Number of groups in each component
                    nobs=nobs,
                    K=K,
                    d=d)

# Run the simulations
source("modules/sims.R")
set.seed(1234)
runs <- run_simulations(runs, test_perc = test_perc)

# Calculate the averages
runs2 <- runs %>%
  mutate(bet_dist = factor(bet_dist)) %>%
  group_by(bet_dist, noise_lev) %>% 
  summarize(avg_nmi = mean(nmi, na.rm = T), 
            avg_n_iter=mean(n_iter, na.rm = T),
            avg_beta_err = mean(beta_err, na.rm = T),
            avg_rmse = mean(rmse, na.rm = T),
            avg_rmse_lm = mean(rmse_lm, na.rm = T),
            avg_rmse_fmr = mean(rmse_fmr, na.rm = T))


xlab = '$\\sigma_k$'
custom_ggplot(runs2, aes(noise_lev, avg_nmi, color=bet_dist), title='Average NMI', xlab=xlab)
ggsave('nmi.pdf')
custom_ggplot(runs2, aes(noise_lev, avg_n_iter, color=bet_dist), title='Average Number of Iterations', xlab=xlab) + scale_y_continuous(trans='log10')
ggsave('niter.pdf')
custom_ggplot(runs2, aes(noise_lev, avg_beta_err, color=bet_dist), title='Average error $\\beta$', xlab=xlab) # + scale_y_continuous(trans='log10')
ggsave('beta_err.pdf')
custom_ggplot(runs2, aes(noise_lev, avg_rmse, color=bet_dist), title='Average RMSE (GMR)', xlab=xlab)
ggsave('gmr_rmse.pdf')
custom_ggplot(runs2, aes(noise_lev, avg_rmse_lm, color=bet_dist), title='Average RMSE (LM)', xlab=xlab)
ggsave('lm_rmse.pdf')
custom_ggplot(runs2, aes(noise_lev, avg_rmse_fmr, color=bet_dist), title='Average RMSE (FMR)', xlab=xlab)#+
ggsave('fmr_rmse.pdf')

