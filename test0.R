library(dplyr)
# library(clue)
library(data.table)
library(ggplot2)
library(latex2exp)
source('modules/GMR_data_gen3.R')
source('modules/fit_GMR.R')
source('modules/network_commons.R')

K <- 3
d <- 10
nobs <- 600       # Total number of observations
G <- 10           # Number of groups in each component
Rtot <- K*G       # Total number of groups
# N <- rep(3e2, K)  # Number of observations for each component
#R <- rep(10, K)  # Number of groups in each component
#Rtot <- sum(R)   # Total number of groups

meanColSSQ <- function(X,Y) mean(colSums((X-Y)^2))
# compute_nmi <- function (label1, label2){
#   #require(clue)
#   cl_agreement(as.cl_hard_partition(label1), as.cl_hard_partition(label2), method = 'NMI') 
# }

# combinations used in simulation
total_num_runs <- 5
runs <- expand.grid(run_id=1:total_num_runs, bet_dist=c(2,4,6), noise_lev=seq(2,10,2))

# Run the simulations
dt <- system.time(
for (r in 1:nrow(runs)) {
      cat(sprintf('Run %4d out of %d\n', r, nrow(runs))) 
      run <- runs[r,]
      #out <- data_gen(K, N, R, run$bet_dist, d, run$noise_lev, VERB=F)
      out <- data_gen(K, nobs, G, run$bet_dist, d, run$noise_lev, VERB=F)
      dat <- out$data
      tru_bets <- out$bets
      tru_obs_labels <- out$tru_labels #  true obs. labels :note the change to label"s" from label
      tru_grp_labels <- out$grp_labels # true group labels
      nr <- tabulate(dat$idx) #as.vector(table(dat$idx))
      fit <- fit_grp_mix_reg(dat, K=K, d=d, n.gr=Rtot, nr=nr, VERB=F)
      
      est_grp_labels <- label_mat2vec(fit$tau) # cluster assignment for groups
      est_obs_labels <- est_grp_labels[dat$idx] # cluster assignment for individual obs.
      est_bets <- do.call(cbind,fit$beta) # estimated beta's
      
      runs[r,"nmi"] <-  compute_mutual_info(tru_obs_labels, est_obs_labels)
      # runs[r,"nmi"] <-  compute_mutual_info(tru_grp_labels, est_grp_labels) 
      runs[r,"n_iter"] <- fit$n.itr
      runs[r,"beta_err"] <- meanColSSQ(tru_bets[ , tru_grp_labels], est_bets[ , est_grp_labels])
}
)["elapsed"]

cat("Total runtime = ", dt, "(s)\n")
# Calculate the averages
runs2 <- runs %>%
  mutate(bet_dist = factor(bet_dist)) %>%
  group_by(bet_dist, noise_lev) %>% 
  summarise(avg_nmi = mean(nmi, na.rm = T), 
            avg_n_iter=mean(n_iter, na.rm = T),
            avg_beta_err = mean(beta_err, na.rm = T)) 


# Plot the results
# p <- ggplot(runs2, aes(x=noise_lev, y=avg_nmi, color=bet_dist)) +
#   geom_point(size=4.5,shape=1,aes(color=bet_dist))+geom_line(lty=2,aes(color=bet_dist)) +
#   theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average NMI') +
#   theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=18),
#         axis.text.y = element_text(size=18)) +
#   scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))
# 
# print(p)
# 
# ggplot(runs2, aes(x=noise_lev, y=avg_n_iter, color=bet_dist)) +
#   geom_point(size=4.5,shape=1,aes(color=bet_dist))+geom_line(lty=2,aes(color=bet_dist))+
#   theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average Number of Iterations')+
#   theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=15),
#         axis.text.y = element_text(size=15))+
#   scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))+scale_y_continuous(trans='log2')#,limits = c(2,120))

custom_ggplot <- function( data, aesth_map, title) {
  ggplot(data, aesth_map) + 
    geom_point(size=4.5,shape=1) +
    geom_line(lty=2) +
    theme_bw() + 
    labs(x=TeX('$\\sigma_k$'), y=TeX(title))+
    theme(text = element_text(size=25),
          panel.grid.major=element_line(colour='gray75'),
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15))+
    scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))#+scale_y_continuous(limits = c(0,1.5))
}
custom_ggplot(runs2, aes(noise_lev, avg_nmi, color=bet_dist), 'Average NMI' )
ggsave('nmi.pdf')
custom_ggplot(runs2, aes(noise_lev, avg_n_iter, color=bet_dist), 'Average Number of Iterations' )
ggsave('niter.pdf')
custom_ggplot(runs2, aes(noise_lev, avg_beta_err, color=bet_dist), 'Average error $\\beta$' )
ggsave('beta_err')
