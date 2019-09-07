library(dplyr)
# library(clue)
source('modules/GMR_data_gen2.R')
source('modules/fit_GMR.R')
source('modules/network_commons.R')

K <- 3
d <- 8
N <- rep(3e2, K)  # Number of observations for each component
R <- rep(10, K)   # Number of groups in each component
Rtot <- sum(R)    # Total number of groups

# compute_nmi <- function (label1, label2){
#   #require(clue)
#   cl_agreement(as.cl_hard_partition(label1), as.cl_hard_partition(label2), method = 'NMI') 
# }

# combinations used in simulation
total_num_runs <- 5
runs <- expand.grid(run_id=1:total_num_runs, bet_dist=c(2,4,6), noise_lev=seq(2,10,2))

# Run the simulations
for (r in 1:nrow(runs)) {
      cat(sprintf('Run %4d out of %d\n', r, nrow(runs))) 
      run <- runs[r,]
      out <- data_gen(K, N, R, run$bet_dist, d, run$noise_lev, VERB=F)
      dat <- out$data
      bets <- out$bets
      tru_label <- out$tru_label
      nr <- as.vector(table(dat$idx))
      fit <- fit_grp_mix_reg(dat, K=K, d=d, n.gr=Rtot, nr=nr, VERB=F)
      
      est_group_label <- apply(fit$tau, 1, which.max) # cluster assignment for groups
      est_label <- est_group_label[dat$idx] # cluster assignment of individual obs.
      # runs[r,"nmi"] <-  compute_nmi(tru_label, est_label)
      runs[r,"nmi"] <-  compute_mutual_info(tru_label, est_label) 
      runs[r,"n_iter"] <- fit$n.itr
}

# Calculate the averages
runs2 <- runs %>%
  mutate(bet_dist = factor(bet_dist)) %>%
  group_by(bet_dist, noise_lev) %>% 
  summarise(avg_nmi = mean(nmi, na.rm = T), avg_n_iter=mean(n_iter, na.rm = T)) 


# Plot the results
ggplot(runs2, aes(x=noise_lev, y=avg_nmi, color=bet_dist)) +
  geom_point(size=4.5,shape=1,aes(color=bet_dist))+geom_line(lty=2,aes(color=bet_dist)) +
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average NMI') +
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))

ggplot(runs2, aes(x=noise_lev, y=avg_n_iter, color=bet_dist)) +
  geom_point(size=4.5,shape=1,aes(color=bet_dist))+geom_line(lty=2,aes(color=bet_dist))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average Number of Iterations')+
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))+
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))+scale_y_continuous(trans='log2')#,limits = c(2,120))

