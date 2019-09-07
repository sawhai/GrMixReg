library(dplyr)
library(clue)
source('modules/GMR_data_gen2.R')
source('modules/fit_GMR.R')

K <- 2
d <- 2
N <- c(3e2, 3e2)  #Number of observations for each component
R <- c(10 ,10) #Number of groups in each component
Rtot <- sum(R)  #Total number of groups

compute_nmi <- function (label1, label2){
  #require(clue)
  cl_agreement(as.cl_hard_partition(label1), as.cl_hard_partition(label2), method = 'NMI') 
}

# combinations used in simulation
runs <- expand.grid(run_id=1:10, bet_dist=c(2.3,4.6,6.9), noise_lev=seq(2,10,2))

      
for (r in 1:nrow(runs)) {
      cat(sprintf('Run %4d out of %d\n', r, nrow(runs))) 
      run <- runs[r,]
      out <- data_gen(K, N, R, run$bet_dist, d, run$noise_lev)
      dat <- out$data
      bets <- out$bets
      tru_label <- out$tru_label
      nr <- as.vector(table(dat$idx))
      fit <- fit_grp_mix_reg(dat, K=2, d=2, n.gr=Rtot, nr=nr, VERB=F)
  
      est_group_label <- apply(fit$tau, 1, which.max) # cluster assignment for groups
      est_label <- est_group_label[dat$idx] # cluster assignment of individual obs.
      runs[r,"nmi"] <-  compute_nmi(tru_label, est_label) 

}

runs2 <- runs %>%
  mutate(bet_dist = factor(bet_dist)) %>%
  group_by(bet_dist, noise_lev) %>% 
  summarise(avg_nmi = mean(nmi, na.rm = T)) 

ggplot(runs2, aes(x=noise_lev, y=avg_nmi, color=bet_dist)) +
  geom_point(size=4.5,shape=1,aes(color=bet_dist))+geom_line(lty=2,aes(color=bet_dist)) +
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average NMI') +
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18)) +
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))
