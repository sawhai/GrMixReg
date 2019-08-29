# This script produces the plot for comparing GMR and MMCL++
rm(list=ls())
library(data.table)
library(ggplot2)
library(latex2exp)
# setwd("/Users/ha/Documents/Gr.EM/")
res1 <- get(load('data/res2_mmcl_k2_d2_n100.RData'))
res2 <- get(load('data/res2_gmr_k2_d2_n100.RData'))
res.combined <- res1[,c('Noise','Beta_dist','NMI','RMSE')]
res.combined2 <- res2[,c('Noise','Beta_dist','NMI','RMSE')]
res.combined$Noise <- as.factor(res.combined$Noise)
res.combined$Beta_dist <- as.factor(res.combined$Beta_dist)
res.combined2$Noise <- as.factor(res.combined2$Noise)
res.combined2$Beta_dist <- as.factor(res.combined2$Beta_dist)
nmi_dtm <- res.combined[,lapply(.SD, mean, na.rm=TRUE), by=.(Noise,Beta_dist),.SDcols=c("NMI",'RMSE')]
nmi_dtg <- res.combined2[,lapply(.SD, mean, na.rm=TRUE), by=.(Noise,Beta_dist),.SDcols=c("NMI",'RMSE')]
nmi_dtm [,Algo:='MMCL++']
nmi_dtg[,Algo:='GMR']
nmi_dt2 <- rbind(nmi_dtm,nmi_dtg)
p1 <- ggplot(nmi_dt2[Beta_dist==4&Noise %in% c(2,4,6,8,10)],aes(Noise,NMI,group=Algo)) + 
  geom_point(size=4.5,shape=1,aes(color=Algo))+
  geom_line(lty=2,aes(color=Algo))+
    theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average NMI')+
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.position="bottom")+ scale_color_manual(name=TeX(' '),values=c('blue3','darkorchid4'))+
  scale_y_continuous(limits = c(0,1))
p1
ggsave(p1,file='nmi_compare_n100_bet4.pdf')

p2 <- ggplot(nmi_dt2[Beta_dist==12 & Noise %in% c(2,4,6,8,10)],aes(Noise,NMI,group=Algo))+geom_point(size=4.5,shape=1,aes(color=Algo))+geom_line(lty=2,aes(color=Algo))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average NMI')+
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.position="bottom")+ scale_color_manual(name=TeX(' '),values=c('blue3','darkorchid4'))+
  scale_y_continuous(limits = c(0,1))
p2
ggsave(p2,file='nmi_compare_n100_bet12.pdf')

p3 <- ggplot(nmi_dt2[Beta_dist==4&Noise %in% c(2,4,6,8,10)],aes(Noise,RMSE,group=Algo))+geom_point(size=4.5,shape=1,aes(color=Algo))+geom_line(lty=2,aes(color=Algo))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average RMSE')+
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.position="bottom")+ scale_color_manual(name=TeX(' '),values=c('blue3','darkorchid4'))
p3
ggsave(p3,file='rmse_compare_n100_bet4.pdf')

p4 <- ggplot(nmi_dt2[Beta_dist==12 & Noise %in% c(2,4,6,8,10)],aes(Noise,RMSE,group=Algo))+geom_point(size=4.5,shape=1,aes(color=Algo))+geom_line(lty=2,aes(color=Algo))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average RMSE')+
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        legend.position="bottom")+ scale_color_manual(name=TeX(' '),values=c('blue3','darkorchid4'))
p4
ggsave(p4,file='rmse_compare_n100_bet12.pdf')
