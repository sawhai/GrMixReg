# Plotting the result (the "RData" file that contains the results needs to be loaded)
rm(list=ls())
setwd("/Users/ha/Documents/Gr.EM//")
res <- get(load('res_k2_d2_n800.RData'))
# load('res_k2_d4_n100.RData')
# load('res_k2_d4_n100.RData')
# load('res_k2_d4_n100.RData')
# load('res_k2_d4_n1600.RData')
library('data.table')
library('ggplot2')
library(latex2exp)
res <- res[complete.cases(res), ]
levels(res$Beta_dist)[levels(res$Beta_dist) == '0'] <- NA
res$Noise <- as.factor(res$Noise)
res$Beta_dist <- as.factor(res$Beta_dist)
res <- res[Noise %in% c(2,4,6,8,10),]
nmi_dt2 <- res[,lapply(.SD, mean, na.rm=TRUE), by=.(Noise,Beta_dist),.SDcols=c("NMI", "Beta_err",'RMSE','num.itr')]
levels(nmi_dt2$Beta_dist) <- c('4','8','12')

p1 <- ggplot(nmi_dt2,aes(Noise,NMI,group=Beta_dist))+geom_point(size=4.5,shape=1,aes(color=Beta_dist))+geom_line(lty=2,aes(color=Beta_dist))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average NMI')+
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18))+
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))#+
  scale_y_continuous(limits = c(0.2,1))
p1

ggsave(p1,width = 6.5, height = 6.5,file='nmi_k4_d4_n200.pdf')
p2 <- ggplot(nmi_dt2,aes(Noise,Beta_err,group=Beta_dist))+geom_point(size=4.5,shape=1,aes(color=Beta_dist))+geom_line(lty=2,aes(color=Beta_dist))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y=TeX('Average error $\\beta$'))+
  theme(text = element_text(size=25),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))+
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))#+scale_y_continuous(limits = c(0,1.5))
p2
ggsave(p2,width = 6.5, height = 6.5,file='beta_k4_d4_n200.pdf')
p3 <- ggplot(nmi_dt2,aes(Noise,num.itr,group=Beta_dist))+geom_point(size=4.5,shape=1,aes(color=Beta_dist))+geom_line(lty=2,aes(color=Beta_dist))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average Number of Iterations')+
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))+
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))+scale_y_continuous(trans='log2',limits = c(2,120))
p3
ggsave(p3,file='itr_k2_d2_n800_log2.pdf')
p4 <- ggplot(nmi_dt2,aes(Noise,RMSE,group=Beta_dist))+geom_point(size=4.5,shape=1,aes(color=Beta_dist))+geom_line(lty=2,aes(color=Beta_dist))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average RMSE')+
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18))+
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))

p4
ggsave(p4,width = 6.5, height = 6.5,file='rmse_k4_d4_n200.pdf')


# #Box plot
q1 <- ggplot(res,aes(Noise,NMI,fill=Beta_dist))+
  geom_boxplot(aes(fill=factor(Beta_dist)))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='NMI')+
  theme(text = element_text(size=13),panel.grid.major=element_line(colour='gray75'))+
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))+guides(fill=guide_legend(title=expression(d^2)))
q1
ggsave('nmi_box_k2_d4_n100.png')
#
q2 <- ggplot(res,aes(Noise,Beta_err,fill=Beta_dist))+
  geom_boxplot(aes(fill=factor(Beta_dist)))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Beta Error')+
  theme(text = element_text(size=13),panel.grid.major=element_line(colour='gray75'))+
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))+guides(fill=guide_legend(title=expression(d^2)))
q2
ggsave('beta_box_k2_d4_n100.png')
q3 <- ggplot(res,aes(Noise,num.itr,fill=Beta_dist))+
  geom_boxplot(aes(fill=factor(Beta_dist)))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Number of Iterations')+
  theme(text = element_text(size=13),panel.grid.major=element_line(colour='gray75'))+
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))+guides(fill=guide_legend(title=expression(d^2)))
q3
ggsave('itr_box_k2_d4_n100.png')
q4 <- ggplot(res,aes(Noise,RMSE,fill=Beta_dist))+
  geom_boxplot(aes(fill=factor(Beta_dist)))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='RMSE')+
  theme(text = element_text(size=13),panel.grid.major=element_line(colour='gray75'))+
  scale_color_manual(name=TeX('$\\delta_\\beta$'),values=c(1:3))+guides(fill=guide_legend(title=expression(d^2)))
q4
ggsave('rmse_box_k2_d4_n100.png')

