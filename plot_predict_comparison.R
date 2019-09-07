rm(list=ls())
# setwd("/Users/ha/Documents/Gr.EM/")
res <- get(load('pred_compare.RData'))

library('data.table')
library('ggplot2')
library(latex2exp)
res <- res[complete.cases(res), ]

res <- res[Beta_dist %in% c(8),]
res$Noise <- as.factor(res$Noise)
res$Beta_dist <- as.factor(res$Beta_dist)
nmi_dt2 <- res[,lapply(.SD, mean, na.rm=TRUE), by=.(Noise,Beta_dist),.SDcols=c('RMSE','RMSE_nr')]
nmi_dt2 <- melt(nmi_dt2,id.vars = c('Noise',"Beta_dist"))
nmi_dt2[,Algorithm:='GMR']
nmi_dt2[1:5,Algorithm:='FMR']

p1 <- ggplot(nmi_dt2,aes(Noise,value,group=Algorithm))+geom_point(size=4.5,shape=1,aes(color=Algorithm))+geom_line(lty=2,aes(color=Algorithm))+
  theme_bw()+labs(x=TeX('$\\sigma_k$'), y='Average RMSE')+
  theme(text = element_text(size=20),panel.grid.major=element_line(colour='gray75'),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),legend.position="bottom")+
  scale_color_manual(name=TeX(' '),values=c('blue3','darkorchid4'))
p1

