## Figure supplementary 4

library(itsadug)


##delta pi~shannon
load('/Users/naima/projet_HMP_gut_nandita_sept20/5juillet_final_version/Review_elife_mai_2022/poyet_juillet_final/pi/gam.pi2.RData')
summary(gam.pi2)


pdf('Fig_S4.pdf',width=3,height = 5,pointsize = 8)

plot_smooth(gam.pi2, view="alpha_div_tp1", plot_all='lag',cond=list(lag=c( 50,130,135,140,150,155,160,180,250)),
            col=c("lightblue1" ,'deepskyblue' ,'steelblue3','blue' ,'blue4','gold', 'darkorange' ,'tomato1' ,'#FF0000'),
            lwd=5,sim.ci=T,shade=T,se=F,
                       legend_plot_all = T,ann=F,bty='l',print.summary = F,rug = F,hide.label = T) 
title(ylab='Predicted log10 polymorphism change', cex.lab=1,line=2.5)
title(xlab='Shannon diversity tp1', cex.lab=1,line=2.5)


dev.off()
