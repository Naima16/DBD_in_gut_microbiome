### figure suppl 1
### Ploymorphism rate ~ diversity at different taxonomic levels (GAMs output)
### Naima Madi dec 2022

library(itsadug)
library(ggplot2)



#### shannon as the predictor
load('~/phyl.1.RData')
load('~/clas.1.RData')
load('~/fam.1.RData')
load('~/genus.1.RData')
load('~/order.1.RData')

## Richness as the predictor 
load('~/phyl.rich.1.RData')
load('~/clas.rich.1.RData')
load('~/order.rich.1.RData')
load('~/fam.rich.1.RData')
load('~/gen.rich.1.RData')


pdf('Figure_S1.pdf',width=6.5,height=3.2,pointsize = 10)  ## width=2 un peu trop etroit, 3.5 bonne fig
par(mfrow=c(2,5),mar=c(4,2.5,1,1))
plot_smooth(phyl.1, view="phyla_shannon", cex.main=1,rug=F,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',lwd=1.5,
            legend_plot_all = F, h0=NA,transform=plogis)
title(main='A1', cex.main = 1.2,line=0)
title(ylab='Polymorphism rate in a focal species', cex.lab=1.2,line=2.5)
title(xlab='Shannon diversity (phyla)', cex.lab=1.2,line = 2.5)
text(x=1,y=0.007,substitute(paste(italic('R'^2),'=0.189',', ',italic('P'),'=0.528')),cex=0.89)


plot_smooth(clas.1, view="class_shannon", cex.main=1,rug=F,
            print.summary = F,hide.label = T, #col='#0070e4ff'
            cex.axis=0.7,ann=F,bty='l',lwd=1.5,
            legend_plot_all = F, h0=NA,transform=plogis)
title(main='B1', cex.main = 1.2,line=0)
title(ylab='Polymorphism rate in a focal species', cex.lab=1.2,line=2.5)
title(xlab='Shannon diversity (class)', cex.lab=1.2,line = 2.5)
text(x=1,y=0.007,substitute(paste(italic('R'^2),'=0.189',', ',italic('P'),'=0.518')),cex=0.89)


plot_smooth(order.1, view="order_shannon", cex.main=1,rug=F,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',lwd=1.5,
            legend_plot_all = F, h0=NA,transform=plogis)
title(main='C1', cex.main = 1.2,line=0)
title(ylab='Polymorphism rate in a focal species', cex.lab=1.2,line=2.5)
title(xlab='Shannon diversity (order)', cex.lab=1.2,line = 2.5)
text(x=1,y=0.007,substitute(paste(italic('R'^2),'=0.189',', ',italic('P'),'=0.395')),cex=0.89)

plot_smooth(fam.1 , view="family_shannon", cex.main=1,rug=F,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',lwd=1.5,
            legend_plot_all = F, h0=NA,transform=plogis)
title(main='D1', cex.main = 1.2,line=0)
title(ylab='Polymorphism rate in a focal species', cex.lab=1.2,line=2.5)
title(xlab='Shannon diversity (family)', cex.lab=1.2,line = 2.5)
text(x=1.3,y=0.009,substitute(paste(italic('R'^2),'=0.192',', ',italic('P'),'=0.031')),cex=0.89)

plot_smooth(genus.1 , view="genus_shannon", cex.main=1,rug=F,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',lwd=1.5,
            legend_plot_all = F, h0=NA,transform=plogis)
title(main='E1', cex.main = 1.2,line=0)
title(ylab='Polymorphism rate in a focal species', cex.lab=1.2,line=2.5)
title(xlab='Shannon diversity (genus)', cex.lab=1.2,line = 2.5)
text(x=1.5,y=0.01,substitute(paste(italic('R'^2),'=0.193',', ',italic('P'),'=0.011')),cex=0.89)

plot_smooth(phyl.rich.1, view="phyla_nb", cex.main=1,rug=F,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',lwd=1.5,
            legend_plot_all = F, h0=NA,transform=plogis)
title(main='A2', cex.main = 1.2,line=0)
title(ylab='Polymorphism rate in a focal species', cex.lab=1.2,line=2.5)
title(xlab='Richness (phyla)', cex.lab=1.2,line = 2.5)
text(x=8.5,y=0.01,substitute(paste(italic('R'^2),'=0.196',', ',italic('P'),'=0.006')),cex=0.89)

plot_smooth(clas.rich.1, view="class_nb", cex.main=1,rug=F,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',lwd=1.5,
            legend_plot_all = F, h0=NA,transform=plogis)
title(main='B2', cex.main = 1.2,line=0)
title(ylab='Polymorphism rate in a focal species', cex.lab=1.2,line=2.5)
title(xlab='Richness (class)', cex.lab=1.2,line = 2.5)
text(x=10,y=0.01,substitute(paste(italic('R'^2),'=0.194',', ',italic('P'),'=0.011')),cex=0.89)

plot_smooth(order.rich.1, view="order_nb", cex.main=1,rug=F,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',lwd=1.5,
            legend_plot_all = F, h0=NA,transform=plogis)
title(main='C2', cex.main = 1.2,line=0)
title(ylab='Polymorphism rate in a focal species', cex.lab=1.2,line=2.5)
title(xlab='Richness (order)', cex.lab=1.2,line = 2.5)
text(x=20,y=0.008,substitute(paste(italic('R'^2),'=0.19',', ',italic('P'),'=0.018')),cex=0.89)

plot_smooth(fam.rich.1 , view="family_nb", cex.main=1,rug=F,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',lwd=1.5,
            legend_plot_all = F, h0=NA,transform=plogis)
title(main='D2', cex.main = 1.2,line=0)
title(ylab='Polymorphism rate in a focal species', cex.lab=1.2,line=2.5)
title(xlab='Richness (family)', cex.lab=1.2,line = 2.5)
text(x=30,y=0.009,substitute(paste(italic('R'^2),'=0.191',', ',italic('P'),'=0.002')),cex=0.89)

plot_smooth(gen.rich.1 , view="genus_nb", cex.main=1,rug=F,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',lwd=1.5,
            legend_plot_all = F, h0=NA,transform=plogis)
title(main='E2', cex.main = 1.2,line=0)
title(ylab='Polymorphism rate in a focal species', cex.lab=1.2,line=2.5)
title(xlab='Richness (genus)', cex.lab=1.2,line = 2.5)
text(x=80,y=0.009,substitute(paste(italic('R'^2),'=0.937',', ',italic('P'),'=0.004')),cex=0.89)
dev.off()

