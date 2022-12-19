### Figure S2 : plot GAM outputs with polymorhism rate at non synonymous sites as a function of diversity at hogh taxonomic levels
### Naima Madi dec 2022

library(itsadug)


## load GAM models with polymrphism rate (at non synonymous sites) as a function of Shannon diversity
load('~/shan_tax/Nphyl.1.sept.RData')
load('~/Nclas.1.sept.RData')
load('~/Norder.1.sept.RData')
load('~/Nfam.1.sept.RData')
load('~/Ngenus.1.sept.RData')

## load GAM models with polymrphism rate (at non synonymous sites) as a function of richness
load('~/Nphyl.rich.1.sept.RData')
load('~/Nclas.rich.1.sept.RData')
load('~/rich_tax/Norder.rich.1.sept.RData')
load('~/Nfam.rich.1.sept.RData')
load('~/Ngenus.rich.1.sept.Rdata')


pdf('Figure_S2.pdf',width=6.5,height=3.2,pointsize = 10)  ## width=2 un peu trop etroit, 3.5 bonne fig
par(mfrow=c(2,5),mar=c(4,2.5,1,1))
plot_smooth(Nphyl.1.sept, view="phyla_shannon", cex.main=1,rug=F,
            print.summary = F,hide.label = T,lwd=1.5,
            cex.axis=0.7,ann=F,bty='l',legend_plot_all = F, h0=NA,
            transform=plogis)
title(main='A1', cex.main = 1.2,line=0)
title(xlab='Phylum Sahnnon', cex.lab=1.2,line=2.5)
text(1,0.0009,substitute(paste(italic('R'^2),'=0.241',', ',italic('P'),'=0.284')),cex=0.89)


plot_smooth(Nclas.1.sept, view="class_shannon", cex.main=1,rug=F,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',legend_plot_all = F, lwd=1.5,
            h0=NA,transform=plogis)
title(xlab='Class Shannon', cex.lab=1.2,line=2.5)
title(main='B1', cex.main = 1.2,line=0)
text(1,0.0016,substitute(paste(italic('R'^2),'=0.24',', ',italic('P'),'= 0.286')),cex=0.89)

plot_smooth(Norder.1.sept, view="order_shannon", cex.main=1,rug=F, lwd=1.5,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',legend_plot_all = F, h0=NA,transform=plogis)
title(xlab='Order Shannon', cex.lab=1.2,line=2.5)
title(main='C1', cex.main = 1.2,line=0)
text(1,0.0009,substitute(paste(italic('R'^2),'= 0.242',', ',italic('P'),'=0.138')),cex=0.89)

plot_smooth(Nfam.1.sept, view="family_shannon", cex.main=1,rug=F,lwd=1.5,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',legend_plot_all = F, h0=NA,transform=plogis)
text(1.5,0.0009,substitute(paste(italic('R'^2),'= 0.244',', ',italic('P'),'=0.842')),cex=0.89)
title(xlab='Family Shannon', cex.lab=1.2,line=2.5)
title(main='D1', cex.main = 1.2,line=0)

plot_smooth(Ngenus.1.sept, view="genus_shannon", cex.main=1,rug=F,
            print.summary = F,hide.label = T,lwd=1.5,
            cex.axis=0.7,ann=F,bty='l',legend_plot_all = F, h0=NA,transform=plogis)
text(2.2,0.0009,substitute(paste(italic('R'^2),'= 0.247',', ',italic('P'),'=0.871')),cex=0.89)
title(xlab='Genus Shannon', cex.lab=1.2,line=2.5)
title(main='E1', cex.main = 1.2,line=0)

plot_smooth(Nphyl.rich.1.sept, view="phyla_nb", cex.main=1,rug=F,
            print.summary = F,hide.label = T,lwd=1.5,
            cex.axis=0.7,ann=F,bty='l',legend_plot_all = F, h0=NA,transform=plogis)
title(main='A2', cex.lab=1.2,line=2.5)
title(xlab='Phylum richness', cex.lab=1.2,line=2.5)
text(11,0.00105,substitute(paste(italic('R'^2),'=0.246',', ',italic('P'),'=3e-04')),cex=0.89)


plot_smooth(Nclas.rich.1.sept, view="class_nb", cex.main=1,rug=F,
            print.summary = F,hide.label = T,lwd=1.5,
            cex.axis=0.7,ann=F,bty='l',legend_plot_all = F, h0=NA,transform=plogis)
title(main='B2', cex.main = 1.2,line=0)
title(xlab='Class richness', cex.lab=1.2,line=2.5)
text(11,0.001,substitute(paste(italic('R'^2),'=0.244',', ',italic('P'),'=0.017')),cex=0.89)

plot_smooth(Norder.rich.1.sept, view="order_nb", cex.main=1,rug=F,
            print.summary = F,hide.label = T,lwd=1.5,
            cex.axis=0.7,ann=F,bty='l',legend_plot_all = F, h0=NA,transform=plogis)
title(main='C2', cex.main = 1.2,line=0)
title(xlab='Order richness', cex.lab=1.2,line=2.5)
text(22,0.001,substitute(paste(italic('R'^2),'= 0.244',', ',italic('P'),'=6.11e-04')),cex=0.89)

plot_smooth(Nfam.rich.1.sept, view="family_nb", cex.main=1,rug=F,
            print.summary = F,hide.label = T,lwd=1.5,
            cex.axis=0.7,ann=F,bty='l',legend_plot_all = F, h0=NA,transform=plogis)
title(main='D2', cex.main = 1.2,line=0)
title(ylab='')
text(40,0.001,substitute(paste(italic('R'^2),'= 0.243',', ',italic('P'),'=0.241')),cex=0.89)
title(xlab='Family richness', cex.lab=1.2,line=2.5)

plot_smooth(Ngenus.rich.1.sept, view="genus_nb", cex.main=1,rug=F,
            print.summary = F,hide.label = T,lwd=1.5,
            cex.axis=0.7,ann=F,bty='l',legend_plot_all = F, h0=NA,transform=plogis)
title(main='E2', cex.main = 1.2,line=0)
text(90,0.001,substitute(paste(italic('R'^2),'= 0.243',', ',italic('P'),'=0.122')),cex=0.89)
title(xlab='Genus richness', cex.lab=1.2,line=2.5)

dev.off()


