### script for figure 5 
### Poyet data mixed models outputs
### Naima Madi dec 2022

library(itsadug)
library(sjPlot)

##delta pi~shannon div gam
load('~/gam.pi2.RData')
summary(gam.pi2)


## gene loss~shannon div glmm
load('~/loss_glmmTMB2.RData')
summary(loss_glmmTMB2)


## gene gain~shannon div glmm
load('~/gain_shan.RData')
summary(gain_shan)


#### gene gain ~ richness
load('~/gain_rich.RData')
summary(gain_rich)


## pi ~ Shannon
pdf('pi_shannon_quartile.pdf',width = 2,height = 3.2,pointsize = 8)
plot_smooth(gam.pi2, view="alpha_div_tp1", plot_all='lag',cond=list(lag=c(50 , 130 , 250)),
            col=c('slategray3','steelblue','navyblue'),lwd=5,sim.ci=T,shade=T,se=F,
            legend_plot_all = F,ann=F,bty='l',print.summary = F,rug = F,hide.label = T) #,transform=exp)#,transform=plogis)
title(ylab='Predicted log10 polymorphism change', cex.lab=1,line=2.5)
title(xlab='Shannon diversity tp1', cex.lab=1,line=2.5)
dev.off()

##gene loss~Shannon
pdf('genesloss_shannon_quartile_v2.pdf',width = 2,height = 3,pointsize = 6)
plot_model(
  loss_glmmTMB2, type = "pred", terms = c("alpha_div_tp1", "lag"),mdrt.values='quart',
  show.legend =F,
  show.values = F,transform=exp,
  value.offset = .4,
  ci.lvl = NA,
  value.size = 4,
  dot.size = 3,
  vline.color = "blue",title='Gene loss',
  axis.title = c('Shannon diversity in tp1','Predicted number of gene lost'),
  width = 1.5,colors =c('slategray3','steelblue','navyblue'),line.size=1.7)+
  theme_sjplot2(base_size = 7, base_family = "")
dev.off()


## gene gain ~ Shannon
pdf('genesgain_shannon_quartile_v1.pdf',width = 2,height = 3,pointsize = 6)
plot_model(
  gain_shan,
  type = "pred", terms = c("alpha_div_tp1", "lag"),
  mdrt.values='quart',show.legend =F,
  show.values = TRUE,
  transform=exp,
  value.offset = .4,
  value.size = 2,
  dot.size = 3,
  ci.lvl = NA,
  vline.color = "blue",
  width = 1.5,
  colors =c('slategray3','steelblue','navyblue'),line.size=1.7,
  axis.title = c('Shannon diversity in tp1','Predicted number of genes gained'))+
  theme_sjplot2(base_size = 7, base_family = "")
dev.off()


## gene gain ~ richness
pdf('genesgain_rich_quartile_v2.pdf',width = 2,height = 3,pointsize = 6)
plot_model(
  gain_rich, type = "pred", terms = c("richness_rare_tp1", "lag"),mdrt.values='quart',
  show.legend =F,
  show.values = F,transform=exp,
  value.offset = .4,
  ci.lvl = NA,
  value.size = 4,
  dot.size = 3,
  vline.color = "blue",title='Gene gain',
  axis.title = c('Species richness in tp1','Predicted number of genes gained'),
  width = 1.5,colors =c('slategray3','steelblue','navyblue'),line.size=1.7)+
  theme_sjplot2(base_size = 7, base_family = "")
dev.off()

## merge plots with inkscape
