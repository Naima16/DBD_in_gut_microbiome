#### This script plot figure 2. 1: Scatter plot showing the relationship between polymorphism rate 
#### and community diversity in the 9 most prevalent species in the raw data, 
#### 2: Gam outputs

Garud.data=read.table('~/HMP1-2_polymorphism_rates_alpha_divs.csv',sep=',',header=T)
dim(Garud.data)
head(Garud.data)
colnames(Garud.data)
# [1] "sample_id"               "subject_id"              "species_name"            "polymorphism_rate"      
# [5] "Shannon_alpha_diversity"

Garud.data$sample_id=as.factor(Garud.data$sample_id)
Garud.data$species_name=as.factor(Garud.data$species_name)
Garud.data$subject_id=as.factor(Garud.data$subject_id)
Garud.data[grepl('c',Garud.data$sample_id),'sample_id']


sample_data_raref=read.table('~/sample_covariates_rarefied_20m.txt',sep='\t',header=T)
colnames(sample_data_raref)
#"sample_id"         "subject_id"        "shannon_diversity" "species_richness"  "total_reads_MIDAS" "total_reads_orig" "

dim(sample_data_raref)


Garud.data.2=merge(Garud.data,sample_data_raref[,c('sample_id','species_richness')],by='sample_id')
colnames(Garud.data.2)

sp.df=data.frame(table(Garud.data.2$species_name))
head(sp.df)
sp.ord=sp.df[with(sp.df,order(-Freq)),]


####################
#########################3
###   biplot 9 most prevalent

c=as.data.frame(table(Garud.data.2$species_name))

c.ord=as.data.frame(c[with(c,order(-Freq)),])

mydata1=Garud.data.2[which(Garud.data.2$species_name %in% c.ord[1:9,'Var1']),]

library(scales)

Plot2 <- ggplot(mydata1,aes(Shannon_alpha_diversity,polymorphism_rate))+
  geom_point(color='gray41',size=0.07) + 
  xlab("Shannon diversity") + 
  ylab("Polymorphism rate") + 
  geom_smooth(method="lm",se=FALSE,fullrange=T,color='gray36',size=0.8)+
  theme_bw()+
  theme(axis.text=element_text(size=4),
        axis.title.y=element_blank(),
        axis.title=element_text(size=7), #element_blank(),
        strip.text = element_text(size=4),
        panel.spacing = unit(0.025, "cm"),
        strip.background = element_blank(), 
        strip.placement = "top"
  ) +
  facet_wrap( ~ species_name, scales = "free",ncol=3) +
  scale_y_continuous(trans=log10_trans())
Plot2

library(scales)
Plot1 <- ggplot(mydata1,aes(species_richness,polymorphism_rate))+
  geom_point(color='gray41',size=0.07) + 
  xlab("Species richness") + 
  ylab("Polymorphism rate") + 
  geom_smooth(method="lm",se=FALSE,fullrange=T,color='gray36',size=0.8)+
  theme_bw()+
  theme(axis.text=element_text(size=4),
        axis.title.y=element_blank(), 
        axis.title=element_text(size=7), #element_blank(),
        strip.text = element_text(size=4),
        panel.spacing = unit(0.025, "cm"),
        strip.background = element_blank(), strip.placement = "top"
  )+
  facet_wrap( ~ species_name, scales = "free",ncol=3) +
  scale_y_continuous(trans=log10_trans())
Plot1

pdf('top9.pdf',width=5.5 ,height = 3.5,pointsize = 7)
p0=plot_grid(Plot2,Plot1,
             ncol=2,nrow=1,labels = c('A','B'),label_size = 6)
p0
dev.off()

##gam plot
library(itsadug)
load('~/gam1.pi.fs.RData')
load('~/gam.pi.rich.all.RData')
load('~/gam.pi.rarefRich.cov.RData')

pdf('gam_2.pdf',width=3.5 ,height = 2,pointsize = 7)
par(mfrow=c(1,3))

plot_smooth(gam1.pi.fs, view="Shannon_alpha_diversity", cex.main=1,rug=F,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',legend_plot_all = F, h0=NA,transform=plogis)
title(main='C', cex.main = 1.4,line=0)
title(ylab='Polymorphism rate in a focal species', cex.lab=1.3,line=2.1)
title(xlab='Shannon diversity',cex.lab=1.3,line = 2.1)
text(1.9,0.011,substitute(paste(italic('R'^2),'=0.198',', ',italic('P'),'=0.031')),cex=1)


plot_smooth(gam.pi.rich.all, view="species_richness", cex.main=1,rug=F,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',legend_plot_all = F, h0=NA,transform=plogis)
title(main='D', cex.main = 1.4,line=0)
title(ylab='')

title(xlab='Species richness (all data)',cex.lab=1.3,line = 2.1)
text(154,0.0097,substitute(paste(italic('R'^2),'=0.191',', ',italic('P'),'=0.017')),cex=1)


plot_smooth(gam.pi.rarefRich.cov, view="species_richness", cex.main=1,rug=F,
            print.summary = F,hide.label = T,
            cex.axis=0.7,ann=F,bty='l',legend_plot_all = F, h0=NA,transform=plogis)
title(main='E', cex.main = 1.4,line=0)
title(ylab='')
title(xlab='Species richness (rarefied data)',cex.lab=1.3,line = 2.1)
text(120,0.0082,substitute(paste(italic('R'^2),'=0.19',', ',italic('P'),'=2.63e-04')),cex=1)

dev.off()

## merge the two plots in inkscape
