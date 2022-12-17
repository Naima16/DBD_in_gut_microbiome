### This script plot the fig 4 in the paper 

library(cowplot)
library(ggplot2)

genes=read.table('~/HMP_gene_changes_full.csv',sep=',',header=T)
colnames(genes)
# [1] "sample_tp1"         "sample_tp2"         "subject"            "tp1"                "tp2"               
# [6] "species"            "num_gene_gains"     "num_gene_losses"    "alpha_div_tp1"      "alpha_div_tp2"     
# [11] "alpha_div_rare_tp1" "alpha_div_rare_tp2" "richness_tp1"       "richness_tp2"       "richness_rare_tp1" 
# [16] "richness_rare_tp2"  "polymorphism_tp1"   "polymorphism_tp2"  
#   


list_sp=c()
spec_names=unique(genes$species)
for (i in spec_names)
{
  cc=genes[genes$species==i,]
  if (all(cc$num_genes_changes == 0))
    list_sp <- c(list_sp, i)
}

###########################################
## plot the same top 9 species as in fig2 and fig3
load('~/top9.RData')
df.all3=genes[genes$species %in% top9,]

Plot_los <- ggplot(df.all3,aes(alpha_div_tp1,log10(num_gene_losses+1)))+
  geom_point(fill='gray41',size=0.9,pch=21,col='black') + 
  xlab("Shannon diversity tp1") + 
  ylab("Log10 gene loss") + 
  geom_smooth(method="lm",se=FALSE,fullrange=T,color='gray36',size=0.8)+
  theme_bw()+
  theme(axis.text=element_text(size=4),
        axis.title.y=element_blank(),
        axis.title=element_text(size=7), #element_blank(),
        strip.text = element_text(size=4),
        panel.spacing = unit(0.025, "cm"),
        strip.background = element_blank(), strip.placement = "top"
  ) +
  facet_wrap( ~ species, scales = "free",ncol=3) 
  Plot_los


Plot_los_rich <- ggplot(df.all3,aes(richness_rare_tp1,log10(num_gene_losses+1)))+
  geom_point(fill='gray41',size=0.9,pch=21,col='black') + 
  xlab("Species richness tp1") + 
  ylab("Log10 gene loss") + 
  geom_smooth(method="lm",se=FALSE,fullrange=T,color='gray36',size=0.8)+
  theme_bw()+
  theme(axis.text=element_text(size=4),
        axis.title.y=element_blank(),
        axis.title=element_text(size=7), #element_blank(),
        strip.text = element_text(size=4),
        panel.spacing = unit(0.025, "cm"),
        strip.background = element_blank(), strip.placement = "top"
  ) +
  facet_wrap( ~ species, scales = "free",ncol=3) 
Plot_los_rich


pdf('top9_gene_loss.pdf',width=5.5 ,height = 3.5,pointsize = 7)
plot_grid(Plot_los,Plot_los_rich,
             ncol=2,nrow=1,labels = c('A','B'),label_size = 6)
dev.off()


### glmm outputs

load('~/loss.nb2.shan.RData')
load('~/new_data_dec2021/raref.r1.RData')
load('~/loss.nb2.rich.sing.RData')


p1=as_grob(effect_plot(loss.nb2.shan,'alpha_div_tp1',main.title= NULL,
                       rug=FALSE,colors = 'black', y.label='Gene loss in a focal species',
                       line.thickness=0.7,
                       x.label='Shannon diversity tp1',interval = T,data=datsc,outcome.scale='link')+
             theme_classic()+
             theme(axis.title = element_text(size = rel(0.5)),
                   axis.text = element_text(size = rel(0.3)),
                   axis.title.y=element_blank())+ #+
             annotate('text',x = -2, y = 1.5,size=2, hjust = 0, 
                       label = paste(" \nP =",'0.027'),fontface = "plain"))  

p2=as_grob(effect_plot(loss.nb2.rich.sing,'richness_tp1',main.title= NULL,
                       rug=FALSE,y.label='Strain count',
                       line.thickness=0.7,colors ='black',
                       x.label='Species richness at tp1 (all data)',interval = T,data=datsc,outcome.scale='link')+
             theme_classic()+
             theme(axis.title = element_text(size = rel(0.5)),
                   axis.text = element_text(size = rel(0.3)),
                   axis.title.y=element_blank())+
             annotate('text',x = -2, y = 1.5,size=2, hjust = 0, 
                      label = paste(" \nP =",'0.034'),fontface = "plain",colors = "navyblue"))  

p3=as_grob(effect_plot(raref.r1,'richness_rare_tp1',main.title= NULL,
                       rug=FALSE,y.label='Strain count',
                       line.thickness=0.7,
                       x.label='Species richness at tp1 (rarefied data)',interval = T,colors ='black',data=datsc,outcome.scale='link' )+
             theme_classic()+
             theme(axis.title = element_text(size = rel(0.5)),
                   axis.text = element_text(size = rel(0.3)),
                   axis.title.y=element_blank())+
             annotate('text',x = -2, y = 1.5,size=2, hjust = 0, 
                      label = paste(" \nP =",'0.01'),fontface = "plain"))  


pdf('Strain_glmm.pdf',width=3,height=1.8,pointsize = 7)
plot_grid(p1,p2,p3,ncol=3,nrow=1,labels = c('C', 'D','E'),label_size = 7)
dev.off()


