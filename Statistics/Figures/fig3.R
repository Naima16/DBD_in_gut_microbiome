#### This script plot figure 3
#### Naima Madi dec 2022

library(cowplot)
library(jtools)
library(ggplot2)
library(forcats)
require(wesanderson)
require(tidyverse)
library(cowplot)

load('~/nb1.rich.RData')
load('~/r3.RData')
load('~/r3_raref.RData')

##plot GLMM outputs

p1=as_grob(effect_plot(r3,'shannon_diversity',main.title= NULL,
                       rug=FALSE,colors = 'black', y.label='Strain count in a focal species',
                       line.thickness=0.7,
                       x.label='Shannon diversity',interval = T)+
                       theme_classic()+
                       theme(axis.title = element_text(size = rel(0.5)),
                       axis.text = element_text(size = rel(0.3)),
                       axis.title.y=element_blank())+ #+
                       annotate('text',x = -2, y = 0.55,size=2, hjust = 0, 
                       label = paste(" \nP =",'3.579e-07'),fontface = "plain"))  #  4.83e-07

p2=as_grob(effect_plot(nb1.rich,'species_richness',main.title= NULL,
                       rug=FALSE,y.label='Strain count',
                       line.thickness=0.7,colors ='black',
                       x.label='Species richness (all data)',interval = T)+
                       theme_classic()+
                       theme(axis.title = element_text(size = rel(0.5)),
                       axis.text = element_text(size = rel(0.3)),
                       axis.title.y=element_blank())+
                       annotate('text',x = -2, y = 0.55,size=2, hjust = 0, 
                       label = paste(" \nP =",'1.675e-06'),fontface = "plain",colors = "navyblue"))  # 3.56e-06

p3=as_grob(effect_plot(r3_raref,'species_richness',main.title= NULL,
                       rug=FALSE,y.label='Strain count',
                       line.thickness=0.7,
                       x.label='Species richness (rarefied data)',interval = T,colors ='black' )+
                       theme_classic()+
                       theme(axis.title = element_text(size = rel(0.5)),
                       axis.text = element_text(size = rel(0.3)),
                       axis.title.y=element_blank())+
                       annotate('text',x = -2, y = 0.55,size=2, hjust = 0, 
                       label = paste(" \nP =",'1.188e-05'),fontface = "plain"))  #1.188e-05  1.42e-05


pdf('Strain_glmm.pdf',width=3,height=1.8,pointsize = 7)
plot_grid(p1,p2,p3,ncol=3,nrow=1,labels = c('C', 'D','E'),label_size = 7)
dev.off()


## Scatter plot for strain number~community diversity for the 9 most prevalent species

gene_strain_df=read.table('~/gene_counts_strain_num_per_sample_species.csv',sep=',',header=T)
length(unique(gene_strain_df$sample_id))  ## 469
gene_strain_df$sample_id=as.character(gene_strain_df$sample_id)

for (i in 1:dim(gene_strain_df)[1] )
{
  if (grepl('c',gene_strain_df[i,'sample_id']))  
  {
    name=gene_strain_df[i,'sample_id']
    gene_strain_df[i,'sample_id']=strsplit(name, "c")[[1]][1]
    
  }
}

strain_df=gene_strain_df[,c('sample_id','species_name','strain_number')]
coverage=read.table('~/depths_geq_20_TF_copie.csv',sep=',',header=T,row.names = 1, stringsAsFactors = FALSE, quote = "", comment.char = "",check.names = F)

###filter out samples and species < X20 coverage
for (i in 1:dim(strain_df)[1])
{
  my_sample=strain_df[i,'sample_id']
  my_Sp=strain_df[i,'species_name']
  if (my_Sp %in% rownames(coverage) & my_sample %in% colnames(coverage))
  {
    if (coverage[rownames(coverage)==my_Sp,colnames(coverage)==my_sample]=='False')
      strain_df[i,'strain_number']=0
  }
}

##filter out samples/species with no  enough coverage
strain_df1=strain_df[strain_df$strain_number>0,]

sample_data.1=read.table('~/sample_covariates_rarefied_20m.txt',sep='\t',header=T)

####### 
df.strain2=merge(strain_df1,sample_data1[,c('sample_id','subject_id','shannon_diversity','species_richness')],by='sample_id')
colnames(df.strain2)[2]='species_id'
colnames(df.strain2)[3]='strain_nb'

## plot exactly the same species as polymorphism rate ~ diversity plot
load('~/top9.RData')
to_plot_strain=df.strain2[df.strain2$species_id %in% top9,]

Plot_shan <- ggplot(to_plot_strain,aes(shannon_diversity,strain_nb))+
  geom_point(fill='gray41',size=0.9,pch=21,col='black') + 
  xlab("Shannon diversity") + 
  ylab("Strain count") + 
  geom_smooth(method="lm",se=FALSE,fullrange=T,color='gray36',size=0.8)+
  theme_bw()+
  theme(axis.text=element_text(size=4),
        axis.title.y=element_blank(),
        axis.title=element_text(size=7), #element_blank(),
        strip.text = element_text(size=4),
        panel.spacing = unit(0.025, "cm"),
        strip.background = element_blank(), strip.placement = "top"
  ) +
  facet_wrap( ~ species_id, scales = "free",ncol=3) 
Plot_shan


Plot_div <- ggplot(to_plot_strain,aes(species_richness,strain_nb))+
  geom_point(fill='gray41',size=0.9,pch=21,col='black') + 
  xlab("Shannon diversity") + 
  ylab("Strain count") + 
  geom_smooth(method="lm",se=FALSE,fullrange=T,color='gray36',size=0.8)+
  theme_bw()+
  theme(axis.text=element_text(size=4),
        axis.title.y=element_blank(),
        axis.title=element_text(size=7), #element_blank(),
        strip.text = element_text(size=4),
        panel.spacing = unit(0.025, "cm"),
        strip.background = element_blank(), strip.placement = "top"
  )+
  facet_wrap( ~ species_id, scales = "free",ncol=3) 
Plot_div

pdf('top9_strain_diversity.pdf',width=5.5 ,height = 3.5,pointsize = 8)
plot_grid(Plot_shan,Plot_div,
             ncol=2,nrow=1,labels = c('A','B'),label_size = 6)
dev.off()

###merge the two plots (glmm output and the scatter plot) with inkscape.
