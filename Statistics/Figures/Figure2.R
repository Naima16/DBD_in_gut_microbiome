
polymorphism.df=read.table('HMP1-2_polymorphism_rates_alpha_divs.csv',sep=',',header=T)
dim(polymorphism.df)
head(polymorphism.df)
colnames(polymorphism.df)
# [1] "sample_id"               "subject_id"              "species_name"            "polymorphism_rate"      
# [5] "Shannon_alpha_diversity"

polymorphism.df$sample_id=as.factor(polymorphism.df$sample_id)
polymorphism.df$species_name=as.factor(polymorphism.df$species_name)
polymorphism.df$subject_id=as.factor(polymorphism.df$subject_id)
polymorphism.df[grepl('c',polymorphism.df$sample_id),'sample_id']


sample_data_raref=read.table('sample_covariates_rarefied_20m.txt',sep='\t',header=T)
colnames(sample_data_raref)
# [1] "sample_id"         "subject_id"        "shannon_diversity" "species_richness"  "total_reads_MIDAS"
# [6] "total_reads_orig"
dim(sample_data_raref)


colnames(sample_data_raref)[4]='species_richness_raref'

polymorphism.df.2=merge(polymorphism.df,sample_data_raref[,c('sample_id','species_richness_raref')],by='sample_id')
colnames(polymorphism.df.2)


length(unique(polymorphism.df.2$sample_id))  ## 467

sp.df=data.frame(table(polymorphism.df.2$species_name))
head(sp.df)
sp.ord=sp.df[with(sp.df,order(-Freq)),]
head(sp.ord)
dim(sp.df) ##62


######### all the data in 1 plot
library(ggpubr)
###plot all species in the same plot


###########################################
######################################
###plot only significant correltions
prev.sp=as.data.frame(table(polymorphism.df.2$species_name))
length(unique(polymorphism.df.2$species_name))  #69
prev.sp.ord=prev.sp[order(prev.sp$Freq),]
to.drop=prev.sp.ord[prev.sp.ord$Freq<4,'Var1']  #1
##drop species with not enough observations
df.all3=polymorphism.df.2[!polymorphism.df.2$species_name %in% to.drop,]
dim(df.all3)
dim(polymorphism.df.2)

length(unique(df.all3$species_name))  #68
length(unique(polymorphism.df.2$species_name))  ##69

 library(dplyr)

corfun<-function(x, y) {
  corr=(cor.test(x, y, method="spearman"))
}

length(unique(df.all3$species_name))
#[1] 68

stat.shan1=ddply(df.all3, .(species_name), summarise,
                 pval=corfun(Shannon_alpha_diversity,polymorphism_rate)$p.value,
                 est=corfun(Shannon_alpha_diversity,polymorphism_rate)$estimate) 

sig.pval.pi=stat.shan1[stat.shan1$pval<=0.05,'species_name']
length(sig.pval.pi) ##15
# [1] Alistipes_onderdonkii_55464        Bacteroides_fragilis_54507        
# [3] Bacteroides_thetaiotaomicron_56941 Bacteroides_uniformis_57318       
# [5] Bacteroides_vulgatus_57955         Bacteroides_xylanisolvens_57185   
# [7] Burkholderiales_bacterium_56577    Faecalibacterium_cf_62236         
# [9] Faecalibacterium_prausnitzii_62201 Lachnospiraceae_bacterium_51870   
# [11] Oscillibacter_sp_60799             Oscillospiraceae_bacterium_54867  
# [13] Parabacteroides_distasonis_56985   Parabacteroides_merdae_56972      
# [15] Roseburia_inulinivorans_61943 

stat.shan1[stat.shan1$pval<0.05 & stat.shan1$est<0 ,'species_name']
#[1] Lachnospiraceae_bacterium_51870

stat.shan1[stat.shan1$pval<0.05 & stat.shan1$est>0 ,] # 14out of 15

to.plot.shan=df.all3[df.all3$species_name%in% sig.pval.pi ,]
length(unique(to.plot.shan$species_name)) #15


library(ggplot2)

library(RColorBrewer)
nb.cols <- 15
myColors.pi.shan <- colorRampPalette(brewer.pal(4, "Spectral"))(nb.cols)

shan_pi=unique(to.plot.shan$species_name)
names(myColors.pi.shan) <- shan_pi

myColors.pi.shan[names(myColors.pi.shan)=='Lachnospiraceae_bacterium_51870']='black'

colScale.pi.shan <- scale_colour_manual(name = "shan_pi",values = myColors.pi.shan)


sp.shan <- ggplot(to.plot.shan,aes(x=Shannon_alpha_diversity,y=log10(polymorphism_rate),color=species_name))+
  geom_point(size=0.5,show.legend = F) + xlab("Shannon") + ylab("Log polymorphism rate") + 
  guides( size = FALSE)+
  geom_smooth(data=to.plot.shan,method=lm, formula=y~x,se=FALSE,size=0.5) +  
  
  colScale.pi.shan+
  theme_bw()+
  
  theme(
        legend.position = 'right',
        legend.title = element_blank(),
        legend.text = element_text (size=3.5),
        legend.direction= 'vertical',
        #legend.spacing.y = unit(0.01, 'cm'),
        legend.key.size = unit(0.1, "cm"),
        axis.text=element_text(size=4),
        axis.title=element_text(size=7)) 
 

sp.shan2=sp.shan+guides(col=guide_legend(ncol=1),label.vjust=0.5)
sp.shan2

ggsave('shan_pi_v1.pdf',width=7,height = 5)


################

stat.rich2=ddply(df.all3, .(species_name), summarise,
      pval=corfun(species_richness_raref, polymorphism_rate)$p.value,
      est=corfun(species_richness_raref, polymorphism_rate)$estimate) 

stat.rich2[stat.rich2$pval<0.05 ,] 
dim(stat.rich2[stat.rich2$pval<0.05 ,]) #18

stat.rich2[stat.rich2$pval<0.05 & stat.rich2$est<0 ,'species_name']
#0

#stat.rich2[stat.rich2$est>0,'species_name']
sig.pval.rich=stat.rich2[stat.rich2$pval<0.05,'species_name']
#sig.pval.rich=sig.pval.rich[!is.na(sig.pval.rich)]
to.plot.rich=df.all3[df.all3$species_name %in% sig.pval.rich,]

####
library(RColorBrewer)
nb.cols <- 18
myColors.pi.rich <- colorRampPalette(brewer.pal(4, "Spectral"))(nb.cols)

#myColors <- brewer.pal(25,"Spectral")
names(myColors.pi.rich) <- sig.pval.rich

for (i in 1:length(myColors.pi.rich))
{
  name_i=names(myColors.pi.rich)[i]
  if (name_i %in% sig.pval.pi)
  {
    myColors.pi.rich[names(myColors.pi.rich)==name_i]=myColors.pi.shan[names(myColors.pi.shan)==name_i]
  }
}



colScale.pi.rich <- scale_colour_manual(name = "sig.pval.rich",values = myColors.pi.rich)


to.plot.rich$species_name=as.factor(to.plot.rich$species_name)

sp.rich <- ggplot(to.plot.rich,aes(x=species_richness_raref, y=log(polymorphism_rate),color=species_name))+
  geom_point(size=0.5,show.legend = F) + xlab("Species richness") +  
  guides( size = FALSE)+
  geom_smooth(data=to.plot.rich,method=lm, formula=y~x,se=FALSE,size=0.5) +
  colScale.pi.rich+
  theme_bw()+
  theme(legend.position = 'right',
        legend.title = element_blank(),
        legend.text = element_text (size=3.5),
        legend.direction= 'vertical',
        legend.key.size = unit(0.1, "cm"),
        axis.text=element_text(size=4),
        axis.title=element_text(size=7),
        
        axis.title.y=element_blank())
sp.rich2=sp.rich+guides(col=guide_legend(ncol=1),label.vjust=0.5)

sp.rich2
ggsave('rich_pi_v1.pdf',width=7,height = 5)



ggarrange(sp.shan2, sp.rich2, widths = c(2,2),labels = c("A", "B"))
ggsave('rich_shan_pi_v1.pdf',width=10,height=4)

####### plot significant
###
mydata=df.all3[df.all3$species_name %in%sig.pval.pi,]

Plot2 <- ggplot(mydata,aes(x=Shannon_alpha_diversity,y=polymorphism_rate))+
  geom_point(color='blue2',size=0.5) + 
  xlab("Shannon") + 
  ylab("Log polymorphism rate") + 
  geom_smooth(method="loess",se=FALSE,fullrange=T,color='#FF0000',size=1)+
  #geom_smooth(method="loess",se=FALSE,fullrange=T,color='green2',size=0.4)+
  
  theme_bw()+
  theme(axis.text=element_text(size=4),
        axis.title=element_text(size=7), #element_blank(),
        strip.text = element_text(size=4),
        panel.spacing = unit(0.025, "cm"),
        strip.background = element_blank(), strip.placement = "outside")+
  facet_wrap( ~ species_name, scales = "free",ncol=4) +
  # stat_cor(size=1.4,label.x.npc = "left",
  #          label.y.npc = "top",method='spearman') 

  stat_cor(size=1.4,label.x.npc = "left",
         label.y.npc = "top",method='spearman',color = 'red',r.accuracy=0.01,p.accuracy=0.001)+
  #panel.background = element_rect(fill = 'lightblue', color = 'purple'))+
  
  scale_y_continuous(trans=log10_trans())
Plot2

ggsave('sig_pi_shannon_hmp_v4.pdf',width=5,height = 5 )  #with log10

ggsave('sig_pi_shannon_hmp_v3.pdf',width=5,height = 5 )

##plot significant richness
mydata.rich=df.all3[df.all3$species_name %in%sig.pval.rich,]

Plot2.rich <- ggplot(mydata.rich,aes(x=species_richness_raref,y=polymorphism_rate))+
#Plot2.rich <- ggplot(mydata,aes(x=Shannon_alpha_diversity,y=polymorphism_rate))+
  geom_point(color='blue2',size=0.5) + 
  xlab("Shannon") + 
  ylab("Log polymorphism rate") + 
  geom_smooth(method="loess",se=FALSE,fullrange=T,color='#FF0000',size=1)+
  #geom_smooth(method="loess",se=FALSE,fullrange=T,color='green2',size=0.4)+
  
  theme_bw()+
  theme(axis.text=element_text(size=4),
        axis.title=element_text(size=7), #element_blank(),
        strip.text = element_text(size=4),
        panel.spacing = unit(0.025, "cm"),
        strip.background = element_blank(), strip.placement = "outside")+
  facet_wrap( ~ species_name, scales = "free",ncol=4) +
  # stat_cor(size=1.4,label.x.npc = "left",
  #          label.y.npc = "top",method='spearman') 
  
  stat_cor(size=1.4,label.x.npc = "left",
           label.y.npc = "top",method='spearman',r.accuracy=0.01,p.accuracy=0.001)+
  #panel.background = element_rect(fill = 'lightblue', color = 'purple'))+
  
  scale_y_continuous(trans=log10_trans())
Plot2.rich

#ggsave('sig_pi_shannon_hmp_v4.pdf',width=5,height = 5 )  #with log10

ggsave('sig_pi_richness_hmp_v4.pdf',width=5,height = 5.5 )  #with log10

#ggsave('sig_pi_richness_hmp_v3.pdf',width=5,height = 5 )



###biplot 10 most prevalent

c=as.data.frame(table(df.all3$species_name))
dim(c)
length(unique(df.all3$species_name))  ##68

c.ord=as.data.frame(c[with(c,order(-Freq)),])

mydata1=df.all3[which(df.all3$species_name %in% c.ord[1:10,'Var1']),]
dim(mydata1)



Plot2 <- ggplot(mydata1,aes(Shannon_alpha_diversity,polymorphism_rate))+
  geom_point(color='blue2',size=0.07) + 
  xlab("Shannon") + 
  ylab("Log polymorphism rate") + 
  geom_smooth(method="loess",se=FALSE,fullrange=T,color='gray53',size=0.8)+
  #geom_smooth(method="loess",se=FALSE,fullrange=T,color='green2',size=0.4)+
  
  theme_bw()+
  
  theme(axis.text=element_text(size=4),
        axis.title=element_text(size=7), #element_blank(),
        strip.text = element_text(size=4),
        panel.spacing = unit(0.025, "cm"),
        strip.background = element_blank(), strip.placement = "top",
  )+
  facet_wrap( ~ species_name, scales = "free",ncol=2) +
  stat_cor(size=1.4,label.x.npc = "left",
           label.y.npc = "top",method='spearman',color = 'red',r.accuracy=0.01,p.accuracy=0.001)+
  #panel.background = element_rect(fill = 'lightblue', color = 'purple'))+
  
  scale_y_continuous(trans=log10_trans())
Plot2
library(scales)
Plot1 <- ggplot(mydata1,aes(species_richness_raref,polymorphism_rate))+
  geom_point(color='blue2',size=0.07) + 
  xlab("Species richness") + 
  ylab("Log polymorphism rate") + 
  geom_smooth(method="loess",se=FALSE,fullrange=T,color='gray53',size=0.8)+
  #geom_smooth(method="loess",se=FALSE,fullrange=T,color='green2',size=0.4)+
  
  theme_bw()+
 
  theme(axis.text=element_text(size=4),
        axis.title=element_text(size=7), #element_blank(),
        strip.text = element_text(size=4),
        panel.spacing = unit(0.025, "cm"),
        strip.background = element_blank(), strip.placement = "top",
        )+
  facet_wrap( ~ species_name, scales = "free",ncol=2) +
  stat_cor(size=1.4,label.x.npc = "left",
           label.y.npc = "top",method='spearman',color = 'red',r.accuracy=0.01,p.accuracy=0.001)+
        #panel.background = element_rect(fill = 'lightblue', color = 'purple'))+
   
  scale_y_continuous(trans=log10_trans())
Plot1
pdf('biplot_1otop_v1.pdf',width=4.5 ,height = 4.5,pointsize = 7)
#### 
grid.arrange(Plot2,Plot1,ncol=2)
dev.off()


####


################
#################

###strain number

gene_strain_df=read.table('gene_counts_strain_num_per_sample_species.csv',sep=',',header=T)
dim(gene_strain_df)
head(gene_strain_df)

colnames(gene_strain_df)
gene_strain_df$sample_id=as.character(gene_strain_df$sample_id)
###enlève le c à la fin des noms de certains samples
for (i in 1:dim(gene_strain_df)[1] )
{
  if (grepl('c',gene_strain_df[i,'sample_id']))  
  {
    name=gene_strain_df[i,'sample_id']
    gene_strain_df[i,'sample_id']=strsplit(name, "c")[[1]][1]
    
  }
}

gene_strain_df[grepl('c',gene_strain_df$sample_id),'sample_id']
colnames(gene_strain_df)
dim(gene_strain_df)
gene_strain_df$sample_id=as.factor(gene_strain_df$sample_id)
length(unique(gene_strain_df$sample_id)) ##467
length(unique(gene_strain_df$species_name)) ##190
#########

###sample data from strain not raref
sample_data=read.table('sample_covariates.txt',sep='\t',header=T)
colnames(sample_data)
# "sample_id"                     "subject_id"               
# [3] "shannon_diversity"         "species_richness"         
# [5] "total_reads_MIDAS"         "total_reads_orig"         
# [7] "total_avg_marker_coverage" "total_med_marker_coverage"


sample_data[grepl('c',sample_data$sample_id),'sample_id']

dim(sample_data)
length(unique(sample_data$sample_id))  ##469

sample_data[!sample_data$sample_id %in% gene_strain_df$sample_id,'sample_id']
###700038918 700165665

##raref data

sample_data_raref=read.table('sample_covariates_rarefied_20m.txt',sep='\t',header=T)
colnames(sample_data_raref)
range(sample_data_raref$total_reads_orig) #70628 20000000
small.samples=sample_data_raref[sample_data_raref$total_reads_orig<20000000,'sample_id']

plot(sample_data$shannon_diversity~sample_data_raref$shannon_diversity)

dim(sample_data) #469
dim(sample_data_raref) #469

##richness is on rarefied data
saple_Data_comb=merge(sample_data[,c('sample_id','subject_id','shannon_diversity','total_reads_orig')],sample_data_raref[c('sample_id','species_richness')],by='sample_id')

#######
df.strain=merge(gene_strain_df,saple_Data_comb[,c('sample_id','subject_id','shannon_diversity','species_richness','total_reads_orig')],by='sample_id')


colnames(df.strain)[2]='species_id'
colnames(df.strain)[10]='strain_nb'
colnames(df.strain)[13]='species_richness_raref'
colnames(df.strain)

dim(df.strain)
length(unique(df.strain$species_id))  #194
length(unique(df.strain$sample_id))  #469
length(unique(df.strain$subject_id))  ###249
length(unique(sample_data$sample_id))  ##469


df.strain$species_id=as.factor(df.strain$species_id)  
df.strain$sample_id=as.factor(df.strain$sample_id)
df.strain$subject_id=as.factor(df.strain$subject_id)



##plot binary relation in df.strain

prev.sp=as.data.frame(table(df.strain$species_id))
prev.sp.ord=prev.sp[order(prev.sp$Freq),]
to.drop=prev.sp.ord[prev.sp.ord$Freq<4,'Var1']
##drop species with not enough observations
df.all3=df.strain[!df.strain$species_id %in% to.drop,]



length(unique(df.all3$species_id))  #134
length(unique(df.strain$species_id))  ##194

library(dplyr)


corfun<-function(x, y) {
  corr=(cor.test(x, y, method="pearson"))
}
length(unique(df.all3$species_id)) #134

stat.shan1=ddply(df.all3, .(species_id), summarise,
                 pval=corfun(shannon_diversity,strain_nb)$p.value,
                 est=corfun(shannon_diversity,strain_nb)$estimate) 

stat.shan1=stat.shan1[!is.na(stat.shan1$pval),]
stat.shan1[stat.shan1$pval<0.05,]
sig.pval.strain1=stat.shan1[stat.shan1$pval<0.05,'species_id']
length(sig.pval.strain1) #23
#sig.pval=sig.pval[!is.na(sig.pval)]
to.plot=df.all3[df.all3$species_id %in% sig.pval.strain1 ,]
length(unique(to.plot$species_id)) ##23

###plot significant

##plot significant richness

Plot3 <- ggplot(to.plot,aes(x=shannon_diversity,y=strain_nb))+
  #Plot2.rich <- ggplot(mydata,aes(x=Shannon_alpha_diversity,y=polymorphism_rate))+
  geom_point(color='blue2',size=0.5) + 
  xlab("Shannon") + 
  ylab("Strain number") + 
  geom_smooth(method="lm",se=FALSE,fullrange=T,color='#FF0000',size=1)+
  #geom_smooth(method="loess",se=FALSE,fullrange=T,color='green2',size=0.4)+
  
  theme_bw()+
  theme(axis.text=element_text(size=4),
        axis.title=element_text(size=7), #element_blank(),
        strip.text = element_text(size=4),
        panel.spacing = unit(0.025, "cm"),
        strip.background = element_blank(), strip.placement = "outside")+
  facet_wrap( ~ species_name, scales = "free",ncol=4) +
  # stat_cor(size=1.4,label.x.npc = "left",
  #          label.y.npc = "top",method='spearman') 
  
  stat_cor(size=1.4,label.x.npc = "left",
           label.y.npc = "top",method='spearman',r.accuracy=0.01,p.accuracy=0.001)+
  #panel.background = element_rect(fill = 'lightblue', color = 'purple'))+
  
  scale_y_continuous(trans=log10_trans())
Plot3

#ggsave('sig_pi_shannon_hmp_v4.pdf',width=5,height = 5 )  #with log10

ggsave('sig_pi_richness_hmp_v4.pdf',width=5,height = 5.5 )  #with log10



##########
###########
############
##################
stat.shan1[stat.shan1$pval<0.05 & stat.shan1$est>0, 'species_id'] ##21
stat.shan1[stat.shan1$pval<0.05 & stat.shan1$est<0, 'species_id'] ##
###[1] Bacteroides_uniformis_57318 Ruminococcus_gnavus_57638 

to.plot$species_id=as.factor(to.plot$species_id)

sig.pval.strain1[sig.pval.strain1 %in% sig.pval.rich]
# [1] Alistipes_finegoldii_56071         Alistipes_onderdonkii_55464       
# [3] Bacteroides_cellulosilyticus_58046 Bacteroides_fragilis_54507        
# [5] Bacteroides_uniformis_57318        Bacteroides_vulgatus_57955        
# [7] Burkholderiales_bacterium_56577    Parabacteroides_merdae_56972  

sig.pval.strain1[sig.pval.strain1 %in% sig.pval.pi] 
# [1] Alistipes_onderdonkii_55464        Bacteroides_fragilis_54507        
# [3] Bacteroides_uniformis_57318        Bacteroides_vulgatus_57955        
# [5] Burkholderiales_bacterium_56577    Faecalibacterium_prausnitzii_62201
# [7] Oscillibacter_sp_60799             Parabacteroides_merdae_56972 

##les negatives
sig.pval.rich[sig.pval.rich=='Ruminococcus_gnavus_57638'] ##non
sig.pval.pi[sig.pval.pi=='Ruminococcus_gnavus_57638']  ##0
##Bacteroides_uniformis_57318


nb.cols <- 23
myColors.strain.shan <- colorRampPalette(brewer.pal(4, "Spectral"))(nb.cols)

#myColors <- brewer.pal(25,"Spectral")
names(myColors.strain.shan) <- sig.pval.strain1

for (i in 1:length(myColors.strain.shan))
{
  name_i=names(myColors.strain.shan)[i]
  if (name_i %in% sig.pval.rich)
  {
     myColors.strain.shan[names(myColors.strain.shan)==name_i]=myColors.pi.rich[names(myColors.pi.rich)==name_i]
  }
     
}

for (i in 1:length(myColors.strain.shan))
{
  name_i=names(myColors.strain.shan)[i]
  if (name_i %in% sig.pval.pi)
  {
    myColors.strain.shan[names(myColors.strain.shan)==name_i]=myColors.pi.shan[names(myColors.pi.shan)==name_i]
  }
  
}

#myColors.strain.shan[names(myColors) sig.pval.strain1[sig.pval.strain1%in%sig.pval.rich]

names(myColors.strain.shan)

myColors.strain.shan[names(myColors.strain.shan)=='Ruminococcus_gnavus_57638']## "#2B83BA" (bleu)
myColors.strain.shan[names(myColors.strain.shan)=='Bacteroides_uniformis_57318']## "#2B83BA" (brique)


colScale.strain.shan <- scale_colour_manual(name = "sig.pval.strain1",values = myColors.strain.shan)

sp.shan <- ggplot(to.plot,aes(x=shannon_diversity,y=strain_nb,color=species_id))+
  geom_point(size=1,show.legend = F) + xlab("Shannon") + ylab("Strains count") + 
  guides( size = FALSE)+
  geom_smooth(data=to.plot,method=lm, formula=y~x,se=FALSE,size=0.5) +
  colScale.strain.shan+
  theme_bw()+
  
  theme(legend.position = 'right',
        legend.title = element_blank(),
        legend.text = element_text (size=3.5),
        legend.direction= 'vertical',
        legend.key.size = unit(0.1, "cm"),
        axis.text=element_text(size=4),
        axis.title=element_text(size=7))
        
        
#guides(fill=guide_legend(ncol=1))
#scale_fill_continuous(guide = guide_legend(ncol=1))

sp.shan1=sp.shan+guides(col=guide_legend(ncol=1),label.vjust=0.5)
sp.shan1
sp.shan1.1=sp.shan1+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

ggsave('shan_strain_significant_correl_point_v1.pdf',width=7,height = 5)

##
library(dplyr)

stat.rich1=ddply(df.all3, .(species_id), summarise,
                 pval=corfun(species_richness_raref,strain_nb)$p.value,
                 est=corfun(species_richness_raref,strain_nb)$estimate) 
stat.rich1=stat.rich1[!is.na(stat.rich1$pval),]

pos.rich=stat.rich1[stat.rich1$est>0 & stat.rich1$pval < 0.05, ]  ##3
dim(stat.rich1[stat.rich1$pval < 0.05, ]) #16

sig.pval.strain.rich=stat.rich1[stat.rich1$pval<0.05,'species_id']

to.plot=df.all3[df.all3$species_id %in% sig.pval.strain.rich ,]
length(unique(to.plot$species_id)) ##16


to.plot$species_id=as.factor(to.plot$species_id)

length(unique(to.plot$species_id))

stat.rich1[stat.rich1$pval<0.05 & stat.rich1$est>0 ,'species_id']
#[1] Alistipes_putredinis_61533     Bacteroides_clarus_62282       Phascolarctobacterium_sp_59818

#library(RColorBrewer)
nb.cols <- 16
myColors.strain.rich <- colorRampPalette(brewer.pal(8, "Spectral"))(nb.cols)

names(myColors.strain.rich) <- sig.pval.strain.rich


for (i in 1:length(myColors.strain.rich))
{
  name_i=names(myColors.strain.rich)[i]
  if (name_i %in% sig.pval.rich)
  {
    myColors.strain.rich[names(myColors.strain.rich)==name_i]=myColors.pi.rich[names(myColors.pi.rich)==name_i]
  }
  
}

for (i in 1:length(myColors.strain.rich))
{
  name_i=names(myColors.strain.rich)[i]
  if (name_i %in% sig.pval.pi)
  {
    myColors.strain.rich[names(myColors.strain.rich)==name_i]=myColors.pi.shan[names(myColors.pi.shan)==name_i]
  }
  
}

#colScale.strain.shan <- scale_colour_manual(name = "sig.pval.strain1",values = myColors.strain.shan)

for (i in 1:length(myColors.strain.rich))
{
  name_i=names(myColors.strain.rich)[i]
  if (name_i %in% sig.pval.strain1)
  {
    myColors.strain.rich[names(myColors.strain.rich)==name_i]=myColors.strain.shan[names(myColors.strain.shan)==name_i]
  }
  
}


colScale.strain.rich <- scale_colour_manual(name = "sig.pval.strain.rich",values = myColors.strain.rich)


sp.rich1 <- ggplot(to.plot,aes(x=species_richness_raref,y=strain_nb,color=species_id))+
  geom_point(size=1,show.legend = F) + xlab("Species richness")  + 
  guides( size = FALSE)+
  geom_smooth(data=to.plot,method=lm, formula=y~x,se=FALSE,size=0.5) +
  colScale.strain.rich+
  theme_bw()+
  theme(legend.position = 'right',
        legend.title = element_blank(),
        legend.text = element_text (size=3.5),
        legend.direction= 'vertical',
        legend.key.size = unit(0.1, "cm"),
        axis.text=element_text(size=4),
        axis.title=element_text(size=7),
        
        axis.title.y=element_blank())
sp.rich1


#ggsave('rich_strain_significant_correl_point_v1.pdf',width=7,height = 5)


library(gridExtra)

ggarrange(sp.shan1, sp.rich1, widths = c(2,2),labels = c("A", "B"))
ggsave('rich_shan_strain_significant_correl_point_v1.pdf',width=7,height=3)



sp.shan2.1=sp.shan2+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
sp.rich2.1=sp.rich2+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
sp.shan1.1=sp.shan1+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
sp.rich1.1=sp.rich1+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

library(cowplot)

pdf('figure2.pdf',width=6.5,heigh=6) 
plot_grid(sp.shan2.1, sp.rich2.1,sp.shan1.1, sp.rich1.1,labels = c('A', 'B','C','D'),label_size=7)
dev.off()
