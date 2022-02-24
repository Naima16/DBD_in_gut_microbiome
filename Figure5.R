library(dplyr)
library(plyr)
library(ggpubr)
pathway.df=read.table('sample_species_pathway_richness.csv',sep=',',header=T)
colnames(pathway.df)
# [1] "sample_id"                      "species"                        "pathway_richness"              
# [4] "community_pathway_richness"     "community_pathway_richness_set"
dim(pathway.df) #23301     5

df.path.cl=pathway.df[pathway.df$species != 'unclassified',] ##469  
dim(df.path.cl) #22832     5

sample_data=read.table('sample_covariates.txt',sep='\t',header=T)
colnames(sample_data)


df.path=merge(df.path.cl,sample_data,by='sample_id')

df.path$species=as.factor(df.path$species)  
df.path$sample_id=as.factor(df.path$sample_id)
df.path$subject_id=as.factor(df.path$subject_id)


library(plyr)
corfun<-function(x, y) {
  corr=(cor.test(x, y, method="spearman"))
}

 
prev.sp=as.data.frame(table(df.path$species))
prev.sp.ord=prev.sp[order(prev.sp$Freq),]
to.drop=prev.sp.ord[prev.sp.ord$Freq<4,'Var1']
length(to.drop) #130
length(unique(df.path$species)) #368

# ##drop species with not enough observations
df.all3.path=df.path[!df.path$species %in% to.drop,]

stat.shan1.path=ddply(df.all3.path, .(species), summarise,
                 pval=corfun(community_pathway_richness_set,pathway_richness)$p.value,
                 est=corfun(community_pathway_richness_set,pathway_richness)$estimate)
stat.shan1.path=stat.shan1.path[!is.na(stat.shan1.path$pval),]
dim(stat.shan1.path[stat.shan1.path$pval<0.05,]) #107
dim(stat.shan1.path[stat.shan1.path$pval<0.05 & stat.shan1.path$est<0,]) #95
dim(stat.shan1.path[stat.shan1.path$pval<0.05 & stat.shan1.path$est>0,]) #12

#
sig.shan.path=stat.shan1.path[stat.shan1.path$pval<0.05,]
sig.shan.path=pos.shan.path[!is.na(pos.shan.path$pval),]
# 
to.plot.path=df.all3.path[df.all3.path$species %in% sig.shan.path$species  ,]
length(unique(to.plot.path$species)) #107


library(ggplot2)

library(RColorBrewer)

nb.cols <- 107
#myColors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
myColors <- colorRampPalette(brewer.pal(10, "Spectral"))(nb.cols)

#myColors <- brewer.pal(25,"polychrome")
ps.corel=sig.shan.path$species
ps.corel=ps.corel[!is.na(ps.corel)]
names(myColors) <- ps.corel


myColors[names(myColors)=='g__Escherichia.s__Escherichia_coli']='#00A08A'
myColors[names(myColors)=='g__Enterobacter.s__Enterobacter_cloacae']="#F2AD00"
myColors[names(myColors)=='g__Klebsiella.s__Klebsiella_pneumoniae']='#FF0000'
#myColors[names(myColors)=='g__Haemophilus.s__Haemophilus_parainfluenzae']='#F2AD00'

colScale <- scale_colour_manual(name = 'sp.neg.1',values = myColors)

sp.shan.neg <- ggplot(to.plot.path,aes(x=community_pathway_richness_set,y=pathway_richness,color=species))+
  geom_point(aes(color=species),size=0.5)  +
   colScale+
  ylab("Number of pathways in a focal species") + 
  xlab("Community pathway richness") + 
  #guides( size = FALSE)+
  geom_smooth(data=to.plot.path,method=lm, formula=y~x,se=FALSE,size=1) +

  theme_bw()+
  theme(legend.position = 'none',
        axis.text=element_text(size=4),
        axis.title=element_text(size=7), 
        strip.text = element_text(size=4))
sp.shan.neg
ggsave('plot4_pathway.pdf',width = 5,heigh=7)




##### genes 
gen.community=read.table('gene_counts_combined.csv',header=T,sep=',')
dim(gen.community) #13105    11
colnames(gen.community)
# [1] "species_name"               "sample_id"                  "pres_genes"                
# [4] "total_genes"                "community_pres_genes"       "pres_acc95_genes"          
# [7] "total_acc95_genes"          "community_pres_acc95_genes" "pres_acc90_genes"          
# [10] "total_acc90_genes"          "community_pres_acc90_genes"

unique(gen.community$species_name)


prev.sp=as.data.frame(table(gen.community$species_name))
prev.sp.ord=prev.sp[order(prev.sp$Freq),]
to.drop=prev.sp.ord[prev.sp.ord$Freq<4,'Var1']
length(to.drop) #44
length(unique(gen.community$species_name)) #151
##drop species with not enough observations
df.all3.genes=gen.community[!gen.community$species_name %in% to.drop,]
length(unique(df.all3.genes$species_name))  #114


library(dplyr)


##pres genes
prev.sp=as.data.frame(table(gen.community$species_name))
length(unique(gen.community$species_name))  #69
prev.sp.ord=prev.sp[order(prev.sp$Freq),]
to.drop=prev.sp.ord[prev.sp.ord$Freq<4,'Var1']  #1
##drop species with not enough observations
df.all3.genes2=gen.community[!gen.community$species_name %in% to.drop,]


stat.shan1.genes2=ddply(df.all3.genes2, .(species_name), summarise,
                 pval=corfun(pres_genes,community_pres_genes)$p.value,
                 est=corfun(pres_genes,community_pres_genes)$estimate) 

sig.pval.genes2=stat.shan1.genes2[stat.shan1.genes2$pval<0.05,'species_name']
length(sig.pval.genes2) ##36

#stat.shan1.genes2[stat.shan1.genes2$pval<0.05 & stat.shan1.genes2$est<0 ,'species_name'] ##0
######[1] Alistipes_sp_59510       Eubacterium_hallii_61477


to.plot.shan.genes2=df.all3.genes2[df.all3.genes2$species_name%in% sig.pval.genes2 ,]
length(unique(to.plot.shan.genes2$species_name)) #36


library(ggplot2)

library(RColorBrewer)
nb.cols <- 36
myColors <- colorRampPalette(brewer.pal(4, "Spectral"))(nb.cols)

names(myColors) <- unique(to.plot.shan.genes2$species_name)

colScale <- scale_colour_manual(name = "sig.pval",values = myColors)

qbh1 <- ggplot(to.plot.shan.genes2,aes(x=community_pres_genes,y=pres_genes,color=species_name))+
  geom_point(size=0.5,show.legend = F) +
  guides( size = FALSE)+
  geom_smooth(data=to.plot.shan.genes2,method=lm, formula=y~x,se=FALSE,size=1) +
  xlab("Community gene richness") + 
  ylab("Number of genes in a focal species") + 
   colScale+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text=element_text(size=4),
        axis.title=element_text(size=7), 
        strip.text = element_text(size=4))
qbh1


library(cowplot)

pdf('Figure5.pdf',width=5,height=3.5,pointsize = 7)
#grid.arrange(grobs = list(qbh1,acc.plot2,pathway),layout_matrix = lay1)
plot_grid(qbh1,sp.shan.neg,ncol=2,labels=c('A','B'),label_size = 8,rel_widths=c(1.2,1.4))
dev.off()

