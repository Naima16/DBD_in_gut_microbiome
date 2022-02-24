
genes=read.table('HMP_gene_changes_full.csv',sep=',',header=T)
dim(genes)  #1339   14
colnames(genes)


# [1] "sample_tp1"         "sample_tp2"         "subject"            "tp1"                "tp2"               
# [6] "species"            "num_gene_gains"     "num_gene_losses"    "alpha_div_tp1"      "alpha_div_tp2"     
# [11] "alpha_div_rare_tp1" "alpha_div_rare_tp2" "richness_tp1"       "richness_tp2"       "richness_rare_tp1" 
# [16] "richness_rare_tp2"  "polymorphism_tp1"   "polymorphism_tp2"  
#   

length(unique(genes$sample_tp1)) #210
length(unique(genes$sample_tp2)) #211
length(unique(genes$subject))  #154
length(unique(genes$species))   #54

colnames(genes)
list_sp=c()
spec_names=unique(genes$species)
for (i in spec_names)
{
  cc=genes[genes$species==i,]
  if (all(cc$num_genes_changes == 0))
    list_sp <- c(list_sp, i)
  
}

# cc=genes[genes$qp_status=='False',]
# > length(unique(cc$species))
# [1] 33
# > length(unique(cc$sample_tp1))
# [1] 164

hist(genes$num_gene_gains)
range(genes$num_gene_gains)  #0 425



###########################################
######################################
###plot only significant correltions
prev.sp=as.data.frame(table(genes$species))
length(unique(genes$species))  #54
prev.sp.ord=prev.sp[order(prev.sp$Freq),]
to.drop=prev.sp.ord[prev.sp.ord$Freq<4,'Var1']  #10

##drop species with not enough observations
df.all3=genes[!genes$species %in% to.drop,]
dim(df.all3)
dim(genes)

sum_loss=df.all3 %>%                                        # Specify data frame
  group_by(species) %>%                         # Specify group indicator
  summarise_at(vars(num_gene_losses),              # Specify column
               list(los = sum)) 

sum_loss_sp=sum_loss[sum_loss$los==0,'species']
sum_loss_sp1=as.list(sum_loss_sp)
genes.3.los=df.all3[!df.all3$species %in% sum_loss_sp$species, ]
dim(genes.3.los)
length(unique(genes.3.los$species)) #33
length(unique(df.all3$species)) #44
dim(df.all3)


length(unique(df.all3$species))  #44
length(unique(genes$species))  ##54

 library(dplyr)

corfun<-function(x, y) {
  corr=(cor.test(x, y, method="pearson"))
}

length(unique(df.all3$species))
#[1] 68

stat.shan1=ddply(genes.3.los, .(species), summarise,
                 pval=corfun(alpha_div_tp1,num_gene_losses)$p.value,
                 est=corfun(alpha_div_tp1,num_gene_losses)$estimate) 
##NA are species with only 0

sig.pval.pi=stat.shan1[stat.shan1$pval<=0.05,'species']
#[1] Bacteroides_uniformis_57318 Bacteroides_vulgatus_57955  Dialister_invisus_61905  

stat.shan1[stat.shan1$pval<0.05 & stat.shan1$est<0 ,'species']
#[1] Dialister_invisus_61905
#sig.pval.pi=sig.pval.pi[!is.na(sig.pval.pi)]
stat.shan1[stat.shan1$pval<0.05 & stat.shan1$est>0 ,] #2

to.plot.shan=df.all3[df.all3$species%in% sig.pval.pi ,]
length(unique(to.plot.shan$species)) #11


library(ggplot2)

library(RColorBrewer)
nb.cols <- 33
myColors <- colorRampPalette(brewer.pal(4, "Spectral"))(nb.cols)

names(myColors) <- unique(genes.3.los$species)

myColors[names(myColors)=='Bacteroides_uniformis_57318']='cyan'
myColors[names(myColors)=='Bacteroides_vulgatus_57955']='#FF0000'
myColors[names(myColors)=='Dialister_invisus_61905']='black'
myColors[names(myColors)=='Alistipes_putredinis_61533']='magenta'
myColors[names(myColors)=='Bacteroides_stercoris_56735']='blue'
myColors[names(myColors)=='Parabacteroides_merdae_56972']='green'


colScale <- scale_colour_manual(name = "sig.pval.pi",values = myColors)


sp.shan_p <- ggplot(genes.3.los,aes(x=alpha_div_tp1,y=log10(num_gene_losses+1),color=species))+
  geom_point(size=1.5,show.legend = F) + xlab("Shannon diversity") + ylab("Log number of genes lost in a focal species") + 
  guides( size = FALSE)+
  #geom_smooth(data=genes.3.los,method=lm, formula=y~x,se=FALSE,size=0.8) +
 
  geom_line(stat="smooth",data=genes.3.los,method = 'lm', formula = y~x,
            size = 0.8,se=FALSE,
            alpha = 0.5)+
  annotate("text", x = 1.5, y = 2.9, label = "B.uniformis",size=1.4)+
  annotate("text", x = 1.5, y = 2.8, label = "B.vulgatus",size=1.4)+
  annotate("text", x = 1.539, y = 2.7, label = "Dialister.invisus",size=1.4)+
  annotate("text", x = 1.47, y = 3,  label = 'Significant',size=1.4,fontface = "bold")+
 
  colScale+
  theme_bw()+
  
  theme(
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text (size=4),
        legend.direction= 'vertical',
        #legend.spacing.y = unit(0.01, 'cm'),
        legend.key.size = unit(0.1, "cm"),
        axis.text=element_text(size=5),
        axis.title=element_text(size=7)) 

p1=sp.shan_p +guides(col=guide_legend(ncol=1),label.vjust=0.5)
p1

ggsave('geneslos_los.pdf',width=7,height = 5)


stat.rich1=ddply(genes.3.los, .(species), summarise,
                 pval=corfun(richness_rare_tp1,num_gene_losses)$p.value,
                 est=corfun(richness_rare_tp1,num_gene_losses)$estimate) 
##NA are species with only 0

sig.rich=stat.rich1[stat.rich1$pval<=0.05,'species']

stat.rich1[stat.rich1$pval<0.05 & stat.rich1$est<0 ,'species']
#0

stat.rich1[stat.rich1$pval<0.05 & stat.shan1$est>0 ,] #all 5

to.plot.rch=df.all3[df.all3$species%in% sig.rich ,]
length(unique(to.plot.rch$species)) #5


library(ggplot2)


sp.shan <- ggplot(genes.3.los,aes(x=richness_rare_tp1,y=log10(num_gene_losses+1),color=species))+
  geom_point(size=1.5,show.legend = F) + xlab("Species richness")  + 
  guides( size = FALSE)+
  #geom_smooth(data=genes.3.los,method=lm, formula=y~x,se=FALSE,size=0.8,alpha=0.5) +
  geom_line(stat="smooth",data=genes.3.los,method = 'lm', formula = y~x,
            size = 0.8,se=FALSE,
            alpha = 0.5)+
  

  colScale+
  theme_bw()+
  annotate("text", x = 52, y = 2.9, label = "B.uniformis",size=1.4)+
  annotate("text", x = 52, y = 2.8, label = "B.vulgatus",size=1.4)+
  annotate("text", x = 52, y = 2.7, label = "B.stercoris",size=1.4)+
  annotate("text", x = 59, y = 2.6, label = "Alistipes.putredinis",size=1.4)+
  annotate("text", x = 64, y = 2.5, label = "Parabacteroides.merdae",size=1.4)+
  annotate("text", x = 50, y = 3,  label = 'Significant',size=1.4,fontface = "bold")+
  
  theme(
    legend.position = 'right',
    legend.title = element_blank(),
    legend.text = element_text (size=4),
    #legend.direction= 'horizontal',
    #legend.spacing.y = unit(0.01, 'cm'),
    legend.key.size = unit(0.1, "cm"),
    axis.text=element_text(size=5),
    axis.title=element_text(size=7),
    axis.title.y=element_blank()) 

sp.shan3=sp.shan+guides(col=guide_legend(ncol=1),label.vjust=0.5)
sp.shan3

ggsave('geneslos_rich.pdf',width=7,height = 5)

# 
# 
p1.1=p1+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
sp.shan3.1=sp.shan3+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

library(cowplot)
pdf('fig3.pdf',width=6.5,heigh=2.7) ##7,3
plot_grid(p1, sp.shan3,labels = c('A', 'B'),label_size=7,align='h',rel_widths=c(1,1.5))
dev.off()


###10 most prevalent
prev.sp.ord=prev.sp[order(-prev.sp$Freq),]
prev.sp.ord[1:10,]

library(ggplot2)

###10 species
#c.ord=as.data.frame(c[with(c,order(-Freq)),])
mydata1=genes.3.los[which(genes.3.los$species %in% prev.sp.ord[1:10,'Var1']),]
dim(mydata1)

rich <- ggplot(mydata1,aes(x=richness_rare_tp1,y=log10(num_gene_losses+1)))+
  geom_point(color='blue2',size=0.9) + 
  xlab("Species richness") + 
  ylab("Log gene loss") + 
  geom_smooth(method="lm",se=FALSE,fullrange=T,color='red',size=0.8)+
  #geom_smooth(method="loess",se=FALSE,fullrange=T,color='green2',size=0.4)+
  
  theme_bw()+
  theme(axis.text=element_text(size=4),
        axis.title=element_text(size=7), #element_blank(),
        strip.text = element_text(size=3.5),
        panel.spacing = unit(0.025, "cm"),
        strip.background = element_blank(), strip.placement = "outside")+
  facet_wrap( ~ species, scales = "free",ncol=2) 
  #stat_cor(data=mydata1,aes(x=richness_rare_tp1,y=num_gene_losses),size=1.4,label.x.npc = "left",
      #     label.y.npc = "top") 
rich

shan <- ggplot(mydata1,aes(x=alpha_div_tp1,y=log10(num_gene_losses+1)))+
  geom_point(color='blue2',size=0.9) + 
  xlab("Shannon") + 
  ylab("Log gene loss") + 
  geom_smooth(method="lm",se=FALSE,fullrange=T,color='red',size=0.8)+
  #geom_smooth(method="loess",se=FALSE,fullrange=T,color='green2',size=0.4)+
  
  theme_bw()+
  theme(axis.text=element_text(size=4),
        axis.title=element_text(size=7), #element_blank(),
        strip.text = element_text(size=3.5),
        panel.spacing = unit(0.025, "cm"),
        strip.background = element_blank(), strip.placement = "outside")+
  facet_wrap( ~ species, scales = "free",ncol=2) 
 # stat_cor(data=mydata1,size=1.4,label.x.npc = "left",
       #    label.y.npc = "top") 
shan

pdf('gene_loss_top10.pdf',width=4 ,height = 4,pointsize = 6)
#### 
grid.arrange(shan,rich,ncol=2)
dev.off()

