
pi_poyet=read.table('/Users/naima/projet_HMP_gut_nandita_sept20/5juillet_final_version/december_2021_Longitudinal_data/Poyet/Poyet_delta_polymorphism.csv',sep=',',header=T)
colnames(pi_poyet)
# [1] "sample_tp1"         "sample_tp2"         "subject"            "tp1"                "tp2"               
# [6] "species"            "delta_polymorphism" "polymorphism_tp1"   "polymorphism_tp2"   "alpha_div_tp1"     
# [11] "alpha_div_tp2"      "alpha_div_rare_tp1" "alpha_div_rare_tp2" "richness_tp1"       "richness_tp2"      
# [16] "richness_rare_tp1"  "richness_rare_tp2"
dim(pi_poyet)

pi_am=pi_poyet[pi_poyet$subject=='am' & pi_poyet$species=='Bacteroides_vulgatus_57955',]
length(unique(pi_am$species))  #1
length(unique(pi_am$sample_tp1))   ##204


range(pi_am$tp1) #1:230

###am data de moi
am_access_date=read.table('/Users/naima/Projet2_zhao/S01_am_analyses/unrarefied_am_24janv/11janv/am_accession_date.csv',sep=",",header=T)
dim(am_access_date)  ##206-6
head(am_access_date)

library(schoolmath)
am_access_date[is.decimal(am_access_date$jour),]

for (i in 1:dim(pi_am)[1])
{
  sampl_name=as.character(pi_am[i,'sample_tp1'])
  pi_am[i,'days']=am_access_date[am_access_date$Run==sampl_name,'jour']
  pi_am[i,'collec_date']=am_access_date[am_access_date$Run==sampl_name,'collec_date']
}

###


range(pi_am$delta_polymorphism)
#-0.005272594  0.006089037


#########
cor.snp=cor.test(pi_am$delta_polymorphism,pi_am$alpha_div_tp1)


library(ggplot2)


range(pi_am$delta_polymorphism)

cor_pi=cor.test(pi_am$alpha_div_tp1,pi_am$delta_polymorphism,method="pearson")
cor_pi

p.1 <- ggplot(pi_am, aes(alpha_div_tp1,delta_polymorphism))
shanon_pi <- p.1 +
  geom_point(size=0.6,color="#FF0000") + ##blue2
  theme_bw()+
  xlab('Shannon diversity tp1')+ylab('Delta polymorphism')+
  #ggtitle("Shannon ~ SNV count")+
  theme(plot.title = element_text(size = 8),
        axis.title=element_text(size=7),
        #axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1,size =5),
        axis.text.y = element_text(size=5))+

  stat_smooth(method = "lm",color ='black',se=F,size=0.8) + ##FF0000
  annotate(geom="text", x = 1.8, y = 0.006, 
           label=paste(" \nr =",signif(cor_pi$estimate, 2),
                       " \nPval =",signif(cor_pi$p.value, 3)), colour="black",
           size=2, fontface="plain",hjust = 0)  #, family="Courier"

shanon_pi
ggsave('shanon_deltapi.pdf')

cor_pi2=cor.test(pi_am$richness_rare_tp1,pi_am$delta_polymorphism,method='pearson')
cor_pi2

pi.rich <- ggplot(pi_am, aes(richness_rare_tp1,delta_polymorphism))
rich_pi <- pi.rich + geom_point(size=0.6,color="#FF0000")  +
  theme_bw()+
  xlab('Species richness tp1')+ylab('Delta polymorphism ')+
  #ggtitle("Richness ~ SNV count")+
  theme(plot.title = element_text(size = 8),
        axis.title=element_text(size=7),
        #axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1,size =5),
        axis.text.y = element_text(size=5))+
  stat_smooth(method = "lm",color ='black',se=F,size=0.8) +
  annotate(geom="text", x = 55, y = 0.006, 
           label=paste(" \nr =",signif(cor_pi2$estimate, 2),
                       " \nPval =",signif(cor_pi2$p.value, 3)), colour="black",
           size=2, fontface="plain",hjust = 0)  #, family="Courier"
rich_pi
ggsave('richness_deltapi.pdf')

pi_am$sample_tp1=as.factor(pi_am$sample_tp1)
hist(pi_am$delta_polymorphism)

am.pi.1=lm(delta_polymorphism~alpha_div_tp1,data=pi_am)
summary(am.pi.1)


###genes
gene_poyet=read.table('/Users/naima/projet_HMP_gut_nandita_sept20/5juillet_final_version/december_2021_Longitudinal_data/Poyet/Poyet_gene_changes_full.csv',sep=',',header=T)
#gene_poyet=read.table('/Users/naima/projet_HMP_gut_nandita_sept20/5juillet_final_version/longitudinal_data_12nov/poyrt/Poyet_gene_changes_full.csv',sep=',',header=T)
dim(gene_poyet)

gene_am=gene_poyet[gene_poyet$subject=='am' & gene_poyet$species=='Bacteroides_vulgatus_57955',]


length(unique(gene_am$species))  #1
length(unique(gene_am$sample_tp1))   ##188

for (i in 1:dim(gene_am)[1])
{
  sampl_name=as.character(gene_am[i,'sample_tp1'])
  gene_am[i,'days']=am_access_date[am_access_date$Run==sampl_name,'jour']
  gene_am[i,'collec_date']=am_access_date[am_access_date$Run==sampl_name,'collec_date']
}


cor_var=cor.test(gene_am$alpha_div_tp1,gene_am$num_gene_gains)

p.var <- ggplot(gene_am, aes(alpha_div_tp1,num_gene_gains))
gain.p <- p.var + geom_point(size=1,color="#FF0000")  +
  stat_smooth(method = "lm",color ='black',se=F,size=0.8) +
  xlab('Shannon diversity tp1')+
  ylab('Number of genes gained')+
  theme_bw()+
  annotate(geom="text", x = 1.8, y = 1.8, 
           label=paste(" \nr =",signif(cor_var$estimate, 2),
                       " \nPval =",signif(cor_var$p.value, 3)), colour="black",
           size=2, fontface="plain",hjust = 0) + #, family="Courier"
        theme(plot.title = element_text(size = 8),
        axis.title=element_text(size=7),
        #axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1,size =5),
        axis.text.y = element_text(size=5))
  


gain.p
ggsave('shannon_gain.pdf')

cor_var2=cor.test(gene_am$richness_rare_tp1,gene_am$num_gene_gains)

p.var2 <- ggplot(gene_am, aes(richness_rare_tp1,num_gene_gains))
gain.p2 <- p.var2 + geom_point(size=1,color="#FF0000") +
  xlab('Species richness tp1')+ylab('Number of genes gained')+
  # ggtitle("Species richness ~ gene gain")+
  theme(plot.title = element_text(size = 10),
        axis.title=element_text(size=10))+
  stat_smooth(method = "lm",color ='black',se=F,size=0.8) +
  theme_bw()+
  
  annotate(geom="text", x = 58, y = 1.8, 
           label=paste(" \nr =",signif(cor_var2$estimate, 2),
                       " \nPval =",signif(cor_var2$p.value, 3)), colour="black",
           size=2, fontface="plain",hjust = 0) +
  theme(plot.title = element_text(size = 8),
        axis.title=element_text(size=7),
        #axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1,size =5),
        axis.text.y = element_text(size=5))
                
gain.p2
ggsave('richness_gain.pdf')



##loss gene



cor_var=cor.test(gene_am$alpha_div_tp1,gene_am$num_gene_losses)

p.var <- ggplot(gene_am, aes(alpha_div_tp1,num_gene_losses))
loss.p1 <- p.var +  geom_point(size=1,color="#FF0000")  +
  theme_bw()+
  xlab('Shannon diversity tp1')+ylab('Number of genes lost')+
  theme(plot.title = element_text(size = 8),
        axis.title=element_text(size=7),
        #axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1,size =5),
        axis.text.y = element_text(size=5))+
  stat_smooth(method = "lm",color ='black',se=F,size=0.8) +
  annotate(geom="text", x = 1.8, y = 3.5, label=paste(" \nr =",signif(cor_var$estimate, 2),
                                                      " \nPval =",signif(cor_var$p.value, 3)), colour="black",
           size=2, fontface="plain",hjust = 0) 
  
loss.p1
ggsave('shannon_loss.pdf')

cor_var2=cor.test(gene_am$richness_rare_tp1,gene_am$num_gene_losses)

#gene_am$collec_date
p.var2 <- ggplot(gene_am, aes(richness_rare_tp1,num_gene_losses))
loss.p2 <- p.var2 + geom_point(size=1,color="#FF0000")  +
  theme_bw()+
  xlab('Species richness tp1')+
  ylab('Number of genes lost')+
  theme(plot.title = element_text(size = 8),
        axis.title=element_text(size=7),
        #axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1,size =5),
        axis.text.y = element_text(size=5))+
  stat_smooth(method = "lm",color ='black',se=F,size=0.8) +
  annotate(geom="text", x=58, y=3.5, label=paste(" \nr =",signif(cor_var2$estimate, 2),
                                                 " \nPval =",signif(cor_var2$p.value, 3)), colour="black",
           size=2, fontface="plain",hjust = 0)  +
  theme(axis.title=element_text(size=8),
        axis.text.y=element_text(size=5),
        axis.text.x = element_text(hjust = 1,size=5))

loss.p2

ggsave('richness_loss.pdf')





p1=shanon_pi+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p2=loss.p1+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p3=gain.p+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p4=rich_pi+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p5=loss.p2+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p6=gain.p2+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))



pdf('Figure4.pdf',width=5.7,height=2.6,pointsize = 8)
plot_grid(shanon_pi,loss.p1,gain.p,
          ncol=3,nrow=1,labels = c('A','B','C'),label_size = 7,rel_widths = c(1.2,1.1,1.1))
dev.off()

pdf('FigureS11.pdf',width=5.7,height=2.6,pointsize = 8)
plot_grid(rich_pi,loss.p2,gain.p2,
          ncol=3,nrow=1,labels = c('A','B','C'),label_size = 7,rel_widths = c(1.2,1.1,1.1))
dev.off()
