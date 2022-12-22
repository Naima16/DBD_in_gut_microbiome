#### gam with polymorphism rate as a function of Shannon diversity
#### Naima Madi dec 2022


library(bbmle) ##AICtab
library(mgcv)  ##gam
library(tidyverse) ##reduce

### polymorphism rates at synonymous sites / species from MIdas
Garud.data=read.table('~/HMP1-2_polymorphism_rates_alpha_divs.csv',sep=',',header=T)
dim(Garud.data)
head(Garud.data)
colnames(Garud.data)
#[1] "sample_id"               "subject_id"              "species_name"            "polymorphism_rate"       "Shannon_alpha_diversity"

Garud.data$sample_id=as.factor(Garud.data$sample_id)
Garud.data$species_name=as.factor(Garud.data$species_name)
Garud.data$subject_id=as.factor(Garud.data$subject_id)


sample_data=read.table('~/sample_covariates.txt',sep='\t',header=T)
colnames(sample_data)
# [1] "sample_id"                 "subject_id"                "shannon_diversity"         "species_richness"          "total_reads_MIDAS"        
# [6] "total_reads_orig"          "total_avg_marker_coverage" "total_med_marker_coverage"

Garud.data1=merge(Garud.data,sample_data[,c('sample_id','shannon_diversity','species_richness','total_reads_orig')],by='sample_id')

Garud.data1$species_name=as.factor(Garud.data1$species_name)  
Garud.data1$sample_id=as.factor(Garud.data1$sample_id)
Garud.data1$subject_id=as.factor(Garud.data1$subject_id)



gam1.pi.fs <- gam(polymorphism_rate ~ s(Shannon_alpha_diversity)+ s(total_reads_orig)+
                 s(species_name,Shannon_alpha_diversity,bs="fs")+
                 s(sample_id,bs="re")+
                 s(subject_id,bs="re"),
                 data=Garud.data1,family=betar(link="logit"))
summary(gam1.pi.fs)
save(gam1.pi.fs,file='gam1.pi.fs.RData')


#### diagnostic plots
pdf('gam1.pi.fs_check.pdf')
par(mfrow=c(2,2))
gam.check(gam1.pi.fs,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()

##polymorphism rate as a function of richness on all data
gam.pi.rich.all <- gam(polymorphism_rate ~ s(species_richness)+s(total_reads_orig)+
                         s(species_name,species_richness,bs="fs")+
                         s(sample_id,bs="re")+s(subject_id,bs="re"),
                       data=Garud.data1,family=betar(link="logit"))
save(gam.pi.rich.all,file='gam.pi.rich.all.RData')

pdf('gam.pi.rich.all_check.pdf')
par(mfrow=c(2,2))
gam.check(gam.pi.rich.all ,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()

### sample data with rarefied richness
sample_data_1=read.table('~/sample_covariates_rarefied_20m.txt',sep='\t',header=T)
colnames(sample_data_1)
# [1] "sample_id"         "subject_id"        "shannon_diversity" "species_richness"  "total_reads_MIDAS"
# [6] "total_reads_orig"


# merge polymorphism with rarefied richness and total reads
sampl.total_reads=sample_data[,c('sample_id','total_reads_orig')]
sampl.raref=sample_data_1[,c('sample_id','species_richness')]

df_list <- list(Garud.data,sampl.raref, sampl.total_reads)      

#merge all data frames together
Garud.data.raref=df_list %>% reduce(inner_join, by='sample_id')  #full_join
colnames(Garud.data.raref)
# [1] "sample_id"               "subject_id"              "species_name"            "polymorphism_rate"       "Shannon_alpha_diversity"
# [6] "species_richness"        "total_reads_orig"

gam.pi.rarefRich.cov <- gam(polymorphism_rate ~ s(species_richness)+
                              s(total_reads_orig)+
                              s(species_name,species_richness,bs="fs")+ 
                              s(sample_id,bs="re")+s(subject_id,bs="re"),
                            data=Garud.data.raref,family=betar(link="logit"))

save(gam.pi.rarefRich.cov,file='gam.pi.rarefRich.cov.RData')

pdf('gam.pi.rarefRich.cov_check.pdf')
par(mfrow=c(2,2))
gam.check(gam.pi.rarefRich.cov ,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()

##############
##### gams with polyorphism rate at synonymous sites as a function of Shannon diversity at higher taxonomic levels

## taxonomic annotation from GTDBK pipeline
taxo=read.table('~/gtdbtk.bac120.summary.tsv',sep='\t',header=T)
dim(taxo) 
colnames(taxo)

taxo.df=taxo[,c('user_genome','classification')]

dim(taxo.df)  #1258    2

for (i in 1:dim(taxo.df)[1] ) {
  taxo.df[i,'genus']=unlist(strsplit(unlist(strsplit(as.character(taxo.df[i,]$classification),";"))[6],"__"))[2]
  taxo.df[i,'family']=unlist(strsplit(unlist(strsplit(as.character(taxo.df[i,]$classification),";"))[5],"__"))[2]
  taxo.df[i,'order']=unlist(strsplit(unlist(strsplit(as.character(taxo.df[i,]$classification),";"))[4],"__"))[2]
  taxo.df[i,'class']=unlist(strsplit(unlist(strsplit(as.character(taxo.df[i,]$classification),";"))[3],"__"))[2]
  taxo.df[i,'phylum']=unlist(strsplit(unlist(strsplit(as.character(taxo.df[i,]$classification),";"))[2],"__"))[2]
  taxo.df[i,'GTDBK_species']=unlist(strsplit(unlist(strsplit(as.character(taxo.df[i,]$classification),";"))[7],"__"))[2]
  taxo.df[i,'species_name']=unlist(strsplit(as.character(taxo.df[i,]$user_genome), "[.]"))[1]
}
dim(taxo.df)  #1258    9

taxo.df[is.na(taxo.df$GTDBK_species),]  #aucun

df.all=merge(Garud.data1,taxo.df[,c('species_name','GTDBK_species','genus','family','order','class','phylum')],by='species_name')

df.all$genus=as.factor(df.all$genus)
df.all$family=as.factor(df.all$family)
df.all$order=as.factor(df.all$order)
df.all$class=as.factor(df.all$class)
df.all$phylum=as.factor(df.all$phylum)
df.all$species_name=as.factor(df.all$species_name)
df.all$GTDBK_species=as.factor(df.all$GTDBK_species)
df.all$sample_id=as.factor(df.all$sample_id)


##abundances from midas
midas_abund=read.table('/Users/naima/projet_HMP_gut_nandita_sept20/31janvier_MidasAbundance_files/relative_abundance.txt',sep='\t',header=T,check.names = F)
dim(midas_abund)  #5952  470

rownames(midas_abund)=midas_abund$species_id
midas_abund=midas_abund[,-1]

## community data
commun_data=t(midas_abund)
dim(commun_data)

##convert to num
commun_data_n=data.matrix(as.data.frame(commun_data))  ##sample*species
dim(commun_data_n) #469 5952

#remove species present in 0 sample
commun_data_n.nonzero=commun_data_n[,-(which(colSums(commun_data_n)==0))]

library(vegan)
shannon_midas=diversity(commun_data_n.nonzero,index='shannon')  ## la même que quand j'enlève pas les 0
shannon_midas=diversity(commun_data_n,index='shannon')
head(shannon_midas)

shannon_midas=as.data.frame(shannon_midas)

library(reshape2)
midas_abund$species_name=rownames(midas_abund)

midas_abund_long=melt(midas_abund,id='species_name')
dim(midas_abund_long)  ##2791488       3
colnames(midas_abund_long)[2] = 'sample_id'
colnames(midas_abund_long)[3] = 'abundance'


midas_taxo=merge(midas_abund_long,taxo.df,by='species_name')
sample_list=unique(midas_abund_long$sample_id)
sample_tax_shannon=data.frame()

i=1

for (i in 1:length(sample_list)) {
  
  sample_name=as.character(sample_list[i])
  sample_tax_shannon[i,'sample_id']=sample_name
  
  my_data=midas_taxo[midas_taxo$sample_id==sample_name & midas_taxo$abundance !=0,]
  t1=unique(my_data$phylum)
  rich_phyla=data.frame()
  j=1
  for (k in t1)
  { 
    rich_phyla[j,'phylum_id']=as.character(k)
    rich_phyla[j,'rich']=sum(my_data[my_data$phylum==k,'abundance'])
    j=j+1
  }
  
  phylum_div=diversity(rich_phyla$rich,index='shannon')
  
  t1=unique(my_data$class)
  rich_class=data.frame()
  j=1
  for (k in t1)
  { 
    rich_class[j,'class_id']=as.character(k)
    rich_class[j,'rich']=sum(my_data[my_data$class==k,'abundance'])
    j=j+1
  }
  class_div=diversity(rich_class$rich,index='shannon')
  
  t1=unique(my_data$order)
  rich_order=data.frame()
  j=1
  for (k in t1)
  { 
    rich_order[j,'order_id']=as.character(k)
    rich_order[j,'rich']=sum(my_data[my_data$order==k,'abundance'])
    j=j+1
  }
  order_div=diversity(rich_order$rich,index='shannon')
  
  t1=unique(my_data$family)
  rich_family=data.frame()
  j=1
  for (k in t1)
  { 
    rich_family[j,'family_id']=as.character(k)
    rich_family[j,'rich']=sum(my_data[my_data$family==k,'abundance'])
    j=j+1
  }
  family_div=diversity(rich_family$rich,index='shannon')
  
  t1=unique(my_data$genus)
  rich_genus=data.frame()
  j=1
  for (k in t1)
  { 
    rich_genus[j,'genus_id']=as.character(k)
    rich_genus[j,'rich']=sum(my_data[my_data$genus==k,'abundance'])
    j=j+1
  }
  genus_div=diversity(rich_genus$rich,index='shannon')
  
  t1=unique(my_data$GTDBK_species)
  rich_GTDBK_species=data.frame()
  j=1
  for (k in t1)
  { 
    rich_GTDBK_species[j,'GTDBK_species']=as.character(k)
    rich_GTDBK_species[j,'rich']=sum(my_data[my_data$GTDBK_species==k,'abundance'])
    j=j+1
  }
  sp_div=diversity(rich_GTDBK_species$rich,index='shannon')
 
  t=as.numeric(my_data$abundance)  
  sp_midas_div=diversity(t,index='shannon')
  
  sample_tax_shannon[i,'phyla_shannon']=phylum_div
  sample_tax_shannon[i,'class_shannon']=class_div
  sample_tax_shannon[i,'order_shannon']=order_div
  sample_tax_shannon[i,'family_shannon']=family_div
  sample_tax_shannon[i,'genus_shannon']=genus_div
  sample_tax_shannon[i,'gtdbk_sp_shannon']=sp_div
  sample_tax_shannon[i,'Midas_sp_shannon']=sp_midas_div
}

df.all2=merge(df.all[,c('species_name','sample_id','subject_id','total_reads_orig','shannon_diversity','species_richness','polymorphism_rate'),],sample_tax_shannon,by='sample_id')


########
###  gams

##phylum level
phyl.1 <- gam(polymorphism_rate ~ s(species_name,phyla_shannon,bs="fs")+
                s(sample_id,bs="re")+
                s(subject_id,bs="re")+
                s(phyla_shannon)+
                s(total_reads_orig),
              data=df.all2,
              family=betar(link="logit"))
save(phyl.1,file='phyl.1.RData')
pdf('phyl.1_check.pdf')
par(mfrow=c(2,2))
gam.check(phyl.1,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()

##class
clas.1 <- gam(polymorphism_rate ~ s(species_name,class_shannon,bs="fs")+
                s(sample_id,bs="re")+
                s(subject_id,bs="re")+
                s(class_shannon)+
                s(total_reads_orig),
              data=df.all2,
              family=betar(link="logit"))
save(clas.1,file='clas.1.RData')

pdf('clas.1_check.pdf')
par(mfrow=c(2,2))
gam.check(clas.1,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()

order.1 <- gam(polymorphism_rate ~ s(species_name,order_shannon,bs="fs")+
                 s(sample_id,bs="re")+
                 s(subject_id,bs="re")+
                 s(order_shannon)+
                 s(total_reads_orig),
               data=df.all2,
               family=betar(link="logit"))
save(order.1,file='order.1.RData')
pdf('order.1_check.pdf')
par(mfrow=c(2,2))
gam.check(order.1,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()

fam.1 <- gam(polymorphism_rate ~ s(species_name,family_shannon,bs="fs")+
               s(sample_id,bs="re")+
               s(subject_id,bs="re")+
               s(family_shannon)+
               s(total_reads_orig),
             data=df.all2,
             family=betar(link="logit"))
save(fam.1,file='fam.1.RData')
pdf('fam.1_check.pdf')
par(mfrow=c(2,2))
gam.check(fam.1,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()

genus.1 <- gam(polymorphism_rate ~ s(species_name,genus_shannon,bs="fs")+
                 s(sample_id,bs="re")+
                 s(subject_id,bs="re")+
                 s(genus_shannon)+
                 s(total_reads_orig),
               data=df.all2,
               family=betar(link="logit"))
save(genus.1,file='genus.1.RData')
pdf('genus.1_check.pdf')
par(mfrow=c(2,2))
gam.check(genus.1,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()


### Gams with polymorphism rate as a function of richness at different taxonomic levels
### estimate richness in higher taxonomic levels

i=1

for (i in 1:length(sample_list)) {
  sample_name=as.character(sample_list[i])
  sample_tax_richness[i,'sample_id']=sample_name
  
  my_data=midas_taxo[midas_taxo$sample_id==sample_name & midas_taxo$abundance !=0,]
  my_data=midas_taxo[midas_taxo$sample_id==sample_name & midas_taxo$abundance !=0,]
  sample_tax_richness[i,'phyla_nb']=length(unique(my_data$phylum))
  sample_tax_richness[i,'class_nb']=length(unique(my_data$class))
  sample_tax_richness[i,'order_nb']=length(unique(my_data$order))
  sample_tax_richness[i,'family_nb']=length(unique(my_data$family))
  sample_tax_richness[i,'genus_nb']=length(unique(my_data$genus))
  sample_tax_richness[i,'species_nb']=length(unique(my_data$species_name))
  sample_tax_richness[i,'gtdbk_species_nb']=length(unique(my_data$GTDBK_species))
  
}

df.all2=merge(df.all[,c('species_name','sample_id','subject_id','total_reads_orig','shannon_diversity','species_richness','polymorphism_rate'),],sample_tax_richness,by='sample_id')

phyl.rich.1 <- gam(polymorphism_rate ~ s(species_name,phyla_nb,bs="fs",k=5)+
                     s(sample_id,bs="re",k=5)+
                     s(subject_id,bs="re",k=5)+
                     s(phyla_nb,k=5)+
                     s(total_reads_orig,k=5),
                   data=df.all2,
                   family=betar(link="logit"))
save(phyl.rich.1,file='phyl.rich.1.RData')

pdf('phyl.rich.1_check.pdf')
par(mfrow=c(2,2))
gam.check(phyl.rich.1,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()


clas.rich.1 <- gam(polymorphism_rate ~ s(species_name,class_nb,bs="fs")+
                     s(sample_id,bs="re")+
                     s(subject_id,bs="re")+
                     s(class_nb)+
                     s(total_reads_orig),
                   data=df.all2,
                   family=betar(link="logit"))
save(clas.rich.1,file='clas.rich.1.RData')

pdf('clas.rich.1_check.pdf')
par(mfrow=c(2,2))
gam.check(clas.rich.1,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()

order.rich.1 <- gam(polymorphism_rate ~ s(species_name,order_nb,bs="fs")+
                      s(sample_id,bs="re")+
                      s(subject_id,bs="re")+
                      s(order_nb)+
                      s(total_reads_orig),
                    data=df.all2,
                    family=betar(link="logit"))
save(order.rich.1,file='order.rich.1.RData')

pdf('order.rich.1_check.pdf')
par(mfrow=c(2,2))
gam.check(order.rich.1,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()

fam.rich.1 <- gam(polymorphism_rate ~ s(species_name,family_nb,bs="fs")+
                    s(sample_id,bs="re")+
                    s(subject_id,bs="re")+
                    s(family_nb)+
                    s(total_reads_orig),
                  data=df.all2,
                  family=betar(link="logit"))
save(fam.rich.1,file='fam.rich.1.RData')

pdf('fam.rich.1_check.pdf')
par(mfrow=c(2,2))
gam.check(fam.rich.1,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()

genus.rich.1 <- gam(polymorphism_rate ~ s(species_name,genus_nb,bs="fs")+
                      s(sample_id,bs="re")+
                      s(subject_id,bs="re")+
                      s(genus_nb)+
                      s(total_reads_orig),
                    data=df.all2,
                    family=betar(link="logit"))
save(genus.rich.1,file='genus.rich.1.Rdata')

pdf('genus.1_check.pdf')
par(mfrow=c(2,2))
gam.check(gen.rich.1,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()




