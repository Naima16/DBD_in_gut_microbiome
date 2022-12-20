#### glmms with strain count as a function of diversity
#### Naima Madi dec 2022

library(glmmTMB)
library(performance)
library(DHARMa)
library(bbmle) ## AICtab

##  strain count table / species
gene_strain_df=read.table('~/gene_counts_strain_num_per_sample_species.csv',sep=',',header=T,check.names = F)
gene_strain_df$sample_id=as.character(gene_strain_df$sample_id)

###  enlève le c à la fin des noms de certains samples
for (i in 1:dim(gene_strain_df)[1] )
{
  if (grepl('c',gene_strain_df[i,'sample_id']))  
  {
    name=gene_strain_df[i,'sample_id']
    gene_strain_df[i,'sample_id']=strsplit(name, "c")[[1]][1]
    
  }
}

strain_df=gene_strain_df[,c('sample_id','species_name','strain_number')]

##coverage 
coverage=read.table('~/depths_geq_20_TF.csv',sep=',',header=T,row.names = 1, stringsAsFactors = FALSE, quote = "", comment.char = "",check.names = F)
colnames(coverage)

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
range(gene_strain_df$strain_number) #1-4

########
# coverage:Attached is a table which contains entries indicating whether or not species had greater than 100 sites with more 
# than 20x coverage in a given sample
sample_data=read.table('~/sample_covariates.txt',sep='\t',header=T)

df.strain=merge(strain_df1,sample_data[,c('sample_id','subject_id','shannon_diversity','species_richness','total_reads_orig')],by='sample_id')
names(df.strain)[names(df.strain)=="species_name"] <- "species_id"
names(df.strain)[names(df.strain)=="strain_number"] <- "strain_nb"

df.strain$species_id=as.factor(df.strain$species_id)  
df.strain$sample_id=as.factor(df.strain$sample_id)
df.strain$subject_id=as.factor(df.strain$subject_id)

### standardize the predixtors to mean 0 and variance 1 
datsc <- df.strain

pvar='shannon_diversity'
datsc[pvar] <- lapply(df.strain[pvar],scale)

pvar1='total_reads_orig'
datsc[pvar1] <- lapply(df.strain[pvar1],scale)

p1.shannon=glmmTMB(strain_nb ~ shannon_diversity  + total_reads_orig + 
                     (shannon_diversity|species_id) + (1|sample_id) + (1|subject_id),
                     data=datsc,family = truncated_poisson(link = "log"))
summary(p1.shannon)


nb1.shanon=glmmTMB(strain_nb ~ shannon_diversity +total_reads_orig  + 
                   (shannon_diversity|species_id) + (1|sample_id) +(1|subject_id),
                   data=datsc,family = truncated_nbinom1(link = "log"),
                   control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nb1.shanon)

nb1.shanon.sing=glmmTMB(strain_nb ~ shannon_diversity +total_reads_orig  + 
                        (shannon_diversity|species_id)  +(1|subject_id),
                        data=datsc,family = truncated_nbinom1(link = "log"),
                        control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nb1.shanon.sing)

nb2.shanon=glmmTMB(strain_nb ~ shannon_diversity +total_reads_orig + 
                   (shannon_diversity|species_id) + (1|sample_id) +(1|subject_id),
                   data=datsc,family = truncated_nbinom2(link = "log"),
                   control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nb2.shanon)

library(bbmle)
AICtab(p1.shannon,nb1.shanon,nb2.shanon)

check_overdispersion(p1.shannon)
# Overdispersion detected. 

## test random effects significance
nb2.shanon=glmmTMB(strain_nb ~ shannon_diversity +total_reads_orig + (shannon_diversity|species_id) + (1|sample_id) +(1|subject_id),
                   data=datsc,family = truncated_nbinom2(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

r1=glmmTMB(strain_nb ~ shannon_diversity +total_reads_orig + (shannon_diversity|species_id)+(1|sample_id),
           data=datsc,family = truncated_nbinom2(link = "log"))
summary(r1)

r2=glmmTMB(strain_nb ~ shannon_diversity +total_reads_orig + (shannon_diversity|species_id)  ,
           data=datsc,family = truncated_nbinom2(link = "log")) #,

r3=glmmTMB(strain_nb ~ shannon_diversity +total_reads_orig + (1|species_id)  ,
           data=datsc,family = truncated_nbinom2(link = "log"),
           control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(r3)

r4=glmmTMB(strain_nb ~ shannon_diversity +total_reads_orig   ,
           data=datsc,family = truncated_nbinom2(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nb2.shanon,r1,r2,r3,r4)

save(r3,file='r3.RData')

## test fixed effects significance
drop1(r3,test='Chisq')

### check diagnostic plots
sims_sp <- simulateResiduals(r3)
png('r3_dharma.png')
plot(sims_sp,quantreg = T)
dev.off()

## test significance / null model
nul=glmmTMB(strain_nb ~ 1 + (1|species_id)  ,
            data=datsc,family = truncated_nbinom2(link = "log"),
            control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(r3,nul)

#### R2
r2(r3)


###############################
####### strain number as a function of richness estimated on all data
##############################

datsc <- df.strain

pvar0='species_richness'
datsc[pvar0] <- lapply(df.strain[pvar0],scale)

pvar2='total_reads_orig'
datsc[pvar2] <- lapply(df.strain[pvar2],scale)

#model with richness on all data, truncated poisson 
poi1.rich=glmmTMB(strain_nb ~ species_richness + total_reads_orig + (1|species_id) + 
                  (1|sample_id) +(1|subject_id),
                  data=datsc,family = truncated_poisson(link = "log"))
summary(poi1.rich)

#truncated negative binom1
nb1.rich=glmmTMB(strain_nb ~ species_richness + total_reads_orig + 
                 (species_richness|species_id) + (1|sample_id) +(1|subject_id),
                 data=datsc,family = truncated_nbinom1(link = "log"),
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nb1.rich)

#truncated negative binom2
nb2.rich=glmmTMB(strain_nb ~ species_richness + total_reads_orig + 
                 (1|species_id) + (1|sample_id) +(1|subject_id),
                 data=datsc,family = truncated_nbinom2(link = "log"),
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nb2.rich)

AICtab(poi1.rich,nb1.rich,nb2.rich)

check_overdispersion(poi1.rich)

## test random effects
nb1.rich=glmmTMB(strain_nb ~ species_richness + total_reads_orig + (species_richness|species_id) + (1|sample_id) +(1|subject_id),
                 data=datsc,family = truncated_nbinom1(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nb1.rich)
r1=glmmTMB(strain_nb ~ species_richness + total_reads_orig + (species_richness|species_id) + (1|sample_id) ,
           data=datsc,family = truncated_nbinom1(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r2=glmmTMB(strain_nb ~ species_richness + total_reads_orig + (species_richness|species_id) ,
           data=datsc,family = truncated_nbinom1(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r3=glmmTMB(strain_nb ~ species_richness + total_reads_orig + (1|species_id) ,
           data=datsc,family = truncated_nbinom1(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r4=glmmTMB(strain_nb ~ species_richness + total_reads_orig  ,
           data=datsc,family = truncated_nbinom1(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
anova(nb1.rich,r1,r2,r3,r4)

save(nb1.rich,file='nb1.rich.RData')

sim1=simulateResiduals(nb1.rich)
png('nb1.rich_dharma.png')
plot(sim1,quantreg = T)
dev.off()

### test fixed effects
drop1(nb1.rich,.~species_richness + total_reads_orig,test='Chisq')


## test significance / null model
nul=glmmTMB(strain_nb ~ 1 + (species_richness|species_id) + (1|sample_id) +(1|subject_id),
            data=datsc,family = truncated_nbinom1(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nb1.rich,nul)

r2(nb1.rich)

###################
###################
### glmm with strain number as a function of rarefied richness
##################
## add raref richness

sample_data.raref=read.table('~/sample_covariates_rarefied_20m.txt',sep='\t',header=T)
colnames(sample_data.raref)

##[1] "sample_id"         "subject_id"        "shannon_diversity" "species_richness"  "total_reads_MIDAS" "total_reads_orig" 

df.strain2=merge(strain_df1,sample_data.raref[,c('sample_id','subject_id','species_richness')],by='sample_id')

##add total reads 
df.strain=merge(df.strain2,sample_data[,c('sample_id','total_reads_orig')],by='sample_id')

datsc <- df.strain

pvar='species_richness'
pvar2='total_reads_orig'

datsc[pvar] <- lapply(df.strain[pvar],scale)
datsc[pvar2] <- lapply(df.strain[pvar2],scale)

###### truncated poisson
p1.cov=glmmTMB(strain_nb ~ species_richness + total_reads_orig+ (species_richness|species_id) + (1|sample_id) +(1|subject_id),
               data=datsc,family = truncated_poisson(link = "log"))
summary(p1.cov)

###### truncated negative binomial 1
nb1.cov=glmmTMB(strain_nb ~ species_richness +total_reads_orig  + (species_richness|species_id) + (1|sample_id) +(1|subject_id),
                data=datsc,family = truncated_nbinom1(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nb1.cov)

###### truncated negative binomial 2
nb2.cov=glmmTMB(strain_nb ~ species_richness +total_reads_orig  + (species_richness|species_id) + (1|sample_id) +(1|subject_id),
                data=datsc,family = truncated_nbinom2(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nb2.cov)


AICtab(p1.cov,nb1.cov,nb2.cov)

check_overdispersion(p1.cov)  #Overdispersion detected.

## test random effects (truncated neg bin 2 model)
nb2.cov=glmmTMB(strain_nb ~ species_richness +total_reads_orig  + (species_richness|species_id) + (1|sample_id) +(1|subject_id),
                data=datsc,family = truncated_nbinom2(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

r1=glmmTMB(strain_nb ~ species_richness +total_reads_orig  + (species_richness|species_id) + (1|sample_id) ,
           data=datsc,family = truncated_nbinom2(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

r2=glmmTMB(strain_nb ~ species_richness +total_reads_orig  + (species_richness|species_id)  ,
           data=datsc,family = truncated_nbinom2(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

r3=glmmTMB(strain_nb ~ species_richness +total_reads_orig  + (1|species_id)  ,
           data=datsc,family = truncated_nbinom2(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

r4=glmmTMB(strain_nb ~ species_richness +total_reads_orig    ,
           data=datsc,family = truncated_nbinom2(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nb2.cov,r1,r2,r3,r4)

## test fixed effects
drop1(r3,test='Chisq')

sims_sp <- simulateResiduals(r3)
png('r3_dharma.png')
plot(sims_sp,quantreg = T)
dev.off()

null.r3=glmmTMB(strain_nb ~ 1  + (1|species_id)  ,
                data=datsc,family = truncated_nbinom2(link = "log"),
                control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(r3,null.r3)

r2(r3)
save(r3,file='r3.RData')

##########################################################
######### strain count as a function of shannon estimated at higher taxonomic levels
##################################################
#########################################################
###################

df.strain=merge(strain_df1,sample_data[,c('sample_id','subject_id','shannon_diversity','species_richness','total_reads_orig')],by='sample_id')
names(df.strain)[names(df.strain)=="species_name"] <- "species_id"
names(df.strain)[names(df.strain)=="strain_number"] <- "strain_nb"


df.strain$species_id=as.factor(df.strain$species_id)  
df.strain$sample_id=as.factor(df.strain$sample_id)
df.strain$subject_id=as.factor(df.strain$subject_id)


## taxonomy from GTDNK pipeline
taxo_1=read.table('~/part1/gtdbtk.bac120.summary.tsv',sep='\t',header=T)
taxo_2=read.table('~/part2/gtdbtk.bac120.summary.tsv',sep='\t',header=T)

taxo=rbind(taxo_1,taxo_2)

taxo.df=taxo[,c('user_genome','classification')]

for (i in 1:dim(taxo.df)[1] ) {
   taxo.df[i,'genus']=unlist(strsplit(unlist(strsplit(as.character(taxo.df[i,]$classification),";"))[6],"__"))[2]
  taxo.df[i,'family']=unlist(strsplit(unlist(strsplit(as.character(taxo.df[i,]$classification),";"))[5],"__"))[2]
  taxo.df[i,'order']=unlist(strsplit(unlist(strsplit(as.character(taxo.df[i,]$classification),";"))[4],"__"))[2]
  taxo.df[i,'class']=unlist(strsplit(unlist(strsplit(as.character(taxo.df[i,]$classification),";"))[3],"__"))[2]
  taxo.df[i,'phylum']=unlist(strsplit(unlist(strsplit(as.character(taxo.df[i,]$classification),";"))[2],"__"))[2]
  taxo.df[i,'GTDBK_species']=unlist(strsplit(unlist(strsplit(as.character(taxo.df[i,]$classification),";"))[7],"__"))[2]
  taxo.df[i,'species_id']=unlist(strsplit(as.character(taxo.df[i,]$user_genome), "[.]"))[1]
}

df.all=merge(df.strain,taxo.df[,c('species_id','GTDBK_species','genus','family','order','class','phylum')],by='species_id')

df.all$genus=as.factor(df.all$genus)
df.all$family=as.factor(df.all$family)
df.all$order=as.factor(df.all$order)
df.all$class=as.factor(df.all$class)
df.all$phylum=as.factor(df.all$phylum)
df.all$species_id=as.factor(df.all$species_id)
df.all$GTDBK_species=as.factor(df.all$GTDBK_species)
df.all$sample_id=as.factor(df.all$sample_id)


midas_abund=read.table('~/relative_abundance.txt',sep='\t',header=T,check.names = F)
dim(midas_abund)  #5952  470

rownames(midas_abund)=midas_abund$species_id
midas_abund=midas_abund[,-1]

commun_data=t(midas_abund)

##convert to num
commun_data_n=data.matrix(as.data.frame(commun_data))  ##sample*species

#remove species present in 0 sample
commun_data_n.nonzero=commun_data_n[,-(which(colSums(commun_data_n)==0))]


library(vegan)
shannon_midas=diversity(commun_data_n.nonzero,index='shannon')  ## la même que quand j'enlève pas les 0
shannon_midas=diversity(commun_data_n,index='shannon')
#####
shannon_midas=as.data.frame(shannon_midas)


library(reshape2)
midas_abund$species_id=rownames(midas_abund)
midas_abund_long=melt(midas_abund,id='species_id')
colnames(midas_abund_long)[2] = 'sample_id'
colnames(midas_abund_long)[3] = 'abundance'


#######
midas_taxo=merge(midas_abund_long,taxo.df,by='species_id')

sample_list=unique(midas_abund_long$sample_id)
sample_tax_shannon=data.frame()

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
  t=as.numeric(my_data$abundance)  ## ça  donne le meme resultat que ds la matrice de daisy
  sp_midas_div=diversity(t,index='shannon')
  
  
  sample_tax_shannon[i,'phyla_shannon']=phylum_div
  sample_tax_shannon[i,'class_shannon']=class_div
  sample_tax_shannon[i,'order_shannon']=order_div
  sample_tax_shannon[i,'family_shannon']=family_div
  sample_tax_shannon[i,'genus_shannon']=genus_div
  sample_tax_shannon[i,'gtdbk_sp_shannon']=sp_div
  sample_tax_shannon[i,'Midas_sp_shannon']=sp_midas_div
}



df.all2=merge(df.all[,c('species_id','sample_id','subject_id','total_reads_orig','shannon_diversity','species_richness','strain_nb'),],sample_tax_shannon,by='sample_id')


datsc2=df.all2

pvar1='phyla_shannon'
pvar2='class_shannon'
pvar3='order_shannon'
pvar4='family_shannon'
pvar5='genus_shannon'
pvar6='gtdbk_sp_shannon'
pvar7='total_reads_orig'
pvar8='shannon_diversity'
pvar9='species_richness'


datsc2[pvar1] <- lapply(df.all2[pvar1],scale)
datsc2[pvar2] <- lapply(df.all2[pvar2],scale)
datsc2[pvar3] <- lapply(df.all2[pvar3],scale)
datsc2[pvar4] <- lapply(df.all2[pvar4],scale)
datsc2[pvar5] <- lapply(df.all2[pvar5],scale)
datsc2[pvar6] <- lapply(df.all2[pvar6],scale)
datsc2[pvar7] <- lapply(df.all2[pvar7],scale)
datsc2[pvar8] <- lapply(df.all2[pvar8],scale)
datsc2[pvar9] <- lapply(df.all2[pvar9],scale)




##phyla
poi.phyla=glmmTMB(strain_nb ~ phyla_shannon + total_reads_orig +(phyla_shannon|species_id) + (1|sample_id) + (1|subject_id) ,
                  data=datsc2,family = truncated_poisson(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(poi.phyla)

##nbinom1
nbinom1.phyla=glmmTMB(strain_nb ~ phyla_shannon + total_reads_orig +  (phyla_shannon|species_id)+(1|sample_id) + (1|subject_id),
                      data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom1.phyla) 


##nbinom2
nbinom2.phyla=glmmTMB(strain_nb ~ phyla_shannon + total_reads_orig +  (phyla_shannon|species_id)+(1|sample_id) + (1|subject_id),
                      data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

summary(nbinom2.phyla)


AICtab(poi.phyla,nbinom1.phyla,nbinom2.phyla)

check_overdispersion(poi.phyla)

###overdispersion so we keep nbinom1
nbinom2.phyla=glmmTMB(strain_nb ~ phyla_shannon + total_reads_orig +  (phyla_shannon|species_id)+(1|sample_id) + (1|subject_id),
                      data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))


## test random effects

m1=glmmTMB(strain_nb ~ phyla_shannon + total_reads_orig +  (phyla_shannon|species_id)+ (1|subject_id) ,
           data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m2=glmmTMB(strain_nb ~ phyla_shannon + total_reads_orig +  (phyla_shannon|species_id) ,
           data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m3=glmmTMB(strain_nb ~ phyla_shannon + total_reads_orig +  (1|species_id) ,
           data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m4=glmmTMB(strain_nb ~ phyla_shannon + total_reads_orig  ,
           data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nbinom2.phyla,m1,m2,m3,m4)

## test fixed  effects
drop1(nbinom2.phyla,.~phyla_shannon + total_reads_orig,test='Chisq')  ##pbm de convergence

save(nbinom2.phyla,file='nbinom2.phyla.RData')


library(DHARMa)
sims_nb <- simulateResiduals(nbinom2.phyla)

png('nbinom2.phyla_dharma.png')
plot(sims_nb,quantreg = T)
dev.off()

r2(nbinom2.phyla)

##  heteroscedasticity
##  add dispformula
nbinom2.phyla.disp=glmmTMB(strain_nb ~ phyla_shannon + total_reads_orig +  (phyla_shannon|species_id)+(1|sample_id) + (1|subject_id),
                           data=datsc2,dispformula=~phyla_shannon,family = truncated_nbinom2(link='log'),
                           control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom2.phyla.disp) #ne converge pas

sims_nb <- simulateResiduals(nbinom2.phyla.disp)

png('nbinom2.phyla.disp_dharma.png')
plot(sims_nb,quantreg = T)
dev.off()  ##has not changed with dispformula.


### class
poi.class=glmmTMB(strain_nb ~ class_shannon + total_reads_orig +  (class_shannon|species_id)+(1|sample_id)+(1|subject_id) ,
                  data=datsc2,family = truncated_poisson(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(poi.class)

nbinom1.class=glmmTMB(strain_nb ~ class_shannon + total_reads_orig +  (class_shannon|species_id)+(1|sample_id)+(1|subject_id) ,
                      data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom1.class)

nbinom2.class=glmmTMB(strain_nb ~ class_shannon + total_reads_orig +  (class_shannon|species_id)+(1|sample_id)+(1|subject_id) ,
                      data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom2.class)

AICtab(poi.class,nbinom1.class,nbinom2.class)


check_overdispersion(poi.class)
# Overdispersion detected.


## test random effects
r1.class=glmmTMB(strain_nb ~ class_shannon + total_reads_orig +  (1|species_id)+(1|subject_id) +(1|sample_id),
                 data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r2.class=glmmTMB(strain_nb ~ class_shannon + total_reads_orig +  (1|species_id)+(1|subject_id) ,
                 data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r3.class=glmmTMB(strain_nb ~ class_shannon + total_reads_orig +  (1|species_id) ,
                 data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r4.class=glmmTMB(strain_nb ~ class_shannon + total_reads_orig ,
                 data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nbinom2.class,r1.class,r2.class,r3.class,r4.class)

summary(nbinom2.class)

drop1(nbinom2.class,.~class_shannon + total_reads_orig ,test='Chisq')

save(nbinom2.class,file='nbinom2.class.RData')

sims_clas <- simulateResiduals(nbinom2.class)

png('nbinom2.class_dharma.png')
plot(sims_clas,quantreg = T)
dev.off()

null=glmmTMB(strain_nb ~ 1 +  (class_shannon|species_id)+(1|sample_id)+(1|subject_id) ,
             data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
anova(nbinom2.class,null)

r2(nbinom2.class)

###order
poi.order=glmmTMB(strain_nb ~ order_shannon + total_reads_orig +  (order_shannon|species_id)+(1|sample_id)+(1|subject_id) ,
                  data=datsc2,family = truncated_poisson(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(poi.order) ## sing

nbinom1.order=glmmTMB(strain_nb ~ order_shannon + total_reads_orig +  (order_shannon|species_id)+(1|sample_id)+(1|subject_id) ,
                      data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom1.order)

nbinom2.order=glmmTMB(strain_nb ~ order_shannon + total_reads_orig +  (order_shannon|species_id)+(1|sample_id)+(1|subject_id) ,
                      data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom2.order)

AICtab(poi.order,nbinom1.order,nbinom2.order)

check_overdispersion(poi.order)
# Overdispersion detected.

summary(nbinom2.order)


## test random effects
r1.order=glmmTMB(strain_nb ~ order_shannon + total_reads_orig +  (order_shannon|species_id)+(1|subject_id) ,
                 data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r2.order=glmmTMB(strain_nb ~ order_shannon + total_reads_orig +  (order_shannon|species_id) ,
                 data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r3.order=glmmTMB(strain_nb ~ order_shannon + total_reads_orig +  (1|species_id),
                 data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r4.order=glmmTMB(strain_nb ~ order_shannon + total_reads_orig ,
                 data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nbinom2.order,r1.order,r2.order,r3.order,r4.order)

drop1(r2.order,.~order_shannon + total_reads_orig,test='Chisq')  


save(r2.order,file='r2.order.RData')

sims_clas <- simulateResiduals(r2.order)
png('r2.order_dharma.png')
plot(sims_clas,quantreg = T)
dev.off()

null=glmmTMB(strain_nb ~ 1 +  (order_shannon|species_id) ,
             data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(r2.order,null)

r2(r2.order)

### family
poi.fam=glmmTMB(strain_nb ~ family_shannon + total_reads_orig +  (family_shannon|species_id)+(1|sample_id) +(1|subject_id),
                data=datsc2,family = truncated_poisson(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(poi.fam)

nbinom1.fam=glmmTMB(strain_nb ~ family_shannon + total_reads_orig +  (family_shannon|species_id)+(1|sample_id) +(1|subject_id),
                    data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom1.fam)

nbinom2.fam=glmmTMB(strain_nb ~ family_shannon + total_reads_orig +  (family_shannon|species_id)+(1|sample_id) +(1|subject_id),
                    data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom2.fam)

AICtab(poi.fam,nbinom1.fam,nbinom2.fam)

check_overdispersion(poi.fam)
# Overdispersion detected.

## test random effects
nbinom2.fam=glmmTMB(strain_nb ~ family_shannon + total_reads_orig +  (family_shannon|species_id)+(1|sample_id) +(1|subject_id),
                    data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

r1.fam=glmmTMB(strain_nb ~ family_shannon + total_reads_orig +  (family_shannon|species_id) +(1|subject_id),
               data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

r2.fam=glmmTMB(strain_nb ~ family_shannon + total_reads_orig +  (family_shannon|species_id),
               data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r3.fam=glmmTMB(strain_nb ~ family_shannon + total_reads_orig+  (1|species_id),
               data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r4.fam=glmmTMB(strain_nb ~ family_shannon + total_reads_orig,
               data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nbinom2.fam,r1.fam,r2.fam,r3.fam,r4.fam)

### test fixed effects
drop1(nbinom2.fam,.~family_shannon + total_reads_orig,test='Chisq')  

null=glmmTMB(strain_nb ~ 1 +  (family_shannon|species_id)+(1|sample_id) +(1|subject_id),
           data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nbinom2.fam,null)

save(nbinom2.fam,file='nbinom2.fam.RData')


sims_fam <- simulateResiduals(nbinom2.fam)
png('nbinom2.fam_dharma.png')
plot(sims_fam,quantreg = T)
dev.off()

r2(nbinom2.fam)

### genus
poi.gen=glmmTMB(strain_nb ~ genus_shannon + total_reads_orig +  (genus_shannon|species_id)+(1|sample_id)+(1|subject_id) ,
                data=datsc2,family = truncated_poisson(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(poi.gen)

nbinom1.gen=glmmTMB(strain_nb ~ genus_shannon + total_reads_orig +  (genus_shannon|species_id)+(1|sample_id)+(1|subject_id) ,
                    data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom1.gen)

nbinom2.gen=glmmTMB(strain_nb ~ genus_shannon + total_reads_orig +  (genus_shannon|species_id)+(1|sample_id)+(1|subject_id) ,
                    data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom2.gen)

AICtab(poi.gen,nbinom1.gen,nbinom2.gen)

check_overdispersion(poi.gen)
# Overdispersion detected.


## test random effects
r1.gen=glmmTMB(strain_nb ~ genus_shannon + total_reads_orig +  (genus_shannon|species_id)+(1|subject_id) ,
               data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r2.gen=glmmTMB(strain_nb ~ genus_shannon + total_reads_orig +  (genus_shannon|species_id) ,
               data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r3.gen=glmmTMB(strain_nb ~ genus_shannon + total_reads_orig +  (1|species_id),
               data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r4.gen=glmmTMB(strain_nb ~ genus_shannon + total_reads_orig ,
               data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nbinom2.gen,r1.gen,r2.gen,r3.gen,r4.gen)

## test fixed effects
drop1(nbinom2.gen,.~genus_shannon + total_reads_orig,test='Chisq')

summary(nbinom2.gen)

save(nbinom2.gen,file='nbinom2.gen.RData')


null=glmmTMB(strain_nb ~ 1 +  (genus_shannon|species_id)+(1|sample_id)+(1|subject_id) ,
             data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
anova(nbinom2.gen,null)

sims_fam <- simulateResiduals(nbinom2.gen)
png('nbinom2.gen_dharma.png')
plot(sims_fam,quantreg = T)
dev.off()

r2(nbinom2.gen)

##############
#################
#####################
## strain count as a function of richness at high taxonomic levels
#####################
#################
##############

sample_list=unique(midas_abund_long$sample_id)
length(sample_list)
sample_tax_richness=data.frame()

for (i in 1:length(sample_list)) {
  
  sample_name=as.character(sample_list[i])
  sample_tax_richness[i,'sample_id']=sample_name
  my_data=midas_taxo[midas_taxo$sample_id==sample_name & midas_taxo$abundance !=0,]
  sample_tax_richness[i,'phyla_nb']=length(unique(my_data$phylum))
  sample_tax_richness[i,'class_nb']=length(unique(my_data$class))
  sample_tax_richness[i,'order_nb']=length(unique(my_data$order))
  sample_tax_richness[i,'family_nb']=length(unique(my_data$family))
  sample_tax_richness[i,'genus_nb']=length(unique(my_data$genus))
  sample_tax_richness[i,'species_nb']=length(unique(my_data$species_id))
  sample_tax_richness[i,'gtdbk_species_nb']=length(unique(my_data$GTDBK_species))
}


df.all2=merge(df.all[,c('species_id','sample_id','subject_id','total_reads_orig','shannon_diversity','species_richness','strain_nb'),],sample_tax_richness,by='sample_id')
df.richness=unique(df.all2[,c('species_id','sample_id','gtdbk_species_nb','species_richness')])


datsc2=df.all2

pvar1='phyla_nb'
pvar2='class_nb'
pvar3='order_nb'
pvar4='family_nb'
pvar5='genus_nb'
pvar6='gtdbk_species_nb'
pvar7='species_richness'
pvar8='total_reads_orig'

datsc2[pvar1] <- lapply(df.all2[pvar1],scale)
datsc2[pvar2] <- lapply(df.all2[pvar2],scale)
datsc2[pvar3] <- lapply(df.all2[pvar3],scale)
datsc2[pvar4] <- lapply(df.all2[pvar4],scale)
datsc2[pvar5] <- lapply(df.all2[pvar5],scale)
datsc2[pvar6] <- lapply(df.all2[pvar6],scale)
datsc2[pvar7] <- lapply(df.all2[pvar7],scale)
datsc2[pvar8] <- lapply(df.all2[pvar8],scale)



##  phyla
## poisson
poi.phyla=glmmTMB(strain_nb ~ phyla_nb + total_reads_orig +  (phyla_nb|species_id) + (1|sample_id) + (1|subject_id) ,
                  data=datsc2,family = truncated_poisson(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(poi.phyla)

##nbinom1
nbinom1.phyla=glmmTMB(strain_nb ~ phyla_nb + total_reads_orig +  (phyla_nb|species_id)+(1|sample_id) + (1|subject_id),
                      data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom1.phyla) 

##nbinom2
nbinom2.phyla=glmmTMB(strain_nb ~ phyla_nb + total_reads_orig +  (phyla_nb|species_id)+(1|sample_id) + (1|subject_id),
                      data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

summary(nbinom2.phyla)

library(bbmle)
AICtab(poi.phyla,nbinom1.phyla,nbinom2.phyla)

check_overdispersion(poi.phyla)
# Overdispersion detected.

### overdispersion so we keep nbinom2
nbinom2.phyla=glmmTMB(strain_nb ~ phyla_nb + total_reads_orig +  (phyla_nb|species_id)+(1|sample_id) + (1|subject_id),
                      data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

## test random effects
m1=glmmTMB(strain_nb ~ phyla_nb + total_reads_orig +  (phyla_nb|species_id)+ (1|subject_id) ,
           data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m2=glmmTMB(strain_nb ~ phyla_nb + total_reads_orig +  (phyla_nb|species_id) ,
           data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m3=glmmTMB(strain_nb ~ phyla_nb + total_reads_orig +  (1|species_id) ,
           data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

m4=glmmTMB(strain_nb ~ phyla_nb + total_reads_orig  ,
           data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nbinom2.phyla,m1,m2,m3,m4)
## test fixed  effects
drop1(m3,.~phyla_nb + total_reads_orig,test='Chisq')  ##pbm de convergence

save(m3,file='m3.RData')

library(DHARMa)
sims_nb <- simulateResiduals(m3)

png('m3_dharma.png')
plot(sims_nb,quantreg = T)
dev.off()

##little heteroscedasticity
##add dispformula
m3.disp=glmmTMB(strain_nb ~ phyla_nb + total_reads_orig +  (1|species_id),
                data=datsc2,dispformula=~phyla_nb,family = truncated_nbinom2(link='log'),
                control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(m3.disp)
anova(m3.disp,m3)  ### disp non significatif

sm <- simulateResiduals(m3.disp)
plot(sm,quantreg = T)  # same plot as the model without  disp

### null model test
m3.nul=glmmTMB(strain_nb ~ 1 +  (1|species_id) ,
               data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(m3,m3.nul)
r2(m3)
#######
### class
poi.class=glmmTMB(strain_nb ~ class_nb + total_reads_orig +  (class_nb|species_id)+(1|sample_id)++(1|subject_id) ,
                  data=datsc2,family = truncated_poisson(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(poi.class)

nbinom1.class=glmmTMB(strain_nb ~ class_nb + total_reads_orig +  (class_nb|species_id)+(1|sample_id)++(1|subject_id) ,
                      data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom1.class)
nbinom1.class.sing=glmmTMB(strain_nb ~ class_nb + total_reads_orig +  (1|species_id)+(1|sample_id)+(1|subject_id) ,
                           data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom1.class.sing)
nbinom2.class=glmmTMB(strain_nb ~ class_nb + total_reads_orig +  (class_nb|species_id)+(1|sample_id)++(1|subject_id) ,
                      data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom2.class)
nbinom2.class.sing=glmmTMB(strain_nb ~ class_nb + total_reads_orig +  (1|species_id)+(1|sample_id)+(1|subject_id) ,
                           data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom2.class.sing)
AICtab(poi.class,nbinom1.class.sing,nbinom2.class.sing)

check_overdispersion(poi.class)
#Overdispersion detected.

##random
r1.class=glmmTMB(strain_nb ~ class_nb + total_reads_orig +  (1|species_id)+(1|subject_id) ,
                 data=datsc2,family = truncated_nbinom1(link='log'),
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r2.class=glmmTMB(strain_nb ~ class_nb + total_reads_orig +  (1|species_id) ,
                 data=datsc2,family = truncated_nbinom1(link='log'),
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r3.class=glmmTMB(strain_nb ~ class_nb + total_reads_orig  ,
                 data=datsc2,family = truncated_nbinom1(link='log'),
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nbinom1.class.sing,r1.class,r2.class,r3.class)

drop1(nbinom1.class.sing,.~ class_nb + total_reads_orig,test='Chisq')  ##pbm de convergence

save(nbinom1.class.sing,file='nbinom1.class.sing.RData')
sims_clas <- simulateResiduals(nbinom1.class.sing)
png('nbinom1.class.sing_dharma.png')
plot(sims_clas,quantreg = T)
dev.off()

null=glmmTMB(strain_nb ~ 1 +  (1|species_id)+(1|sample_id)++(1|subject_id) ,
             data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
anova(nbinom1.class.sing,null)
r2(nbinom1.class.sing)


###order
poi.order=glmmTMB(strain_nb ~ order_nb + total_reads_orig +  (order_nb|species_id)+(1|sample_id)+(1|subject_id) ,
                  data=datsc2,family = truncated_poisson(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

nbinom1.order=glmmTMB(strain_nb ~ order_nb + total_reads_orig +  (order_nb|species_id)+(1|sample_id)+(1|subject_id) ,
                      data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom1.order)
nbinom1.order.sing=glmmTMB(strain_nb ~ order_nb + total_reads_orig +  (1|species_id)+(1|sample_id)+(1|subject_id) ,
                           data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

nbinom2.order=glmmTMB(strain_nb ~ order_nb + total_reads_orig +  (order_nb|species_id)+(1|sample_id)+(1|subject_id) ,
                      data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(nbinom2.order)

AICtab(poi.order,nbinom1.order.sing,nbinom2.order)
check_overdispersion(poi.order) #Overdispersion detected

summary(nbinom2.order)

## test random effects
r1.order=glmmTMB(strain_nb ~ order_nb + total_reads_orig +  (order_nb|species_id)+(1|subject_id) ,
                 data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r2.order=glmmTMB(strain_nb ~ order_nb + total_reads_orig +  (order_nb|species_id) ,
                 data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r3.order=glmmTMB(strain_nb ~ order_nb + total_reads_orig + (1|species_id),
                 data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

r4.order=glmmTMB(strain_nb ~ order_nb + total_reads_orig,
                 data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nbinom2.order,r1.order,r2.order,r3.order,r4.order)

save(nbinom2.order,file='nbinom2.order.RData')

null=glmmTMB(strain_nb ~ 1 +  (order_nb|species_id)+(1|sample_id)+(1|subject_id) ,
             data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nbinom2.order,null)
sims_clas <- simulateResiduals(nbinom2.order)
png('nbinom2.order_dharma.png')
plot(sims_clas,quantreg = T)
dev.off()
r2(nbinom2.order)

### family
poi.fam=glmmTMB(strain_nb ~ family_nb + total_reads_orig +  (family_nb|species_id)+(1|sample_id) +(1|subject_id),
                data=datsc2,family = truncated_poisson(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

nbinom1.fam=glmmTMB(strain_nb ~ family_nb + total_reads_orig +  (family_nb|species_id)+(1|sample_id) +(1|subject_id),
                    data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

nbinom2.fam=glmmTMB(strain_nb ~ family_nb + total_reads_orig +  (family_nb|species_id)+(1|sample_id) +(1|subject_id),
                    data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
AICtab(poi.fam,nbinom1.fam,nbinom2.fam)

check_overdispersion(poi.fam)

# Overdispersion detected.
summary(nbinom1.fam)

## test random effects
r1.fam=glmmTMB(strain_nb ~ family_nb + total_reads_orig +  (family_nb|species_id)+(1|subject_id),
               data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

r2.fam=glmmTMB(strain_nb ~ family_nb + total_reads_orig +  (family_nb|species_id),
               data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r3.fam=glmmTMB(strain_nb ~ family_nb + total_reads_orig +  (1|species_id),
               data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r4.fam=glmmTMB(strain_nb ~ family_nb + total_reads_orig ,
               data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
anova(nbinom1.fam,r1.fam,r2.fam,r3.fam,r4.fam)

save(nbinom1.fam,file='nbinom1.fam.RData')
null=glmmTMB(strain_nb ~ 1 +  (family_nb|species_id)+(1|sample_id) +(1|subject_id),
             data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nbinom1.fam,null)

sims_fam <- simulateResiduals(nbinom1.fam)
png('nbinom1.fam_dharma.png')
plot(sims_fam,quantreg = T)
dev.off()

r2(nbinom1.fam)

### genus
poi.gen=glmmTMB(strain_nb ~ genus_nb + total_reads_orig +  (genus_nb|species_id)+(1|sample_id)+(1|subject_id) ,
                data=datsc2,family = truncated_poisson(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

nbinom1.gen=glmmTMB(strain_nb ~ genus_nb + total_reads_orig +  (genus_nb|species_id)+(1|sample_id)+(1|subject_id) ,
                    data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

nbinom2.gen=glmmTMB(strain_nb ~ genus_nb + total_reads_orig +  (genus_nb|species_id)+(1|sample_id)+(1|subject_id) ,
                    data=datsc2,family = truncated_nbinom2(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

AICtab(poi.gen,nbinom1.gen,nbinom2.gen)

check_overdispersion(poi.gen)
# Overdispersion detected.

## test random effects
r1.gen=glmmTMB(strain_nb ~ genus_nb + total_reads_orig +  (genus_nb|species_id)+(1|subject_id) ,
               data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r2.gen=glmmTMB(strain_nb ~ genus_nb + total_reads_orig +  (genus_nb|species_id) ,
               data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r3.gen=glmmTMB(strain_nb ~ genus_nb + total_reads_orig +  (1|species_id) ,
               data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
r4.gen=glmmTMB(strain_nb ~ genus_nb + total_reads_orig  ,
               data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
anova(nbinom1.gen,r1.gen,r2.gen,r3.gen,r4.gen)

save(nbinom1.gen,file='nbinom1.gen.RData')

null=glmmTMB(strain_nb ~ 1 +  (genus_nb|species_id)+(1|sample_id)+(1|subject_id) ,
             data=datsc2,family = truncated_nbinom1(link='log'),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
anova(nbinom1.gen,null)

sims_fam <- simulateResiduals(nbinom1.gen)
png('nbinom1.gen_dharma.png')
plot(sims_fam,quantreg = T)
dev.off()

r2(nbinom1.gen)

