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
#nb2.cov=glmmTMB(strain_nb ~ species_richness +total_reads_orig  + (species_richness|species_id) + (1|sample_id) +(1|subject_id),
 #               data=datsc,family = truncated_nbinom2(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

r1_raref=glmmTMB(strain_nb ~ species_richness +total_reads_orig  + (species_richness|species_id) + (1|sample_id) ,
           data=datsc,family = truncated_nbinom2(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

r2_raref=glmmTMB(strain_nb ~ species_richness +total_reads_orig  + (species_richness|species_id)  ,
           data=datsc,family = truncated_nbinom2(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

r3_raref=glmmTMB(strain_nb ~ species_richness +total_reads_orig  + (1|species_id)  ,
           data=datsc,family = truncated_nbinom2(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

r4_raref=glmmTMB(strain_nb ~ species_richness +total_reads_orig    ,
           data=datsc,family = truncated_nbinom2(link = "log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(nb2.cov,r1_raref,r2_raref,r3_raref,r4_raref)

## test fixed effects
drop1(r3_raref,test='Chisq')

sims_sp <- simulateResiduals(r3_raref)
png('r3_raref_dharma.png')
plot(sims_sp,quantreg = T)
dev.off()

null.r3=glmmTMB(strain_nb ~ 1  + (1|species_id)  ,
                data=datsc,family = truncated_nbinom2(link = "log"),
                control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(r3_raref,null.r3)

r2(r3_raref)
save(r3_raref,file='r3_raref.RData')




