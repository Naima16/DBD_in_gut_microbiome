### HMP time series
### GAMs with polymorphism change , gene gains and loss as a function of diversity

library(mgcv)
library(glmmTMB)
library(bbmle)
library(performance)

pi=read.table('~/HMP_delta_polymorphism.csv',sep=',',header=T)
colnames(pi)
# [1] "sample_tp1"         "sample_tp2"         "subject"            "tp1"                "tp2"                "species"
# [7] "delta_polymorphism" "polymorphism_tp1"   "polymorphism_tp2"   "alpha_div_tp1"      "alpha_div_tp2"      "alpha_div_rare_tp1"
# [13] "alpha_div_rare_tp2" "richness_tp1"       "richness_tp2"       "richness_rare_tp1"  "richness_rare_tp2"

sample_data=read.table('~/sample_covariates.txt',sep='\t',header=T)

names(sample_data)[names(sample_data)=="sample_id"] <- "sample_tp1"
names(sample_data)[names(sample_data)=="total_reads_orig"] <- "total_reads_orig1"

pi.1=merge(pi,sample_data[,c('sample_tp1','total_reads_orig1')],by='sample_tp1')

pi.1$subject=as.factor(pi.1$subject)
pi.1$sample_tp1=as.factor(pi.1$sample_tp1)
pi.1$sample_tp2=as.factor(pi.1$sample_tp2)
pi.1$species=as.factor(pi.1$species)

pi.1$log.delta_pi <- log10(1+pi.1$delta_polymorphism)

## la reponse est normale et la relation delta_pi~shannon et non lineaire
library(mgcv)
pi0.log <- gam(log.delta_pi ~ s(total_reads_orig1) + 
             s(alpha_div_tp1) + 
             s(species, alpha_div_tp1, bs = "fs")+
             s(sample_tp1, bs = "re") + 
             s(subject, bs = "re"), 
             data=pi.1)

pdf('pi0.log_check.pdf')
par(mfrow=c(2,2))
gam.check(pi0.log,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()

save(pi0.log,file='pi0.log.RData')

## richness is the predictor
## compare between truncated-poisson,neg binom1 and neg binom2 
## and keep the lowest AIC

pi0.rich.log <- gam(log.delta_pi ~ s(total_reads_orig1) + 
                 s(richness_tp1) + 
                 s(species, richness_tp1, bs = "fs")+
                 s(sample_tp1, bs = "re") + 
                 s(subject, bs = "re") 
               , data=pi.1)

summary(pi0.rich.log)

pdf('pi0.rich.log_check.pdf')
par(mfrow=c(2,2))
gam.check(pi0.rich.log,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()

save(pi0.rich.log,file='pi0.rich.log.RData')

###################################
## raref rich is the predictor
pi0.Rarefrich.log <- gam(log.delta_pi ~   s(richness_rare_tp1) + 
                       s(total_reads_orig1)+
                      s(species, richness_rare_tp1, bs = "fs")+
                      s(sample_tp1, bs = "re") + 
                      s(subject, bs = "re") 
                    , data=pi.2)
summary(pi0.Rarefrich.log)


pdf('pi0.Rarefrich_check.pdf')
par(mfrow=c(2,2))
gam.check(pi0.Rarefrich,pch=19,cex=.3,rl.col=2, rep.col="steelblue",old.style=FALSE,pages=1)
dev.off()

save(pi0.Rarefrich.log,file='pi0.Rarefrich.log.RData')

######
############
## glmms with gene gaines as a function of Shannon 
############
######
genes=read.table('~/HMP_gene_changes_full.csv',sep=',',header=T)
colnames(genes)
# [1] "sample_tp1"         "sample_tp2"         "subject"            "tp1"                "tp2"                "species"           
# [7] "num_gene_gains"     "num_gene_losses"    "alpha_div_tp1"      "alpha_div_tp2"      "alpha_div_rare_tp1" "alpha_div_rare_tp2"
# [13] "richness_tp1"       "richness_tp2"       "richness_rare_tp1"  "richness_rare_tp2"  "polymorphism_tp1"   "polymorphism_tp2"  
#   

genes.1=merge(genes,sample_data[,c('sample_tp1','total_reads_orig1')],by='sample_tp1')

genes.1$subject=as.factor(genes.1$subject)
genes.1$sample_tp1=as.factor(genes.1$sample_tp1)
genes.1$sample_tp2=as.factor(genes.1$sample_tp2)
genes.1$species=as.factor(genes.1$species)


datsc <- genes.1

pvar='alpha_div_tp1'
pvar1='total_reads_orig1'
pvar2='richness_tp1'
pvar3='richness_rare_tp1'

datsc[pvar] <- lapply(genes.1[pvar],scale)
datsc[pvar1] <- lapply(genes.1[pvar1],scale)
datsc[pvar2] <- lapply(genes.1[pvar2],scale)
datsc[pvar3] <- lapply(genes.1[pvar3],scale)

gain.1 <- glmmTMB(num_gene_gains ~ alpha_div_tp1 + total_reads_orig1 + 
                  (alpha_div_tp1|species)+(1|subject)+(1|sample_tp1),
                  data=datsc,family = poisson(link = "log"))
summary(gain.1)

gain.nb1.shan <- glmmTMB(num_gene_gains ~ alpha_div_tp1 + total_reads_orig1 + 
                           (alpha_div_tp1|species)+(1|subject)+(1|sample_tp1),
                         data=datsc,family = nbinom1(link = "log"))
summary(gain.nb1.shan)

gain.nb2.shan <- glmmTMB(num_gene_gains ~ alpha_div_tp1 + total_reads_orig1 + 
                           (alpha_div_tp1|species)+(1|subject)+(1|sample_tp1),
                         data=datsc,family = nbinom2(link = "log"))
summary(gain.nb2.shan)

AICtab(gain.nb2.shan,gain.nb1.shan,gain.1)

## test random effects
gain.nb1.shan <- glmmTMB(num_gene_gains ~ alpha_div_tp1 + total_reads_orig1 + 
                         (alpha_div_tp1|species)+(1|subject)+(1|sample_tp1),
                         data=datsc,family = nbinom1(link = "log"))
gain.r.1=glmmTMB(num_gene_gains ~ alpha_div_tp1 + total_reads_orig1 + 
                 (alpha_div_tp1|species)+(1|subject),
                 data=datsc,family = nbinom1(link = "log"))

gain.r.2=glmmTMB(num_gene_gains ~ alpha_div_tp1 + total_reads_orig1 + 
                 (alpha_div_tp1|species),
                 data=datsc,family = nbinom1(link = "log"))

gain.r.3=glmmTMB(num_gene_gains ~ alpha_div_tp1 + total_reads_orig1 + 
                 (1|species),
                 data=datsc,family = nbinom1(link = "log"))
gain.r.4=glmmTMB(num_gene_gains ~ alpha_div_tp1 + total_reads_orig1,
                 data=datsc,family = nbinom1(link = "log"))

anova(gain.nb1.shan,gain.r.1,gain.r.2,gain.r.3,gain.r.4)

## test fixed effects
drop1(gain.nb1.shan,test='Chisq')

save(gain.nb1.shan,file='gain.nb1.shan.RData')
summary(gain.nb1.shan)

## genes gain rich
gain.poi.rich <- glmmTMB(num_gene_gains ~ richness_tp1 + total_reads_orig1 + 
                         (richness_tp1|species)+(1|subject)+(1|sample_tp1),
                         data=datsc,family = poisson(link = "log"))
summary(gain.poi.rich)
gain.nb1.rich <- glmmTMB(num_gene_gains ~ richness_tp1 + total_reads_orig1 + 
                         (richness_tp1|species)+(1|subject)+(1|sample_tp1),
                         data=datsc,family = nbinom1(link = "log"))
summary(gain.nb1.rich)
gain.nb1.rich.sing <- glmmTMB(num_gene_gains ~ richness_tp1 + total_reads_orig1 + 
                              (1|species)+(1|subject)+(1|sample_tp1),
                              data=datsc,family = nbinom1(link = "log"))
gain.nb2.rich <- glmmTMB(num_gene_gains ~ richness_tp1 + total_reads_orig1 + 
                           (richness_tp1|species)+(1|subject)+(1|sample_tp1),
                         data=datsc,family = nbinom2(link = "log"))
summary(gain.nb2.rich)
gain.nb2.rich.sing <- glmmTMB(num_gene_gains ~ richness_tp1 + total_reads_orig1 + 
                                (1|species)+(1|subject)+(1|sample_tp1),
                              data=datsc,family = nbinom2(link = "log"))
AICtab(gain.poi.rich,gain.nb1.rich,gain.nb2.rich)

## test random effects
# gain.nb1.rich.sing <- glmmTMB(num_gene_gains ~ richness_tp1 + total_reads_orig1 + 
#                                 (1|species)+(1|subject)+(1|sample_tp1),
#                               data=datsc,family = nbinom1(link = "log"))
rand1=glmmTMB(num_gene_gains ~ richness_tp1 + total_reads_orig1 + 
              (1|species)+(1|subject),
              data=datsc,family = nbinom1(link = "log"))

rand2=glmmTMB(num_gene_gains ~ richness_tp1 + total_reads_orig1 + 
              (1|species),
              data=datsc,family = nbinom1(link = "log"))

rand3=glmmTMB(num_gene_gains ~ richness_tp1 + total_reads_orig1,
              data=datsc,family = nbinom1(link = "log"))

anova(gain.nb1.rich,gain.nb1.rich.sing,rand1,rand2,rand3)

drop1(gain.nb1.rich.sing,test='Chisq')

save(gain.nb1.rich.sing,file='gain.nb1.rich.sing.RData')

summary(gain.nb1.rich.sing)

sims_0 <- simulateResiduals(gain.nb1.rich.sing)

png('gain.nb1.rich.sing_dharma.png')
plot(sims_0,quantreg = T)
dev.off()

r2(gain.nb1.rich.sing)


##raref richness
gain.poi.Rarefrich <- glmmTMB(num_gene_gains ~ richness_rare_tp1  + 
                              (richness_rare_tp1|species)+(1|subject)+(1|sample_tp1),
                              data=datsc,family = poisson(link = "log"))

gain.nb1.Rarefrich <- glmmTMB(num_gene_gains ~ richness_rare_tp1+ 
                              (richness_rare_tp1|species)+(1|subject)+(1|sample_tp1),
                              data=datsc,family = nbinom1(link = "log"))
summary(gain.nb1.Rarefrich)

gain.nb1.Rarefrich.sing <- glmmTMB(num_gene_gains ~ richness_rare_tp1+ 
                              (1|species)+(1|subject)+(1|sample_tp1),
                              data=datsc,family = nbinom1(link = "log"))


gain.nb2.Rarefrich <- glmmTMB(num_gene_gains ~ richness_rare_tp1 + 
                            (richness_rare_tp1|species)+(1|subject)+(1|sample_tp1),
                            data=datsc,family = nbinom2(link = "log"))
summary(gain.nb2.Rarefrich)
gain.nb2.Rarefrich.sing <- glmmTMB(num_gene_gains ~ richness_rare_tp1 + 
                                 (1|species)+(1|subject)+(1|sample_tp1),
                                 data=datsc,family = nbinom2(link = "log"))

AICtab(gain.poi.Rarefrich,gain.nb1.Rarefrich.sing,gain.nb2.Rarefrich.sing)

summary(gain.nb1.Rarefrich.sing)

## test random effects
gain.nb1.Rarefrich.sing <- glmmTMB(num_gene_gains ~ richness_rare_tp1+ 
                                   (1|species)+(1|subject)+(1|sample_tp1),
                                   data=datsc,family = nbinom1(link = "log"))


raref.r1=glmmTMB(num_gene_gains ~ richness_rare_tp1+ 
                 (1|species)+(1|subject),
                 data=datsc,family = nbinom1(link = "log"))

raref.r2=glmmTMB(num_gene_gains ~ richness_rare_tp1+ 
                 (1|species),
                 data=datsc,family = nbinom1(link = "log"))

raref.r3=glmmTMB(num_gene_gains ~ richness_rare_tp1,
                 data=datsc,family = nbinom1(link = "log"))

anova(gain.nb1.Rarefrich.sing,raref.r1,raref.r2,raref.r3)

## test fixed effects
drop1(gain.nb1.Rarefrich.sing,test='Chisq')

save(gain.nb1.Rarefrich.sing,file='gain.nb1.Rarefrich.sing.RData')

sims_0 <- simulateResiduals(gain.nb1.Rarefrich.sing)

png('gain.nb1.Rarefrich.sing_dharma.png')
plot(sims_0,quantreg = T)
dev.off()

r2(gain.nb1.Rarefrich.sing)


##################
##############
######
### glmm with gene loss as a function of diversity

loss.poi.shan <- glmmTMB(num_gene_losses ~ alpha_div_tp1 + total_reads_orig1 + 
                         (alpha_div_tp1|species)+(1|subject)+(1|sample_tp1),
                         data=datsc,family = poisson(link = "log"))
summary(loss.poi.shan)

loss.nb1.shan <- glmmTMB(num_gene_losses ~ alpha_div_tp1 + total_reads_orig1 + 
                         (alpha_div_tp1|species)+(1|subject)+(1|sample_tp1),
                         data=datsc,family = nbinom1(link = "log"))
summary(loss.nb1.shan)

loss.nb2.shan <- glmmTMB(num_gene_losses ~ alpha_div_tp1 + total_reads_orig1 + 
                         (alpha_div_tp1|species)+(1|subject)+(1|sample_tp1),
                         data=datsc,family = nbinom2(link = "log"))
summary(loss.nb2.shan)

AICtab(loss.nb2.shan,loss.nb1.shan,loss.poi.shan)

summary(loss.nb2.shan)

## test random effects 
loss.nb2.shan <- glmmTMB(num_gene_losses ~ alpha_div_tp1 + total_reads_orig1 + 
                         (alpha_div_tp1|species)+(1|subject)+(1|sample_tp1),
                         data=datsc,family = nbinom2(link = "log"))
loss.r1=glmmTMB(num_gene_losses ~ alpha_div_tp1 + total_reads_orig1 + 
                (alpha_div_tp1|species)+(1|subject),
                data=datsc,family = nbinom2(link = "log"))

loss.r2=glmmTMB(num_gene_losses ~ alpha_div_tp1 + total_reads_orig1 + 
                (alpha_div_tp1|species),
                data=datsc,family = nbinom2(link = "log"))
loss.r3=glmmTMB(num_gene_losses ~ alpha_div_tp1 + total_reads_orig1 + 
                (1|species),
                data=datsc,family = nbinom2(link = "log"))
loss.r4=glmmTMB(num_gene_losses ~ alpha_div_tp1 + total_reads_orig1 , 
                data=datsc,family = nbinom2(link = "log"))

anova(loss.nb2.shan,loss.r1,loss.r2,loss.r3,loss.r4)

drop1(loss.nb2.shan,test='Chisq')

save(loss.nb2.shan,file='loss.nb2.shan.RData')

library(DHARMa)
sims_0 <- simulateResiduals(loss.nb2.shan)

png('loss.nb2.shan_dharma.png')
plot(sims_0,quantreg = T)
dev.off()

r2(loss.nb2.shan)

### compare to null model
loss.null.1=glmmTMB(num_gene_losses ~ 1 + 
                      (alpha_div_tp1|species)+(1|subject)+(1|sample_tp1),
                    data=datsc,family = nbinom2(link = "log"))
anova(loss.nb2.shan,loss.null.1)


## glmms with gene loss as a function of richness
loss.poi.rich <- glmmTMB(num_gene_losses ~ richness_tp1 + total_reads_orig1 + 
                         (richness_tp1|species)+(1|subject)+(1|sample_tp1),
                         data=datsc,family = poisson(link = "log"))
summary(loss.poi.rich)
loss.nb1.rich <- glmmTMB(num_gene_losses ~ richness_tp1 + total_reads_orig1 + 
                         (richness_tp1|species)+(1|subject)+(1|sample_tp1),
                         data=datsc,family = nbinom1(link = "log"))
summary(loss.nb1.rich)
loss.nb2.rich <- glmmTMB(num_gene_losses ~ richness_tp1 + total_reads_orig1 + 
                           (richness_tp1|species)+(1|subject)+(1|sample_tp1),
                         data=datsc,family = nbinom2(link = "log"))
summary(loss.nb2.rich)
loss.nb2.rich.sing <- glmmTMB(num_gene_losses ~ richness_tp1 + total_reads_orig1 + 
                                (1|species)+(1|subject)+(1|sample_tp1),
                              data=datsc,family = nbinom2(link = "log"))
AICtab(loss.poi.rich,loss.nb1.rich,loss.nb2.rich.sing)

summary(loss.nb2.rich.sing)

save(loss.nb2.rich.sing,file='loss.nb2.rich.sing.RData')


## test random effects
r1=glmmTMB(num_gene_losses ~ richness_tp1 + total_reads_orig1 + 
           (1|species)+(1|subject),
           data=datsc,family = nbinom2(link = "log"))

r2=glmmTMB(num_gene_losses ~ richness_tp1 + total_reads_orig1 + 
           (1|species),
           data=datsc,family = nbinom2(link = "log"))

r3=glmmTMB(num_gene_losses ~ richness_tp1 + total_reads_orig1 ,
           data=datsc,family = nbinom2(link = "log"))

anova(loss.nb2.rich.sing,r1,r2,r3)

## test fixed effects
drop1(loss.nb2.rich.sing,test='Chisq')

## test / null model
null.rich=glmmTMB(num_gene_losses ~ 1 + 
                  (1|species)+(1|subject)+(1|sample_tp1),
                  data=datsc,family = nbinom2(link = "log"))

anova(loss.nb2.rich.sing,null.rich)

r2(loss.nb2.rich.sing)

sims_0 <- simulateResiduals(loss.nb2.rich.sing)
png('loss.nb2.rich.sing_dharma.png')
plot(sims_0,quantreg = T)
dev.off()



#### gene loss as a function of rarefed richness

loss.poi.Rarefrich <- glmmTMB(num_gene_losses ~ richness_rare_tp1  + 
                              (richness_rare_tp1|species)+(1|subject)+(1|sample_tp1),
                              data=datsc,family = poisson(link = "log"))
summary(loss.poi.Rarefrich)
loss.nb1.Rarefrich <- glmmTMB(num_gene_losses ~ richness_rare_tp1  + 
                              (richness_rare_tp1|species)+(1|subject)+(1|sample_tp1),
                              data=datsc,family = nbinom1(link = "log"))

loss.nb2.Rarefrich <- glmmTMB(num_gene_losses ~ richness_rare_tp1  + 
                              (richness_rare_tp1|species)+(1|subject)+(1|sample_tp1),
                              data=datsc,family = nbinom2(link = "log"))
summary(loss.nb2.rich)
loss.nb2.Rarefrich.sing <- glmmTMB(num_gene_losses ~ richness_rare_tp1  + 
                                  (1|species)+(1|subject)+(1|sample_tp1),
                                  data=datsc,family = nbinom2(link = "log"))
summary(loss.nb2.Rarefrich.sing)

AICtab(loss.poi.Rarefrich,loss.nb1.Rarefrich,loss.nb2.Rarefrich.sing)

## test random effects
loss.nb2.Rarefrich.sing <- glmmTMB(num_gene_losses ~ richness_rare_tp1  + 
                                   (1|species)+(1|subject)+(1|sample_tp1),
                                   data=datsc,family = nbinom2(link = "log"))

raref.r1= glmmTMB(num_gene_losses ~ richness_rare_tp1  + 
                  (1|species)+(1|subject),
                  data=datsc,family = nbinom2(link = "log"))


raref.r2= glmmTMB(num_gene_losses ~ richness_rare_tp1  + 
                  (1|species),
                  data=datsc,family = nbinom2(link = "log"))

raref.r3= glmmTMB(num_gene_losses ~ richness_rare_tp1  ,
                  data=datsc,family = nbinom2(link = "log"))

anova(loss.nb2.Rarefrich.sing,raref.r1,raref.r2,raref.r3)

summary(raref.r1)

## test fixed effects
drop1(raref.r1,test='Chisq')

save(raref.r1,file='raref.r1.RData')

## test null model
raref.null=glmmTMB(num_gene_losses ~ 1  + 
                     (1|species)+(1|subject),
                   data=datsc,family = nbinom2(link = "log"))
anova(raref.r1,raref.null)

r2(raref.r1)

sims_0 <- simulateResiduals(raref.r1)

png('raref.r1_dharma.png')
plot(sims_0,quantreg = T)
dev.off()
