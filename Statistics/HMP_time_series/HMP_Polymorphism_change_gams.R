### HMP time series
### GAMs with polymorphism change as a function of diversity

library(mgcv)

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
