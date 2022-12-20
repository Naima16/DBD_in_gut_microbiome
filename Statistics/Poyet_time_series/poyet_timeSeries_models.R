##############
### Gams and glmms with polymorphism change, gene loss and gain as a function of diversity (Poyet time series)
##############

library(mgcv)
library(bbmle)
library(glmmTMB)
library(performance)
library(DHARMa)

gene_poyet=read.table('~/Poyet_gene_changes_full.csv',sep=',',header=T)

coverage=read.table('~/Poyet_read_counts.txt',sep=',',header=T)
names(coverage)[names(coverage)=="run_accession"] <- "sample_tp1"

gene_all1=merge(gene_poyet,coverage,by='sample_tp1')

gene_all1$delta_pi=gene_all1$polymorphism_tp2-gene_all1$polymorphism_tp1

gene_all1$log.delta_pi <- log10(1+gene_all1$delta_pi)

gene_all1$tp2=as.numeric(gene_all1$tp2)
gene_all1$tp1=as.numeric(gene_all1$tp1)
gene_all1$lag=gene_all1$tp2-gene_all1$tp1


gam.pi2<- gam(log.delta_pi ~ ti(alpha_div_tp1,lag)+
                s(read_count)+
                s(species,bs='re')+
                s(subject,bs='re') +
                s(sample_tp1,bs='re')+
                s(sample_tp2,bs='re'),
                data=gene_all1,
                method='ML')

save(gam.pi2,file='gam.pi2.RData')
summary(gam.pi2)

pdf('check_pi2.pdf')
gam.check(gam.pi2)
dev.off()

### Polymorphism change as a function of richness
gam.rich.Logpi1<- gam(log.delta_pi ~ ti(richness_rare_tp1,lag)+
                        s(read_count)+
                        s(species,bs='re')+
                        s(subject,bs='re') +
                        s(sample_tp1,bs='re')+
                        s(sample_tp2,bs='re'),
                      data=gene_all1,select=T,
                      method='ML')
save(gam.rich.Logpi1,file='gam.rich.Logpi1.RData')
summary(gam.rich.Logpi1)
save(gam.rich.Logpi1,file='gam.rich.Logpi1.RData')

pdf('check_rich_Logpi1.pdf')
gam.check(gam.rich.Logpi1)
dev.off()

### gene gain ~ Shannon
datsc=gene_all1

pvar1='alpha_div_tp1'
pvar2='lag'
pvar3='read_count'

datsc[pvar1] <- lapply(gene_all1[pvar1],scale)
datsc[pvar2] <- lapply(gene_all1[pvar2],scale)
datsc[pvar3] <- lapply(gene_all1[pvar3],scale)

## gene gain
gain_shan <- glmmTMB(num_gene_gains ~ alpha_div_tp1*lag + read_count+ (1|species) + (1|subject)+
                       (1|sample_tp1)+ (1|sample_tp2),
                     family = nbinom2(),
                     data=datsc)
save(gain_shan,file='gain_shan.RData')
summary(gain_shan)

## test random effects
gain_shan.r1 <- glmmTMB(num_gene_gains ~ alpha_div_tp1*lag + read_count+ (1|species) + (1|subject)+
                        (1|sample_tp1),
                        family = nbinom2(),
                        data=datsc)


gain_shan.r2 <- glmmTMB(num_gene_gains ~ alpha_div_tp1*lag + read_count+ (1|species) + (1|subject),
                        family = nbinom2(),
                        data=datsc)


gain_shan.r3 <- glmmTMB(num_gene_gains ~ alpha_div_tp1*lag + read_count+ (1|species) ,
                        family = nbinom2(),
                        data=datsc)

gain_shan.r4 <- glmmTMB(num_gene_gains ~ alpha_div_tp1*lag + read_count ,
                        family = nbinom2(),
                        data=datsc)


anova(gain_shan,gain_shan.r1,gain_shan.r2,gain_shan.r3,gain_shan.r4)

drop1(gain_shan,.~alpha_div_tp1*lag + read_count,test='Chisq')
## null model test
gain.shan.nul=glmmTMB(num_gene_gains ~ 1+ (1|species) + (1|subject)+
                      (1|sample_tp1)+ (1|sample_tp2),
                      family = nbinom2(),
                      data=datsc)

anova(gain_shan,gain.shan.nul)

r2(gain_shan)

#################
## gene gain ~ richness
#################

datsc1=gene_all1

pvar1='richness_rare_tp1'
pvar2='lag'
pvar3='read_count'

datsc1[pvar1] <- lapply(gene_all1[pvar1],scale)
datsc1[pvar2] <- lapply(gene_all1[pvar2],scale)
datsc1[pvar3] <- lapply(gene_all1[pvar3],scale)

gain_rich <- glmmTMB(num_gene_gains  ~ richness_rare_tp1*lag + read_count+ (1|species) + (1|subject)+
                    (1|sample_tp1)+
                    (1|sample_tp2),
                     family = nbinom2(),
                     data=datsc1)


gain_rich.r1 <- glmmTMB(num_gene_gains  ~ richness_rare_tp1*lag + read_count+ (1|species) + (1|subject)+
                        (1|sample_tp1),
                        family = nbinom2(),
                        data=datsc1)
gain_rich.r2 <- glmmTMB(num_gene_gains  ~ richness_rare_tp1*lag + read_count+ (1|species) + (1|subject),
                        family = nbinom2(),
                        data=datsc1)

gain_rich.r3 <- glmmTMB(num_gene_gains  ~ richness_rare_tp1*lag + read_count+ (1|species) ,
                        family = nbinom2(),
                        data=datsc1)

gain_rich.r4 <- glmmTMB(num_gene_gains  ~ richness_rare_tp1*lag + read_count ,
                        family = nbinom2(),
                        data=datsc1)

anova(gain_rich,gain_rich.r1,gain_rich.r2,gain_rich.r3,gain_rich.r4)

summary(gain_rich)
save(gain_rich,file='gain_rich.RData')

## test fixed effects
drop1(gain_rich,.~richness_rare_tp1*lag + read_count,test='Chisq')

r2(gain_rich)

## test null model
null.rich.gain <- glmmTMB(num_gene_gains  ~ 1 + (1|species) + (1|subject)+
                          (1|sample_tp1)+
                          (1|sample_tp2),
                          family = nbinom2(),
                          data=datsc1)

anova(gain_rich,null.rich.gain)

#################
## gene loss ~ Shannon
#################

datsc=gene_all1

pvar1='alpha_div_tp1'
pvar2='lag'
pvar3='read_count'

datsc[pvar1] <- lapply(gene_all1[pvar1],scale)
datsc[pvar2] <- lapply(gene_all1[pvar2],scale)
datsc[pvar3] <- lapply(gene_all1[pvar3],scale)

loss_glmmTMB1 <- glmmTMB(num_gene_losses  ~ alpha_div_tp1*lag + read_count+ (1|species) + (1|subject)+
                         (1|sample_tp1)+
                         (1|sample_tp2),
                         family = nbinom1(),
                         data=datsc)

loss_glmmTMB2 <- glmmTMB(num_gene_losses ~ alpha_div_tp1*lag + read_count+ (1|species) + (1|subject)+
                         (1|sample_tp1)+ (1|sample_tp2),
                         family = nbinom2(),
                         data=datsc)

loss_glmmTMB.pois <- glmmTMB(num_gene_losses ~ alpha_div_tp1*lag + read_count+ (1|species) + (1|subject)+
                             (1|sample_tp1)+
                             (1|sample_tp2),
                             family = poisson(),
                             data=datsc)
AICtab(loss_glmmTMB.pois,loss_glmmTMB1,loss_glmmTMB2)

r2(loss_glmmTMB2)

summary(loss_glmmTMB2)

###### test random effects
loss_glmmTMB2.r1 <- glmmTMB(num_gene_losses ~ alpha_div_tp1*lag + read_count+ (1|species) + (1|subject)+
                            (1|sample_tp1),
                            family = nbinom2(),
                            data=datsc)
loss_glmmTMB2.r2 <- glmmTMB(num_gene_losses ~ alpha_div_tp1*lag + read_count+ (1|species) + (1|subject),
                            family = nbinom2(),
                            data=datsc)
loss_glmmTMB2.r3 <- glmmTMB(num_gene_losses ~ alpha_div_tp1*lag + read_count+ (1|species) ,
                            family = nbinom2(),
                            data=datsc)
loss_glmmTMB2.r4 <- glmmTMB(num_gene_losses ~ alpha_div_tp1*lag + read_count ,
                            family = nbinom2(),
                            data=datsc)

anova(loss_glmmTMB2.r1,loss_glmmTMB2.r2,loss_glmmTMB2.r3,loss_glmmTMB2.r4)

#### test fixed effects
drop1(loss_glmmTMB2,.~ alpha_div_tp1*lag + read_count ,test='Chisq')

# test null model
null.los=glmmTMB(num_gene_losses ~ 1 + (1|species) + (1|subject)+
                 (1|sample_tp1)+ (1|sample_tp2),
                 family = nbinom2(),
                 data=datsc)
anova(loss_glmmTMB2,null.los)

sims_beta <- simulateResiduals(loss_glmmTMB2)
pdf('dharma_loss_glmmTMB2.1.pdf')
plot(sims_beta,quantreg = T)
dev.off()

###################
######### gene loss ~ richness
###################

datsc1=gene_all1

pvar1='richness_rare_tp1'
pvar2='lag'
pvar3='read_count'

datsc1[pvar1] <- lapply(gene_all1[pvar1],scale)
datsc1[pvar2] <- lapply(gene_all1[pvar2],scale)
datsc1[pvar3] <- lapply(gene_all1[pvar3],scale)


loss_rich_glmmTMB1 <- glmmTMB(num_gene_losses  ~ richness_rare_tp1*lag + read_count+ (1|species) + (1|subject)+
                              (1|sample_tp1)+
                              (1|sample_tp2),
                              family = nbinom1(),
                              data=datsc1)
loss_rich_glmmTMB2 <- glmmTMB(num_gene_losses ~ richness_rare_tp1*lag + read_count+ (1|species) + (1|subject)+
                              (1|sample_tp1)+ (1|sample_tp2),
                              family = nbinom2(),
                              data=datsc1)
loss_rich_glmmTMB.pois <- glmmTMB(num_gene_losses ~ richness_rare_tp1*lag + read_count+ (1|species) + (1|subject)+
                                  (1|sample_tp1)+
                                  (1|sample_tp2),
                                  family = poisson(),
                                  data=datsc1)

AICtab(loss_rich_glmmTMB.pois,loss_rich_glmmTMB2,loss_rich_glmmTMB1)

summary(loss_rich_glmmTMB2)

drop1(loss_rich_glmmTMB2,~richness_rare_tp1*lag + read_count,test='Chisq')

## test random efects
loss_rich_glmmTMB2.r1 <- glmmTMB(num_gene_losses ~ richness_rare_tp1*lag + 
                                       (1|species) + (1|subject)+
                                       (1|sample_tp1),
                                     family = nbinom2(),
                                     data=datsc1)
loss_rich_glmmTMB2.r2 <- glmmTMB(num_gene_losses ~ richness_rare_tp1*lag + 
                                       (1|species) + (1|subject),
                                     family = nbinom2(),
                                     data=datsc1)
loss_rich_glmmTMB2.r3 <- glmmTMB(num_gene_losses ~ richness_rare_tp1*lag + 
                                       (1|species) ,
                                     family = nbinom2(),
                                     data=datsc1)

loss_rich_glmmTMB2.r4 <- glmmTMB(num_gene_losses ~ richness_rare_tp1*lag ,
                                     family = nbinom2(),
                                     data=datsc1)
anova(loss_rich_glmmTMB2,loss_rich_glmmTMB2.r1,loss_rich_glmmTMB2.r2,loss_rich_glmmTMB2.r3,loss_rich_glmmTMB2.r4)


null.los.rich=glmmTMB(num_gene_losses ~ 1 + 
                        (1|species) + (1|subject)+
                        (1|sample_tp1)+ (1|sample_tp2),
                      family = nbinom2(),
                      data=datsc1)
anova(loss_rich_glmmTMB2,null.los.rich)

r2(loss_rich_glmmTMB2)
save(loss_rich_glmmTMB2,file='loss_rich_glmmTMB2.RData')

