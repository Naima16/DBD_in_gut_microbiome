##### GLMMs with strain count in a focal species as a function of diversity at high taxonomic levels
##### Naima Madi dec 2022

library(cowplot)  ## as.grob
library(jtools)   ## effect_plot
library(ggplot2)


## load glmm outputs : strain count ~ shannon at high taxonomic levels
load('/Users/naima/projet_HMP_gut_nandita_sept20/5juillet_final_version/Review_elife_mai_2022/strain_higher_taxo/nbinom2.class.RData')
load('/Users/naima/projet_HMP_gut_nandita_sept20/5juillet_final_version/Review_elife_mai_2022/strain_higher_taxo/nbinom2.fam.RData')
load('/Users/naima/projet_HMP_gut_nandita_sept20/5juillet_final_version/Review_elife_mai_2022/strain_higher_taxo/nbinom2.gen.RData')
load('/Users/naima/projet_HMP_gut_nandita_sept20/5juillet_final_version/Review_elife_mai_2022/strain_higher_taxo/nbinom2.phyla.RData')
load('/Users/naima/projet_HMP_gut_nandita_sept20/5juillet_final_version/Review_elife_mai_2022/strain_higher_taxo/r2.order.RData')

## load glmm outputs : strain count ~ richness at high taxonomic levels
load('/Users/naima/projet_HMP_gut_nandita_sept20/5juillet_final_version/Review_elife_mai_2022/strain_higherTaxo_richness/m3.RData')
load('/Users/naima/projet_HMP_gut_nandita_sept20/5juillet_final_version/Review_elife_mai_2022/strain_higherTaxo_richness/nbinom1.class.sing.RData')
load('/Users/naima/projet_HMP_gut_nandita_sept20/5juillet_final_version/Review_elife_mai_2022/strain_higherTaxo_richness/nbinom1.fam.RData')
load('/Users/naima/projet_HMP_gut_nandita_sept20/5juillet_final_version/Review_elife_mai_2022/strain_higherTaxo_richness/nbinom1.gen.RData')
load('/Users/naima/projet_HMP_gut_nandita_sept20/5juillet_final_version/Review_elife_mai_2022/strain_higherTaxo_richness/nbinom2.order.RData')


strain.phyl.shan=as_grob(effect_plot(nbinom2.phyla,'phyla_shannon',main.title= NULL,
                                     rug=FALSE,colors=pal2, col=1,y.label='Strain count',
                                     line.thickness=2,
                                     x.label='Shannon (phylum)',interval = T)+
                           theme_bw()+
                           theme(axis.title = element_text(size = rel(0.7)),
                                 axis.text = element_text(size = rel(0.5)),
                                 axis.title.y = element_blank())+
                           annotate('text',x = -2, y = 0.55,size=2, hjust = 0, 
                                    label = paste(" \nP =",'0.12190'),fontface = "plain"))

strain.clas.shan=as_grob(effect_plot(nbinom2.class,'class_shannon',main.title= NULL,
                                     rug=FALSE,colors=pal2, col=1,
                                     line.thickness=2,
                                     x.label='Shannon (class)',interval = T)+
                           theme_bw()+
                           theme(axis.title = element_text(size = rel(0.7)),
                                 axis.text = element_text(size = rel(0.5)),
                                 axis.title.y = element_blank())+
                           annotate('text',x = -2, y = 0.55,size=2, hjust = 0, 
                                    label = paste(" \nP =",'0.10333'),fontface = "plain"))

strain.ord.shan=as_grob(effect_plot(r2.order,'order_shannon',main.title= NULL,
                                    rug=FALSE,colors=pal2, col=1,
                                    line.thickness=2,
                                    x.label='Shannon (order)',interval = T)+
                          theme_bw()+
                          theme(axis.title = element_text(size = rel(0.7)),
                                axis.text = element_text(size = rel(0.5)),
                                axis.title.y = element_blank())+
                          annotate('text',x = -2, y = 0.55,size=2, hjust = 0, 
                                   label = paste(" \nP =",'0.0117'),fontface = "plain"))



strain.fam.shan=as_grob(effect_plot(nbinom2.fam,'family_shannon',main.title= NULL,
                                    rug=FALSE,colors=pal2, col=1,
                                    line.thickness=2,
                                    x.label='Shannon (family)',interval = T)+
                          theme_bw()+
                          theme(axis.title = element_text(size = rel(0.7)),
                                axis.text = element_text(size = rel(0.5)),
                                axis.title.y = element_blank())+
                          annotate('text',x = -2, y = 0.55,size=2, hjust = 0, 
                                   label = paste(" \nP =",'0.12952'),fontface = "plain"))

strain.gen.shan=as_grob(effect_plot(nbinom2.gen,'genus_shannon',main.title= NULL,
                                    rug=FALSE,colors=pal2, col=1,
                                    line.thickness=2,
                                    x.label='Shannon (genus)',interval = T)+
                          theme_bw()+
                          theme(axis.title = element_text(size = rel(0.7)),
                                axis.text = element_text(size = rel(0.5)),
                                axis.title.y = element_blank())+
                          annotate('text',x = -3, y = 0.6,size=2, hjust = 0, 
                                   label = paste(" \nP =",'0.00455'),fontface = "plain"))


strain.ph.rich=as_grob(effect_plot(m3,'phyla_nb',main.title= NULL,
                                   rug=FALSE,colors=pal2, col=1,#y.label='Strain count',
                                   line.thickness=2,
                                   x.label='Richness (phyla)',interval = T)+
                         theme_bw()+
                         theme(axis.title = element_text(size = rel(0.7)),
                               axis.text = element_text(size = rel(0.5)),
                               axis.title.y = element_blank())+
                         annotate('text',x = -5, y = 0.7,size=2, hjust = 0, 
                                  label = paste(" \nP =",'2.24e-05'),fontface = "plain"))

strain.cl.rich=as_grob(effect_plot(nbinom1.class.sing,'class_nb',main.title= NULL,
                                   rug=FALSE,colors=pal2, col=1,y.label='Strain count',
                                   line.thickness=2,
                                   x.label='Richness (class)',interval = T)+
                         theme_bw()+
                         theme(axis.title = element_text(size = rel(0.7)),
                               axis.text = element_text(size = rel(0.5)),
                               axis.title.y = element_blank())+
                         annotate('text',x = -5, y = 0.64,size=2, hjust = 0, 
                                  label = paste(" \nP =",'0.000221'),fontface = "plain"))

strain.ord.rich=as_grob(effect_plot(nbinom2.order,'order_nb',main.title= NULL,
                                    rug=FALSE,colors=pal2, col=1,y.label='Strain count',
                                    line.thickness=2,
                                    x.label='Richness (order)',interval = T)+
                          theme_bw()+
                          theme(axis.title = element_text(size = rel(0.7)),
                                axis.text = element_text(size = rel(0.5)),
                                axis.title.y = element_blank())+
                          annotate('text',x = -4, y = 0.77,size=2, hjust = 0, 
                                   label = paste(" \nP =",'2.67e-09'),fontface = "plain"))


strain.fam.rich=as_grob(effect_plot(nbinom1.fam,'family_nb',main.title= NULL,
                                    rug=FALSE,colors=pal2, col=1,y.label='Strain count',
                                    line.thickness=2,
                                    x.label='Richness (family)',interval = T)+
                          theme_bw()+
                          theme(axis.title = element_text(size = rel(0.7)),
                                axis.text = element_text(size = rel(0.5)),
                                axis.title.y = element_blank())+
                          annotate('text',x = -4, y = 0.65,size=2, hjust = 0, 
                                   label = paste(" \nP =",'2.66e-11'),fontface = "plain"))

strain.gen.rich=as_grob(effect_plot(nbinom1.gen,'genus_nb',main.title= NULL,
                                    rug=FALSE,colors=pal2, col=1,y.label='Strain count',
                                    line.thickness=2,
                                    x.label='Richness (genus)',interval = T)+
                          theme_bw()+
                          theme(axis.title = element_text(size = rel(0.7)),
                                axis.text = element_text(size = rel(0.5)),
                                axis.title.y = element_blank())+
                          annotate('text',x = -3, y = 0.7,size=2, hjust = 0, 
                                   label = paste(" \nP =",'6.09e-05'),fontface = "plain")) 


pdf('Figure_S3.pdf',width=7,height=3.5,pointsize = 10)
plot_grid(strain.phyl.shan,strain.clas.shan,strain.ord.shan,strain.fam.shan,strain.gen.shan,
          strain.ph.rich,strain.cl.rich,strain.ord.rich,strain.fam.rich,strain.gen.rich,
          ncol=5,nrow=2,labels = c('A1', 'B1','C1','D1','E1','A2', 'B2','C2','D2','E2'),label_size = 6)
dev.off()

