library(ggplot2)
library(wesanderson)
library(tidyverse)
library(effsize)
library(gridExtra)
library(RColorBrewer)
library(visdat)
# stability plots, plots of the data with the naive and deconfounded mean, 
# plots of the missingness

#' Examine the stability of z-statistics under different seeds:
results.all=NULL
for (i in 1:200) {
  load(paste0('./Results/noImputation/ic30_pc85_glm_gam_drtmle_seed',i,'.RData'))
  results.all = rbind(results.all,results.df)
}
rm(results.df)

# second set of seeds:
results.all2=NULL
for (i in 201:400) {
  load(paste0('./Results/noImputation/ic30_pc85_glm_gam_drtmle_seed',i,'.RData'))
  results.all2 = rbind(results.all2,results.df)
}
rm(results.df)

p=ggplot(data=results.all[order(results.all$z.stat.diff.naive),], aes(x=z.stat.diff.naive, y=z.stat.diff.SL, color=EdgeName))+
  geom_point() + theme(legend.position='none')+geom_abline(slope=1,intercept=0)+
  scale_color_manual(values = colorRampPalette(brewer.pal(n=8,name='Accent'))(153))

#+scale_color_manual(values = colorRampPalette(brewer.pal(n=12,name='Set3'))(153))

#library(wesanderson)
#+scale_color_manual(values = colorRampPalette(wes_palette(n=30,name='Zissou1', type='continuous'))(153))

pdf(file='./Application_Figures/Stability_Assess.pdf')
p
dev.off()

# Edges of interest from Lombardo:
# edgeList = c('r.ic2.ic17','r.ic8.ic17','r.ic17.ic19')
# results.all[results.all$EdgeName%in%edgeList,]
# three key edges are not significant

# average the z-statistics:
results.ave = results.all%>%
  group_by(EdgeName)%>%
  summarize(mean.ASD.naive=mean(mean.ASD.naive), mean.TD.naive=mean(mean.TD.naive),mean.diff.naive=mean(mean.diff.naive),
            z.stat.ASD.naive=mean(z.stat.ASD.naive),z.stat.TD.naive=mean(z.stat.TD.naive),z.stat.diff.naive=mean(z.stat.diff.naive),
            mean.ASD.SL=mean(mean.ASD.SL),mean.TD.SL=mean(mean.TD.SL),mean.diff.SL=mean(mean.diff.SL),
            z.stat.ASD.SL=mean(z.stat.ASD.SL),z.stat.TD.SL=mean(z.stat.TD.SL),z.stat.diff.SL=mean(z.stat.diff.SL))

results.ave2 = results.all2%>%
  group_by(EdgeName)%>%
  summarize(mean.ASD.naive=mean(mean.ASD.naive), mean.TD.naive=mean(mean.TD.naive),mean.diff.naive=mean(mean.diff.naive),
            z.stat.ASD.naive=mean(z.stat.ASD.naive),z.stat.TD.naive=mean(z.stat.TD.naive),z.stat.diff.naive=mean(z.stat.diff.naive),
            mean.ASD.SL=mean(mean.ASD.SL),mean.TD.SL=mean(mean.TD.SL),mean.diff.SL=mean(mean.diff.SL),
            z.stat.ASD.SL=mean(z.stat.ASD.SL),z.stat.TD.SL=mean(z.stat.TD.SL),z.stat.diff.SL=mean(z.stat.diff.SL))

# compare two estimates:
plot(results.ave$z.stat.diff.SL~results.ave2$z.stat.diff.SL)
head(cbind(results.ave$z.stat.diff.SL,results.ave2$z.stat.diff.SL))

#' Correlation of z statistics from two sets of 200 seeds across all edges (this result appears in the manuscript):
cor(results.ave$z.stat.diff.SL,results.ave2$z.stat.diff.SL)


###' Compare the selected edges at FDR=0.20 and 0.05 in the two sets of seeds:
results.ave$p.SL = 2*(1-pnorm(abs(results.ave$z.stat.diff.SL)))
results.ave$p.SL.fdr = p.adjust(results.ave$p.SL,method='BH')

results.ave2$p.SL = 2*(1-pnorm(abs(results.ave2$z.stat.diff.SL)))
results.ave2$p.SL.fdr = p.adjust(results.ave2$p.SL,method='BH')

#' These results appear in the manuscript:

sum(results.ave$p.SL.fdr<0.05)
sum(results.ave2$p.SL.fdr<0.05)
results.ave[results.ave$p.SL.fdr<0.05,c('EdgeName')]
results.ave2[results.ave2$p.SL.fdr<0.05,c('EdgeName')]
# same seven edges at fdr=0.05

sum(results.ave$p.SL.fdr<0.20)
sum(results.ave2$p.SL.fdr<0.20)
edges.20.set1 = results.ave[results.ave$p.SL.fdr<0.2,c('EdgeName','z.stat.diff.SL','p.SL.fdr')]
edges.20.set2 = results.ave2[results.ave2$p.SL.fdr<0.2,c('EdgeName','z.stat.diff.SL','p.SL.fdr')]

cbind(edges.20.set1$EdgeName,edges.20.set2$EdgeName)
all(edges.20.set1$EdgeName==edges.20.set2$EdgeName)

###########################
###########################
# Create results.ave from all seeds:
results.all.both = rbind(results.all,results.all2)
# standard errors:
results.all.both$se.diff.naive = results.all.both$mean.diff.naive/results.all.both$z.stat.diff.naive
results.all.both$se.diff.SL = results.all.both$mean.diff.SL/results.all.both$z.stat.diff.SL

results.ave = results.all.both%>%
  group_by(EdgeID, EdgeName)%>%
  summarize(mean.ASD.naive=mean(mean.ASD.naive), mean.TD.naive=mean(mean.TD.naive),mean.diff.naive=mean(mean.diff.naive),
            z.stat.ASD.naive=mean(z.stat.ASD.naive),z.stat.TD.naive=mean(z.stat.TD.naive),z.stat.diff.naive=mean(z.stat.diff.naive),
            mean.ASD.SL=mean(mean.ASD.SL),mean.TD.SL=mean(mean.TD.SL),mean.diff.SL=mean(mean.diff.SL),
            z.stat.ASD.SL=mean(z.stat.ASD.SL),z.stat.TD.SL=mean(z.stat.TD.SL),z.stat.diff.SL=mean(z.stat.diff.SL),
            se.diff.naive=mean(se.diff.naive),se.diff.SL=mean(se.diff.SL))



results.ave$p.naive = 2*(1-pnorm(abs(results.ave$z.stat.diff.naive)))
results.ave$p.naive.fdr = p.adjust(results.ave$p.naive,method='BH')
results.ave$p.SL = 2*(1-pnorm(abs(results.ave$z.stat.diff.SL)))
results.ave$p.SL.fdr = p.adjust(results.ave$p.SL,method='BH')

#' Number of edges indicated by DRTMLE at FDR=.20 when results averaged across all 400 seeds
sum(results.ave$p.SL.fdr<0.20)


#' Number of edges indicated by DRTMLE at FDR=.05 when results averaged across all 400 seeds
sum(results.ave$p.SL.fdr<0.05)


# Changes in ASD mean:
temp = results.ave[,c("EdgeName","mean.ASD.SL","mean.ASD.naive")]
temp$Change=temp$mean.ASD.SL-temp$mean.ASD.naive
results.ASD = temp%>%pivot_longer(cols=c(2,3), names_to = 'Method',values_to='mean.ASD')

p0 = results.ASD%>%ggplot(aes(x=Method, y=mean.ASD, group=EdgeName,color=Change))+
  geom_point()+geom_line()+ylim(c(-0.2,0.25))+ylab('Mean FC')+
  scale_x_discrete(labels=c("mean.ASD.naive" = "Naive", "mean.ASD.SL" = "DRTMLE"))+
  ggtitle('A) ASD')+scale_color_gradient2(limits=c(-0.005,0.005),oob=scales::squish)

# Changed in TD mean: 
temp = results.ave[,c("EdgeName","mean.TD.naive","mean.TD.SL")]
temp$Change=temp$mean.TD.SL-temp$mean.TD.naive

results.TD = temp%>%pivot_longer(cols=c(2,3), names_to = 'Method',values_to='mean.TD')
p1=results.TD%>%ggplot(aes(x=Method, y=mean.TD, group=EdgeName,color=Change))+
  geom_line()+geom_point()+ylim(c(-0.2,0.25))+ggtitle('B) TD')+
  scale_x_discrete(labels=c("mean.TD.naive" = "Naive", "mean.TD.SL" = "DRTMLE"))+
  ylab('Mean FC')+scale_color_gradient2(limits=c(-0.005,0.005))

# Difference in mean between ASD and TD: 
results.ave$diff.ASD.TD.naive = results.ave$mean.ASD.naive - results.ave$mean.TD.naive
results.ave$diff.ASD.TD.SL = results.ave$mean.ASD.SL - results.ave$mean.TD.SL

temp = results.ave[,c("EdgeName","diff.ASD.TD.naive","diff.ASD.TD.SL")]
temp$Change=temp$diff.ASD.TD.SL-temp$diff.ASD.TD.naive
results.diff = temp%>%pivot_longer(cols=c(2,3), names_to = 'Method',values_to='difference')

p2=results.diff%>%ggplot(aes(x=Method, y=difference, group=EdgeName, color=Change))+
  geom_line()+geom_point()+ylim(c(-0.015,0.015))+ggtitle('C) ASD-TD')+
  scale_x_discrete(labels=c("diff.ASD.TD.naive" = "Naive", "diff.ASD.TD.SL" = "DRTMLE"))+
  ylab('Mean ASD - TD')+scale_color_gradient2()

#pdf(file='~/Dropbox/Apps/Overleaf/MotionSelectionBias_rsfMRI/Figures/DeconfoundedGroupMeans.pdf',width=8,height=5)
grid.arrange(p0,p1,p2,ncol=3)
#dev.off()


#' Examine histograms of uncorrected p values from both approaches
par(mfrow=c(1,2))
p0=ggplot(results.ave,aes(x=p.naive))+geom_histogram()+ggtitle('A) P-values from naive')+ylim(0,20)+xlab('Naive')
p1=ggplot(results.ave,aes(x=p.SL))+geom_histogram()+ggtitle('B) P-values from DRTMLE')+ylim(0,20)+xlab('DRTMLE')

#pdf(file='~/Dropbox/Apps/Overleaf/MotionSelectionBias_rsfMRI/Figures/Pvalues_ASD_vs_TD.pdf',width=6,height=3)
grid.arrange(p0,p1,ncol=2)
#dev.off()

#' We observe more clustering of p values near 0 for DRTMLE

#' Examine histograms of standard errors for both approaches
par(mfrow=c(1,2))
hist(results.ave$se.diff.naive)
hist(results.ave$se.diff.SL)
mean(results.ave$se.diff.SL<results.ave$se.diff.naive)
#' Only slightly greater than 50% have smaller SEs

hist(results.ave$se.diff.SL/results.ave$se.diff.naive,main='Relative efficiency of SL to naive')

#' Calculate Cohen's D in naive estimates

#' Create plots to visualize the effect of drtmle on the deconfounded mean

#' Load propensities
propensities.all=NULL
for (seed in 1:400) {
  load(paste0('./Data/noImputation/DataWithPropensities_seed',seed,'.RData'))
  dat3$seed = rep(seed, nrow(dat3))
  dat3 <- select(dat3, c(ID, seed, Delta.KKI, propensities.glm, propensities.gam, propensities.SL, CompletePredictorCases))
  propensities.all = rbind(propensities.all, dat3)
}
rm(dat3)

seedtib <- tibble(propensities.all)

#' Calculate superlearner AUC for each seed
seedNest <- seedtib %>% 
  select(c(ID, seed, Delta.KKI, propensities.glm, propensities.gam, propensities.SL, CompletePredictorCases)) %>% 
  group_by(seed) %>% 
  filter(CompletePredictorCases==TRUE) %>% 
  tidyr::nest() %>% 
  mutate(auc = map(data, ~unlist(ROCit::rocit(score=.x$propensities.SL, class=.x$Delta.KKI, method = 'nonparametric')[6]))) %>% 
  unnest(auc)

# AUC for the first 5 seeds
seedNest$auc[1:5]

#' Summarize superlearner AUC across all seeds
auc.summary <- seedNest %>% 
  ungroup() %>% 
  summarise(min.SL.auc = min(auc),
            max.SL.auc = max(auc),
            mean.SL.auc = mean(auc))

auc.summary

# re-load seed 1 to get additional information that is the same for all seeds 
load('./Data/noImputation/DataWithPropensities_seed1.RData')

#' Compare to AUC for glm
summary(ROCit::rocit(score=dat3$propensities.glm, class=dat3$Delta.KKI, method = 'nonparametric'))

#' Compare to AUC for gam
summary(ROCit::rocit(score=dat3$propensities.gam, class=dat3$Delta.KKI, method = 'nonparametric'))

#' Find min propensity for each seed
propensities.summary  <- propensities.all %>%
  group_by(seed) %>% 
  summarise(min.SL.propensity = min(propensities.SL, na.rm = TRUE))

#min superlearner propensity for the first 5 seeds
propensities.summary$min.SL.propensity[1:5]

#' What is the average smallest propensity across seeds?
mean(propensities.summary$min.SL.propensity)
  
#' Average propensities across seeds for each participant
propensities.ave = propensities.all%>%
  group_by(ID)%>%
  summarize(mean.SL.propensity = mean(propensities.SL))

# merge average propensities with modified partial correlations
dat3 <- merge(dat3, propensities.ave)


nEdges=153
idx.pass.cc = dat3$KKI_criteria=='Pass' & !is.na(dat3$propensities.SL) & !is.na(dat3$r.ic1.ic2)
idx.pass.cc.asd = idx.pass.cc & dat3$PrimaryDiagnosis=='Autism'
idx.pass.cc.td = idx.pass.cc & dat3$PrimaryDiagnosis=='None'

#' check cohen's d calculation:
n1 = sum(idx.pass.cc.asd)
n2 = sum(idx.pass.cc.td)
(mean(dat3[idx.pass.cc.asd,'r.ic2.ic27'])-mean(dat3[idx.pass.cc.td,'r.ic2.ic27']))/sqrt(((n1-1)*var(dat3[idx.pass.cc.asd,'r.ic2.ic27'])+(n2-1)*var(dat3[idx.pass.cc.td,'r.ic2.ic27']))/(n1+n2-2))
cohen.d(dat3[idx.pass.cc.asd,'r.ic2.ic27'],dat3[idx.pass.cc.td,'r.ic2.ic27'])
# they match


startEdgeidx=which(names(dat3)=='r.ic1.ic2')

cohen.d.df = data.frame('EdgeID'=numeric(nEdges),'EdgeName'=numeric(nEdges),'cohensd'=numeric(nEdges),'adj.cohensd'=numeric(nEdges))

#' NOTE: Adjusted cohen's d is a crude adjustment to cohen's d based on the new means estimated from drtmle. 
#' It does not update the variance estimates, which is a novel problem. 
#' The adjusted cohen's d is not used in the manuscript,
#' but I did it here to gain insight into the magnitude of the change in mean:
for (edgeidx in 1:nEdges) {
  cohen.d.df[edgeidx,'EdgeID'] = edgeidx
  dat3.edgeidx = startEdgeidx+edgeidx-1 
  cohen.d.df[edgeidx,'EdgeName'] = names(dat3)[dat3.edgeidx]
  temp = cohen.d(dat3[idx.pass.cc.asd,dat3.edgeidx],dat3[idx.pass.cc.td,dat3.edgeidx])
  cohen.d.df[edgeidx,'cohensd']=temp$estimate
  pooled.sd = sqrt(((n1-1)*var(dat3[idx.pass.cc.asd,dat3.edgeidx])+(n2-1)*var(dat3[idx.pass.cc.td,dat3.edgeidx]))/(n1+n2-2))
  cohen.d.df[edgeidx,'adj.cohensd']=results.ave[edgeidx,'mean.diff.SL']/pooled.sd
}

cohen.d.df
# check that the indexing is correct (that the correct drtmle
# estimate was used):
with(cohen.d.df,cor(cohensd,adj.cohensd))
# high correlation (~0.90) so looks like calcs are correct

max(abs(cohen.d.df$adj.cohensd))
#' We see increased effect size but this is a crude estimate that only adjusts the means

cohen.d.df[cohen.d.df$EdgeName=='r.ic2.ic27',]

#' Proportion of edge effect sizes that "increased":
mean(abs(cohen.d.df$cohensd)<abs(cohen.d.df$adj.cohensd))
sum(abs(cohen.d.df$cohensd)<abs(cohen.d.df$adj.cohensd))

#' There was a tendency to increase
#' In the absence of any changes in SEs, we see the effect sizes increase, indicating a correction in bias. 


#############################################
#' Ordered edges with smallest p-values in drtmle:
(list.edges = results.ave$EdgeName[results.ave$p.SL.fdr<0.20])
list.pvalues = results.ave$p.SL[results.ave$p.SL.fdr<0.20]
list.edges = list.edges[order(list.pvalues)]
list.pvalues = sort(list.pvalues)
# naive cohen's d of selected edges:
cohen.d.df[cohen.d.df$EdgeName%in%list.edges,]


#' These values appear in the manuscript:
max(abs(cohen.d.df$cohensd))
#' average of the absolute value of Cohen's D among selected edges:
mean(abs(cohen.d.df$cohensd[cohen.d.df$EdgeName%in%list.edges]))

### this is how much the effect sizes would change due to mean only,
# do not use this, but it provides insight into mean versus variance:
max(abs(cohen.d.df$adj.cohensd))
# average of the absolute value of Cohen's D among selected edges:
mean(abs(cohen.d.df$adj.cohensd[cohen.d.df$EdgeName%in%list.edges]))

#' Examine the relative contributions of changes in means versus changes in standard errors:
examine_means_se_selectededges = results.ave[results.ave$EdgeName%in%list.edges,c('EdgeName','mean.diff.naive',
                                                                                  'mean.diff.SL','se.diff.naive',
                                                                                  'se.diff.SL', "p.SL.fdr", "p.naive.fdr")]

#' number of selected edges in which absolute mean difference increased among selected edges:
sum(abs(examine_means_se_selectededges$mean.diff.SL)>abs(examine_means_se_selectededges$mean.diff.naive))

#' number of  edges selected by DRTMLE but not the naive approach in which absolute mean difference increased:
sum(abs(examine_means_se_selectededges$mean.diff.SL)>abs(examine_means_se_selectededges$mean.diff.naive) &
      examine_means_se_selectededges$p.naive.fdr>.2)

#' number of edges selected by DRTMLE but not by the naive approach in which se decreased:
sum(examine_means_se_selectededges$se.diff.SL<examine_means_se_selectededges$se.diff.naive & 
      examine_means_se_selectededges$p.naive.fdr>.2)

#' number ofedges selected by DRTMLE but not the naive approach in which absolute mean difference increased among selected edges:
sum(examine_means_se_selectededges$se.diff.SL<examine_means_se_selectededges$se.diff.naive)

#' number of selected edges in which both abs mean diff increased and se decreased:
sum(examine_means_se_selectededges$se.diff.SL<examine_means_se_selectededges$se.diff.naive & 
      abs(examine_means_se_selectededges$mean.diff.SL)>abs(examine_means_se_selectededges$mean.diff.naive))

#' number of edges selected by DRTMLE but not by the naive approach in which both abs mean diff increases and se decreased:
sum(examine_means_se_selectededges$se.diff.SL<examine_means_se_selectededges$se.diff.naive & 
      abs(examine_means_se_selectededges$mean.diff.SL)>abs(examine_means_se_selectededges$mean.diff.naive) &
      examine_means_se_selectededges$p.naive.fdr>.2)


#' number of edges in which abs mean difference increase among 153 edges:
sum(abs(results.ave$mean.diff.SL)>abs(results.ave$mean.diff.naive))

#' number of edges in which se decrease among 153 edges:
sum(results.ave$se.diff.SL<results.ave$se.diff.naive)

#' number of edges in which abs mean difference increase and se decreased among 153 edges:
sum(results.ave$se.diff.SL<results.ave$se.diff.naive & abs(results.ave$mean.diff.SL)>abs(results.ave$mean.diff.naive))


examine_means_se_selectededges$abs.mean.diff = abs(examine_means_se_selectededges$mean.diff.SL)-abs(examine_means_se_selectededges$mean.diff.naive)
examine_means_se_selectededges$se.rel.eff = examine_means_se_selectededges$se.diff.SL/examine_means_se_selectededges$se.diff.naive

examine_means_se_selectededges$both_sig=as.factor(examine_means_se_selectededges$p.naive.fdr<=.2)

vizMeanSE = ggplot(examine_means_se_selectededges, aes(x=abs.mean.diff, y=se.rel.eff, 
                                                       shape=both_sig))+
  geom_point(size=4, alpha = .8)+
  scale_shape_discrete(labels = c("naive FDR p <= .2", "naive FDR p > .2"))+
  geom_hline(yintercept = 1, linetype = 2)+
  geom_vline(xintercept = 0, linetype = 2)+
  theme_light()+theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size =12))+
  ylab("DRTMLE/Naive ASD-TD SE")+
  xlab("DRTMLE-Naive (ASD-TD mean)")

pdf('./ReviewerResponse/compare_meansSE.pdf')
vizMeanSE
dev.off()


#####################################
# Visualize propensity scores for a few edges:
# these figures do not appear in the manuscript, as the 
# propensities did not appear to shed insight into drtmle:
pal <- wes_palette("Zissou1", 100, type = "continuous")
#library(colorRamps)
#library(gridExtra)
#pal = matlab.like(100)
set.seed(123)
plot_pcorr_fun_propensities = function(EdgeName,EdgeNamePlot,legend=FALSE) { 
      subdata = dat3[idx.pass.cc,c('PrimaryDiagnosis',EdgeName,'mean.SL.propensity')]
      names(subdata)[2] = 'EdgeName'
      temp = results.ave[results.ave$EdgeName==EdgeName,c('mean.ASD.naive','mean.TD.naive')]
      temp.naive = data.frame('PrimaryDiagnosis'=c('Autism','None'),'EdgeName'=c(temp$mean.ASD.naive,temp$mean.TD.naive))
      temp = results.ave[results.ave$EdgeName==EdgeName,c('mean.ASD.SL','mean.TD.SL')]
      temp.SL = data.frame('PrimaryDiagnosis'=c('Autism','None'),'EdgeName'=c(temp$mean.ASD.SL,temp$mean.TD.SL))
      if(legend==FALSE) {
    outplot = ggplot(subdata, aes(x=PrimaryDiagnosis, y=EdgeName))+
      geom_jitter(shape=16,position=position_jitter(0.2),size=2,aes(color=1/mean.SL.propensity))+
      scale_colour_gradientn(colours=pal,oob = scales::squish,name = expression(1/p[i]),limits=c(1,2.1))+
      xlab("Primary Diagnosis")+ylab(EdgeNamePlot)+geom_point(data=temp.naive,shape=15,size=2,color="red")+
      geom_point(data=temp.SL,shape=15,size=2,color='blue',alpha=0.5)+theme(legend.position="none") 
    } else {
        outplot = ggplot(subdata, aes(x=PrimaryDiagnosis, y=EdgeName))+ 
          geom_jitter(shape=16,position=position_jitter(0.2),size=2,aes(color=1/mean.SL.propensity))+
          scale_colour_gradientn(colours=pal,oob = scales::squish,name = expression(1/p[i]),limits=c(1,2.1))+
          xlab("Primary Diagnosis")+ylab(EdgeNamePlot)+geom_point(data=temp.naive,shape=15,size=2,color="red")+
          geom_point(data=temp.SL,shape=15,size=2,color='blue',alpha=0.5)}
    outplot
  }

#' list the selected edges ordered by their uncorrected pvalues:
list.edges

#' manually enter the desired edges:
gn.p0=plot_pcorr_fun_propensities(EdgeName='r.ic2.ic27',EdgeNamePlot='IC02-IC27')
gn.p1=plot_pcorr_fun_propensities(EdgeName='r.ic14.ic21',EdgeNamePlot='IC14-IC21')
gn.p2=plot_pcorr_fun_propensities(EdgeName='r.ic19.ic26',EdgeNamePlot='IC19-IC26')
gn.p3=plot_pcorr_fun_propensities(EdgeName='r.ic4.ic17',EdgeNamePlot='IC04-IC17')
gn.p4=plot_pcorr_fun_propensities(EdgeName='r.ic14.ic19',EdgeNamePlot='IC14-IC19')
gn.p5=plot_pcorr_fun_propensities(EdgeName='r.ic13.ic26',EdgeNamePlot='IC13-IC26')
gn.p6=plot_pcorr_fun_propensities(EdgeName='r.ic1.ic21',EdgeNamePlot='IC01-IC21')
gn.p7=plot_pcorr_fun_propensities(EdgeName='r.ic17.ic24',EdgeNamePlot='IC17-IC24')
gn.p8=plot_pcorr_fun_propensities(EdgeName='r.ic19.ic21',EdgeNamePlot='IC19-IC21')

pdf(file='./ReviewerResponse/DataForNineComponentsWithPropensities_naive_drtmle.pdf')
grid.arrange(gn.p0,gn.p1,gn.p2,gn.p3,gn.p4,gn.p5,gn.p6,gn.p7,gn.p8,nrow=3)
dev.off()


###################
##################
#' Create a simplified version of the previous plots: DRTMLE diagnosis group means in red, naive group means in blue; modified partial correlations in grey

(slOnly_list.edges = results.ave$EdgeName[results.ave$p.SL.fdr<0.20 & results.ave$p.naive.fdr>0.2])
slOnly_list.pvalues = results.ave$p.SL[results.ave$p.SL.fdr<0.20 & results.ave$p.naive.fdr>0.2]
slOnly_list.edges = list.edges[order(slOnly_list.pvalues)]
slOnly_list.pvalues = sort(slOnly_list.pvalues)

slOnly_list.edges

set.seed(123)
plot_pcorr_fun = function(EdgeName,EdgeNamePlot,legend=FALSE) { 
  subdata = dat3[idx.pass.cc,c('PrimaryDiagnosis',EdgeName,'propensities.SL')]
  names(subdata)[2] = 'EdgeName'
  temp = results.ave[results.ave$EdgeName==EdgeName,c('mean.ASD.naive','mean.TD.naive')]
  temp.naive = data.frame('PrimaryDiagnosis'=c('Autism','None'),'EdgeName'=c(temp$mean.ASD.naive,temp$mean.TD.naive))
  temp = results.ave[results.ave$EdgeName==EdgeName,c('mean.ASD.SL','mean.TD.SL')]
  temp.SL = data.frame('PrimaryDiagnosis'=c('Autism','None'),'EdgeName'=c(temp$mean.ASD.SL,temp$mean.TD.SL))
  if(legend==FALSE) {
    outplot = ggplot(subdata, aes(x=PrimaryDiagnosis, y=EdgeName))+ geom_jitter(shape=16,position=position_jitter(0.2),size=2,color='grey50')+xlab("Primary Diagnosis")+ylab(EdgeNamePlot)+geom_point(data=temp.naive,shape=15,size=2,color="red")+geom_point(data=temp.SL,shape=15,size=2,color='blue',alpha=0.5)+theme(legend.position="none") 
  } else {
    outplot = ggplot(subdata, aes(x=PrimaryDiagnosis, y=EdgeName))+ geom_jitter(shape=16,position=position_jitter(0.2),size=2,color='grey50')+xlab("Primary Diagnosis")+ylab(EdgeNamePlot)+geom_point(data=temp.naive,shape=15,size=2,color="red")+geom_point(data=temp.SL,shape=15,size=2,color='blue',alpha=0.5)}
  outplot
}


gn.p0=plot_pcorr_fun(EdgeName='r.ic14.ic21',EdgeNamePlot='IC14-IC21')
gn.p1=plot_pcorr_fun(EdgeName='r.ic1.ic14',EdgeNamePlot='IC01-IC14')
gn.p2=plot_pcorr_fun(EdgeName='r.ic19.ic26',EdgeNamePlot='IC19-IC26')
gn.p3=plot_pcorr_fun(EdgeName='r.ic17.ic24',EdgeNamePlot='IC17-IC24')
gn.p4=plot_pcorr_fun(EdgeName='r.ic14.ic19',EdgeNamePlot='IC14-IC19')
gn.p5=plot_pcorr_fun(EdgeName='r.ic2.ic27',EdgeNamePlot='IC02-IC27')
gn.p6=plot_pcorr_fun(EdgeName='r.ic13.ic26',EdgeNamePlot='IC13-IC26')
gn.p7=plot_pcorr_fun(EdgeName='r.ic21.ic30',EdgeNamePlot='IC21-IC30')
gn.p8=plot_pcorr_fun(EdgeName='r.ic8.ic22',EdgeNamePlot='IC08-IC22',legend=TRUE)

#' Print Figure S.4: Plots of partial correlations, naive means, and DRTMLE means for each diagnosis group for the nine components with the smallest DRTMLE p-values
pdf(file='./Application_figures/DataForNineComponents_naive_drtmle.pdf')
grid.arrange(gn.p0,gn.p1,gn.p2,gn.p3,gn.p4,gn.p5,gn.p6,gn.p7,gn.p8,nrow=3)
dev.off()

#' Examine means for the edge with the largest change:
results.ave$mag_change.ASD = abs(results.ave$mean.ASD.SL - results.ave$mean.ASD.naive)
results.ave$EdgeName[which.max(results.ave$mag_change.ASD)]
plot_pcorr_fun(EdgeName='r.ic24.ic27',EdgeNamePlot='IC24-IC27')


