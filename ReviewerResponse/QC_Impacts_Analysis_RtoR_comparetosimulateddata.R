# BRisk
#' This code examines some aspects of the real data for comparison with the simulated data. 
#' It uses the revised dataset that has no imputation. 
# This code is modified to be used locally or on a cluster. 
#seedID = 1;

local=TRUE
# Change this for cluster or local:

if(local){
    #setwd('~/Dropbox/QualityControlImpactsFMRI')
    save.input.data = FALSE
    getOption("mc.cores")
    options(mc.cores=8)
    seed=1
  } else {
    setwd('~/risk_share/QualityControlImpacts')
    save.input.data = FALSE
    options(mc.cores=1)
    seed = seedID
}
 
a = .libPaths()
.libPaths(c('/home/benjamin.risk/Rlibs',a))

.libPaths()

getOption("mc.cores")
set.seed(seed, "L'Ecuyer-CMRG")

#+ load-packages, echo=FALSE, warnings=FALSE
tic = proc.time()
library(drtmle)
library(SuperLearner)
library(earth)
library(quadprog)
library(gam)
library(nloptr)
library(future)
library(future.apply)
library(xgboost)
library(ranger)
library(visdat)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(e1071)
library(glmnet)
library(readxl)
library(ROCit)
# KKI: Lenient
# Ciric: Strict

#+ load-data
dat=read.csv('../Data/Master_HeadMotion_wholeGroup_partialCorrelations_ic30_20210803.csv',header=T)

# sort by ID: the simplifies the later indexing, as a merge with 
# propensities sorts by ID:
dat = dat[order(dat$ID),]

table(dat$PrimaryDiagnosis)
table(dat$Race,dat$Sex,useNA='always') # Used in Protection of Human Subjects table in grant

# note: create a "caucasian/other" race variable that allocates the four children without race to caucasian
dat$Race2 = dat$Race
dat$Race2[dat$Race2%in%c('Hispanic','Other','Unknown')]='Caucasian'
table(dat$Race2)


# Add ADHD variable:
dat$ADHD_Secondary = ifelse(dat$ADHD.Subtype=='No dx',0,1)
table(dat$ADHD_Secondary,dat$ADHD.Subtype,useNA='always') # there should be no strange categories. 2/3/2021: Looks good
table(dat$PrimaryDiagnosis[dat$KKI_criteria=='Pass'])
table(dat$PrimaryDiagnosis[dat$Ciric_length=='Pass'])

# subset to signal components:
#ic_class = read_excel('./DeconfoundedFMRI/componentLabels_pca85_ica30.xlsx')
ic_class = read_excel('../componentLabels_pca85_ica30.xlsx')

# create names we want to exclude:
artifacts = c(1:nrow(ic_class))[ic_class$signal==0]
badICs = paste0('ic',artifacts)
allICs = paste0('ic',1:nrow(ic_class))
badEdges=NULL
for (i in 1:nrow(ic_class)) {
  for (j in 1:length(badICs)) {
    badEdges = c(badEdges,paste0(allICs[i],'.',badICs[j]))
  }
}
for (i in 1:length(badICs)) {
  for (j in 1:nrow(ic_class)) {
    badEdges = c(badEdges,paste0(badICs[i],'.',allICs[j]))
  }
}

dat2 = dat[,!names(dat)%in%badEdges]

# Variables containing Fisher z partial correlations 

startEdgeidx = which(names(dat2)=='ic1.ic2')
endEdgeidx = which(names(dat2)=='ic29.ic30')
names(dat2)[startEdgeidx:endEdgeidx] #18*17/2 = nEdges
(nEdges = endEdgeidx - startEdgeidx+1)

# note: check whether all "PASS" have correlations. Should be all FALSE:
table(is.na(dat2[,startEdgeidx]) & dat2$KKI_criteria=='Pass')
# No missing. Cool.

# Motion variables:
# MeanFramewiseDisplacement.KKI
# MaxFramewiseDisplacement.KKI
# FramesWithFDLessThanOrEqualTo250microns
# Here, we control for sex and motion effects. We do not control for age, but 
# age is included in the propensity and outcome models, resulting
# in the age distribution for each group (including fails), which is 
# approximately equal between groups. 

# Create a variable r.ic1.ic2 for each signal pairing that 
# contains the residuals from the three motion variables
temp = dat2[,c('PrimaryDiagnosis','MeanFramewiseDisplacement.KKI','MaxFramewiseDisplacement.KKI',"FramesWithFDLessThanOrEqualTo250microns","Sex","Race2","SES.Family","ic1.ic2")]
completeCases = complete.cases(temp)

#nrow(temp[is.na(temp$ic1.ic2) ,])

t.values.lm=NULL
t.values.naive=NULL

lm.variables=c('MeanFramewiseDisplacement.KKI')
# NOTE: for better interpretation of the plot of mean changes, center all variables:
for (i in c(startEdgeidx:endEdgeidx)) {
  model.temp = lm(dat2[,i]~scale(dat2$MeanFramewiseDisplacement.KKI,center = TRUE, scale=FALSE)+scale(dat2$MaxFramewiseDisplacement.KKI,center=TRUE,scale=FALSE)+scale(dat2$FramesWithFDLessThanOrEqualTo250microns,center=TRUE,scale=FALSE)+dat2$Sex+dat2$Race2+scale(dat2$SES.Family,center=TRUE,scale=FALSE)+dat2$PrimaryDiagnosis)
  # create the motion-adjusted fconn data:
  dat2[completeCases,paste0('r.',names(dat2)[i])]=residuals(model.temp)+coef(model.temp)["(Intercept)"]+coef(model.temp)["dat2$PrimaryDiagnosisNone"]*(dat2$PrimaryDiagnosis[completeCases]=='None')
    
    # Audits:
    # check ordering is same:
    trash.audit = residuals(model.temp)+fitted(model.temp)
    trash = abs(dat2[completeCases,i]-trash.audit)<1e-16
    if(!all(trash)) stop('Ordering mismatch -- check variables in complete cases')     # check that t-statistic on residuals is approximately equal to t-statistic from lm:
    t.values.lm = c(t.values.lm,coef(summary(model.temp))['dat2$PrimaryDiagnosisNone',3])
    t.values.naive =   c(t.values.naive,t.test(dat2[completeCases & dat2$PrimaryDiagnosis=='None',paste0('r.',names(dat2)[i])],dat2[completeCases & dat2$PrimaryDiagnosis=='Autism',paste0('r.',names(dat2)[i])])$statistic)
}

# Make all ADOS in TD = 0
dat2$ADOS.Comparable.Total[dat2$PrimaryDiagnosis=='None'] = 0

# specify learners for super learner. gn is the propensity model and Qbar is the outcome model,
# also used for imputing PANESS: 
my.SL.libs.gn= c("SL.earth","SL.glmnet","SL.gam","SL.glm","SL.ranger","SL.step","SL.step.interaction","SL.xgboost","SL.mean")
my.SL.libs.Qbar= c("SL.earth","SL.glmnet","SL.gam","SL.glm","SL.ranger","SL.ridge","SL.step","SL.step.interaction","SL.svm","SL.xgboost","SL.mean")


#<---------------------------------

# Dataset for propensity model:
gn.variables = c('KKI_criteria','PrimaryDiagnosis','ADHD_Secondary','AgeAtScan','handedness','CurrentlyOnStimulants','PANESS.TotalOverflowNotAccountingForAge','WISC.GAI','DuPaulHome.InattentionRaw','DuPaulHome.HyperactivityRaw','ADOS.Comparable.Total')

# these indices will be used in the outcome model and drtmle as well:
temp.data = dat2[,c(gn.variables)]

# note: include variables from the linear model here to define a consistent set of
# complete predictor cases for the propensity and outcome models
idx.all.cc = complete.cases(temp.data) & complete.cases(dat2[,c('Sex','SES.Family','Race2')])
temp.data = temp.data[idx.all.cc,]
#temp.data$AgeAtScanXdx = (temp.data$PrimaryDiagnosis=='Autism')*temp.data$AgeAtScan
#AgeAtScan
gn.xmat = data.frame(model.matrix(KKI_criteria~.,data=temp.data)[,-1])

# Create a variable equal to one if the observation is used:
# NOTE: due to missing observations in variables in the initial linear
# model, need to include linear model variables in this: 
dat2$CompletePredictorCases = idx.all.cc
sum(dat2$CompletePredictorCases)


#' For RtoR describing simulations section-------------------------------->
#+ compare-sim-real-prop-rejected
Delta.KKI = ifelse(temp.data$KKI_criteria=='Pass',1,0)
mean(Delta.KKI[gn.xmat$PrimaryDiagnosisNone==0])
mean(Delta.KKI[gn.xmat$PrimaryDiagnosisNone==1])

#+ compare-sim-real-ADOS-usability
model.gam = gam(Delta.KKI~s(ADOS.Comparable.Total),data=gn.xmat,family=binomial)
summary(model.gam)

mean(gn.xmat$ADOS.Comparable.Total[gn.xmat$PrimaryDiagnosisNone==0])

newd <- gn.xmat[1, ] # grab any row; we are going to change temperature only
newd$ADOS.Comparable.Total <- mean(gn.xmat$ADOS.Comparable.Total[gn.xmat$PrimaryDiagnosisNone==0]) - 1e-05 # subtract some small number
y1 <- predict(model.gam, newd)
newd$ADOS.Comparable.Total <- mean(gn.xmat$ADOS.Comparable.Total[gn.xmat$PrimaryDiagnosisNone==0]) + 1e-05 # add some small number
y2 <- predict(model.gam, newd)

#' slope at mean ADOS in real data (full sample)
(y2 - y1)/2e-05

#' effect at mean ADOS in real data (full sample)
mean(gn.xmat$ADOS.Comparable.Total[gn.xmat$PrimaryDiagnosisNone==0])*(y2 - y1)/2e-05

newd <- gn.xmat[1, ] # grab any row; we are going to change temperature only
newd$ADOS.Comparable.Total <- 20 - 1e-05 # subtract some small number
y1 <- predict(model.gam, newd)
newd$ADOS.Comparable.Total <- 20 + 1e-05 # add some small number
y2 <- predict(model.gam, newd)
(y2 - y1)/2e-05

## <-----------------------------------------

(propensity.KKI = mcSuperLearner(Y = Delta.KKI, X = gn.xmat, family=binomial(link='logit'),SL.library = my.SL.libs.gn, cvControl = list(V = 10), method='method.CC_nloglik')) # 10-fold CV

propensities.SL = propensity.KKI$SL.predict

# merge propensities back to dataset:
prop.asdtd = data.frame('ID'=dat2$ID[idx.all.cc],'propensities.SL'=propensities.SL,Delta.KKI)
dat3 = merge(dat2,prop.asdtd,all = TRUE) # here, we keep all observations

# this merge changes the ordering of ID. Hence, new index vectors are created,
# and new xmat need to be made even when using the same variables in the 
# propensity and outcome models

# save datasets to be loaded for DRTMLE:
# NOTE: These datasets include the propensities, which change with each seed:
#if (save.input.data) {
#  save(file=paste0('./Data/DataWithPropensities_seed',seed,'.RData'),dat3)
#}

# check that the merge is correct:
table(dat3$KKI_criteria,dat3$Delta.KKI)
# these should be all in agreement (0 on off diagonal)
all(idx.all.cc==!is.na(dat3$propensities.SL))

  # Delta should be correlated with propensities: 
  cor(1*(dat3$KKI_criteria=='Pass'),dat3$propensities.SL,use='pairwise.complete.obs')

#' Create outcome regression datasets:
#' NOTE: using same variables for propensity and outcome model. 
Qn.variables =  gn.variables
#NOTE: KKI_criteria is not used in prediction, but is included as a trick to construct the design matrix

# complete cases pass defined by pass, no missing propensities (behavioral variables), and no missing fconn (driven by the initial linear model with imbalanced variables)
idx.pass.cc = dat3$KKI_criteria=='Pass' & !is.na(dat3$propensities.SL)
# use in naive estimates in for loop:
idx.pass.cc.asd = idx.pass.cc & dat3$PrimaryDiagnosis=='Autism'
idx.pass.cc.td = idx.pass.cc & dat3$PrimaryDiagnosis=='None'

# used in drtmle:
idx.all.cc.asd = idx.all.cc & dat3$PrimaryDiagnosis == 'Autism'
idx.all.cc.td = idx.all.cc & dat3$PrimaryDiagnosis == 'None'

# These datasets necessary to run SuperLearner without error
# (note: ordering of indices differs from the propensity models due to the merge with the propensities)
temp.data = dat3[idx.pass.cc,Qn.variables]
Qn.xmat.fit = data.frame(model.matrix(KKI_criteria~.,data=temp.data))[,-1]

temp.data = dat3[idx.all.cc,c('ID',Qn.variables)]
Qn.xmat.predict = data.frame(model.matrix(KKI_criteria~.,data=temp.data))[,-1]

# Separate ASD and TD datasets are necessary to obtain the DRTMLE estimates:
temp.data = Qn.xmat.predict[Qn.xmat.predict$PrimaryDiagnosisNone==0,]
# check that the order matches idx.all.cc.asd:
all(temp.data$ID==dat3$ID[idx.all.cc.asd])
Qn.xmat.predict.asd = data.frame(model.matrix(numeric(nrow(temp.data))~.,data=temp.data)[,-c(1,2)])

temp.data = Qn.xmat.predict[Qn.xmat.predict$PrimaryDiagnosisNone==1,]
# check that the order matches idx.all.cc.td:
all(temp.data$ID==dat3$ID[idx.all.cc.td])
Qn.xmat.predict.td = data.frame(model.matrix(numeric(nrow(temp.data))~.,data=temp.data))[,-c(1,2)]


## FOR RtoR----->

par(mfrow=c(1,2))
hist(Qn.xmat.predict.asd$ADOS.Comparable.Total,main='a) ASD severity in real data (full sample)', col="#7FAE88")
abline(v=mean(Qn.xmat.predict.asd$ADOS.Comparable.Total),col='red', lwd=3)
#' Mean ADOS severity in the real data (full sample)
mean(Qn.xmat.predict.asd$ADOS.Comparable.Total)
# usable and unusable ADOS:
quantile(Qn.xmat.predict$ADOS.Comparable.Total[Qn.xmat.predict$PrimaryDiagnosisNone==0])

hist(Qn.xmat.fit$ADOS.Comparable.Total[Qn.xmat.fit$PrimaryDiagnosisNone==0], main='B) ASD severity in real data (usable)', col="#5D60AA")
abline(v=mean(Qn.xmat.fit$ADOS.Comparable.Total[Qn.xmat.fit$PrimaryDiagnosisNone==0]),col='red', lwd=3) 
#' Mean ADOS severity in the real data (usable after lenient motion QC)
mean(Qn.xmat.fit$ADOS.Comparable.Total[Qn.xmat.fit$PrimaryDiagnosisNone==0])
# usable ADOS:
quantile(Qn.xmat.fit$ADOS.Comparable.Total[Qn.xmat.fit$PrimaryDiagnosisNone==0])


#' Examine the range of correlations for ADOS and the range of intercepts for ASD:
startEdgeidx=which(names(dat3)=='r.ic1.ic2')

p.value.ados = rep(NA,nEdges)
cor.ados=rep(NA,nEdges)
for (edgeidx in 1:nEdges) {
  dat3.edgeidx = startEdgeidx+edgeidx-1 
  cor.ados[edgeidx]=cor(dat3[idx.pass.cc,dat3.edgeidx],dat3$ADOS.Comparable.Total[idx.pass.cc])
  model.gam.out = gam(dat3[idx.pass.cc,dat3.edgeidx]~s(ADOS.Comparable.Total)+PrimaryDiagnosis,data=dat3[idx.pass.cc,])
  plot(model.gam.out,terms = 's(ADOS.Comparable.Total)')
}

max(cor.ados)
min(cor.ados)

## <---------------------
# note on power:
# the wilcoxon power does not vary with the numbers in each group,
# which is a bummer.
# here is an example of how the group size should matter:
library(pwr)
pwr.t2n.test(n1 = 20, n2=80, d = 0.5, sig.level = 0.05, power = NULL,alternative = c("two.sided"))
pwr.t2n.test(n1 = 50, n2=50, d = 0.5, sig.level = 0.05, power = NULL,alternative = c("two.sided"))

