# brisk
# This code is modified to be used locally or on a cluster. 

local=TRUE
# Change this for cluster or local:
if(local){
    #setwd('~/Dropbox/QualityControlImpactsFMRI')
    save.input.data = TRUE
    getOption("mc.cores")
    options(mc.cores=1)
    seeds = seq(1, 200, by = 8)
  } else {
    setwd('~/risk_share/QualityControlImpacts')
    save.input.data = FALSE
    options(mc.cores=1)
    seed = seedID
}
 

a = .libPaths()
#.libPaths(c('/home/benjamin.risk/Rlibs',a))

.libPaths()

getOption("mc.cores")
#set.seed(seed, "L'Ecuyer-CMRG")

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

dat=read.csv('./Data/Master_HeadMotion_wholeGroup_partialCorrelations_ic30_20210803.csv',header=T)

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
ic_class = read_excel('./componentLabels_pca85_ica30.xlsx')

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

t.values.lm=NULL
t.values.naive=NULL

# Make all ADOS in TD = 0
dat2$ADOS.Comparable.Total[dat2$PrimaryDiagnosis=='None'] = 0

# specify learners for super learner. gn is the propensity model and Qbar is the outcome model,
# also used for imputing PANESS: 
my.SL.libs.gn= c("SL.earth","SL.glmnet","SL.gam","SL.glm","SL.ranger","SL.step","SL.step.interaction","SL.xgboost","SL.mean")
my.SL.libs.Qbar= c("SL.earth","SL.glmnet","SL.gam","SL.glm","SL.ranger","SL.ridge","SL.step","SL.step.interaction","SL.svm","SL.xgboost","SL.mean")

# Learners that produce issues:
#   SL.bartMachine: NA
#   SL.glminteraction: produces a rank-deficient model

# Predict PANESS and other variables with minimal missingness:
# Use YearOfScan for measures with multiple versions
# NOTE: handedness is only missing for cases that are also missing paness, gai, and DuPaul scores

# Not including GAI in PANESS model doesn't change R^2 much and can recover an additional 7 cases who are only missing PANESS and GAI
paness.predictors = c('PrimaryDiagnosis','ADHD_Secondary','AgeAtScan','Sex','SES.Family','Race2','handedness','CurrentlyOnStimulants',
                      'DuPaulHome.InattentionRaw','DuPaulHome.HyperactivityRaw',
                      'ADOS.Comparable.Total')

#PANESS and ADOS are both missing for one case
ados.predictors = c('PrimaryDiagnosis', 'ADHD_Secondary','AgeAtScan','Sex','SES.Family','Race2','handedness','CurrentlyOnStimulants',
                    'DuPaulHome.InattentionRaw','DuPaulHome.HyperactivityRaw',
                    'WISC.GAI')


ses.predictors = c('PrimaryDiagnosis','ADHD_Secondary','AgeAtScan','Sex','Race2','handedness','CurrentlyOnStimulants',
                   'DuPaulHome.InattentionRaw','DuPaulHome.HyperactivityRaw','ADOS.Comparable.Total', 
                   'WISC.GAI')

gai.predictors = c('PrimaryDiagnosis','ADHD_Secondary','AgeAtScan','Sex','SES.Family', 'Race2','handedness','CurrentlyOnStimulants',
                   'DuPaulHome.InattentionRaw','DuPaulHome.HyperactivityRaw',
                   'ADOS.Comparable.Total')

ina.predictors = c('PrimaryDiagnosis','ADHD_Secondary','AgeAtScan','Sex','SES.Family', 'Race2','handedness','CurrentlyOnStimulants', 
                   'ADOS.Comparable.Total', 
                   'WISC.GAI',
                   'PANESS.TotalOverflowNotAccountingForAge')

hyp.predictors = c('PrimaryDiagnosis','ADHD_Secondary','AgeAtScan','Sex','SES.Family', 'Race2','handedness','CurrentlyOnStimulants', 
                   'ADOS.Comparable.Total', 
                   'WISC.GAI',
                   'PANESS.TotalOverflowNotAccountingForAge')

temp.data = dat2[,c('ADOS.Comparable.Total',ados.predictors)]
adosPredict = complete.cases(temp.data[,2:ncol(temp.data)])
adosFit = complete.cases(temp.data)
ADOS=temp.data$ADOS.Comparable.Total[adosFit]

sum(adosPredict) - sum(adosFit) # Imputes values for 1 child

temp2.data = dat2[,c('PANESS.TotalOverflowNotAccountingForAge',paness.predictors)]
completeCasesPredict = complete.cases(temp2.data[,2:ncol(temp2.data)])
completeCasesFit = complete.cases(temp2.data)
PANESS=temp2.data$PANESS.TotalOverflowNotAccountingForAge[completeCasesFit]
# WISC.GAI has some missing values where PANESS is missing.
# These subjects will be dropped

sum(completeCasesPredict) - sum(completeCasesFit) # Imputes values for 26 children

temp3.data = dat2[,c('SES.Family', ses.predictors)]
sesPredict = complete.cases(temp3.data[,2:ncol(temp3.data)])
sesFit = complete.cases(temp3.data)
SES=temp3.data$SES.Family[sesFit]

sum(sesPredict) - sum(sesFit) # Imputes values for 8 children

temp4.data = dat2[,c('WISC.GAI', gai.predictors)]
gaiPredict = complete.cases(temp4.data[,2:ncol(temp4.data)])
gaiFit = complete.cases(temp4.data)
GAI=temp4.data$WISC.GAI[gaiFit]

sum(gaiPredict) - sum(gaiFit) # Imputes values for 16 children


temp5.data = dat2[,c('DuPaulHome.InattentionRaw', ina.predictors)]
inaPredict = complete.cases(temp5.data[,2:ncol(temp5.data)])
inaFit = complete.cases(temp5.data)
INA=temp5.data$DuPaulHome.InattentionRaw[inaFit]

sum(inaPredict) - sum(inaFit) # Imputes values for 5 children

temp6.data = dat2[,c('DuPaulHome.HyperactivityRaw', hyp.predictors)]
hypPredict = complete.cases(temp6.data[,2:ncol(temp6.data)])
hypFit = complete.cases(temp6.data)
HYP=temp6.data$DuPaulHome.HyperactivityRaw[hypFit]

sum(hypPredict) - sum(hypFit) # Imputes values for 5 children


ados.xmat.fit = data.frame(model.matrix(ADOS.Comparable.Total~.,data=temp.data[adosFit,])[,-1])
temp.data$ADOS.Comparable.Total=0
ados.xmat.predict = data.frame(model.matrix(ADOS.Comparable.Total~.,data=temp.data[adosPredict,])[,-1])

paness.xmat.fit = data.frame(model.matrix(PANESS.TotalOverflowNotAccountingForAge~.,data=temp2.data[completeCasesFit,])[,-1])
temp2.data$PANESS.TotalOverflowNotAccountingForAge=0
paness.xmat.predict = data.frame(model.matrix(PANESS.TotalOverflowNotAccountingForAge~.,data=temp2.data[completeCasesPredict,])[,-1])

ses.xmat.fit = data.frame(model.matrix(SES.Family~.,data=temp3.data[sesFit,])[,-1])
temp3.data$SES.Family=0
ses.xmat.predict = data.frame(model.matrix(SES.Family~.,data=temp3.data[sesPredict,])[,-1])

gai.xmat.fit = data.frame(model.matrix(WISC.GAI~.,data=temp4.data[gaiFit,])[,-1])
temp4.data$WISC.GAI=0
gai.xmat.predict = data.frame(model.matrix(WISC.GAI~.,data=temp4.data[gaiPredict,])[,-1])

ina.xmat.fit = data.frame(model.matrix(DuPaulHome.InattentionRaw~., data=temp5.data[inaFit,])[,-1])
temp5.data$DuPaulHome.InattentionRaw=0
ina.xmat.predict = data.frame(model.matrix(DuPaulHome.InattentionRaw~.,data=temp5.data[inaPredict,])[,-1])

hyp.xmat.fit = data.frame(model.matrix(DuPaulHome.HyperactivityRaw~., data=temp6.data[hypFit,])[,-1])
temp6.data$DuPaulHome.HyperactivityRaw=0
hyp.xmat.predict = data.frame(model.matrix(DuPaulHome.HyperactivityRaw~.,data=temp6.data[hypPredict,])[,-1])


for(seed in seeds){
  set.seed(seed, "L'Ecuyer-CMRG")
  
  paness.model = mcSuperLearner(Y = PANESS, X = paness.xmat.fit,SL.library = my.SL.libs.Qbar, family=gaussian(),cvControl = list(V = 10),method = drtmle:::tmp_method.CC_LS) # 10-fold CV
  paness.model
  paness.model$times$everything
  
  #This function produces an error within CV.SuperLearner:
  #method = drtmle:::tmp_method.CC_LS
  
  #set.seed(1, "L'Ecuyer-CMRG")
  #tic = proc.time()
  #paness.model = CV.SuperLearner(Y = PANESS, X = paness.xmat.fit, SL.library = my.SL.libs, family=gaussian(), cvControl = list(V = 10),parallel = 'multicore') # 10-fold CV
  #proc.time() - tic
  #summary(paness.model)
  #plot(paness.model,type='bw')
  
  #appear to be equivalent:
  #set.seed(1, "L'Ecuyer-CMRG")
  #paness.model = mcSuperLearner(Y = PANESS, X = paness.xmat.fit,SL.library = my.SL.libs, family=gaussian(),cvControl = list(V = 10),method = method.CC_LS) # 10-fold CV
  #paness.model
  
  # check fit:
  dat2$Predicted.PANESS[completeCasesPredict]=predict(paness.model,newdata = paness.xmat.predict)[[1]]
  plot(dat2$Predicted.PANESS,dat2$PANESS.TotalOverflowNotAccountingForAge)
  cor(dat2$Predicted.PANESS,dat2$PANESS.TotalOverflowNotAccountingForAge,use='complete.obs')^2 
  
  dat2$iPANESS.TotalOverflowNotAccountingForAge=dat2$PANESS.TotalOverflowNotAccountingForAge
  dat2$iPANESS.TotalOverflowNotAccountingForAge[is.na(dat2$PANESS.TotalOverflowNotAccountingForAge)]=dat2$Predicted.PANESS[is.na(dat2$PANESS.TotalOverflowNotAccountingForAge)]
  
  #check it worked:
  cor(dat2$PANESS.TotalOverflowNotAccountingForAge,dat2$iPANESS.TotalOverflowNotAccountingForAge,use='complete.obs') # should be equal to one, it is, cool!
  sum(is.na(dat2$PANESS.TotalOverflowNotAccountingForAge))
  sum(is.na(dat2$iPANESS.TotalOverflowNotAccountingForAge))
  
  # ADOS
  ados.model = mcSuperLearner(Y = ADOS, X = ados.xmat.fit, SL.library = my.SL.libs.Qbar, family=gaussian(), cvControl = list(V = 10),method = drtmle:::tmp_method.CC_LS) # 10-fold CV
  ados.model
  ados.model$times$everything
  
  # check fit:
  dat2$Predicted.ADOS[adosPredict]=predict(ados.model, newdata = ados.xmat.predict)[[1]]
  with(dplyr::filter(dat2, PrimaryDiagnosis=="Autism"), plot(Predicted.ADOS, ADOS.Comparable.Total))
  with(dplyr::filter(dat2, PrimaryDiagnosis=="Autism"), cor(Predicted.ADOS, ADOS.Comparable.Total, use='complete.obs')^2)
  
  dat2$iADOS.Comparable.Total=dat2$ADOS.Comparable.Total
  dat2$iADOS.Comparable.Total[is.na(dat2$ADOS.Comparable.Total)]=dat2$Predicted.ADOS[is.na(dat2$ADOS.Comparable.Total)]
  
  #check it worked:
  cor(dat2$ADOS.Comparable.Total,dat2$iADOS.Comparable.Total,use='complete.obs') # should be equal to one, it is, cool!
  sum(is.na(dat2$ADOS.Comparable.Total))
  sum(is.na(dat2$iADOS.Comparable.Total))
  
  # SES
  ses.model = mcSuperLearner(Y = SES, X = ses.xmat.fit, SL.library = my.SL.libs.Qbar, family=gaussian(), cvControl = list(V = 10),method = drtmle:::tmp_method.CC_LS) # 10-fold CV
  ses.model
  ses.model$times$everything
  
  # check fit:
  dat2$Predicted.SES[sesPredict]=predict(ses.model,newdata = ses.xmat.predict)[[1]]
  plot(dat2$Predicted.SES,dat2$SES.Family)
  cor(dat2$Predicted.SES,dat2$SES.Family,use='complete.obs')^2 
  
  dat2$iSES.Family=dat2$SES.Family
  dat2$iSES.Family[is.na(dat2$SES.Family)]=dat2$Predicted.SES[is.na(dat2$SES.Family)]
  
  #check it worked:
  cor(dat2$SES.Family,dat2$iSES.Family,use='complete.obs') # should be equal to one, it is, cool!
  sum(is.na(dat2$SES.Family))
  sum(is.na(dat2$iSES.Family))
  
  # GAI
  gai.model = mcSuperLearner(Y = GAI, X = gai.xmat.fit, SL.library = my.SL.libs.Qbar, family=gaussian(), cvControl = list(V = 10),method = drtmle:::tmp_method.CC_LS) # 10-fold CV
  gai.model
  gai.model$times$everything
  
  # check fit:
  dat2$Predicted.GAI[gaiPredict]=predict(gai.model,newdata = gai.xmat.predict)[[1]]
  plot(dat2$Predicted.GAI,dat2$WISC.GAI)
  cor(dat2$Predicted.GAI,dat2$WISC.GAI,use='complete.obs')^2 
  
  dat2$iWISC.GAI=dat2$WISC.GAI
  dat2$iWISC.GAI[is.na(dat2$WISC.GAI)]=dat2$Predicted.GAI[is.na(dat2$WISC.GAI)]
  
  #check it worked:
  cor(dat2$WISC.GAI,dat2$iWISC.GAI,use='complete.obs') # should be equal to one, it is, cool!
  sum(is.na(dat2$WISC.GAI))
  sum(is.na(dat2$iWISC.GAI))
  
  # DuPaul Inattention
  ina.model = mcSuperLearner(Y = INA, X = ina.xmat.fit, SL.library = my.SL.libs.Qbar, family=gaussian(), cvControl = list(V = 10),method = drtmle:::tmp_method.CC_LS) # 10-fold CV
  ina.model
  ina.model$times$everything
  
  # check fit:
  dat2$Predicted.Inattention[inaPredict]=predict(ina.model,newdata = ina.xmat.predict)[[1]]
  plot(dat2$Predicted.Inattention, dat2$DuPaulHome.InattentionRaw)
  cor(dat2$Predicted.Inattention, dat2$DuPaulHome.InattentionRaw, use='complete.obs')^2 
  
  dat2$iDuPaulHome.InattentionRaw=dat2$DuPaulHome.InattentionRaw
  dat2$iDuPaulHome.InattentionRaw[is.na(dat2$DuPaulHome.InattentionRaw)]=dat2$Predicted.Inattention[is.na(dat2$DuPaulHome.InattentionRaw)]
  
  #check it worked:
  cor(dat2$DuPaulHome.InattentionRaw, dat2$iDuPaulHome.InattentionRaw, use='complete.obs') # should be equal to one, it is, cool!
  sum(is.na(dat2$DuPaulHome.InattentionRaw))
  sum(is.na(dat2$iDuPaulHome.InattentionRaw))
  
  # DuPaul Hyperactivity
  hyp.model = mcSuperLearner(Y = HYP, X = hyp.xmat.fit, SL.library = my.SL.libs.Qbar, family=gaussian(), cvControl = list(V = 10),method = drtmle:::tmp_method.CC_LS) # 10-fold CV
  hyp.model
  hyp.model$times$everything
  
  # check fit:
  dat2$Predicted.Hyperactivity[hypPredict]=predict(hyp.model, newdata = hyp.xmat.predict)[[1]]
  plot(dat2$Predicted.Hyperactivity, dat2$DuPaulHome.HyperactivityRaw)
  cor(dat2$Predicted.Hyperactivity, dat2$DuPaulHome.HyperactivityRaw, use='complete.obs')^2 
  
  dat2$iDuPaulHome.HyperactivityRaw=dat2$DuPaulHome.HyperactivityRaw
  dat2$iDuPaulHome.HyperactivityRaw[is.na(dat2$DuPaulHome.HyperactivityRaw)]=dat2$Predicted.Hyperactivity[is.na(dat2$DuPaulHome.HyperactivityRaw)]
  
  #check it worked:
  cor(dat2$DuPaulHome.HyperactivityRaw, dat2$iDuPaulHome.HyperactivityRaw, use='complete.obs') # should be equal to one, it is, cool!
  sum(is.na(dat2$DuPaulHome.HyperactivityRaw))
  sum(is.na(dat2$iDuPaulHome.HyperactivityRaw))
  
  #<---------------------------------
  
  # Dataset for propensity model:
  #gn.variables = c('KKI_criteria','HeadCoil','YearOfScan','PrimaryDiagnosis','ADHD_Secondary','AgeAtScan','Sex','handedness','CurrentlyOnStimulants','iPANESS.TotalOverflowNotAccountingForAge','WISC.GAI','DuPaulHome.InattentionRaw','DuPaulHome.HyperactivityRaw','ADOS.Comparable.Total','ADOS.Comparable.StereotypedBehaviorsRestrictedInterests')
  
  gn.variables = c('KKI_criteria','PrimaryDiagnosis','ADHD_Secondary','AgeAtScan','handedness','CurrentlyOnStimulants','iPANESS.TotalOverflowNotAccountingForAge','iWISC.GAI','iDuPaulHome.InattentionRaw','iDuPaulHome.HyperactivityRaw','iADOS.Comparable.Total')
  
  # these indices will be used in the outcome model and drtmle as well:
  temp.data = dat2[,c(gn.variables)]
  
  # note: include variables from the linear model here to define a consistent set of
  # complete predictor cases for the propensity and outcome models
  idx.all.cc = complete.cases(temp.data) & complete.cases(dat2[,c('Sex','iSES.Family','Race2')])
  temp.data = temp.data[idx.all.cc,]
  #temp.data$AgeAtScanXdx = (temp.data$PrimaryDiagnosis=='Autism')*temp.data$AgeAtScan
  #AgeAtScan
  gn.xmat = data.frame(model.matrix(KKI_criteria~.,data=temp.data)[,-1])
  
  Delta.KKI = ifelse(temp.data$KKI_criteria=='Pass',1,0)
  #corrplot::corrplot(cor(gn.xmat),method='number')
  
  # Create a variable equal to one if the observation is used:
  # NOTE: due to missing observations in variables in the initial linear
  # model, need to include linear model variables in this: 
  dat2$CompletePredictorCases = idx.all.cc
  sum(dat2$CompletePredictorCases)
  
  
  # fit with glm and gam: eventually, compare AUCs
  glm.prop.model = glm(Delta.KKI~as.matrix(gn.xmat),family=binomial)
  propensities.glm = predict(glm.prop.model,type = 'response')
  gam.prop.model = mgcv::gam(Delta.KKI~PrimaryDiagnosisNone+ADHD_Secondary+handednessMixed+handednessRight+CurrentlyOnStimulants+s(AgeAtScan)+
                               s(iDuPaulHome.InattentionRaw)+s(iDuPaulHome.HyperactivityRaw)+
                               s(iWISC.GAI)+s(iPANESS.TotalOverflowNotAccountingForAge)+s(iADOS.Comparable.Total),method='REML',family=binomial,data=gn.xmat)
  propensities.gam=predict(gam.prop.model,type='response')
  
  (propensity.KKI = mcSuperLearner(Y = Delta.KKI, X = gn.xmat, family=binomial(link='logit'),SL.library = my.SL.libs.gn, cvControl = list(V = 10), method='method.CC_nloglik')) # 10-fold CV
  #warning: mgcv and gam packages are both in use. You might see an error because both packages use the same function names
  
  # check stability of propensities:
  min(propensity.KKI$SL.predict[Delta.KKI==1])
  
  # examine if there is a positive probability of usability among all data (positivity):
  min(propensity.KKI$SL.predict)
  #hist(propensity.KKI$SL.predict)
  
  propensities.SL = propensity.KKI$SL.predict
  summary(rocit(score=propensities.glm,class=Delta.KKI,method='nonparametric'))
  summary(rocit(score=propensities.gam,class=Delta.KKI,method='nonparametric'))
  summary(rocit(score=propensities.SL,class=Delta.KKI,method='nonparametric'))
  # For most seeds, much better fit with SuperLearner
  
  # merge propensities back to dataset:
  prop.asdtd = data.frame('ID'=dat2$ID[idx.all.cc],'propensities.glm'=propensities.glm,'propensities.gam'=propensities.gam,'propensities.SL'=propensities.SL,Delta.KKI)
  dat3 = merge(dat2,prop.asdtd,all = TRUE) # here, we keep all observations
  
  # this merge changes the ordering of ID. Hence, new index vectors are created,
  # and new xmat need to be made even when using the same variables in the 
  # propensity and outcome models

    # check that the merge is correct:
  table(dat3$KKI_criteria,dat3$Delta.KKI)
  # these should be all in agreement (0 on off diagonal)
  all(idx.all.cc==!is.na(dat3$propensities.SL))
  
  # Delta should be correlated with propensities: 
  cor(1*(dat3$KKI_criteria=='Pass'),dat3$propensities.SL,use='pairwise.complete.obs')
  
  ### Create outcome regression datasets:
  # NOTE: Currently using same variables for propensity and outcome model. 
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
  
  #save(file=paste0('~/Dropbox/QualityControlImpactsFMRI/Data/PropensitiesXmats_seed',seed,'.RData'),Qn.variables,dat3,idx.all.cc,idx.pass.cc,edgeList,idx.all.cc.asd,idx.all.cc.td)
  
  # These datasets necessary to run SuperLearner without error
  # (note: ordering of indices differs from the propensity models due to the merge with the propensities)
  temp.data = dat3[idx.pass.cc, Qn.variables]
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
  
  
  # define startEdgeIdx using the residuals:
  #startEdgeidx=which(names(dat3)=='r.ic1.ic2')
  
  results.df = data.frame('EdgeID'=numeric(nEdges),'EdgeName'=numeric(nEdges),'mean.ASD.naive'=numeric(nEdges),'mean.TD.naive'=numeric(nEdges),'mean.diff.naive'=numeric(nEdges),'z.stat.ASD.naive'=numeric(nEdges),'z.stat.TD.naive'=numeric(nEdges),'z.stat.diff.naive'=numeric(nEdges),'mean.ASD.SL'=numeric(nEdges),'mean.TD.SL'=numeric(nEdges),'mean.diff.SL'=numeric(nEdges),'z.stat.ASD.SL'=numeric(nEdges),'z.stat.TD.SL'=numeric(nEdges),'z.stat.diff.SL'=numeric(nEdges))
  
  startREdgeidx = length(dat3)
  
  ###########################################
  ###########################
  # OUTCOME MODEL and DRTMLE Estimates:
  for (edgeidx in 1:nEdges) {
    dat3.edgeidx = startEdgeidx+edgeidx-1 
    dat3.redgeidx = startREdgeidx+edgeidx
    
    temp = dat3[,c('PrimaryDiagnosis','MeanFramewiseDisplacement.KKI','MaxFramewiseDisplacement.KKI',"FramesWithFDLessThanOrEqualTo250microns","Sex","Race2","iSES.Family","ic1.ic2")]
    lm.completeCases = complete.cases(temp)
    # Motion variables:
    # MeanFramewiseDisplacement.KKI
    # MaxFramewiseDisplacement.KKI
    # FramesWithFDLessThanOrEqualTo250microns
    # Here, we control for sex, imputed SES, and motion effects. We do not control for age, but 
    # age is included in the propensity and outcome models, resulting
    # in the age distribution for each group (including fails), which is 
    # approximately equal between groups. 
    # NOTE: for better interpretation of the plot of mean changes, center all variables:
    
    model.temp = lm(dat3[,dat3.edgeidx]~scale(dat3$MeanFramewiseDisplacement.KKI,center = TRUE, scale=FALSE)+
                      scale(dat3$MaxFramewiseDisplacement.KKI,center=TRUE,scale=FALSE)+
                      scale(dat3$FramesWithFDLessThanOrEqualTo250microns,center=TRUE,scale=FALSE)+
                      dat3$Sex+dat3$Race2+scale(dat3$iSES.Family,center=TRUE,scale=FALSE)+dat3$PrimaryDiagnosis)
    # create the motion-adjusted fconn data:
    dat3[lm.completeCases,paste0('r.',names(dat3)[dat3.edgeidx])]=residuals(model.temp)+coef(model.temp)["(Intercept)"]+
      coef(model.temp)["dat3$PrimaryDiagnosisNone"]*(dat3$PrimaryDiagnosis[lm.completeCases]=='None')
    
    # Audits:
    # check ordering is same:
    trash.audit = residuals(model.temp)+fitted(model.temp)
    trash = abs(dat3[lm.completeCases,dat3.edgeidx]-trash.audit)<1e-16
    if(!all(trash)) stop('Ordering mismatch -- check variables in complete cases')     # check that t-statistic on residuals is approximately equal to t-statistic from lm:
    t.values.lm = c(t.values.lm,coef(summary(model.temp))['dat3$PrimaryDiagnosisNone',3])
    t.values.naive =   c(t.values.naive, t.test(dat3[lm.completeCases & dat3$PrimaryDiagnosis=='None', paste0('r.',names(dat3)[dat3.edgeidx])], 
                                                dat3[lm.completeCases & dat3$PrimaryDiagnosis=='Autism',paste0('r.',names(dat2)[dat3.edgeidx])])$statistic)
    
    
    results.df[edgeidx,'EdgeID'] = edgeidx
    
    results.df[edgeidx,'EdgeName'] = paste0('r.',names(dat3)[dat3.edgeidx])
    
    # naive estimates:
    results.df[edgeidx,'mean.ASD.naive'] = mean(dat3[idx.pass.cc.asd,dat3.redgeidx])
    results.df[edgeidx,'mean.TD.naive'] = mean(dat3[idx.pass.cc.td,dat3.redgeidx])
    results.df[edgeidx,'mean.diff.naive'] =  mean(dat3[idx.pass.cc.asd,dat3.redgeidx])-mean(dat3[idx.pass.cc.td,dat3.redgeidx])
    results.df[edgeidx,'z.stat.ASD.naive'] = t.test(dat3[idx.pass.cc.asd,dat3.redgeidx])$statistic[[1]]
    results.df[edgeidx,'z.stat.TD.naive'] = t.test(dat3[idx.pass.cc.td,dat3.redgeidx])$statistic[[1]]
    results.df[edgeidx,'z.stat.diff.naive'] = t.test(dat3[idx.pass.cc.asd,dat3.redgeidx],dat3[idx.pass.cc.td,dat3.redgeidx])$statistic[[1]]
    
    outcome.SL = mcSuperLearner(Y = dat3[idx.pass.cc,dat3.redgeidx],X=Qn.xmat.fit,family=gaussian(), SL.library = my.SL.libs.Qbar,cvControl = list(V = 10), method = drtmle:::tmp_method.CC_LS)
    
    Qbar.SL.asd = predict(outcome.SL, newdata = Qn.xmat.predict.asd)[[1]]
    Qbar.SL.td = predict(outcome.SL, newdata = Qn.xmat.predict.td)[[1]]
    # rank deficiency in glm.interaction
    
    mean_fconn_asd.SL <- drtmle(Y = dat3[idx.all.cc.asd, dat3.redgeidx],
                                A = dat3[idx.all.cc.asd, c('Delta.KKI')],
                                W = NULL, # Does not do anything with user-input Qn and gn
                                a_0 = 1, # predict for counterfactual that data are usable
                                Qn = list(Qbar.SL.asd), # pass in fitted values
                                gn = list(dat3$propensities.SL[idx.all.cc.asd]), # pass in fitted values
                                SL_Qr = "SL.npreg",
                                SL_gr = "SL.npreg",
                                maxIter = 1 # between 1 and 3 is probably fine
    )
    
    mean_fconn_td.SL <- drtmle(Y = dat3[idx.all.cc.td, dat3.redgeidx],
                               A = dat3[idx.all.cc.td, c('Delta.KKI')], 
                               W = NULL, a_0 = 1, Qn = list(Qbar.SL.td), 
                               gn = list(dat3$propensities.SL[idx.all.cc.td]), 
                               SL_Qr = "SL.npreg", 
                               SL_gr = "SL.npreg", 
                               maxIter = 1)
    
    results.df[edgeidx,'mean.ASD.SL'] = mean_fconn_asd.SL$drtmle$est
    results.df[edgeidx,'mean.TD.SL'] = mean_fconn_td.SL$drtmle$est
    results.df[edgeidx,'mean.diff.SL'] = mean_fconn_asd.SL$drtmle$est-mean_fconn_td.SL$drtmle$est
    
    results.df[edgeidx,'z.stat.ASD.SL'] = mean_fconn_asd.SL$drtmle$est/sqrt(mean_fconn_asd.SL$drtmle$cov)
    results.df[edgeidx,'z.stat.TD.SL'] = mean_fconn_td.SL$drtmle$est/sqrt(mean_fconn_td.SL$drtmle$cov)
    results.df[edgeidx,'z.stat.diff.SL'] = (mean_fconn_asd.SL$drtmle$est - mean_fconn_td.SL$drtmle$est)/sqrt(mean_fconn_asd.SL$drtmle$cov + mean_fconn_td.SL$drtmle$cov)
    message(paste0('Seed ',seed,': Finished Edge',edgeidx))
  }
  
  # save datasets to be loaded for DRTMLE:
  # NOTE: These datasets include the propensities, which change with each seed:
  if (save.input.data) {
    save(file=paste0('./Data/DataWithPropensities_seed',seed,'.RData'),dat3)
  }
  
  results.df$seed = seed
  save(file=paste0('./Results/ic30_pc85_glm_gam_drtmle_seed',seed,'.RData'),results.df)
  
  proc.time()-tic
  message(paste0('Seed ',seed,' Finished'))
}

