

```r
# BRisk
```

This code examines some aspects of the real data for comparison with the simulated data. 
It uses the revised dataset that has no imputation. 


```r
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
```

```
## [1] "/Library/Frameworks/R.framework/Versions/4.1/Resources/library"
```

```r
getOption("mc.cores")
```

```
## [1] 8
```

```r
set.seed(seed, "L'Ecuyer-CMRG")
```

```
## drtmle: TMLE with doubly robust inference
```

```
## Version: 1.1.0
```

```
## Loading required package: nnls
```

```
## Loading required package: gam
```

```
## Loading required package: splines
```

```
## Loading required package: foreach
```

```
## Loaded gam 1.20
```

```
## Super Learner
```

```
## Version: 2.0-28
```

```
## Package created on 2021-05-04
```

```
## Loading required package: Formula
```

```
## Loading required package: plotmo
```

```
## Loading required package: plotrix
```

```
## Loading required package: TeachingDemos
```

```
## Warning: package 'nloptr' was built under R version 4.1.2
```

```
## Warning: package 'tidyr' was built under R version 4.1.2
```

```
## Loading required package: Matrix
```

```
## 
## Attaching package: 'Matrix'
```

```
## The following objects are masked from 'package:tidyr':
## 
##     expand, pack, unpack
```

```
## Loaded glmnet 4.1-3
```

```r
dat=read.csv('../Data/Master_HeadMotion_wholeGroup_partialCorrelations_ic30_20210803.csv',header=T)

# sort by ID: the simplifies the later indexing, as a merge with 
# propensities sorts by ID:
dat = dat[order(dat$ID),]

table(dat$PrimaryDiagnosis)
```

```
## 
## Autism   None 
##    173    372
```

```r
table(dat$Race,dat$Sex,useNA='always') # Used in Protection of Human Subjects table in grant
```

```
##                   
##                      F   M <NA>
##   African American  14  37    0
##   Asian             10  21    0
##   Biracial          20  41    0
##   Caucasian         95 302    0
##   Hispanic           0   1    0
##   Other              0   1    0
##   Unknown            0   2    0
##   <NA>               0   1    0
```

```r
# note: create a "caucasian/other" race variable that allocates the four children without race to caucasian
dat$Race2 = dat$Race
dat$Race2[dat$Race2%in%c('Hispanic','Other','Unknown')]='Caucasian'
table(dat$Race2)
```

```
## 
## African American            Asian         Biracial        Caucasian 
##               51               31               61              401
```

```r
# Add ADHD variable:
dat$ADHD_Secondary = ifelse(dat$ADHD.Subtype=='No dx',0,1)
table(dat$ADHD_Secondary,dat$ADHD.Subtype,useNA='always') # there should be no strange categories. 2/3/2021: Looks good
```

```
##       
##        Combined Hyperactive/Impulsive Inattentive No dx <NA>
##   0           0                     0           0   425    0
##   1          69                     7          44     0    0
##   <NA>        0                     0           0     0    0
```

```r
table(dat$PrimaryDiagnosis[dat$KKI_criteria=='Pass'])
```

```
## 
## Autism   None 
##    114    308
```

```r
table(dat$PrimaryDiagnosis[dat$Ciric_length=='Pass'])
```

```
## 
## Autism   None 
##     29    151
```

```r
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
```

```
##   [1] "ic1.ic2"   "ic1.ic4"   "ic1.ic8"   "ic1.ic13"  "ic1.ic14"  "ic1.ic15" 
##   [7] "ic1.ic17"  "ic1.ic19"  "ic1.ic21"  "ic1.ic22"  "ic1.ic24"  "ic1.ic25" 
##  [13] "ic1.ic26"  "ic1.ic27"  "ic1.ic28"  "ic1.ic29"  "ic1.ic30"  "ic2.ic4"  
##  [19] "ic2.ic8"   "ic2.ic13"  "ic2.ic14"  "ic2.ic15"  "ic2.ic17"  "ic2.ic19" 
##  [25] "ic2.ic21"  "ic2.ic22"  "ic2.ic24"  "ic2.ic25"  "ic2.ic26"  "ic2.ic27" 
##  [31] "ic2.ic28"  "ic2.ic29"  "ic2.ic30"  "ic4.ic8"   "ic4.ic13"  "ic4.ic14" 
##  [37] "ic4.ic15"  "ic4.ic17"  "ic4.ic19"  "ic4.ic21"  "ic4.ic22"  "ic4.ic24" 
##  [43] "ic4.ic25"  "ic4.ic26"  "ic4.ic27"  "ic4.ic28"  "ic4.ic29"  "ic4.ic30" 
##  [49] "ic8.ic13"  "ic8.ic14"  "ic8.ic15"  "ic8.ic17"  "ic8.ic19"  "ic8.ic21" 
##  [55] "ic8.ic22"  "ic8.ic24"  "ic8.ic25"  "ic8.ic26"  "ic8.ic27"  "ic8.ic28" 
##  [61] "ic8.ic29"  "ic8.ic30"  "ic13.ic14" "ic13.ic15" "ic13.ic17" "ic13.ic19"
##  [67] "ic13.ic21" "ic13.ic22" "ic13.ic24" "ic13.ic25" "ic13.ic26" "ic13.ic27"
##  [73] "ic13.ic28" "ic13.ic29" "ic13.ic30" "ic14.ic15" "ic14.ic17" "ic14.ic19"
##  [79] "ic14.ic21" "ic14.ic22" "ic14.ic24" "ic14.ic25" "ic14.ic26" "ic14.ic27"
##  [85] "ic14.ic28" "ic14.ic29" "ic14.ic30" "ic15.ic17" "ic15.ic19" "ic15.ic21"
##  [91] "ic15.ic22" "ic15.ic24" "ic15.ic25" "ic15.ic26" "ic15.ic27" "ic15.ic28"
##  [97] "ic15.ic29" "ic15.ic30" "ic17.ic19" "ic17.ic21" "ic17.ic22" "ic17.ic24"
## [103] "ic17.ic25" "ic17.ic26" "ic17.ic27" "ic17.ic28" "ic17.ic29" "ic17.ic30"
## [109] "ic19.ic21" "ic19.ic22" "ic19.ic24" "ic19.ic25" "ic19.ic26" "ic19.ic27"
## [115] "ic19.ic28" "ic19.ic29" "ic19.ic30" "ic21.ic22" "ic21.ic24" "ic21.ic25"
## [121] "ic21.ic26" "ic21.ic27" "ic21.ic28" "ic21.ic29" "ic21.ic30" "ic22.ic24"
## [127] "ic22.ic25" "ic22.ic26" "ic22.ic27" "ic22.ic28" "ic22.ic29" "ic22.ic30"
## [133] "ic24.ic25" "ic24.ic26" "ic24.ic27" "ic24.ic28" "ic24.ic29" "ic24.ic30"
## [139] "ic25.ic26" "ic25.ic27" "ic25.ic28" "ic25.ic29" "ic25.ic30" "ic26.ic27"
## [145] "ic26.ic28" "ic26.ic29" "ic26.ic30" "ic27.ic28" "ic27.ic29" "ic27.ic30"
## [151] "ic28.ic29" "ic28.ic30" "ic29.ic30"
```

```r
(nEdges = endEdgeidx - startEdgeidx+1)
```

```
## [1] 153
```

```r
# note: check whether all "PASS" have correlations. Should be all FALSE:
table(is.na(dat2[,startEdgeidx]) & dat2$KKI_criteria=='Pass')
```

```
## 
## FALSE 
##   545
```

```r
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
```

```
## [1] 485
```

For RtoR describing simulations section-------------------------------->


```r
Delta.KKI = ifelse(temp.data$KKI_criteria=='Pass',1,0)
mean(Delta.KKI[gn.xmat$PrimaryDiagnosisNone==0])
```

```
## [1] 0.7153285
```

```r
mean(Delta.KKI[gn.xmat$PrimaryDiagnosisNone==1])
```

```
## [1] 0.8390805
```

```r
model.gam = mgcv::gam(Delta.KKI~s(ADOS.Comparable.Total),data=gn.xmat,family=binomial)
summary(model.gam)
```

```
## 
## Family: binomial 
## Link function: logit 
## 
## Formula:
## Delta.KKI ~ s(ADOS.Comparable.Total)
## 
## Parametric coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   1.4533     0.1182   12.29   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                            edf Ref.df Chi.sq p-value   
## s(ADOS.Comparable.Total) 1.513  1.836   15.3 0.00152 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.0302   Deviance explained = 2.86%
## UBRE = -0.028622  Scale est. = 1         n = 485
```

```r
mean(gn.xmat$ADOS.Comparable.Total[gn.xmat$PrimaryDiagnosisNone==0])
```

```
## [1] 14.34307
```

```r
newd <- gn.xmat[1, ] # grab any row; we are going to change temperature only
newd$ADOS.Comparable.Total <- mean(gn.xmat$ADOS.Comparable.Total[gn.xmat$PrimaryDiagnosisNone==0]) - 1e-05 # subtract some small number
y1 <- predict(model.gam, newd)
newd$ADOS.Comparable.Total <- mean(gn.xmat$ADOS.Comparable.Total[gn.xmat$PrimaryDiagnosisNone==0]) + 1e-05 # add some small number
y2 <- predict(model.gam, newd)
```

slope at mean ADOS in real data (full sample)


```r
(y2 - y1)/2e-05
```

```
##         467 
## -0.07740867
```

effect at mean ADOS in real data (full sample)


```r
mean(gn.xmat$ADOS.Comparable.Total[gn.xmat$PrimaryDiagnosisNone==0])*(y2 - y1)/2e-05
```

```
##       467 
## -1.110278
```

```r
newd <- gn.xmat[1, ] # grab any row; we are going to change temperature only
newd$ADOS.Comparable.Total <- 20 - 1e-05 # subtract some small number
y1 <- predict(model.gam, newd)
newd$ADOS.Comparable.Total <- 20 + 1e-05 # add some small number
y2 <- predict(model.gam, newd)
(y2 - y1)/2e-05
```

```
##        467 
## -0.1000773
```

```r
## <-----------------------------------------

(propensity.KKI = mcSuperLearner(Y = Delta.KKI, X = gn.xmat, family=binomial(link='logit'),SL.library = my.SL.libs.gn, cvControl = list(V = 10), method='method.CC_nloglik')) # 10-fold CV
```

```
## 
## Call:  
## mcSuperLearner(Y = Delta.KKI, X = gn.xmat, family = binomial(link = "logit"),  
##     SL.library = my.SL.libs.gn, method = "method.CC_nloglik", cvControl = list(V = 10)) 
## 
## 
## 
##                             Risk       Coef
## SL.earth_All            485.5807 0.13608363
## SL.glmnet_All           464.0483 0.53460726
## SL.gam_All              478.7078 0.00000000
## SL.glm_All              475.4646 0.00000000
## SL.ranger_All           475.5737 0.13471628
## SL.step_All             477.2543 0.00000000
## SL.step.interaction_All 548.0877 0.00000000
## SL.xgboost_All          545.0375 0.08515944
## SL.mean_All             481.6851 0.10943339
```

```r
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
```

```
##       
##          0   1
##   Fail  95   0
##   Pass   0 390
```

```r
# these should be all in agreement (0 on off diagonal)
all(idx.all.cc==!is.na(dat3$propensities.SL))
```

```
## [1] TRUE
```

```r
# Delta should be correlated with propensities: 
cor(1*(dat3$KKI_criteria=='Pass'),dat3$propensities.SL,use='pairwise.complete.obs')
```

```
## [1] 0.6619975
```

Create outcome regression datasets:
NOTE: using same variables for propensity and outcome model. 


```r
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
```

```
## [1] TRUE
```

```r
Qn.xmat.predict.asd = data.frame(model.matrix(numeric(nrow(temp.data))~.,data=temp.data)[,-c(1,2)])

temp.data = Qn.xmat.predict[Qn.xmat.predict$PrimaryDiagnosisNone==1,]
# check that the order matches idx.all.cc.td:
all(temp.data$ID==dat3$ID[idx.all.cc.td])
```

```
## [1] TRUE
```

```r
Qn.xmat.predict.td = data.frame(model.matrix(numeric(nrow(temp.data))~.,data=temp.data))[,-c(1,2)]


## FOR RtoR----->

par(mfrow=c(1,2))
hist(Qn.xmat.predict.asd$ADOS.Comparable.Total,main='a) ASD severity in real data (full sample)', col="#7FAE88")
abline(v=mean(Qn.xmat.predict.asd$ADOS.Comparable.Total),col='red', lwd=3)
```

<img src="tmp//unnamed-chunk-5-1.png" width="672" />

Mean ADOS severity in the real data (full sample)


```r
mean(Qn.xmat.predict.asd$ADOS.Comparable.Total)
```

```
## [1] 14.34307
```

```r
# usable and unusable ADOS:
quantile(Qn.xmat.predict$ADOS.Comparable.Total[Qn.xmat.predict$PrimaryDiagnosisNone==0])
```

```
##   0%  25%  50%  75% 100% 
##    7   11   14   17   26
```

```r
hist(Qn.xmat.fit$ADOS.Comparable.Total[Qn.xmat.fit$PrimaryDiagnosisNone==0], main='B) ASD severity in real data (usable)', col="#5D60AA")
abline(v=mean(Qn.xmat.fit$ADOS.Comparable.Total[Qn.xmat.fit$PrimaryDiagnosisNone==0]),col='red', lwd=3) 
```

<img src="tmp//unnamed-chunk-6-1.png" width="672" />

Mean ADOS severity in the real data (usable after lenient motion QC)


```r
mean(Qn.xmat.fit$ADOS.Comparable.Total[Qn.xmat.fit$PrimaryDiagnosisNone==0])
```

```
## [1] 13.90816
```

```r
# usable ADOS:
quantile(Qn.xmat.fit$ADOS.Comparable.Total[Qn.xmat.fit$PrimaryDiagnosisNone==0])
```

```
##    0%   25%   50%   75%  100% 
##  7.00 11.00 13.00 16.75 23.00
```

Examine the range of correlations for ADOS and the range of intercepts for ASD:


```r
startEdgeidx=which(names(dat3)=='r.ic1.ic2')

p.value.ados = rep(NA,nEdges)
cor.ados=rep(NA,nEdges)
for (edgeidx in 1:nEdges) {
  dat3.edgeidx = startEdgeidx+edgeidx-1 
  cor.ados[edgeidx]=cor(dat3[idx.pass.cc,dat3.edgeidx],dat3$ADOS.Comparable.Total[idx.pass.cc])
  model.gam.out = gam(dat3[idx.pass.cc,dat3.edgeidx]~s(ADOS.Comparable.Total)+PrimaryDiagnosis,data=dat3[idx.pass.cc,])
  plot(model.gam.out,terms = 's(ADOS.Comparable.Total)')
}
```

<img src="tmp//unnamed-chunk-8-1.png" width="672" /><img src="tmp//unnamed-chunk-8-2.png" width="672" /><img src="tmp//unnamed-chunk-8-3.png" width="672" /><img src="tmp//unnamed-chunk-8-4.png" width="672" /><img src="tmp//unnamed-chunk-8-5.png" width="672" /><img src="tmp//unnamed-chunk-8-6.png" width="672" /><img src="tmp//unnamed-chunk-8-7.png" width="672" /><img src="tmp//unnamed-chunk-8-8.png" width="672" /><img src="tmp//unnamed-chunk-8-9.png" width="672" /><img src="tmp//unnamed-chunk-8-10.png" width="672" /><img src="tmp//unnamed-chunk-8-11.png" width="672" /><img src="tmp//unnamed-chunk-8-12.png" width="672" /><img src="tmp//unnamed-chunk-8-13.png" width="672" /><img src="tmp//unnamed-chunk-8-14.png" width="672" /><img src="tmp//unnamed-chunk-8-15.png" width="672" /><img src="tmp//unnamed-chunk-8-16.png" width="672" /><img src="tmp//unnamed-chunk-8-17.png" width="672" /><img src="tmp//unnamed-chunk-8-18.png" width="672" /><img src="tmp//unnamed-chunk-8-19.png" width="672" /><img src="tmp//unnamed-chunk-8-20.png" width="672" /><img src="tmp//unnamed-chunk-8-21.png" width="672" /><img src="tmp//unnamed-chunk-8-22.png" width="672" /><img src="tmp//unnamed-chunk-8-23.png" width="672" /><img src="tmp//unnamed-chunk-8-24.png" width="672" /><img src="tmp//unnamed-chunk-8-25.png" width="672" /><img src="tmp//unnamed-chunk-8-26.png" width="672" /><img src="tmp//unnamed-chunk-8-27.png" width="672" /><img src="tmp//unnamed-chunk-8-28.png" width="672" /><img src="tmp//unnamed-chunk-8-29.png" width="672" /><img src="tmp//unnamed-chunk-8-30.png" width="672" /><img src="tmp//unnamed-chunk-8-31.png" width="672" /><img src="tmp//unnamed-chunk-8-32.png" width="672" /><img src="tmp//unnamed-chunk-8-33.png" width="672" /><img src="tmp//unnamed-chunk-8-34.png" width="672" /><img src="tmp//unnamed-chunk-8-35.png" width="672" /><img src="tmp//unnamed-chunk-8-36.png" width="672" /><img src="tmp//unnamed-chunk-8-37.png" width="672" /><img src="tmp//unnamed-chunk-8-38.png" width="672" /><img src="tmp//unnamed-chunk-8-39.png" width="672" /><img src="tmp//unnamed-chunk-8-40.png" width="672" /><img src="tmp//unnamed-chunk-8-41.png" width="672" /><img src="tmp//unnamed-chunk-8-42.png" width="672" /><img src="tmp//unnamed-chunk-8-43.png" width="672" /><img src="tmp//unnamed-chunk-8-44.png" width="672" /><img src="tmp//unnamed-chunk-8-45.png" width="672" /><img src="tmp//unnamed-chunk-8-46.png" width="672" /><img src="tmp//unnamed-chunk-8-47.png" width="672" /><img src="tmp//unnamed-chunk-8-48.png" width="672" /><img src="tmp//unnamed-chunk-8-49.png" width="672" /><img src="tmp//unnamed-chunk-8-50.png" width="672" /><img src="tmp//unnamed-chunk-8-51.png" width="672" /><img src="tmp//unnamed-chunk-8-52.png" width="672" /><img src="tmp//unnamed-chunk-8-53.png" width="672" /><img src="tmp//unnamed-chunk-8-54.png" width="672" /><img src="tmp//unnamed-chunk-8-55.png" width="672" /><img src="tmp//unnamed-chunk-8-56.png" width="672" /><img src="tmp//unnamed-chunk-8-57.png" width="672" /><img src="tmp//unnamed-chunk-8-58.png" width="672" /><img src="tmp//unnamed-chunk-8-59.png" width="672" /><img src="tmp//unnamed-chunk-8-60.png" width="672" /><img src="tmp//unnamed-chunk-8-61.png" width="672" /><img src="tmp//unnamed-chunk-8-62.png" width="672" /><img src="tmp//unnamed-chunk-8-63.png" width="672" /><img src="tmp//unnamed-chunk-8-64.png" width="672" /><img src="tmp//unnamed-chunk-8-65.png" width="672" /><img src="tmp//unnamed-chunk-8-66.png" width="672" /><img src="tmp//unnamed-chunk-8-67.png" width="672" /><img src="tmp//unnamed-chunk-8-68.png" width="672" /><img src="tmp//unnamed-chunk-8-69.png" width="672" /><img src="tmp//unnamed-chunk-8-70.png" width="672" /><img src="tmp//unnamed-chunk-8-71.png" width="672" /><img src="tmp//unnamed-chunk-8-72.png" width="672" /><img src="tmp//unnamed-chunk-8-73.png" width="672" /><img src="tmp//unnamed-chunk-8-74.png" width="672" /><img src="tmp//unnamed-chunk-8-75.png" width="672" /><img src="tmp//unnamed-chunk-8-76.png" width="672" /><img src="tmp//unnamed-chunk-8-77.png" width="672" /><img src="tmp//unnamed-chunk-8-78.png" width="672" /><img src="tmp//unnamed-chunk-8-79.png" width="672" /><img src="tmp//unnamed-chunk-8-80.png" width="672" /><img src="tmp//unnamed-chunk-8-81.png" width="672" /><img src="tmp//unnamed-chunk-8-82.png" width="672" /><img src="tmp//unnamed-chunk-8-83.png" width="672" /><img src="tmp//unnamed-chunk-8-84.png" width="672" /><img src="tmp//unnamed-chunk-8-85.png" width="672" /><img src="tmp//unnamed-chunk-8-86.png" width="672" /><img src="tmp//unnamed-chunk-8-87.png" width="672" /><img src="tmp//unnamed-chunk-8-88.png" width="672" /><img src="tmp//unnamed-chunk-8-89.png" width="672" /><img src="tmp//unnamed-chunk-8-90.png" width="672" /><img src="tmp//unnamed-chunk-8-91.png" width="672" /><img src="tmp//unnamed-chunk-8-92.png" width="672" /><img src="tmp//unnamed-chunk-8-93.png" width="672" /><img src="tmp//unnamed-chunk-8-94.png" width="672" /><img src="tmp//unnamed-chunk-8-95.png" width="672" /><img src="tmp//unnamed-chunk-8-96.png" width="672" /><img src="tmp//unnamed-chunk-8-97.png" width="672" /><img src="tmp//unnamed-chunk-8-98.png" width="672" /><img src="tmp//unnamed-chunk-8-99.png" width="672" /><img src="tmp//unnamed-chunk-8-100.png" width="672" /><img src="tmp//unnamed-chunk-8-101.png" width="672" /><img src="tmp//unnamed-chunk-8-102.png" width="672" /><img src="tmp//unnamed-chunk-8-103.png" width="672" /><img src="tmp//unnamed-chunk-8-104.png" width="672" /><img src="tmp//unnamed-chunk-8-105.png" width="672" /><img src="tmp//unnamed-chunk-8-106.png" width="672" /><img src="tmp//unnamed-chunk-8-107.png" width="672" /><img src="tmp//unnamed-chunk-8-108.png" width="672" /><img src="tmp//unnamed-chunk-8-109.png" width="672" /><img src="tmp//unnamed-chunk-8-110.png" width="672" /><img src="tmp//unnamed-chunk-8-111.png" width="672" /><img src="tmp//unnamed-chunk-8-112.png" width="672" /><img src="tmp//unnamed-chunk-8-113.png" width="672" /><img src="tmp//unnamed-chunk-8-114.png" width="672" /><img src="tmp//unnamed-chunk-8-115.png" width="672" /><img src="tmp//unnamed-chunk-8-116.png" width="672" /><img src="tmp//unnamed-chunk-8-117.png" width="672" /><img src="tmp//unnamed-chunk-8-118.png" width="672" /><img src="tmp//unnamed-chunk-8-119.png" width="672" /><img src="tmp//unnamed-chunk-8-120.png" width="672" /><img src="tmp//unnamed-chunk-8-121.png" width="672" /><img src="tmp//unnamed-chunk-8-122.png" width="672" /><img src="tmp//unnamed-chunk-8-123.png" width="672" /><img src="tmp//unnamed-chunk-8-124.png" width="672" /><img src="tmp//unnamed-chunk-8-125.png" width="672" /><img src="tmp//unnamed-chunk-8-126.png" width="672" /><img src="tmp//unnamed-chunk-8-127.png" width="672" /><img src="tmp//unnamed-chunk-8-128.png" width="672" /><img src="tmp//unnamed-chunk-8-129.png" width="672" /><img src="tmp//unnamed-chunk-8-130.png" width="672" /><img src="tmp//unnamed-chunk-8-131.png" width="672" /><img src="tmp//unnamed-chunk-8-132.png" width="672" /><img src="tmp//unnamed-chunk-8-133.png" width="672" /><img src="tmp//unnamed-chunk-8-134.png" width="672" /><img src="tmp//unnamed-chunk-8-135.png" width="672" /><img src="tmp//unnamed-chunk-8-136.png" width="672" /><img src="tmp//unnamed-chunk-8-137.png" width="672" /><img src="tmp//unnamed-chunk-8-138.png" width="672" /><img src="tmp//unnamed-chunk-8-139.png" width="672" /><img src="tmp//unnamed-chunk-8-140.png" width="672" /><img src="tmp//unnamed-chunk-8-141.png" width="672" /><img src="tmp//unnamed-chunk-8-142.png" width="672" /><img src="tmp//unnamed-chunk-8-143.png" width="672" /><img src="tmp//unnamed-chunk-8-144.png" width="672" /><img src="tmp//unnamed-chunk-8-145.png" width="672" /><img src="tmp//unnamed-chunk-8-146.png" width="672" /><img src="tmp//unnamed-chunk-8-147.png" width="672" /><img src="tmp//unnamed-chunk-8-148.png" width="672" /><img src="tmp//unnamed-chunk-8-149.png" width="672" /><img src="tmp//unnamed-chunk-8-150.png" width="672" /><img src="tmp//unnamed-chunk-8-151.png" width="672" /><img src="tmp//unnamed-chunk-8-152.png" width="672" /><img src="tmp//unnamed-chunk-8-153.png" width="672" />

```r
max(cor.ados)
```

```
## [1] 0.1698128
```

```r
min(cor.ados)
```

```
## [1] -0.2142168
```

```r
## <---------------------
# note on power:
# the wilcoxon power does not vary with the numbers in each group,
# which is a bummer.
# here is an example of how the group size should matter:
library(pwr)
pwr.t2n.test(n1 = 20, n2=80, d = 0.5, sig.level = 0.05, power = NULL,alternative = c("two.sided"))
```

```
## 
##      t test power calculation 
## 
##              n1 = 20
##              n2 = 80
##               d = 0.5
##       sig.level = 0.05
##           power = 0.5081857
##     alternative = two.sided
```

```r
pwr.t2n.test(n1 = 50, n2=50, d = 0.5, sig.level = 0.05, power = NULL,alternative = c("two.sided"))
```

```
## 
##      t test power calculation 
## 
##              n1 = 50
##              n2 = 50
##               d = 0.5
##       sig.level = 0.05
##           power = 0.6968934
##     alternative = two.sided
```

