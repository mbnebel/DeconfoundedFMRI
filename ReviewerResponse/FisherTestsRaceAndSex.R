#BRisk

# Reviewer 1 

# These results were used to support the decision to not use any imputation

load('../Data/imputeMotorOverflow/DataWithPropensities_seed1.RData')

# for ADOS, subset to ASD:
(p.motor.overflow = wilcox.test(ADOS.Comparable.Total~is.na(PANESS.TotalOverflowNotAccountingForAge),data=dat3[dat3$PrimaryDiagnosis=='Autism',]))

# pull mean differences from here:
summary(lm(ADOS.Comparable.Total~is.na(PANESS.TotalOverflowNotAccountingForAge),data=dat3[dat3$PrimaryDiagnosis=='Autism',]))
# in the non-imputed data, ADOS is 2.9 higher in the participants missing data

# still present in the remaining 11 children:
sum(is.na(dat3$iPANESS.TotalOverflowNotAccountingForAge))
summary(lm(ADOS.Comparable.Total~is.na(iPANESS.TotalOverflowNotAccountingForAge),data=dat3[dat3$PrimaryDiagnosis=='Autism',]))
# ados is on average 4.12 higher in the children still missing motor overflow after imputation
wilcox.test(ADOS.Comparable.Total~is.na(iPANESS.TotalOverflowNotAccountingForAge),data=dat3[dat3$PrimaryDiagnosis=='Autism',])

#############################
# other notes on missingness not in RtoR:
# create a vector of pvalues for six variables: motor overflow, GAI, SES, inattention, hyperactivity, stimulants: 
p.missing.motor = rep(NA,6)
p.missing.motor[1]=p.motor.overflow$p.value

# Look at other variables in Table S.1:
p.missing.motor[2]=wilcox.test(ADOS.Comparable.Total~is.na(WISC.GAI),data=dat3[dat3$PrimaryDiagnosis=='Autism',])$p.value

p.missing.motor[3]=
wilcox.test(ADOS.Comparable.Total~is.na(SES.Family),data=dat3[dat3$PrimaryDiagnosis=='Autism',])$p.value

p.missing.motor[4]=
wilcox.test(ADOS.Comparable.Total~is.na(DuPaulHome.InattentionRaw),data=dat3[dat3$PrimaryDiagnosis=='Autism',])$p.value

p.missing.motor[5]=wilcox.test(ADOS.Comparable.Total~is.na(DuPaulHome.HyperactivityRaw),data=dat3[dat3$PrimaryDiagnosis=='Autism',])$p.value

p.missing.motor[6]=wilcox.test(ADOS.Comparable.Total~is.na(CurrentlyOnStimulants),data=dat3[dat3$PrimaryDiagnosis=='Autism',])$p.value
p.adjust(p.missing.motor,method='BH')


# no missing handedness with ADOS
wilcox.test(ADOS.Comparable.Total~is.na(handedness),data=dat3[dat3$PrimaryDiagnosis=='Autism',])

# race only missing one observation:
wilcox.test(ADOS.Comparable.Total~is.na(Race2),data=dat3[dat3$PrimaryDiagnosis=='Autism',])



############################
#########################
# GAI:

# create a vector of pvalues for six variables: motor overflow, GAI, SES, inattention, hyperactivity, stimulants: 
p.missing.GAI = rep(NA,5)
p.missing.GAI[1]=wilcox.test(WISC.GAI~is.na(PANESS.TotalOverflowNotAccountingForAge),data=dat3)$p.value

# ADOS is not missing

p.missing.GAI[2]=  wilcox.test(WISC.GAI~is.na(SES.Family),data=dat3)$p.value
# GAI differs by SES.Family
summary(lm(WISC.GAI~is.na(SES.Family),data=dat3))


p.missing.GAI[3]= wilcox.test(WISC.GAI~is.na(DuPaulHome.InattentionRaw),data=dat3)$p.value

p.missing.GAI[4]=wilcox.test(WISC.GAI~is.na(DuPaulHome.HyperactivityRaw),data=dat3)$p.value

p.missing.GAI[5]=wilcox.test(WISC.GAI~is.na(CurrentlyOnStimulants),data=dat3[dat3$PrimaryDiagnosis=='Autism',])$p.value

p.adjust(p.missing.GAI,method='BH')


# no missing handedness with ADOS
wilcox.test(ADOS.Comparable.Total~is.na(handedness),data=dat3[dat3$PrimaryDiagnosis=='Autism',])

# race only missing one observation:
wilcox.test(ADOS.Comparable.Total~is.na(Race2),data=dat3[dat3$PrimaryDiagnosis=='Autism',])



# Reviewer 3 Q Look at proportions excluded by sex and race
# For this, load the non-imputed data:
library(mgcv)

load('../Data/noImputation/DataWithPropensities_seed1.RData')

names(dat3)
with(dat3,table(Delta.KKI,Sex))

# Lenient:
# Sex:
fisher.test(with(dat3,table(Delta.KKI,Sex)))

# race:
fisher.test(with(dat3,table(Delta.KKI,Race)))

# SES:
summary(gam(Delta.KKI~s(SES.Family),data=dat3,family=binomial))


# Strict:
# sex
fisher.test(with(dat3,table(Ciric_length,Sex)))

# race
fisher.test(with(dat3,table(Ciric_length,Race)))

# SES:
dat3$Delta.Ciric = dat3$Ciric_length=="Pass"
summary(gam(Delta.Ciric~s(SES.Family),data=dat3,family=binomial))


####################
###########
# examine whether motion differs by sex:
library(coin)
library(rstatix)
boxplot(dat3$MeanFramewiseDisplacement.KKI~dat3$Sex,dat3[dat3$CompletePredictorCases==1,])


# complete predictor cases:
dat3$Sex = as.factor(dat3$Sex)
wilcox_test(MeanFramewiseDisplacement.KKI~Sex,dat3[dat3$CompletePredictorCases==1,])
wilcox_effsize(formula=MeanFramewiseDisplacement.KKI~Sex,data=dat3[dat3$CompletePredictorCases==1,])

# complete cases lenient pass:
wilcox.test(MeanFramewiseDisplacement.KKI~Sex,dat3[dat3$CompletePredictorCases==1 & dat3$Delta.KKI==1,])

# complete cases strict pass:
# complete cases lenient pass:
wilcox.test(MeanFramewiseDisplacement.KKI~Sex,dat3[dat3$CompletePredictorCases==1 & dat3$Delta.Ciric==1,])





