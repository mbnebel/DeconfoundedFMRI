---
title: "probability of exclusion as a function of the covariates controlling for sex, race, and SES"
author: "MB Nebel"
date: '2022-05-12'
output: html_document
---




```r
library(grid)
library(xtable)
library(dplyr)
library(ggplot2)
library(gridExtra)
#nest/unnest
library(tidyr)
#map function (kind of like a for loop)
library(purrr)
#tidy model summary
library(broom)

#Also uses functions from plyr, mgcv, cowplot
```

#### Load initial data set

```r
load('./Data/noImputation/DataWithPropensities_seed1.RData')

#convert PrimaryDiagnosis to factor
dat3$PrimaryDiagnosis <- factor(dat3$PrimaryDiagnosis, levels = c("Autism", "None"))
levels(dat3$PrimaryDiagnosis)
```

```
## [1] "Autism" "None"
```

```r
dat3$PrimaryDiagnosis <- relevel(dat3$PrimaryDiagnosis, "None")
levels(dat3$PrimaryDiagnosis) = c("TD", "ASD")
```

#### Reshape data to combine motion quality control (QC) levels

```r
# create dummy factor to include all subjects
dat3$noExclusion <- ifelse(dat3$ID > 0, "Pass", "Pass")
dat3$noExclusion <- factor(dat3$noExclusion, levels = c("Pass", "Fail"))

# convert KKI_criteria to factor with reference level 'Pass'
dat3$KKI_criteria <- factor(dat3$KKI_criteria, levels = c("Pass", "Fail"))

# convert Ciric_length to factor with reference level 'Pass' to match
# KKI_criteria
dat3$Ciric_length <- factor(dat3$Ciric_length, levels = c("Pass", "Fail"))

# combine Ciric_length, KKI, and noExclusion exclusion into one variable
allVariables = c("ID", "PrimaryDiagnosis", "AgeAtScan", "Ciric_length", "KKI_criteria",
    "noExclusion", "PANESS.TotalOverflowNotAccountingForAge", "SRS.Score", "WISC.GAI",
    "DuPaulHome.InattentionRaw", "DuPaulHome.HyperactivityRaw", "ADOS.Comparable.Total",
    "CurrentlyOnStimulants", "HeadCoil", "Sex", "ADHD_Secondary", "SES.Family", "Race2",
    "handedness", "CompletePredictorCases")

idVariables = c("ID", "PrimaryDiagnosis", "AgeAtScan", "PANESS.TotalOverflowNotAccountingForAge",
    "SRS.Score", "WISC.GAI", "DuPaulHome.InattentionRaw", "DuPaulHome.HyperactivityRaw",
    "ADOS.Comparable.Total", "CurrentlyOnStimulants", "HeadCoil", "Sex", "ADHD_Secondary",
    "SES.Family", "Race2", "handedness", "CompletePredictorCases")

qcMelt <- reshape2::melt(dat3[, allVariables], id.vars = names(dat3)[which(names(dat3) %in%
    idVariables)], variable.name = "Motion.Exclusion.Level", value.name = "Included")

# rename exclusion levels NOTE: need None to be highest level for
# geom_split_violin
levels(qcMelt$Motion.Exclusion.Level) <- c("Strict", "Lenient", "None")

# convert Included to factor with pass as reference
qcMelt$Included <- factor(qcMelt$Included, levels = c("Pass", "Fail"))

# rename levels of value
levels(qcMelt$Included) <- c("Included", "Excluded")
```

#### Limit initial dataset to complete cases

```r
dat3 <- filter(dat3, CompletePredictorCases==1)

completeCases <- filter(qcMelt, CompletePredictorCases==1)

completeCases$Sex <- relevel(as.factor(completeCases$Sex), "M")
```



```r
phenoVariables <- c("ID", "PrimaryDiagnosis", 
                    "ADOS.Comparable.Total", 
                    "SRS.Score", 
                    "PANESS.TotalOverflowNotAccountingForAge", 
                    "DuPaulHome.InattentionRaw",
                    "DuPaulHome.HyperactivityRaw", 
                    "AgeAtScan", 
                    "WISC.GAI", "SES.Family",
                    "Motion.Exclusion.Level", "Included", "Sex", "Race2")

phenoIDs <- c("ID", "PrimaryDiagnosis", "Motion.Exclusion.Level", "Included", "Sex", "SES.Family", "Race2")

aim1 <- reshape2::melt(completeCases[, phenoVariables],
                         id.vars=names(completeCases)[which(names(completeCases) %in% phenoIDs)])

levels(aim1$variable) <- c("ADOS", "SRS", "Motor Overflow", "Inattention", 
                           "Hyperactivity", "Age", "GAI")

aim1G <- group_by(aim1, PrimaryDiagnosis, Motion.Exclusion.Level, Included, variable)
```



```r
aim1$delta = rep(NA,length=nrow(aim1))
aim1$delta = ifelse(aim1$Included=="Included",1,0)
aim1tib <- tibble(filter(aim1, Motion.Exclusion.Level!="None"))
aim1tib$Motion.Exclusion.Level <- droplevels(aim1tib$Motion.Exclusion.Level)

#initialize dummy variable for Sex and centered SES.Family
aim1tib$dummySex = rep(NA, length=nrow(aim1tib))
aim1tib$SES.Family_centered = rep(NA, length=nrow(aim1tib))

#if male, dummySex = 1; if female, dummySex = -1
aim1tib$dummySex[aim1tib$Sex=="M"]=1
aim1tib$dummySex[aim1tib$Sex=="F"]=-1

#make Race2 a factor and use sum to to zero contrasts
aim1tib$Race2 = as.factor(aim1tib$Race2)
#you put it in the code before the call to mgcv::gam that includes Race2 as a covariate
contrasts(aim1tib$Race2) = contr.sum(4, contrasts=TRUE)

aim1tib$SES.Family_centered = scale(aim1tib$SES.Family, center=TRUE, scale=FALSE)

#model sex separately, needs to be a factor
aim1Nest <- aim1tib %>% 
  group_by(Motion.Exclusion.Level, variable) %>% 
  tidyr::nest()

#nested models
nested_gams <- aim1Nest %>%
  mutate(model = map(data, ~mgcv::gam(1-delta~s(value, k=-1)+dummySex+Race2+s(SES.Family_centered, k=-1), data = na.omit(.x), 
                                      family=binomial(link=logit), method="REML")), 
         coefs = map(model, tidy, conf.int = FALSE),
         Rsq = map_dbl(model, ~summary(.)$r.sq)) %>% 
  unnest(coefs)

#Ben: correct for 7 lenient and 7 strict 
nested_gams_len <-  nested_gams %>% 
  filter(Motion.Exclusion.Level=="Lenient" & term =="s(value)")

nested_gams_len$p.fdr = p.adjust(nested_gams_len$p.value, method = "BH")

nested_gams_strict <-  nested_gams %>% 
  filter(Motion.Exclusion.Level=="Strict" & term =="s(value)")
  
nested_gams_strict$p.fdr = p.adjust(nested_gams_strict$p.value, method = "BH")

#combine to print
nested_gams <- rbind(nested_gams_len, nested_gams_strict)

#list adjusted p values
nested_gams[, c(1:2,6:11)]
```

```
## # A tibble: 14 × 8
## # Groups:   Motion.Exclusion.Level, variable [14]
##    Motion.Exclusion.Lev… variable   edf ref.df statistic p.value     Rsq   p.fdr
##    <fct>                 <fct>    <dbl>  <dbl>     <dbl>   <dbl>   <dbl>   <dbl>
##  1 Lenient               ADOS      1.52   1.84     15.1  1.71e-3 0.0285  2.99e-3
##  2 Lenient               SRS       1.73   2.16     14.6  9.52e-4 0.0400  2.22e-3
##  3 Lenient               Motor O…  1.00   1.00     15.2  9.80e-5 0.0307  6.86e-4
##  4 Lenient               Inatten…  1.46   1.78      7.06 1.68e-2 0.0136  1.68e-2
##  5 Lenient               Hyperac…  1.91   2.39     12.8  3.49e-3 0.0232  4.08e-3
##  6 Lenient               Age       1.00   1.00      9.43 2.14e-3 0.0137  2.99e-3
##  7 Lenient               GAI       1.00   1.00     11.0  9.34e-4 0.0200  2.22e-3
##  8 Strict                ADOS      1.00   1.00     20.3  7.49e-6 0.0366  3.54e-5
##  9 Strict                SRS       1.00   1.00     19.6  1.01e-5 0.0476  3.54e-5
## 10 Strict                Motor O…  1.00   1.00     10.1  1.49e-3 0.0129  2.09e-3
## 11 Strict                Inatten…  1.00   1.00     15.7  7.56e-5 0.0274  1.32e-4
## 12 Strict                Hyperac…  1.62   2.01     22.2  1.64e-5 0.0427  3.83e-5
## 13 Strict                Age       1.85   2.31      8.26 2.54e-2 0.00851 2.54e-2
## 14 Strict                GAI       1.00   1.00      6.52 1.07e-2 0.00515 1.24e-2
```

```r
#max p value for 7 lenient models
max(nested_gams_len$p.fdr)
```

```
## [1] 0.01675631
```

```r
#max p value for 7 strict models
max(nested_gams_strict$p.fdr)
```

```
## [1] 0.02535675
```

```r
#average SES and sex both 0
#0 didn't work for Race2; used most common = Caucasian
nested_gams <- nested_gams %>% 
           mutate(LB = map(data, ~round(min(na.omit(.x$value)))),
                  UB = map(data, ~round(max(na.omit(.x$value)))),
                  range = map2(LB, UB, ~seq(from=.x, to=.y, by=.2)),
                  #sex = map(value, ~rep(0, length=length(value))),
                  #ses = map(value, ~rep(0, length=length(value))),
                  fit.data = map(range, ~data.frame(value = .x, dummySex = rep(0, length=length(.x)),
                                                    Race2 = rep("Caucasian", length=length(.x)), 
                                                    SES.Family_centered = rep(0, length=length(.x)))),
                  logpredict = map2(model, fit.data, ~predict(.x, newdata = .y, type="link", se.fit=TRUE)),
                  fit = map(logpredict, ~plogis(.x$fit)),
                  lCI = map(logpredict, ~plogis(.x$fit-1.96*.x$se.fit)),
                  hCI = map(logpredict, ~plogis(.x$fit+1.96*.x$se.fit)))
```

#### Define theme for Alternate Figure 4 top row

```r
gam_theme = theme(
  axis.title.x=element_text(size=12),
  axis.title.y=element_text(size=12),
  axis.text.x=element_text(size=8),
  axis.text.y=element_text(size=10),
  plot.title = element_text(size = 16),
  plot.caption = element_text(size = 16,hjust = 0),
  legend.title = element_blank(), legend.position ="none")
```


#### Alternate Figure 4a top. Probability of exclusion as a function of ADOS

```r
ados <- nested_gams %>% 
  filter(variable=="ADOS") %>% 
  select("variable", "Motion.Exclusion.Level", "range", "fit", 'lCI', 'hCI') %>% 
  unnest(c(range, fit, lCI, hCI))

p_ados <- ggplot(ados, aes(x=range, y=fit))+
  geom_line(aes(colour = Motion.Exclusion.Level),size=1.2) + ylim(0,1) + theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI, fill = Motion.Exclusion.Level), linetype='blank', alpha=0.2)+  
  scale_color_manual(labels=c('Strict', 'Lenient'), values = c("#f55154","#9FB0CC"))+
  scale_fill_manual(labels=c('Strict', 'Lenient'), values = c("#f55154","#9FB0CC"))+
  labs(x='', y='Probability of Exclusion', fill='Motion Control', colour='Motion Control')+  
  scale_x_continuous(expand = c(0, 0))+ 
  gam_theme+
  ggtitle("ADOS")+
  theme(plot.title = element_text(size = 11, hjust = 0.5))
```

#### Alternate Figure 4b top. Probability of exclusion as a function of SRS

```r
srs <- nested_gams %>% 
  filter(variable=="SRS") %>% 
  select("variable", "Motion.Exclusion.Level", "range", "fit", 'lCI', 'hCI') %>% 
  unnest(c(range, fit, lCI, hCI))

p_srs <- ggplot(srs, aes(x=range, y=fit))+
  geom_line(aes(colour = Motion.Exclusion.Level),size=1.2)+ylim(0,1)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI, fill=Motion.Exclusion.Level), linetype='blank', alpha=0.2)+  
  scale_color_manual(labels=c('Strict', 'Lenient'), values = c("#f55154","#9FB0CC"))+
  scale_fill_manual(labels=c('Strict', 'Lenient'), values = c("#f55154","#9FB0CC"))+
  labs(x='', y='', fill='Motion Control', colour='Motion Control')+  
  scale_x_continuous(expand = c(0, 0))+
  gam_theme+
  ggtitle("SRS")+
  theme(plot.title = element_text(size = 11, hjust = 0.5))+
  theme(axis.title.y = element_blank())
```

#### Alternate Figure 4c top. Probability of exclusion as a function of Inattention

```r
ina <- nested_gams %>% 
  filter(variable=="Inattention") %>% 
  select("variable", "Motion.Exclusion.Level", "range", "fit", 'lCI', 'hCI') %>% 
  unnest(c(range, fit, lCI, hCI))

p_in <- ggplot(ina, aes(x=range, y=fit))+
  geom_line(aes(colour = Motion.Exclusion.Level),size=1.2)+ylim(0,1)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI, fill=Motion.Exclusion.Level), linetype='blank', alpha=0.2)+  
  scale_color_manual(labels=c('Strict', 'Lenient'), values = c("#f55154","#9FB0CC"))+
  scale_fill_manual(labels=c('Strict', 'Lenient'), values = c("#f55154","#9FB0CC"))+
  labs(x='', y='', fill='Motion Control', colour='Motion Control')+  
  scale_x_continuous(expand = c(0, 0))+ 
  gam_theme+
  ggtitle("Inattention")+
  theme(plot.title = element_text(size = 11, hjust = 0.5))+
  theme(axis.title.y = element_blank())
```

#### Alternate Figure 4d top. Probability of exclusion as a function of Hyperactivity

```r
hi <- nested_gams %>% 
  filter(variable=="Hyperactivity") %>% 
  select("variable", "Motion.Exclusion.Level", "range", "fit", 'lCI', 'hCI') %>% 
  unnest(c(range, fit, lCI, hCI))

p_hi <- ggplot(hi, aes(x=range, y=fit))+
  geom_line(aes(colour = Motion.Exclusion.Level),size=1.2)+ylim(0,1)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI, fill=Motion.Exclusion.Level), linetype='blank', alpha=0.2)+  
  scale_color_manual(labels=c('Strict', 'Lenient'), values = c("#f55154","#9FB0CC"))+
  scale_fill_manual(labels=c('Strict', 'Lenient'), values = c("#f55154","#9FB0CC"))+
  labs(x='', y='', fill='Motion Control', colour='Motion Control')+  
  scale_x_continuous(expand = c(0, 0))+ 
  gam_theme+
  ggtitle("Hyperactivity")+
theme(plot.title = element_text(size = 11, hjust = 0.5))+
  theme(axis.title.y = element_blank())
```

#### Alternate Figure 4e top. Probability of exclusion as a function of Motor Overflow

```r
mo <- nested_gams %>% 
  filter(variable=="Motor Overflow") %>% 
  select("variable", "Motion.Exclusion.Level", "range", "fit", 'lCI', 'hCI') %>% 
  unnest(c(range, fit, lCI, hCI))

p_mo <- ggplot(mo, aes(x=range, y=fit))+
  geom_line(aes(colour = Motion.Exclusion.Level),size=1.2)+ylim(0,1)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI, fill=Motion.Exclusion.Level), linetype='blank', alpha=0.2)+  
  scale_color_manual(labels=c('Strict', 'Lenient'), values = c("#f55154","#9FB0CC"))+
  scale_fill_manual(labels=c('Strict', 'Lenient'), values = c("#f55154","#9FB0CC"))+
  labs(x='', y='', fill='Motion Control', colour='Motion Control')+  
  scale_x_continuous(expand = c(0, 0))+ 
  gam_theme+
  ggtitle("Motor Overflow")+
  theme(plot.title = element_text(size = 11, hjust = 0.5))+
  theme(axis.title.y = element_blank())
```

#### Alternate Figure 4f top. Probability of exclusion as a function of Age

```r
age<- nested_gams %>% 
  filter(variable=="Age") %>% 
  select("variable", "Motion.Exclusion.Level", "range", "fit", 'lCI', 'hCI') %>% 
  unnest(c(range, fit, lCI, hCI))

p_age <- ggplot(age, aes(x=range, y=fit))+
  geom_line(aes(colour = Motion.Exclusion.Level),size=1.2)+ylim(0,1)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI, fill=Motion.Exclusion.Level), linetype='blank', alpha=0.2)+  
  scale_color_manual(labels=c('Strict', 'Lenient'), values = c("#f55154","#9FB0CC"))+
  scale_fill_manual(labels=c('Strict', 'Lenient'), values = c("#f55154","#9FB0CC"))+
  labs(x='', y='', fill='Motion Control', colour='Motion Control')+  
  scale_x_continuous(expand = c(0, 0))+ 
  gam_theme+
  ggtitle("Age")+
  theme(plot.title = element_text(size = 11, hjust = 0.5))+
  theme(axis.title.y = element_blank())
```

#### Alternate Figure 4g top. Probability of exclusion as a function of GAI

```r
gai <- nested_gams %>% 
  filter(variable=="GAI") %>% 
  select("variable", "Motion.Exclusion.Level", "range", "fit", 'lCI', 'hCI') %>% 
  unnest(c(range, fit, lCI, hCI))

p_gai <- ggplot(gai, aes(x=range, y=fit))+
  geom_line(aes(colour = Motion.Exclusion.Level),size=1.2)+ylim(0,1)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI, fill=Motion.Exclusion.Level), linetype='blank', alpha=0.2)+  
  scale_color_manual(labels=c('Strict', 'Lenient'), values = c("#f55154","#9FB0CC"))+
  scale_fill_manual(labels=c('Strict', 'Lenient'), values = c("#f55154","#9FB0CC"))+
  labs(x='', y='', fill='Motion Control', colour='Motion Control')+  
  scale_x_continuous(expand = c(0, 0))+ 
  gam_theme+
  ggtitle("GAI")+
  theme(plot.title = element_text(size = 11, hjust = 0.5))+
  theme(axis.title.y = element_blank())

p_legend = cowplot::get_legend(p_gai + guides(color = guide_legend(nrow = 1))+
                                 theme(legend.position = "bottom", legend.text = element_text(size = 11),
                                       legend.key.size=unit(.15, "in")))
```

#### Define theme for density plots of covariates across included and excluded children



```r
ddata <- nested_gams %>% 
  filter(variable=="ADOS") %>% 
  filter(Motion.Exclusion.Level=="Lenient") %>% 
  select("variable", "Motion.Exclusion.Level", "data", "LB", "UB") %>% 
  unnest(data, LB, UB) %>% 
  filter(PrimaryDiagnosis=="ASD")
```

```
## Warning: unnest() has a new interface. See ?unnest for details.
## Try `df %>% unnest(c(data, LB, UB))`, with `mutate()` if needed
```

```r
ddata$PrimaryDiagnosis <- droplevels(ddata$PrimaryDiagnosis)

d_ados=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB[1], ddata$UB[1]), breaks=seq(0, ddata$UB[1], by=5))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0, .09), breaks=seq(0, .08, by=.02))+  
  labs(x='', y='Density')+
  scale_fill_manual(values = c("#FDE599"))+ 
  scale_color_manual(values = c("#E9D38D"))+ 
  den_theme
```

#### Alternate Figure 4b bottom. SRS density

```r
ddata <- nested_gams %>% 
  filter(variable=="SRS") %>% 
  filter(Motion.Exclusion.Level=="Lenient") %>% 
  select("variable", "Motion.Exclusion.Level", "data") %>% 
  unnest(data)

d_srs=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits=c(0,max(srs$range)),breaks = seq(0, 100 , by = 50))+  
  scale_y_continuous(expand = c(0, 0))+  
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```

#### Alternate Figure 4c bottom. Inattention density

```r
ddata <- nested_gams %>% 
  filter(variable=="Inattention") %>% 
  filter(Motion.Exclusion.Level=="Lenient") %>% 
  select("variable", "Motion.Exclusion.Level", "data") %>% 
  unnest(data)

d_in=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE)+  
  scale_x_continuous(expand = c(0, 0))+  
  scale_y_continuous(expand = c(0, 0))+  
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  labs(x='')+  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  theme(axis.title.y = element_blank())+
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```

#### Alternate Figure 4d bottom. Hyperactivity/Impulsivity Density

```r
ddata <- nested_gams %>% 
  filter(variable=="Hyperactivity") %>% 
  filter(Motion.Exclusion.Level=="Lenient") %>% 
  select("variable", "Motion.Exclusion.Level", "data") %>% 
  unnest(data)

d_hi=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE)+  
  scale_x_continuous(expand = c(0, 0))+  
  scale_y_continuous(expand = c(0, 0))+  
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```

#### Alternate Figure 4e bottom. Motor Overflow Density

```r
ddata <- nested_gams %>% 
  filter(variable=="Motor Overflow") %>% 
  filter(Motion.Exclusion.Level=="Lenient") %>% 
  select("variable", "Motion.Exclusion.Level", "data") %>% 
  unnest(data)

d_mo=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE)+  
  scale_x_continuous(expand = c(0, 0))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .07), breaks=seq(0, .06, by=.03))+ 
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```

#### Alternate Figure 4f bottom. Age Density

```r
ddata <- nested_gams %>% 
  filter(variable=="Age") %>% 
  filter(Motion.Exclusion.Level=="Lenient") %>% 
  select("variable", "Motion.Exclusion.Level", "data") %>% 
  unnest(data)

d_age=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits=c(8,13), breaks = seq(8, 13 , by = 1))+  
  scale_y_continuous(expand = c(0, 0), limits=c(0,.29), breaks=seq(0, .25, by = .05))+ 
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```

#### Alternate Figure 4g bottom. GAI Density

```r
ddata <- nested_gams %>% 
  filter(variable=="GAI") %>% 
  filter(Motion.Exclusion.Level=="Lenient") %>% 
  select("variable", "Motion.Exclusion.Level", "data") %>% 
  unnest(data)

d_gai=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE)+  
  scale_x_continuous(expand = c(0, 0))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0, .035), breaks=seq(0., .03, by=.01))+  
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  labs(x='')+  
  theme(axis.title.y = element_blank())

hist_legend = cowplot::get_legend(d_gai + guides(color = guide_legend(nrow = 1))+theme(legend.position = "bottom", legend.text = element_text(size = 11),   legend.key.size=unit(.15, "in")))
```


#### combine gam plots with densities & print

```r
top_row <- cowplot::plot_grid(p_ados, p_srs, p_in, p_hi, p_mo, p_age, p_gai, ncol=7, 
                              rel_widths=c(1.18/7, .97/7, .97/7, .97/7, .97/7, .97/7, .97/7))
second_row <- cowplot::plot_grid(d_ados, d_srs, d_in, d_hi, d_mo, d_age, d_gai, ncol=7, 
                                 rel_widths=c(1.18/7, .97/7, .97/7, .97/7, .97/7, .97/7, .97/7))
```

```
## Warning: Removed 89 rows containing non-finite values (stat_density).
```

```r
png("./ReviewerResponse/probEx_withSexSESRace/fig_probExclusion_allGAM_TD_ASD_cc_sexRaceSES.png",width=10,height=6,units="in",res=200)
cowplot::plot_grid(p_legend, top_row, NULL, 
                   second_row, NULL, hist_legend, 
                   nrow=6, rel_heights=c(.1, 1, -.01, .5, -.05, .1))
dev.off()
```

```
## quartz_off_screen 
##                 2
```

#### Print alternate figure for the report

```r
cowplot::plot_grid(p_legend, top_row, NULL, 
                   second_row, NULL, hist_legend, 
                   nrow=6, rel_heights=c(.1, 1, -.01, .5, -.05, .1))
```

<img src="forReviewers_probEx_withSexSESRace//unnamed-chunk-24-1.png" width="960" />

