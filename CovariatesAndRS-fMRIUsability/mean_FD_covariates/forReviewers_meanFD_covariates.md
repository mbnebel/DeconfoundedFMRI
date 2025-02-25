---
title: "mean FD as a function of the covariates"
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
#Also uses functions from plyr, scales, mgcv, cowplot
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
    "handedness", "CompletePredictorCases", "YearOfScan", "MeanFramewiseDisplacement.KKI",
    "MeanFramewiseDisplacement")
idVariables = c("ID", "PrimaryDiagnosis", "AgeAtScan", "PANESS.TotalOverflowNotAccountingForAge",
    "SRS.Score", "WISC.GAI", "DuPaulHome.InattentionRaw", "DuPaulHome.HyperactivityRaw",
    "ADOS.Comparable.Total", "CurrentlyOnStimulants", "HeadCoil", "Sex", "ADHD_Secondary",
    "SES.Family", "Race2", "handedness", "CompletePredictorCases", "YearOfScan",
    "MeanFramewiseDisplacement.KKI", "MeanFramewiseDisplacement")
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
Motion QC levels:

1. **Strict motion QC** = Ciric_length

In the strict case, scans were excluded if mean FD exceeded .2 mm or they included less than five minutes of data free from frames with FD exceeding .25 mm

2. **Lenient motion QC** = KKI_criteria

In the lenient case, scans were excluded if the participant had less than 5 minutes of continuous data after removing frames in which the participant moved more than the nominal size of a voxel between any two frames (3 mm) or their head rotated 3\degree. This procedure was modeled after common head motion exclusion criteria for task fMRI data, which rely on voxel size to determine thresholds for unacceptable motion.

3. **None** = all participants

#### Limit initial dataset to complete cases


```r
dat3 <- filter(dat3, CompletePredictorCases==1)

completeCases <- filter(qcMelt, CompletePredictorCases==1)

#make M reference level for sex
completeCases$Sex <- relevel(as.factor(completeCases$Sex), "M")
```


```r
#mean FD is different for lenient and strict because some frames at the beginning and/or end of the scan are excluded for scans to pass lenient motion QC
completeCases$MeanFramewiseDisplacement[completeCases$Motion.Exclusion.Level=="Lenient"]=completeCases$MeanFramewiseDisplacement.KKI[completeCases$Motion.Exclusion.Level=="Lenient"]
```


### 3.1.2 mean framewise displacement as a function of phenotype and age (controlling for sex and SES)
#### Covariates

```r
phenoVariables <- c("ID", "PrimaryDiagnosis", 
                    "ADOS.Comparable.Total", 
                    "SRS.Score", 
                    "PANESS.TotalOverflowNotAccountingForAge", 
                    "DuPaulHome.InattentionRaw",
                    "DuPaulHome.HyperactivityRaw", 
                    "AgeAtScan", 
                    "WISC.GAI", "SES.Family",
                    "Motion.Exclusion.Level", "Included", "Sex",
                    "MeanFramewiseDisplacement")

phenoIDs <- c("ID", "PrimaryDiagnosis", "Motion.Exclusion.Level", "Included", "Sex", "SES.Family", "MeanFramewiseDisplacement")

aim1 <- reshape2::melt(completeCases[, phenoVariables],
                         id.vars=names(completeCases)[which(names(completeCases) %in% phenoIDs)])

levels(aim1$variable) <- c("ADOS", "SRS", "Motor Overflow", "Inattention", 
                           "Hyperactivity", "Age", "GAI")

aim1G <- group_by(aim1, PrimaryDiagnosis, Motion.Exclusion.Level, Included, variable)
```

#### Fit univarate GAMs using data from all children (those with usable rsfMRI data and those with unusable rsfMRI data)


```r
#For Motion.Exclusion.Level=None, all children are included
meanFDtib <- tibble(aim1) %>% 
  filter(Included=="Included")

meanFDNest <- meanFDtib %>% 
  group_by(variable, Motion.Exclusion.Level) %>% 
  tidyr::nest()

#nested models
nested_gams_all <- meanFDNest %>%
  mutate(model = map(data, ~mgcv::gam(MeanFramewiseDisplacement ~ s(value, k=-1),
                                      data = na.omit(.x), 
                                      method="REML")),
         coefs = map(model, tidy, conf.int = FALSE),
         Rsq = map_dbl(model, ~summary(.)$r.sq)) %>% 
  unnest(coefs)

nested_gams_all <- nested_gams_all %>% 
           mutate(LB = map(data, ~round(min(na.omit(.x$value)))),
                  UB = map(data, ~round(max(na.omit(.x$value)))),
                  range = map2(LB, UB, ~seq(from=.x, to=.y, by=.1)),
                  fit.data = map(range, ~data.frame(value = .x)),
                  FDpredict = map2(model, fit.data, ~predict(.x, newdata = .y, se.fit=TRUE)),
                  fit = map(FDpredict, ~.x$fit),
                  lCI = map(FDpredict, ~(.x$fit-1.96*.x$se.fit)),
                  hCI = map(FDpredict, ~(.x$fit+1.96*.x$se.fit)))
```

#### Define theme for Figure S2 top row

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

#### Figure S2a top. mean FD as a function of ADOS


```r
nested_gams <- filter(nested_gams_all, Motion.Exclusion.Level=="None")

nested_gams$Motion.Exclusion.Level <- droplevels(nested_gams$Motion.Exclusion.Level)

maxX = nested_gams %>% 
  ungroup(Motion.Exclusion.Level) %>% 
  select("variable", "LB", "UB")

ados <- nested_gams %>% 
  filter(variable=="ADOS") %>% 
  select("Motion.Exclusion.Level", "variable", "range", "fit", "lCI", "hCI", "LB", "UB") %>% 
  unnest(c(range, fit, lCI, hCI, LB, UB))


p_ados <- ggplot(ados, aes(x=range, y=fit))+
  geom_line(aes(),size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), linetype='blank', alpha=0.2)+  
  labs(x='', y='mean FD')+  
  scale_x_continuous(expand = c(0, 0), limits = c(ados$LB[1], ados$UB[1]))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.75), breaks=seq(0, 1.75, by=.25))+
  gam_theme+
  ggtitle("ADOS")+
  theme(plot.title = element_text(size = 11, hjust = 0.5))

#top_title = cowplot::get_title(p_ados + theme(plot.title = element_text(size = 11, hjust = 0.5)))
```

#### Figure S2b top. mean FD as a function of SRS

```r
srs <- nested_gams %>% 
  filter(variable=="SRS") %>% 
  select("Motion.Exclusion.Level", "variable","range", "fit", 'lCI', 'hCI') %>% 
  unnest(c(range, fit, lCI, hCI))

p_srs <- ggplot(srs, aes(x=range, y=fit))+
  geom_line(size=1.2)+
  theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), linetype='blank', alpha=0.2)+  
  labs(x='', y='')+  
  scale_x_continuous(expand = c(0, 0))+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.75), breaks=seq(0, 1.75, by=.25))+
  gam_theme+
  ggtitle("SRS")+
  theme(plot.title = element_text(size = 11, hjust = 0.5))+
  theme(axis.title.y = element_blank())
```

#### Figure S2c top. mean FD as a function of Inattention

```r
ina <- nested_gams %>% 
  filter(variable=="Inattention") %>% 
  select("Motion.Exclusion.Level","variable", "range", "fit", 'lCI', 'hCI') %>% 
  unnest(c(range, fit, lCI, hCI))

p_in <- ggplot(ina, aes(x=range, y=fit))+
  geom_line(size=1.2)+
  theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), linetype='blank', alpha=0.2)+  
  labs(x='', y='')+  
  scale_x_continuous(expand = c(0, 0))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.75), breaks=seq(0, 1.75, by=.25))+
  gam_theme+
  ggtitle("Inattention")+
  theme(plot.title = element_text(size = 11, hjust = 0.5))+
  theme(axis.title.y = element_blank())
```

#### Figure S2d top. mean FD as a function of Hyperactivity

```r
hi <- nested_gams %>% 
  filter(variable=="Hyperactivity") %>% 
  select("Motion.Exclusion.Level", "variable",  "range", "fit", 'lCI', 'hCI') %>% 
  unnest(c(range, fit, lCI, hCI))

p_hi <- ggplot(hi, aes(x=range, y=fit))+
  geom_line(size=1.2)+
  theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), linetype='blank', alpha=0.2)+  
  labs(x='', y='')+  
  scale_x_continuous(expand = c(0, 0))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.75), breaks=seq(0, 1.75, by=.25))+
  gam_theme+
  ggtitle("Hyperactivity")+
  theme(plot.title = element_text(size = 11, hjust = 0.5))+
  theme(axis.title.y = element_blank())
```

#### Figure S2e top. mean FD as a function of Motor Overflow

```r
mo <- nested_gams %>% 
  filter(variable=="Motor Overflow") %>% 
  select("Motion.Exclusion.Level", "variable", "range", "fit", 'lCI', 'hCI') %>% 
  unnest(c(range, fit, lCI, hCI))

p_mo <- ggplot(mo, aes(x=range, y=fit))+
  geom_line(size=1.2)+
  theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), linetype='blank', alpha=0.2)+  
  labs(x='', y='')+  
  scale_x_continuous(expand = c(0, 0))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.75), breaks=seq(0, 1.75, by=.25))+
  gam_theme+
  ggtitle("Motor Overflow")+
  theme(plot.title = element_text(size = 11, hjust = 0.5))+
  theme(axis.title.y = element_blank())
```

#### Figure S2f top. mean FD as a function of Age

```r
age<- nested_gams %>% 
  filter(variable=="Age") %>% 
  select("Motion.Exclusion.Level", "variable", "range", "fit", 'lCI', 'hCI') %>% 
  unnest(c(range, fit, lCI, hCI))

p_age <- ggplot(age, aes(x=range, y=fit))+
  geom_line(size=1.2)+
  theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), linetype='blank', alpha=0.2)+  
  labs(x='', y='')+  
  scale_x_continuous(expand = c(0, 0))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.75), breaks=seq(0, 1.75, by=.25))+
  gam_theme+
  ggtitle("Age")+
  theme(plot.title = element_text(size = 11, hjust = 0.5))+
  theme(axis.title.y = element_blank())
```

#### Figure S2g top. mean FD as a function of GAI

```r
gai <- nested_gams %>% 
  filter(variable=="GAI") %>% 
  select("Motion.Exclusion.Level", "variable", "range", "fit", 'lCI', 'hCI') %>% 
  unnest(c(range, fit, lCI, hCI))

p_gai <- ggplot(gai, aes(x=range, y=fit))+
  geom_line(size=1.2)+
  theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), linetype='blank', alpha=0.2)+  
  labs(x='', y='')+  
  scale_x_continuous(expand = c(0, 0))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.75), breaks=seq(0, 1.75, by=.25))+
  gam_theme+
  ggtitle("GAI")+
  theme(plot.title = element_text(size = 11, hjust = 0.5))+
  theme(axis.title.y = element_blank())

p_legend = cowplot::get_legend(p_gai + guides(color = guide_legend(nrow = 1))+
                                 theme(legend.position = "bottom", legend.text = element_text(size = 11),
                                       legend.key.size=unit(.15, "in")))
```


#### Figure S2 second row: Density plots of covariates used to fit GAMs (across included & excluded children) 
#### Define theme for density plots of covariates across included and excluded children


#### Figure S2a second row ADOS density

```r
ddata <- nested_gams %>% 
  filter(variable=="ADOS") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB", "UB") %>% 
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
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB[1], ddata$UB[1]), breaks=seq(0, ddata$UB[1], by=5))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0, .1), breaks=seq(0, .1, by=.05))+  
  labs(x='', y='Density')+
  scale_fill_manual(values = c("#FDE599"))+ 
  scale_color_manual(values = c("#E9D38D"))+ 
  den_theme
```

#### Figure S2b second row. SRS density

```r
ddata <- nested_gams %>% 
  filter(variable=="SRS") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB", "UB") %>% 
  unnest(c(data, LB, UB))

d_srs=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits=c(0, ddata$UB[1]),breaks = seq(0, 100 , by = 50))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .05), breaks=seq(0, .06, by=.02))+  
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```

#### Figure S2c second row. Inattention density

```r
ddata <- nested_gams %>% 
  filter(variable=="Inattention") %>% 
  select("Motion.Exclusion.Level", "variable", "data") %>% 
  unnest(data)

d_in=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .18), breaks=seq(0, .18, by=.06))+ 
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

#### Figure S2d second row. Hyperactivity/Impulsivity Density

```r
ddata <- nested_gams %>% 
  filter(variable=="Hyperactivity") %>% 
  select("Motion.Exclusion.Level", "variable", "data") %>% 
  unnest(data)

d_hi=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .3), breaks=seq(0, .3, by=.15))+ 
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```

#### Figure S2e second row. Motor Overflow Density

```r
ddata <- nested_gams %>% 
  filter(variable=="Motor Overflow") %>% 
  select("Motion.Exclusion.Level", "variable", "data") %>% 
  unnest(data)

d_mo=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .08), breaks=seq(0, .08, by=.04))+ 
  theme(axis.title.y = element_blank())+
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```

#### Figure S2f second row. Age Density

```r
ddata <- nested_gams %>% 
  filter(variable=="Age") %>% 
  select("Motion.Exclusion.Level", "variable", "data") %>% 
  unnest(data)

d_age=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits=c(8,13), breaks = seq(8, 13 , by = 1))+  
  scale_y_continuous(expand = c(0, 0), limits=c(0,.32), breaks=seq(0, .3, by = .16))+ 
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```

#### Figure S2g second row. GAI Density

```r
ddata <- nested_gams %>% 
  filter(variable=="GAI") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB", "UB") %>% 
  unnest(data, LB, UB)
```

```
## Warning: unnest() has a new interface. See ?unnest for details.
## Try `df %>% unnest(c(data, LB, UB))`, with `mutate()` if needed
```

```r
d_gai=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB[1], ddata$UB))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0, .04), breaks=seq(0., .04, by=.02))+  
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

#### Figure S2a third row. mean FD as a function of ADOS

```r
nested_gams <- filter(nested_gams_all, Motion.Exclusion.Level=="Lenient")
nested_gams$Motion.Exclusion.Level <- droplevels(nested_gams$Motion.Exclusion.Level)

nested_gams <- merge(nested_gams, maxX, by="variable")

temp <- nested_gams %>% 
  filter(variable=="ADOS") %>% 
  select("variable", "range", "fit", "lCI", "hCI", "LB.y", "UB.y") %>% 
  unnest(c(range, fit, lCI, hCI, LB.y, UB.y))


p_ados_l <- ggplot(temp, aes(x=range, y=fit))+
  geom_line(aes(), color = "#9FB0CC", size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), fill = "#9FB0CC", linetype='blank', alpha=0.2)+  
  labs(x='', y='mean FD')+  
  scale_x_continuous(expand = c(0, 0), limits = c(temp$LB.y[1], temp$UB.y[1]))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, .5), breaks=seq(0, .5, by=.25))+
  gam_theme+
  ggtitle("Lenient")+
  theme(plot.title = element_blank())
  

second_title = cowplot::get_title(p_ados_l + theme(plot.title = element_text(size = 11, hjust = 0.5)))
```

#### SRS

```r
temp <- nested_gams %>% 
  filter(variable=="SRS") %>% 
  select("variable", "range", "fit", "lCI", "hCI", "LB.y", "UB.y") %>% 
  unnest(c(range, fit, lCI, hCI, LB.y, UB.y))


p_srs_l <- ggplot(temp, aes(x=range, y=fit))+
  geom_line(aes(), color = "#9FB0CC", size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), fill = "#9FB0CC", linetype='blank', alpha=0.2)+ 
  labs(x='', y='mean FD')+  
  scale_x_continuous(expand = c(0, 0), limits = c(temp$LB.y[1], temp$UB.y[1]))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, .5), breaks=seq(0, .5, by=.25))+
  gam_theme+
  theme(plot.title = element_blank())+
  theme(axis.title.y = element_blank())
```

#### Inattention

```r
temp <- nested_gams %>% 
  filter(variable=="Inattention") %>% 
  select("variable", "range", "fit", "lCI", "hCI", "LB.y", "UB.y") %>% 
  unnest(c(range, fit, lCI, hCI, LB.y, UB.y))


p_in_l <- ggplot(temp, aes(x=range, y=fit))+
  geom_line(aes(), color = "#9FB0CC", size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), fill = "#9FB0CC", linetype='blank', alpha=0.2)+ 
  labs(x='', y='mean FD')+  
  scale_x_continuous(expand = c(0, 0), limits = c(temp$LB.y[1], temp$UB.y[1]))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, .5), breaks=seq(0, .5, by=.25))+
  gam_theme+
  theme(plot.title = element_blank())+
  theme(axis.title.y = element_blank())
```

#### Hyperactivity

```r
temp <- nested_gams %>% 
  filter(variable=="Hyperactivity") %>% 
  select("variable", "range", "fit", "lCI", "hCI", "LB.y", "UB.y") %>% 
  unnest(c(range, fit, lCI, hCI, LB.y, UB.y))


p_hi_l <- ggplot(temp, aes(x=range, y=fit))+
  geom_line(aes(), color = "#9FB0CC", size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), fill = "#9FB0CC", linetype='blank', alpha=0.2)+ 
  labs(x='', y='mean FD')+  
  scale_x_continuous(expand = c(0, 0), limits = c(temp$LB.y[1], temp$UB.y[1]))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, .5), breaks=seq(0, .5, by=.25))+
  gam_theme+
  theme(plot.title = element_blank())+
  theme(axis.title.y = element_blank())
```


#### Motor Overflow

```r
temp <- nested_gams %>% 
  filter(variable=="Motor Overflow") %>% 
  select("variable", "range", "fit", "lCI", "hCI", "LB.y", "UB.y") %>% 
  unnest(c(range, fit, lCI, hCI, LB.y, UB.y))


p_mo_l <- ggplot(temp, aes(x=range, y=fit))+
  geom_line(aes(), color = "#9FB0CC", size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), fill = "#9FB0CC", linetype='blank', alpha=0.2)+ 
  labs(x='', y='mean FD')+  
  scale_x_continuous(expand = c(0, 0), limits = c(temp$LB.y[1], temp$UB.y[1]))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, .5), breaks=seq(0, .5, by=.25))+
  gam_theme+
  theme(plot.title = element_blank())+
  theme(axis.title.y = element_blank())
```

#### Age

```r
temp <- nested_gams %>% 
  filter(variable=="Age") %>% 
  select("variable", "range", "fit", "lCI", "hCI", "LB.y", "UB.y") %>% 
  unnest(c(range, fit, lCI, hCI, LB.y, UB.y))


p_age_l <- ggplot(temp, aes(x=range, y=fit))+
  geom_line(aes(), color = "#9FB0CC", size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), fill = "#9FB0CC", linetype='blank', alpha=0.2)+ 
  labs(x='', y='mean FD')+  
  scale_x_continuous(expand = c(0, 0), limits = c(temp$LB.y[1], temp$UB.y[1]))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, .5), breaks=seq(0, .5, by=.25))+
  gam_theme+
  theme(plot.title = element_blank())+
  theme(axis.title.y = element_blank())
```

#### GAI

```r
temp <- nested_gams %>% 
  filter(variable=="GAI") %>% 
  select("variable", "range", "fit", "lCI", "hCI", "LB.y", "UB.y") %>% 
  unnest(c(range, fit, lCI, hCI, LB.y, UB.y))


p_gai_l <- ggplot(temp, aes(x=range, y=fit))+
  geom_line(aes(), color = "#9FB0CC", size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), fill = "#9FB0CC", linetype='blank', alpha=0.2)+ 
  labs(x='', y='mean FD')+  
  scale_x_continuous(expand = c(0, 0), limits = c(temp$LB.y[1], temp$UB.y[1]))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, .5), breaks=seq(0, .5, by=.25))+
  scale_color_manual("#9FB0CC")+
  scale_fill_manual("#9FB0CC")+
  gam_theme+
  theme(plot.title = element_blank())+
  theme(axis.title.y = element_blank())
```


#### Figure S2a bottom. ADOS density

```r
ddata <- nested_gams %>% 
  filter(variable=="ADOS") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB.y", "UB.y") %>% 
  unnest(c(data, LB.y, UB.y)) %>% 
  filter(PrimaryDiagnosis=="ASD")

ddata$PrimaryDiagnosis <- droplevels(ddata$PrimaryDiagnosis)

d_ados_l=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB.y[1], ddata$UB.y[1]))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0, .1), breaks=seq(0, .1, by=.05))+  
  labs(x='', y='Density')+
  scale_fill_manual(values = c("#FDE599"))+ 
  scale_color_manual(values = c("#E9D38D"))+ 
  den_theme
```

#### Figure S2b bottom. SRS density

```r
ddata <- nested_gams %>% 
  filter(variable=="SRS") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB.y", "UB.y") %>% 
  unnest(c(data, LB.y, UB.y))

d_srs_l=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB.y[1], ddata$UB.y[1]))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .05), breaks=seq(0, .06, by=.02))+  
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```


```r
ddata <- nested_gams %>% 
  filter(variable=="Inattention") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB.y", "UB.y") %>% 
  unnest(c(data, LB.y, UB.y))

d_in_l=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB.y[1], ddata$UB.y[1]))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .18), breaks=seq(0, .18, by=.06))+ 
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```


```r
ddata <- nested_gams %>% 
  filter(variable=="Hyperactivity") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB.y", "UB.y") %>% 
  unnest(c(data, LB.y, UB.y))

d_hi_l=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB.y[1], ddata$UB.y[1]))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .3), breaks=seq(0, .3, by=.15))+   
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```


```r
ddata <- nested_gams %>% 
  filter(variable=="Motor Overflow") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB.y", "UB.y") %>% 
  unnest(c(data, LB.y, UB.y))

d_mo_l=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB.y[1], ddata$UB.y[1]))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .08), breaks=seq(0, .08, by=.04))+  
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```


```r
ddata <- nested_gams %>% 
  filter(variable=="Age") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB.y", "UB.y") %>% 
  unnest(c(data, LB.y, UB.y))

d_age_l=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB.y[1], ddata$UB.y[1]))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .32), breaks=seq(0, .3, by=.16))+ 
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```


```r
ddata <- nested_gams %>% 
  filter(variable=="GAI") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB.y", "UB.y") %>% 
  unnest(c(data, LB.y, UB.y))

d_gai_l=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB.y[1], ddata$UB.y[1]))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .04), breaks=seq(0, .04, by=.02))+  
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```

#### Repeat for strict motion QC

```r
nested_gams <- filter(nested_gams_all, Motion.Exclusion.Level=="Strict")
nested_gams$Motion.Exclusion.Level <- droplevels(nested_gams$Motion.Exclusion.Level)

nested_gams <- merge(nested_gams, maxX, by="variable")

ados <- nested_gams %>% 
  filter(variable=="ADOS") %>% 
  select("variable", "range", "fit", "lCI", "hCI", "LB.y", "UB.y") %>% 
  unnest(c(range, fit, lCI, hCI, LB.y, UB.y))


p_ados_s <- ggplot(ados, aes(x=range, y=fit))+
  geom_line(aes(), color = "#f55154", size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), fill = "#f55154", linetype='blank', alpha=0.2)+  
  labs(x='', y='mean FD')+  
  scale_x_continuous(expand = c(0, 0), limits = c(ados$LB.y[1], ados$UB.y[1]))+ 
  scale_y_continuous(expand = c(0, 0),  limits = c(.08, .18), breaks=seq(.08, .18, by=.02))+ 
  gam_theme+
  ggtitle("Strict")+
  theme(plot.title = element_blank())

third_title = cowplot::get_title(p_ados_s + theme(plot.title = element_text(size = 11, hjust = 0.5)))
```

#### SRS

```r
temp <- nested_gams %>% 
  filter(variable=="SRS") %>% 
  select("variable", "range", "fit", "lCI", "hCI", "LB.y", "UB.y") %>% 
  unnest(c(range, fit, lCI, hCI, LB.y, UB.y))


p_srs_s <- ggplot(temp, aes(x=range, y=fit))+
  geom_line(aes(), color = "#f55154", size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), fill = "#f55154", linetype='blank', alpha=0.2)+  
  labs(x='', y='mean FD')+  
  scale_x_continuous(expand = c(0, 0), limits = c(temp$LB.y[1], temp$UB.y[1]))+ 
  scale_y_continuous(expand = c(0, 0),  limits = c(.08, .18), breaks=seq(.08, .18, by=.02))+ 
  gam_theme+
  theme(plot.title = element_blank())+
  theme(axis.title.y = element_blank())
```

#### Inattention

```r
temp <- nested_gams %>% 
  filter(variable=="Inattention") %>% 
  select("variable", "range", "fit", "lCI", "hCI", "LB.y", "UB.y") %>% 
  unnest(c(range, fit, lCI, hCI, LB.y, UB.y))


p_in_s <- ggplot(temp, aes(x=range, y=fit))+
  geom_line(aes(), color = "#f55154", size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), fill = "#f55154", linetype='blank', alpha=0.2)+  
  labs(x='', y='')+  
  scale_x_continuous(expand = c(0, 0), limits = c(temp$LB.y[1], temp$UB.y[1]))+ 
  scale_y_continuous(expand = c(0, 0),  limits = c(.08, .18), breaks=seq(.08, .18, by=.02))+ 
  gam_theme+
  theme(plot.title = element_blank())+
  theme(axis.title.y = element_blank())
```

#### Hyperactivity

```r
temp <- nested_gams %>% 
  filter(variable=="Hyperactivity") %>% 
  select("variable", "range", "fit", "lCI", "hCI", "LB.y", "UB.y") %>% 
  unnest(c(range, fit, lCI, hCI, LB.y, UB.y))


p_hi_s <- ggplot(temp, aes(x=range, y=fit))+
  geom_line(aes(), color = "#f55154", size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), fill = "#f55154", linetype='blank', alpha=0.2)+  
  labs(x='', y='')+  
  scale_x_continuous(expand = c(0, 0), limits = c(temp$LB.y[1], temp$UB.y[1]))+ 
  scale_y_continuous(expand = c(0, 0),  limits = c(.08, .18), breaks=seq(.08, .18, by=.02))+ 
  gam_theme+
  theme(plot.title = element_blank())+
  theme(axis.title.y = element_blank())
```


#### Motor Overflow

```r
temp <- nested_gams %>% 
  filter(variable=="Motor Overflow") %>% 
  select("variable", "range", "fit", "lCI", "hCI", "LB.y", "UB.y") %>% 
  unnest(c(range, fit, lCI, hCI, LB.y, UB.y))


p_mo_s <- ggplot(temp, aes(x=range, y=fit))+
  geom_line(aes(), color = "#f55154", size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), fill = "#f55154", linetype='blank', alpha=0.2)+  
  labs(x='', y='')+  
  scale_x_continuous(expand = c(0, 0), limits = c(temp$LB.y[1], temp$UB.y[1]))+ 
  scale_y_continuous(expand = c(0, 0),  limits = c(.08, .18), breaks=seq(.08, .18, by=.02))+ 
  gam_theme+
  theme(plot.title = element_blank())+
  theme(axis.title.y = element_blank())
```

#### Age

```r
temp <- nested_gams %>% 
  filter(variable=="Age") %>% 
  select("variable", "range", "fit", "lCI", "hCI", "LB.y", "UB.y") %>% 
  unnest(c(range, fit, lCI, hCI, LB.y, UB.y))


p_age_s <- ggplot(temp, aes(x=range, y=fit))+
  geom_line(aes(), color = "#f55154", size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), fill = "#f55154", linetype='blank', alpha=0.2)+  
  labs(x='', y='mean FD')+  
  scale_x_continuous(expand = c(0, 0), limits = c(temp$LB.y[1], temp$UB.y[1]))+ 
  scale_y_continuous(expand = c(0, 0),  limits = c(.08, .18), breaks=seq(.08, .18, by=.02))+ 
  gam_theme+
  theme(plot.title = element_blank())+
  theme(axis.title.y = element_blank())
```

#### GAI

```r
temp <- nested_gams %>% 
  filter(variable=="GAI") %>% 
  select("variable", "range", "fit", "lCI", "hCI", "LB.y", "UB.y") %>% 
  unnest(c(range, fit, lCI, hCI, LB.y, UB.y))


p_gai_s <- ggplot(temp, aes(x=range, y=fit))+
  geom_line(aes(), color = "#f55154", size=1.2)+theme_bw()+
  geom_ribbon(aes(ymin=lCI, ymax=hCI), fill = "#f55154", linetype='blank', alpha=0.2)+  
  labs(x='', y='')+  
  scale_x_continuous(expand = c(0, 0), limits = c(temp$LB.y[1], temp$UB.y[1]))+ 
  scale_y_continuous(expand = c(0, 0),  limits = c(.08, .18), breaks=seq(.08, .18, by=.02))+ 
  scale_color_manual("#9FB0CC")+
  scale_fill_manual("#9FB0CC")+
  gam_theme+
  theme(plot.title = element_blank())+
  theme(axis.title.y = element_blank())
```

#### Figure S2a bottom. ADOS density

```r
ddata <- nested_gams %>% 
  filter(variable=="ADOS") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB.y", "UB.y") %>% 
  unnest(c(data, LB.y, UB.y)) %>% 
  filter(PrimaryDiagnosis=="ASD")

ddata$PrimaryDiagnosis <- droplevels(ddata$PrimaryDiagnosis)

d_ados_s=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB.y[1], ddata$UB.y[1]))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0, .1), breaks=seq(0, .1, by=.05))+  
  labs(x='', y='Density')+
  scale_fill_manual(values = c("#FDE599"))+ 
  scale_color_manual(values = c("#E9D38D"))+ 
  den_theme
```

#### Figure S2b bottom. SRS density

```r
ddata <- nested_gams %>% 
  filter(variable=="SRS") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB.y", "UB.y") %>% 
  unnest(c(data, LB.y, UB.y))

d_srs_s=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB.y[1], ddata$UB.y[1]))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .05), breaks=seq(0, .06, by=.02))+  
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```


```r
ddata <- nested_gams %>% 
  filter(variable=="Inattention") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB.y", "UB.y") %>% 
  unnest(c(data, LB.y, UB.y))

d_in_s=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB.y[1], ddata$UB.y[1]))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .18), breaks=seq(0, .18, by=.06))+ 
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```


```r
ddata <- nested_gams %>% 
  filter(variable=="Hyperactivity") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB.y", "UB.y") %>% 
  unnest(c(data, LB.y, UB.y))

d_hi_s=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB.y[1], ddata$UB.y[1]))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .3), breaks=seq(0, .3, by=.15))+   
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```


```r
ddata <- nested_gams %>% 
  filter(variable=="Motor Overflow") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB.y", "UB.y") %>% 
  unnest(c(data, LB.y, UB.y))

d_mo_s=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+ 
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB.y[1], ddata$UB.y[1]))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .08), breaks=seq(0, .08, by=.04))+  
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```


```r
ddata <- nested_gams %>% 
  filter(variable=="Age") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB.y", "UB.y") %>% 
  unnest(c(data, LB.y, UB.y))

d_age_s=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB.y[1], ddata$UB.y[1]))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .32), breaks=seq(0, .3, by=.16))+ 
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```


```r
ddata <- nested_gams %>% 
  filter(variable=="GAI") %>% 
  select("Motion.Exclusion.Level", "variable", "data", "LB.y", "UB.y") %>% 
  unnest(c(data, LB.y, UB.y))

d_gai_s=ggplot(ddata, aes(x=value, fill=PrimaryDiagnosis, color=PrimaryDiagnosis))+  
  geom_density(alpha=0.5, inherit.aes=TRUE, trim=TRUE)+  
  scale_x_continuous(expand = c(0, 0), limits = c(ddata$LB.y[1], ddata$UB.y[1]))+  
  scale_y_continuous(expand = c(0, 0),  limits = c(0, .04), breaks=seq(0, .04, by=.02))+  
  scale_fill_manual(labels=c('TD','ASD'), values = c("#009E73", "#FDE599"))+
  scale_color_manual(labels=c('TD','ASD'), values = c("#05634a", "#E9D38D"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(size=.5), panel.border = element_blank())+  
  den_theme+
  theme(axis.title.y = element_blank())+
  labs(x='')
```

#### combine gam plots with densities & print

```r
top_row <- cowplot::plot_grid(p_ados, p_srs, p_in, p_hi, p_mo, p_age, p_gai, ncol=7, 
                              rel_widths=c(1.18/7, .97/6, .97/6, .97/6, .97/6, .97/6, .97/6))
second_row <- cowplot::plot_grid(d_ados, d_srs, d_in, d_hi, d_mo, d_age, d_gai, ncol=7, 
                                 rel_widths=c(1.18/7, .97/6, .97/6, .97/6, .97/6, .97/6, .97/6))
```

```
## Warning: Removed 89 rows containing non-finite values (stat_density).
```

```r
third_row <- cowplot::plot_grid(p_ados_l, p_srs_l, p_in_l, p_hi_l, p_mo_l, p_age_l, p_gai_l, ncol=7, 
                                 rel_widths=c(1.18/7, .97/6, .97/6, .97/6, .97/6, .97/6, .97/6))

fourth_row <- cowplot::plot_grid(d_ados_l, d_srs_l, d_in_l, d_hi_l, d_mo_l, d_age_l, d_gai_l, ncol=7, 
                                 rel_widths=c(1.18/7, .97/6, .97/6, .97/6, .97/6, .97/6, .97/6))
```

```
## Warning: Removed 74 rows containing non-finite values (stat_density).
```

```r
fifth_row <- cowplot::plot_grid(p_ados_s, p_srs_s, p_in_s, p_hi_s, p_mo_s, p_age_s, p_gai_s, ncol=7, 
                                 rel_widths=c(1.18/7, .97/6, .97/6, .97/6, .97/6, .97/6, .97/6))

sixth_row <- cowplot::plot_grid(d_ados_s, d_srs_s, d_in_s, d_hi_s, d_mo_s, d_age_s, d_gai_s, ncol=7, 
                                 rel_widths=c(1.18/7, .97/6, .97/6, .97/6, .97/6, .97/6, .97/6))
```

```
## Warning: Removed 37 rows containing non-finite values (stat_density).
```

```r
png("./CovariatesAndRS-fMRIUsability/fig_meanFD_allGAM_TD_ASD_cc_3levels.png",width=11,height=11,units="in",res=200)
cowplot::plot_grid(top_row, NULL, 
                   second_row, NULL, 
                   third_row, NULL, 
                   fourth_row, NULL, 
                   fifth_row, NULL,
                   sixth_row, NULL,
                   hist_legend, nrow=13, rel_heights=c( 1, -.02, 
                                                        .5, -.05, 
                                                        1, -.02, 
                                                        .5, -.05,
                                                        1, -.02,
                                                       .5, -.07, .1))
dev.off()
```

```
## quartz_off_screen 
##                 2
```


```r
nested_gams_none <- nested_gams_all %>% 
  filter(Motion.Exclusion.Level=="None")

nested_gams_none$p.fdr = p.adjust(nested_gams_none$p.value, method = "BH")

nested_gams_len <-  nested_gams_all %>% 
  filter(Motion.Exclusion.Level=="Lenient")

nested_gams_len$p.fdr = p.adjust(nested_gams_len$p.value, method = "BH")

nested_gams_strict <-  nested_gams_all %>% 
  filter(Motion.Exclusion.Level=="Strict")
  
nested_gams_strict$p.fdr = p.adjust(nested_gams_strict$p.value, method = "BH")
#combine to print
nested_gams_table <- rbind(nested_gams_none, nested_gams_len, nested_gams_strict)

nested_gams_table[nested_gams_table$Motion.Exclusion.Level=="None", c(1:2,6:10, 19)]
```

```
## # A tibble: 7 × 8
## # Groups:   variable, Motion.Exclusion.Level [7]
##   Motion.Exclusion.Level variable    edf ref.df statistic p.value    Rsq   p.fdr
##   <fct>                  <fct>     <dbl>  <dbl>     <dbl>   <dbl>  <dbl>   <dbl>
## 1 None                   ADOS       2.35   2.85     23.2  0       0.121  0      
## 2 None                   SRS        5.00   6.12      8.56 0       0.118  0      
## 3 None                   Motor Ov…  1.00   1.00     46.2  0       0.0867 0      
## 4 None                   Inattent…  3.45   4.26      9.91 0       0.0813 0      
## 5 None                   Hyperact…  1.44   1.75     16.3  1.55e-6 0.0600 2.17e-6
## 6 None                   Age        1.00   1.00     17.2  3.98e-5 0.0329 4.65e-5
## 7 None                   GAI        1.30   1.55      7.83 4.11e-3 0.0205 4.11e-3
```

```r
nested_gams_table[nested_gams_table$Motion.Exclusion.Level=="Lenient", c(1:2,6:10, 19)]
```

```
## # A tibble: 7 × 8
## # Groups:   variable, Motion.Exclusion.Level [7]
##   Motion.Exclusion.Lev… variable   edf ref.df statistic p.value      Rsq   p.fdr
##   <fct>                 <fct>    <dbl>  <dbl>     <dbl>   <dbl>    <dbl>   <dbl>
## 1 Lenient               ADOS      1.00   1.00    22.6   3.07e-6  0.0526  1.07e-5
## 2 Lenient               SRS       1.00   1.00    26.0   8.64e-7  0.0735  6.05e-6
## 3 Lenient               Motor O…  1.00   1.00    19.0   1.70e-5  0.0442  3.98e-5
## 4 Lenient               Inatten…  1.00   1.00    16.4   6.12e-5  0.0382  8.56e-5
## 5 Lenient               Hyperac…  2.03   2.53     8.77  5.42e-5  0.0540  8.56e-5
## 6 Lenient               Age       2.95   3.68     5.57  3.84e-4  0.0492  4.48e-4
## 7 Lenient               GAI       1.05   1.10     0.558 4.96e-1 -0.00105 4.96e-1
```

```r
nested_gams_table[nested_gams_table$Motion.Exclusion.Level=="Strict", c(1:2,6:10, 19)]
```

```
## # A tibble: 7 × 8
## # Groups:   variable, Motion.Exclusion.Level [7]
##   Motion.Exclusion.Level variable    edf ref.df statistic p.value      Rsq p.fdr
##   <fct>                  <fct>     <dbl>  <dbl>     <dbl>   <dbl>    <dbl> <dbl>
## 1 Strict                 ADOS       1.00   1.00  0.000235  0.991  -0.00613 0.991
## 2 Strict                 SRS        1.14   1.27  0.306     0.572  -0.00250 0.801
## 3 Strict                 Motor Ov…  3.25   4.06  1.35      0.273   0.0294  0.801
## 4 Strict                 Inattent…  1.00   1.00  0.785     0.377  -0.00131 0.801
## 5 Strict                 Hyperact…  1.00   1.00  0.00517   0.946  -0.00610 0.991
## 6 Strict                 Age        1.32   1.58  0.881     0.517   0.00121 0.801
## 7 Strict                 GAI        1.00   1.00  3.73      0.0550  0.0164  0.385
```

#### print figure for report

```r
cowplot::plot_grid(top_row, NULL, 
                   second_row, NULL, 
                   third_row, NULL, 
                   fourth_row, NULL, 
                   fifth_row, NULL,
                   sixth_row, NULL,
                   hist_legend, nrow=13, rel_heights=c( 1, -.02, 
                                                        .5, -.05, 
                                                        1, -.02, 
                                                        .5, -.05,
                                                        1, -.02,
                                                       .5, -.07, .1))
```

<img src="forReviewers_meanFD_covariates//unnamed-chunk-47-1.png" width="1056" />


**Figure. mean FD changes with phenotype and age.** Univariate analysis of mean framewise displacement as a function of participant characteristics}. From left to right: Autism Diagnostic Observation Schedule (ADOS) total scores, social responsiveness scale (SRS) scores, inattentive symptoms, hyperactive/impulsive symptoms, total motor overflow, age, and general ability index (GAI) for children with usable and unusable rsfMRI data (top row, black lines), children with usable data under lenient motion QC (slate blue lines), and children with usable data under strict motion QC (red lines). Variable distributions for each diagnosis group are displayed across the bottom panel (TD=typically developing, green; ASD=autism spectrum disorder, yellow).

