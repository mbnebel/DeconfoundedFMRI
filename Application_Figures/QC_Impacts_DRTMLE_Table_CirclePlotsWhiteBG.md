---
title: "Tables and Circle Plots comparing Naive and DRTMLE group differences"
author: "MB Nebel"
date: '2022-05-02'
output: html_document
---




```r
library(tidyverse)
#for tableby
library(arsenal)
library(xtable)
library(readxl)
library(circlize)
library(png)
#to remove background from .pngs
library(colordistance)
library(magick)
```

#### Define Circle plot function:
Based on [Dr. Mowinckel's](https://github.com/Athanasiamo) [tutorial](https://drmowinckels.io/blog/2018-05-25-circluar-plots-in-r-and-adding-images/) 
and adapted by [Daniel Lidstone](https://github.com/lidstone) & [Ben Risk](https://github.com/benjaminrisk) 

```r
#### Function to Make Transparent Colors 
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

# varname: a variable name corresponding to a column of results.ave
# thresh: threshold to use for displaying chords

myChordDiagram = function(data,varname,alpha=0.20,alphaname,title) {
  # alpha: pvalue threshold
  # alphaname: column of data with pvalues
  
  #
  
  ## Set transparent color, used to mask unwanted chords
  trcol <- t_col("white", perc = 100, name = "white")
  
  ## Set colors for negative and positive values
  poscol <- t_col("DodgerBlue4", perc = 15, name = "lt.blue") #positive correlation color
  negcol <- t_col("red", perc = 15, name = "lt.red") #negative correlation color
  
  ## Set names for networks (netLab)
  icInfo <- read_excel("componentLabels_pca85_ica30.xlsx")
  netLab <-icInfo$cognitive.domain[icInfo$signal==1]
  
  DATA = data.frame(images = list.files("ic_pngs_noBG",full.names = T), stringsAsFactors = F) %>% 
    mutate(
      names = gsub("[a-zA-Z]|[[:punct:]]","", images),
      values = sample(0:100, size=nrow(.), replace = T)
    )
  head(DATA)
  
  # Assign group in the original data frame
  DATA$Group = netLab %>% as.factor()
  
  nGroups = length(unique(DATA$Group))
  
  #### arrange circle plot by group
  oDATA = DATA %>% arrange(Group)
  
  #### trying to avoid track text changing direction within a network
  oDATA$NewOrder = rep(NA, nrow(oDATA))
  oDATA$NewOrder[oDATA$Group=="Attn"]=7
  oDATA$NewOrder[oDATA$Group=="Control"]=1
  oDATA$NewOrder[oDATA$Group=="SalVenAtt"]=3
  oDATA$NewOrder[oDATA$Group=="CB"]=4
  oDATA$NewOrder[oDATA$Group=="Default"]=2
  oDATA$NewOrder[oDATA$Group=="SomMot"]=5
  oDATA$NewOrder[oDATA$Group=="Vis"]=6
  
  oDATA = oDATA %>% arrange(NewOrder)
  
   # Get colours, turn the vector into characters, and then use the factor numbers of DATA$NewOrder to assign colour.
  
  #cColors = wesanderson::wes_palette("FantasticFox1", nGroups, type = 'continuous')
  cColors = viridis::plasma(nGroups)
  
  oDATA$gColor = cColors[1:nGroups] %>% as.character %>% 
    .[as.numeric(oDATA$NewOrder)]
  
  #### Define small gap between ICs and big gap between networks
  # sectorGap[1] is the gap after the first sector
  # sectorGap[end] should be bigGap
  
  smallGap = 1
  bigGap = 7
  sectorGap = rep(bigGap, nrow(DATA))
  
  for (ig in 1:nrow(oDATA)-1){
    sectorGap[ig] = ifelse(oDATA$NewOrder[ig+1]==oDATA$NewOrder[ig], smallGap, bigGap)
  }
  
  ############################################
  # Create the weighted adjacency matrix:
  
  # Construct the "from" edge:
  ic.from =   substr(data$EdgeName,1,2)
  ic.from = ifelse(substr(ic.from,2,2)=='-',paste0('0',substr(ic.from,1,1)),ic.from)
  # visually check it looks good:
  cbind(data$EdgeName,ic.from)
  
  # Construct the "to" edge
  ic.to = substr(data$EdgeName,nchar(data$EdgeName)-1,nchar(data$EdgeName))
  ic.to = ifelse(substr(ic.to,1,1)=='-',paste0('0',substr(ic.to,2,2)),ic.to)
  # visually check it looks good:
  cbind(data$EdgeName,ic.to)
  
  # modify the THIRD variable to create chords for different variables:
  PAIRS = data.frame(ic.from,ic.to,data[,varname])
  
  #### Index specific pairs by a threshold
  
  neg_ix = data[,alphaname] < alpha & PAIRS[,3]<0
  pos_ix = data[,alphaname] < alpha & PAIRS[,3]>0
  
  cols = rep(trcol,nrow(PAIRS))
  cols[neg_ix]=negcol
  cols[pos_ix]=poscol
  
  ## Set plot background to white
  par(bg = 'white')
  
  ## 
  circos.par(gap.after = sectorGap)
  
  ##plot updated Chord diagram with new colors
  myPlot = chordDiagram(PAIRS, annotationTrack = c("grid"), 
                        annotationTrackHeight=c(0.05, 0.01), 
                        preAllocateTracks = list(track.height = 0.32), 
                        order=oDATA$names, grid.col = setNames(oDATA$gColor, oDATA$names), col=cols)
  title(title,col.main="black",cex.main=1.5)
  
  as_radians = function(x) x*pi/180
  
  ## Loop through images and plot brain maps
  u=0
  for(si in get.all.sector.index()){
    xplot=get.cell.meta.data("xplot",si)
    u=u+1
    
    # Small workaround because coordinate 0 should be 360
    if(xplot[1] == 0) xplot[1] = 360
    
    
    x=.86*cos(as_radians((xplot[2]+xplot[1])/2))
    y=.86*sin(as_radians((xplot[2]+xplot[1])/2))
    
    oDATA$images[grep(si, oDATA$images)] %>% 
      readPNG() %>% 
      rasterImage(x-0.09, y-0.09, x+0.09, y+0.09)
  }
  
  #### set labels to be black
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    
    circos.text(mean(xlim), ylim[1], sector.name, facing = "inside", 
                niceFacing = TRUE, adj = c(0.5, 0), col= "black")
  }, bg.border = NA)
  
  #    myPlot
}
```

#### Load results from two sets of 200 seeds

```r
results.all=NULL
for (i in 1:200) {
  load(paste0('./Results/noImputation/ic30_pc85_glm_gam_drtmle_seed',i,'.RData'))
  results.all = rbind(results.all,results.df)
}
rm(results.df)

# second set of seeds:
#Ben doesn't know why seed 365 failed
results.all2=NULL
for (i in 201:400) {
  load(paste0('./Results/noImputation/ic30_pc85_glm_gam_drtmle_seed',i,'.RData'))
  results.all2 = rbind(results.all2,results.df)
}
rm(results.df)
```

#### Compare results from two sets of 200 seeds

```r
results.ave = results.all%>%
  group_by(EdgeName)%>%
  summarize(mean.ASD.naive=mean(mean.ASD.naive), 
            mean.TD.naive=mean(mean.TD.naive),
            mean.diff.naive=mean(mean.diff.naive),
            z.stat.ASD.naive=mean(z.stat.ASD.naive),
            z.stat.TD.naive=mean(z.stat.TD.naive),
            z.stat.diff.naive=mean(z.stat.diff.naive),
            mean.ASD.SL=mean(mean.ASD.SL),
            mean.TD.SL=mean(mean.TD.SL),
            mean.diff.SL=mean(mean.diff.SL),
            z.stat.ASD.SL=mean(z.stat.ASD.SL),
            z.stat.TD.SL=mean(z.stat.TD.SL),
            z.stat.diff.SL=mean(z.stat.diff.SL),
            SL.naive.mean.diff=mean(abs(mean.diff.SL)-abs(mean.diff.naive)))

results.ave2 = results.all2%>%
  group_by(EdgeName)%>%
  summarize(mean.ASD.naive=mean(mean.ASD.naive), 
            mean.TD.naive=mean(mean.TD.naive),
            mean.diff.naive=mean(mean.diff.naive),
            z.stat.ASD.naive=mean(z.stat.ASD.naive),
            z.stat.TD.naive=mean(z.stat.TD.naive),
            z.stat.diff.naive=mean(z.stat.diff.naive),
            mean.ASD.SL=mean(mean.ASD.SL),
            mean.TD.SL=mean(mean.TD.SL),
            mean.diff.SL=mean(mean.diff.SL),
            z.stat.ASD.SL=mean(z.stat.ASD.SL),
            z.stat.TD.SL=mean(z.stat.TD.SL),
            z.stat.diff.SL=mean(z.stat.diff.SL),
            SL.naive.mean.diff=mean(abs(mean.diff.SL)-abs(mean.diff.naive)))

# compare two estimates:
plot(results.ave$z.stat.diff.SL~results.ave2$z.stat.diff.SL)
```

<img src="QC_Impacts_DRTMLE_Table_CirclePlotsWhiteBG//average-two-sets-1.png" width="672" />

```r
head(cbind(results.ave$z.stat.diff.SL,results.ave2$z.stat.diff.SL))
```

```
##            [,1]       [,2]
## [1,]  1.5813133  1.5659236
## [2,] -2.4105445 -2.4171422
## [3,] -0.3030594 -0.3970811
## [4,] -0.9387843 -0.8936481
## [5,]  1.4080784  1.2850823
## [6,] -0.1019924 -0.1256482
```

```r
# this result appears in the manuscript:
seedCor = cor(results.ave$z.stat.diff.SL,results.ave2$z.stat.diff.SL)
```

The correlation between the average DRTMLE ASD-TD Z-statistics across the 153 edges was **0.99877296268611**. 


#### Calculate FDR adjusted p values

```r
results.ave$p.SL = 2*(1-pnorm(abs(results.ave$z.stat.diff.SL)))
results.ave$p.SL.fdr = p.adjust(results.ave$p.SL,method='BH')

results.ave2$p.SL = 2*(1-pnorm(abs(results.ave2$z.stat.diff.SL)))
results.ave2$p.SL.fdr = p.adjust(results.ave2$p.SL,method='BH')
```

#### Compare the selected edges at FDR=0.05 in the two sets of seeds:

```r
results.ave[results.ave$p.SL.fdr<0.05,c('EdgeName')]
```

```
## # A tibble: 6 × 1
##   EdgeName   
##   <chr>      
## 1 r.ic13.ic26
## 2 r.ic14.ic19
## 3 r.ic14.ic21
## 4 r.ic19.ic26
## 5 r.ic2.ic27 
## 6 r.ic4.ic17
```

```r
results.ave2[results.ave2$p.SL.fdr<0.05,c('EdgeName')]
```

```
## # A tibble: 6 × 1
##   EdgeName   
##   <chr>      
## 1 r.ic13.ic26
## 2 r.ic14.ic19
## 3 r.ic14.ic21
## 4 r.ic19.ic26
## 5 r.ic2.ic27 
## 6 r.ic4.ic17
```

6 were selected using DRTMLE and a false discovery rate FDR=.05 in the first set of seeds and the same 6 in the second set.


#### Compare the selected edges at FDR=0.20 in the two sets of seeds:

```r
results.ave[results.ave$p.SL.fdr<0.2,c('EdgeName')]
```

```
## # A tibble: 25 × 1
##    EdgeName   
##    <chr>      
##  1 r.ic1.ic14 
##  2 r.ic1.ic21 
##  3 r.ic1.ic24 
##  4 r.ic13.ic22
##  5 r.ic13.ic24
##  6 r.ic13.ic25
##  7 r.ic13.ic26
##  8 r.ic14.ic19
##  9 r.ic14.ic21
## 10 r.ic14.ic24
## # … with 15 more rows
```

```r
results.ave2[results.ave2$p.SL.fdr<0.2,c('EdgeName')]
```

```
## # A tibble: 25 × 1
##    EdgeName   
##    <chr>      
##  1 r.ic1.ic14 
##  2 r.ic1.ic21 
##  3 r.ic1.ic24 
##  4 r.ic13.ic22
##  5 r.ic13.ic24
##  6 r.ic13.ic25
##  7 r.ic13.ic26
##  8 r.ic14.ic19
##  9 r.ic14.ic21
## 10 r.ic14.ic24
## # … with 15 more rows
```

25 were selected using DRTMLE and a false discovery rate FDR=.20 in the first set of seeds and the same 25 in the second set.

#### Average results from all seeds

```r
#SL = superlearner

results.ave = rbind(results.all, results.all2) %>%
  group_by(EdgeName) %>%
  summarize(mean.ASD.naive=round(mean(mean.ASD.naive), digits = 3),
            mean.TD.naive=round(mean(mean.TD.naive), digits = 3),
            mean.diff.naive=mean(mean.diff.naive), 
            z.stat.ASD.naive=mean(z.stat.ASD.naive),
            z.stat.TD.naive=mean(z.stat.TD.naive),
            z.stat.diff.naive=round(mean(z.stat.diff.naive), digits = 3),
            mean.ASD.SL=round(mean(mean.ASD.SL), digits = 3), 
            mean.TD.SL=round(mean(mean.TD.SL), digits = 3),
            mean.diff.SL=round(mean(mean.diff.SL), digits = 3), 
            z.stat.ASD.SL=mean(z.stat.ASD.SL),
            z.stat.TD.SL=mean(z.stat.TD.SL), 
            z.stat.diff.SL=round(mean(z.stat.diff.SL), digits = 2),
            SL.naive.z.stat.diff=mean(abs(z.stat.diff.SL)-abs(z.stat.diff.naive)))
```



#### Add network names

```r
icInfo <- read_excel("componentLabels_pca85_ica30.xlsx")

results.ave <- results.ave %>% 
  mutate(EdgeName = gsub("r.|ic", "", EdgeName)) %>% 
  mutate(EdgeName = gsub("[.]", "-", EdgeName)) %>% 
  mutate(ic.from = gsub("[-]", "", substr(EdgeName, 1, 2))) %>% 
  mutate(ic.to = gsub("[-]", "", substr(EdgeName, nchar(EdgeName)-1,nchar(EdgeName)))) %>% 
  mutate(net1 = icInfo$cognitive.domain[as.numeric(ic.from)]) %>% 
  mutate(net2 = icInfo$cognitive.domain[as.numeric(ic.to)])
```

#### Calculate uncorrected p values for results averaged across 400 seeds for both methods and plot histograms

```r
#calculate p values
results.ave$p.naive = 2*(1-pnorm(abs(results.ave$z.stat.diff.naive)))
results.ave$p.SL = 2*(1-pnorm(abs(results.ave$z.stat.diff.SL)))

pVariables = c("EdgeName", "p.naive", "p.SL")

pvalues <- reshape2::melt(results.ave[, pVariables],
                         id.vars="EdgeName",
                         variable.name = "method",
                         value.name = "p")

method.labs <- c("naive", "DRTMLE")
names(method.labs) <- c("p.naive", "p.SL")

ggplot(pvalues, aes(x=p))+
  geom_histogram(position = "identity", alpha=0.5, inherit.aes=TRUE, binwidth = .05)+
  scale_x_continuous(expand = c(0, 0))+  
  scale_y_continuous(expand = c(0, 0))+ 
  facet_grid(method~., labeller = labeller(method = method.labs))
```

<img src="QC_Impacts_DRTMLE_Table_CirclePlotsWhiteBG//plot-uncorrected-p-values-1.png" width="672" />

#### FDR-adjusted p values

```r
results.ave$p.naive.fdr = round(p.adjust(results.ave$p.naive,method='BH'), digits = 4)
results.ave$p.SL.fdr = round(p.adjust(results.ave$p.SL,method='BH'), digits = 4)
```

#### Determine edges showing a significant group difference FDR=.05 averaged across all 400 seeds

```r
fdr05.edges.SL = which(results.ave$p.SL.fdr<.05) 
fdr05.edges.naive = which(results.ave$p.naive.fdr<.05) 

all(fdr05.edges.naive %in% fdr05.edges.SL)
```

```
## [1] TRUE
```

All edges indicated by naive approach are also indicated by DRTMLE at the FDR=0.05 level

#### Summary of network pairs showing a significant group difference using DRTMLE, FDR=.05

```r
tabfdr05<- tableby(net1 ~ net2,
                 data=filter(results.ave , p.SL.fdr<.05))
summary(tabfdr05,  
        title='Network pairs showing a group difference using DRTMLE FDR=.05',digits=1, digits.p=4,digits.pct=1, numeric.simplify=TRUE, total=TRUE, test=FALSE)
```



Table: Network pairs showing a group difference using DRTMLE FDR=.05

|                          | Attn (N=1) | CB (N=3)  | SomMot (N=1) | Vis (N=1)  | Total (N=6) |
|:-------------------------|:----------:|:---------:|:------------:|:----------:|:-----------:|
|**net2**                  |            |           |              |            |             |
|&nbsp;&nbsp;&nbsp;Attn    |  0 (0.0%)  | 1 (33.3%) |   0 (0.0%)   |  0 (0.0%)  |  1 (16.7%)  |
|&nbsp;&nbsp;&nbsp;Control |  0 (0.0%)  | 1 (33.3%) |   0 (0.0%)   | 1 (100.0%) |  2 (33.3%)  |
|&nbsp;&nbsp;&nbsp;Default | 1 (100.0%) | 1 (33.3%) |  1 (100.0%)  |  0 (0.0%)  |  3 (50.0%)  |


#### Remove background from .pngs of IC maps

```r
#create directory to save pngs with transparent background
dir.create("./ic_pngs_noBG")
```

```
## Warning in dir.create("./ic_pngs_noBG"): './ic_pngs_noBG' already exists
```

```r
DATA = data.frame(images = list.files("ic_pngs",full.names = T), stringsAsFactors = F) %>% 
    mutate(
      names = gsub("[a-zA-Z]|[[:punct:]]","", images),
      values = sample(0:100, size=nrow(.), replace = T)
    )
  head(DATA)
```

```
##             images names values
## 1 ic_pngs/ic01.png    01     44
## 2 ic_pngs/ic02.png    02     54
## 3 ic_pngs/ic04.png    04     65
## 4 ic_pngs/ic08.png    08     57
## 5 ic_pngs/ic13.png    13     16
## 6 ic_pngs/ic14.png    14     22
```

```r
for(iimg in seq(1:nrow(DATA))) {
 img_data <- magick::image_read(DATA$images[iimg])
 new_img <- magick::image_transparent(img_data, 'black')
 fname <- sub('ic_pngs', 'ic_pngs_noBG', DATA$images[iimg])
 magick::image_write(new_img, fname)
}
```


#### Naive ASD-TD Z-Statistic FDR=0.05

```r
png("Application_Figures/Naive_zstat_groupdifference_fdr05_noImpute.png",width=5,height=5,units="in",res=200)

myChordDiagram(data=results.ave,varname='z.stat.diff.naive',alpha=0.05,alphaname='p.naive.fdr', title='')
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
circos.clear()
```

#### DRTMLE ASD-TD Z-Statistic FDR=0.05

```r
png("Application_Figures/DRTMLE_zstat_groupdifference_fdr05_noImpute.png",width=5,height=5,units="in",res=200)
myChordDiagram(data=results.ave,varname='z.stat.diff.SL',alpha=0.05,alphaname='p.SL.fdr',title='')
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
circos.clear()
```

#### Print a and b together for report
<img src="QC_Impacts_DRTMLE_Table_CirclePlotsWhiteBG//circle-fdr-05-1.png" width="672" /><img src="QC_Impacts_DRTMLE_Table_CirclePlotsWhiteBG//circle-fdr-05-2.png" width="672" />

**Figure S3. The DRTMLE deconfounded group difference revealed more extensive differences than the naïve approach.** Z-statistics for autism spectrum disorder (ASD) versus typically developing (TD) using a) the na\"ive test and b) using DRTMLE. Connections are thresholded using a false discovery rate (FDR) of  0.05. Blue lines indicate ASD>TD (0 in naive, 3 in DRTMLE). Red lines indicate ASD<TD (1 in naive,  3 in DRTMLE). Brain regions contributing to each independent component are illustrated and components are grouped by functional assignment. Navy nodes: control. Blue violet: default mode. Purple: salience/ventral attention. Magenta: pontomedullary/cerebellar. Coral: somatomotor. Orange: visual. Yellow: dorsal attention.

#### Determine edges showing a significant group difference FDR=.20

```r
fdr20.edges.SL = which(results.ave$p.SL.fdr<=.20) 
fdr20.edges.naive = which(results.ave$p.naive.fdr<=.20) 

all(fdr20.edges.naive %in% fdr20.edges.SL)
```

```
## [1] TRUE
```
All edges indicated by naive approach are also indicated by DRTMLE at the FDR=0.2 level

#### Summary of network pairs showing a significant group difference using the naive approach, FDR=.2

```r
tabfdr2n<- tableby(net1 ~ net2,
                 data=filter(results.ave , p.naive.fdr<=.2))
summary(tabfdr2n,  
        title='Network pairs showing a group difference using the naive approach, FDR=.2',digits=1, digits.p=4,digits.pct=1, numeric.simplify=TRUE, total=TRUE, test=FALSE)
```



Table: Network pairs showing a group difference using the naive approach, FDR=.2

|                          | Attn (N=2) | CB (N=3)  | Default (N=1) | SomMot (N=1) | Vis (N=1)  | Total (N=8) |
|:-------------------------|:----------:|:---------:|:-------------:|:------------:|:----------:|:-----------:|
|**net2**                  |            |           |               |              |            |             |
|&nbsp;&nbsp;&nbsp;Attn    |  0 (0.0%)  | 1 (33.3%) |   0 (0.0%)    |   0 (0.0%)   |  0 (0.0%)  |  1 (12.5%)  |
|&nbsp;&nbsp;&nbsp;Control |  0 (0.0%)  | 1 (33.3%) |  1 (100.0%)   |   0 (0.0%)   | 1 (100.0%) |  3 (37.5%)  |
|&nbsp;&nbsp;&nbsp;Default | 2 (100.0%) | 1 (33.3%) |   0 (0.0%)    |  1 (100.0%)  |  0 (0.0%)  |  4 (50.0%)  |

#### Summary of network pairs showing a significant group difference using DRTMLE, FDR=.2

```r
tabfdr2<- tableby(net1 ~ net2,
                 data=filter(results.ave , p.SL.fdr<=.2))
summary(tabfdr2,  
        title='Network pairs showing a group difference using DRTMLE, FDR=.2',digits=1, digits.p=4,digits.pct=1, numeric.simplify=TRUE, total=TRUE, test=FALSE)
```



Table: Network pairs showing a group difference using DRTMLE, FDR=.2

|                            | Attn (N=4) | CB (N=5)  | Control (N=1) | Default (N=4) | SomMot (N=9) | Vis (N=2) | Total (N=25) |
|:---------------------------|:----------:|:---------:|:-------------:|:-------------:|:------------:|:---------:|:------------:|
|**net2**                    |            |           |               |               |              |           |              |
|&nbsp;&nbsp;&nbsp;Attn      |  0 (0.0%)  | 1 (20.0%) |   0 (0.0%)    |   0 (0.0%)    |   0 (0.0%)   | 0 (0.0%)  |   1 (4.0%)   |
|&nbsp;&nbsp;&nbsp;CB        |  0 (0.0%)  | 0 (0.0%)  |   0 (0.0%)    |   0 (0.0%)    |  1 (11.1%)   | 0 (0.0%)  |   1 (4.0%)   |
|&nbsp;&nbsp;&nbsp;Control   | 2 (50.0%)  | 3 (60.0%) |   0 (0.0%)    |   2 (50.0%)   |  4 (44.4%)   | 1 (50.0%) |  12 (48.0%)  |
|&nbsp;&nbsp;&nbsp;Default   | 2 (50.0%)  | 1 (20.0%) |   0 (0.0%)    |   1 (25.0%)   |  2 (22.2%)   | 0 (0.0%)  |  6 (24.0%)   |
|&nbsp;&nbsp;&nbsp;SalVenAtt |  0 (0.0%)  | 0 (0.0%)  |  1 (100.0%)   |   1 (25.0%)   |  2 (22.2%)   | 1 (50.0%) |  5 (20.0%)   |

#### Create table displaying group Z-statistics from DRTMLE for edges showing a group difference at FDR=.20

```r
# group means for edges showing a group difference using DRTMLE:
results.ave$NetworkPair <- stringr::str_c(results.ave$net1, '-', results.ave$net2)

#results.sigEdges = results.ave[fdr20.edges.SL,c(1, 23, 2:3, 7:9,13:14, 21:22)]
results.sigEdges = results.ave[fdr20.edges.SL, c("EdgeName", "NetworkPair",
                                                 "mean.ASD.naive", "mean.TD.naive", "z.stat.diff.naive",
                                                 "mean.ASD.SL", "mean.TD.SL", "z.stat.diff.SL",
                                                 "p.naive.fdr", "p.SL.fdr")]

results.sigEdges <- results.sigEdges[order(results.sigEdges$p.naive.fdr), ]

knitr::kable(results.sigEdges, format="markdown", 
             col.names = c("Edge", "Networks", 
                           "naive ASD mean", "naive TD mean",  "naive ASD-TD Z",
                           "DRTMLE ASD mean", "DRTMLE TD mean", "DRTMLE ASD-TD Z",
                           "naive FDR-adjusted p", 
                           "DRTMLE FDR-adjusted p"))
```



|Edge  |Networks          | naive ASD mean| naive TD mean| naive ASD-TD Z| DRTMLE ASD mean| DRTMLE TD mean| DRTMLE ASD-TD Z| naive FDR-adjusted p| DRTMLE FDR-adjusted p|
|:-----|:-----------------|--------------:|-------------:|--------------:|---------------:|--------------:|---------------:|--------------------:|---------------------:|
|2-27  |Vis-Control       |         -0.038|        -0.024|         -4.088|          -0.040|         -0.024|           -4.67|               0.0067|                0.0005|
|14-19 |CB-Attn           |          0.068|         0.078|         -3.351|           0.068|          0.078|           -3.18|               0.0616|                0.0451|
|13-26 |SomMot-Default    |         -0.031|        -0.039|          2.805|          -0.030|         -0.039|            3.12|               0.1283|                0.0461|
|14-21 |CB-Control        |          0.044|         0.053|         -2.862|           0.041|          0.052|           -3.36|               0.1283|                0.0451|
|19-25 |Attn-Default      |         -0.085|        -0.094|          2.867|          -0.086|         -0.094|            2.45|               0.1283|                0.1681|
|19-26 |Attn-Default      |          0.037|         0.028|          2.997|           0.037|          0.027|            3.29|               0.1283|                0.0451|
|4-17  |CB-Default        |         -0.067|        -0.076|          2.671|          -0.067|         -0.076|            3.19|               0.1653|                0.0451|
|17-24 |Default-Control   |         -0.005|         0.003|         -2.561|          -0.006|          0.003|           -2.94|               0.1996|                0.0628|
|1-14  |SomMot-CB         |         -0.002|         0.007|         -2.225|          -0.003|          0.006|           -2.41|               0.3990|                0.1685|
|14-24 |CB-Control        |          0.015|         0.021|         -2.099|           0.014|          0.021|           -2.30|               0.4567|                0.1685|
|21-30 |Control-SalVenAtt |          0.063|         0.069|         -2.115|           0.063|          0.069|           -2.33|               0.4567|                0.1685|
|13-25 |SomMot-Default    |         -0.097|        -0.102|          1.869|          -0.094|         -0.102|            2.66|               0.4696|                0.1087|
|15-30 |SomMot-SalVenAtt  |          0.113|         0.117|         -1.830|           0.112|          0.117|           -2.29|               0.4696|                0.1685|
|19-24 |Attn-Control      |         -0.160|        -0.165|          2.037|          -0.159|         -0.165|            2.25|               0.4696|                0.1713|
|4-21  |CB-Control        |         -0.041|        -0.034|         -1.918|          -0.043|         -0.034|           -2.17|               0.4696|                0.1836|
|8-22  |Vis-SalVenAtt     |          0.019|         0.024|         -1.782|           0.017|          0.025|           -2.48|               0.4696|                0.1675|
|26-30 |Default-SalVenAtt |         -0.088|        -0.082|         -1.720|          -0.090|         -0.082|           -2.29|               0.5027|                0.1685|
|1-24  |SomMot-Control    |         -0.027|        -0.032|          1.695|          -0.023|         -0.033|            2.74|               0.5104|                0.0940|
|1-21  |SomMot-Control    |          0.003|        -0.002|          1.668|           0.008|         -0.002|            3.02|               0.5208|                0.0552|
|13-22 |SomMot-SalVenAtt  |         -0.030|        -0.035|          1.564|          -0.029|         -0.035|            2.35|               0.5223|                0.1685|
|17-27 |Default-Control   |          0.034|         0.028|          1.557|           0.036|          0.028|            2.17|               0.5223|                0.1836|
|13-24 |SomMot-Control    |         -0.091|        -0.095|          1.117|          -0.088|         -0.095|            2.24|               0.6793|                0.1713|
|15-21 |SomMot-Control    |          0.039|         0.042|         -1.098|           0.038|          0.044|           -2.23|               0.6793|                0.1713|
|19-21 |Attn-Control      |         -0.009|        -0.013|          1.086|          -0.006|         -0.015|            2.75|               0.6793|                0.0940|
|17-25 |Default-Default   |          0.033|         0.032|          0.425|           0.041|          0.031|            2.31|               0.9164|                0.1685|

DRTMLE indicates a significant deconfounded group difference for the first **6** edges at the FDR=.05 level and an additional **19** edges at the FDR=.20 level.

Mean functional connectivity group estimates further from zero reflect stronger functional connectivity regardless of sign; positive scores reflect positive partial correlations, or more integrated intrinsic activity between independent components (ICs). Negative scores reflect negative partial correlations, or more segregated intrinsic activity between ICs.

#### Create version for paper

```r
#tab <- knitr::kable(results.sigEdges, format="latex", linesep = '')
#print(tab, type="latex")

knitr::kable(results.sigEdges, "latex", linesep = "", vline = "",
             col.names = c("Edge", "Network Pair", 
                           "naive ASD mean", "naive TD mean", "naive ASD-TD Z",
                           "DRTMLE ASD mean", "DRTMLE TD mean",  "DRTMLE ASD-TD Z",
                           "naive FDR-adjusted p", 
                           "DRTMLE FDR-adjusted p")) %>%
  #kableExtra::column_spec(9:10, width = "1 cm") %>% 
  kableExtra::kable_styling(latex_options = c("scale_down")) %>% 
  kableExtra::as_image(width = 7, file = "Application_Figures/tableS3.pdf")
```

```
## Warning in kableExtra::as_image(., width = 7, file = "Application_Figures/
## tableS3.pdf"): You need to install magick in order to use width/height in
## as_image.
```

#### Naive ASD-TD Z-Statistic FDR=0.20

```r
png("Application_Figures/Naive_zstat_groupdifference_fdr20_noImpute.png",width=5,height=5,units="in",res=200)
myChordDiagram(data=results.ave,varname='z.stat.diff.naive',alpha=0.20,alphaname='p.naive.fdr', title='')
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
circos.clear()
```

#### DRTMLE ASD-TD Z-Statistic FDR=0.20

```r
png("Application_Figures/DRTMLE_zstat_groupdifference_fdr20_noImpute.png",width=5,height=5,units="in",res=200)
#myChordDiagram(data=results.ave,varname='z.stat.diff.SL',alpha=0.20,alphaname='p.SL.fdr',title='B. DRTMLE Z-Statistic, ASD-TD')
myChordDiagram(data=results.ave,varname='z.stat.diff.SL',alpha=0.20,alphaname='p.SL.fdr',title='')
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
circos.clear()
```

#### Print a and b together for .html
<img src="QC_Impacts_DRTMLE_Table_CirclePlotsWhiteBG//circle-fdr-20-1.png" width="672" /><img src="QC_Impacts_DRTMLE_Table_CirclePlotsWhiteBG//circle-fdr-20-2.png" width="672" />


**Figure 5. The DRTMLE deconfounded group difference revealed more extensive differences than the naïve approach.** Z-statistics for autism spectrum disorder (ASD) versus typically developing (TD) using a) the na\"ive test and b) using DRTMLE. Connections are thresholded using a false discovery rate (FDR) of  0.20. Blue lines indicate ASD>TD (4 in naive, 13 in DRTMLE). Red lines indicate ASD<TD (4 in naive,  12 in DRTMLE). Brain regions contributing to each independent component are illustrated and components are grouped by functional assignment. Navy nodes: control. Blue violet: default mode. Purple: salience/ventral attention. Magenta: pontomedullary/cerebellar. Coral: somatomotor. Orange: visual. Yellow: dorsal attention.


#### [Daniel Lidstone's](https://github.com/lidstone) analysis of edges indicated by DRTMLE to show a group difference at FDR=0.20


```r
SM_ix = grep('SomMot',c(icInfo$cognitive.domain))
Vis_ix = grep('Vis',c(icInfo$cognitive.domain))
CB_ix = grep('CB',c(icInfo$cognitive.domain))
Con_ix = grep('Control',c(icInfo$cognitive.domain))
Def_ix = grep('Default',c(icInfo$cognitive.domain))
Sal_ix = grep('SalVenAtt',c(icInfo$cognitive.domain))
DorAtt_ix = grep('Attn',c(icInfo$cognitive.domain))

naivefdr = results.ave[abs(results.ave$p.naive.fdr)<0.20,c('EdgeName','p.naive.fdr','z.stat.diff.naive')]
slfdr = results.ave[abs(results.ave$p.SL.fdr)<0.20,c('EdgeName','p.SL.fdr','z.stat.diff.SL')]

naivefdr_edgeNum = as.numeric(unlist((stringr::str_extract_all(naivefdr$EdgeName, "\\d+"))))
slfdr_edgeNum = as.numeric(unlist((stringr::str_extract_all(slfdr$EdgeName, "\\d+"))))

naive_SM = length((which(naivefdr_edgeNum%in%SM_ix)))
naive_Vis = length((which(naivefdr_edgeNum%in%Vis_ix)))
naive_CB = length((which(naivefdr_edgeNum%in%CB_ix)))
naive_Con = length((which(naivefdr_edgeNum%in%Con_ix)))
naive_Def = length((which(naivefdr_edgeNum%in%Def_ix)))
naive_Sal = length((which(naivefdr_edgeNum%in%Sal_ix)))
naive_DorAtt = length((which(naivefdr_edgeNum%in%DorAtt_ix)))


sl_SM = length((which(slfdr_edgeNum%in%SM_ix)))
sl_Vis = length((which(slfdr_edgeNum%in%Vis_ix)))
sl_CB = length((which(slfdr_edgeNum%in%CB_ix)))
sl_Con = length((which(slfdr_edgeNum%in%Con_ix)))
sl_Def = length((which(slfdr_edgeNum%in%Def_ix)))
sl_Sal = length((which(slfdr_edgeNum%in%Sal_ix)))
sl_DorAtt = length((which(slfdr_edgeNum%in%DorAtt_ix)))


## Bar Plot
df <- data.frame(Network = c("SMN","Visual","Cerebellum","Control","Default","Ventral Attn","Dorsal Attn"), Naive = c(naive_SM,naive_Vis,naive_CB,naive_Con,naive_Def,naive_Sal,naive_DorAtt), DRTMLE = c(sl_SM,sl_Vis,sl_CB,sl_Con,sl_Def,sl_Sal,sl_DorAtt))
data.m <- reshape2::melt(df, id.vars='Network')

p <- ggplot(data=data.m, aes(x=Network, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=value), vjust=-0.1, color="black",
            position = position_dodge(0.9), size=5)+
  scale_fill_brewer(palette="Paired")+
  theme_minimal()

p + labs(x="Network", y = "Number of Edges")+
  scale_fill_manual(values=c('black','lightgray'))+
  theme_classic()+theme(legend.title = element_blank(),legend.text = element_text(size=14),axis.text=element_text(size=12),
                        axis.title=element_text(size=14,face="bold"))+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

<img src="QC_Impacts_DRTMLE_Table_CirclePlotsWhiteBG//lidstone-network-pairs-1.png" width="672" />

