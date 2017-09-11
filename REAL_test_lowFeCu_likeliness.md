# Linear Model_STATS-TABLE - is low FeCu response like lowFe or lowCu?
Anna A. Hippmann  
September 6, 2017  



##Libraries used

```r
suppressPackageStartupMessages(library(lsmeans))
suppressPackageStartupMessages(library(phia))
suppressPackageStartupMessages(library(visreg))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
library(ggplot2)
library(knitr)
```
<a id="BackUP"></a>

  + [Growthrate (dd-1)](#Growthrate.dd)
  + [Growthrate (d)](#Growthrate.specific.d.1)
  + [Growthrate (percent of umax)](#Growthrate.Percent..u.umax.)
  + [Fv/Fm](#FvFm.old)
  + [Sigma](#Sig.old)
  + [PQ Size](#PQ_Siz.old)
  + [Chl a per Cell Vol (fg/fL)](#Chla.per.cell.vol.fg.fL)
  + [Chl a per Cell (pg/cell)](#Chla.per.cell.pg.cell)
  + [Cell Size (um)](#cell.size..um)
  + [Cell vol (fL)](#cell.volume.fl.cell)
  + [Cell Surface area (um^2)](#cell.SA.um2)
  + [Cell Surface area / Vol - SA/V (um^2/fL)](#cell.SA.Vol.ratio)
  + [FeDFB uptake per cell (zmol/cell*h)](#FeDFB.zmol.cell.1.h.)
  + [FeDFB uptake per cell surface area (zmol/um^2*h))](#FeDFB.zmol.um.2.h.)
  + [14C uptake - at 155 uEinstein (per Chla)](#X14C.per.Chla.at.155uE)
  + [14C uptake - alpha (per Chl a)](#X14C.per.Chla.alpha)
  + [14C uptake - ek (per Chl a)](#X14C.per.Chla.ek)
  + [14C uptake - pmax (per Chl a)](#X14C.per.Chla.pmax)


  + [AOX - Alternative Oxidase (percent)](#AOXactivity)
  + [Gross Production per Chl a (mol O2 / mol Chla * h)](#GrossPchla.mol.O2.mol.Chla..h.)
  + [Gross Production per cell (mol O2 / cell * h)](#GrossPcell..mol.O2.cell..h.)
  + [Respiration per Chl a (mol O2 / mol Chla * h)](#Resp.mol.O2.mol.Chla.h)
  + [Respiration per cell (mol O2 / cell * h)](#Resp.umol.O2.cell.h)
  + [Gross Production per cell vol (mol O2 / cell vol * h)](#Gross.P.CellVol..mol.O2.L.)
  + [Total Cellular Protein (pg/cell)](#pg.prot.cell)
  
  + [Combined Table for physiological parameters](#Comb_Table_Phys)
  
* [FRRF DATA](#FRRF)
  + [sigma (A^2/RC)](#sig.)  
  + [Fq'/Fv' @ growth irradiance (A.U.)](#F.q.F.v)  
  + [Fv'/Fm'@ growth irradiance (A.U.)](#F.v.F.m)
  + [Fq'/Fm'@ growth irradiance (A.U.)](#F.q.F.m)
  + [ETR @ growth irradiance (mol e / RC * sec)](#ETR)
  + [NPQ (nsv) @ growth irradiance](#NPQ)
  + [ETR PE curve - alpha](#ETR.alpha.JP)
  + [ETR PE curve - ek)](#ETR.ek.JP)
  + [ETR PE curve - pmax)](#ETR.pmax.JP)
  + [Conversion Factor](#Converse_corr)
  
   + [Combined Table for FRRF parameters](#Comb_Table_FRRF)
  
  
  

```r
#I used the original "ALL_PhysiologicalData_2015_04.txt" and added columns "Treatment" and "Merged" to it, then saved the ensuing file under "ALL_PhysiologicalData_2015_04_withMerged_Variable.txt", so I can use it for all my plots [2017_09: added last protein/cell values]

mydata <- read.delim("ALL_PhysiologicalData_2015_04_withMerged_Variable.txt", sep="\t", header=T)


FRRF <- read.delim("Input_Data/STATS/FRRF_modified_TO03lowFe_Ave.txt", sep="\t", header=T)

#FRRF <- read.delim("Input_Data/ALL_Phys_Barplots/Compiled_ALL_Raw_R_2016_07_FRRF.txt", sep="\t", header=T)
#mean.df <- read.delim("Physiological_Data_Mean_Values.txt", sep="\t", header=T)

#mean.new <- read.delim("Input_Data/ALL_Phys_Barplots/ALL_Phys_both_TO03_TO05_mean_stderror.txt", sep="\t", header=T)


#str(mydata)
treatment <- ( c("TO03_lowCu", "TO03_lowFe", "TO03_lowFeCu", "TO05_lowCu", "TO05_lowFe", "TO05_lowFeCu"))

treatment <- as.data.frame(treatment)
```

<a id="Growthrate.dd"></a>

#### Growthrate.dd lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_Growthrate.dd <- lm(data=mydata, Growthrate.dd.1~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_Growthrate.dd_lowCu <- testInteractions(lm_Growthrate.dd, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_Growthrate.dd_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high  0.24333  1   0.08882  6.8189  0.037803 *  
## TO 1005 : high  0.58000  1   0.50460 38.7409 4.866e-05 ***
## TO 1003 :  low -0.13000  1   0.02535  1.9463  0.182056    
## TO 1005 :  low  0.41667  1   0.26042 19.9936  0.001157 ** 
## Residuals               16   0.20840                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Growthrate.dd_lowCu <- phia_Growthrate.dd_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "03_Growthrate.dd",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_Growthrate.dd_lowCu), DF_here = phia_Growthrate.dd_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_Growthrate.dd_lowFe <- testInteractions(lm_Growthrate.dd, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_Growthrate.dd_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq      F    Pr(>F)    
## TO 1003 : high 0.79333  1   0.94407 72.481 9.805e-07 ***
## TO 1005 : high 0.62667  1   0.58907 45.226 1.463e-05 ***
## TO 1003 :  low 0.42000  1   0.26460 20.315 0.0003582 ***
## TO 1005 :  low 0.46333  1   0.32202 24.723 0.0002768 ***
## Residuals              16   0.20840                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Growthrate.dd_lowFe <- phia_Growthrate.dd_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "03_Growthrate.dd",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_Growthrate.dd_lowFe), DF_here = phia_Growthrate.dd_lowFe[5,2]) %>% 
     slice (3:4)  





#Making Growthrate.dd data.frame
kable(( Growthrate.dd <- bind_rows(phia_Growthrate.dd_lowCu, phia_Growthrate.dd_lowFe)), format = "markdown")
```



|      Value| Df| Sum of Sq|         F|    Pr(>F)|Analysis         |Data             |test_intrctn     |tested         | DF_here|
|----------:|--:|---------:|---------:|---------:|:----------------|:----------------|:----------------|:--------------|-------:|
| -0.1300000|  1| 0.0253500|  1.946257| 0.1820559|phia_interaction |03_Growthrate.dd |lowFeCu vs lowCu |TO 1003 :  low |      16|
|  0.4166667|  1| 0.2604167| 19.993602| 0.0011571|phia_interaction |03_Growthrate.dd |lowFeCu vs lowCu |TO 1005 :  low |      16|
|  0.4200000|  1| 0.2646000| 20.314779| 0.0003582|phia_interaction |03_Growthrate.dd |lowFeCu vs lowFe |TO 1003 :  low |      16|
|  0.4633333|  1| 0.3220167| 24.722969| 0.0002768|phia_interaction |03_Growthrate.dd |lowFeCu vs lowFe |TO 1005 :  low |      16|

<a id="Growthrate.specific.d.1"></a>

#### Growthrate.specific.d.1 lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_Growthrate.specific.d.1 <- lm(data=mydata, Growthrate.specific.d.1~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_Growthrate.specific.d.1_lowCu <- testInteractions(lm_Growthrate.specific.d.1, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_Growthrate.specific.d.1_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high  0.16667  1  0.041667  6.7659  0.038593 *  
## TO 1005 : high  0.40000  1  0.240000 38.9716 4.702e-05 ***
## TO 1003 :  low -0.08333  1  0.010417  1.6915  0.211828    
## TO 1005 :  low  0.28667  1  0.123267 20.0162  0.001151 ** 
## Residuals               16  0.098533                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Growthrate.specific.d.1_lowCu <- phia_Growthrate.specific.d.1_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "04_Growthrate.specific.d.1",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_Growthrate.specific.d.1_lowCu), DF_here = phia_Growthrate.specific.d.1_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_Growthrate.specific.d.1_lowFe <- testInteractions(lm_Growthrate.specific.d.1, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_Growthrate.specific.d.1_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq      F    Pr(>F)    
## TO 1003 : high 0.54667  1   0.44827 72.790 9.532e-07 ***
## TO 1005 : high 0.43667  1   0.28602 46.444 1.246e-05 ***
## TO 1003 :  low 0.29667  1   0.13202 21.437 0.0002781 ***
## TO 1005 :  low 0.32333  1   0.15682 25.464 0.0002385 ***
## Residuals              16   0.09853                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Growthrate.specific.d.1_lowFe <- phia_Growthrate.specific.d.1_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "04_Growthrate.specific.d.1",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_Growthrate.specific.d.1_lowFe), DF_here = phia_Growthrate.specific.d.1_lowFe[5,2]) %>% 
     slice (3:4)  





#Making Growthrate.specific.d.1 data.frame
kable(( Growthrate.specific.d.1 <- bind_rows(phia_Growthrate.specific.d.1_lowCu, phia_Growthrate.specific.d.1_lowFe)), format = "markdown")
```



|      Value| Df| Sum of Sq|         F|    Pr(>F)|Analysis         |Data                       |test_intrctn     |tested         | DF_here|
|----------:|--:|---------:|---------:|---------:|:----------------|:--------------------------|:----------------|:--------------|-------:|
| -0.0833333|  1| 0.0104167|  1.691475| 0.2118281|phia_interaction |04_Growthrate.specific.d.1 |lowFeCu vs lowCu |TO 1003 :  low |      16|
|  0.2866667|  1| 0.1232667| 20.016238| 0.0011510|phia_interaction |04_Growthrate.specific.d.1 |lowFeCu vs lowCu |TO 1005 :  low |      16|
|  0.2966667|  1| 0.1320167| 21.437077| 0.0002781|phia_interaction |04_Growthrate.specific.d.1 |lowFeCu vs lowFe |TO 1003 :  low |      16|
|  0.3233333|  1| 0.1568167| 25.464141| 0.0002385|phia_interaction |04_Growthrate.specific.d.1 |lowFeCu vs lowFe |TO 1005 :  low |      16|



<a id="Growthrate.Percent..u.umax."></a>

#### Growthrate.Percent..u.umax. lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_Growthrate.Percent..u.umax. <- lm(data=mydata, Growthrate.Percent..u.umax.~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_Growthrate.Percent..u.umax._lowCu <- testInteractions(lm_Growthrate.Percent..u.umax., fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_Growthrate.Percent..u.umax._lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                 Value Df Sum of Sq       F   Pr(>F)    
## TO 1003 : high 16.003  1    384.16  7.9604 0.024572 *  
## TO 1005 : high 31.867  1   1523.23 31.5635 0.000154 ***
## TO 1003 :  low -8.310  1    103.58  2.1464 0.162277    
## TO 1005 :  low 23.000  1    793.50 16.4425 0.002759 ** 
## Residuals             16    772.15                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Growthrate.Percent..u.umax._lowCu <- phia_Growthrate.Percent..u.umax._lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "05_Growthrate.Percent..u.umax.",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_Growthrate.Percent..u.umax._lowCu), DF_here = phia_Growthrate.Percent..u.umax._lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_Growthrate.Percent..u.umax._lowFe <- testInteractions(lm_Growthrate.Percent..u.umax., fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_Growthrate.Percent..u.umax._lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                 Value Df Sum of Sq      F    Pr(>F)    
## TO 1003 : high 52.187  1    4085.2 84.651 3.457e-07 ***
## TO 1005 : high 34.363  1    1771.3 36.703 4.978e-05 ***
## TO 1003 :  low 27.873  1    1165.4 24.148 0.0003114 ***
## TO 1005 :  low 25.497  1     975.1 20.206 0.0003672 ***
## Residuals             16     772.1                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Growthrate.Percent..u.umax._lowFe <- phia_Growthrate.Percent..u.umax._lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "05_Growthrate.Percent..u.umax.",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_Growthrate.Percent..u.umax._lowFe), DF_here = phia_Growthrate.Percent..u.umax._lowFe[5,2]) %>% 
     slice (3:4)  





#Making Growthrate.Percent..u.umax. data.frame
kable(( Growthrate.Percent..u.umax. <- bind_rows(phia_Growthrate.Percent..u.umax._lowCu, phia_Growthrate.Percent..u.umax._lowFe)), format = "markdown")
```



|    Value| Df| Sum of Sq|         F|    Pr(>F)|Analysis         |Data                           |test_intrctn     |tested         | DF_here|
|--------:|--:|---------:|---------:|---------:|:----------------|:------------------------------|:----------------|:--------------|-------:|
| -8.31000|  1|  103.5841|  2.146417| 0.1622768|phia_interaction |05_Growthrate.Percent..u.umax. |lowFeCu vs lowCu |TO 1003 :  low |      16|
| 23.00000|  1|  793.5000| 16.442498| 0.0027586|phia_interaction |05_Growthrate.Percent..u.umax. |lowFeCu vs lowCu |TO 1005 :  low |      16|
| 27.87333|  1| 1165.3841| 24.148487| 0.0003114|phia_interaction |05_Growthrate.Percent..u.umax. |lowFeCu vs lowFe |TO 1003 :  low |      16|
| 25.49667|  1|  975.1200| 20.205934| 0.0003672|phia_interaction |05_Growthrate.Percent..u.umax. |lowFeCu vs lowFe |TO 1005 :  low |      16|





<a id="FvFm.old"></a>

#### FvFm.old lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_FvFm.old <- lm(data=mydata, FvFm.old~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_FvFm.old_lowCu <- testInteractions(lm_FvFm.old, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_FvFm.old_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high 0.24167  1  0.070083 36.6077 8.885e-05 ***
## TO 1005 : high 0.11000  1  0.018150  9.4806   0.02291 *  
## TO 1003 :  low 0.01000  1  0.000150  0.0784   0.78337    
## TO 1005 :  low 0.10667  1  0.017067  8.9147   0.02291 *  
## Residuals              15  0.028717                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_FvFm.old_lowCu <- phia_FvFm.old_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "12_FvFm.old",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_FvFm.old_lowCu), DF_here = phia_FvFm.old_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_FvFm.old_lowFe <- testInteractions(lm_FvFm.old, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_FvFm.old_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high 0.305000  1  0.111630 58.3093 6.082e-06 ***
## TO 1005 : high 0.036667  1  0.002017  1.0534    0.6420    
## TO 1003 :  low 0.073333  1  0.008067  4.2136    0.1739    
## TO 1005 :  low 0.033333  1  0.001667  0.8706    0.6420    
## Residuals               15  0.028717                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_FvFm.old_lowFe <- phia_FvFm.old_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "12_FvFm.old",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_FvFm.old_lowFe), DF_here = phia_FvFm.old_lowFe[5,2]) %>% 
     slice (3:4)  





#Making FvFm.old data.frame
kable(( FvFm.old <- bind_rows(phia_FvFm.old_lowCu, phia_FvFm.old_lowFe)), format = "markdown")
```



|     Value| Df| Sum of Sq|         F|    Pr(>F)|Analysis         |Data        |test_intrctn     |tested         | DF_here|
|---------:|--:|---------:|---------:|---------:|:----------------|:-----------|:----------------|:--------------|-------:|
| 0.0100000|  1| 0.0001500| 0.0783517| 0.7833668|phia_interaction |12_FvFm.old |lowFeCu vs lowCu |TO 1003 :  low |      15|
| 0.1066667|  1| 0.0170667| 8.9146837| 0.0229111|phia_interaction |12_FvFm.old |lowFeCu vs lowCu |TO 1005 :  low |      15|
| 0.0733333|  1| 0.0080667| 4.2135810| 0.1739240|phia_interaction |12_FvFm.old |lowFeCu vs lowFe |TO 1003 :  low |      15|
| 0.0333333|  1| 0.0016667| 0.8705746| 0.6419878|phia_interaction |12_FvFm.old |lowFeCu vs lowFe |TO 1005 :  low |      15|




<a id="Sig.old"></a>

#### Sig.old lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_Sig.old <- lm(data=mydata, Sig.old~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_Sig.old_lowCu <- testInteractions(lm_Sig.old, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_Sig.old_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high  -22.000  1     580.8  0.6146 0.448252   
## TO 1005 : high -115.500  1   13340.3 14.1171 0.008203 **
## TO 1003 :  low  -32.667  1    1600.7  1.6939 0.435046   
## TO 1005 :  low -133.500  1   21386.7 22.6321 0.001865 **
## Residuals               12   11339.7                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Sig.old_lowCu <- phia_Sig.old_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "13_Sig.old",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_Sig.old_lowCu), DF_here = phia_Sig.old_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_Sig.old_lowFe <- testInteractions(lm_Sig.old, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_Sig.old_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high -180.00  1     38880 41.1441 9.993e-05 ***
## TO 1005 : high    7.00  1        49  0.0519         1    
## TO 1003 :  low -190.67  1     54531 57.7061 2.544e-05 ***
## TO 1005 :  low  -11.00  1       145  0.1537         1    
## Residuals              12     11340                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Sig.old_lowFe <- phia_Sig.old_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "13_Sig.old",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_Sig.old_lowFe), DF_here = phia_Sig.old_lowFe[5,2]) %>% 
     slice (3:4)  





#Making Sig.old data.frame
kable(( Sig.old <- bind_rows(phia_Sig.old_lowCu, phia_Sig.old_lowFe)), format = "markdown")
```



|      Value| Df| Sum of Sq|          F|    Pr(>F)|Analysis         |Data       |test_intrctn     |tested         | DF_here|
|----------:|--:|---------:|----------:|---------:|:----------------|:----------|:----------------|:--------------|-------:|
|  -32.66667|  1|  1600.667|  1.6938770| 0.4350461|phia_interaction |13_Sig.old |lowFeCu vs lowCu |TO 1003 :  low |      12|
| -133.50000|  1| 21386.700| 22.6320938| 0.0018650|phia_interaction |13_Sig.old |lowFeCu vs lowCu |TO 1005 :  low |      12|
| -190.66667|  1| 54530.667| 57.7061054| 0.0000254|phia_interaction |13_Sig.old |lowFeCu vs lowFe |TO 1003 :  low |      12|
|  -11.00000|  1|   145.200|  0.1536553| 1.0000000|phia_interaction |13_Sig.old |lowFeCu vs lowFe |TO 1005 :  low |      12|




<a id="PQ_Siz.old"></a>

#### PQ_Siz.old lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_PQ_Siz.old <- lm(data=mydata, PQ_Siz.old~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_PQ_Siz.old_lowCu <- testInteractions(lm_PQ_Siz.old, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_PQ_Siz.old_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq      F Pr(>F)
## TO 1003 : high 1.45000  1    2.5230 1.3315 0.8130
## TO 1005 : high 1.40000  1    1.9600 1.0343 0.8130
## TO 1003 :  low 0.33333  1    0.1667 0.0880 0.8130
## TO 1005 :  low 2.19667  1    5.7904 3.0558 0.4239
## Residuals              12   22.7389
```

```r
phia_PQ_Siz.old_lowCu <- phia_PQ_Siz.old_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "14_PQ_Siz.old",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_PQ_Siz.old_lowCu), DF_here = phia_PQ_Siz.old_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_PQ_Siz.old_lowFe <- testInteractions(lm_PQ_Siz.old, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_PQ_Siz.old_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq      F Pr(>F)
## TO 1003 : high -1.6500  1    3.2670 1.7241 0.6412
## TO 1005 : high  0.5500  1    0.3025 0.1596 0.6965
## TO 1003 :  low -2.7667  1   11.4817 6.0592 0.1198
## TO 1005 :  low  1.3467  1    2.1762 1.1485 0.6412
## Residuals              12   22.7389
```

```r
phia_PQ_Siz.old_lowFe <- phia_PQ_Siz.old_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "14_PQ_Siz.old",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_PQ_Siz.old_lowFe), DF_here = phia_PQ_Siz.old_lowFe[5,2]) %>% 
     slice (3:4)  





#Making PQ_Siz.old data.frame
kable(( PQ_Siz.old <- bind_rows(phia_PQ_Siz.old_lowCu, phia_PQ_Siz.old_lowFe)), format = "markdown")
```



|      Value| Df|  Sum of Sq|         F|    Pr(>F)|Analysis         |Data          |test_intrctn     |tested         | DF_here|
|----------:|--:|----------:|---------:|---------:|:----------------|:-------------|:----------------|:--------------|-------:|
|  0.3333333|  1|  0.1666667| 0.0879549| 0.8130291|phia_interaction |14_PQ_Siz.old |lowFeCu vs lowCu |TO 1003 :  low |      12|
|  2.1966667|  1|  5.7904133| 3.0557704| 0.4238525|phia_interaction |14_PQ_Siz.old |lowFeCu vs lowCu |TO 1005 :  low |      12|
| -2.7666667|  1| 11.4816667| 6.0592112| 0.1198142|phia_interaction |14_PQ_Siz.old |lowFeCu vs lowFe |TO 1003 :  low |      12|
|  1.3466667|  1|  2.1762133| 1.1484514| 0.6411668|phia_interaction |14_PQ_Siz.old |lowFeCu vs lowFe |TO 1005 :  low |      12|




<a id="Chla.per.cell.vol.fg.fL"></a>

#### Chla.per.cell.vol.fg.fL lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_Chla.per.cell.vol.fg.fL <- lm(data=mydata, Chla.per.cell.vol.fg.fL~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_Chla.per.cell.vol.fg.fL_lowCu <- testInteractions(lm_Chla.per.cell.vol.fg.fL, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_Chla.per.cell.vol.fg.fL_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high -0.2733  1    0.1121  0.2819  0.602739    
## TO 1005 : high  1.9557  1    5.7369 14.4321  0.004729 ** 
## TO 1003 :  low  4.3810  1   28.7897 72.4247 9.856e-07 ***
## TO 1005 :  low  0.8391  1    1.0561  2.6567  0.245276    
## Residuals              16    6.3602                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Chla.per.cell.vol.fg.fL_lowCu <- phia_Chla.per.cell.vol.fg.fL_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "11_Chla.per.cell.vol.fg.fL",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_Chla.per.cell.vol.fg.fL_lowCu), DF_here = phia_Chla.per.cell.vol.fg.fL_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_Chla.per.cell.vol.fg.fL_lowFe <- testInteractions(lm_Chla.per.cell.vol.fg.fL, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_Chla.per.cell.vol.fg.fL_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high -3.06667  1   14.1067 35.4873 8.035e-05 ***
## TO 1005 : high  0.48600  1    0.3543  0.8913   0.47665    
## TO 1003 :  low  1.58767  1    3.7810  9.5117   0.02134 *  
## TO 1005 :  low -0.63059  1    0.5965  1.5005   0.47665    
## Residuals               16    6.3602                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Chla.per.cell.vol.fg.fL_lowFe <- phia_Chla.per.cell.vol.fg.fL_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "11_Chla.per.cell.vol.fg.fL",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_Chla.per.cell.vol.fg.fL_lowFe), DF_here = phia_Chla.per.cell.vol.fg.fL_lowFe[5,2]) %>% 
     slice (3:4)  





#Making Chla.per.cell.vol.fg.fL data.frame
kable(( Chla.per.cell.vol.fg.fL <- bind_rows(phia_Chla.per.cell.vol.fg.fL_lowCu, phia_Chla.per.cell.vol.fg.fL_lowFe)), format = "markdown")
```



|      Value| Df|  Sum of Sq|         F|    Pr(>F)|Analysis         |Data                       |test_intrctn     |tested         | DF_here|
|----------:|--:|----------:|---------:|---------:|:----------------|:--------------------------|:----------------|:--------------|-------:|
|  4.3810000|  1| 28.7897415| 72.424712| 0.0000010|phia_interaction |11_Chla.per.cell.vol.fg.fL |lowFeCu vs lowCu |TO 1003 :  low |      16|
|  0.8390809|  1|  1.0560851|  2.656733| 0.2452762|phia_interaction |11_Chla.per.cell.vol.fg.fL |lowFeCu vs lowCu |TO 1005 :  low |      16|
|  1.5876667|  1|  3.7810282|  9.511717| 0.0213381|phia_interaction |11_Chla.per.cell.vol.fg.fL |lowFeCu vs lowFe |TO 1003 :  low |      16|
| -0.6305858|  1|  0.5964577|  1.500474| 0.4766535|phia_interaction |11_Chla.per.cell.vol.fg.fL |lowFeCu vs lowFe |TO 1005 :  low |      16|




<a id="Chla.per.cell.pg.cell"></a>

#### Chla.per.cell.pg.cell lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_Chla.per.cell.pg.cell <- lm(data=mydata, Chla.per.cell.pg.cell~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_Chla.per.cell.pg.cell_lowCu <- testInteractions(lm_Chla.per.cell.pg.cell, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_Chla.per.cell.pg.cell_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                    Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high -0.071000  1  0.007561  2.4221 0.170038   
## TO 1005 : high  0.206667  1  0.064067 20.5221 0.001366 **
## TO 1003 :  low  0.178667  1  0.047883 15.3380 0.003692 **
## TO 1005 :  low  0.083757  1  0.010523  3.3707 0.170038   
## Residuals                16  0.049949                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Chla.per.cell.pg.cell_lowCu <- phia_Chla.per.cell.pg.cell_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "10_Chla.per.cell.pg.cell",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_Chla.per.cell.pg.cell_lowCu), DF_here = phia_Chla.per.cell.pg.cell_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_Chla.per.cell.pg.cell_lowFe <- testInteractions(lm_Chla.per.cell.pg.cell, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_Chla.per.cell.pg.cell_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                    Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high -0.043333  1  0.002817  0.9022 0.712616   
## TO 1005 : high  0.121667  1  0.022204  7.1125 0.050627 . 
## TO 1003 :  low  0.206333  1  0.063860 20.4560 0.001387 **
## TO 1005 :  low -0.001243  1  0.000002  0.0007 0.978595   
## Residuals                16  0.049949                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Chla.per.cell.pg.cell_lowFe <- phia_Chla.per.cell.pg.cell_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "10_Chla.per.cell.pg.cell",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_Chla.per.cell.pg.cell_lowFe), DF_here = phia_Chla.per.cell.pg.cell_lowFe[5,2]) %>% 
     slice (3:4)  





#Making Chla.per.cell.pg.cell data.frame
kable(( Chla.per.cell.pg.cell <- bind_rows(phia_Chla.per.cell.pg.cell_lowCu, phia_Chla.per.cell.pg.cell_lowFe)), format = "markdown")
```



|      Value| Df| Sum of Sq|          F|    Pr(>F)|Analysis         |Data                     |test_intrctn     |tested         | DF_here|
|----------:|--:|---------:|----------:|---------:|:----------------|:------------------------|:----------------|:--------------|-------:|
|  0.1786667|  1| 0.0478827| 15.3379888| 0.0036918|phia_interaction |10_Chla.per.cell.pg.cell |lowFeCu vs lowCu |TO 1003 :  low |      16|
|  0.0837567|  1| 0.0105228|  3.3707051| 0.1700384|phia_interaction |10_Chla.per.cell.pg.cell |lowFeCu vs lowCu |TO 1005 :  low |      16|
|  0.2063333|  1| 0.0638602| 20.4559727| 0.0013871|phia_interaction |10_Chla.per.cell.pg.cell |lowFeCu vs lowFe |TO 1003 :  low |      16|
| -0.0012433|  1| 0.0000023|  0.0007427| 0.9785954|phia_interaction |10_Chla.per.cell.pg.cell |lowFeCu vs lowFe |TO 1005 :  low |      16|



<a id="cell.size..um"></a>

#### cell.size..um lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_cell.size..um <- lm(data=mydata, cell.size..um~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_cell.size..um_lowCu <- testInteractions(lm_cell.size..um, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_cell.size..um_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F  Pr(>F)  
## TO 1003 : high -0.13333  1   0.02667  0.8245 0.75466  
## TO 1005 : high  0.49333  1   0.36507 11.2879 0.01594 *
## TO 1003 :  low  0.02333  1   0.00082  0.0253 0.87573  
## TO 1005 :  low  0.26333  1   0.10402  3.2162 0.27550  
## Residuals               16   0.51746                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_cell.size..um_lowCu <- phia_cell.size..um_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "06_cell.size..um",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_cell.size..um_lowCu), DF_here = phia_cell.size..um_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_cell.size..um_lowFe <- testInteractions(lm_cell.size..um, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_cell.size..um_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high 0.50000  1   0.37500 11.5950 0.007241 **
## TO 1005 : high 0.59667  1   0.53402 16.5118 0.002710 **
## TO 1003 :  low 0.65667  1   0.64682 19.9996 0.001541 **
## TO 1005 :  low 0.36667  1   0.20167  6.2355 0.023812 * 
## Residuals              16   0.51746                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_cell.size..um_lowFe <- phia_cell.size..um_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "06_cell.size..um",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_cell.size..um_lowFe), DF_here = phia_cell.size..um_lowFe[5,2]) %>% 
     slice (3:4)  





#Making cell.size..um data.frame
kable(( cell.size..um <- bind_rows(phia_cell.size..um_lowCu, phia_cell.size..um_lowFe)), format = "markdown")
```



|     Value| Df| Sum of Sq|          F|    Pr(>F)|Analysis         |Data             |test_intrctn     |tested         | DF_here|
|---------:|--:|---------:|----------:|---------:|:----------------|:----------------|:----------------|:--------------|-------:|
| 0.0233333|  1| 0.0008167|  0.0252514| 0.8757308|phia_interaction |06_cell.size..um |lowFeCu vs lowCu |TO 1003 :  low |      16|
| 0.2633333|  1| 0.1040167|  3.2161980| 0.2755028|phia_interaction |06_cell.size..um |lowFeCu vs lowCu |TO 1005 :  low |      16|
| 0.6566667|  1| 0.6468167| 19.9995877| 0.0015406|phia_interaction |06_cell.size..um |lowFeCu vs lowFe |TO 1003 :  low |      16|
| 0.3666667|  1| 0.2016667|  6.2355384| 0.0238124|phia_interaction |06_cell.size..um |lowFeCu vs lowFe |TO 1005 :  low |      16|




<a id="cell.volume.fl.cell"></a>

#### cell.volume.fl.cell lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_cell.volume.fl.cell <- lm(data=mydata, cell.volume.fl.cell~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_cell.volume.fl.cell_lowCu <- testInteractions(lm_cell.volume.fl.cell, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_cell.volume.fl.cell_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F Pr(>F)  
## TO 1003 : high -5.7500  1     49.59  0.9048 0.7113  
## TO 1005 : high 21.2433  1    676.92 12.3499 0.0115 *
## TO 1003 :  low  0.7667  1      0.88  0.0161 0.9007  
## TO 1005 :  low  9.9867  1    149.60  2.7294 0.3540  
## Residuals              16    876.99                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_cell.volume.fl.cell_lowCu <- phia_cell.volume.fl.cell_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "07_cell.volume.fl.cell",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_cell.volume.fl.cell_lowCu), DF_here = phia_cell.volume.fl.cell_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_cell.volume.fl.cell_lowFe <- testInteractions(lm_cell.volume.fl.cell, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_cell.volume.fl.cell_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                 Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high 19.520  1    571.55 10.4275 0.010492 * 
## TO 1005 : high 24.157  1    875.32 15.9696 0.003121 **
## TO 1003 :  low 26.037  1   1016.86 18.5519 0.002170 **
## TO 1005 :  low 12.900  1    249.62  4.5541 0.048664 * 
## Residuals             16    876.99                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_cell.volume.fl.cell_lowFe <- phia_cell.volume.fl.cell_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "07_cell.volume.fl.cell",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_cell.volume.fl.cell_lowFe), DF_here = phia_cell.volume.fl.cell_lowFe[5,2]) %>% 
     slice (3:4)  





#Making cell.volume.fl.cell data.frame
kable(( cell.volume.fl.cell <- bind_rows(phia_cell.volume.fl.cell_lowCu, phia_cell.volume.fl.cell_lowFe)), format = "markdown")
```



|      Value| Df|    Sum of Sq|          F|    Pr(>F)|Analysis         |Data                   |test_intrctn     |tested         | DF_here|
|----------:|--:|------------:|----------:|---------:|:----------------|:----------------------|:----------------|:--------------|-------:|
|  0.7666667|  1|    0.8816667|  0.0160854| 0.9006565|phia_interaction |07_cell.volume.fl.cell |lowFeCu vs lowCu |TO 1003 :  low |      16|
|  9.9866667|  1|  149.6002667|  2.7293535| 0.3540156|phia_interaction |07_cell.volume.fl.cell |lowFeCu vs lowCu |TO 1005 :  low |      16|
| 26.0366667|  1| 1016.8620167| 18.5519449| 0.0021699|phia_interaction |07_cell.volume.fl.cell |lowFeCu vs lowFe |TO 1003 :  low |      16|
| 12.9000000|  1|  249.6150000|  4.5540532| 0.0486643|phia_interaction |07_cell.volume.fl.cell |lowFeCu vs lowFe |TO 1005 :  low |      16|



<a id="cell.SA.um2"></a>

#### cell.SA.um2 lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_cell.SA.um2 <- lm(data=mydata, cell.SA.um2~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_cell.SA.um2_lowCu <- testInteractions(lm_cell.SA.um2, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_cell.SA.um2_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F  Pr(>F)  
## TO 1003 : high -4.4202  1     29.31  0.8782 0.72524  
## TO 1005 : high 16.1893  1    393.14 11.7807 0.01368 *
## TO 1003 :  low  0.6973  1      0.73  0.0219 0.88431  
## TO 1005 :  low  8.0844  1     98.04  2.9377 0.31750  
## Residuals              16    533.94                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_cell.SA.um2_lowCu <- phia_cell.SA.um2_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "40_08_cell.SA.um2",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_cell.SA.um2_lowCu), DF_here = phia_cell.SA.um2_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_cell.SA.um2_lowFe <- testInteractions(lm_cell.SA.um2, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_cell.SA.um2_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                 Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high 15.627  1    366.30 10.9765 0.008793 **
## TO 1005 : high 19.024  1    542.90 16.2684 0.002886 **
## TO 1003 :  low 20.744  1    645.50 19.3429 0.001796 **
## TO 1005 :  low 10.920  1    178.86  5.3597 0.034214 * 
## Residuals             16    533.94                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_cell.SA.um2_lowFe <- phia_cell.SA.um2_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "40_08_cell.SA.um2",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_cell.SA.um2_lowFe), DF_here = phia_cell.SA.um2_lowFe[5,2]) %>% 
     slice (3:4)  





#Making cell.SA.um2 data.frame
kable(( cell.SA.um2 <- bind_rows(phia_cell.SA.um2_lowCu, phia_cell.SA.um2_lowFe)), format = "markdown")
```



|      Value| Df|   Sum of Sq|          F|    Pr(>F)|Analysis         |Data              |test_intrctn     |tested         | DF_here|
|----------:|--:|-----------:|----------:|---------:|:----------------|:-----------------|:----------------|:--------------|-------:|
|  0.6973288|  1|   0.7294013|  0.0218571| 0.8843150|phia_interaction |40_08_cell.SA.um2 |lowFeCu vs lowCu |TO 1003 :  low |      16|
|  8.0844028|  1|  98.0363529|  2.9377385| 0.3175044|phia_interaction |40_08_cell.SA.um2 |lowFeCu vs lowCu |TO 1005 :  low |      16|
| 20.7444599|  1| 645.4989242| 19.3428964| 0.0017962|phia_interaction |40_08_cell.SA.um2 |lowFeCu vs lowFe |TO 1003 :  low |      16|
| 10.9196902|  1| 178.8594500|  5.3596678| 0.0342139|phia_interaction |40_08_cell.SA.um2 |lowFeCu vs lowFe |TO 1005 :  low |      16|




<a id="cell.SA.Vol.ratio"></a>

#### cell.SA.Vol.ratio lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_cell.SA.Vol.ratio <- lm(data=mydata, cell.SA.Vol.ratio~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_cell.SA.Vol.ratio_lowCu <- testInteractions(lm_cell.SA.Vol.ratio, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_cell.SA.Vol.ratio_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                    Value Df Sum of Sq       F  Pr(>F)  
## TO 1003 : high  0.027490  1 0.0011335  0.6614 0.85597  
## TO 1005 : high -0.110274  1 0.0182405 10.6438 0.01956 *
## TO 1003 :  low -0.004758  1 0.0000340  0.0198 0.88982  
## TO 1005 :  low -0.068830  1 0.0071063  4.1467 0.17585  
## Residuals                16 0.0274195                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_cell.SA.Vol.ratio_lowCu <- phia_cell.SA.Vol.ratio_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "09_cell.SA.Vol.ratio",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_cell.SA.Vol.ratio_lowCu), DF_here = phia_cell.SA.Vol.ratio_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_cell.SA.Vol.ratio_lowFe <- testInteractions(lm_cell.SA.Vol.ratio, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_cell.SA.Vol.ratio_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                    Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high -0.123008  1  0.022697 13.2441 0.004417 **
## TO 1005 : high -0.138688  1  0.028852 16.8357 0.002493 **
## TO 1003 :  low -0.155256  1  0.036157 21.0983 0.001200 **
## TO 1005 :  low -0.097244  1  0.014185  8.2771 0.010951 * 
## Residuals                16  0.027419                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_cell.SA.Vol.ratio_lowFe <- phia_cell.SA.Vol.ratio_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "09_cell.SA.Vol.ratio",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_cell.SA.Vol.ratio_lowFe), DF_here = phia_cell.SA.Vol.ratio_lowFe[5,2]) %>% 
     slice (3:4)  





#Making cell.SA.Vol.ratio data.frame
kable(( cell.SA.Vol.ratio <- bind_rows(phia_cell.SA.Vol.ratio_lowCu, phia_cell.SA.Vol.ratio_lowFe)), format = "markdown")
```



|      Value| Df| Sum of Sq|          F|    Pr(>F)|Analysis         |Data                 |test_intrctn     |tested         | DF_here|
|----------:|--:|---------:|----------:|---------:|:----------------|:--------------------|:----------------|:--------------|-------:|
| -0.0047575|  1| 0.0000340|  0.0198115| 0.8898215|phia_interaction |09_cell.SA.Vol.ratio |lowFeCu vs lowCu |TO 1003 :  low |      16|
| -0.0688297|  1| 0.0071063|  4.1467116| 0.1758528|phia_interaction |09_cell.SA.Vol.ratio |lowFeCu vs lowCu |TO 1005 :  low |      16|
| -0.1552558|  1| 0.0361565| 21.0983112| 0.0011995|phia_interaction |09_cell.SA.Vol.ratio |lowFeCu vs lowFe |TO 1003 :  low |      16|
| -0.0972440|  1| 0.0141846|  8.2770937| 0.0109513|phia_interaction |09_cell.SA.Vol.ratio |lowFeCu vs lowFe |TO 1005 :  low |      16|

<a id="FeDFB.zmol.cell.1.h."></a>

#### FeDFB.zmol.cell.1.h. lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_FeDFB.zmol.cell.1.h. <- lm(data=mydata, FeDFB.zmol.cell.1.h.~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_FeDFB.zmol.cell.1.h._lowCu <- testInteractions(lm_FeDFB.zmol.cell.1.h., fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_FeDFB.zmol.cell.1.h._lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high  258.66  1    100360  2.8887 0.1113020    
## TO 1005 : high -478.16  1    274360  7.8971 0.0277984 *  
## TO 1003 :  low -666.73  1    666800 19.1929 0.0018822 ** 
## TO 1005 :  low -838.68  1    844061 24.2951 0.0008878 ***
## Residuals              14    486389                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_FeDFB.zmol.cell.1.h._lowCu <- phia_FeDFB.zmol.cell.1.h._lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "41_FeDFB.zmol.cell.1.h.",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_FeDFB.zmol.cell.1.h._lowCu), DF_here = phia_FeDFB.zmol.cell.1.h._lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_FeDFB.zmol.cell.1.h._lowFe <- testInteractions(lm_FeDFB.zmol.cell.1.h., fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_FeDFB.zmol.cell.1.h._lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high  292.63  1    128448  3.6972 0.150174   
## TO 1005 : high  412.16  1    203851  5.8676 0.088730 . 
## TO 1003 :  low -632.77  1    600590 17.2871 0.003867 **
## TO 1005 :  low   51.64  1      3200  0.0921 0.765992   
## Residuals              14    486389                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_FeDFB.zmol.cell.1.h._lowFe <- phia_FeDFB.zmol.cell.1.h._lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "41_FeDFB.zmol.cell.1.h.",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_FeDFB.zmol.cell.1.h._lowFe), DF_here = phia_FeDFB.zmol.cell.1.h._lowFe[5,2]) %>% 
     slice (3:4)  





#Making FeDFB.zmol.cell.1.h. data.frame
kable(( FeDFB.zmol.cell.1.h. <- bind_rows(phia_FeDFB.zmol.cell.1.h._lowCu, phia_FeDFB.zmol.cell.1.h._lowFe)), format = "markdown")
```



|      Value| Df|  Sum of Sq|          F|    Pr(>F)|Analysis         |Data                    |test_intrctn     |tested         | DF_here|
|----------:|--:|----------:|----------:|---------:|:----------------|:-----------------------|:----------------|:--------------|-------:|
| -666.73333|  1| 666800.007| 19.1928569| 0.0018822|phia_interaction |41_FeDFB.zmol.cell.1.h. |lowFeCu vs lowCu |TO 1003 :  low |      14|
| -838.68002|  1| 844061.019| 24.2950543| 0.0008878|phia_interaction |41_FeDFB.zmol.cell.1.h. |lowFeCu vs lowCu |TO 1005 :  low |      14|
| -632.76667|  1| 600590.482| 17.2871132| 0.0038672|phia_interaction |41_FeDFB.zmol.cell.1.h. |lowFeCu vs lowFe |TO 1003 :  low |      14|
|   51.63658|  1|   3199.604|  0.0920959| 0.7659925|phia_interaction |41_FeDFB.zmol.cell.1.h. |lowFeCu vs lowFe |TO 1005 :  low |      14|



<a id="FeDFB.zmol.um.2.h."></a>

#### FeDFB.zmol.um.2.h. lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_FeDFB.zmol.um.2.h. <- lm(data=mydata, FeDFB.zmol.um.2.h.~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_FeDFB.zmol.um.2.h._lowCu <- testInteractions(lm_FeDFB.zmol.um.2.h., fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_FeDFB.zmol.um.2.h._lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high   3.0000  1    13.500  1.5766 0.2298039    
## TO 1005 : high  -8.8351  1    93.671 10.9394 0.0103693 *  
## TO 1003 :  low  -9.6333  1   139.202 16.2567 0.0037085 ** 
## TO 1005 :  low -13.7485  1   226.825 26.4897 0.0005934 ***
## Residuals               14   119.878                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_FeDFB.zmol.um.2.h._lowCu <- phia_FeDFB.zmol.um.2.h._lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "29_FeDFB.zmol.um.2.h.",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_FeDFB.zmol.um.2.h._lowCu), DF_here = phia_FeDFB.zmol.um.2.h._lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_FeDFB.zmol.um.2.h._lowFe <- testInteractions(lm_FeDFB.zmol.um.2.h., fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_FeDFB.zmol.um.2.h._lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F  Pr(>F)   
## TO 1003 : high  3.1667  1    15.042  1.7566 0.61880   
## TO 1005 : high  3.4551  1    14.326  1.6730 0.61880   
## TO 1003 :  low -9.4667  1   134.427 15.6990 0.00567 **
## TO 1005 :  low -1.4582  1     2.552  0.2980 0.61880   
## Residuals              14   119.878                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_FeDFB.zmol.um.2.h._lowFe <- phia_FeDFB.zmol.um.2.h._lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "29_FeDFB.zmol.um.2.h.",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_FeDFB.zmol.um.2.h._lowFe), DF_here = phia_FeDFB.zmol.um.2.h._lowFe[5,2]) %>% 
     slice (3:4)  





#Making FeDFB.zmol.um.2.h. data.frame
kable(( FeDFB.zmol.um.2.h. <- bind_rows(phia_FeDFB.zmol.um.2.h._lowCu, phia_FeDFB.zmol.um.2.h._lowFe)), format = "markdown")
```



|      Value| Df|  Sum of Sq|          F|    Pr(>F)|Analysis         |Data                  |test_intrctn     |tested         | DF_here|
|----------:|--:|----------:|----------:|---------:|:----------------|:---------------------|:----------------|:--------------|-------:|
|  -9.633333|  1| 139.201667| 16.2566755| 0.0037085|phia_interaction |29_FeDFB.zmol.um.2.h. |lowFeCu vs lowCu |TO 1003 :  low |      14|
| -13.748473|  1| 226.824622| 26.4897279| 0.0005934|phia_interaction |29_FeDFB.zmol.um.2.h. |lowFeCu vs lowCu |TO 1005 :  low |      14|
|  -9.466667|  1| 134.426667| 15.6990268| 0.0056695|phia_interaction |29_FeDFB.zmol.um.2.h. |lowFeCu vs lowFe |TO 1003 :  low |      14|
|  -1.458227|  1|   2.551711|  0.2980017| 0.6188015|phia_interaction |29_FeDFB.zmol.um.2.h. |lowFeCu vs lowFe |TO 1005 :  low |      14|



<a id="X14C.per.Chla.at.155uE"></a>

#### X14C.per.Chla.at.155uE lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_X14C.per.Chla.at.155uE <- lm(data=mydata, X14C.per.Chla.at.155uE~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_X14C.per.Chla.at.155uE_lowCu <- testInteractions(lm_X14C.per.Chla.at.155uE, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_X14C.per.Chla.at.155uE_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq      F Pr(>F)
## TO 1003 : high  0.17667  1    0.0468 0.1783 0.8043
## TO 1005 : high -0.41667  1    0.2604 0.9920 0.8043
## TO 1003 :  low -0.48000  1    0.3456 1.3164 0.8043
## TO 1005 :  low  0.94333  1    1.3348 5.0845 0.1540
## Residuals               16    4.2004
```

```r
phia_X14C.per.Chla.at.155uE_lowCu <- phia_X14C.per.Chla.at.155uE_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "28_X14C.per.Chla.at.155uE",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_X14C.per.Chla.at.155uE_lowCu), DF_here = phia_X14C.per.Chla.at.155uE_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_X14C.per.Chla.at.155uE_lowFe <- testInteractions(lm_X14C.per.Chla.at.155uE, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_X14C.per.Chla.at.155uE_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F  Pr(>F)  
## TO 1003 : high  1.39333  1    2.9121 11.0925 0.01695 *
## TO 1005 : high -0.77333  1    0.8971  3.4171 0.24928  
## TO 1003 :  low  0.73667  1    0.8140  3.1007 0.24928  
## TO 1005 :  low  0.58667  1    0.5163  1.9665 0.24928  
## Residuals               16    4.2004                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_X14C.per.Chla.at.155uE_lowFe <- phia_X14C.per.Chla.at.155uE_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "28_X14C.per.Chla.at.155uE",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_X14C.per.Chla.at.155uE_lowFe), DF_here = phia_X14C.per.Chla.at.155uE_lowFe[5,2]) %>% 
     slice (3:4)  





#Making X14C.per.Chla.at.155uE data.frame
kable(( X14C.per.Chla.at.155uE <- bind_rows(phia_X14C.per.Chla.at.155uE_lowCu, phia_X14C.per.Chla.at.155uE_lowFe)), format = "markdown")
```



|      Value| Df| Sum of Sq|        F|    Pr(>F)|Analysis         |Data                      |test_intrctn     |tested         | DF_here|
|----------:|--:|---------:|--------:|---------:|:----------------|:-------------------------|:----------------|:--------------|-------:|
| -0.4800000|  1| 0.3456000| 1.316446| 0.8042959|phia_interaction |28_X14C.per.Chla.at.155uE |lowFeCu vs lowCu |TO 1003 :  low |      16|
|  0.9433333|  1| 1.3348167| 5.084532| 0.1540156|phia_interaction |28_X14C.per.Chla.at.155uE |lowFeCu vs lowCu |TO 1005 :  low |      16|
|  0.7366667|  1| 0.8140167| 3.100721| 0.2492769|phia_interaction |28_X14C.per.Chla.at.155uE |lowFeCu vs lowFe |TO 1003 :  low |      16|
|  0.5866667|  1| 0.5162667| 1.966543| 0.2492769|phia_interaction |28_X14C.per.Chla.at.155uE |lowFeCu vs lowFe |TO 1005 :  low |      16|



<a id="X14C.per.Chla.alpha"></a>

#### X14C.per.Chla.alpha lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_X14C.per.Chla.alpha <- lm(data=mydata, X14C.per.Chla.alpha~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_X14C.per.Chla.alpha_lowCu <- testInteractions(lm_X14C.per.Chla.alpha, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_X14C.per.Chla.alpha_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                     Value Df  Sum of Sq      F Pr(>F)
## TO 1003 : high -0.0007496  1 8.4300e-07 0.0506 0.9746
## TO 1005 : high -0.0032710  1 1.6049e-05 0.9645 0.9746
## TO 1003 :  low -0.0033831  1 1.7168e-05 1.0317 0.9746
## TO 1005 :  low  0.0059094  1 5.2381e-05 3.1479 0.3802
## Residuals                 16 2.6624e-04
```

```r
phia_X14C.per.Chla.alpha_lowCu <- phia_X14C.per.Chla.alpha_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "25_X14C.per.Chla.alpha",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_X14C.per.Chla.alpha_lowCu), DF_here = phia_X14C.per.Chla.alpha_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_X14C.per.Chla.alpha_lowFe <- testInteractions(lm_X14C.per.Chla.alpha, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_X14C.per.Chla.alpha_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                     Value Df  Sum of Sq      F Pr(>F)  
## TO 1003 : high  0.0094690  1 1.3449e-04 8.0825 0.0470 *
## TO 1005 : high -0.0064355  1 6.2124e-05 3.7334 0.1706  
## TO 1003 :  low  0.0068355  1 7.0085e-05 4.2119 0.1706  
## TO 1005 :  low  0.0027448  1 1.1301e-05 0.6792 0.4220  
## Residuals                 16 2.6624e-04                
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_X14C.per.Chla.alpha_lowFe <- phia_X14C.per.Chla.alpha_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "25_X14C.per.Chla.alpha",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_X14C.per.Chla.alpha_lowFe), DF_here = phia_X14C.per.Chla.alpha_lowFe[5,2]) %>% 
     slice (3:4)  





#Making X14C.per.Chla.alpha data.frame
kable(( X14C.per.Chla.alpha <- bind_rows(phia_X14C.per.Chla.alpha_lowCu, phia_X14C.per.Chla.alpha_lowFe)), format = "markdown")
```



|      Value| Df| Sum of Sq|         F|    Pr(>F)|Analysis         |Data                   |test_intrctn     |tested         | DF_here|
|----------:|--:|---------:|---------:|---------:|:----------------|:----------------------|:----------------|:--------------|-------:|
| -0.0033831|  1|  1.72e-05| 1.0317440| 0.9745922|phia_interaction |25_X14C.per.Chla.alpha |lowFeCu vs lowCu |TO 1003 :  low |      16|
|  0.0059094|  1|  5.24e-05| 3.1478791| 0.3802144|phia_interaction |25_X14C.per.Chla.alpha |lowFeCu vs lowCu |TO 1005 :  low |      16|
|  0.0068355|  1|  7.01e-05| 4.2118584| 0.1706296|phia_interaction |25_X14C.per.Chla.alpha |lowFeCu vs lowFe |TO 1003 :  low |      16|
|  0.0027448|  1|  1.13e-05| 0.6791589| 0.4219921|phia_interaction |25_X14C.per.Chla.alpha |lowFeCu vs lowFe |TO 1005 :  low |      16|



<a id="X14C.per.Chla.ek"></a>

#### X14C.per.Chla.ek lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_X14C.per.Chla.ek <- lm(data=mydata, X14C.per.Chla.ek~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_X14C.per.Chla.ek_lowCu <- testInteractions(lm_X14C.per.Chla.ek, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_X14C.per.Chla.ek_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high  89.459  1   12004.4 26.2884 0.0004054 ***
## TO 1005 : high  20.054  1     603.3  1.3211 0.2672881    
## TO 1003 :  low -29.026  1    1263.8  2.7676 0.2313016    
## TO 1005 :  low  60.212  1    5438.3 11.9093 0.0098603 ** 
## Residuals              16    7306.3                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_X14C.per.Chla.ek_lowCu <- phia_X14C.per.Chla.ek_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "26_X14C.per.Chla.ek",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_X14C.per.Chla.ek_lowCu), DF_here = phia_X14C.per.Chla.ek_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_X14C.per.Chla.ek_lowFe <- testInteractions(lm_X14C.per.Chla.ek, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_X14C.per.Chla.ek_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high  79.971  1    9593.1 21.0077 0.001224 **
## TO 1005 : high  26.013  1    1015.0  2.2228 0.155434   
## TO 1003 :  low -38.515  1    2225.1  4.8727 0.084470 . 
## TO 1005 :  low  66.171  1    6568.0 14.3831 0.004793 **
## Residuals              16    7306.3                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_X14C.per.Chla.ek_lowFe <- phia_X14C.per.Chla.ek_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "26_X14C.per.Chla.ek",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_X14C.per.Chla.ek_lowFe), DF_here = phia_X14C.per.Chla.ek_lowFe[5,2]) %>% 
     slice (3:4)  





#Making X14C.per.Chla.ek data.frame
kable(( X14C.per.Chla.ek <- bind_rows(phia_X14C.per.Chla.ek_lowCu, phia_X14C.per.Chla.ek_lowFe)), format = "markdown")
```



|     Value| Df| Sum of Sq|         F|    Pr(>F)|Analysis         |Data                |test_intrctn     |tested         | DF_here|
|---------:|--:|---------:|---------:|---------:|:----------------|:-------------------|:----------------|:--------------|-------:|
| -29.02647|  1|  1263.804|  2.767590| 0.2313016|phia_interaction |26_X14C.per.Chla.ek |lowFeCu vs lowCu |TO 1003 :  low |      16|
|  60.21241|  1|  5438.302| 11.909274| 0.0098603|phia_interaction |26_X14C.per.Chla.ek |lowFeCu vs lowCu |TO 1005 :  low |      16|
| -38.51467|  1|  2225.070|  4.872654| 0.0844697|phia_interaction |26_X14C.per.Chla.ek |lowFeCu vs lowFe |TO 1003 :  low |      16|
|  66.17141|  1|  6567.982| 14.383147| 0.0047933|phia_interaction |26_X14C.per.Chla.ek |lowFeCu vs lowFe |TO 1005 :  low |      16|



<a id="X14C.per.Chla.pmax"></a>

#### X14C.per.Chla.pmax lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_X14C.per.Chla.pmax <- lm(data=mydata, X14C.per.Chla.pmax~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_X14C.per.Chla.pmax_lowCu <- testInteractions(lm_X14C.per.Chla.pmax, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_X14C.per.Chla.pmax_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq      F Pr(>F)  
## TO 1003 : high  1.30308  1    2.5470 3.6461 0.2229  
## TO 1005 : high -0.54506  1    0.4456 0.6379 0.4570  
## TO 1003 :  low -0.85450  1    1.0953 1.5679 0.4570  
## TO 1005 :  low  2.01792  1    6.1080 8.7437 0.0371 *
## Residuals               16   11.1770                
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_X14C.per.Chla.pmax_lowCu <- phia_X14C.per.Chla.pmax_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "27_X14C.per.Chla.pmax",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_X14C.per.Chla.pmax_lowCu), DF_here = phia_X14C.per.Chla.pmax_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_X14C.per.Chla.pmax_lowFe <- testInteractions(lm_X14C.per.Chla.pmax, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_X14C.per.Chla.pmax_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high  2.92435  1   12.8277 18.3630 0.002272 **
## TO 1005 : high -0.93730  1    1.3178  1.8865 0.377086   
## TO 1003 :  low  0.76677  1    0.8819  1.2624 0.377086   
## TO 1005 :  low  1.62568  1    3.9642  5.6748 0.089885 . 
## Residuals               16   11.1770                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_X14C.per.Chla.pmax_lowFe <- phia_X14C.per.Chla.pmax_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "27_X14C.per.Chla.pmax",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_X14C.per.Chla.pmax_lowFe), DF_here = phia_X14C.per.Chla.pmax_lowFe[5,2]) %>% 
     slice (3:4)  





#Making X14C.per.Chla.pmax data.frame
kable(( X14C.per.Chla.pmax <- bind_rows(phia_X14C.per.Chla.pmax_lowCu, phia_X14C.per.Chla.pmax_lowFe)), format = "markdown")
```



|      Value| Df| Sum of Sq|        F|    Pr(>F)|Analysis         |Data                  |test_intrctn     |tested         | DF_here|
|----------:|--:|---------:|--------:|---------:|:----------------|:---------------------|:----------------|:--------------|-------:|
| -0.8545003|  1| 1.0952562| 1.567868| 0.4570240|phia_interaction |27_X14C.per.Chla.pmax |lowFeCu vs lowCu |TO 1003 :  low |      16|
|  2.0179243|  1| 6.1080275| 8.743690| 0.0370997|phia_interaction |27_X14C.per.Chla.pmax |lowFeCu vs lowCu |TO 1005 :  low |      16|
|  0.7667661|  1| 0.8818955| 1.262441| 0.3770861|phia_interaction |27_X14C.per.Chla.pmax |lowFeCu vs lowFe |TO 1003 :  low |      16|
|  1.6256754|  1| 3.9642309| 5.674828| 0.0898849|phia_interaction |27_X14C.per.Chla.pmax |lowFeCu vs lowFe |TO 1005 :  low |      16|



<a id="AOXactivity"></a>

#### AOXactivity lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_AOXactivity <- lm(data=mydata, AOXactivity~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_AOXactivity_lowCu <- testInteractions(lm_AOXactivity, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_AOXactivity_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq      F Pr(>F)
## TO 1003 : high -19.7100  1     466.2 1.7194 0.8327
## TO 1005 : high   1.3683  1       2.2 0.0083 0.9517
## TO 1003 :  low  -9.8500  1     145.5 0.5368 0.9517
## TO 1005 :  low -17.7400  1     472.1 1.7411 0.8327
## Residuals               14    3795.8
```

```r
phia_AOXactivity_lowCu <- phia_AOXactivity_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "34_AOXactivity",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_AOXactivity_lowCu), DF_here = phia_AOXactivity_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_AOXactivity_lowFe <- testInteractions(lm_AOXactivity, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_AOXactivity_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq      F Pr(>F)
## TO 1003 : high  -4.7600  1      34.0 0.1254 1.0000
## TO 1005 : high  -3.1233  1      14.6 0.0540 1.0000
## TO 1003 :  low   5.1000  1      31.2 0.1151 1.0000
## TO 1005 :  low -22.2317  1     593.1 2.1875 0.6451
## Residuals               14    3795.8
```

```r
phia_AOXactivity_lowFe <- phia_AOXactivity_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "34_AOXactivity",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_AOXactivity_lowFe), DF_here = phia_AOXactivity_lowFe[5,2]) %>% 
     slice (3:4)  





#Making AOXactivity data.frame
kable(( AOXactivity <- bind_rows(phia_AOXactivity_lowCu, phia_AOXactivity_lowFe)), format = "markdown")
```



|     Value| Df| Sum of Sq|         F|    Pr(>F)|Analysis         |Data           |test_intrctn     |tested         | DF_here|
|---------:|--:|---------:|---------:|---------:|:----------------|:--------------|:----------------|:--------------|-------:|
|  -9.85000|  1|  145.5337| 0.5367732| 0.9517230|phia_interaction |34_AOXactivity |lowFeCu vs lowCu |TO 1003 :  low |      14|
| -17.74000|  1|  472.0614| 1.7411075| 0.8326847|phia_interaction |34_AOXactivity |lowFeCu vs lowCu |TO 1005 :  low |      14|
|   5.10000|  1|   31.2120| 0.1151194| 1.0000000|phia_interaction |34_AOXactivity |lowFeCu vs lowFe |TO 1003 :  low |      14|
| -22.23167|  1|  593.0964| 2.1875217| 0.6451168|phia_interaction |34_AOXactivity |lowFeCu vs lowFe |TO 1005 :  low |      14|





<a id="GrossPchla.mol.O2.mol.Chla..h."></a>

#### GrossPchla.mol.O2.mol.Chla..h. lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_GrossPchla.mol.O2.mol.Chla..h. <- lm(data=mydata, GrossPchla.mol.O2.mol.Chla..h.~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_GrossPchla.mol.O2.mol.Chla..h._lowCu <- testInteractions(lm_GrossPchla.mol.O2.mol.Chla..h., fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_GrossPchla.mol.O2.mol.Chla..h._lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F  Pr(>F)  
## TO 1003 : high  271.280  1     88311 12.5120 0.01313 *
## TO 1005 : high   77.228  1      7157  1.0140 0.66207  
## TO 1003 :  low -157.193  1     37065  5.2513 0.11385  
## TO 1005 :  low  -57.167  1      4902  0.6945 0.66207  
## Residuals               14     98814                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_GrossPchla.mol.O2.mol.Chla..h._lowCu <- phia_GrossPchla.mol.O2.mol.Chla..h._lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "30_GrossPchla.mol.O2.mol.Chla..h.",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_GrossPchla.mol.O2.mol.Chla..h._lowCu), DF_here = phia_GrossPchla.mol.O2.mol.Chla..h._lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_GrossPchla.mol.O2.mol.Chla..h._lowFe <- testInteractions(lm_GrossPchla.mol.O2.mol.Chla..h., fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_GrossPchla.mol.O2.mol.Chla..h._lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                 Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high 380.79  1    217505 30.8163 0.0002856 ***
## TO 1005 : high  66.93  1      6719  0.9519 1.0000000    
## TO 1003 :  low -47.68  1      2728  0.3865 1.0000000    
## TO 1005 :  low -67.47  1      5462  0.7739 1.0000000    
## Residuals             14     98814                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_GrossPchla.mol.O2.mol.Chla..h._lowFe <- phia_GrossPchla.mol.O2.mol.Chla..h._lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "30_GrossPchla.mol.O2.mol.Chla..h.",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_GrossPchla.mol.O2.mol.Chla..h._lowFe), DF_here = phia_GrossPchla.mol.O2.mol.Chla..h._lowFe[5,2]) %>% 
     slice (3:4)  





#Making GrossPchla.mol.O2.mol.Chla..h. data.frame
kable(( GrossPchla.mol.O2.mol.Chla..h. <- bind_rows(phia_GrossPchla.mol.O2.mol.Chla..h._lowCu, phia_GrossPchla.mol.O2.mol.Chla..h._lowFe)), format = "markdown")
```



|      Value| Df| Sum of Sq|         F|    Pr(>F)|Analysis         |Data                              |test_intrctn     |tested         | DF_here|
|----------:|--:|---------:|---------:|---------:|:----------------|:---------------------------------|:----------------|:--------------|-------:|
| -157.19333|  1| 37064.616| 5.2513403| 0.1138532|phia_interaction |30_GrossPchla.mol.O2.mol.Chla..h. |lowFeCu vs lowCu |TO 1003 :  low |      14|
|  -57.16679|  1|  4902.063| 0.6945277| 0.6620708|phia_interaction |30_GrossPchla.mol.O2.mol.Chla..h. |lowFeCu vs lowCu |TO 1005 :  low |      14|
|  -47.68000|  1|  2728.059| 0.3865133| 1.0000000|phia_interaction |30_GrossPchla.mol.O2.mol.Chla..h. |lowFeCu vs lowFe |TO 1003 :  low |      14|
|  -67.46755|  1|  5462.244| 0.7738945| 1.0000000|phia_interaction |30_GrossPchla.mol.O2.mol.Chla..h. |lowFeCu vs lowFe |TO 1005 :  low |      14|



<a id="GrossPcell..mol.O2.cell..h."></a>

#### GrossPcell..mol.O2.cell..h. lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
mydata2 <- mydata %>%
  mutate(GrossPcell..mol.O2.cell..h._times1000000 = GrossPcell..mol.O2.cell..h.*1000000)

lm_GrossPcell..mol.O2.cell..h. <- lm(data=mydata2, GrossPcell..mol.O2.cell..h._times1000000~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_GrossPcell..mol.O2.cell..h._lowCu <- testInteractions(lm_GrossPcell..mol.O2.cell..h., fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_GrossPcell..mol.O2.cell..h._lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                    Value Df Sum of Sq       F  Pr(>F)   
## TO 1003 : high  0.083167  1 0.0083000 11.7372 0.01229 * 
## TO 1005 : high  0.093117  1 0.0104049 14.7137 0.00727 **
## TO 1003 :  low -0.011967  1 0.0002148  0.3038 1.00000   
## TO 1005 :  low  0.011800  1 0.0002089  0.2954 1.00000   
## Residuals                14 0.0099002                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_GrossPcell..mol.O2.cell..h._lowCu <- phia_GrossPcell..mol.O2.cell..h._lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "42_GrossPcell..mol.O2.cell..h.",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_GrossPcell..mol.O2.cell..h._lowCu), DF_here = phia_GrossPcell..mol.O2.cell..h._lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_GrossPcell..mol.O2.cell..h._lowFe <- testInteractions(lm_GrossPcell..mol.O2.cell..h., fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_GrossPcell..mol.O2.cell..h._lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                    Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high  0.132700  1 0.0264139 37.3524 0.0001075 ***
## TO 1005 : high  0.069233  1 0.0071899 10.1673 0.0197022 *  
## TO 1003 :  low  0.037567  1 0.0016935  2.3948 0.2880819    
## TO 1005 :  low -0.012083  1 0.0001752  0.2478 0.6263801    
## Residuals                14 0.0099002                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_GrossPcell..mol.O2.cell..h._lowFe <- phia_GrossPcell..mol.O2.cell..h._lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "42_GrossPcell..mol.O2.cell..h.",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_GrossPcell..mol.O2.cell..h._lowFe), DF_here = phia_GrossPcell..mol.O2.cell..h._lowFe[5,2]) %>% 
     slice (3:4)  





#Making GrossPcell..mol.O2.cell..h. data.frame
kable(( GrossPcell..mol.O2.cell..h. <- bind_rows(phia_GrossPcell..mol.O2.cell..h._lowCu, phia_GrossPcell..mol.O2.cell..h._lowFe)), format = "markdown")
```



|      Value| Df| Sum of Sq|         F|    Pr(>F)|Analysis         |Data                           |test_intrctn     |tested         | DF_here|
|----------:|--:|---------:|---------:|---------:|:----------------|:------------------------------|:----------------|:--------------|-------:|
| -0.0119667|  1| 0.0002148| 0.3037549| 1.0000000|phia_interaction |42_GrossPcell..mol.O2.cell..h. |lowFeCu vs lowCu |TO 1003 :  low |      14|
|  0.0118000|  1| 0.0002089| 0.2953527| 1.0000000|phia_interaction |42_GrossPcell..mol.O2.cell..h. |lowFeCu vs lowCu |TO 1005 :  low |      14|
|  0.0375667|  1| 0.0016935| 2.3948161| 0.2880819|phia_interaction |42_GrossPcell..mol.O2.cell..h. |lowFeCu vs lowFe |TO 1003 :  low |      14|
| -0.0120833|  1| 0.0001752| 0.2477652| 0.6263801|phia_interaction |42_GrossPcell..mol.O2.cell..h. |lowFeCu vs lowFe |TO 1005 :  low |      14|


<a id="Resp.mol.O2.mol.Chla.h"></a>

#### Resp.mol.O2.mol.Chla.h lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_Resp.mol.O2.mol.Chla.h <- lm(data=mydata, Resp.mol.O2.mol.Chla.h~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_Resp.mol.O2.mol.Chla.h_lowCu <- testInteractions(lm_Resp.mol.O2.mol.Chla.h, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_Resp.mol.O2.mol.Chla.h_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F Pr(>F)  
## TO 1003 : high -22.613  1     613.6  2.3238 0.4490  
## TO 1005 : high -16.698  1     334.6  1.2672 0.5585  
## TO 1003 :  low -47.053  1    3321.0 12.5769 0.0129 *
## TO 1005 :  low  -4.535  1      30.8  0.1168 0.7376  
## Residuals              14    3696.8                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Resp.mol.O2.mol.Chla.h_lowCu <- phia_Resp.mol.O2.mol.Chla.h_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "31_Resp.mol.O2.mol.Chla.h",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_Resp.mol.O2.mol.Chla.h_lowCu), DF_here = phia_Resp.mol.O2.mol.Chla.h_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_Resp.mol.O2.mol.Chla.h_lowFe <- testInteractions(lm_Resp.mol.O2.mol.Chla.h, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_Resp.mol.O2.mol.Chla.h_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq      F Pr(>F)
## TO 1003 : high   9.5994  1     138.2 0.5235      1
## TO 1005 : high -11.6384  1     203.2 0.7694      1
## TO 1003 :  low -14.8410  1     264.3 1.0009      1
## TO 1005 :  low   0.5250  1       0.3 0.0013      1
## Residuals               14    3696.8
```

```r
phia_Resp.mol.O2.mol.Chla.h_lowFe <- phia_Resp.mol.O2.mol.Chla.h_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "31_Resp.mol.O2.mol.Chla.h",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_Resp.mol.O2.mol.Chla.h_lowFe), DF_here = phia_Resp.mol.O2.mol.Chla.h_lowFe[5,2]) %>% 
     slice (3:4)  





#Making Resp.mol.O2.mol.Chla.h data.frame
kable(( Resp.mol.O2.mol.Chla.h <- bind_rows(phia_Resp.mol.O2.mol.Chla.h_lowCu, phia_Resp.mol.O2.mol.Chla.h_lowFe)), format = "markdown")
```



|       Value| Df|    Sum of Sq|          F|    Pr(>F)|Analysis         |Data                      |test_intrctn     |tested         | DF_here|
|-----------:|--:|------------:|----------:|---------:|:----------------|:-------------------------|:----------------|:--------------|-------:|
| -47.0533333|  1| 3321.0242667| 12.5768710| 0.0128959|phia_interaction |31_Resp.mol.O2.mol.Chla.h |lowFeCu vs lowCu |TO 1003 :  low |      14|
|  -4.5349586|  1|   30.8487745|  0.1168257| 0.7375795|phia_interaction |31_Resp.mol.O2.mol.Chla.h |lowFeCu vs lowCu |TO 1005 :  low |      14|
| -14.8409831|  1|  264.3057354|  1.0009379| 1.0000000|phia_interaction |31_Resp.mol.O2.mol.Chla.h |lowFeCu vs lowFe |TO 1003 :  low |      14|
|   0.5250447|  1|    0.3308064|  0.0012528| 1.0000000|phia_interaction |31_Resp.mol.O2.mol.Chla.h |lowFeCu vs lowFe |TO 1005 :  low |      14|



<a id="Resp.umol.O2.cell.h"></a>

#### Resp.umol.O2.cell.h lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
mydata2 <- mydata %>% 
  mutate(Resp.umol.O2.cell.h_times1000000 = Resp.umol.O2.cell.h * 1000000)

lm_Resp.umol.O2.cell.h <- lm(data=mydata2, Resp.umol.O2.cell.h_times1000000~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_Resp.umol.O2.cell.h_lowCu <- testInteractions(lm_Resp.umol.O2.cell.h, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_Resp.umol.O2.cell.h_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                     Value Df  Sum of Sq      F Pr(>F)  
## TO 1003 : high -0.0116367  1 0.00016249 7.1893 0.0716 .
## TO 1005 : high  0.0072833  1 0.00006366 2.8164 0.3464  
## TO 1003 :  low -0.0042033  1 0.00002650 1.1725 0.4089  
## TO 1005 :  low  0.0051667  1 0.00004004 1.7716 0.4089  
## Residuals                 14 0.00031643                
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Resp.umol.O2.cell.h_lowCu <- phia_Resp.umol.O2.cell.h_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "32_Resp.umol.O2.cell.h",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_Resp.umol.O2.cell.h_lowCu), DF_here = phia_Resp.umol.O2.cell.h_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_Resp.umol.O2.cell.h_lowFe <- testInteractions(lm_Resp.umol.O2.cell.h, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_Resp.umol.O2.cell.h_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                     Value Df  Sum of Sq      F Pr(>F)
## TO 1003 : high  0.0024667  1 0.00000913 0.4038 1.0000
## TO 1005 : high  0.0019700  1 0.00000582 0.2576 1.0000
## TO 1003 :  low  0.0099000  1 0.00011761 5.2035 0.1548
## TO 1005 :  low -0.0001467  1 0.00000003 0.0011 1.0000
## Residuals                 14 0.00031643
```

```r
phia_Resp.umol.O2.cell.h_lowFe <- phia_Resp.umol.O2.cell.h_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "32_Resp.umol.O2.cell.h",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_Resp.umol.O2.cell.h_lowFe), DF_here = phia_Resp.umol.O2.cell.h_lowFe[5,2]) %>% 
     slice (3:4)  





#Making Resp.umol.O2.cell.h data.frame
kable(( Resp.umol.O2.cell.h <- bind_rows(phia_Resp.umol.O2.cell.h_lowCu, phia_Resp.umol.O2.cell.h_lowFe)), format = "markdown")
```



|      Value| Df| Sum of Sq|         F|    Pr(>F)|Analysis         |Data                   |test_intrctn     |tested         | DF_here|
|----------:|--:|---------:|---------:|---------:|:----------------|:----------------------|:----------------|:--------------|-------:|
| -0.0042033|  1| 0.0000265| 1.1725358| 0.4089161|phia_interaction |32_Resp.umol.O2.cell.h |lowFeCu vs lowCu |TO 1003 :  low |      14|
|  0.0051667|  1| 0.0000400| 1.7715741| 0.4089161|phia_interaction |32_Resp.umol.O2.cell.h |lowFeCu vs lowCu |TO 1005 :  low |      14|
|  0.0099000|  1| 0.0001176| 5.2035390| 0.1548387|phia_interaction |32_Resp.umol.O2.cell.h |lowFeCu vs lowFe |TO 1003 :  low |      14|
| -0.0001467|  1| 0.0000000| 0.0011421| 1.0000000|phia_interaction |32_Resp.umol.O2.cell.h |lowFeCu vs lowFe |TO 1005 :  low |      14|



<a id="Gross.P.CellVol..mol.O2.L."></a>

#### Gross.P.CellVol..mol.O2.L. lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_Gross.P.CellVol..mol.O2.L. <- lm(data=mydata, Gross.P.CellVol..mol.O2.L.~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_Gross.P.CellVol..mol.O2.L._lowCu <- testInteractions(lm_Gross.P.CellVol..mol.O2.L., fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_Gross.P.CellVol..mol.O2.L._lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq      F Pr(>F)
## TO 1003 : high  0.60333  1   0.43681 3.5735 0.2388
## TO 1005 : high  0.73000  1   0.63948 5.2314 0.1531
## TO 1003 :  low -0.23333  1   0.08167 0.6681 0.8548
## TO 1005 :  low  0.17333  1   0.04507 0.3687 0.8548
## Residuals               14   1.71133
```

```r
phia_Gross.P.CellVol..mol.O2.L._lowCu <- phia_Gross.P.CellVol..mol.O2.L._lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "44_Gross.P.CellVol..mol.O2.L.",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_Gross.P.CellVol..mol.O2.L._lowCu), DF_here = phia_Gross.P.CellVol..mol.O2.L._lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_Gross.P.CellVol..mol.O2.L._lowFe <- testInteractions(lm_Gross.P.CellVol..mol.O2.L., fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_Gross.P.CellVol..mol.O2.L._lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high  1.20000  1   2.16000 17.6704 0.003537 **
## TO 1005 : high  0.27333  1   0.11207  0.9168 0.822187   
## TO 1003 :  low  0.36333  1   0.15841  1.2959 0.822187   
## TO 1005 :  low -0.28333  1   0.09633  0.7881 0.822187   
## Residuals               14   1.71133                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Gross.P.CellVol..mol.O2.L._lowFe <- phia_Gross.P.CellVol..mol.O2.L._lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "44_Gross.P.CellVol..mol.O2.L.",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_Gross.P.CellVol..mol.O2.L._lowFe), DF_here = phia_Gross.P.CellVol..mol.O2.L._lowFe[5,2]) %>% 
     slice (3:4)  





#Making Gross.P.CellVol..mol.O2.L. data.frame
kable(( Gross.P.CellVol..mol.O2.L. <- bind_rows(phia_Gross.P.CellVol..mol.O2.L._lowCu, phia_Gross.P.CellVol..mol.O2.L._lowFe)), format = "markdown")
```



|      Value| Df| Sum of Sq|         F|    Pr(>F)|Analysis         |Data                          |test_intrctn     |tested         | DF_here|
|----------:|--:|---------:|---------:|---------:|:----------------|:-----------------------------|:----------------|:--------------|-------:|
| -0.2333333|  1| 0.0816667| 0.6680951| 0.8548166|phia_interaction |44_Gross.P.CellVol..mol.O2.L. |lowFeCu vs lowCu |TO 1003 :  low |      14|
|  0.1733333|  1| 0.0450667| 0.3686794| 0.8548166|phia_interaction |44_Gross.P.CellVol..mol.O2.L. |lowFeCu vs lowCu |TO 1005 :  low |      14|
|  0.3633333|  1| 0.1584133| 1.2959408| 0.8221874|phia_interaction |44_Gross.P.CellVol..mol.O2.L. |lowFeCu vs lowFe |TO 1003 :  low |      14|
| -0.2833333|  1| 0.0963333| 0.7880795| 0.8221874|phia_interaction |44_Gross.P.CellVol..mol.O2.L. |lowFeCu vs lowFe |TO 1005 :  low |      14|


<a id="pg.prot.cell"></a>

#### pg.prot.cell lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_pg.prot.cell <- lm(data=mydata, pg.prot.cell~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_pg.prot.cell_lowCu <- testInteractions(lm_pg.prot.cell, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_pg.prot.cell_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high -10.6033  1    168.65  9.5753 0.0158403 *  
## TO 1005 : high  23.1783  1    644.68 36.6031 0.0001195 ***
## TO 1003 :  low -16.0167  1    384.80 21.8479 0.0010749 ** 
## TO 1005 :  low   7.9855  1     76.52  4.3446 0.0559242 .  
## Residuals               14    246.58                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_pg.prot.cell_lowCu <- phia_pg.prot.cell_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "35_pg.prot.cell",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_pg.prot.cell_lowCu), DF_here = phia_pg.prot.cell_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_pg.prot.cell_lowFe <- testInteractions(lm_pg.prot.cell, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_pg.prot.cell_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high -2.8300  1    12.013  0.6821 0.845443   
## TO 1005 : high 15.7900  1   249.324 14.1559 0.008404 **
## TO 1003 :  low -8.2433  1   101.929  5.7872 0.091607 . 
## TO 1005 :  low  0.5972  1     0.535  0.0304 0.864147   
## Residuals              14   246.578                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_pg.prot.cell_lowFe <- phia_pg.prot.cell_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "35_pg.prot.cell",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_pg.prot.cell_lowFe), DF_here = phia_pg.prot.cell_lowFe[5,2]) %>% 
     slice (3:4)  





#Making pg.prot.cell data.frame
kable(( pg.prot.cell <- bind_rows(phia_pg.prot.cell_lowCu, phia_pg.prot.cell_lowFe)), format = "markdown")
```



|       Value| Df|   Sum of Sq|          F|    Pr(>F)|Analysis         |Data            |test_intrctn     |tested         | DF_here|
|-----------:|--:|-----------:|----------:|---------:|:----------------|:---------------|:----------------|:--------------|-------:|
| -16.0166667|  1| 384.8004167| 21.8478978| 0.0010749|phia_interaction |35_pg.prot.cell |lowFeCu vs lowCu |TO 1003 :  low |      14|
|   7.9854555|  1|  76.5209990|  4.3446496| 0.0559242|phia_interaction |35_pg.prot.cell |lowFeCu vs lowCu |TO 1005 :  low |      14|
|  -8.2433333|  1| 101.9288167|  5.7872348| 0.0916068|phia_interaction |35_pg.prot.cell |lowFeCu vs lowFe |TO 1003 :  low |      14|
|   0.5971628|  1|   0.5349051|  0.0303704| 0.8641474|phia_interaction |35_pg.prot.cell |lowFeCu vs lowFe |TO 1005 :  low |      14|




<a id="pg.prot.fl.cell"></a>

#### pg.prot.fl.cell lowFeCu like lowFe or lowCu ?

[Back Up](#BackUP)


```r
lm_pg.prot.fl.cell <- lm(data=mydata, pg.prot.fl.cell~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowFeCu vs lowCu
phia_pg.prot.fl.cell_lowCu <- testInteractions(lm_pg.prot.fl.cell, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowCu vs control
phia_pg.prot.fl.cell_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high -0.12000  1  0.021600  5.4225 0.0353735 *  
## TO 1005 : high  0.22671  1  0.061678 15.4838 0.0044860 ** 
## TO 1003 :  low -0.29667  1  0.132017 33.1419 0.0001986 ***
## TO 1005 :  low  0.15610  1  0.029240  7.3404 0.0338885 *  
## Residuals               14  0.055767                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_pg.prot.fl.cell_lowCu <- phia_pg.prot.fl.cell_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "36_pg.prot.fl.cell",test_intrctn = "lowFeCu vs lowCu", tested = rownames(phia_pg.prot.fl.cell_lowCu), DF_here = phia_pg.prot.fl.cell_lowCu[5,2]) %>% 
     slice (3:4)  



#lowFeCu vs lowFe
phia_pg.prot.fl.cell_lowFe <- testInteractions(lm_pg.prot.fl.cell, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowFe vs control
phia_pg.prot.fl.cell_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                    Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high -0.113333  1  0.019267  4.8368 0.1354870    
## TO 1005 : high  0.050000  1  0.002500  0.6276 0.8829032    
## TO 1003 :  low -0.290000  1  0.126150 31.6692 0.0002494 ***
## TO 1005 :  low -0.020614  1  0.000637  0.1600 0.8829032    
## Residuals                14  0.055767                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_pg.prot.fl.cell_lowFe <- phia_pg.prot.fl.cell_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "36_pg.prot.fl.cell",test_intrctn = "lowFeCu vs lowFe", tested = rownames(phia_pg.prot.fl.cell_lowFe), DF_here = phia_pg.prot.fl.cell_lowFe[5,2]) %>% 
     slice (3:4)  





#Making pg.prot.fl.cell data.frame
kable(( pg.prot.fl.cell <- bind_rows(phia_pg.prot.fl.cell_lowCu, phia_pg.prot.fl.cell_lowFe)), format = "markdown")
```



|      Value| Df| Sum of Sq|          F|    Pr(>F)|Analysis         |Data               |test_intrctn     |tested         | DF_here|
|----------:|--:|---------:|----------:|---------:|:----------------|:------------------|:----------------|:--------------|-------:|
| -0.2966667|  1| 0.1320167| 33.1419434| 0.0001986|phia_interaction |36_pg.prot.fl.cell |lowFeCu vs lowCu |TO 1003 :  low |      14|
|  0.1560975|  1| 0.0292397|  7.3404479| 0.0338885|phia_interaction |36_pg.prot.fl.cell |lowFeCu vs lowCu |TO 1005 :  low |      14|
| -0.2900000|  1| 0.1261500| 31.6691541| 0.0002494|phia_interaction |36_pg.prot.fl.cell |lowFeCu vs lowFe |TO 1003 :  low |      14|
| -0.0206142|  1| 0.0006374|  0.1600202| 0.8829032|phia_interaction |36_pg.prot.fl.cell |lowFeCu vs lowFe |TO 1005 :  low |      14|








<a id="Comb_Table_Phys"></a>

# Combining all Tables - physiology without FRRF

[Back Up](#BackUP)


```r
physiology <- bind_rows(Growthrate.dd, Growthrate.Percent..u.umax., Growthrate.specific.d.1, cell.SA.um2, cell.SA.Vol.ratio, cell.size..um, cell.volume.fl.cell, Chla.per.cell.pg.cell, Chla.per.cell.vol.fg.fL, FvFm.old, Sig.old, PQ_Siz.old, X14C.per.Chla.alpha, X14C.per.Chla.at.155uE, X14C.per.Chla.ek, X14C.per.Chla.pmax, FeDFB.zmol.cell.1.h., FeDFB.zmol.um.2.h., Gross.P.CellVol..mol.O2.L., GrossPcell..mol.O2.cell..h., GrossPchla.mol.O2.mol.Chla..h., Resp.mol.O2.mol.Chla.h, Resp.umol.O2.cell.h, AOXactivity, pg.prot.cell, pg.prot.fl.cell)


physiology <- physiology %>% 
  mutate(order = 1:104) %>% 
  unite (col = Comb, tested, test_intrctn, remove=FALSE) %>% 
  select ( Comb, Data, starts_with("Pr"))

physiology_wide <- physiology %>% 
  spread(Comb, 'Pr(>F)')

write.table(physiology_wide, "Output_Data/STATS_Comparison_Table_Physiology_lowFeCu_vs_lowFe_or_lowCu.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
```

<a id="FRRF"></a>

#FRRF DATA

[Back Up](#BackUP)

<a id="sig."></a>

#### sig. 

[Back Up](#BackUP)


```r
lm_sig. <- lm(data=FRRF, sig.~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowCu vs control
phia_sig._lowCu <- testInteractions(lm_sig., fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowCu vs control
phia_sig._lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F Pr(>F)  
## TO 1003 : high -80.956  1    9830.7 10.0167 0.0360 *
## TO 1005 : high   9.390  1      88.2  0.0898 1.0000  
## TO 1003 :  low  49.185  1    1814.4  1.8487 0.6035  
## TO 1005 :  low  -2.412  1       7.0  0.0071 1.0000  
## Residuals              11   10795.7                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_sig._lowCu <- phia_sig._lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "sig.",test_intrctn = "lowCu", tested = rownames(phia_sig._lowCu), DF_here = phia_sig._lowCu[5,2]) %>% 
     slice (1:2)  



#lowFe vs control
phia_sig._lowFe <- testInteractions(lm_sig., fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowFe vs control
phia_sig._lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high  -15.885  1       189  0.1928 0.6690686    
## TO 1005 : high -175.915  1     30946 31.5316 0.0004701 ***
## TO 1003 :  low  114.256  1     19581 19.9520 0.0019039 ** 
## TO 1005 :  low -187.717  1     42285 43.0851 0.0001623 ***
## Residuals               11     10796                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_sig._lowFe <- phia_sig._lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "sig.",test_intrctn = "lowFe", tested = rownames(phia_sig._lowFe), DF_here = phia_sig._lowFe[5,2]) %>% 
     slice (1:2)  



#lowFeCu vs control
#TO03
TO03_lowFeCu_linearHypothesis <-  linearHypothesis(lm_sig., "Fe.levellow+Cu.levellow+Fe.levellow:Cu.levellow")

TO03_lowFeCu_linearHypothesis <- TO03_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "sig.",test_intrctn = "lowFeCu", tested = "TO 1003 : high", DF_here = TO03_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)

#TO05
TO05_lowFeCu_linearHypothesis <-  linearHypothesis(lm_sig., c("Fe.levellow + Cu.levellow + SpeciesTO 1005:Fe.levellow + SpeciesTO 1005:Cu.levellow +
  Fe.levellow:Cu.levellow + SpeciesTO 1005:Fe.levellow:Cu.levellow"))

TO05_lowFeCu_linearHypothesis <- TO05_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "sig.",test_intrctn = "lowFeCu", tested = "TO 1005 : high", DF_here = TO05_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)


#Strain Comparison: TO03 vs TO05 
Strain_Comparison <- testInteractions(lm_sig., fixed=c("Cu.level", "Fe.level"), across="Species")
Strain_Comparison <- Strain_Comparison%>% 
  mutate(Analysis = "phia_interaction", Data = "sig.",test_intrctn = c("control", "lowCu", "lowFe", "lowFeCu", "DF"), tested = "TO03 vs TO05", DF_here = Strain_Comparison[5,2]) %>% 
  slice(1:4)


#Making growthrate data.frame
kable((sig. <- bind_rows(phia_sig._lowCu, phia_sig._lowFe, TO03_lowFeCu_linearHypothesis, TO05_lowFeCu_linearHypothesis, Strain_Comparison)), format = "markdown")
```



|      Value| Df|   Sum of Sq|          F|    Pr(>F)|Analysis         |Data |test_intrctn |tested         | DF_here| Res.Df|      RSS|
|----------:|--:|-----------:|----------:|---------:|:----------------|:----|:------------|:--------------|-------:|------:|--------:|
|  -80.95556|  1|  9830.70297| 10.0167065| 0.0359989|phia_interaction |sig. |lowCu        |TO 1003 : high |      11|     NA|       NA|
|    9.39000|  1|    88.17210|  0.0898404| 1.0000000|phia_interaction |sig. |lowCu        |TO 1005 : high |      11|     NA|       NA|
|  -15.88519|  1|   189.25433|  0.1928352| 0.6690686|phia_interaction |sig. |lowFe        |TO 1003 : high |      11|     NA|       NA|
| -175.91500|  1| 30946.08722| 31.5316082| 0.0004701|phia_interaction |sig. |lowFe        |TO 1005 : high |      11|     NA|       NA|
|         NA|  1|  1663.33500|  1.6948064| 0.2195707|phia_interaction |sig. |lowFeCu      |TO 1003 : high |      11|     11| 10795.74|
|         NA|  1| 38160.48005| 38.8825022| 0.0000639|phia_interaction |sig. |lowFeCu      |TO 1005 : high |      11|     11| 10795.74|
|   -1.91000|  1|     4.37772|  0.0044605| 0.9479495|phia_interaction |sig. |control      |TO03 vs TO05   |      11|     NA|       NA|
|   88.43556|  1|  9385.01699|  9.5625879| 0.0204863|phia_interaction |sig. |lowCu        |TO03 vs TO05   |      11|     NA|       NA|
| -161.93981|  1| 17483.00241| 17.8137927| 0.0043045|phia_interaction |sig. |lowFe        |TO03 vs TO05   |      11|     NA|       NA|
| -213.53667|  1| 68396.86202| 69.6909771| 0.0000174|phia_interaction |sig. |lowFeCu      |TO03 vs TO05   |      11|     NA|       NA|


<a id="F.q.F.v"></a>

#### F.q.F.v 

[Back Up](#BackUP)


```r
lm_F.q.F.v <- lm(data=FRRF, F.q.F.v~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowCu vs control
phia_F.q.F.v_lowCu <- testInteractions(lm_F.q.F.v, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowCu vs control
phia_F.q.F.v_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                    Value Df Sum of Sq       F  Pr(>F)  
## TO 1003 : high  0.123761  1 0.0229750 10.4816 0.03163 *
## TO 1005 : high  0.029303  1 0.0008587  0.3917 1.00000  
## TO 1003 :  low -0.112572  1 0.0095044  4.3361 0.18435  
## TO 1005 :  low -0.002199  1 0.0000058  0.0026 1.00000  
## Residuals                11 0.0241113                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_F.q.F.v_lowCu <- phia_F.q.F.v_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "F.q.F.v",test_intrctn = "lowCu", tested = rownames(phia_F.q.F.v_lowCu), DF_here = phia_F.q.F.v_lowCu[5,2]) %>% 
     slice (1:2)  



#lowFe vs control
phia_F.q.F.v_lowFe <- testInteractions(lm_F.q.F.v, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowFe vs control
phia_F.q.F.v_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                    Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high  0.005594  1  0.000023  0.0107 0.9194469    
## TO 1005 : high  0.129216  1  0.016697  7.6174 0.0556717 .  
## TO 1003 :  low -0.230739  1  0.079861 36.4338 0.0003391 ***
## TO 1005 :  low  0.097714  1  0.011458  5.2271 0.0861151 .  
## Residuals                11  0.024111                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_F.q.F.v_lowFe <- phia_F.q.F.v_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "F.q.F.v",test_intrctn = "lowFe", tested = rownames(phia_F.q.F.v_lowFe), DF_here = phia_F.q.F.v_lowFe[5,2]) %>% 
     slice (1:2)  



#lowFeCu vs control
#TO03
TO03_lowFeCu_linearHypothesis <-  linearHypothesis(lm_F.q.F.v, "Fe.levellow+Cu.levellow+Fe.levellow:Cu.levellow")

TO03_lowFeCu_linearHypothesis <- TO03_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "F.q.F.v",test_intrctn = "lowFeCu", tested = "TO 1003 : high", DF_here = TO03_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)

#TO05
TO05_lowFeCu_linearHypothesis <-  linearHypothesis(lm_F.q.F.v, c("Fe.levellow + Cu.levellow + SpeciesTO 1005:Fe.levellow + SpeciesTO 1005:Cu.levellow +
  Fe.levellow:Cu.levellow + SpeciesTO 1005:Fe.levellow:Cu.levellow"))

TO05_lowFeCu_linearHypothesis <- TO05_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "F.q.F.v",test_intrctn = "lowFeCu", tested = "TO 1005 : high", DF_here = TO05_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)


#Strain Comparison: TO03 vs TO05 
Strain_Comparison <- testInteractions(lm_F.q.F.v, fixed=c("Cu.level", "Fe.level"), across="Species")
Strain_Comparison <- Strain_Comparison%>% 
  mutate(Analysis = "phia_interaction", Data = "F.q.F.v",test_intrctn = c("control", "lowCu", "lowFe", "lowFeCu", "DF"), tested = "TO03 vs TO05", DF_here = Strain_Comparison[5,2]) %>% 
  slice(1:4)


#Making growthrate data.frame
kable((F.q.F.v <- bind_rows(phia_F.q.F.v_lowCu, phia_F.q.F.v_lowFe, TO03_lowFeCu_linearHypothesis, TO05_lowFeCu_linearHypothesis, Strain_Comparison)), format = "markdown")
```



|      Value| Df| Sum of Sq|          F|    Pr(>F)|Analysis         |Data    |test_intrctn |tested         | DF_here| Res.Df|       RSS|
|----------:|--:|---------:|----------:|---------:|:----------------|:-------|:------------|:--------------|-------:|------:|---------:|
|  0.1237606|  1| 0.0229750| 10.4815893| 0.0316346|phia_interaction |F.q.F.v |lowCu        |TO 1003 : high |      11|     NA|        NA|
|  0.0293034|  1| 0.0008587|  0.3917484| 1.0000000|phia_interaction |F.q.F.v |lowCu        |TO 1005 : high |      11|     NA|        NA|
|  0.0055941|  1| 0.0000235|  0.0107075| 0.9194469|phia_interaction |F.q.F.v |lowFe        |TO 1003 : high |      11|     NA|        NA|
|  0.1292162|  1| 0.0166968|  7.6173761| 0.0556717|phia_interaction |F.q.F.v |lowFe        |TO 1005 : high |      11|     NA|        NA|
|         NA|  1| 0.0171666|  7.8316723| 0.0173208|phia_interaction |F.q.F.v |lowFeCu      |TO 1003 : high |      11|     11| 0.0241113|
|         NA|  1| 0.0193600|  8.8323510| 0.0127015|phia_interaction |F.q.F.v |lowFeCu      |TO 1005 : high |      11|     11| 0.0241113|
|  0.0136439|  1| 0.0002234|  0.1019130| 0.7555264|phia_interaction |F.q.F.v |control      |TO03 vs TO05   |      11|     NA|        NA|
| -0.0808133|  1| 0.0078369|  3.5753447| 0.1705172|phia_interaction |F.q.F.v |lowCu        |TO03 vs TO05   |      11|     NA|        NA|
|  0.1372661|  1| 0.0125613|  5.7306827| 0.1068475|phia_interaction |F.q.F.v |lowFe        |TO03 vs TO05   |      11|     NA|        NA|
|  0.2476393|  1| 0.0919878| 41.9663908| 0.0001825|phia_interaction |F.q.F.v |lowFeCu      |TO03 vs TO05   |      11|     NA|        NA|

<a id="F.v.F.m"></a>

#### F.v.F.m 

[Back Up](#BackUP)


```r
lm_F.v.F.m <- lm(data=FRRF, F.v.F.m~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowCu vs control
phia_F.v.F.m_lowCu <- testInteractions(lm_F.v.F.m, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowCu vs control
phia_F.v.F.m_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq        F   Pr(>F)    
## TO 1003 : high  0.32118  1  0.154733 137.5812 5.88e-07 ***
## TO 1005 : high -0.02573  1  0.000662   0.5887 0.918184    
## TO 1003 :  low  0.16563  1  0.020575  18.2945 0.003914 ** 
## TO 1005 :  low  0.00493  1  0.000029   0.0259 0.918184    
## Residuals               11  0.012371                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_F.v.F.m_lowCu <- phia_F.v.F.m_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "F.v.F.m",test_intrctn = "lowCu", tested = rownames(phia_F.v.F.m_lowCu), DF_here = phia_F.v.F.m_lowCu[5,2]) %>% 
     slice (1:2)  



#lowFe vs control
phia_F.v.F.m_lowFe <- testInteractions(lm_F.v.F.m, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowFe vs control
phia_F.v.F.m_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq      F    Pr(>F)    
## TO 1003 : high 0.243404  1  0.044434 39.509 0.0002382 ***
## TO 1005 : high 0.129469  1  0.016762 14.904 0.0053012 ** 
## TO 1003 :  low 0.087857  1  0.011578 10.295 0.0083267 ** 
## TO 1005 :  low 0.160124  1  0.030768 27.357 0.0008428 ***
## Residuals               11  0.012371                     
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_F.v.F.m_lowFe <- phia_F.v.F.m_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "F.v.F.m",test_intrctn = "lowFe", tested = rownames(phia_F.v.F.m_lowFe), DF_here = phia_F.v.F.m_lowFe[5,2]) %>% 
     slice (1:2)  



#lowFeCu vs control
#TO03
TO03_lowFeCu_linearHypothesis <-  linearHypothesis(lm_F.v.F.m, "Fe.levellow+Cu.levellow+Fe.levellow:Cu.levellow")

TO03_lowFeCu_linearHypothesis <- TO03_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "F.v.F.m",test_intrctn = "lowFeCu", tested = "TO 1003 : high", DF_here = TO03_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)

#TO05
TO05_lowFeCu_linearHypothesis <-  linearHypothesis(lm_F.v.F.m, c("Fe.levellow + Cu.levellow + SpeciesTO 1005:Fe.levellow + SpeciesTO 1005:Cu.levellow +
  Fe.levellow:Cu.levellow + SpeciesTO 1005:Fe.levellow:Cu.levellow"))

TO05_lowFeCu_linearHypothesis <- TO05_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "F.v.F.m",test_intrctn = "lowFeCu", tested = "TO 1005 : high", DF_here = TO05_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)


#Strain Comparison: TO03 vs TO05 
Strain_Comparison <- testInteractions(lm_F.v.F.m, fixed=c("Cu.level", "Fe.level"), across="Species")
Strain_Comparison <- Strain_Comparison%>% 
  mutate(Analysis = "phia_interaction", Data = "F.v.F.m",test_intrctn = c("control", "lowCu", "lowFe", "lowFeCu", "DF"), tested = "TO03 vs TO05", DF_here = Strain_Comparison[5,2]) %>% 
  slice(1:4)


#Making growthrate data.frame
kable((F.v.F.m <- bind_rows(phia_F.v.F.m_lowCu, phia_F.v.F.m_lowFe, TO03_lowFeCu_linearHypothesis, TO05_lowFeCu_linearHypothesis, Strain_Comparison)), format = "markdown")
```



|      Value| Df| Sum of Sq|           F|    Pr(>F)|Analysis         |Data    |test_intrctn |tested         | DF_here| Res.Df|       RSS|
|----------:|--:|---------:|-----------:|---------:|:----------------|:-------|:------------|:--------------|-------:|------:|---------:|
|  0.3211778|  1| 0.1547328| 137.5812459| 0.0000006|phia_interaction |F.v.F.m |lowCu        |TO 1003 : high |      11|     NA|        NA|
| -0.0257301|  1| 0.0006620|   0.5886528| 0.9181842|phia_interaction |F.v.F.m |lowCu        |TO 1005 : high |      11|     NA|        NA|
|  0.2434043|  1| 0.0444343|  39.5088774| 0.0002382|phia_interaction |F.v.F.m |lowFe        |TO 1003 : high |      11|     NA|        NA|
|  0.1294689|  1| 0.0167622|  14.9041665| 0.0053012|phia_interaction |F.v.F.m |lowFe        |TO 1005 : high |      11|     NA|        NA|
|         NA|  1| 0.2509646| 223.1461192| 0.0000000|phia_interaction |F.v.F.m |lowFeCu      |TO 1003 : high |      11|     11| 0.0123713|
|         NA|  1| 0.0216742|  19.2716644| 0.0010811|phia_interaction |F.v.F.m |lowFeCu      |TO 1005 : high |      11|     11| 0.0123713|
|  0.0111939|  1| 0.0001504|   0.1336972| 0.7215632|phia_interaction |F.v.F.m |control      |TO03 vs TO05   |      11|     NA|        NA|
| -0.3357140|  1| 0.1352447| 120.2532997| 0.0000012|phia_interaction |F.v.F.m |lowCu        |TO03 vs TO05   |      11|     NA|        NA|
| -0.1027415|  1| 0.0070372|   6.2571642| 0.0588627|phia_interaction |F.v.F.m |lowFe        |TO03 vs TO05   |      11|     NA|        NA|
| -0.2634470|  1| 0.1041065|  92.5666718| 0.0000033|phia_interaction |F.v.F.m |lowFeCu      |TO03 vs TO05   |      11|     NA|        NA|

<a id="F.q.F.m"></a>

#### F.q.F.m 

[Back Up](#BackUP)


```r
lm_F.q.F.m <- lm(data=FRRF, F.q.F.m~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowCu vs control
phia_F.q.F.m_lowCu <- testInteractions(lm_F.q.F.m, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowCu vs control
phia_F.q.F.m_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high 0.264513  1  0.104950 62.5928 2.904e-05 ***
## TO 1005 : high 0.000586  1  0.000000  0.0002    1.0000    
## TO 1003 :  low 0.090596  1  0.006156  3.6713    0.2451    
## TO 1005 :  low 0.005497  1  0.000036  0.0216    1.0000    
## Residuals               11  0.018444                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_F.q.F.m_lowCu <- phia_F.q.F.m_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "F.q.F.m",test_intrctn = "lowCu", tested = rownames(phia_F.q.F.m_lowCu), DF_here = phia_F.q.F.m_lowCu[5,2]) %>% 
     slice (1:2)  



#lowFe vs control
phia_F.q.F.m_lowFe <- testInteractions(lm_F.q.F.m, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowFe vs control
phia_F.q.F.m_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high 0.177555  1 0.0236442 14.1015 0.009543 **
## TO 1005 : high 0.146639  1 0.0215029 12.8244 0.009543 **
## TO 1003 :  low 0.003638  1 0.0000199  0.0118 0.915309   
## TO 1005 :  low 0.151549  1 0.0275605 16.4372 0.007607 **
## Residuals               11 0.0184439                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_F.q.F.m_lowFe <- phia_F.q.F.m_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "F.q.F.m",test_intrctn = "lowFe", tested = rownames(phia_F.q.F.m_lowFe), DF_here = phia_F.q.F.m_lowFe[5,2]) %>% 
     slice (1:2)  



#lowFeCu vs control
#TO03
TO03_lowFeCu_linearHypothesis <-  linearHypothesis(lm_F.q.F.m, "Fe.levellow+Cu.levellow+Fe.levellow:Cu.levellow")

TO03_lowFeCu_linearHypothesis <- TO03_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "F.q.F.m",test_intrctn = "lowFeCu", tested = "TO 1003 : high", DF_here = TO03_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)

#TO05
TO05_lowFeCu_linearHypothesis <-  linearHypothesis(lm_F.q.F.m, c("Fe.levellow + Cu.levellow + SpeciesTO 1005:Fe.levellow + SpeciesTO 1005:Cu.levellow +
  Fe.levellow:Cu.levellow + SpeciesTO 1005:Fe.levellow:Cu.levellow"))

TO05_lowFeCu_linearHypothesis <- TO05_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "F.q.F.m",test_intrctn = "lowFeCu", tested = "TO 1005 : high", DF_here = TO05_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)


#Strain Comparison: TO03 vs TO05 
Strain_Comparison <- testInteractions(lm_F.q.F.m, fixed=c("Cu.level", "Fe.level"), across="Species")
Strain_Comparison <- Strain_Comparison%>% 
  mutate(Analysis = "phia_interaction", Data = "F.q.F.m",test_intrctn = c("control", "lowCu", "lowFe", "lowFeCu", "DF"), tested = "TO03 vs TO05", DF_here = Strain_Comparison[5,2]) %>% 
  slice(1:4)


#Making growthrate data.frame
kable((F.q.F.m <- bind_rows(phia_F.q.F.m_lowCu, phia_F.q.F.m_lowFe, TO03_lowFeCu_linearHypothesis, TO05_lowFeCu_linearHypothesis, Strain_Comparison)), format = "markdown")
```



|      Value| Df| Sum of Sq|          F|    Pr(>F)|Analysis         |Data    |test_intrctn |tested         | DF_here| Res.Df|       RSS|
|----------:|--:|---------:|----------:|---------:|:----------------|:-------|:------------|:--------------|-------:|------:|---------:|
|  0.2645127|  1| 0.1049505| 62.5927693| 0.0000290|phia_interaction |F.q.F.m |lowCu        |TO 1003 : high |      11|     NA|        NA|
|  0.0005862|  1| 0.0000003|  0.0002050| 1.0000000|phia_interaction |F.q.F.m |lowCu        |TO 1005 : high |      11|     NA|        NA|
|  0.1775545|  1| 0.0236442| 14.1014719| 0.0095433|phia_interaction |F.q.F.m |lowFe        |TO 1003 : high |      11|     NA|        NA|
|  0.1466385|  1| 0.0215029| 12.8243679| 0.0095433|phia_interaction |F.q.F.m |lowFe        |TO 1005 : high |      11|     NA|        NA|
|         NA|  1| 0.1078573| 64.3263983| 0.0000064|phia_interaction |F.q.F.m |lowFeCu      |TO 1003 : high |      11|     11| 0.0184439|
|         NA|  1| 0.0277741| 16.5645706| 0.0018516|phia_interaction |F.q.F.m |lowFeCu      |TO 1005 : high |      11|     11| 0.0184439|
|  0.0168036|  1| 0.0003388|  0.2020810| 1.0000000|phia_interaction |F.q.F.m |control      |TO03 vs TO05   |      11|     NA|        NA|
| -0.2471229|  1| 0.0732837| 43.7066029| 0.0001522|phia_interaction |F.q.F.m |lowCu        |TO03 vs TO05   |      11|     NA|        NA|
| -0.0141124|  1| 0.0001328|  0.0791862| 1.0000000|phia_interaction |F.q.F.m |lowFe        |TO03 vs TO05   |      11|     NA|        NA|
| -0.0992121|  1| 0.0147645|  8.8056192| 0.0384124|phia_interaction |F.q.F.m |lowFeCu      |TO03 vs TO05   |      11|     NA|        NA|



<a id="ETR"></a>

#### ETR 

[Back Up](#BackUP)


```r
lm_ETR <- lm(data=FRRF, ETR~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowCu vs control
phia_ETR_lowCu <- testInteractions(lm_ETR, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowCu vs control
phia_ETR_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq      F Pr(>F)
## TO 1003 : high  24.6804  1     913.7 1.0678      1
## TO 1005 : high  22.4300  1     503.1 0.5880      1
## TO 1003 :  low -30.2287  1     685.3 0.8009      1
## TO 1005 :  low  -6.8233  1      55.9 0.0653      1
## Residuals               11    9412.1
```

```r
phia_ETR_lowCu <- phia_ETR_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR",test_intrctn = "lowCu", tested = rownames(phia_ETR_lowCu), DF_here = phia_ETR_lowCu[5,2]) %>% 
     slice (1:2)  



#lowFe vs control
phia_ETR_lowFe <- testInteractions(lm_ETR, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowFe vs control
phia_ETR_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq      F Pr(>F)
## TO 1003 : high  -2.774  1       5.8 0.0067 1.0000
## TO 1005 : high -18.780  1     352.7 0.4122 1.0000
## TO 1003 :  low -57.683  1    4991.0 5.8330 0.1372
## TO 1005 :  low -48.033  1    2768.6 3.2357 0.2985
## Residuals              11    9412.1
```

```r
phia_ETR_lowFe <- phia_ETR_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR",test_intrctn = "lowFe", tested = rownames(phia_ETR_lowFe), DF_here = phia_ETR_lowFe[5,2]) %>% 
     slice (1:2)  



#lowFeCu vs control
#TO03
TO03_lowFeCu_linearHypothesis <-  linearHypothesis(lm_ETR, "Fe.levellow+Cu.levellow+Fe.levellow:Cu.levellow")

TO03_lowFeCu_linearHypothesis <- TO03_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR",test_intrctn = "lowFeCu", tested = "TO 1003 : high", DF_here = TO03_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)

#TO05
TO05_lowFeCu_linearHypothesis <-  linearHypothesis(lm_ETR, c("Fe.levellow + Cu.levellow + SpeciesTO 1005:Fe.levellow + SpeciesTO 1005:Cu.levellow +
  Fe.levellow:Cu.levellow + SpeciesTO 1005:Fe.levellow:Cu.levellow"))

TO05_lowFeCu_linearHypothesis <- TO05_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR",test_intrctn = "lowFeCu", tested = "TO 1005 : high", DF_here = TO05_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)


#Strain Comparison: TO03 vs TO05 
Strain_Comparison <- testInteractions(lm_ETR, fixed=c("Cu.level", "Fe.level"), across="Species")
Strain_Comparison <- Strain_Comparison%>% 
  mutate(Analysis = "phia_interaction", Data = "ETR",test_intrctn = c("control", "lowCu", "lowFe", "lowFeCu", "DF"), tested = "TO03 vs TO05", DF_here = Strain_Comparison[5,2]) %>% 
  slice(1:4)


#Making growthrate data.frame
kable((ETR <- bind_rows(phia_ETR_lowCu, phia_ETR_lowFe, TO03_lowFeCu_linearHypothesis, TO05_lowFeCu_linearHypothesis, Strain_Comparison)), format = "markdown")
```



|      Value| Df|   Sum of Sq|         F|    Pr(>F)|Analysis         |Data |test_intrctn |tested         | DF_here| Res.Df|      RSS|
|----------:|--:|-----------:|---------:|---------:|:----------------|:----|:------------|:--------------|-------:|------:|--------:|
|  24.680438|  1|  913.686027| 1.0678297| 1.0000000|phia_interaction |ETR  |lowCu        |TO 1003 : high |      11|     NA|       NA|
|  22.430000|  1|  503.104900| 0.5879814| 1.0000000|phia_interaction |ETR  |lowCu        |TO 1005 : high |      11|     NA|       NA|
|  -2.774109|  1|    5.771762| 0.0067455| 1.0000000|phia_interaction |ETR  |lowFe        |TO 1003 : high |      11|     NA|       NA|
| -18.780000|  1|  352.688400| 0.4121888| 1.0000000|phia_interaction |ETR  |lowFe        |TO 1005 : high |      11|     NA|       NA|
|         NA|  1| 1633.773852| 1.9094002| 0.1944477|phia_interaction |ETR  |lowFeCu      |TO 1003 : high |      11|     11| 9412.125|
|         NA|  1|  786.636813| 0.9193466| 0.3582474|phia_interaction |ETR  |lowFeCu      |TO 1005 : high |      11|     11| 9412.125|
|   7.455299|  1|   66.697779| 0.0779500| 1.0000000|phia_interaction |ETR  |control      |TO03 vs TO05   |      11|     NA|       NA|
|   5.204861|  1|   32.508693| 0.0379931| 1.0000000|phia_interaction |ETR  |lowCu        |TO03 vs TO05   |      11|     NA|       NA|
|  -8.550592|  1|   48.741746| 0.0569647| 1.0000000|phia_interaction |ETR  |lowFe        |TO03 vs TO05   |      11|     NA|       NA|
|  14.854732|  1|  330.994579| 0.3868351| 1.0000000|phia_interaction |ETR  |lowFeCu      |TO03 vs TO05   |      11|     NA|       NA|

<a id="NPQ"></a>

#### NPQ 

[Back Up](#BackUP)


```r
lm_NPQ <- lm(data=FRRF, NPQ~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowCu vs control
phia_NPQ_lowCu <- testInteractions(lm_NPQ, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowCu vs control
phia_NPQ_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                    Value Df Sum of Sq       F    Pr(>F)    
## TO 1003 : high -0.021093  1  0.000667  0.4810 1.0000000    
## TO 1005 : high  0.030000  1  0.000900  0.6486 1.0000000    
## TO 1003 :  low -0.229647  1  0.039553 28.5035 0.0009513 ***
## TO 1005 :  low  0.021667  1  0.000563  0.4060 1.0000000    
## Residuals                11  0.015264                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_NPQ_lowCu <- phia_NPQ_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "NPQ",test_intrctn = "lowCu", tested = rownames(phia_NPQ_lowCu), DF_here = phia_NPQ_lowCu[5,2]) %>% 
     slice (1:2)  



#lowFe vs control
phia_NPQ_lowFe <- testInteractions(lm_NPQ, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowFe vs control
phia_NPQ_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq        F    Pr(>F)    
## TO 1003 : high -0.12537  1  0.011788   8.4950   0.04222 *  
## TO 1005 : high -0.01500  1  0.000225   0.1621   1.00000    
## TO 1003 :  low -0.33392  1  0.167257 120.5317 1.155e-06 ***
## TO 1005 :  low -0.02333  1  0.000653   0.4708   1.00000    
## Residuals               11  0.015264                       
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_NPQ_lowFe <- phia_NPQ_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "NPQ",test_intrctn = "lowFe", tested = rownames(phia_NPQ_lowFe), DF_here = phia_NPQ_lowFe[5,2]) %>% 
     slice (1:2)  



#lowFeCu vs control
#TO03
TO03_lowFeCu_linearHypothesis <-  linearHypothesis(lm_NPQ, "Fe.levellow+Cu.levellow+Fe.levellow:Cu.levellow")

TO03_lowFeCu_linearHypothesis <- TO03_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "NPQ",test_intrctn = "lowFeCu", tested = "TO 1003 : high", DF_here = TO03_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)

#TO05
TO05_lowFeCu_linearHypothesis <-  linearHypothesis(lm_NPQ, c("Fe.levellow + Cu.levellow + SpeciesTO 1005:Fe.levellow + SpeciesTO 1005:Cu.levellow +
  Fe.levellow:Cu.levellow + SpeciesTO 1005:Fe.levellow:Cu.levellow"))

TO05_lowFeCu_linearHypothesis <- TO05_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "NPQ",test_intrctn = "lowFeCu", tested = "TO 1005 : high", DF_here = TO05_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)


#Strain Comparison: TO03 vs TO05 
Strain_Comparison <- testInteractions(lm_NPQ, fixed=c("Cu.level", "Fe.level"), across="Species")
Strain_Comparison <- Strain_Comparison%>% 
  mutate(Analysis = "phia_interaction", Data = "NPQ",test_intrctn = c("control", "lowCu", "lowFe", "lowFeCu", "DF"), tested = "TO03 vs TO05", DF_here = Strain_Comparison[5,2]) %>% 
  slice(1:4)


#Making growthrate data.frame
kable((NPQ <- bind_rows(phia_NPQ_lowCu, phia_NPQ_lowFe, TO03_lowFeCu_linearHypothesis, TO05_lowFeCu_linearHypothesis, Strain_Comparison)), format = "markdown")
```



|      Value| Df| Sum of Sq|           F|    Pr(>F)|Analysis         |Data |test_intrctn |tested         | DF_here| Res.Df|       RSS|
|----------:|--:|---------:|-----------:|---------:|:----------------|:----|:------------|:--------------|-------:|------:|---------:|
| -0.0210934|  1| 0.0006674|   0.4809525| 1.0000000|phia_interaction |NPQ  |lowCu        |TO 1003 : high |      11|     NA|        NA|
|  0.0300000|  1| 0.0009000|   0.6485724| 1.0000000|phia_interaction |NPQ  |lowCu        |TO 1005 : high |      11|     NA|        NA|
| -0.1253701|  1| 0.0117883|   8.4950391| 0.0422164|phia_interaction |NPQ  |lowFe        |TO 1003 : high |      11|     NA|        NA|
| -0.0150000|  1| 0.0002250|   0.1621431| 1.0000000|phia_interaction |NPQ  |lowFe        |TO 1005 : high |      11|     NA|        NA|
|         NA|  1| 0.1890556| 136.2402383| 0.0000002|phia_interaction |NPQ  |lowFeCu      |TO 1003 : high |      11|     11| 0.0152643|
|         NA|  1| 0.0000533|   0.0384339| 0.8481491|phia_interaction |NPQ  |lowFeCu      |TO 1005 : high |      11|     11| 0.0152643|
| -0.0306939|  1| 0.0011305|   0.8147083| 0.7721409|phia_interaction |NPQ  |control      |TO03 vs TO05   |      11|     NA|        NA|
|  0.0203995|  1| 0.0004994|   0.3598627| 0.7721409|phia_interaction |NPQ  |lowCu        |TO03 vs TO05   |      11|     NA|        NA|
|  0.0796762|  1| 0.0042322|   3.0498745| 0.3257086|phia_interaction |NPQ  |lowFe        |TO03 vs TO05   |      11|     NA|        NA|
|  0.3309897|  1| 0.1643313| 118.4230236| 0.0000013|phia_interaction |NPQ  |lowFeCu      |TO03 vs TO05   |      11|     NA|        NA|

<a id="ETR.alpha.JP"></a>

#### ETR.alpha.JP 

[Back Up](#BackUP)


```r
lm_ETR.alpha.JP <- lm(data=FRRF, ETR.alpha.JP~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowCu vs control
phia_ETR.alpha.JP_lowCu <- testInteractions(lm_ETR.alpha.JP, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowCu vs control
phia_ETR.alpha.JP_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high -0.54801  1   0.45048 24.9128 0.001632 **
## TO 1005 : high -0.16807  1   0.02825  1.5622 0.711830   
## TO 1003 :  low  0.02820  1   0.00060  0.0330 1.000000   
## TO 1005 :  low  0.01488  1   0.00027  0.0147 1.000000   
## Residuals               11   0.19890                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_ETR.alpha.JP_lowCu <- phia_ETR.alpha.JP_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.alpha.JP",test_intrctn = "lowCu", tested = rownames(phia_ETR.alpha.JP_lowCu), DF_here = phia_ETR.alpha.JP_lowCu[5,2]) %>% 
     slice (1:2)  



#lowFe vs control
phia_ETR.alpha.JP_lowFe <- testInteractions(lm_ETR.alpha.JP, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowFe vs control
phia_ETR.alpha.JP_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high -0.25991  1   0.05066  2.8019 0.122319   
## TO 1005 : high -0.63607  1   0.40459 22.3751 0.002477 **
## TO 1003 :  low  0.31630  1   0.15007  8.2994 0.029896 * 
## TO 1005 :  low -0.45312  1   0.24639 13.6260 0.010665 * 
## Residuals               11   0.19890                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_ETR.alpha.JP_lowFe <- phia_ETR.alpha.JP_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.alpha.JP",test_intrctn = "lowFe", tested = rownames(phia_ETR.alpha.JP_lowFe), DF_here = phia_ETR.alpha.JP_lowFe[5,2]) %>% 
     slice (1:2)  



#lowFeCu vs control
#TO03
TO03_lowFeCu_linearHypothesis <-  linearHypothesis(lm_ETR.alpha.JP, "Fe.levellow+Cu.levellow+Fe.levellow:Cu.levellow")

TO03_lowFeCu_linearHypothesis <- TO03_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.alpha.JP",test_intrctn = "lowFeCu", tested = "TO 1003 : high", DF_here = TO03_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)

#TO05
TO05_lowFeCu_linearHypothesis <-  linearHypothesis(lm_ETR.alpha.JP, c("Fe.levellow + Cu.levellow + SpeciesTO 1005:Fe.levellow + SpeciesTO 1005:Cu.levellow +
  Fe.levellow:Cu.levellow + SpeciesTO 1005:Fe.levellow:Cu.levellow"))

TO05_lowFeCu_linearHypothesis <- TO05_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.alpha.JP",test_intrctn = "lowFeCu", tested = "TO 1005 : high", DF_here = TO05_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)


#Strain Comparison: TO03 vs TO05 
Strain_Comparison <- testInteractions(lm_ETR.alpha.JP, fixed=c("Cu.level", "Fe.level"), across="Species")
Strain_Comparison <- Strain_Comparison%>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.alpha.JP",test_intrctn = c("control", "lowCu", "lowFe", "lowFeCu", "DF"), tested = "TO03 vs TO05", DF_here = Strain_Comparison[5,2]) %>% 
  slice(1:4)


#Making growthrate data.frame
kable((ETR.alpha.JP <- bind_rows(phia_ETR.alpha.JP_lowCu, phia_ETR.alpha.JP_lowFe, TO03_lowFeCu_linearHypothesis, TO05_lowFeCu_linearHypothesis, Strain_Comparison)), format = "markdown")
```



|      Value| Df| Sum of Sq|         F|    Pr(>F)|Analysis         |Data         |test_intrctn |tested         | DF_here| Res.Df|       RSS|
|----------:|--:|---------:|---------:|---------:|:----------------|:------------|:------------|:--------------|-------:|------:|---------:|
| -0.5480130|  1| 0.4504773| 24.912839| 0.0016324|phia_interaction |ETR.alpha.JP |lowCu        |TO 1003 : high |      11|     NA|        NA|
| -0.1680730|  1| 0.0282485|  1.562234| 0.7118297|phia_interaction |ETR.alpha.JP |lowCu        |TO 1005 : high |      11|     NA|        NA|
| -0.2599080|  1| 0.0506641|  2.801889| 0.1223186|phia_interaction |ETR.alpha.JP |lowFe        |TO 1003 : high |      11|     NA|        NA|
| -0.6360734|  1| 0.4045894| 22.375090| 0.0024765|phia_interaction |ETR.alpha.JP |lowFe        |TO 1005 : high |      11|     NA|        NA|
|         NA|  1| 0.0805351|  4.453848| 0.0585348|phia_interaction |ETR.alpha.JP |lowFeCu      |TO 1003 : high |      11|     11| 0.1989035|
|         NA|  1| 0.4630639| 25.608916| 0.0003660|phia_interaction |ETR.alpha.JP |lowFeCu      |TO 1005 : high |      11|     11| 0.1989035|
|  0.1998192|  1| 0.0479132|  2.649756| 0.3360430|phia_interaction |ETR.alpha.JP |control      |TO03 vs TO05   |      11|     NA|        NA|
|  0.5797591|  1| 0.4033448| 22.306259| 0.0025058|phia_interaction |ETR.alpha.JP |lowCu        |TO03 vs TO05   |      11|     NA|        NA|
| -0.1763462|  1| 0.0207320|  1.146546| 0.3360430|phia_interaction |ETR.alpha.JP |lowFe        |TO03 vs TO05   |      11|     NA|        NA|
| -0.1896674|  1| 0.0539606|  2.984192| 0.3360430|phia_interaction |ETR.alpha.JP |lowFeCu      |TO03 vs TO05   |      11|     NA|        NA|


<a id="ETR.ek.JP"></a>

#### ETR.ek.JP 

[Back Up](#BackUP)


```r
lm_ETR.ek.JP <- lm(data=FRRF, ETR.ek.JP~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowCu vs control
phia_ETR.ek.JP_lowCu <- testInteractions(lm_ETR.ek.JP, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowCu vs control
phia_ETR.ek.JP_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F  Pr(>F)  
## TO 1003 : high  35.454  1   1885.52 10.7064 0.02975 *
## TO 1005 : high  27.947  1    781.04  4.4349 0.17697  
## TO 1003 :  low -22.742  1    387.89  2.2025 0.33173  
## TO 1005 :  low  -2.608  1      8.16  0.0464 0.83347  
## Residuals              11   1937.22                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_ETR.ek.JP_lowCu <- phia_ETR.ek.JP_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.ek.JP",test_intrctn = "lowCu", tested = rownames(phia_ETR.ek.JP_lowCu), DF_here = phia_ETR.ek.JP_lowCu[5,2]) %>% 
     slice (1:2)  



#lowFe vs control
phia_ETR.ek.JP_lowFe <- testInteractions(lm_ETR.ek.JP, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowFe vs control
phia_ETR.ek.JP_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq       F  Pr(>F)   
## TO 1003 : high   6.356  1      30.3  0.1721 1.00000   
## TO 1005 : high  27.401  1     750.8  4.2633 0.19002   
## TO 1003 :  low -51.840  1    4031.1 22.8893 0.00227 **
## TO 1005 :  low  -3.154  1      11.9  0.0678 1.00000   
## Residuals              11    1937.2                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_ETR.ek.JP_lowFe <- phia_ETR.ek.JP_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.ek.JP",test_intrctn = "lowFe", tested = rownames(phia_ETR.ek.JP_lowFe), DF_here = phia_ETR.ek.JP_lowFe[5,2]) %>% 
     slice (1:2)  



#lowFeCu vs control
#TO03
TO03_lowFeCu_linearHypothesis <-  linearHypothesis(lm_ETR.ek.JP, "Fe.levellow+Cu.levellow+Fe.levellow:Cu.levellow")

TO03_lowFeCu_linearHypothesis <- TO03_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.ek.JP",test_intrctn = "lowFeCu", tested = "TO 1003 : high", DF_here = TO03_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)

#TO05
TO05_lowFeCu_linearHypothesis <-  linearHypothesis(lm_ETR.ek.JP, c("Fe.levellow + Cu.levellow + SpeciesTO 1005:Fe.levellow + SpeciesTO 1005:Cu.levellow +
  Fe.levellow:Cu.levellow + SpeciesTO 1005:Fe.levellow:Cu.levellow"))

TO05_lowFeCu_linearHypothesis <- TO05_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.ek.JP",test_intrctn = "lowFeCu", tested = "TO 1005 : high", DF_here = TO05_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)


#Strain Comparison: TO03 vs TO05 
Strain_Comparison <- testInteractions(lm_ETR.ek.JP, fixed=c("Cu.level", "Fe.level"), across="Species")
Strain_Comparison <- Strain_Comparison%>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.ek.JP",test_intrctn = c("control", "lowCu", "lowFe", "lowFeCu", "DF"), tested = "TO03 vs TO05", DF_here = Strain_Comparison[5,2]) %>% 
  slice(1:4)


#Making growthrate data.frame
kable((ETR.ek.JP <- bind_rows(phia_ETR.ek.JP_lowCu, phia_ETR.ek.JP_lowFe, TO03_lowFeCu_linearHypothesis, TO05_lowFeCu_linearHypothesis, Strain_Comparison)), format = "markdown")
```



|      Value| Df|  Sum of Sq|          F|    Pr(>F)|Analysis         |Data      |test_intrctn |tested         | DF_here| Res.Df|      RSS|
|----------:|--:|----------:|----------:|---------:|:----------------|:---------|:------------|:--------------|-------:|------:|--------:|
|  35.454397|  1| 1885.52139| 10.7064217| 0.0297513|phia_interaction |ETR.ek.JP |lowCu        |TO 1003 : high |      11|     NA|       NA|
|  27.947134|  1|  781.04230|  4.4349368| 0.1769745|phia_interaction |ETR.ek.JP |lowCu        |TO 1005 : high |      11|     NA|       NA|
|   6.356288|  1|   30.30180|  0.1720606| 1.0000000|phia_interaction |ETR.ek.JP |lowFe        |TO 1003 : high |      11|     NA|       NA|
|  27.401147|  1|  750.82285|  4.2633438| 0.1900208|phia_interaction |ETR.ek.JP |lowFe        |TO 1005 : high |      11|     NA|       NA|
|         NA|  1|  402.72849|  2.2867845| 0.1586655|phia_interaction |ETR.ek.JP |lowFeCu      |TO 1003 : high |      11|     11| 1937.224|
|         NA|  1|  737.62387|  4.1883970| 0.0653621|phia_interaction |ETR.ek.JP |lowFeCu      |TO 1005 : high |      11|     11| 1937.224|
| -11.752807|  1|  165.75417|  0.9411901| 0.7056413|phia_interaction |ETR.ek.JP |control      |TO03 vs TO05   |      11|     NA|       NA|
| -19.260070|  1|  445.14036|  2.5276088| 0.4205308|phia_interaction |ETR.ek.JP |lowCu        |TO03 vs TO05   |      11|     NA|       NA|
|   9.292051|  1|   57.56148|  0.3268472| 0.7056413|phia_interaction |ETR.ek.JP |lowFe        |TO03 vs TO05   |      11|     NA|       NA|
|  29.425598|  1| 1298.79871|  7.3748762| 0.0803449|phia_interaction |ETR.ek.JP |lowFeCu      |TO03 vs TO05   |      11|     NA|       NA|

<a id="ETR.pmax.JP"></a>

#### ETR.pmax.JP 

[Back Up](#BackUP)


```r
lm_ETR.pmax.JP <- lm(data=FRRF, ETR.pmax.JP~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowCu vs control
phia_ETR.pmax.JP_lowCu <- testInteractions(lm_ETR.pmax.JP, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowCu vs control
phia_ETR.pmax.JP_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                  Value Df Sum of Sq      F Pr(>F)
## TO 1003 : high  68.051  1      6946 2.2377 0.6512
## TO 1005 : high  70.935  1      5032 1.6209 0.6596
## TO 1003 :  low -83.700  1      5254 1.6926 0.6596
## TO 1005 :  low  -8.079  1        78 0.0252 0.8767
## Residuals              11     34147
```

```r
phia_ETR.pmax.JP_lowCu <- phia_ETR.pmax.JP_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.pmax.JP",test_intrctn = "lowCu", tested = rownames(phia_ETR.pmax.JP_lowCu), DF_here = phia_ETR.pmax.JP_lowCu[5,2]) %>% 
     slice (1:2)  



#lowFe vs control
phia_ETR.pmax.JP_lowFe <- testInteractions(lm_ETR.pmax.JP, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowFe vs control
phia_ETR.pmax.JP_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F  Pr(>F)  
## TO 1003 : high   -7.825  1        46  0.0148 1.00000  
## TO 1005 : high   10.126  1       103  0.0330 1.00000  
## TO 1003 :  low -159.576  1     38197 12.3044 0.01961 *
## TO 1005 :  low  -68.888  1      5695  1.8344 0.60832  
## Residuals               11     34147                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_ETR.pmax.JP_lowFe <- phia_ETR.pmax.JP_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.pmax.JP",test_intrctn = "lowFe", tested = rownames(phia_ETR.pmax.JP_lowFe), DF_here = phia_ETR.pmax.JP_lowFe[5,2]) %>% 
     slice (1:2)  



#lowFeCu vs control
#TO03
TO03_lowFeCu_linearHypothesis <-  linearHypothesis(lm_ETR.pmax.JP, "Fe.levellow+Cu.levellow+Fe.levellow:Cu.levellow")

TO03_lowFeCu_linearHypothesis <- TO03_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.pmax.JP",test_intrctn = "lowFeCu", tested = "TO 1003 : high", DF_here = TO03_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)

#TO05
TO05_lowFeCu_linearHypothesis <-  linearHypothesis(lm_ETR.pmax.JP, c("Fe.levellow + Cu.levellow + SpeciesTO 1005:Fe.levellow + SpeciesTO 1005:Cu.levellow +
  Fe.levellow:Cu.levellow + SpeciesTO 1005:Fe.levellow:Cu.levellow"))

TO05_lowFeCu_linearHypothesis <- TO05_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.pmax.JP",test_intrctn = "lowFeCu", tested = "TO 1005 : high", DF_here = TO05_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)


#Strain Comparison: TO03 vs TO05 
Strain_Comparison <- testInteractions(lm_ETR.pmax.JP, fixed=c("Cu.level", "Fe.level"), across="Species")
Strain_Comparison <- Strain_Comparison%>% 
  mutate(Analysis = "phia_interaction", Data = "ETR.pmax.JP",test_intrctn = c("control", "lowCu", "lowFe", "lowFeCu", "DF"), tested = "TO03 vs TO05", DF_here = Strain_Comparison[5,2]) %>% 
  slice(1:4)


#Making growthrate data.frame
kable((ETR.pmax.JP <- bind_rows(phia_ETR.pmax.JP_lowCu, phia_ETR.pmax.JP_lowFe, TO03_lowFeCu_linearHypothesis, TO05_lowFeCu_linearHypothesis, Strain_Comparison)), format = "markdown")
```



|     Value| Df|    Sum of Sq|         F|    Pr(>F)|Analysis         |Data        |test_intrctn |tested         | DF_here| Res.Df|      RSS|
|---------:|--:|------------:|---------:|---------:|:----------------|:-----------|:------------|:--------------|-------:|------:|--------:|
| 68.051245|  1|  6946.457912| 2.2376774| 0.6512451|phia_interaction |ETR.pmax.JP |lowCu        |TO 1003 : high |      11|     NA|       NA|
| 70.934942|  1|  5031.766039| 1.6208936| 0.6595594|phia_interaction |ETR.pmax.JP |lowCu        |TO 1005 : high |      11|     NA|       NA|
| -7.824533|  1|    45.917483| 0.0147915| 1.0000000|phia_interaction |ETR.pmax.JP |lowFe        |TO 1003 : high |      11|     NA|       NA|
| 10.126267|  1|   102.541291| 0.0330318| 1.0000000|phia_interaction |ETR.pmax.JP |lowFe        |TO 1005 : high |      11|     NA|       NA|
|        NA|  1| 12565.195274| 4.0476534| 0.0693767|phia_interaction |ETR.pmax.JP |lowFeCu      |TO 1003 : high |      11|     11| 34147.48|
|        NA|  1|     5.029519| 0.0016202| 0.9686140|phia_interaction |ETR.pmax.JP |lowFeCu      |TO 1005 : high |      11|     11| 34147.48|
| -8.384718|  1|    84.364205| 0.0271764| 1.0000000|phia_interaction |ETR.pmax.JP |control      |TO03 vs TO05   |      11|     NA|       NA|
| -5.501021|  1|    36.313481| 0.0116977| 1.0000000|phia_interaction |ETR.pmax.JP |lowCu        |TO03 vs TO05   |      11|     NA|       NA|
|  9.566081|  1|    61.006610| 0.0196522| 1.0000000|phia_interaction |ETR.pmax.JP |lowFe        |TO03 vs TO05   |      11|     NA|       NA|
| 85.187382|  1| 10885.335206| 3.5065165| 0.3517289|phia_interaction |ETR.pmax.JP |lowFeCu      |TO03 vs TO05   |      11|     NA|       NA|


<a id="Converse_corr"></a>

#### Converse_corr 

[Back Up](#BackUP)


```r
lm_Converse_corr <- lm(data=FRRF, Converse_corr~(Species*Fe.level*Cu.level)) #this makes the statistical model

#lowCu vs control
phia_Converse_corr_lowCu <- testInteractions(lm_Converse_corr, fixed=c("Species", "Fe.level"), across="Cu.level") # in first two rows, this tests lowCu vs control
phia_Converse_corr_lowCu
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq       F   Pr(>F)   
## TO 1003 : high -214.608  1     55268  7.2562 0.050001 . 
## TO 1005 : high    5.567  1        46  0.0061 0.938761   
## TO 1003 :  low -301.585  1    136430 17.9119 0.002897 **
## TO 1005 :  low -105.927  1     16831  2.2097 0.315720   
## Residuals               15    114251                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
phia_Converse_corr_lowCu <- phia_Converse_corr_lowCu %>% 
  mutate(Analysis = "phia_interaction", Data = "Converse_corr",test_intrctn = "lowCu", tested = rownames(phia_Converse_corr_lowCu), DF_here = phia_Converse_corr_lowCu[5,2]) %>% 
     slice (1:2)  



#lowFe vs control
phia_Converse_corr_lowFe <- testInteractions(lm_Converse_corr, fixed=c("Species", "Cu.level"), across="Fe.level") # in first two rows, this tests lowFe vs control
phia_Converse_corr_lowFe
```

```
## F Test: 
## P-value adjustment method: holm
##                   Value Df Sum of Sq      F Pr(>F)
## TO 1003 : high  115.739  1     20093 2.6381 0.3852
## TO 1005 : high  -14.942  1       335 0.0440 1.0000
## TO 1003 :  low   28.763  1       993 0.1303 1.0000
## TO 1005 :  low -126.436  1     23979 3.1482 0.3852
## Residuals               15    114251
```

```r
phia_Converse_corr_lowFe <- phia_Converse_corr_lowFe %>% 
  mutate(Analysis = "phia_interaction", Data = "Converse_corr",test_intrctn = "lowFe", tested = rownames(phia_Converse_corr_lowFe), DF_here = phia_Converse_corr_lowFe[5,2]) %>% 
     slice (1:2)  



#lowFeCu vs control
#TO03
TO03_lowFeCu_linearHypothesis <-  linearHypothesis(lm_Converse_corr, "Fe.levellow+Cu.levellow+Fe.levellow:Cu.levellow")

TO03_lowFeCu_linearHypothesis <- TO03_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "Converse_corr",test_intrctn = "lowFeCu", tested = "TO 1003 : high", DF_here = TO03_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)

#TO05
TO05_lowFeCu_linearHypothesis <-  linearHypothesis(lm_Converse_corr, c("Fe.levellow + Cu.levellow + SpeciesTO 1005:Fe.levellow + SpeciesTO 1005:Cu.levellow +
  Fe.levellow:Cu.levellow + SpeciesTO 1005:Fe.levellow:Cu.levellow"))

TO05_lowFeCu_linearHypothesis <- TO05_lowFeCu_linearHypothesis %>% 
  mutate(Analysis = "phia_interaction", Data = "Converse_corr",test_intrctn = "lowFeCu", tested = "TO 1005 : high", DF_here = TO05_lowFeCu_linearHypothesis[2,1]) %>% 
  slice (2)


#Strain Comparison: TO03 vs TO05 
Strain_Comparison <- testInteractions(lm_Converse_corr, fixed=c("Cu.level", "Fe.level"), across="Species")
Strain_Comparison <- Strain_Comparison%>% 
  mutate(Analysis = "phia_interaction", Data = "Converse_corr",test_intrctn = c("control", "lowCu", "lowFe", "lowFeCu", "DF"), tested = "TO03 vs TO05", DF_here = Strain_Comparison[5,2]) %>% 
  slice(1:4)


#Making growthrate data.frame
kable((Converse_corr <- bind_rows(phia_Converse_corr_lowCu, phia_Converse_corr_lowFe, TO03_lowFeCu_linearHypothesis, TO05_lowFeCu_linearHypothesis, Strain_Comparison)), format = "markdown")
```



|      Value| Df|    Sum of Sq|          F|    Pr(>F)|Analysis         |Data          |test_intrctn |tested         | DF_here| Res.Df|      RSS|
|----------:|--:|------------:|----------:|---------:|:----------------|:-------------|:------------|:--------------|-------:|------:|--------:|
| -214.60825|  1|  55268.04259|  7.2561623| 0.0500013|phia_interaction |Converse_corr |lowCu        |TO 1003 : high |      15|     NA|       NA|
|    5.56714|  1|     46.48957|  0.0061036| 0.9387607|phia_interaction |Converse_corr |lowCu        |TO 1005 : high |      15|     NA|       NA|
|  115.73935|  1|  20093.39566|  2.6380696| 0.3852095|phia_interaction |Converse_corr |lowFe        |TO 1003 : high |      15|     NA|       NA|
|  -14.94246|  1|    334.91580|  0.0439712| 1.0000000|phia_interaction |Converse_corr |lowFe        |TO 1005 : high |      15|     NA|       NA|
|         NA|  1|  51807.74625|  6.8018587| 0.0197799|phia_interaction |Converse_corr |lowFeCu      |TO 1003 : high |      15|     15| 114250.6|
|         NA|  1|  21914.04844|  2.8771038| 0.1104961|phia_interaction |Converse_corr |lowFeCu      |TO 1005 : high |      15|     15| 114250.6|
|   74.68698|  1|   8367.21824|  1.0985353| 0.6223662|phia_interaction |Converse_corr |control      |TO03 vs TO05   |      15|     NA|       NA|
|  294.86238|  1| 104332.58515| 13.6978648| 0.0085378|phia_interaction |Converse_corr |lowCu        |TO03 vs TO05   |      15|     NA|       NA|
|  -55.99483|  1|   4703.13139|  0.6174759| 0.6223662|phia_interaction |Converse_corr |lowFe        |TO03 vs TO05   |      15|     NA|       NA|
|  139.66313|  1|  29258.68640|  3.8413841| 0.2065534|phia_interaction |Converse_corr |lowFeCu      |TO03 vs TO05   |      15|     NA|       NA|



<a id="Comb_Table_FRRF"></a>

# Combining all Tables - FRRF

[Back Up](#BackUP)


```r
frrf_table <- bind_rows(sig., F.q.F.v, F.v.F.m, F.q.F.m, ETR, ETR.alpha.JP, ETR.ek.JP, ETR.pmax.JP, NPQ, Converse_corr)


frrf_table <- frrf_table %>% 
  mutate(order = 1:100) %>% 
  unite (col = Comb, tested, test_intrctn, remove=FALSE) %>% 
  select ( Comb, Data, starts_with("Pr"))

kable((frrf_table_wide <- frrf_table %>% 
  spread(Comb, 'Pr(>F)')), format = "markdown")
```



|Data          | TO 1003 : high_lowCu| TO 1003 : high_lowFe| TO 1003 : high_lowFeCu| TO 1005 : high_lowCu| TO 1005 : high_lowFe| TO 1005 : high_lowFeCu| TO03 vs TO05_control| TO03 vs TO05_lowCu| TO03 vs TO05_lowFe| TO03 vs TO05_lowFeCu|
|:-------------|--------------------:|--------------------:|----------------------:|--------------------:|--------------------:|----------------------:|--------------------:|------------------:|------------------:|--------------------:|
|Converse_corr |            0.0500013|            0.3852095|              0.0197799|            0.9387607|            1.0000000|              0.1104961|            0.6223662|          0.0085378|          0.6223662|            0.2065534|
|ETR           |            1.0000000|            1.0000000|              0.1944477|            1.0000000|            1.0000000|              0.3582474|            1.0000000|          1.0000000|          1.0000000|            1.0000000|
|ETR.alpha.JP  |            0.0016324|            0.1223186|              0.0585348|            0.7118297|            0.0024765|              0.0003660|            0.3360430|          0.0025058|          0.3360430|            0.3360430|
|ETR.ek.JP     |            0.0297513|            1.0000000|              0.1586655|            0.1769745|            0.1900208|              0.0653621|            0.7056413|          0.4205308|          0.7056413|            0.0803449|
|ETR.pmax.JP   |            0.6512451|            1.0000000|              0.0693767|            0.6595594|            1.0000000|              0.9686140|            1.0000000|          1.0000000|          1.0000000|            0.3517289|
|F.q.F.m       |            0.0000290|            0.0095433|              0.0000064|            1.0000000|            0.0095433|              0.0018516|            1.0000000|          0.0001522|          1.0000000|            0.0384124|
|F.q.F.v       |            0.0316346|            0.9194469|              0.0173208|            1.0000000|            0.0556717|              0.0127015|            0.7555264|          0.1705172|          0.1068475|            0.0001825|
|F.v.F.m       |            0.0000006|            0.0002382|              0.0000000|            0.9181842|            0.0053012|              0.0010811|            0.7215632|          0.0000012|          0.0588627|            0.0000033|
|NPQ           |            1.0000000|            0.0422164|              0.0000002|            1.0000000|            1.0000000|              0.8481491|            0.7721409|          0.7721409|          0.3257086|            0.0000013|
|sig.          |            0.0359989|            0.6690686|              0.2195707|            1.0000000|            0.0004701|              0.0000639|            0.9479495|          0.0204863|          0.0043045|            0.0000174|

```r
write.table(frrf_table_wide, "Output_Data/STATS_Comparison_Table_FRRF.txt", sep = "\t", col.names = TRUE, row.names = FALSE)
```
