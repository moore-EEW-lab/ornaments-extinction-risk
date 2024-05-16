setwd('/Users/michaelmoore/desktop/Working Directory')


library(lme4)
library(lmerTest)
library(ggplot2)
library(phytools)
library(plyr)
library(metafor)
library(lsmeans)
library(visreg)
library(car)
library(nlme)
library(MuMIn)
library(ggplot2)
library(phylosignal)
library(phylobase)
library(viridis)
library(raster)
library(sp)
library(spgwr)


mw <- read.csv('assemblage.comp.csv')

## subset to grid cells that have 5 or more observations
mw.eff <- subset(mw, effort > 4)
dim(mw.eff) # 2181    9

## z transform so all of our variables are on a common scale. some variables bounded at zero and have zeros in there, so add small positive value to make ln transformation work
mw.eff$z.eff <- (log(mw.eff$effort + 0.001) - mean(log(mw.eff$effort + 0.001)))/sd(log(mw.eff$effort + 0.001))
mw.eff$z.imperv <- (log(mw.eff$impervious_percent + 0.001) - mean(log(mw.eff$impervious_percent + 0.001)))/sd(log(mw.eff$impervious_percent + 0.001))
mw.eff$z.hist.temp <- (log(mw.eff$avmaxt_61to65 + 0.001) - mean(log(mw.eff$avmaxt_61to65 + 0.001)))/sd(log(mw.eff$avmaxt_61to65 + 0.001))
mw.eff$z.warm <-(mw.eff$change_maxt - mean(mw.eff$change_maxt))/sd(mw.eff$change_maxt)
mw.eff$z.ag <- (log(mw.eff$agricul + 0.001) - mean(log(mw.eff$agricul + 0.001)))/sd(log(mw.eff$agricul + 0.001))

## naive model with no acconting for spatial autocorrelation
mod00 <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag, data = mw.eff)
summary(mod00)

### test what's the best correlation structure for accounting for spatial autocorrelation. Be warned, these models were taking 10-30 mins each on my computer
mod01a <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag, correlation = corExp(form = ~ Longitude + Latitude, nugget = TRUE), data = mw.eff)
mod01b <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag, correlation = corExp(value = 0.3018767 , form = ~ Longitude + Latitude, nugget = FALSE), data = mw.eff)
mod01c <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag, correlation = corGaus(form = ~ Longitude + Latitude, nugget = TRUE), data = mw.eff)
mod01d <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag, correlation = corGaus(form =~ Longitude + Latitude, nugget = FALSE), data = mw.eff)
mod01e <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag, correlation = corSpher(form =~ Longitude + Latitude, nugget = TRUE), data = mw.eff)
mod01f <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag, correlation = corSpher(form =~ Longitude + Latitude, nugget = FALSE), data = mw.eff)
mod01g <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag, correlation = corLin(form =~ Longitude + Latitude, nugget = TRUE), data = mw.eff)
mod01h <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag, correlation = corLin(form =~ Longitude + Latitude, nugget = FALSE), data = mw.eff)
mod01i <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag, correlation = corRatio(form =~ Longitude + Latitude, nugget = TRUE), data = mw.eff)
mod01j <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag, correlation = corRatio(form =~ Longitude + Latitude, nugget = FALSE), data = mw.eff)
AIC(mod00, mod01a, mod01b, mod01c, mod01d, mod01e, mod01f, mod01g, mod01h, mod01i, mod01j)

       # df       AIC
# mod00   7 -1536.403
# mod01a  9 -1695.050
# mod01b  8 -1649.137
# mod01c  9 -1688.482
# mod01d  8 -1596.453
# mod01e  9 -1690.906
# mod01f  8 -1587.710
# mod01g  9 -1578.759
# mod01h  8 -1580.759
# mod01i  9 -1694.705
# mod01j  8 -1667.106


## exponential model w/o a nugget is best. second best is rational quadratic model w/o nugget. 
## let's take a look at those models as well as the spatially uninformed one

# uninformed model first
summary(mod00)

# Generalized least squares fit by REML
  # Model: rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag 
  # Data: mw.eff 
        # AIC      BIC   logLik
  # -1536.403 -1496.61 775.2017

# Coefficients:
                 # Value   Std.Error   t-value p-value
# (Intercept)  0.6426127 0.003589658 179.01782  0.0000
# z.eff       -0.0222326 0.003699710  -6.00927  0.0000
# z.imperv     0.0230529 0.003725395   6.18805  0.0000
# z.hist.temp -0.0227259 0.003614146  -6.28803  0.0000
# z.warm      -0.0044414 0.003640923  -1.21984  0.2227
# z.ag         0.0113540 0.003691545   3.07567  0.0021

 # Correlation: 
            # (Intr) z.eff  z.mprv z.hst. z.warm
# z.eff        0.000                            
# z.imperv     0.000 -0.186                     
# z.hist.temp  0.000  0.087  0.053              
# z.warm       0.000  0.046  0.071 -0.011       
# z.ag         0.000  0.080  0.154  0.022  0.151

# Standardized residuals:
       # Min         Q1        Med         Q3        Max 
# -4.0631480 -0.5862187  0.0634300  0.6276560  2.8175238 

# Residual standard error: 0.1676412 
# Degrees of freedom: 2181 total; 2175 residual

### best model - exponential correlation structure w/o a nugget
summary(mod01a)

# Generalized least squares fit by REML
  # Model: rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag 
  # Data: mw.eff 
       # AIC       BIC   logLik
  # -1695.05 -1643.887 856.5251

# Correlation Structure: Exponential spatial correlation
 # Formula: ~Longitude + Latitude 
 # Parameter estimate(s):
    # range    nugget 
# 0.3044056 0.6881985 

# Coefficients:
                 # Value   Std.Error  t-value p-value
# (Intercept)  0.6535963 0.007397911 88.34877  0.0000
# z.eff       -0.0062047 0.004008948 -1.54772  0.1218
# z.imperv     0.0168648 0.004337143  3.88845  0.0001
# z.hist.temp -0.0242021 0.006310733 -3.83508  0.0001
# z.warm       0.0042357 0.006662272  0.63577  0.5250
# z.ag         0.0200595 0.004608926  4.35231  0.0000

 # Correlation: 
            # (Intr) z.eff  z.mprv z.hst. z.warm
# z.eff        0.134                            
# z.imperv     0.097 -0.158                     
# z.hist.temp -0.140  0.010  0.000              
# z.warm      -0.058  0.033  0.024  0.014       
# z.ag        -0.012  0.098  0.036  0.004  0.104

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -4.04589786 -0.67552648 -0.02954363  0.58067227  2.87750702 

# Residual standard error: 0.1689713 
# Degrees of freedom: 2181 total; 2175 residual



### second best model - rational quadratic correlation structure
summary(mod01i)

# Generalized least squares fit by REML
  # Model: rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag 
  # Data: mw.eff 
        # AIC       BIC   logLik
  # -1694.705 -1643.542 856.3523

# Correlation Structure: Rational quadratic spatial correlation
 # Formula: ~Longitude + Latitude 
 # Parameter estimate(s):
    # range    nugget 
# 0.2585759 0.7418707 

# Coefficients:
                 # Value   Std.Error  t-value p-value
# (Intercept)  0.6541781 0.008728940 74.94359  0.0000
# z.eff       -0.0061242 0.004015253 -1.52522  0.1273
# z.imperv     0.0166262 0.004351139  3.82111  0.0001
# z.hist.temp -0.0234575 0.006806027 -3.44658  0.0006
# z.warm       0.0040788 0.006941924  0.58756  0.5569
# z.ag         0.0197147 0.004624354  4.26323  0.0000

 # Correlation: 
            # (Intr) z.eff  z.mprv z.hst. z.warm
# z.eff        0.116                            
# z.imperv     0.089 -0.158                     
# z.hist.temp -0.143  0.007 -0.015              
# z.warm      -0.087  0.034  0.017  0.033       
# z.ag        -0.009  0.094  0.035  0.007  0.099

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -4.04689570 -0.68283331 -0.03153505  0.57674920  2.86288401 

# Residual standard error: 0.1689114 
# Degrees of freedom: 2181 total; 2175 residual


## set up a bunch of maximum likelihood models so we can do likelihood ratio tests. Since the exponential correlation structure model is best, let's use that

mod02a <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag, correlation = corExp(value = c(0.3044056, 0.6881985), form =~ Longitude + Latitude, nugget = TRUE), data = mw.eff, method = 'ML')
mod02b <- gls(rg_per_orn_spp ~ z.eff  + z.hist.temp + z.warm + z.ag, correlation = corExp(value = c(0.3044056, 0.6881985), form =~ Longitude + Latitude, nugget = TRUE), data = mw.eff, method = 'ML')
mod02c <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.warm + z.ag, correlation = corExp(value = c(0.3044056, 0.6881985), form =~ Longitude + Latitude, nugget = TRUE), data = mw.eff, method = 'ML')
mod02d <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.ag, correlation = corExp(value = c(0.3044056, 0.6881985), form =~ Longitude + Latitude, nugget = TRUE), data = mw.eff, method = 'ML')
mod02e <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm, correlation = corExp(value = c(0.3044056, 0.6881985), form =~ Longitude + Latitude, nugget = TRUE), data = mw.eff, method = 'ML')


anova(mod02a, mod02b) # test of effect of ISA
       # Model df       AIC       BIC   logLik   Test  L.Ratio p-value
# mod02a     1  9 -1746.861 -1695.673 882.4304                        
# mod02b     2  8 -1734.128 -1688.627 875.0638 1 vs 2 14.73326   1e-04

anova(mod02a, mod02c) # test of effect of historical climate
       # Model df       AIC       BIC   logLik   Test  L.Ratio p-value
# mod02a     1  9 -1746.861 -1695.673 882.4304                        
# mod02c     2  8 -1735.528 -1690.028 875.7643 1 vs 2 13.33225   3e-04

anova(mod02a, mod02d) # test of effect of warming
       # Model df       AIC       BIC   logLik   Test   L.Ratio p-value
# mod02a     1  9 -1746.861 -1695.673 882.4304                         
# mod02d     2  8 -1748.455 -1702.955 882.2276 1 vs 2 0.4055265  0.5242

anova(mod02a, mod02e) # test of effect of agriculture
       # Model df       AIC       BIC   logLik   Test  L.Ratio p-value
# mod02a     1  9 -1746.861 -1695.673 882.4304                        
# mod02e     2  8 -1729.942 -1684.441 872.9708 1 vs 2 18.91923  <.0001



##### let's look at the parameter estimates using the same model but this time only including grid cells w/ 7 or more observations. give it a couple starting values to speed up the computation time
mod01a.sub7 <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag, correlation = corExp(value = c(0.3044056, 0.6881985), form =~ Longitude + Latitude, nugget = TRUE), data = subset(mw.eff, effort >6))
summary(mod01a.sub7)

# Generalized least squares fit by REML
  # Model: rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag 
  # Data: subset(mw.eff, effort > 6) 
        # AIC       BIC  logLik
  # -1708.778 -1659.639 863.389

# Correlation Structure: Exponential spatial correlation
 # Formula: ~Longitude + Latitude 
 # Parameter estimate(s):
    # range    nugget 
# 0.3291419 0.6910116 

# Coefficients:
                 # Value   Std.Error  t-value p-value
# (Intercept)  0.6557998 0.007368350 89.00226  0.0000
# z.eff       -0.0115608 0.004184843 -2.76253  0.0058
# z.imperv     0.0163271 0.004328867  3.77168  0.0002
# z.hist.temp -0.0260805 0.006423659 -4.06008  0.0001
# z.warm      -0.0018320 0.006706265 -0.27318  0.7847
# z.ag         0.0146789 0.004540287  3.23303  0.0012

 # Correlation: 
            # (Intr) z.eff  z.mprv z.hst. z.warm
# z.eff       -0.005                            
# z.imperv     0.075 -0.137                     
# z.hist.temp -0.134  0.004 -0.010              
# z.warm      -0.046  0.024  0.031  0.024       
# z.ag         0.006  0.074  0.073 -0.005  0.107

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -4.47708509 -0.67486505 -0.03440255  0.58536598  2.94267983 

# Residual standard error: 0.1523419 
# Degrees of freedom: 1743 total; 1737 residual

### let's do the same analysis but this time including grid cells that have 3 or more observations. Also need to create the z-transformed variables again
mw.sub <- subset(mw, effort > 2 & !is.na(rg_per_orn_spp))
mw.sub$z.eff <- (log(mw.sub$effort + 0.001) - mean(log(mw.sub$effort + 0.001)))/sd(log(mw.sub$effort + 0.001))
mw.sub$z.imperv <- (log(mw.sub$impervious_percent + 0.001) - mean(log(mw.sub$impervious_percent + 0.001)))/sd(log(mw.sub$impervious_percent + 0.001))
mw.sub$z.hist.temp <- (log(mw.sub$avmaxt_61to65 + 0.001) - mean(log(mw.sub$avmaxt_61to65 + 0.001)))/sd(log(mw.sub$avmaxt_61to65 + 0.001))
mw.sub$z.warm <-(mw.sub$change_maxt - mean(mw.sub$change_maxt))/sd(mw.sub$change_maxt)
mw.sub$z.ag <- (log(mw.sub$agricul + 0.001) - mean(log(mw.sub$agricul + 0.001)))/sd(log(mw.sub$agricul + 0.001))


mod01a.sub3 <- gls(rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag, correlation = corExp(value = c(0.3044056, 0.6881985), form =~ Longitude + Latitude, nugget = TRUE), data = mw.sub)
summary(mod01a.sub3)
# Generalized least squares fit by REML
  # Model: rg_per_orn_spp ~ z.eff + z.imperv + z.hist.temp + z.warm + z.ag 
  # Data: mw.sub 
        # AIC       BIC   logLik
  # -973.3554 -919.1464 495.6777

# Correlation Structure: Exponential spatial correlation
 # Formula: ~Longitude + Latitude 
 # Parameter estimate(s):
    # range    nugget 
# 0.2889785 0.7618888 

# Coefficients:
                 # Value   Std.Error  t-value p-value
# (Intercept)  0.6590020 0.007369425 89.42380  0.0000
# z.eff       -0.0073483 0.004375450 -1.67945  0.0932
# z.imperv     0.0157432 0.004580225  3.43722  0.0006
# z.hist.temp -0.0288331 0.006535644 -4.41166  0.0000
# z.warm       0.0054446 0.006780674  0.80296  0.4221
# z.ag         0.0188250 0.004840186  3.88931  0.0001

 # Correlation: 
            # (Intr) z.eff  z.mprv z.hst. z.warm
# z.eff        0.121                            
# z.imperv     0.085 -0.204                     
# z.hist.temp -0.112  0.021 -0.002              
# z.warm      -0.034  0.033  0.014  0.031       
# z.ag        -0.014  0.141 -0.004  0.015  0.103

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -3.48240436 -0.63520903 -0.03164409  0.51971823  2.25987075 

# Residual standard error: 0.2101594 
# Degrees of freedom: 3057 total; 3051 residual


## the parameter estimates are pretty similar across the three different slices of the data, which is good.
