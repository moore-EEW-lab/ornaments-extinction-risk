setwd('/Users/michaelmoore/desktop/Working Directory')

library(car)
library(plyr)
library(lme4)
library(nlme)
library(phytools)
library(MuMIn)
library(ggplot2)
library(metafor)


pop.dat <- read.csv('pop.ext.csv')
phylo <- read.newick('pop.ext.tre')
rownames(pop.dat) <- pop.dat$binom


### meta-analytic models so that we can account for measurement error
occ.mod <- rma.uni(temporal.change ~ m.wing.color + z.size + z.MAT, vi = (pop.dat$SEM)^2, data = pop.dat, method = 'REML')
summary(occ.mod)

# Mixed-Effects Model (k = 19; tau^2 estimator: REML)

  # logLik  deviance       AIC       BIC      AICc 
# -46.2364   92.4728  102.4728  106.0131  109.1395   

# tau^2 (estimated amount of residual heterogeneity):     0.0972 (SE = 0.0756)
# tau (square root of estimated tau^2 value):             0.3117
# I^2 (residual heterogeneity / unaccounted variability): 57.95%
# H^2 (unaccounted variability / sampling variability):   2.38
# R^2 (amount of heterogeneity accounted for):            55.14%

# Test for Residual Heterogeneity:
# QE(df = 15) = 66.9820, p-val < .0001

# Test of Moderators (coefficients 2:4):
# QM(df = 3) = 13.2263, p-val = 0.0042

# Model Results:

               # estimate      se     zval    pval    ci.lb   ci.ub 
# intrcpt         -0.1367  0.1437  -0.9514  0.3414  -0.4184  0.1449      
# m.wing.colory    0.8797  0.2659   3.3085  0.0009   0.3586  1.4009  *** 
# z.size           0.2487  0.1021   2.4369  0.0148   0.0487  0.4487    * 
# z.MAT            0.2957  0.1602   1.8457  0.0649  -0.0183  0.6097    . 

# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## hypothesis tests
mod01a <- rma.uni(temporal.change ~ m.wing.color + z.size + z.MAT , vi = (pop.dat$SEM)^2, data = pop.dat, method = 'ML')
mod01b <- rma.uni(temporal.change ~ 1 + z.size + z.MAT , vi = (pop.dat$SEM)^2, data = pop.dat, method = 'ML')
mod01c <- rma.uni(temporal.change ~ m.wing.color + 1 + z.MAT , vi = (pop.dat$SEM)^2, data = pop.dat, method = 'ML')
mod01d <- rma.uni(temporal.change ~ m.wing.color + z.size + 1, vi = (pop.dat$SEM)^2, data = pop.dat, method = 'ML')


anova(mod01a, mod01b) # test of wing color
        # df      AIC      BIC     AICc   logLik     LRT   pval       QE  tau^2      R^2 
# Full     5 103.5418 108.2640 108.1572 -46.7709                 66.9820 0.0482          
# Reduced  4 112.3411 116.1189 115.1983 -52.1706 10.7993 0.0010 107.9676 0.1692 71.5019% 


anova(mod01a, mod01c) # test of body size
        # df      AIC      BIC     AICc   logLik    LRT   pval      QE  tau^2      R^2 
# Full     5 103.5418 108.2640 108.1572 -46.7709               66.9820 0.0482          
# Reduced  4 108.3091 112.0869 111.1663 -50.1546 6.7673 0.0093 84.8115 0.1172 58.8586% 

anova(mod01a, mod01d) # test of range-wide temperature
        # df      AIC      BIC     AICc   logLik    LRT   pval      QE  tau^2      R^2 
# Full     5 103.5418 108.2640 108.1572 -46.7709               66.9820 0.0482          
# Reduced  4 105.6897 109.4674 108.5468 -48.8448 4.1479 0.0417 74.9256 0.0731 34.0790% 




### test if including a phylogeny is better
pop.dat$species <- pop.dat$binom # need to make a separate column called "species" to align the data to the phylogeny within a metafor-based function
tree <- vcv(phylo, cor = TRUE)

mod.phylo.ml <- rma.mv(temporal.change ~ m.wing.color + z.size + z.MAT, random = ~1|species, R = list(species = tree), V = (pop.dat$SEM)^2, data = pop.dat, method = 'ML')
AICc(mod.phylo.ml) # 113.8939. model with phylo more complex and worse
AICc(mod01a) # 108.1572



  
  