setwd('/Users/michaelmoore/desktop/Working Directory')

library(phytools)
library(car)
library(MASS)
library(nlme)
library(ggplot2)
library(lme4)
library(lmerTest)
library(lsmeans)
library(ggplot2)
library(plyr)
library(letsR)
library(phylolm)
library(phylosignal)
library(geiger)
library(MuMIn)
library(emmeans)
library(phylobase)

r.dat <- read.csv('species.ext.csv') 
rownames(r.dat) <- r.dat$binom
r.tree <- read.newick('species.ext.tre')


## let's start by taking a pgls approach
mod00a <- gls(log(risk) ~ m.wing.color + z.MAT, data = r.dat, correlation = corPagel(1, phy = r.tree, fixed = TRUE, form = ~binom))
mod00b <- gls(log(risk) ~ m.wing.color + z.MAT, data = r.dat, correlation = corPagel(1, phy = r.tree, fixed = FALSE, form = ~binom))
mod00c <- gls(log(risk) ~ m.wing.color + z.MAT, data = r.dat, correlation = corPagel(0.5, phy = r.tree, fixed = FALSE, form = ~binom))
mod00d <- gls(log(risk) ~ m.wing.color + z.MAT, data = r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom))
AICc(mod00a, mod00b, mod00c, mod00d) 

       # df       AICc
# mod00a  4   24.23801
# mod00b  5 -228.60426
# mod00c  5 -271.16569
# mod00d  4 -230.66237

## different AICc for different starting values for Pagel's lambda indicates model isn't converging well. Star phylogeny MUCH better supported than Brownian Motion

summary(mod00d)

# Generalized least squares fit by REML
  # Model: log(risk) ~ m.wing.color + z.MAT 
  # Data: big.dat 
        # AIC       BIC   logLik
  # -230.7931 -215.8727 119.3965

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
# lambda 
     # 0 

# Coefficients:
                    # Value   Std.Error   t-value p-value
# (Intercept)    0.04366090 0.011320971  3.856639  0.0001
# m.wing.colory -0.04268040 0.018956774 -2.251459  0.0251
# z.MAT          0.02148287 0.009096616  2.361633  0.0188

 # Correlation: 
              # (Intr) m.wng.
# m.wing.colory -0.598       
# z.MAT         -0.020  0.034

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -0.50084120 -0.30184389 -0.17507497 -0.03225269  8.26127089 

# Residual standard error: 0.1600695 
# Degrees of freedom: 311 total; 308 residual


## hypothesis testing
mod01a <- gls(log(risk) ~ m.wing.color + z.MAT, data =r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom), method = 'ML')
mod01b <- gls(log(risk) ~ m.wing.color + 1, data = r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom), method = 'ML')
mod01c <- gls(log(risk) ~ 1 + z.MAT, data = r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom), method = 'ML')

anova(mod01a, mod01b) # test effect of range-wide temp
       # Model df       AIC       BIC   logLik   Test  L.Ratio p-value
# mod01a     1  4 -252.0302 -237.0710 130.0151                        
# mod01b     2  3 -248.4489 -237.2295 127.2245 1 vs 2 5.581251  0.0182

anova(mod01a, mod01c) # test effect of male wing color
       # Model df       AIC      BIC   logLik   Test  L.Ratio p-value
# mod01a     1  4 -252.0302 -237.071 130.0151                        
# mod01c     2  3 -248.9534 -237.734 127.4767 1 vs 2 5.076779  0.0242


### body size is thought to maybe be involved. let's see if we should have included it. 
bs.mod00a <- gls(log(risk) ~ m.wing.color + z.MAT + z.size, data = r.dat, correlation = corPagel(1, phy = r.tree, fixed = TRUE, form = ~binom))
bs.mod00b <- gls(log(risk) ~ m.wing.color + z.MAT + z.size, data = r.dat, correlation = corPagel(1, phy = r.tree, fixed = FALSE, form = ~binom))
bs.mod00c <- gls(log(risk) ~ m.wing.color + z.MAT + z.size, data = r.dat, correlation = corPagel(0.5, phy = r.tree, fixed = FALSE, form = ~binom))
bs.mod00d <- gls(log(risk) ~ m.wing.color + z.MAT + z.size, data = r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom))
AICc(bs.mod00a, bs.mod00b, bs.mod00c, bs.mod00d) 
          # df       AICc
# bs.mod00a  5   29.06336
# bs.mod00b  6 -219.34738
# bs.mod00c  6 -219.34738
# bs.mod00d  5 -221.38211

# pagel's lambda converges this time, but star phylogeny still best supported

# let's test whether body size has a significant effect
bs.mod01a <- gls(log(risk) ~ m.wing.color + z.MAT + z.size, data = r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom), method = 'ML')
bs.mod01b <- gls(log(risk) ~ m.wing.color + z.MAT + 1, data = r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom), method = 'ML')

anova(bs.mod01a, bs.mod01b) # test effect of body size

          # Model df       AIC       BIC   logLik   Test   L.Ratio p-value
# bs.mod01a     1  5 -250.1504 -231.4514 130.0752                         
# bs.mod01b     2  4 -252.0302 -237.0710 130.0151 1 vs 2 0.1202266  0.7288

# not even close.


### since there's no effect of phylogeny, we can double check the qualitative results of these analyses using a generalized linear model with a quasipoisson error distribution.

mod02 <- glm(risk ~ m.wing.color + z.MAT, data = r.dat, family = 'quasipoisson')
summary(mod02)

# Call:
# glm(formula = risk ~ m.wing.color + z.MAT, family = "quasipoisson", 
    # data = r.dat)

# Deviance Residuals: 
     # Min        1Q    Median        3Q       Max  
# -0.13488  -0.07965  -0.04644  -0.00793   2.11505  

# Coefficients:
              # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    0.07082    0.01866   3.796 0.000177 ***
# m.wing.colory -0.06990    0.03195  -2.188 0.029447 *  
# z.MAT          0.03595    0.01528   2.353 0.019272 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# (Dispersion parameter for quasipoisson family taken to be 0.07466565)

    # Null deviance: 16.789  on 310  degrees of freedom
# Residual deviance: 15.988  on 308  degrees of freedom
# AIC: NA

# Number of Fisher Scoring iterations: 4

mod02b <- glm(risk ~ m.wing.color + 1, data = r.dat, family = 'quasipoisson')
mod02c <- glm(risk ~ 1 + z.MAT, data = r.dat, family = 'quasipoisson') 


anova(mod02, mod02b, test = 'F') # test of effect of range-wide temp
# Analysis of Deviance Table

# Model 1: risk ~ m.wing.color + z.MAT
# Model 2: risk ~ m.wing.color + 1
  # Resid. Df Resid. Dev Df Deviance      F  Pr(>F)  
# 1       308     15.988                             
# 2       309     16.403 -1 -0.41491 5.5569 0.01903 *

anova(mod02, mod02c, test = 'F') # test of effect of male wing color
# Analysis of Deviance Table

# Model 1: risk ~ m.wing.color + z.MAT
# Model 2: risk ~ 1 + z.MAT
  # Resid. Df Resid. Dev Df Deviance      F  Pr(>F)  
# 1       308     15.988                             
# 2       309     16.348 -1 -0.35993 4.8206 0.02887 *


### should we have included body size in this model?
bs.mod02 <- glm(risk ~ m.wing.color + z.MAT + z.size, data = r.dat, family = 'quasipoisson')
anova(bs.mod02, mod02, test = 'F') # test of effect of body size

# Analysis of Deviance Table

# Model 1: risk ~ m.wing.color + z.MAT + z.size
# Model 2: risk ~ m.wing.color + z.MAT
  # Resid. Df Resid. Dev Df    Deviance  F Pr(>F)
# 1       307     15.988                         
# 2       308     15.988 -1 -1.0914e-06  0  0.997

## not even close.

###### let's look at population trends next

pt.dat <- read.csv('species.ext.pt.csv')
rownames(pt.dat) <- pt.dat$binom
pt.tree <- read.newick('species.ext.pt.tre')

## start with a pgls approach
pt.mod00a <- gls(log(trend.num) ~ m.wing.color + z.MAT, data = pt.dat, correlation = corPagel(1, phy = pt.tree, fixed = TRUE, form = ~binom))
pt.mod00b <- gls(log(trend.num) ~ m.wing.color + z.MAT, data = pt.dat, correlation = corPagel(1, phy = pt.tree, fixed = FALSE, form = ~binom))
pt.mod00c <- gls(log(trend.num) ~ m.wing.color + z.MAT, data = pt.dat, correlation = corPagel(0.5, phy = pt.tree, fixed = FALSE, form = ~binom))
pt.mod00d <- gls(log(trend.num) ~ m.wing.color + z.MAT, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom))
AICc(pt.mod00a, pt.mod00b, pt.mod00c, pt.mod00d)

          # df      AICc
# pt.mod00a  4 -304.6028
# pt.mod00b  5 -509.2304
# pt.mod00c  5 -509.2304
# pt.mod00d  4 -510.9810

## star phylogeny is best

summary(pt.mod00d)
# Generalized least squares fit by REML
  # Model: log(trend.num) ~ m.wing.color + z.MAT 
  # Data: pt.dat 
        # AIC       BIC logLik
  # -511.1199 -496.4404 259.56

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
# lambda 
     # 0 

# Coefficients:
                   # Value   Std.Error  t-value p-value
# (Intercept)    0.6804030 0.007100327 95.82699  0.0000
# m.wing.colory  0.0275252 0.011642786  2.36415  0.0187
# z.MAT         -0.0028771 0.005808082 -0.49535  0.6207

 # Correlation: 
              # (Intr) m.wng.
# m.wing.colory -0.610       
# z.MAT         -0.034  0.038

# Standardized residuals:
       # Min         Q1        Med         Q3        Max 
# -7.0612452 -0.1381344  0.1071588  0.1411674  4.3471030 

# Residual standard error: 0.0962569 
# Degrees of freedom: 293 total; 290 residual

pt.mod01a <- gls(log(trend.num) ~ m.wing.color + z.MAT, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom), method = "ML")
pt.mod01b <- gls(log(trend.num) ~ m.wing.color + 1, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom), method = "ML")
pt.mod01c <- gls(log(trend.num) ~ 1 + z.MAT, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom), method = "ML")

anova(pt.mod01a, pt.mod01b) # test effect of range-wide temp
          # Model df       AIC       BIC   logLik   Test   L.Ratio p-value
# pt.mod01a     1  4 -535.1880 -520.4673 271.5940                         
# pt.mod01b     2  3 -536.9402 -525.8996 271.4701 1 vs 2 0.2478097  0.6186

anova(pt.mod01a, pt.mod01c) # test effect of male wing color
          # Model df       AIC       BIC   logLik   Test  L.Ratio p-value
# pt.mod01a     1  4 -535.1880 -520.4673 271.5940                        
# pt.mod01c     2  3 -531.5947 -520.5542 268.7973 1 vs 2 5.593274   0.018




# should we be including body size in these models?
bs.pt.mod00a <- gls(log(trend.num) ~ m.wing.color + z.MAT + z.size, data = pt.dat, correlation = corPagel(1, phy = pt.tree, fixed = TRUE, form = ~binom))
bs.pt.mod00b <- gls(log(trend.num) ~ m.wing.color + z.MAT + z.size, data = pt.dat, correlation = corPagel(1, phy = pt.tree, fixed = FALSE, form = ~binom))
bs.pt.mod00c <- gls(log(trend.num) ~ m.wing.color + z.MAT + z.size, data = pt.dat, correlation = corPagel(0.5, phy = pt.tree, fixed = FALSE, form = ~binom))
bs.pt.mod00d <- gls(log(trend.num) ~ m.wing.color + z.MAT + z.size, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom))
AICc(bs.pt.mod00a, bs.pt.mod00b, bs.pt.mod00c, bs.pt.mod00d)
             # df      AICc
# bs.pt.mod00a  5 -296.6282
# bs.pt.mod00b  6 -499.6318
# bs.pt.mod00c  6 -549.3972
# bs.pt.mod00d  5 -501.2994

# pagel's lambda model not converging. star phylo is best

bs.pt.mod01a <- gls(log(trend.num) ~ m.wing.color + z.MAT + z.size, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom), method = 'ML')
bs.pt.mod01b <- gls(log(trend.num) ~ m.wing.color + z.MAT + 1, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom), method = 'ML')
anova(bs.pt.mod01a, bs.pt.mod01b)
             # Model df       AIC       BIC   logLik   Test   L.Ratio p-value
# bs.pt.mod01a     1  5 -533.9059 -515.5050 271.9529                         
# bs.pt.mod01b     2  4 -535.1880 -520.4673 271.5940 1 vs 2 0.7178947  0.3968

## no significant effect of body size, so don't need to worry about including it

### generalized linear model approach
pt.mod02 <- glm(trend.num ~ m.wing.color + z.MAT, data = pt.dat, family = 'quasipoisson')
summary(pt.mod02)

# Call:
# glm(formula = trend.num ~ m.wing.color + z.MAT, family = "quasipoisson", 
    # data = big.dat2)

# Deviance Residuals: 
     # Min        1Q    Median        3Q       Max  
# -0.77427  -0.02515   0.01055   0.01200   0.67025  

# Coefficients:
               # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.6849276  0.0064400 106.355   <2e-16 ***
# m.wing.colory 0.0264290  0.0104732   2.523   0.0122 *  
# z.MAT         0.0008045  0.0052369   0.154   0.8780    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# (Dispersion parameter for quasipoisson family taken to be 0.01511968)

    # Null deviance: 4.6171  on 292  degrees of freedom
# Residual deviance: 4.5209  on 290  degrees of freedom
# AIC: NA

# Number of Fisher Scoring iterations: 4

pt.mod02b <- glm(trend.num ~ m.wing.color + 1, data = pt.dat, family = 'quasipoisson')
pt.mod02c <- glm(trend.num ~ 1 + z.MAT, data = pt.dat, family = 'quasipoisson')

anova(pt.mod02, pt.mod02b, test = 'F') # test effect of range-wide temp

# Analysis of Deviance Table

# Model 1: trend.num ~ m.wing.color + z.MAT
# Model 2: trend.num ~ m.wing.color + 1
  # Resid. Df Resid. Dev Df    Deviance      F Pr(>F)
# 1       290     4.5209                             
# 2       291     4.5213 -1 -0.00035682 0.0236  0.878

anova(pt.mod02, pt.mod02c, test = 'F') # test effect of wing color

# Analysis of Deviance Table

# Model 1: trend.num ~ m.wing.color + z.MAT
# Model 2: trend.num ~ 1 + z.MAT
  # Resid. Df Resid. Dev Df  Deviance      F  Pr(>F)  
# 1       290     4.5209                              
# 2       291     4.6170 -1 -0.096073 6.3542 0.01225 *

## should we include body size in this analysis?
bs.pt.mod02 <- glm(trend.num ~ m.wing.color + z.MAT +z.size, data = pt.dat, family = 'quasipoisson')
anova(bs.pt.mod02, pt.mod02, test = 'F') # test if adding body size improves the model

# Analysis of Deviance Table

# Model 1: trend.num ~ m.wing.color + z.MAT + z.size
# Model 2: trend.num ~ m.wing.color + z.MAT
  # Resid. Df Resid. Dev Df  Deviance      F Pr(>F)
# 1       289     4.5090                           
# 2       290     4.5209 -1 -0.011993 0.7928  0.374

# no need to include body size as a covariate


