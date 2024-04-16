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


## let's start by taking a pgls approach - no log transformation
nl.mod00a <- gls(risk ~ m.wing.color + z.MAT, data = r.dat, correlation = corPagel(1, phy = r.tree, fixed = TRUE, form = ~binom))
nl.mod00b <- gls(risk ~ m.wing.color + z.MAT, data = r.dat, correlation = corPagel(1, phy = r.tree, fixed = FALSE, form = ~binom))
nl.mod00c <- gls(risk ~ m.wing.color + z.MAT, data = r.dat, correlation = corPagel(0.5, phy = r.tree, fixed = FALSE, form = ~binom))
nl.mod00d <- gls(risk ~ m.wing.color + z.MAT, data = r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom))
AICc(nl.mod00a, nl.mod00b, nl.mod00c, nl.mod00d) 
          # df     AICc
# nl.mod00a  4 356.6808
# nl.mod00b  5 131.8691
# nl.mod00c  5 131.8691
# nl.mod00d  4 129.8954


summary(nl.mod00d) 
# Generalized least squares fit by REML
  # Model: risk ~ m.wing.color + z.MAT 
  # Data: r.dat 
       # AIC     BIC    logLik
  # 129.7646 144.685 -60.88232

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
# lambda 
     # 0 

# Coefficients:
                   # Value  Std.Error  t-value p-value
# (Intercept)    1.0740541 0.02032757 52.83732  0.0000
# m.wing.colory -0.0723496 0.03403816 -2.12555  0.0343
# z.MAT          0.0373439 0.01633359  2.28632  0.0229

 # Correlation: 
              # (Intr) m.wng.
# m.wing.colory -0.598       
# z.MAT         -0.020  0.034

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -0.47846104 -0.28580941 -0.16839506 -0.03122421 10.05768348 

# Residual standard error: 0.2874156 
# Degrees of freedom: 311 total; 308 residual

summary(nl.mod00b) # estimate for pagel's lambda
# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
      # lambda 
# -0.005527406 

##
nl.mod01a <- gls(risk ~ m.wing.color + z.MAT, data =r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom), method = 'ML')
nl.mod01b <- gls(risk ~ m.wing.color + 1, data = r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom), method = 'ML')
nl.mod01c <- gls(risk ~ 1 + z.MAT, data = r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom), method = 'ML')

anova(nl.mod01a, nl.mod01b) # test effect of range-wide temp
          # Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# nl.mod01a     1  4 112.0395 126.9987 -52.01975                        
# nl.mod01b     2  3 115.2734 126.4928 -54.63670 1 vs 2 5.233905  0.0222

anova(nl.mod01a, nl.mod01c) # test effect of male wing color
          # Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# nl.mod01a     1  4 112.0395 126.9987 -52.01975                        
# nl.mod01c     2  3 114.5683 125.7877 -54.28415 1 vs 2 4.528813  0.0333

### test for effect of body size
nl.bs.mod00a <- gls(risk ~ m.wing.color + z.MAT + z.size, data = r.dat, correlation = corPagel(1, phy = r.tree, fixed = TRUE, form = ~binom))
nl.bs.mod00b <- gls(risk ~ m.wing.color + z.MAT + z.size, data = r.dat, correlation = corPagel(1, phy = r.tree, fixed = FALSE, form = ~binom))
nl.bs.mod00c <- gls(risk ~ m.wing.color + z.MAT + z.size, data = r.dat, correlation = corPagel(0.5, phy = r.tree, fixed = FALSE, form = ~binom))
nl.bs.mod00d <- gls(risk ~ m.wing.color + z.MAT + z.size, data = r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom))
AICc(nl.bs.mod00a, nl.bs.mod00b, nl.bs.mod00c, nl.bs.mod00d) 
             # df     AICc
# nl.bs.mod00a  5 360.4783
# nl.bs.mod00b  6 140.1708
# nl.bs.mod00c  6 140.1708
# nl.bs.mod00d  5 138.1236


nl.bs.mod01a <- gls(risk ~ m.wing.color + z.MAT + z.size, data = r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom), method = 'ML')
nl.bs.mod01b <- gls(risk ~ m.wing.color + z.MAT + 1, data = r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom), method = 'ML')
anova(nl.bs.mod01a, nl.bs.mod01b)
             # Model df      AIC      BIC    logLik   Test      L.Ratio p-value
# nl.bs.mod01a     1  5 114.0394 132.7384 -52.01972                            
# nl.bs.mod01b     2  4 112.0395 126.9987 -52.01975 1 vs 2 4.900926e-05  0.9944

# no evidence that body size needs to go in the model to control for anything



# Number of Fisher Scoring iterations: 4
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


## start with a pgls approach - no log transformation
nl.pt.mod00a <- gls(trend.num ~ m.wing.color + z.MAT, data = pt.dat, correlation = corPagel(1, phy = pt.tree, fixed = TRUE, form = ~binom))
nl.pt.mod00b <- gls(trend.num ~ m.wing.color + z.MAT, data = pt.dat, correlation = corPagel(1, phy = pt.tree, fixed = FALSE, form = ~binom))
nl.pt.mod00c <- gls(trend.num ~ m.wing.color + z.MAT, data = pt.dat, correlation = corPagel(0.5, phy = pt.tree, fixed = FALSE, form = ~binom))
nl.pt.mod00d <- gls(trend.num ~ m.wing.color + z.MAT, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom))
AICc(nl.pt.mod00a, nl.pt.mod00b, nl.pt.mod00c, nl.pt.mod00d)
             # df        AICc
# nl.pt.mod00a  4    4.893916
# nl.pt.mod00b  5 -164.786770
# nl.pt.mod00c  5 -164.786770
# nl.pt.mod00d  4 -166.824397

summary(nl.pt.mod00d)

# Generalized least squares fit by REML
  # Model: trend.num ~ m.wing.color + z.MAT 
  # Data: pt.dat 
        # AIC       BIC   logLik
  # -166.9633 -152.2838 87.48164

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
# lambda 
     # 0 

# Coefficients:
                  # Value  Std.Error   t-value p-value
# (Intercept)   1.9836280 0.01285219 154.34161  0.0000
# m.wing.colory 0.0531245 0.02107443   2.52081  0.0122
# z.MAT         0.0016149 0.01051312   0.15361  0.8780

 # Correlation: 
              # (Intr) m.wng.
# m.wing.colory -0.610       
# z.MAT         -0.034  0.038

# Standardized residuals:
        # Min          Q1         Med          Q3         Max 
# -5.65774757 -0.20546511  0.08535792  0.09715863  5.83266452 

# Residual standard error: 0.1742331 
# Degrees of freedom: 293 total; 290 residual

summary(nl.pt.mod00c) # estimate for pagel's lambda

# Correlation Structure: corPagel
 # Formula: ~binom 
 # Parameter estimate(s):
      # lambda 
# -0.004016886 

# likelihood ratio tests for effects
nl.pt.mod01a <- gls(trend.num ~ m.wing.color + z.MAT, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom), method = "ML")
nl.pt.mod01b <- gls(trend.num ~ m.wing.color + 1, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom), method = "ML")
nl.pt.mod01c <- gls(trend.num ~ 1 + z.MAT, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom), method = "ML")

anova(nl.pt.mod01a, nl.pt.mod01b) # test effect of range-wide temp
             # Model df       AIC       BIC   logLik   Test    L.Ratio p-value
# nl.pt.mod01a     1  4 -187.4711 -172.7504 97.73555                          
# nl.pt.mod01b     2  3 -189.4473 -178.4067 97.72364 1 vs 2 0.02383773  0.8773

anova(nl.pt.mod01a, nl.pt.mod01c) # test effect of male wing color
             # Model df       AIC       BIC   logLik   Test  L.Ratio p-value
# nl.pt.mod01a     1  4 -187.4711 -172.7504 97.73555                        
# nl.pt.mod01c     2  3 -183.1202 -172.0797 94.56012 1 vs 2 6.350868  0.0117

##
nl.bs.pt.mod00a <- gls(trend.num ~ m.wing.color + z.MAT + z.size, data = pt.dat, correlation = corPagel(1, phy = pt.tree, fixed = TRUE, form = ~binom))
nl.bs.pt.mod00b <- gls(trend.num ~ m.wing.color + z.MAT + z.size, data = pt.dat, correlation = corPagel(1, phy = pt.tree, fixed = FALSE, form = ~binom))
nl.bs.pt.mod00c <- gls(trend.num ~ m.wing.color + z.MAT + z.size, data = pt.dat, correlation = corPagel(0.5, phy = pt.tree, fixed = FALSE, form = ~binom))
nl.bs.pt.mod00d <- gls(trend.num ~ m.wing.color + z.MAT + z.size, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom))
AICc(nl.bs.pt.mod00a, nl.bs.pt.mod00b, nl.bs.pt.mod00c, nl.bs.pt.mod00d)

                # df      AICc
# nl.bs.pt.mod00a  5   11.7040
# nl.bs.pt.mod00b  6 -156.4154
# nl.bs.pt.mod00c  6 -156.4154
# nl.bs.pt.mod00d  5 -158.4077


# is body size mucking things up?
nl.bs.pt.mod01a <- gls(trend.num ~ m.wing.color + z.MAT + z.size, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom), method = 'ML')
nl.bs.pt.mod01b <- gls(trend.num ~ m.wing.color + z.MAT + 1, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom), method = 'ML')
anova(nl.bs.pt.mod01a, nl.bs.pt.mod01b)
                # Model df       AIC       BIC   logLik   Test   L.Ratio p-value
# nl.bs.pt.mod01a     1  5 -186.2682 -167.8674 98.13412                         
# nl.bs.pt.mod01b     2  4 -187.4711 -172.7504 97.73555 1 vs 2 0.7971255   0.372


## use better fitting model

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


#### Figures ####

# Fig 2a - range wide
nl.mod01a <- gls(risk ~ m.wing.color + z.MAT, data =r.dat, correlation = corPagel(0, phy = r.tree, fixed = TRUE, form = ~binom), method = 'ML')
library(emmeans)
gls.out <- emmeans(nl.mod01a, specs = 'm.wing.color')

gls.dat <- data.frame(gls.out)


bgrd01 =
  theme(axis.text = element_text(color="Black"),
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),
        axis.title.y = element_text(size = 12), 
        axis.text.y = element_text(size = 8),
        panel.background = element_rect(fill = "White"),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        axis.line.x = element_line(linetype = 'solid', color = 'black', size = 0.8),
        axis.line.y = element_line(linetype = 'solid', color = 'black', size = 0.8),
        legend.position = c(1,1), legend.justification = c(1,1), 
        legend.title = element_blank(), legend.text = element_text(size = 9),
        legend.key = element_rect(fill = 'white'),
        plot.margin = margin(10,10,10,10))


risk_gls_plot <-
  ggplot(data = gls.dat, aes(x = m.wing.color, y = exp(emmean), group = m.wing.color)) +
  geom_linerange(aes(x = m.wing.color, ymin = exp(emmean - SE), ymax = exp(emmean + SE), color = m.wing.color), size = 1) +
  geom_point(shape = 22, size = 5, stroke = 1, aes(fill = m.wing.color, color = m.wing.color)) +
  scale_color_manual(values = c('#8497B0', '#714D21'), guide = "none") +
  scale_fill_manual(values = c('#DAE3F3', '#EBC06B'), guide = "none") +
  scale_x_discrete(labels = c('Non-Ornamented\nSpecies\n(n = 200)', 'Ornamented\nSpecies\n(n = 111)')) +
  labs(y = 'IUCN Extinction Risk') +
  bgrd01 + 
  theme(axis.text.x = element_text(size = 10))
risk_gls_plot

#svg('risk_gls_plot.svg', height = 3, width = 3)
#print(risk_gls_plot)
#dev.off()




# Fig 2b - pop trends
nl.pt.mod01a <- gls(trend.num ~ m.wing.color + z.MAT, data = pt.dat, correlation = corPagel(0, phy = pt.tree, fixed = TRUE, form = ~binom), method = "ML")

pop.tab <- emmeans(nl.pt.mod01a, specs = 'm.wing.color')


pop.dat <- data.frame(pop.tab)
pop.dat$est <- pop.dat$emmean
pop.dat$est.l.sem <- pop.dat$emmean - pop.dat$SE
pop.dat$est.u.sem <- pop.dat$emmean + pop.dat$SE
pop.dat$lcl <- pop.dat$lower.CL
pop.dat$ucl <- pop.dat$upper.CL

gls_trend_plot <-
  ggplot(data = pop.dat, aes(x = m.wing.color, y = est, group = m.wing.color)) +
  geom_linerange(aes(x = m.wing.color, ymin = est.l.sem, ymax = est.u.sem, color = m.wing.color), size = 1) +
  geom_point(shape = 22, size = 5, stroke = 1, aes(fill = m.wing.color, color = m.wing.color)) +
  scale_color_manual(values = c('#8497B0', '#714D21'), guide = "none") +
  scale_fill_manual(values = c('#DAE3F3', '#EBC06B'), guide = "none") +
  geom_hline(yintercept = 2, linetype = 32, alpha = 1, color = '#999999') +
  scale_y_continuous(limits = c(1.96,2.056)) +
  scale_x_discrete(labels = c('Non-Ornamented\nSpecies\n(n = 184)', 'Ornamented\nSpecies\n(n = 109)')) +
  ylab('IUCN Demographic Trend') +
  bgrd01+ 
  theme(axis.text.x = element_text(size = 10))
gls_trend_plot

#svg('gls_trend_plot.svg', height = 3, width = 3)
#print(gls_trend_plot)
#dev.off()

#trend.gls.fig <- gls.trend.plot + bgrd01

#png('supp.ptrend.gls.png', height = 5, width = 2.75, units = 'in', res = 650)
#print(trend.gls.fig)
#dev.off()



library(cowplot)
iucn_plots <- plot_grid(risk_gls_plot, gls_trend_plot, align = "v", label_size = 12, nrow = 1, ncol = 2)
iucn_plots

svg('iucn_plots.svg', height = 3, width = 6.5)
print(iucn_plots)
dev.off()
