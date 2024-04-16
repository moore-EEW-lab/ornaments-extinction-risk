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


#### make figures
 
## fig 1a 
occ.odes <- read.csv("pop.ext.csv")

head(occ.odes)

sort.occ.odes <- occ.odes[order(-occ.odes$upper),]
sort.occ.odes$ID <- seq(from = 1, to = 19, by = 1)
sort.occ.odes$m.wing.color <- factor(sort.occ.odes$m.wing.color, levels = c('n', 'y'))


bgrd01 =
theme(axis.text = element_text(color="Black"),
axis.title.x = element_blank(), 
axis.text.x = element_text(size = 8, angle = 55, hjust = 1),
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




sort.occ.odes$sig <- ifelse(sort.occ.odes$lower > 0, "pos", ifelse(sort.occ.odes$upper < 0, "neg", "non"))

pop_occ_plot <-
  ggplot() +
  geom_hline(yintercept = 0, linetype = 32, alpha = 1, color = '#999999') +
  geom_linerange(data = sort.occ.odes, aes(x = binom, ymin=lower, ymax = upper, color = m.wing.color), size = 2) +
  scale_color_manual(values = c('#DAE3F3', '#EBC06B'), labels = c('Non-Ornamented', 'Ornamented')) +
  new_scale_color() +
  geom_point(data = sort.occ.odes, aes(x = binom, y = temporal.change, color = m.wing.color), shape = 15, size = 1.5) +
  scale_color_manual(values = c('#8497B0', '#714D21'), labels = c('Non-Ornamented', 'Ornamented')) +
  geom_point(data = sort.occ.odes, aes(x = binom, y = -32, shape = as.factor(sig)), size = 4) +
  scale_shape_manual(values = c("-","","+"), guide = "none") +
  ylab('Annual change in probability of\npopulation occurance (logit)') + 
  scale_x_discrete(limits = c('Tramea_lacerata',          'Libellula_luctuosa',       'Plathemis_lydia',          'Pantala_hymenaea',       'Erythemis_collocata',      'Libellula_saturata',       'Aeshna_interrupta',        'Anax_junius',             'Libellula_forensis',       'Pachydiplax_longipennis',  'Gomphus_kurilis',          'Libellula_pulchella',      'Libellula_quadrimaculata', 'Rhionaeschna_californica', 'Cordulegaster_dorsalis', 'Sympetrum_corruptum',     'Sympetrum_illotum',        'Sympetrum_pallipes',       'Libellula_nodisticta'), labels = c('Tramea lacerata',          'Libellula luctuosa',       'Plathemis lydia',          'Pantala hymenaea',       'Erythemis collocata',      'Libellula saturata',       'Aeshna interrupta',        'Anax junius',             'Libellula forensis',       'Pachydiplax longipennis',  'Gomphus kurilis',          'Libellula pulchella',      'Libellula quadrimaculata', 'Rhionaeschna californica', 'Cordulegaster dorsalis', 'Sympetrum corruptum',     'Sympetrum illotum',        'Sympetrum pallipes',       'Libellula nodisticta')) + 
  annotate("text", x = 13, y = 50, hjust = 0, label = "Ornamented", size = 10/.pt) +
  annotate("text", x = 13, y = 20, hjust = 0, label = "Non-Ornamented", size = 10/.pt) +
  annotate("segment", x = 12.5, xend = 12.5, y = 17, yend = 23, size = 2, color = "#DAE3F3") +
  annotate("segment", x = 12.5, xend = 12.5, y = 47, yend = 53, size = 2, color = "#EBC06B") +
  annotate("point", x = 12.5, y = 20, shape = 15, size = 1.5, color = "#8497B0") +
  annotate("point", x = 12.5, y = 50, shape = 15, size = 1.5, color = "#714D21") +
  #annotate("text", x = 10, y = 45, hjust = 0, label = "Ornamented") +
  bgrd01 +
  theme(legend.position="none", axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1))

pop_occ_plot

svg('pop_occ_blank_plot.svg', height = 4.28, width = 4.75)
print(pop_occ_plot)
dev.off()




#######

## fig 1b

#parameter estimates from meta-analytic model
m.wing.color <- c('n', 'y')
est <- c(-0.0640, 0.6194)
sem <- c(0.1384, 0.2051)
lcl <- c(-0.3352,0.2174)
ucl <- c(0.2073, 1.0214)
temp.dat <- data.frame(m.wing.color, est, sem, lcl, ucl)



temp_change_plot <-
  ggplot(temp.dat, aes(x = m.wing.color, y = est, group = m.wing.color)) +
  geom_hline(yintercept = 0, linetype = 32, alpha = 1, color = '#999999') +
  geom_linerange(aes(x = m.wing.color, ymin = lcl, ymax = ucl, color = m.wing.color), size = 1) +
  geom_point(shape = 22, size = 5, stroke = 1, aes(fill = m.wing.color, color = m.wing.color)) +
  scale_color_manual(values = c('#8497B0', '#714D21'), guide = "none") +
  scale_fill_manual(values = c('#DAE3F3', '#EBC06B'), guide = "none") +
  scale_x_discrete(labels = c('Non-Ornamented\nSpecies', 'Ornamented\nSpecies')) +
  scale_y_continuous(breaks = c(-0.3, 0, 0.3, 0.6, 0.9)) +
  ylab('Annual change in probability of\npopulation occurance (logit)') + 
  bgrd01 +
  theme(axis.text.x = element_text(size = 10,angle = 50, vjust = 1, hjust = 1))

temp_change_plot



temp_change_plots <- plot_grid(pop_occ_plot, temp_change_plot, labels = c("(a)","(b)"), label_size = 12, nrow = 1, ncol = 2, rel_widths = c(0.7,0.3))
temp_change_plots


svg('temp_change_plots.svg', height = 4, width = 6.5)
print(temp_change_plots)
dev.off()

png('pop-level.annu.change.png', height = 4, width = 2.4, units = 'in', res = 650)
print(temp.change.fig)
dev.off()  
  