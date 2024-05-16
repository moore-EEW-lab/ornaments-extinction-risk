library(car)
library(ggplot2)
library(MASS)
library(ggplot2)
library(viridis)


## file with predicted valuees for impervious surface, heating, agricultural land cover, and historic temp. The spatially informed models can take multiple hours to run depending on computational power, so separate files with model outputs were created to facilitate figure design
imp.fit.dat <- read.csv('ss.ext.assemb.imperv.model.csv')
warm.fit.dat <- read.csv('ss.ext.assemb.warming.model.csv')
ag.fit.dat <- read.csv('ss.ext.assemb.ag.model.csv')
hist.fit.dat <- read.csv('ss.ext.assemb.hist.model.csv')

## raw data for density plots
mw <- read.csv('assemblage.comp.csv')
mw.eff <- subset(mw, effort >4)



### Figure 3c - impervious surface

# middle panel
bgrd01 =
theme(axis.text = element_text(color="Black"),
axis.title.x = element_text(face = "plain", size = 14, margin = margin(t = 10)), 
axis.text.x = element_text(size = 10),
axis.title.y = element_text(face = "plain", size = 14, margin = margin(r = 10)), 
axis.text.y = element_text(size = 10),
panel.background = element_rect(fill = "White"),
panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
axis.line.x = element_line(linetype = 'solid', color = 'black', size = 0.8),
axis.line.y = element_line(linetype = 'solid', color = 'black', size = 0.8),
legend.key = element_rect(fill = 'white')) 

summary(imp.fit.dat$imperv)
imperv.plot <-
ggplot() +
geom_ribbon(data = imp.fit.dat, aes(x = imperv, ymax = upper.CI, ymin = lower.CI), alpha = 0.3, fill = 'black', linetype = 0) +
geom_line(data = imp.fit.dat, aes(x = imperv, y = fitted.prop), color = 'black', linetype = 1, size = 1.7) +
scale_y_continuous(breaks = c(0.55, 0.61, 0.67, 0.73)) +
scale_x_continuous(breaks = c(0,25,50,75), limits = c(0, 75)) +
  coord_cartesian(ylim = c(0.53, 0.738)) +
labs(x = 'Impervious Surface Area (%)', y = '% Ornamented Spp in Assemblage') +
bgrd01
imperv.plot

## top panel
imperv.dens.plot <-
ggplot(data = mw.eff, aes(x = impervious_percent, y = ..scaled..)) +
geom_density(fill = 'black', alpha = 0.3) +
scale_y_continuous(breaks = c(0, 0.33, 0.66, 0.99)) +
  scale_x_continuous(breaks = c(0,25,50,75), limits = c(0, 75)) +
  ylab(label = "Rel. Freq.") +
   bgrd01 + 
  theme(legend.position = "none", axis.title.x = element_blank())

imperv_plots <- plot_grid(imperv.dens.plot, imperv.plot, ncol = 1, nrow = 2, rel_heights = c(0.3,0.7), align = "h")
imperv_plots
pdf('imperv_plots.pdf', height = 5, width = 4)
print(imperv_plots)
dev.off()
svg('imperv_plots.svg', height = 5, width = 4)
print(imperv_plots)
dev.off()




## figure 3a - historical temp



## middle panel
clim.temp.plot <-
ggplot(data = hist.fit.dat) +
geom_ribbon(aes(x = temp, ymax = upper.CI, ymin = lower.CI), alpha = 0.3, fill = 'purple', linetype = 0) +
geom_line(aes(x = temp, y = fitted.prop), color = 'purple', linetype = 1, size = 1.7) +
scale_y_continuous(breaks = c(0.55, 0.61, 0.67, 0.73)) +
scale_x_continuous(breaks = c(12, 14, 16,18,20,22, 24)) +
coord_cartesian(ylim = c(0.53, 0.738)) +
labs(x = 'Historical Max Temperature (°C)', y = '% Ornamented Spp in Assemblage') +
bgrd01
clim.temp.plot

clim.dens.plot <- 
ggplot(data = mw.eff, aes(x = avmaxt_61to65, y = ..scaled..)) +
geom_density(color = 'purple', fill = 'purple', alpha = 0.3) +
scale_y_continuous(breaks = c(0, 0.33, 0.66, 0.99)) +
scale_x_continuous(breaks = c(12, 14, 16,18,20,22, 24)) +
  ylab(label = "Rel. Freq.") +
  bgrd01 +
  theme(legend.position = "none", axis.title.x = element_blank())
clim.dens.plot

clim_max_plots <- plot_grid(clim.dens.plot , clim.temp.plot, ncol = 1, nrow = 2, rel_heights = c(0.3,0.7), align = "h")
clim_max_plots
pdf('clim_max_plots.pdf', height = 5, width = 4)
print(clim_max_plots)
dev.off()
svg('clim_max_plots.svg', height = 5, width = 4)
print(clim_max_plots)
dev.off()




## Figure 3b - warming


#middle panel
warm.plot <-
ggplot(data = warm.fit.dat) +
geom_ribbon(aes(x = warm, ymax = upper.CI, ymin = lower.CI), alpha = 0.3, fill = 'red4', linetype = 0) +
geom_line(aes(x = warm, y = fitted.prop), color = 'red4', linetype = 1, size = 1.7) +
scale_y_continuous(breaks = c(0.55, 0.61, 0.67, 0.73)) +
scale_x_continuous(breaks = c(0, 0.55, 1.10)) +
coord_cartesian(ylim = c(0.53, 0.738)) +
labs(x = 'Change in Max Temperature (°C)', y = '% Ornamented Spp in Assemblage') +
bgrd01
warm.plot 

# top panel
warm.dens.plot <- 
ggplot(data = mw.eff, aes(x = change_maxt, y = ..scaled..)) +
geom_density(color = 'red4', fill = 'red4', alpha = 0.3) +
scale_y_continuous(breaks = c(0, 0.33, 0.66, 0.99)) +
scale_x_continuous(breaks = c(0, 0.55, 1.10)) +
  ylab(label = "Rel. Freq.") +
  bgrd01 +
  theme(legend.position = "none", axis.title.x = element_blank())


warm_plots <- plot_grid(warm.dens.plot , warm.plot, ncol = 1, nrow = 2, rel_heights = c(0.3,0.7), align = "h")
warm_plots
pdf('warm_plots.pdf', height = 5, width = 4)
print(warm_plots)
dev.off()
svg('warm_plots.svg', height = 5, width = 4)
print(warm_plots)
dev.off()


##### Figure 3d - agricultural land area

#middle panel
ag.plot <-
  ggplot(data = ag.fit.dat) +
  geom_ribbon(aes(x = ag, ymax = upper.CI, ymin = lower.CI), alpha = 0.3, fill = 'forestgreen', linetype = 0) +
  geom_line(aes(x = ag, y = fitted.prop), color = 'forestgreen', linetype = 1, size = 1.7) +
  scale_y_continuous(breaks = c(0.55, 0.61, 0.67, 0.73))+
  scale_x_continuous(breaks = c(0,30,60,90), limits = c(0, 96)) +
  coord_cartesian(ylim = c(0.53, 0.738)) +
  labs(x = 'Agricultural Land Area (%)', y = '% Ornamented Spp in Assemblage') +
  bgrd01
ag.plot

# top panel
ag.dens.plot <- 
  ggplot(data = mw.eff, aes(x = agricul*100, y = ..scaled..)) +
  geom_density(color = 'forestgreen', fill = 'forestgreen', alpha = 0.3) +
  scale_y_continuous(breaks = c(0, 0.33, 0.66, 0.99)) +
  scale_x_continuous(breaks = c(0,30,60,90), limits = c(0, 96)) +
  ylab(label = "Rel. Freq.") +
  bgrd01 + 
  theme(legend.position = "none", axis.title.x = element_blank())

ag.dens.plot

ag_plots <- plot_grid(ag.dens.plot, ag.plot, ncol = 1, nrow = 2, rel_heights = c(0.3,0.7), align = "h")
ag_plots
pdf('ag_plots.pdf', height = 5, width = 4)
print(ag_plots)
dev.off()
svg('ag_plots.svg', height = 5, width = 4)
print(ag_plots)
dev.off()



