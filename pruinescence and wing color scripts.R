
dat <- read.csv('m.orn.pruin.csv')


libs <- subset(dat, Family == 'Libellulidae')
# chisq.test(table(libs$m.wing.color, libs$pruinosity))
	# Pearson's Chi-squared test with Yates' continuity correction

# data:  table(libs$m.wing.color, libs$pruinosity)
# X-squared = 1.4473, df = 1, p-value = 0.229


gomphs <- subset(dat, Family == 'Gomphidae')
chisq.test(table(gomphs$m.wing.color, gomphs$pruinosity))
	# Pearson's Chi-squared test with Yates' continuity correction

# data:  table(gomphs$m.wing.color, gomphs$pruinosity)
# X-squared = 1.528e-31, df = 1, p-value = 1

