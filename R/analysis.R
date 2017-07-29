# Calibrating a multi-species coalescent tree of mouse lemurs to
# geological times of divergence. This tutorial reproduces the
# results from:

# Yoder et al. (2016) "Geogenetic patterns in mouse lemurs (genus 
# Microcebus) reveal the ghosts of Madagascar’s forests past" 
# Proceedings of the National Academy of Sciences, 113: 8049–8056.

# Detail tutorial explanation can be found in bit.ly/mousies

rm(list=ls()) # clean up the workspace

# *******************************************************************
# Constructing priors on 'g' (genartion time) and
# 'u' (per-generation mutation rate)
# *******************************************************************
# Generation time is 3.0 to 4.5 years (see paper above).
# We build a gamma prior on the generation time, 'g' with 95% CI
# roughly within 3 to 4.5.

# Recall the following relationship between mean, 'm', and variance, 
# 'v', and parameters 'a' and 'b' in the gamma distribution:
# m = a/b -> a = bm; v = a/b^2 -> v = bm/b^2 -> v = m/b
# Thus b = m/v and a = m^2/v.

# We want the mean of the prior of 'g' to be 3.75, midway between 3 to
# 4.5. The range of the prior is 4.5 - 3 = 1.5 years, so we want the
# standard deviation to be 1.5 / 4 = 0.375, thus 
# the variance is (1.5 / 4)^2.
m <- 3.75; v <- (4.5-3)^2/16
a <- m^2/v; b <- m/v
# 'a' is 100, and 'b' is 26.67

# This is how the prior on 'g' looks like:
curve(dgamma(x, a, b), from=0, to=6, xlab="Generation time, g (y)", ylab="Gamma(g | a, b)", las=1, n=5e2)

# We reapeat the procedure above to obtain a prior on the
# per-generation mutation rate 'u', set to be between 0.54 to 1.2
# substitutions per site per 10^8 years.
m.r <- (1.2+.54)/2; v.r <- (1.2-.54)^2/16
a.r <- m.r^2/v.r; b.r <- m.r/v.r
# 'a.r' is 27.80 and 'b.r' is 31.96

# This is how the prior on 'u' looks like:
curve(dgamma(x, a.r, b.r), from=0, to=2, xlab="Per-generation mutation rate, u (x 10^-8)", ylab="Gamma(u | a, b)", las=1, n=5e2)

# File mcmc.txt has the output of running bpp under the A00 analysis
# (fixed species tree, no species delimitation) on 80,662 RADseq
# fragments from six mouse lemur species.
# Make sure R's workding directory is set to 'mousies/R' for this
# command to work:
m1 <- read.table("mcmc.txt", head=TRUE)  # contains 20,000 MCMC samples

names(m1)
# The output contains 5 divergence times (tau's) and 5 ancestral
# scaled population sizes (theta's)

# By using the priors on 'g' and 'u' we will convert the tau's into
# geological divergence times (in years) and the theta's into
# effective population sizes (as numbers of individuals)

# For example, this is how posterior distribution of the root's tau 
# (age of the root in substitutions per site) looks like before
# re-calibrating it to geological time:
plot(density(m1$tau_7OLMXRB), xlim=c(0,0.002), xlab="tau (in substitutions per site)", main="tau_7OLMXRB")

# To obtain the calibrated times and population sizes, we simply
# obtain 20,000 samples from the priors on 'g' and 'u'.

# Recal that the per-yer mutation rate is r=u/g. Thus the calibrated
# times are t = tau/r. Also recall that theta = 4Nu, thus N = theta/(4u).
# So we simply use the sampled values of 'g' and 'u' to recalibrate all
# the tau's and theta's in the m1 dataframe:

n <- nrow(m1)  # 20,000
set.seed(123357)  # We set the set so that the analysis is reproducible
gi <- rgamma(n, a, b)               # sample from prior on 'g'
ri <- rgamma(n, a.r, b.r) * 1e-8    # sample from prior on 'u'

# Column indices for thau's and theta's
tau.i <- 7:11; theta.i <- 2:6

# Obtain population sizes (Ne) and geological times in years (ty):
Ne <- m1[,theta.i] / (4*ri)   # N = theta / (4*u)
ty <- m1[,tau.i] * gi / ri    # t = tau * g / u

# Voilá! Ne and ty contain our posterior estimates of population sizes
# and geological divergence times!

# For example, this is how the posterior distribution of the root's Ne
# and age look like:
plot(density(Ne$theta_7OLMXRB, from=0, to=6e5), xlab = "Root's Ne (number of individuals)", main="Effective size of the root's ancestral population")
plot(density(ty$tau_7OLMXRB/1e3, from=0, to=1.2e3), xlab = "Root's age (thousands of years ago)", main="Root's age")

# Calculate posterior means and 95% credibility-intervals:
N.m <- apply(Ne, 2, mean)
N.95 <- apply(Ne, 2, quantile, prob=c(.025, .975))
t.m <- apply(ty, 2, mean)
t.95 <- apply(ty, 2, quantile, prob=c(.025, .975))

# print out a table of posterior means and CI's for Ne and ty for all
# the populations:
pop.names <- c("7OLMXRB", "8LMXRB", "9LM", "10XRB", "11RB")
Ne.df <- cbind(N.m, t(N.95)); row.names(Ne.df) <-  pop.names
t.df <-  cbind(t.m, t(t.95)); row.names(t.df)  <-  pop.names
Ne.df; t.df

# Note that these results will differ slightly from those in the PNAS
# paper. That's because for the PNAS paper I combined 8 MCMC runs, each
# with 20,000 samples (thus 160,000 samples in total) to calculate the
# posterior. Here we only used samples from one run.