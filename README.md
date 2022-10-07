# Calibrating coalescent trees to geological time
(Note this tutorial is superseded by https://dosreislab.github.io/2018/08/31/bppr-times.html)

A tutorial on calibrating coalescent trees to geological divergence times using next-gen estimates of the per-generation mutation rate.  The tutorial uses as example a mouse lemur phylogeny estimated under the multispecies coalescent using RADseq data. The phylogeny is then calibrated to geological time by using next-gen estimates of de novo mutation rate in humans and mice, and estimates of mouse lemur generation times. A detailed explanation is given in https://dosreislab.github.io/2016/07/03/coalescent-calibration.html.

File `R/mcmc.txt` contains an MCMC sample (obtained with the program BPP) from the posterior distribution of tau (coalescent branch lengths) and theta (ancestral scaled population sizes) for a fixed tree of six mouse lemur species.

File `R/analysis.R` contains an R script that uses `R/mcmc.txt` as input to convert the coalescent branch lengths (tau) into geological times of divergence, and to convert the scaled population sizes (theta) into actual effective population sizes (i.e., in numbe of individuals).

Directory `bpp/` contains the necessary files to run an A00 analysis with BPP to obtain the MCMC sample for tau and theta.

This reproduces the divergence time analysis in Yoder et al. (PNAS, 113: 8049)

***

![](mousies.png)
