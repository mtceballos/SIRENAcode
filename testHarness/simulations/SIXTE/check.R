#!/usr/bin/Rscript
#
setwd("/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/tessimLPA1")
pulseLength <- 1024
samprate <- 156250
nrmfctr <- read.table("nrmfctr.txt",header=F)[,1]

filter <- read.table("filter1024.txt",header=F)[,1]
pulse <- read.table("template.txt",header=F)[,1]
#P(t) * OF(t) ) * 2*pulse_length/(samp_freq*nrmfctr
energy <- sum(pulse * filter) *2*pulseLength/(samprate*nrmfctr[1])

cat(energy)


