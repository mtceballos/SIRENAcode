# plot histogram from reconstructed energies
#===========================================

# READ data
#--------------
rm(list=ls())
library("MASS")
library(stats)
HEADAS="/home/ceballos/sw/heasoft/x86_64-pc-linux-gnu-libc2.28"
setwd("/dataj6/ceballos/INSTRUMEN/EURECA/realData/carlos/20191115_zaragoza_kickoffRTI2018")
dataname <- "evt_PulsosSIRENA_HR.fits+1"
colname <- "SIGNAL"
command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                 "fdump wrap=yes infile=",dataname," columns=",colname, " rows='-' prhead=no ",
                 "showcol=yes showunit=no showrow=no outfile=pulse.txt clobber=yes", sep="")
system(command)
data <- read.table(file='pulse.txt', header=TRUE )[,1]

# Plot data (ka1,ka2)  5.75-5.90 (ka1=5.89875 keV (16.2%); ka2=5.88765 (8.2%))
#-----------------------------------------------------------------------------
dataKa <- data[data>5.75 & data<5.90]
histo<-hist(dataKa, xlab='Reconstructed PH', breaks = 100, xlim=c(5.7,5.9),
            main="Reconstructed PHs")
rug(dataKa)

# fit a single gaussian and calibrate energies (assuming 5.9 keV)
#----------------------------------------------------------------
ftg0 <- fitdistr(dataKa, "normal")
histo<-hist(dataKa, xlab='Reconstructed PH', breaks=100,xlim=c(5.7,5.9), freq = FALSE,
            main="Ka1+Ka2")
meanG1 <- ftg0$estimate[1]
dataKaCalib <- dataKa / ftg0$estimate[1]*5.90

# Use calibrated data & fit single gaussian
#---------------------------------------------
histo<-hist(dataKaCalib, xlab='Reconstructed PH', breaks = 100, xlim=c(5.8,6.0))
rug(dataKaCalib)
ftg <- fitdistr(dataKaCalib, "normal")
histo<-hist(dataKaCalib, xlab='Calibrated energies', breaks=100,xlim=c(5.8,6.0), freq = FALSE,
            main="Ka1+Ka2 calibrated")
meanG1 <- ftg$estimate[1]
sdG1 <- ftg$estimate[2]
fwhmG1 <- sdG1*2.35*1000 #eV
curve(dnorm(x,mean=meanG1, sd=sdG1), add=TRUE)
text(5.85,25, "Single Gaussian")
text(5.85,20, paste("Mean=",sprintf("%.2f",meanG1),"keV"))
text(5.85,16, paste("FWHM=",sprintf("%.2f",fwhmG1), "eV"))


# fit 2 Gaussians
#----------------
x <- histo$mids
y <- histo$density
df <- data.frame(x, y)
fit2G <- nls(y ~ (C1 * exp(-(x-mean1)**2/(2 * sigma1**2)) +
                      C2 * exp(-(x-mean2)**2/(2 * sigma2**2))), data=df,
                      start=list(C1=25, mean1=5.89, sigma1=0.01,
                      C2=30, mean2=5.9, sigma2=0.01), algorithm="port")  
newx <- seq(5.80,6.0,length.out = 100)
fit2G_prediction <- predict(fit2G, newdata=list(x=newx))
lines(newx, fit2G_prediction, col="red")
meanG2_1 <- coef(fit2G)[2] # keV 
fwhmG2_1 <- coef(fit2G)[3]*2.35*1000 #eV
meanG2_2 <- coef(fit2G)[5] # keV
fwhmG2_2 <- coef(fit2G)[6]*2.35*1000 #eV
text(5.98,25, "Double Gaussian", col="red")
text(5.98,20, paste("Mean(Ka1)=",sprintf("%.2f",meanG2_1),"keV"), col="red")
text(5.98,18, paste("FWHM(Ka1)=",sprintf("%.2f",fwhmG2_1), "eV"), col="red")
text(5.98,14, paste("Mean(Ka2)=",sprintf("%.2f",meanG2_2),"keV"), col="red")
text(5.98,12, paste("FWHM(Ka2)=",sprintf("%.2f",fwhmG2_2), "eV"), col="red")


