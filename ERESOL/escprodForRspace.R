# escprod.R
# Maite Ceballos (IFCA)
# 03/2019
#
# script to calculate scalar product (data*optimalFilter) at different offsets
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

display <- "W" # ('P' for pdf or 'W' for window)
setwd("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/checkR/")

library(Hmisc)
# General variables
#-------------------
HEADAS="/home/ceballos/sw/heasoft/x86_64-pc-linux-gnu-libc2.28"
outPDFbasename <-"scProd"
energyStr <- "7keV"
#energyStr <- "0.2keV_0.2keV"
samprate<-156250
noffs <- 100 # number of 'offsets' in which scalar product will be calculated

#Read filter in Time Domain
#--------------------------
dirfilterR <- "/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/xifusimLPA75um/GLOBAL/R/"
filtername <- paste(dirfilterR,"library6keV_PL8192_20000p_jitter_bbfb.fits+2",sep="") #FIXFILTT
colname <- "T8192" # filter to be used from the lib
command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                 "fdump wrap=yes infile=",filtername," columns=",colname, " rows=-",
                 " prhead=no showcol=yes showunit=no showrow=no outfile=filterR.txt clobber=yes", sep="")
system(command)
filterR <- read.table(file='filterR.txt', header=TRUE )[,1]

dirfilterADC <- "/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/xifusimLPA75um/GLOBAL/ADC/"
filtername <- paste(dirfilterADC,"library6keV_PL8192_20000p_jitter_bbfb.fits+2",sep="") #FIXFILTT
colname <- "T8192" # filter to be used from the lib
command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                 "fdump wrap=yes infile=",filtername," columns=",colname, " rows=-",
                 " prhead=no showcol=yes showunit=no showrow=no outfile=filterADC.txt clobber=yes", sep="")
system(command)
filterADC <- read.table(file='filterADC.txt', header=TRUE )[,1]

#Read Data Pulse
#--------------------------
dataname <- paste(dirfilterR,"library6keV_PL8192_20000p_jitter_bbfb.fits+1",sep="") #FIXFILTT
colname <- "PULSE"
recordNum <- 1
startSample <-1 
command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                 "fdump wrap=yes infile=",dataname," columns=",colname, " rows=",recordNum,
                 " prhead=no showcol=yes showunit=no showrow=no outfile=pulseR.txt clobber=yes", sep="")
system(command)
dimdata <- 8192
dataTotalR <- read.table(file='pulseR.txt', header=TRUE )[,1]

dataname <- paste(dirfilterADC,"library6keV_PL8192_20000p_jitter_bbfb.fits+1",sep="") #FIXFILTT
colname <- "PULSE"
recordNum <- 1
startSample <-1 
command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                 "fdump wrap=yes infile=",dataname," columns=",colname, " rows=",recordNum,
                 " prhead=no showcol=yes showunit=no showrow=no outfile=pulseADC.txt clobber=yes", sep="")
system(command)
dimdata <- 8192
dataTotalADC <- read.table(file='pulseADC.txt', header=TRUE )[,1]

# do scalar product
finSample <- 8192
dataR <- numeric(8192)
dataR[startSample:finSample] <-dataTotalR[startSample:finSample]
dataADC <- numeric(8192)
dataADC[startSample:finSample] <-dataTotalADC[startSample:finSample]
par(mfrow=c(1,2))
plot(dataADC)
abline(v=256,col="red",lty=2)
plot(dataR)

plot(filterADC)
plot(filterR)



# calculate scalar product in 'noff' offsets
#-------------------------------------------
flen <- length(filter)  #samples
flensecs <- flen/samprate #seconds
scprodR <- numeric(noffs)
reconER <- numeric(noffs)
scprodR <- sum(dataR * filterR)
reconER <- scprodR/flen

scprodADC <- numeric(noffs)
reconEADC <- numeric(noffs)
scprodADC <- sum(dataADC * filterADC)
reconEADC <- scprodADC/flen

cat("Escalar prods(R,ADC)=",scprodR,scprodADC,"\n")
cat("Energ(keV)(R,ADC)=",reconER,reconEADC,"\n")
