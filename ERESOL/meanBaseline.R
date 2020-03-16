# meanBaseline.R
# Maite Ceballos (IFCA)
# 01/2020
#
# script to calculate mean of baseline data in different records
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

setwd("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/")

library(Hmisc)
# General variables
#-------------------
HEADAS="/home/ceballos/sw/heasoft/x86_64-pc-linux-gnu-libc2.28"
energyStr <- "7keV"
samprate<-156250

#Read data
#---------
dirdata <- "/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/xifusimLPA75um/gainScale/"
dataname <- paste(dirdata,"sep40000sam_5000p_",energyStr,"_jitter_bbfb.fits+1",sep="")
colname <- "ADC"

for (nr in 1:50){
    recordNum <- nr
    command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                 "fdump wrap=yes infile=",dataname," columns=",colname, " rows=",recordNum,
                 " prhead=no showcol=yes showunit=no showrow=no outfile=pulse.txt clobber=yes", sep="")
    system(command)
    dataTotal <- read.table(file='pulse.txt', header=TRUE )[,1]
    baseline <- dataTotal[100:900]
    cat("Mean baseline=", mean(baseline), "\n")
}
