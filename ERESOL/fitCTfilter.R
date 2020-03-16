# fitCTfilter.R
# Maite Ceballos (IFCA)
# 01/2020
#
# script to calculate the constant which fits the central part of an optimal filter
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
rm(list=ls())
display <- "W" # ('P' for pdf or 'W' for window)

library(Hmisc)
# General variables
#-------------------
HEADAS="/home/ceballos/sw/heasoft/x86_64-pc-linux-gnu-libc2.28"
samprate<-156250
pB=75
flen <- c(8192, 4096, 2048)
# limits for different filters
#              8192  4096  2048
init_fit <-  c(1000, 1000, 1000) + pB
fin_fit <-   c(7000, 3500, 1500) + pB
init_repl <- c(1000, 1000, 1000) + pB
fin_repl <-  c(7000, 3500, 1500)

for (i in 1:length(flen)){

    #Read filter in Time Domain
    #--------------------------
    dirfilter <- "/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/xifusimLPA75um/GLOBAL/ADC/"
    filtername <- paste(dirfilter,"library6keV_PL8192_20000p_jitter_bbfb.fits+2",sep="") #FIXFILTT
    #dirfilter <- "/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/xifusimLPA75um/GLOBAL/R/"
    #filtername <- paste(dirfilter,"library6keV_PL8192_20000p_jitter_bbfb_pB75.fits+2",sep="") #FIXFILTT
    
    colname <- paste("T", flen[i],sep="") # filter to be used from the lib
    command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                     "fdump wrap=yes infile=",filtername," columns=",colname, " rows=-",
                     " prhead=no showcol=yes showunit=no showrow=no outfile=filter.txt clobber=yes", sep="")
    system(command)
    filter <- read.table(file='filter.txt', header=TRUE )[,1]
    totalsum <- sum(filter)
    cat("\nWorking with filter length=", flen[i],"\n")
    cat("----------------------------------\n")
    cat("Filter sum=", totalsum, "\n")
    #cat("0-pad sum=", sum(filter[1:4096]), "\n")
    
    
    plot(seq(1,flen[i]),filter[1:flen[i]], col="green",cex=0.2)
    fitFilter <-lm(filter[init_fit[i]:fin_fit[i]]~1)
    summary(fitFilter)
    abline(h=fitFilter$coefficients[1], col="red")
    
    newfilter<-filter
    newfilter[init_repl[i]:fin_repl[i]] <- fitFilter$coefficients[1]
    lines(1:flen[i], newfilter)
    cat("New Filter sum=", sum(newfilter), "\n")
    cat("Fitted constant=", fitFilter$coefficients[1],"\n")
}
