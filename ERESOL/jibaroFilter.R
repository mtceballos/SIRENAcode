# jibaroFilter.R
# Maite Ceballos (IFCA)
# 01/2020
#
# script to calculate reduced filters using as base filter the largest filter 8192
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
rm(list=ls())
display <- "W" # ('P' for pdf or 'W' for window)
library("RColorBrewer")
library(FITSio)
library(Hmisc)
library(reticulate)

# General variables
#-------------------
HEADAS="/home/ceballos/sw/heasoft/x86_64-pc-linux-gnu-libc2.28"
samprate<-156250
fmaxlen=8192

# Read pulse data
# ------------------
dirdata <- "/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/xifusimLPA75um/gainScale/"
dataname <- paste(dirdata,"sep40000sam_5000p_6keV_jitter_bbfb_xifusim040.fits+1",sep="")

colname <- "ADC"
recordNum <- 1
command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                 "fdump wrap=yes infile=",dataname," columns=",colname, " rows=",recordNum,
                 " prhead=no showcol=yes showunit=no showrow=no outfile=pulse.txt clobber=yes", sep="")
system(command)
dataTotal <- read.table(file='pulse.txt', header=TRUE )[,1]


#Read (original) filters in Time Domain
#----------------------------------------
dirfilter <- "/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/xifusimLPA75um/GLOBAL/ADC/"
#filterlib <- paste(dirfilter,"library6keV_PL8192_20000p_jitter_bbfb_040.fits+2",sep="") #FIXFILTT
filterlib <- paste(dirfilter,"library6keV_PL8192_20000p_jitter_bbfb.fits+2",sep="") #FIXFILTT
columns <- c("T8192","T4096","T2048","T1024","T512","T256","T128")
ncols <- length(columns)

filters <- array(data=NA,dim=c(ncols,fmaxlen)) # FWHM of Erecons
for (i in 1:ncols){
    colname <- columns[i]
    command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                 "fdump wrap=yes infile=",filterlib," columns=",colname, " rows=-",
                 " prhead=no showcol=yes showunit=no showrow=no outfile=filter.txt clobber=yes", sep="")
    system(command)
    ff <- read.table(file='filter.txt', header=TRUE )[,1]
    filters[i,1:length(ff)] <- ff
    
}
# Plot filters
#-------------
x <- 1:fmaxlen
plot(x, 1:fmaxlen, type="n", ylim=c(-50,50), xlab="Samples", ylab="Filter value", 
     main="Filters in Time Domain")

colors <- brewer.pal(n = ncols, name = "Set2")
for (i in 1:ncols){
    lines(1:fmaxlen, filters[i,], col=colors[i])
}
legend("topright",legend=columns, col=colors, lty=rep(1,ncols), cex=0.7, bty="n")

# plot filters (rescaled) differences
#-------------------------------------
plot(x, x, type="n", xlim=c(1,1024), ylim=c(-1,0.5), xlab="Samples", ylab="Filter-Filter8192 (rescaled)",
     main="Re-scales differences of filters")
lines(x, filters[5,]-filters[1,]*max(filters[5,],na.rm=TRUE)/max(filters[1,],na.rm=TRUE), col=colors[4]) #1024
lines(x, filters[4,]-filters[1,]*max(filters[4,],na.rm=TRUE)/max(filters[1,],na.rm=TRUE), col=colors[3]) #1024
lines(x, filters[3,]-filters[1,]*max(filters[3,],na.rm=TRUE)/max(filters[1,],na.rm=TRUE), col=colors[2]) #2048
lines(x, filters[2,]-filters[1,]*max(filters[2,],na.rm=TRUE)/max(filters[1,],na.rm=TRUE), col=colors[1]) #4096
abline(h=0, col="grey", lty=2)
lines(x, dataTotal[1001:(length(dataTotal)-1000)]/(2.*max(dataTotal)))
legend("topright",legend=paste(columns[2:ncols],"-T8192 resc",sep=""), col=colors, lty=rep(1,ncols), cex=0.7, bty="n")

# plot filters[i] vs filters[1] and initial points (black crosses, used to fit relations)
#-----------------------------------------------------------------------------------------
plot(filters[1,], filters[2,], col=colors[2], xlim=c(-30,20), ylim=c(-10,10),lty=1, cex=0.8,
     xlab="Filter8192", ylab="Shorter filter") 
points(filters[1,], filters[3,], col=colors[3]) #2048 vs 8192
points(filters[1,], filters[4,], col=colors[4]) #1024 vs 8192
points(filters[1,], filters[5,], col=colors[5]) #512 vs 8192
points(filters[1,], filters[6,], col=colors[6]) #256 vs 8192
points(filters[1,], filters[7,], col=colors[7]) #128 vs 8192
#points(filters[1,1:1000], filters[2,1:1000], col="black", pch=3) #4096 vs 8192
#points(filters[1,1:500], filters[3,1:500], col="black", pch=3) #2048 vs 8192
#points(filters[1,1:500], filters[4,1:500], col="black", pch=3) #1024 vs 8192
points(filters[1,2:200], filters[5,2:200], col="black", pch=3) #512 vs 8192
#points(filters[1,13:40], filters[6,13:40], col="black", pch=3) #256 vs 8192
#points(filters[1,13:30], filters[7,13:30], col="black", pch=3) #128 vs 8192
abline(h=0, col="grey", lty=2)
abline(v=0, col="grey", lty=2)
legend("topleft",legend=c(columns[2:ncols], "Initial pts used for fit"), 
       col=c(colors[2:ncols], "black"), pch=c(rep(1,ncols-1),3), cex=0.7, bty="n")







# PLOT fitted filters vs original filters and 
#         RECONSTRUCT LIBRARY REPLACING ORIGINAL FILTER BY FITTED FILTER (optional)
#=====================================================================================
createNewLib <- 0
if (createNewLib){
    newfilterlib <- paste(dirfilter,"library6keV_PL8192_20000p_jitter_bbfb_fit.fits",sep="")
}

#           8192 4096 2048 1024 512 256 128
initfit <- c(NA,   1,   1,   1,   2, 13, 13 )
finfit <-  c(NA, 500, 500, 500, 200, 40, 30 )
tailsam <- c(NA, 200, 200, 200, 100, 0, 0)

a0 <- numeric(ncols)
b <- numeric(ncols)
flen <- numeric(ncols)
#filterFIT <- array(data=NA,dim=c(ncols,8192)) # FWHM of Erecons
for (i in 2:ncols){
    flen[i] <- as.numeric(substr(columns[i],2,5))
    cat("Fitting filter of ", flen[i], "samples\n")
    # fit initial points: filters[i] ~ filters[1]
    fitLine <- lm(filters[i,initfit[i]:finfit[i]]~filters[1,initfit[i]:finfit[i]])
    summary(fitLine) 
    a0[i] <- fitLine$coefficients[1]
    b[i] <- fitLine$coefficients[2]
    
    # add fitted tail to fitted filter
    filterFIT <- a0[i] + b[i]*filters[1,1:flen[i]]
    tailsam2 <- flen[i] - tailsam[i]+1
    if (tailsam[i] > 0){
        filterFIT[tailsam2:flen[i]] <- a0[i] + b[i]*tail(filters[1,], n=(tailsam[i]))
    }
    cat("Sum(fitted filter for", flen[i], ") =",sum(filterFIT),"\n")
    
    # plot original filter vs fitted filter
    #plot(x, x, type="n", xlim=c(1,flen[i]), ylim=range(filters[i,], na.rm=TRUE),
    #     xlab="Samples", ylab="Filter value", main="Filter (fitted and original)")
    plot(x, x, type="n", xlim=c(1,flen[i]), ylim=range(filterFIT, na.rm=TRUE),
         xlab="Samples", ylab="Filter value", main="Filter (fitted and original)")
    lines(x, filters[i,], col="black")
    lines(x[1:flen[i]], filterFIT, col=colors[i])
    points(which.min(filters[i,]), min(filters[i,], na.rm=TRUE), pch=3, col="black")
    points(which.min(filterFIT), min(filterFIT), pch=4, col=colors[i])
    points(flen[i],filters[i,flen[i]], pch=3, col="black")
    points(flen[i],filterFIT[flen[i]], pch=4, col=colors[i])
    legend("bottomright", legend=c(paste("Original filter",flen[i]), paste("fitted filter",flen[i])),
           bty="n", col=c("black", colors[i]), lty=c(1,1))
    
    if (createNewLib){
        write.matrix(filterFIT, "pp.csv")
        #pycom <- paste("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/set_fit_filter.py --infile=",
        #               newfilterlib, " --ext=2 --colname='", columns[i], "' --arrayfile=pp.csv", 
        #               " --init=1 --fin=", flen[i], sep="")
        #source_python(pycom)
        source_python("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/set_fit_filter.py")
        replaceCol(infile=newfilterlib, ext=2, colname=columns[i], 
                arrayfile="pp.csv", init=1, fin=flen[i])
    }
    if(readline(prompt = "Press <Enter> to continue...(q to quit):") == "q") {
        break
    }
}

# Plot fitting coefficients (a0,b)
# --------------------------------
plot(flen[2:ncols],a0[2:ncols], xlab="Filter number", ylab="a0")
plot(flen[2:ncols],b[2:ncols], xlab="Filter number", ylab="b")
                       
