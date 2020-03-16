# escprod.R
# Maite Ceballos (IFCA)
# 03/2019
#
# script to calculate scalar product (data*optimalFilter) at different offsets
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

display <- "W" # ('P' for pdf or 'W' for window)
setwd("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/Gamil/")

library(Hmisc)
# General variables
#-------------------
HEADAS="/home/ceballos/sw/heasoft/x86_64-pc-linux-gnu-libc2.28"
outPDFbasename <-"scProd"
energyStr <- "7keV"
#energyStr <- "0.2keV_0.2keV"
samprate<-156250
noffs <- 100 # number of 'offsets' in which scalar product will be calculated

#Read data
#---------
dirdata <- "/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/xifusimLPA75um/gainScale/"
dataname <- paste(dirdata,"sep40000sam_5000p_",energyStr,"_jitter_bbfb.fits+1",sep="")
#dataname <- paste(dirdata,"sep00030sam_100p_",energyStr,"_jitter_dcmt100.fits+7",sep="")
colname <- "ADC"
recordNum <- 1
startSample <-950 # pulse is supposed to start ~ sample 1000, so go 50 samples before
finSample <- 10192 #10192 final sample of the stream to be analized (1256 for a 256 filter)
command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                 "fdump wrap=yes infile=",dataname," columns=",colname, " rows=",recordNum,
                 " prhead=no showcol=yes showunit=no showrow=no outfile=pulse.txt clobber=yes", sep="")
system(command)
dimdata <- 10192
dataTotal <- read.table(file='pulse.txt', header=TRUE )[,1]
# data contains prebuffer samples of baseline + pulseFragment + postbuffer
prebuffer <- dataTotal[10:100]
pulse <-dataTotal[startSample:finSample]
postbuffer <- length(dataTotal)-finSample
data <- c(pulse, rep(prebuffer, len=postbuffer))

#Read filter in Time Domain
#--------------------------
dirfilter <- "/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/xifusimLPA75um/GLOBAL/ADC/"
filtername <- paste(dirfilter,"library6keV_PL8192_20000p_jitter_bbfb.fits+2",sep="") #FIXFILTT
#filtername <- paste(dirfilter,"library6keV_PL8192_20000p_jitter_dcmt100.fits+2",sep="") #FIXFILTT
#dirfilter <- "/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/simsDBbbfb/"
#filtername <- paste(dirfilter,"library6keV_PL8192_20000p_bbfb.fits+2",sep="") #FIXFILTT
colname <- "T8192" # filter to be used from the lib
command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                 "fdump wrap=yes infile=",filtername," columns=",colname, " rows=-",
                 " prhead=no showcol=yes showunit=no showrow=no outfile=filter.txt clobber=yes", sep="")
system(command)
filter <- read.table(file='filter.txt', header=TRUE )[,1]
totalsum <- sum(filter)

# calculate scalar product in 'noff' offsets
#-------------------------------------------
flen <- length(filter)  #samples
flensecs <- flen/samprate #seconds
scprod <- numeric(noffs)
reconE <- numeric(noffs)
for (i in 1:noffs){
    ifin <- flen + (i-1)
    scprod[i] <- sum(data[i:ifin] * filter[1:flen])
    reconE[i] <- scprod[i]/flen
    #reconE[i] <- scprod[i]*2*flensecs
    
}
# apply normalising factor since filter in Time Domain comes from Frequency Domain
# check that reconstructed energy is ok (off by a factor due to the energy scale corr)
#--------------------------------------------------------------------------------------
MaxreconE <- max(reconE)
cat("Max escprod=",max(scprod),"\n")
cat("Reconstructed E ~ ", MaxreconE,"\n")
maxoff<-which.max(scprod)
cat("Offset of max prod=",maxoff,"\n")
maxdata<-which.max(data)
cat("Max Pulse data=",maxdata,"\n")
start<-min(which(data>1.1*data[1]))-1
cat("Start Pulse data=",start,"\n")
#plot scalar product
ifin <- flen+ (maxoff-1)

par(mfrow=c(2,1))
plot(seq(1,flen),cumsum(data[maxoff:ifin]*filter[1:flen]), pch=1,cex=0.2, 
     xlim=c(100,9000), ylim=c(5e7,6e7),log="x")
plot(seq(1,flen),data[maxoff:ifin]*filter[1:flen], ylim=c(-2000,1000),pch=1,cex=0.2, xlim=c(7000,8200))
lines(seq(1,flen),data[maxoff:ifin]-1800,col="blue")
lines(seq(1,flen),filter[1:flen]*1000, col="green")
fitFilter <-lm(filter[1000:7000]~seq(1000,7000))
summary(fitFilter)
abline(h=fitFilter$coefficients[1]*1000, col="red")
grid()
#plot(seq(1,flen),data[maxoff:ifin]*filter[1:flen], ylim=c(-2000,1000),xlim=c(100,9000), log="x",pch=1,cex=0.2)
#lines(seq(1,flen),data[maxoff:ifin]-1800,col="blue")
#lines(seq(1,flen),filter[1:flen]*1000, col="green")


# Plot Pulse/data vs samples  && Re-scaled Scalar product vs. data samples
#--------------------------------------------------------------------------
plotlen<- 2*flen
#plotlen <- finSample+100

maxy <- max(data)
for (i in 45:noffs){
    pdffile=paste(outPDFbasename,"_",energyStr,"_filter",flen,"_",sprintf("%02d",i),".pdf",sep="")
    if(display == "P" || display == "p"){
        pdf(pdffile, width=10, height=7)
        cat("Creating pdf:", pdffile,"\n")
    }
    par(mfrow=c(1,2))
    # fig1: pulse data & filter
    plot(1:plotlen,data[1:plotlen], col="blue",ylim=c(200,maxy),cex=0.5, xlab="Samples", 
         ylab="Pulse data and re-scaled filter", cex.lab=1.2)
    
    # fig1: filter (rescaled)
    lines(i:(i+(plotlen-1)),filter[1:plotlen]/max(filter)*maxy, col="red")
    abline(v=start,lty=2,col="darkgray",lw=3)
    abline(v=i,lty=2,col="red",lw=1)
    legend("topright",c(paste(energyStr," pulse",sep=""),"Rescaled 6keV filter", 
                        "Pulse Start"), pch=c(19,NA,NA),lty=c(NA,1,2),
           col=c("blue","red","darkgray"), bty="n",cex=0.8)
    
    # fig2: scalar product vs data samples 
    plot(1:noffs,reconE[1:noffs]/max(reconE[1:noffs])*maxy, 
           pch=19, cex=0.5, ylim=c(800,maxy),cex.lab=1.2,
            xlab="Samples", ylab="Scalar product (normalized to data)")
    # plot diamond showing each value of scalar product as filter runs over data
    points(i,reconE[i]/max(reconE[1:noffs])*maxy, pch=5, col="magenta",cex=2)
    E<-sprintf("%.2f",reconE[i])
    text(i-20,reconE[i]/max(reconE[1:noffs])*maxy, paste("E=",E,sep=""))
    # plot start of data pulse
    points(1:length(data),data,col="blue", cex=0.5)
    lines(1:length(data),data,col="blue")
    # plot line marking maximum value of product
    abline(v=start,lty=2,col="darkgray",lw=3)
    legend("bottomright",c("Scalar product",paste(energyStr, " pulse", sep=""), 
                        "Pulse start"), pch=c(19,NA,NA),lty=c(NA,1,2),
           col=c("black","blue","darkgray"), bty="n",cex=0.8)
    par(mgp=c(1,0.5,0))
    subplot(plot(filter,col="red", type="l",
                 xlab="Sample", ylab="Filter",cex.lab=0.7, cex.axis=0.7),
            size=c(0.8,0.8), 20, maxy/4.)
    par(mgp=c(3,1,0))
    
    if(readline(prompt = "Press <Enter> to continue...(q to quit):") == "q") {
        dev.off()
        break
    }
    if(display == "P" || display == "p"){
        cat("Closing PDF dislay\n")
        dev.off()
    }
}
