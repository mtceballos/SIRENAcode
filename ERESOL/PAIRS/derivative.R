#
# BAGPLOTS: PACKAGE Documented on derivative.ipynb (also GitHub)
#
rm(list=ls())
library(FITSio)
library(aplpack)
source("~/R/Rfunctions/derivPulseFits.r")
#source("~/.Rprofile")
setwd("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/")
# --------------------------------------------------------------------------------------
#  For input: 
#
npulses <- 1000      # Number of pulses at each energy
nPairs <- 500    
samprate <- 156250
#samprateStr="" # to name files with full samprate
#pulseLength<- 4096   # pulse length
samprate <- samprate/2.
samprateStr="_samprate2" # to name files with 1/2 samprate
pulseLength<- 2048   # pulse length

energies <- c("0.2", "0.5", "1", "2", "3", "4", "5", "6", "7", "8") # pulses energies
nenergies <- length(energies)
Esec <- "1" # keV : energy of secondary pulses

pdf(paste("baselineLPA2/derivativePairsStudy_Esec",Esec,"keV",samprateStr,".pdf",sep=""),width=10, height=7,version="1.4")
#energies <- c("0.2", "0.5")
#separations <-c("00005", "00010", "00020", "00045", "00060", "00100", "00200", "00250","00300","00400", "00800")
separations <- c("00002","00005","00010","00020","00022","00030","00045","00050","00060",
                 "00100","00125","00150","00200","00250","00300","00400","00500","00800","01000")
nseps <- length(separations)
separationsToPlot <-c("00200","00250","00300")
separationsToPlot <-c("00100","00125","00200")
seps.ms <- as.numeric(separations)/samprate*1E3 #separations in ms
saveFileFull <- paste("derivateCalibPulsesFull",samprateStr,".dat",sep="")  # File with the array with the derivatives (R object )
saveFileMean4 <- paste("derivateCalibPulsesMean4",samprateStr,".dat",sep="")  # File with the array with the mean 4 derivatives (R object )
pairFileFull <- paste("derivatePairPulsesFull_sec",Esec,"keV",samprateStr,".dat",sep="")  # File with the array with the derivatives (R object ) for pairs
pairFileMean4 <- paste("derivatePairPulsesMean4_sec",Esec,"keV",samprateStr,".dat",sep="")  # File with the array with the mean 4 derivatives (R object ) for pairs
headas <- "/home/ceballos/sw/heasoft/x86_64-unknown-linux-gnu-libc2.24"
repeatCal <- FALSE
# --------------------------------------------------------------------------------------
derivArray <- array(NA,dim=c(npulses,length(energies),pulseLength-1))
derivArrayPair <- array(NA,dim=c(nPairs,length(energies),length(separations),pulseLength-1))
derivArrayMean4samples <- matrix(NA,nrow=npulses, ncol=length(energies))
derivArrayMean4samplesPair <- array(NA, dim=c(nPairs, length(energies), length(separations)))

# Calculate or read DERIVATIVE of SINGLE pulses
if(file.exists(saveFileFull)){
    cat("Looking for data in ",saveFileFull,"\n")
    repeatCal <- FALSE
    load(saveFileFull)
    load(saveFileMean4)
    # check that energies and number of pulses are ok
    if(dim(derivArray)[1] != npulses || 
       dim(derivArray)[2] != length(energies) ||
       dim(derivArray)[3] != (pulseLength-1)  ||
       !identical(colnames(derivArray), energies)) repeatCal <- TRUE
    if(repeatCal) cat("File ", saveFileFull," exists but it is not the required one: it will be calculated again\n")
}
if(!file.exists(saveFileFull) || repeatCal){
    cat("Repeating file ", saveFileFull, " calculation\n")
    
    # save also the mean of the first 4 samples of the derivative
    colnames(derivArray) <- energies
    rownames(derivArray) <- paste("pulse",1:npulses,sep="")
    for (ie in 1:length(energies)){
        cat("Working with energy=",energies[ie],"\n")
        fitsFileLarge <- paste("tessimLPA2shunt/sep40000sam_20000p_",energies[ie],"keV",samprateStr,".fits",sep="")
        fitsFile <- paste("tessimLPA2shunt/sep40000sam_1000p_",energies[ie],"keV",samprateStr,".fits",sep="")
        if(!file.exists(fitsFile)){
            # use larger one
            stopifnot((file.exists(fitsFileLarge) || file.exists(fitsFile)))
            fitsFile <- fitsFileLarge
        }
        istart <- 1000 # Start pulses at sample=istart
        if (energies == "0.2" || energies == "0.5") istart = 999
        derivateListForE <- derivPulseFits(fitsFile,fitsExt = 1, fitsCol = "ADC", npulses = npulses, 
                                           pulseLength = pulseLength,startSample = istart, headas=headas)
        for (ip in 1:npulses){
            derivArray[ip,ie,] <- derivateListForE[[ip]]
            derivArrayMean4samples[ip,ie] <- mean(derivateListForE[[ip]][1:4])
        } # for each pulse
    } # foreach calib energy
    # Save R object (it can be loaded again using "load(saveFile)" where saveFile is the output filename)
    save(derivArray, file=saveFileFull)
    save(derivArrayMean4samples, file=saveFileMean4) 
}

# Get reconstructed energies of SINGLE pulses ---- TAKE CARE of SAMPRATE ----
EkeVrecons <- matrix(NA,nrow=npulses, ncol=length(energies))
BiasEkeVrecons <- matrix(NA,nrow=npulses, ncol=length(energies))
minErecons <- numeric()
maxErecons <- numeric()
for (ie in 1:length(energies)){
    cat("Working with energy=",energies[ie],"\n")
    # get reconstructed energies
    eventsFile <- paste("eresolLPA2shunt/events_sep40000sam_20000p_SIRENA",pulseLength,"_pL",pulseLength,"_",energies[ie],
                        "keV_F0F_fixedlib1OF_OPTFILT",samprateStr,"_NTRIG.fits",sep="")
    zz <- file(description = eventsFile, open = "rb")
    header0 <- readFITSheader(zz, fixHdr = 'remove') # read primary header
    header <- readFITSheader(zz, fixHdr = 'remove') # read extension header
    evtTable <- readFITSbintable(zz, header)
    close(zz)
    idcol <- which(evtTable$colNames == "SIGNAL")
    pulse <- list()
    EkeVrecons[,ie] <- evtTable$col[[idcol]][1:npulses]
    BiasEkeVrecons[,ie] <- EkeVrecons[,ie] - as.numeric(energies[ie])
    minErecons[ie] <- min(EkeVrecons[,ie])
    maxErecons[ie] <- max(EkeVrecons[,ie])
} # foreach calib energy

# Calculate or read DERIVATIVE of PAIRS of pulses
if(file.exists(pairFileFull)){
    cat("Looking for data in ",pairFileFull,"\n")
    repeatcal <- FALSE
    load(pairFileFull)
    load(pairFileMean4)
    # check that energies and number of pulses are ok
    if(dim(derivArrayPair)[1] < nPairs || 
       dim(derivArrayPair)[2] != length(energies) ||
       dim(derivArrayPair)[3] < length(separations)  ||
       dim(derivArrayPair)[4] < (pulseLength-1)  ||
       !identical(colnames(derivArrayPair), energies) ||
       !identical(dimnames(derivArrayPair)[[3]], separations)) repeatCal <- TRUE
    if(repeatCal) cat("File ", pairFileFull," exists but it is not the required one: it will be calculated again\n")
}
if(!file.exists(pairFileFull) || repeatCal){
    cat("Repeating file ", pairFileFull, " calculation\n")
    
    # save also the mean of the first 4 samples of the derivative
    colnames(derivArrayPair) <- energies
    rownames(derivArrayPair) <- paste("pulse",1:nPairs,sep="")
    dimnames(derivArrayPair)[[3]] <- separations
    for (ie in 1:length(energies)){
        for (is in 1:length(separations)){
            cat("Working with Primary energy in Pair=",energies[ie]," at separation ",separations[is],"\n")
            fitsFile <- paste("tessimLPA2shunt/sep",separations[is],"sam_2000p_",
                              energies[ie],"keV_",Esec,"keV",samprateStr,".fits",sep="")
            istart <- 1000 # Start pulses at sample=istart
            if (energies == "0.2" || energies == "0.5") istart = 999
            derivateListForE <- derivPulseFits(fitsFile,fitsExt = 1, fitsCol = "ADC", npulses = nPairs, 
                                           pulseLength = pulseLength,startSample = istart, headas=headas)
            for (ip in 1:nPairs){
                derivArrayPair[ip,ie,is,] <- derivateListForE[[ip]]
                derivArrayMean4samplesPair[ip,ie,is] <- mean(derivateListForE[[ip]][1:4])
            } # for each pulse
        } # foreach separation    
    } # foreach calib energy
    # Save R object (it can be loaded again using "load(pairFileFull)" where pairFileFull is the output filename)
    save(derivArrayPair, file=pairFileFull)
    save(derivArrayMean4samplesPair, file=pairFileMean4) 
}
EkeVreconsPairs <- array(NA, dim=c(nPairs, length(energies), length(separations)))
BiasEkeVreconsPairs <- array(NA, dim=c(nPairs, length(energies), length(separations)))
# Get reconstructed energies of PAIRS of pulses
for (ie in 1:length(energies)){
    for (is in 1:length(separations)){
        if(separations[is] == "00002"){ # too short filter: cannot calculate energy
            EkeVreconsPairs[,ie,is] <- rep(NaN,nPairs)
            BiasEkeVreconsPairs[,ie,is] <- rep(NaN,nPairs)
            next
        }
        # get reconstructed energies
        eventsFile <- paste("eresolLPA2shunt/nodetSP/events_sep",separations[is],"sam_2000p_SIRENA",pulseLength,"_pL",
                            pulseLength,"_",energies[ie],"keV_",Esec,"keV_F0F_fixedlib1OF_OPTFILT",
                            samprateStr,"_NTRIG.fits",sep="")
        zz <- file(description = eventsFile, open = "rb")
        header0 <- readFITSheader(zz, fixHdr = 'remove') # read primary header
        header <- readFITSheader(zz, fixHdr = 'remove') # read extension header
        evtTable <- readFITSbintable(zz, header)
        close(zz)
        idcol <- which(evtTable$colNames == "SIGNAL")
        pulse <- list()
        EkeVreconsPairs[,ie,is] <- evtTable$col[[idcol]][1:nPairs]
        BiasEkeVreconsPairs[,ie,is] <- EkeVreconsPairs[,ie,is] - as.numeric(energies[ie])
    } # foreach separation    
} # foreach calib energy

#
# Draw Bagplots
#
par(mfrow=c(2,3))
for (ie in 1:length(energies)){
    cat("Plotting figure for Eprim=",energies[ie],"keV\n")
    xmin<-min(derivArrayMean4samples[,ie],derivArrayMean4samplesPair[,ie,2:nseps])
    xmax<- max(derivArrayMean4samples[,ie],derivArrayMean4samplesPair[,ie,2:nseps])
    ymin <- min(EkeVrecons[,ie],EkeVreconsPairs[,ie,which(separations %in% separationsToPlot)])
    ymax <- max(EkeVrecons[,ie],EkeVreconsPairs[,ie,which(separations %in% separationsToPlot)])
    #ymin <- min(1E3*BiasEkeVrecons[,ie],1E3*BiasEkeVreconsPairs[,ie,which(separations %in% separationsToPlot)])
    #ymax <- max(1E3*BiasEkeVrecons[,ie],1E3*BiasEkeVreconsPairs[,ie,which(separations %in% separationsToPlot)])
    
    # Yellow-ish bagplot (single pulses)
    bagplot(derivArrayMean4samples[,ie],EkeVrecons[,ie],xlab="<4 derivative samples>", ylab="Reconstructed Energy (keV)",
            main=(bquote(paste(E[prim],"=",.(energies[ie])," keV  ",E[sec],"=",.(Esec)," keV",sep=""))), 
            xlim=c(floor(xmin),ceiling(xmax*1.005)), ylim=c(ymin,ymax),
            col.loophull="cornsilk",col.looppoints="peachpuff",col.baghull="orange",show.outlier = FALSE)
    #bagplot(derivArrayMean4samples[,ie],1E3*BiasEkeVrecons[,ie],xlab="<4 derivative samples>", 
    #        ylab="Bias in Reconstructed Energy (eV)",
    #        main=(bquote(paste(E[prim],"=",.(energies[ie])," keV  ",E[sec],"=",.(Esec)," keV",sep=""))), 
    #        xlim=c(floor(xmin),ceiling(xmax*1.005)), ylim=c(ymin,ymax),
    #        col.loophull="cornsilk",col.looppoints="peachpuff",col.baghull="orange",show.outlier = FALSE)
    # blue-ish bagplots (pairs of pulses)
    for(is in 1:length(separations)){
        if(!separations[is] %in% separationsToPlot) next
        cat("............Plotting bagplot for sep=",separations[is],"samples\n")
        bag <- try(compute.bagplot(derivArrayMean4samplesPair[,ie,is],EkeVreconsPairs[,ie,is]))
        #bag <- try(compute.bagplot(derivArrayMean4samplesPair[,ie,is],1E3*BiasEkeVreconsPairs[,ie,is]))
        if(class(bag) == "try-error") next
        plot.bagplot(bag,add=TRUE, transparency=TRUE)
        #text(xmax*1.001,ymax,"Sep(sam)", cex=0.8 )
        text(xmax*1.001,ymax,"Separation", cex=0.8 )
        #text(xmax*1.001,bag$center[2],separations[is], cex=0.7)
        text(xmax*1.001,bag$center[2],paste(separations[is],"sam/",seps.ms[is],"ms",sep=""), cex=0.7)
    }
}
#lines(colMeans(derivArrayMean4samples), minErecons, col="green")
#lines(colMeans(derivArrayMean4samples), maxErecons, col="red")

#
# Trying to plot 10,10 windows with plots in the diagonal...
#

if(0){

#for (ie in 1:length(energies)){
screens<-list()
nw <- 10
nw1 <-nw+1
xmins <- seq(0.,1.,length.out=nw+1)[1:nw]
xmaxs <- seq(0.,1.,length.out=nw+1)[2:nw1]
ymins <- seq(0.,1.,length.out=nw+1)[1:nw]
ymaxs <- seq(0.,1.,length.out=nw+1)[2:nw1]
screensMat <-rbind()
for (ir in 1:nw){
    for (ic in 1:nw){
        winnum <- (ir-1)*nw+ic
        winname <- paste("sc",winnum,sep="")
        screens[[winname]] <- c(xmins[ic],xmaxs[ic],ymins[ir],ymaxs[ir])
        screensMat <- rbind(screensMat,c(xmins[ic],xmaxs[ic],ymins[ir],ymaxs[ir]))
    }
}
#split.screen(screensMat)
#par(mfrow=c(nenergies,nenergies))
nw <- nenergies
nw <- 5
nw1 <- nw+1
screensToPlot <- screen.num <- 1 + nw1*(1:nw-1)
par(mfrow=c(nw,nw))
for (is in 1:(nw*nw)){
    cat("Window=",is,"\n")
    #screen.num <- 1 + nw1*(ie-1)
    #screen(screen.num)
    #cat("Using screen ",screen.num," in ",screensMat[ie,],"\n")
    xmin<-min(derivArrayMean4samples[,ie],derivArrayMean4samplesPair[,ie,])
    xmax<- max(derivArrayMean4samples[,ie],derivArrayMean4samplesPair[,ie,])
    ymin <- min(EkeVrecons[,ie],EkeVreconsPairs[,ie,which(separations %in% separationsToPlot)])
    ymax <- max(EkeVrecons[,ie],EkeVreconsPairs[,ie,which(separations %in% separationsToPlot)])
    #bagplot(derivArrayMean4samples[,ie],EkeVrecons[,ie],xlab="<4 derivative samples>", ylab="Reconstructed Energy (keV)",
    #        main=paste("Primary pulses E=",energies[ie]," keV",sep=""), 
    #        col.loophull="cornsilk",col.looppoints="peachpuff",col.baghull="orange",show.outlier = FALSE)
    if(is %in% screensToPlot){    
        cat("plotting in is=",is)
        plot(1:is,1:is,cex=0.5)
    }else{
        plot.new()
    }
}
#close.screen(all.screens = TRUE)
}
dev.off()

