########################################################################################
#
# Calculate derivative of single (sep=40000 samples) pulses
#
########################################################################################
#rm(list = ls())
library(FITSio)
#source("~/Rfunctions/DerivPulseFits.r")
#source("~/.Rprofile")
setwd("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/")

# --------------------------------------------------------------------------------------
#  For input: 
#
npulses <- 1000      # Number of pulses at each energy
pulseLength<- 4096   # pulse length
energies <- c("0.2", "0.5", "1", "2", "3", "4", "5", "6", "7", "8") # pulses energies
saveFile <- "derivateCalibPulses.dat"  # File with the array with the derivatives (R object )
headas <- "/home/ceballos/sw/heasoft/x86_64-unknown-linux-gnu-libc2.24"
# --------------------------------------------------------------------------------------

# Create output array :
#                    Ener1                      Ener2                 ...       
#   pulse1   p1E1Der[1:pulseLength-1]  p1E2Der[1:pulseLength-1]       ... 
#   pulse2   p2E1Der[1:pulseLength-1]  p2E2Der[1:pulseLength-1]       ... 
#     ..
#   pulseN   pNE1Der[1:pulseLength-1]  pNE2Der[1:pulseLength-1]       ... 
#
#  i.e: derivArray[2,6,] is the vector with the derivative of the 2nd pulse at 6th energy
#  i.e: derivArray[2,6,1] is the 1st sample of the derivative of the 2nd pulse at 6th energy

derivArray <- array(NA,dim=c(npulses,length(energies),pulseLength-1))
# save also the mean of the first 4 samples of the derivative
derivArrayMean4samples <- matrix(NA,nrow=npulses, ncol=length(energies))
EkeVrecons <- matrix(NA,nrow=npulses, ncol=length(energies))

plot(1:npulses, 1:npulses, xlim=c(0,1500),ylim=c(0,6), type = "n", xlab="Mean of 4 samples in derivative",
     ylab="Reconstructed Energy (keV)")

colnames(derivArray) <- energies
rownames(derivArray) <- paste("pulse",1:npulses,sep="")

minErecons <- numeric()
maxErecons <- numeric()
for (ie in 1:length(energies)){
    cat("Working with energy=",energies[ie],"\n")
    fitsFileLarge <- paste("tessimLPA2shunt/sep40000sam_20000p_",energies[ie],"keV.fits",sep="")
    fitsFile <- paste("tessimLPA2shunt/sep40000sam_1000p_",energies[ie],"keV.fits",sep="")
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
    }
    
    # get reconstructed energies
    eventsFile <- paste("eresolLPA2shunt/events_sep40000sam_20000p_SIRENA4096_pL4096_",energies[ie],
                        "keV_F0F_fixedlib1OF_OPTFILT_NTRIG.fits",sep="")
    zz <- file(description = eventsFile, open = "rb")
    header0 <- readFITSheader(zz, fixHdr = 'remove') # read primary header
    header <- readFITSheader(zz, fixHdr = 'remove') # read extension header
    evtTable <- readFITSbintable(zz, header)
    close(zz)
    idcol <- which(evtTable$colNames == "SIGNAL")
    pulse <- list()
    EkeVrecons[,ie] <- evtTable$col[[idcol]][1:npulses]
    minErecons[ie] <- min(EkeVrecons)
    maxErecons[ie] <- max(EkeVrecons)
    
    points(derivArrayMean4samples[,ie], EkeVrecons[,ie], col="red", pch=20,cex=0.1)
}
lines(colMeans(derivArrayMean4samples), minErecons, col="green")
lines(colMeans(derivArrayMean4samples), maxErecons, col="red")

# Save R object (it can be loaded again using "load(saveFile)" where saveFile is the output filename)
save(derivArray, file=saveFile) 



