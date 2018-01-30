#
# BAGPLOTS: PACKAGE Documented on derivative.ipynb (also GitHub)
#
rm(list=ls())
library(FITSio)
library(aplpack)
setwd("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/")
# --------------------------------------------------------------------------------------
#  For input: 
#
npulses <- 1000      # Number of pulses at each energy
nPairs <- 500    
filterLengths <- c(4096,512,256)
samprate <- 156250
samprateStr="" # to name files with full samprate
pulseLength<- 4096   # pulse length
fEnergy="6" #keV
#samprate <- samprate/2.
#samprateStr="_samprate2" # to name files with 1/2 samprate
#pulseLength<- 2048   # pulse length

energies <- c("0.2", "0.5", "1", "2", "3", "4", "5", "6", "7", "8") # pulses energies
nenergies <- length(energies)
Esec <- "0.2" # keV : energy of secondary pulses

pdf(paste("baselineLPA2/derivativePairsStudy_Esec",Esec,"keV",samprateStr,".pdf",sep=""),width=10, height=7,version="1.4")
#energies <- c("0.2", "0.5")

separations <-c("00090","00100","00110","00120","00130","00140","00150",
                "00250","00260","00270","00280","00290",
                "00300","00310","00320","00330","00340","00350","00360",
                "00500","00510")

nseps <- length(separations)
seps.ms <- as.numeric(separations)/samprate*1E3 #separations in ms

# --------------------------------------------------------------------------------------

derivArrayMean4samples <- array(data=NA,dim=c(npulses,length(energies),length(filterLengths)))
derivArrayMean4samplesPair <- array(data=NA, dim=c(nPairs, length(energies), length(separations),length(filterLengths)))

#=====================================================================
#                                                                    #
#                              SINGLE PULSES                         #
#                                                                    #
#=====================================================================


# Get reconstructed energies AND <4-samples> derivative of SINGLE pulses
#=======================================================================
EkeVrecons <- array(data=NA,dim=c(npulses,length(energies),length(filterLengths)))
BiasEkeVrecons <- array(data=NA,dim=c(npulses,length(energies),length(filterLengths)))
minErecons <- matrix(NA,nrow = length(energies), ncol = length(filterLengths))
maxErecons <- matrix(NA,nrow = length(energies), ncol = length(filterLengths))

for (ie in 1:length(energies)){
    cat("Working with SINGLES energy=",energies[ie],"\n")
    for (ifl in 1:length(filterLengths)){
        fl <- filterLengths[ifl]
        # get reconstructed energies
        eventsFile <- paste("eresolLPA2shunt/nodetSP/events_sep40000sam_20000p_SIRENA",pulseLength,"_pL",pulseLength,
                            "_",energies[ie],"keV_NTRIG_F0F_fixedlib",fEnergy,"OF_OPTFILT",fl,samprateStr,".fits",sep="")
        zz <- file(description = eventsFile, open = "rb")
        header0 <- readFITSheader(zz, fixHdr = 'remove') # read primary header
        header <- readFITSheader(zz, fixHdr = 'remove') # read extension header
        evtTable <- readFITSbintable(zz, header)
        close(zz)
        idcol <- which(evtTable$colNames == "SIGNAL")
        EkeVrecons[,ie,ifl] <- evtTable$col[[idcol]][1:npulses]
        BiasEkeVrecons[,ie,ifl] <- EkeVrecons[,ie,ifl] - as.numeric(energies[ie])
        minErecons[ie,ifl] <- min(EkeVrecons[,ie,ifl])
        maxErecons[ie,ifl] <- max(EkeVrecons[,ie,ifl])
        # get derivative
        idcol <- which(evtTable$colNames == "AVG4SD")
        derivArrayMean4samples[,ie,ifl] <-  evtTable$col[[idcol]][1:npulses]
    } # foreach filter length
} # foreach calib energy
cat("SINGLE pulses (yellow bagplot) finished\n")
#=====================================================================
#                                                                    #
#                              PAIR PULSES                           #
#                                                                    #
#=====================================================================

# Get reconstructed energies AND <4-samples> derivative of PAIR pulses
#=======================================================================
EkeVreconsPairs <- array(NA, dim=c(nPairs, length(energies), length(separations),length(filterLengths)))
BiasEkeVreconsPairs <- array(NA, dim=c(nPairs, length(energies), length(separations),length(filterLengths)))
# Get reconstructed energies of PAIRS of pulses
for (ie in 1:length(energies)){
    cat("Working with PAIRS: EnergyPrim=",energies[ie],"EnergySec=",Esec,"\n")
    for (is in 1:length(separations)){
        for (ifl in 1:length(filterLengths)){
            fl <- filterLengths[ifl]
            if(separations[is] == "00002"){ # too short filter: cannot calculate energy
                EkeVreconsPairs[,ie,is,ifl] <- rep(NaN,nPairs)
                BiasEkeVreconsPairs[,ie,is,ifl] <- rep(NaN,nPairs)
                next
            }
            # get reconstructed energies
            eventsFile <- paste("eresolLPA2shunt/nodetSP/events_sep",separations[is],"sam_2000p_SIRENA",pulseLength,"_pL",
                                pulseLength,"_",energies[ie],"keV_",Esec,"keV_NTRIG_F0F_fixedlib",fEnergy,"OF_OPTFILT",fl,
                                samprateStr,".fits",sep="")
            zz <- file(description = eventsFile, open = "rb")
            header0 <- readFITSheader(zz, fixHdr = 'remove') # read primary header
            header <- readFITSheader(zz, fixHdr = 'remove') # read extension header
            evtTable <- readFITSbintable(zz, header)
            close(zz)
            idcol <- which(evtTable$colNames == "SIGNAL")
            EkeVreconsPairs[,ie,is,ifl] <- evtTable$col[[idcol]][1:nPairs]
            BiasEkeVreconsPairs[,ie,is,ifl] <- EkeVreconsPairs[,ie,is,ifl] - as.numeric(energies[ie])
            # get derivative
            idcol <- which(evtTable$colNames == "AVG4SD")
            derivArrayMean4samplesPair[,ie,is,ifl] <- evtTable$col[[idcol]][1:nPairs]
            
        } # foreach filter length    
    } # foreach separation    
} # foreach calib energy
cat("PAIR pulses (yellow bagplot) finished\n")

#
# Draw Bagplots
#
par(mfrow=c(2,3))
# DEFINE array of separations to Plot (min,max)
#                                            Eprim    Esec           flength      min max
separationsToPlot <- array(data=NA, dim=c(nenergies,nenergies,length(filterLengths),2), 
                           dimnames = list(energies,energies,filterLengths, c("min","max")))
# Esec=0.2 keV
separationsToPlot[    ,1,1,1]  <-"00250"
separationsToPlot[1:6 ,1,1,2]  <-"00360"
separationsToPlot[7:10,1,1,2] <-"00330"
separationsToPlot[    ,1,2,1] <-"00075"
separationsToPlot[    ,1,2,2] <-"00200"
separationsToPlot[    ,1,3,1] <-"00039"
separationsToPlot[    ,1,3,2] <-"00260"
# Esec=0.5 keV
separationsToPlot[    ,2,1,1] <-"00260"
separationsToPlot[1:8 ,2,1,2] <-"00310"
separationsToPlot[9:10,2,1,2] <-"00300"
separationsToPlot[    ,2,2,1] <-"00500"
separationsToPlot[    ,2,2,2] <-"00520"
separationsToPlot[    ,2,3,1] <-"00220"
separationsToPlot[    ,2,3,2] <-"00260"
# Esec=1 keV
separationsToPlot[,3,1,1] <-"00270"
separationsToPlot[,3,1,2] <-"00300"
separationsToPlot[,3,2,1] <-"00500"
separationsToPlot[,3,2,2] <-"00520"
separationsToPlot[,3,3,1] <-"00220"
separationsToPlot[,3,3,2] <-"00260"
# Esec=2 keV
separationsToPlot[,4,1,1] <-"00270"
separationsToPlot[,4,1,2] <-"00290"
separationsToPlot[,4,2,1] <-"00500"
separationsToPlot[,4,2,2] <-"00520"
separationsToPlot[,4,3,1] <-"00220"
separationsToPlot[,4,3,2] <-"00260"
# Esec=3 keV
separationsToPlot[,5,1,1] <-"00270"
separationsToPlot[,5,1,2] <-"00290"
separationsToPlot[,5,2,1] <-"00510"
separationsToPlot[,5,2,2] <-"00520"
separationsToPlot[,5,3,1] <-"00230"
separationsToPlot[,5,3,2] <-"00260"
# Esec=4 keV
separationsToPlot[    ,6,1,1] <-"00270"
separationsToPlot[    ,6,1,2] <-"00290"
separationsToPlot[1:3 ,6,2,1] <-"00500"
separationsToPlot[4:10,6,2,1] <-"00510"
separationsToPlot[    ,6,2,2] <-"00520"
separationsToPlot[1   ,6,3,1] <-"00220"
separationsToPlot[2:10,6,3,1] <-"00230"
separationsToPlot[    ,6,3,2] <-"00260"

# Esec=5 keV
separationsToPlot[    ,7,1,1] <-"00270"
separationsToPlot[    ,7,1,2] <-"00280"
separationsToPlot[    ,7,2,1] <-"00510"
separationsToPlot[    ,7,2,2] <-"00520"
separationsToPlot[1   ,7,3,1] <-"00220"
separationsToPlot[2:10,7,3,1] <-"00230"
separationsToPlot[    ,7,3,2] <-"00260"

# Esec=6 keV
separationsToPlot[    ,8,1,1] <-"00270"
separationsToPlot[    ,8,1,2] <-"00280"
separationsToPlot[    ,8,2,1] <-"00510"
separationsToPlot[    ,8,2,2] <-"00520"
separationsToPlot[1   ,8,3,1] <-"00220"
separationsToPlot[2:10,8,3,1] <-"00230"
separationsToPlot[    ,8,3,2] <-"00260"

# Esec=7 keV
separationsToPlot[1:8 ,9,1,1] <-"00270"
separationsToPlot[9:10,9,1,1] <-"00260"
separationsToPlot[    ,9,1,2] <-"00280"
separationsToPlot[    ,9,2,1] <-"00510"
separationsToPlot[    ,9,2,2] <-"00520"
separationsToPlot[1   ,9,3,1] <-"00220"
separationsToPlot[2:10,9,3,1] <-"00230"
separationsToPlot[    ,9,3,2] <-"00260"

# Esec=8 keV
separationsToPlot[    ,10,1,1] <-"00260"
separationsToPlot[    ,10,1,2] <-"00280"
separationsToPlot[    ,10,2,1] <-"00510"
separationsToPlot[    ,10,2,2] <-"00520"
separationsToPlot[1   ,10,3,1] <-"00220"
separationsToPlot[2:10,10,3,1] <-"00230"
separationsToPlot[    ,10,3,2] <-"00260"

indexEsec <- which(energies == Esec)
for (ie in 1:length(energies)){
    for (ifl in 1:length(filterLengths)){
        fLength=filterLengths[ifl]
        minSep <- as.numeric(separationsToPlot[ie,indexEsec,ifl,1])
        maxSep <- as.numeric(separationsToPlot[ie,indexEsec,ifl,2])
        cat("Plotting figure for Eprim=",energies[ie],"keV, Esec=",Esec, "and Filter Length=",fLength,"\n")
        xmin<-min(derivArrayMean4samples[,ie,ifl],derivArrayMean4samplesPair[,ie,,ifl])
        xmax<- max(derivArrayMean4samples[,ie,ifl],derivArrayMean4samplesPair[,ie,,ifl])
        ymin <- min(EkeVrecons[,ie,ifl],EkeVreconsPairs[,ie,which(as.numeric(separations) %in% minSep:maxSep),ifl])
        ymax <- max(EkeVrecons[,ie,ifl],EkeVreconsPairs[,ie,which(as.numeric(separations) %in% minSep:maxSep),ifl])
        
        # Yellow-ish bagplot (single pulses)
        bagplot(derivArrayMean4samples[,ie,ifl],EkeVrecons[,ie,ifl],xlab="<4 derivative samples>", 
                ylab="Reconstructed Energy (keV)", 
                main=(bquote(paste(E[prim],"=",.(energies[ie])," keV  ",
                                   E[sec],"=",.(Esec)," keV  ", "FLength=",.(fLength),
                                   "  FEnergy=",.(fEnergy)," keV",sep=""))), cex.main=0.8,
                xlim=c(floor(xmin),ceiling(xmax*1.005)), ylim=c(ymin,ymax),
                col.loophull="cornsilk",col.looppoints="peachpuff",col.baghull="orange",show.outlier = FALSE)
                # blue-ish bagplots (pairs of pulses)
        for(is in 1:length(separations)){
            if(!as.numeric(separations[is]) %in% minSep:maxSep) next
            cat("............Plotting bagplot for sep=",separations[is],"samples\n")
            bag <- try(compute.bagplot(derivArrayMean4samplesPair[,ie,is,ifl],EkeVreconsPairs[,ie,is,ifl]))
    
            if(class(bag) == "try-error") next
            plot.bagplot(bag,add=TRUE, transparency=TRUE)
    
            text(xmax,ymax,"Separation(sam)", cex=0.8 )
    
            #text(xmax*1.001,bag$center[2],paste(separations[is],"sam/",seps.ms[is],"ms",sep=""), cex=0.7)
            text(xmax*1.001,bag$center[2],paste(separations[is],"sam",sep=""), cex=0.7)
        } #separations
        abline(h=max(EkeVrecons[,ie,ifl]), col="orange",lty=2)
        abline(h=min(EkeVrecons[,ie,ifl]), col="orange",lty=2)
        
    } # filter lengths
}
dev.off()

