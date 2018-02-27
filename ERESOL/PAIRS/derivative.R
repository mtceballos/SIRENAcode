#
# BAGPLOTS: PACKAGE Documented on derivative.ipynb (also GitHub)
#
#rm(list=ls())
library(FITSio)
library(aplpack)
library(Hmisc)
setwd("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/")
# --------------------------------------------------------------------------------------
#  For input: 
#
npulses <- 1000      # Number of pulses at each energy
nPairs <- 500    
samprate <- 156250
samprateStr="" # to name files with full samprate
samprateStr="_samprate2" # to name files with 1/2 samprate
if (samprateStr == "_samprate2"){
    samprate <- samprate/2.    
    filterLengths <- c(4096,256,128)
    pulseLength<- 4096   # pulse length
    xmax <- c(400, 400, 150 )
    separations <- sprintf("%05d",
        sort(c(seq(15,35,5),seq(40,190,2),200,250,300,400,500)))
}else{
    filterLengths <- c(8192,512,256)
    pulseLength<- 8192   # pulse length
    xmax <- c(800, 800, 300 )
    separations <- sprintf("%05d",
        sort(c(seq(30,300,5),52,57,62,seq(310,500,5),510,seq(520,800,20),900,1000,2000)))
}
fEnergy="6" #keV
energies <- c("0.2", "0.5", "1", "2", "3", "4", "5", "6", "7", "8") # pulses energies
nenergies <- length(energies)
Esec <- "0.2" # keV : energy of secondary pulses
nseps <- length(separations)
seps.ms <- as.numeric(separations)/samprate*1E3 #separations in ms

pdf(paste("baselineLPA2/derivativePairsStudy_Esec",Esec,"keV",samprateStr,".pdf",sep=""),width=10, height=7,version="1.4")

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
            header0 <- readFITSheader(zz, fixHdr = 'none') # read primary header
            header <- readFITSheader(zz, fixHdr = 'none') # read extension header
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

par(mfrow=c(2,3))
#
# DEFINE array of separations to Plot (min,max)
#                                            Eprim    Esec           flength      min max
separationsToPlot <- array(data="", dim=c(nenergies,length(filterLengths)), 
                           dimnames = list(energies,filterLengths))
#
# PLOTTING
#
indexEsec <- which(energies == Esec)
lcols <- rainbow(length(energies))

for (ie in 1:length(energies)){
    #
    # Plot variation of Reconstructed energy with separation for every filter length
    #
    bagplots.separations <- list()
    
    for (ifl in 1:length(filterLengths)){
        fLength=filterLengths[ifl]
        meanErecon <- numeric(length(separations))
        maxErecon <- numeric(length(separations))
        minErecon <- numeric(length(separations))
        # define orange band
        bandMin <- min(EkeVrecons[,ie,ifl])
        bandMax <- max(EkeVrecons[,ie,ifl])
        
        for(is in 1:length(separations)){
            meanErecon[is] <- mean(EkeVreconsPairs[,ie,is,ifl])
            maxErecon[is] <- max(EkeVreconsPairs[,ie,is,ifl])
            minErecon[is] <- min(EkeVreconsPairs[,ie,is,ifl])
            # look for minimum separation where points enter the orange band
            if((meanErecon[is]<=bandMax && meanErecon[is]>=bandMin ) || 
               (maxErecon[is] <=bandMax && maxErecon[is] >=bandMin ) ||
               (minErecon[is] <=bandMax && minErecon[is] >=bandMin ) ||
               (minErecon[is] <=bandMin && maxErecon[is] >=bandMax )){
                if (separationsToPlot[ie,ifl] == ""){
                    separationsToPlot[ie,ifl] <- separations[is]
                }else{
                    separationsToPlot[ie,ifl] <- paste(separationsToPlot[ie,ifl],",",separations[is],sep="")
                }
                #cat("Adding",separations[is]," ")
            }
        }
        fileBPfail <- paste("baselineLPA2/BPfail_",Esec,"keV",samprateStr,".dat",sep="")
        save(separationsToPlot,file=fileBPfail)
        bgpseps <-unlist(strsplit(separationsToPlot[ie,ifl],","))
        bgpseps <- bgpseps[bgpseps != ""]
        bagplots.separations[[ifl]] <-bgpseps
        
        x1 = 0
        x2 = fLength-1
        if(length(bagplots.separations[[ifl]]) >0){
            x1 <- min(0,as.numeric(bagplots.separations[[ifl]]))
            x2 <- min(fLength-1,max(as.numeric(bagplots.separations[[ifl]])))
        }
        errbar(x=as.numeric(separations),y=meanErecon, yplus=maxErecon, yminus=minErecon,
               col="darkmagenta", pch=1, typ="b",cex=0.7, errbar.col="cornflowerblue",
               xlab="Pair separation (samples)", ylab="<Erecon>", xlim=c(0,xmax[ifl]))
        title(main=(bquote(paste(E[prim],"=",.(energies[ie])," keV  ",
                                 E[sec],"=",.(Esec)," keV  ", "FLength=",.(fLength),
                                 "  FEnergy=",.(fEnergy)," keV",sep=""))),cex.main=0.8)
        abline(h=bandMax, col="orange",lty=2)
        abline(h=bandMin, col="orange",lty=2)
        abline(v=fLength, col="black",lty=2)
        abline(v=x1, col="orange", lty=2)
        abline(v=x2, col="orange", lty=2)
        text(max(0,(x1-20)),(bandMax+0.01), x1,cex=0.7)
        text(min(xmax[ifl],(x2+20)),(bandMax+0.01), x2,cex=0.7)
    }
    #
    # Draw Bagplots  and reconstruction curves
    #
    for (ifl in 1:length(filterLengths)){
         fLength=filterLengths[ifl]
         #minSep <- as.numeric(separationsToPlot[ie,indexEsec,ifl,1])
         #maxSep <- as.numeric(separationsToPlot[ie,indexEsec,ifl,2])
         cat("Plotting figure for Eprim=",energies[ie],"keV, Esec=",Esec, "and Filter Length=",fLength,"\n")
         xmin2<-min(derivArrayMean4samples[,ie,ifl],derivArrayMean4samplesPair[,ie,,ifl])
         xmax2<- max(derivArrayMean4samples[,ie,ifl],derivArrayMean4samplesPair[,ie,,ifl])
    
         # use as ylimits those of the most extreme bagplots among those to be plotted
         ymin2 <- min(EkeVrecons[,ie,ifl],
                     EkeVreconsPairs[,ie,which(separations %in% bagplots.separations[[ifl]]),ifl])
         ymax2 <- max(EkeVrecons[,ie,ifl],
                     EkeVreconsPairs[,ie,which(separations %in% bagplots.separations[[ifl]]),ifl])
    
         # Yellow-ish bagplot (single pulses)
         bagplot(derivArrayMean4samples[,ie,ifl],EkeVrecons[,ie,ifl],xlab="<4 derivative samples>",
                 ylab="Reconstructed Energy (keV)",
                 main=(bquote(paste(E[prim],"=",.(energies[ie])," keV  ",
                                    E[sec],"=",.(Esec)," keV  ", "FLength=",.(fLength),
                                    "  FEnergy=",.(fEnergy)," keV",sep=""))), cex.main=0.8,
                 xlim=c(floor(xmin2),ceiling(xmax2*1.005)), ylim=c(ymin2,ymax2),
                 col.loophull="cornsilk",col.looppoints="peachpuff",col.baghull="orange",show.outlier = FALSE)
                 # blue-ish bagplots (pairs of pulses)
         for(is in 1:length(separations)){
             if(!separations[is] %in% bagplots.separations[[ifl]]) next
             if(as.numeric(separations[is]) >= fLength-1) next
             cat("............Plotting bagplot for sep=",separations[is],"samples\n")
             bag <- try(compute.bagplot(derivArrayMean4samplesPair[,ie,is,ifl],EkeVreconsPairs[,ie,is,ifl]))
    
             if(class(bag) == "try-error") next
             plot.bagplot(bag,add=TRUE, transparency=TRUE)
    
             text(xmax2,ymax2,"Separation(sam)", cex=0.8 )
    
             #text(xmax*1.001,bag$center[2],paste(separations[is],"sam/",seps.ms[is],"ms",sep=""), cex=0.7)
             text(xmax2*1.001,bag$center[2],separations[is], cex=0.7)
         } #separations
         abline(h=max(EkeVrecons[,ie,ifl]), col="orange",lty=2)
         abline(h=min(EkeVrecons[,ie,ifl]), col="orange",lty=2)

     } # filter lengths
}
dev.off()

