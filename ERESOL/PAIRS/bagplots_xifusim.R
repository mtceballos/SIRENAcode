#
# BAGPLOTS: plots the bagplots in 3 columns (3 filters) for all the combinations Eprim/Esec
# 
# Simulated files (xifusim) must be previously run (bagplots.csh):
#      - singles
#      - pairs as if secondary pulse has been missed (in noDetSP folder)
#
#
#rm(list=ls())
read=0
source("~/R/Rfunctions/drawLogPlotBox.r")
library(scales)
dcmt <- 100
domain="T" #or 'F0F'
thresholdProb=0.01 # 1% - to be used to decide if separation is in singles area
args = commandArgs(trailingOnly = TRUE)
# test if there is one argument: if not, return an error
if (length(args)==0) {
    # default secondary energy
    args[1] = "0.2"
}
#cat("args=",args[1])

library(FITSio)
library(aplpack)
library(Hmisc)
setwd("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/")
# --------------------------------------------------------------------------------------
#  For input: 
#
#npulses <- 1000      # Number of pulses (singles or pairs) at each energy
#nPairs <- 1000
nSinglePulses <- 5000
nPairs <- 200
npulses <- nPairs*2
samprate <- 156250
samprateStr="" # to name files with full samprate
#samprateStr="_samprate2" # to name files with 1/2 samprate

filterLengths <- c(8192,512,256)
pulseLength<- 8192   # pulse length
xmax <- c(8200, 550, 260 )
separationsBGplots <- sprintf("%05d", sort(c(seq(30,250,40),235,seq(345,470,10),
                                             seq(245,255,5),seq(260,510,40),
                                             seq(510,515,5),seq(520,8130,100), 
                                             seq(7900,8151,50),seq(8190,8200,5))))

if (samprateStr == "_samprate2"){
    samprate <- samprate/2.    
    filterLengths <- filterLengths/2
    pulseLength<- pulseLength/2   # pulse length
    xmax <- c(4200, 260, 150 )
    #separationsBGplots <- sprintf("%05d",
    #                       sort(c(15,seq(19,25,2),seq(16,34,2),35,seq(40,256,2),300,
    #                        seq(320,340,10),seq(360,400,20),450,500,750,800,seq(1000,4000,250),
    #                        seq(4010,4196,10),4200)))
    separationsBGplots <- sprintf("%05d", round(as.numeric(separationsBGplots)/2))
}

fEnergy="6" # filter energy in keV
#energies <- c("0.2", "0.5", "1", "2", "3", "4", "5", "6", "7", "8") # primary pulses energies
Eprims <- c("0.2", "0.5", "1", "4", "6", "8") # primary pulses energies
nprims <- length(Eprims)
Esecs <- c("0.2", "0.5", "0.8", "1", "1.3", "2", "2.3", "3", "3.3",
           "4", "4.3", "5", "5.3", "6", "6.3", "7", "7.3", "8")# keV : energy of secondary pulses
nsecs <- length(Esecs)
nseps <- length(separationsBGplots)
nfilters <- length(filterLengths)
seps.ms <- as.numeric(separationsBGplots)/samprate*1E3 #separations in ms
minsep <- min(as.numeric(separationsBGplots))
maxsep <- max(as.numeric(separationsBGplots))
# --------------------------------------------------------------------------------------

if(read){
    arrayLowResEstimE <- array(data=NA,dim=c(nSinglePulses,nprims,nfilters))
    arrayLowResEstimEpair <- array(data=NA, dim=c(nPairs,nprims,nsecs,nseps,nfilters))
    arrayPhiPair <- array(data=NA, dim=c(nPairs,nprims,nsecs,nseps,nfilters))

#=====================================================================
#                                                                    #
#                              SINGLE PULSES                         #
#                                                                    #
#=====================================================================


# Get reconstructed energies AND Low-Res Energy estimation of SINGLE pulses
#============================================================================
EkeVrecons <- array(data=NA,dim=c(nSinglePulses,nprims,nfilters))

for (ip in 1:nprims){
    cat("Working with SINGLES energy=",Eprims[ip],"\n")
    for (ifl in 1:nfilters){
        fl <- filterLengths[ifl]
        # get reconstructed energies
        eventsFile <- paste("eresolLPA75um/nodetSP/events_sep40000sam_5000p_SIRENA",
                            pulseLength,"_pL",fl,"_",Eprims[ip],
                            "keV_STC_",domain,"_fixedlib",fEnergy,"OF_OPTFILT",pulseLength,
                            samprateStr,"_jitter_dcmt",dcmt,"_HR.fits",sep="")
        zz <- file(description = eventsFile, open = "rb")
        header0 <- readFITSheader(zz, fixHdr = 'remove') # read primary header
        header <- readFITSheader(zz, fixHdr = 'remove') # read extension header
        evtTable <- readFITSbintable(zz, header)
        close(zz)
        idcol <- which(evtTable$colNames == "SIGNAL")
        EkeVrecons[,ip,ifl] <- evtTable$col[[idcol]][1:nSinglePulses]
        # get low resolution (4-samples length filter) E estimation
        idcol <- which(evtTable$colNames == "ELOWRES")
        arrayLowResEstimE[,ip,ifl] <-  evtTable$col[[idcol]][1:nSinglePulses]
    } # foreach filter length
} # foreach calib energy
cat("SINGLE pulses (yellow bagplot) finished\n")
#=====================================================================
#                                                                    #
#                              PAIR PULSES                           #
#                                                                    #
#=====================================================================

# Get reconstructed energies AND LowResEstimation of PAIR pulses
#=======================================================================
EkeVreconsPairs <- array(NA, dim=c(nPairs, nprims, nsecs, nseps, nfilters))
# Get reconstructed energies of PAIRS of pulses
for (ip in 1:nprims){
    Eprim <- Eprims[ip]
    for (is in 1:nsecs){
        Esec <- Esecs[is]
        cat("Working with PAIRS: EnergyPrim=",Eprim,"EnergySec=",Esec,"\n")
        for (iss in 1:nseps){
            for (ifl in 1:length(filterLengths)){
                fl <- filterLengths[ifl]
                if(separationsBGplots[iss] == "00002"){ # too short filter: cannot calculate energy
                    EkeVreconsPairs[,ip,is,iss,ifl] <- rep(NaN,nPairs)
                    next
                }
                # get reconstructed energies
                eventsFile <- paste("eresolLPA75um/nodetSP/events_sep",separationsBGplots[iss],
                                "sam_",npulses,"p_SIRENA",pulseLength,"_pL", fl,
                                "_",Eprim,"keV_",Esec,"keV_STC_",domain,"_fixedlib",fEnergy,
                                "OF_OPTFILT",pulseLength,samprateStr,"_jitter_dcmt",
                                dcmt,".fits",sep="")
                #cat("  Reading file:", eventsFile,"\n")
                zz <- file(description = eventsFile, open = "rb")
                header0 <- readFITSheader(zz, fixHdr = 'none') # read primary header
                header <- readFITSheader(zz, fixHdr = 'none') # read extension header
                evtTable <- readFITSbintable(zz, header)
                close(zz)
                idcol <- which(evtTable$colNames == "SIGNAL")
                EkeVreconsPairs[,ip,is,iss,ifl] <- evtTable$col[[idcol]][1:nPairs]
                # get low resolution (4-samples length filter) E estimation
                idcol <- which(evtTable$colNames == "ELOWRES")
                arrayLowResEstimEpair[,ip,is,iss,ifl] <- evtTable$col[[idcol]][1:nPairs]
                idcol <- which(evtTable$colNames == "PHI")
                arrayPhiPair[,ip,is,iss,ifl] <- evtTable$col[[idcol]][1:nPairs]
            } # foreach filter length    
        } # foreach separation    
    } # foreach secondary energy
} # foreach primary energy
cat("PAIR pulses (bluish bagplots) finished\n")
} #finish reading

#
# PLOTTING
#
#pdf(paste("baselineLPA75um/bagplotsStudy",samprateStr,".pdf",sep=""),width=10, height=7,version="1.4")
pdf(paste("baselineLPA75um/bagplotsStudy",samprateStr,".pdf",sep=""),width=7, height=10,version="1.4")

separationsInSingles <- array(NA, dim=c(nseps, nprims, nsecs, nfilters))
for (ip in 1:nprims){
    Eprim <- Eprims[ip]
    cat("plotting for Eprim=", Eprim,"\n")
    #
    # Plot variation of Reconstructed energy with separation for every filter length
    # and secondary energy
    #layout(matrix(c(1,3,5,2,4,6),2,3,byrow = TRUE), widths=c(1,1,1), heights = c(2,1))
    layout(matrix(c(1,5,9,2,6,10,3,7,11,4,8,12),4,3,byrow = TRUE), 
           widths=c(1,1,1,1), heights = c(1,1,1))
    for (is in 1:nsecs){
        Esec <- Esecs[is]
        for (ifl in 1:length(filterLengths)){
            fLength=filterLengths[ifl]
            # define orange band
            bandMin <- min(EkeVrecons[,ip,ifl])
            bandMax <- max(EkeVrecons[,ip,ifl])
        
            #
            # Draw reconstruction curves
            #
            minE<-max(1E-5,min(EkeVreconsPairs[,ip,is,,ifl],na.rm=TRUE))
            maxE<-max(EkeVreconsPairs[,ip,is,,ifl],na.rm=TRUE)
            
            plot(seq(1,10),seq(1,10),log="y",type="n",
                 xlim=c(minsep,min(fLength+0.3*fLength,maxsep)), ylim=c(min(minE,bandMin), maxE),
                 xlab="Pair separation (samples)", ylab="Reconstructed Pulse Height")
            inprob <- as.numeric()
            for (iss in 1:length(separationsBGplots)){
                # not NA reconstructed energies
                ypoints <- EkeVreconsPairs[!is.na(EkeVreconsPairs[,ip,is,iss,ifl]),ip,is,iss,ifl]
                # points inside orange interval
                inpoints <- ypoints[(ypoints<=bandMax & ypoints>=bandMin)]
                inprob <- append(inprob,length(inpoints)/length(ypoints)) #prob inside interval
            
                #draw one point for each pulse
                points(x=rep(as.numeric(separationsBGplots[iss]),length(ypoints)),
                    y=ypoints, pch=19, cex=0.2,col="darkmagenta")
                #cat("inpoints=",length(inpoints),"ypoints=",length(ypoints),"\n")
                # label position of mean of points
                text(as.numeric(separationsBGplots[iss]),mean(ypoints),labels="--",
                        col="blue")
                #cat("mean=",mean(EkeVreconsPairs[,ie,is,ifl]))
            }
            title(main=(bquote(paste(E[prim],"=",.(Eprim)," keV  ",
                                 E[sec],"=",.(Esec)," keV  ", "FLength=",.(fLength),
                                 "  FEnergy=",.(fEnergy)," keV",sep=""))),cex.main=0.8)
            abline(h=bandMax, col="orange",lty=2)
            abline(h=bandMin, col="orange",lty=2)
            abline(v=fLength, col="black",lty=2)
        
            #
            # Draw probability curves
            #
            plot(as.numeric(separationsBGplots),inprob,pch=1, cex=0.5,
                 xlim=c(minsep,min(fLength+100,maxsep)),lty=1,type = 'b',
                xlab="Pair separation (samples)", ylab="Probability of being in singles region")
            abline(v=fLength, col="black",lty=2)
            abline(h=thresholdProb, col="orange",lty=1)
            
            # separations where the probability of having double pulses inside singles zone
            # is larger than thresholdProb
            badseps <-separationsBGplots[which(inprob>thresholdProb&
                                                        as.numeric(separationsBGplots)<=fLength)]
            if(length(badseps)>0){
                separationsInSingles[1:length(badseps),ip,is,ifl]<-badseps
            }
            points(separationsInSingles, rep(0.1,length(separationsInSingles)),
                   pch=4,cex=0.5,col="orange")
            #
            # Draw reconstruction vs low-resolution estimation of energy
            #
            minlow <- min(arrayLowResEstimE[,ip,ifl], na.rm = TRUE)
            maxlow <- max(arrayLowResEstimE[,ip,ifl], na.rm = TRUE)
            ypos1 <- numeric()
            ypos2 <- numeric()
            iss1 <- numeric()
            iss2 <- numeric()
            iss1 <-append(iss1,1)
            iss2 <-append(iss2,2)
            ypos1 <- append(ypos1,mean(EkeVreconsPairs[,ip,is,1,ifl],na.rm=TRUE))
            ypos2 <- append(ypos2,mean(EkeVreconsPairs[,ip,is,2,ifl],na.rm=TRUE))
            
            plot(seq(1,10),seq(1,10),log="y",type="n",
                 xlim=c(minlow,maxlow), ylim=c(minE, maxE),
                 xlab="Low-resolution PH estimation", ylab="Reconstructed PH")
            for (iss in 1:length(separationsBGplots)){
                points(arrayLowResEstimEpair[,ip,is,iss,ifl], EkeVreconsPairs[,ip,is,iss,ifl],
                       pch=1, cex=0.3, col = alpha("blue", 0.4))
                ypos <- mean(EkeVreconsPairs[,ip,is,iss,ifl],na.rm=TRUE)
                if(iss<3) next
                if(! separationsBGplots[iss] %in% separationsInSingles) next
                if (iss %% 2 == 0 ){
                    #&&  abs(ypos - ypos2[whichClosest(ypos2, ypos)])>(maxE-minE)*0.03){
                    ypos2 <- append(ypos2,ypos)
                    iss2 <- append(iss2,iss)
                    #cat("Adding label", separationsBGplots[iss],"\n")
                }else if(iss %% 2 == 1 ){
                         #&& abs(ypos - ypos1[whichClosest(ypos1, ypos)])>(maxE-minE)*0.009){
                    ypos1 <- append(ypos1,ypos)
                    iss1 <- append(iss1,iss)
                    #cat("Adding label", separationsBGplots[iss],"\n")
                }
            }
            points(arrayLowResEstimE[,ip,ifl], EkeVrecons[,ip,ifl],
                   pch=1,cex=0.3, col = alpha("orange", 0.4))
            # print separation labels so that they do not overplot
            text(x=runif(length(ypos1),max=maxlow, min=minlow),
                 y=ypos1,labels=paste(separationsBGplots[iss1]),cex=0.5)
            text(x=runif(length(ypos2),max=maxlow, min=minlow),
                 y=ypos2,labels=paste(separationsBGplots[iss2]),cex=0.5)
            #
            # print reconstruction vs offset
            #
            minPHI <- min(arrayPhiPair[,ip,is,,ifl], na.rm = TRUE)
            maxPHI <- max(arrayPhiPair[,ip,is,,ifl], na.rm = TRUE)
            plot(seq(1,10),seq(1,10),log="y",type="n",
                 xlim=c(minPHI,maxPHI), ylim=c(minE, maxE),
                 xlab="Arrival Offset (samples)", ylab="Reconstructed PH")
            for (iss in 1:length(separationsBGplots)){
                points(arrayPhiPair[,ip,is,iss,ifl], EkeVreconsPairs[,ip,is,iss,ifl],
                       pch=1, cex=0.3, col = alpha("blue", 0.4))
            }
            
        } # foreach ifl (curves)
        
    } # foreach secondary energy
} # foreach primary energy

# save interesting quantities for detectionMaps y frequencyPairs:
fileBPfail <- paste("baselineLPA75um/BPfail", samprateStr,".dat",sep="")
save(separationsInSingles,file=fileBPfail)
dev.off()

