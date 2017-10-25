#
# Calculate Frequency of pairs given a spectrum
#
rm(list=ls())
library(FITSio)
library(Hmisc)
source("~/R/Rfunctions/drawLogPlotBox.r")
setwd("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/baselineLPA2")
sepRange <- c(200,300) # separation range for failure in rejectability plots (see bagplots)
samprate <- 156250 # Hz-1 - sampling rate
XT <- "" #crosstalk for filenames
EkeVrecons <- numeric()
evt.times  <- numeric()
evt.phids  <- numeric()
evt.pixids <- numeric()
pix.phids  <- numeric()
pdf(paste("e2e/FrequencyPairs",XT,".pdf",sep=""),width=10, height=7,version="1.4")
calibEnergies <- c("0.2","0.5","1","2","3","4","5","6","7","8") 
EprimSakai <- c("0.2","1","2","4","6","8") 
nen <- length(calibEnergies)
nSec <- nen
nPrim <- length(EprimSakai)
separations <- c(4,10,20,40,44,60,90,100,120,200,250,300,400,500,600,800,1000,1600,2000)
nseps <- length(separations)
detectionModes <- c("AD","A1")
#flux210_mCrab=0.000000000021147 = 2.1147E-11 erg/cm2/s 2-10kev   ---> 1 mCrab 
fluxes.mcrab <- c("0.0001", "0.0005", "0.001", "0.005", "0.01", "0.036", "0.13", "0.46",
                  "0.60", "0.80", "1.67", "6.", "21.54", "77.43", "278.26", "1000.") # mCrab
nFluxes <- length(fluxes.mcrab)
filters <- c("_filter_", "_") # Be Filter & w/o filter
legendTitles<-c("Be Filter","NO Be Filter")
legendPos <- c("topleft","top")
legendPch <- c(4,1)
nfilters <- length(filters)
nmodes <- length(detectionModes)
# Bagplots Failure areas (as a funtion of secondary energy)
BP.failRange <- list("0.2"=c(200,300),"0.5"=c(200,300),
                     "1"=c(240,260),"2"=c(240,260),"3"=c(240,260),
                     "4"=c(235,245),"5"=c(235,245),"6"=c(235,245),
                     "7"=c(235,245),"8"=c(235,245))

percentMiss   <- array(data=NA, dim=c(nFluxes,nfilters,nmodes),dimnames = list(rep("",nFluxes), c("BeFilter", "NoFilter"), c("AD","A1")))
percentBPfail <- array(data=NA, dim=c(nFluxes,nfilters,nmodes),dimnames = list(rep("",nFluxes), c("BeFilter", "NoFilter"), c("AD","A1")))
probDetMatrix <- array(data=NA, dim=c(nSec,nseps,nmodes,nPrim), dimnames = list(calibEnergies, separations, detectionModes, EprimSakai))
# Read probability matrices
for (i1 in 1:length(EprimSakai)){
    for(idet in 1:nmodes){
        matrixFile <- paste("../imageMatrix_",detectionModes[idet],"_",EprimSakai[i1],"keV_old.mat",sep="")
        load(matrixFile)
        for(i2 in 1:nSec){
            for(is in 1:nseps){
                probDetMatrix[i2,is,idet,i1] <- mat_detectedPulses[i2,is]                
            }
        }
    }
}

# Start processing
#==================
for (ifi in 1:nfilters){
    cat("Working with filter ", filters[ifi],"\n")
    filter <- filters[ifi]
    for (ifl in 1:nFluxes){
        if(fluxes.mcrab[ifl] < 0.5) {
            # read evt file
            evtFile <- paste("e2e/crabSpec",fluxes.mcrab[ifl],"mCrab",XT,".fits",sep="")
            pixFile <- paste("e2e/crabSpec",fluxes.mcrab[ifl],"mCrab",XT,".piximpact",sep="")
        }else{
            
            evtFile <- paste("e2e/crabSpec",fluxes.mcrab[ifl],"mCrab",filter,"35mm",XT,".fits",sep="")
            pixFile <- paste("e2e/crabSpec",fluxes.mcrab[ifl],"mCrab",filter,"35mm",XT,".piximpact",sep="")
        }
        
        cat("Working with ", evtFile,"\n")
        zz <- file(description = evtFile, open = "rb")
        header0 <- readFITSheader(zz, fixHdr = 'remove') # read primary header
        header <- readFITSheader(zz, fixHdr = 'remove') # read extension header
        zz.hdr <- parseHdr(header)
        evtTable <- readFITSbintable(zz, header)
        close(zz)
        # read simulated & triggered events
        ntrigs <- as.numeric(zz.hdr[which(zz.hdr=="NIMP")+1])
        nsims <-as.numeric(zz.hdr[which(zz.hdr=="NESTOT")+1])
        cat("   reading table in EVT file\n")
        idcol <- which(evtTable$colNames == "SIGNAL")
        EkeVrecons <- evtTable$col[[idcol]][1:ntrigs]
        idcol <- which(evtTable$colNames == "TIME")
        evt.times <- evtTable$col[[idcol]][1:ntrigs]
        idcol <- which(evtTable$colNames == "PH_ID")
        evt.phids <- evtTable$col[[idcol]][1:ntrigs]
        idcol <- which(evtTable$colNames == "PIXID")
        evt.pixids <- evtTable$col[[idcol]][1:ntrigs]
        
        # Effectivity of detections
        if(any(EkeVrecons == 0.0)){
            cat("   reading piximpact file to get sim energy of '0.000' photons\n")
            pp <- file(description = pixFile, open = "rb")
            header0 <- readFITSheader(pp, fixHdr = 'remove') # read primary header
            header <- readFITSheader(pp, fixHdr = 'remove') # read extension header
            pp.hdr <- parseHdr(header)
            evtTable <- readFITSbintable(pp, header)
            close(pp)
            idcol <- which(evtTable$colNames == "ENERGY")
            EkeVsim <- evtTable$col[[idcol]][1:nsims]
            idcol <- which(evtTable$colNames == "PH_ID")
            pix.phids <- evtTable$col[[idcol]][1:nsims]
        }
        for(idet in 1:nmodes){
            detMod <- detectionModes[idet]
            cat(" Detection mode: ",detectionModes[idet],"\n")    
            cat("   Photons loop (nsims=",nsims,"ntrigs=",ntrigs,")\n")
            probMiss   <- rep(0,ntrigs)
            probBPfail <- rep(0,ntrigs)
            for (ip in 1:ntrigs){
                #if(ip == 1 || ip %% 1000 == 0) cat("   Working with photon: ", ip)
                # save simulated energy instead of 0.00 energy
                if(EkeVrecons[ip] == 0.0) EkeVrecons[ip] <- EkeVsim[pix.phids == evt.phids[ip]]
                if(ip==1) next
                BPfails <- FALSE
                BPfailsConserv <- FALSE
                diffSamples <- (evt.times[ip]-evt.times[ip-1])*samprate
                
                # look for the closest energy of photons Prim & Sec in energies list and closest separation in matrices
                idclosPrim <- whichClosest(EprimSakai,EkeVrecons[ip-1])
                #closestEnergyPrim <- calibEnergies[idclosPrim]
                closestEnergyPrim <- EprimSakai[idclosPrim]
                idclosSec <- whichClosest(calibEnergies,EkeVrecons[ip])
                closestEnergySec <- calibEnergies[idclosSec]
                idclosSep <- whichClosest(separations,diffSamples)
                closestSep <- separations[idclosSep]
        
                # BPA1AD: Pulses not rejected by bagplots and not detected by A1/AD
                #===============================================================
                # if in same pixel than previous && in required temporal range &&
                #    energy prim/sec in range of interest && in same pixel:
                #  add current & previous pulse to "affected" to take into account 
                #  also primaries in the count (later discard repetitions)
                if((diffSamples>=BP.failRange[[closestEnergySec]][1] && diffSamples<=BP.failRange[[closestEnergySec]][2]))
                    BPfails <- TRUE
                if((diffSamples>=BP.failRange[[calibEnergies[1]]][1] && diffSamples<=BP.failRange[[calibEnergies[1]]][2]))
                    BPfailsConserv <- TRUE
            
                # ip = secondary
                # ip -1 = primary
                #======================================================================
                # prob=140%     
                #     Prim (50pulses)      Secondaries (50 pulses)
                #         (50)ok                 (10) ok -> 10/50=20%  (well)         
                #                            (40+40 very close) -> 40/50=80% (bad)
                # prob=70%     
                #     Prim (50pulses)      Secondaries (50 pulses)
                #         (50)ok             (20) ok -> 20/50=40%  (well)         
                #                            (30) missing-> 30/50=60% (bad)
                # prob=30%     
                #     Prim (50pulses)      Secondaries (50 pulses)
                #         (30)ok             (0) ok -> 0/50=0%  (well)         
                #         ->30/50=60%        (50) missing-> 50/50=100% (bad)    
                
                idxPrim <- paste(detMod,ip-1,sep="_")
                idxSec  <- paste(detMod,ip,sep="_")
                probDet <- probDetMatrix[idclosSec,idclosSep,idet,idclosPrim] # 0-100%
                # Probs in 0-1 scale
                if(probDet>=50){
                    probDetPrim <- 1.
                    probMissPrim <- 0.                       # primary
                    probMissSec <- abs(probDet-100)/50.   # secondary
                    probDetSec <- 1-probMissSec 
                }else if (probDet < 50){
                    probDetPrim <- probDet/50.
                    probMissPrim <- 1.-probDetPrim      # primary
                    probDetSec <- 0.
                    probMissSec <- 1. #secondary
                }
                #probDetection[[idxSec]] <- append(probDetection[[idxSec]], probDetSec)
                #probDetection[[idxPrim]]<- append(probDetection[[idxPrim]], probDetPrim)
                probMiss[ip] <- probMissSec
                if(BPfails) probBPfail[ip] <- probMissSec # if not in 'bad area' probBPfail=0
            } #photons (given a ctrate and detection mode)
            
            #for (pulseid in probDetection$names){
            #    totalProbDetection[[pulseid]] <- prod(probDetection[[pulseid]])
            #}
            percentMiss[ifl,ifi,idet] <- sum(probMiss)/nsims # 0-1
            percentBPfail[ifl,ifi,idet] <- sum(probBPfail)/nsims # 0-1
        } #detMod
    } #each flux
}#each filter
par(mfrow=c(1,2))

# PLOTTING (UN)DETECTION
##########################
drawLogPlotBox(xlimits=c(1E-4,1E3),ylimits=c(0,100), logxy="x",xlabel="Intensity (mCrab)
               1 mCrab=2.11E-11 erg/cm2/s 2-10kev", 
               ylabel=expression("Fraction of photons (%)"), naxes=c(T,T,F,F))
title(main="Undetected Fraction of Photons")
abline(v=0.5,col="grey",lty=2)
for (ifi in 1:nfilters){
    for(idet in 1:nmodes){
        coldet<- "blue"
        if(detectionModes[idet] == "A1") coldet <- "red"
        points(fluxes.mcrab,percentMiss[,ifi,idet]*100,col=coldet,type = "b",pch=legendPch[ifi],cex=0.8) 
    }
    #legend(legendPos[ifi],legend=c("A1","AD"), col=c("red","blue"),
    #       pch=c(legendPch[ifi],legendPch[ifi]),cex=0.7,bty="n", title=legendTitles[ifi])
}
legend("topleft",legend=c(paste(detectionModes[1],legendTitles[1]),
                          paste(detectionModes[2],legendTitles[1]),
                          paste(detectionModes[1],legendTitles[2]),
                          paste(detectionModes[2],legendTitles[2])),        
                col=c("blue","red","blue","red"),
                pch=c(legendPch[1],legendPch[1],legendPch[2],legendPch[2]),
                cex=0.7,bty="n")

text(1E-3,80,"PSF=athena_psf_onaxis_20150602.fits",cex=0.4)
text(1E2,80,"PSF=\nathena_ladapt_defocus_35mm_kev_psf_20161013.fits",cex=0.4)
# PLOTTING (UN)DETECTION * BAGPLTS flagging
##############################################
drawLogPlotBox(xlimits=c(1E-4,1E3),ylimits=c(0,5E-3), logxy="x",xlabel="Intensity (mCrab)", 
               ylabel=expression("Fraction of photons (%)"), naxes=c(T,T,F,F))
title(main="Undetected Fraction of Photons \n(not spotted by derivative)",cex=0.8)
for (ifi in 1:nfilters){
    for(idet in 1:nmodes){
        coldet<- "blue"
        if(detectionModes[idet] == "A1") coldet <- "red"
        points(fluxes.mcrab,percentBPfail[,ifi,idet]*100,col=coldet,lty=1,type = "b",pch=legendPch[ifi],cex=0.8)
    }
    #legend(legendPos[ifi],legend=c("A1","AD"), 
    #       col=c("red","blue"), pch=c(legendPch[ifi],legendPch[ifi]),cex=0.7,bty="n",
    #       title=legendTitles[ifi])
}
legend("topleft",legend=c(paste(detectionModes[1],legendTitles[1]),
                          paste(detectionModes[2],legendTitles[1]),
                          paste(detectionModes[1],legendTitles[2]),
                          paste(detectionModes[2],legendTitles[2])),        
                          col=c("blue","red","blue","red"),
                          pch=c(legendPch[1],legendPch[1],legendPch[2],legendPch[2]),
                          cex=0.7,bty="n")
abline(v=0.5,col="grey",lty=2)
text(1E-3,4E-3,"PSF=athena_psf_onaxis_20150602.fits",cex=0.4)
text(1E2,4E-3,"PSF=\nathena_ladapt_defocus_35mm_kev_psf_20161013.fits",cex=0.4)

dev.off()
