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
nen <- length(calibEnergies)
#srcFluxes <- c("1E-11", "2.15E-11", "4.64E-11", "1E-10", "2.15E-10", "4.64E-10", 
#                "1E-9", "2.15E-9", "4.64E-9", "1E-8", "2.15E-8")
#ctr.mcrab<-numeric(length(srcFluxes))
fluxes.mcrab <- c("0.0001","0.0005","0.001","0.005","0.01", "0.036", "0.13", "0.46", 
                  "0.80", "1.67", "6.", "21.54", "77.43", "278.26", "1000.") # flux in mCrab
nFluxes <- length(fluxes.mcrab)
filters <- c("_filter_", "_") # Be Filter & w/o filter
legendTitles<-c("Be Filter","NO Be Filter")
legendPos <- c("topright","right")
legendPch <- c(4,1)
nfilters <- length(filters)

# Bagplots Failure areas (as a funtion of secondary energy)
BP.failRange <- list("0.2"=c(200,300),"0.5"=c(200,300),
                     "1"=c(240,260),"2"=c(240,260),"3"=c(240,260),
                     "4"=c(235,245),"5"=c(235,245),"6"=c(235,245),
                     "7"=c(235,245),"8"=c(235,245))
# limiting values for detection for different primary energies (for 0.2keV secondaries) 
minSampAD_lt_02 <- list("0.2"=45, "0.5"=5, "1"=5, "2"=5, "3"=5, "4"=5, "5"=5, "6"=5, "7"=5, "8"=5) # using lower threshold
minSampAD_ht_02 <- list("0.2"=45, "0.5"=5, "1"=5, "2"=5, "3"=100, "4"=200, "5"=200, "6"=200, "7"=200, "8"=200) # using lower threshold
minSampA1_lt_02 <- list("0.2"=20, "0.5"=45, "1"=200, "2"=200, "3"=250, "4"=250, "5"=300, "6"=300, "7"=400, "8"=400)
minSampA1_ht_02 <- list("0.2"=200, "0.5"=250, "1"=300, "2"=400, "3"=400, "4"=400, "5"=800, "6"=800, "7"=800, "8"=800)

# limiting values for detection for different primary energies (for 1keV secondaries) - better detection
minSampAD_lt_1 <- list("0.2"=60, "0.5"=20, "1"=5, "2"=5, "3"=5, "4"=5, "5"=5, "6"=5, "7"=5, "8"=5) # using lower threshold
minSampAD_ht_1 <- list("0.2"=60, "0.5"=20, "1"=5, "2"=5, "3"=5, "4"=5, "5"=5, "6"=5, "7"=5, "8"=5) # using lower threshold
minSampA1_lt_1 <- list("0.2"=45, "0.5"=45, "1"=45, "2"=45, "3"=45, "4"=45, "5"=45, "6"=45, "7"=45, "8"=45)
minSampA1_ht_1 <- list("0.2"=5, "0.5"=10, "1"=20, "2"=20, "3"=20, "4"=20, "5"=20, "6"=20, "7"=20, "8"=200)

percent.AD_lt <- matrix(nrow=nFluxes,ncol = nfilters, dimnames=list(rep("",nFluxes), c("BeFilter", "NoFilter")))
percent.AD_ht <- matrix(nrow=nFluxes,ncol = nfilters, dimnames=list(rep("",nFluxes), c("BeFilter", "NoFilter")))
percent.A1_lt <- matrix(nrow=nFluxes,ncol = nfilters, dimnames=list(rep("",nFluxes), c("BeFilter", "NoFilter")))
percent.A1_ht <- matrix(nrow=nFluxes,ncol = nfilters, dimnames=list(rep("",nFluxes), c("BeFilter", "NoFilter")))
percent.BPAD_lt <- matrix(nrow=nFluxes,ncol = nfilters, dimnames=list(rep("",nFluxes), c("BeFilter", "NoFilter")))
percent.BPAD_ht <- matrix(nrow=nFluxes,ncol = nfilters, dimnames=list(rep("",nFluxes), c("BeFilter", "NoFilter")))
percent.BPA1_lt <- matrix(nrow=nFluxes,ncol = nfilters, dimnames=list(rep("",nFluxes), c("BeFilter", "NoFilter")))
percent.BPA1_ht <- matrix(nrow=nFluxes,ncol = nfilters, dimnames=list(rep("",nFluxes), c("BeFilter", "NoFilter")))
percent.BPA1_ht_cnsrv <- matrix(nrow=nFluxes,ncol = nfilters, dimnames=list(rep("",nFluxes), c("BeFilter", "NoFilter")))

for (ifi in 1:nfilters){
    filter <- filters[ifi]
    for (i in 1:nFluxes){
        if(fluxes.mcrab[i] < 0.5) {
            # read evt file
            evtFile <- paste("e2e/crabSpec",fluxes.mcrab[i],"mCrab",XT,".fits",sep="")
            pixFile <- paste("e2e/crabSpec",fluxes.mcrab[i],"mCrab",XT,".piximpact",sep="")
        }else{
            
            evtFile <- paste("e2e/crabSpec",fluxes.mcrab[i],"mCrab",filter,"35mm",XT,".fits",sep="")
            pixFile <- paste("e2e/crabSpec",fluxes.mcrab[i],"mCrab",filter,"35mm",XT,".piximpact",sep="")
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
        missedAD_lt_ids <- numeric()
        missedAD_ht_ids <- numeric()
        missedA1_lt_ids <- numeric()
        missedA1_ht_ids <- numeric()
        BPAD_lt_ids <- numeric()
        BPAD_ht_ids <- numeric()
        BPA1_lt_ids <- numeric()
        BPA1_ht_ids <- numeric()
        BPA1_ht_cnsrv_ids <- numeric()
        cat("   Photons loop (nsims=",nsims,"ntrigs=",ntrigs,")\n")
        for (ip in 1:ntrigs){
            # save simulated energy instead of 0.00 energy
            if(EkeVrecons[ip] == 0.0) EkeVrecons[ip] <- EkeVsim[pix.phids == evt.phids[ip]]
            if(ip==1) next
            BPfails <- FALSE
            BPfailsConserv <- FALSE
            diffSamples <- (evt.times[ip]-evt.times[ip-1])*samprate
            
            # look for the closest energy of photons Prim & Sec in calib list
            idclosPrim <- whichClosest(calibEnergies,EkeVrecons[ip-1])
            closestEnergyPrim <- calibEnergies[idclosPrim]
            idclosSec <- whichClosest(calibEnergies,EkeVrecons[ip])
            closestEnergySec <- calibEnergies[idclosSec]
        
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

            # missedAD_lw/ht: Pulses missed by AD using lower and higher threshold 
            # (for secondaries in (0.2-1)keV and (1-8)keV
            #======================================================================
            if(EkeVrecons[ip]<1.){ #secondary
                if(diffSamples<minSampAD_lt_02[[closestEnergyPrim]] && (evt.pixids[ip]==evt.pixids[ip-1])){
                    missedAD_lt_ids <- append(missedAD_lt_ids, c(evt.phids[ip],evt.phids[ip-1]))
                    if(BPfails)  BPAD_lt_ids <- append(BPAD_lt_ids, c(evt.phids[ip],evt.phids[ip-1]))
                }
                if(diffSamples<minSampAD_ht_02[[closestEnergyPrim]] && (evt.pixids[ip]==evt.pixids[ip-1])){
                    missedAD_ht_ids <- append(missedAD_ht_ids, c(evt.phids[ip],evt.phids[ip-1]))
                    #cat("Missed AD_ht for ip=",ip,"-- Diff=",diffSamples,"\n")
                    if(BPfails)  BPAD_ht_ids <- append(BPAD_ht_ids, c(evt.phids[ip],evt.phids[ip-1]))
                }
            }else{
                if(diffSamples<minSampAD_lt_1[[closestEnergyPrim]] && (evt.pixids[ip]==evt.pixids[ip-1])){
                    missedAD_lt_ids <- append(missedAD_lt_ids, c(evt.phids[ip],evt.phids[ip-1]))
                    if(BPfails)  BPAD_lt_ids <- append(BPAD_lt_ids, c(evt.phids[ip],evt.phids[ip-1]))
                }
                if(diffSamples<minSampAD_ht_1[[closestEnergyPrim]] && (evt.pixids[ip]==evt.pixids[ip-1])){
                    missedAD_ht_ids <- append(missedAD_ht_ids, c(evt.phids[ip],evt.phids[ip-1]))
                    #cat("Missed AD_ht for ip=",ip,"-- Diff=",diffSamples,"\n")
                    if(BPfails)  BPAD_ht_ids <- append(BPAD_ht_ids, c(evt.phids[ip],evt.phids[ip-1]))
                }
            }
        
        
            # missedA1_lw/th: Pulses missed by A1 using lower and higher threshold
            #======================================================================
            if(EkeVrecons[ip]<1.){ #secondary
                if(diffSamples<minSampA1_lt_02[[closestEnergyPrim]] && (evt.pixids[ip]==evt.pixids[ip-1])){
                    missedA1_lt_ids <- append(missedA1_lt_ids, c(evt.phids[ip],evt.phids[ip-1]))
                    if(BPfails)  BPA1_lt_ids <- append(BPA1_lt_ids, c(evt.phids[ip],evt.phids[ip-1]))
                }
                if(diffSamples<minSampA1_ht_02[[closestEnergyPrim]] && (evt.pixids[ip]==evt.pixids[ip-1])){
                    missedA1_ht_ids <- append(missedA1_ht_ids, c(evt.phids[ip],evt.phids[ip-1]))    
                    if(BPfails)  BPA1_ht_ids <- append(BPA1_ht_ids, c(evt.phids[ip],evt.phids[ip-1]))
                    if(BPfailsConserv)  BPA1_ht_cnsrv_ids <- append(BPA1_ht_cnsrv_ids, c(evt.phids[ip],evt.phids[ip-1]))
                }
            }else{
                if(diffSamples<minSampA1_lt_1[[closestEnergyPrim]] && (evt.pixids[ip]==evt.pixids[ip-1])){
                    missedA1_lt_ids <- append(missedA1_lt_ids, c(evt.phids[ip],evt.phids[ip-1]))
                    if(BPfails)  BPA1_lt_ids <- append(BPA1_lt_ids, c(evt.phids[ip],evt.phids[ip-1]))
                }
                if(diffSamples<minSampA1_ht_1[[closestEnergyPrim]] && (evt.pixids[ip]==evt.pixids[ip-1])){
                    missedA1_ht_ids <- append(missedA1_ht_ids, c(evt.phids[ip],evt.phids[ip-1]))    
                    if(BPfails)  BPA1_ht_ids <- append(BPA1_ht_ids, c(evt.phids[ip],evt.phids[ip-1]))
                    if(BPfailsConserv)  BPA1_ht_cnsrv_ids <- append(BPA1_ht_cnsrv_ids, c(evt.phids[ip],evt.phids[ip-1]))
                }
            }
        } # trig photons
        cat("   percentages calculation\n")
        # percentage of photons missed by AD/A1 and out-of-scope of BP
        percent.BPA1_lt[i,ifi] <- 100*length(unique(BPA1_lt_ids))/nsims # % fraction of affected photons
        percent.BPA1_ht[i,ifi] <- 100*length(unique(BPA1_ht_ids))/nsims # % fraction of affected photons
        percent.BPA1_ht_cnsrv[i,ifi] <- 100*length(unique(BPA1_ht_cnsrv_ids))/nsims # % fraction of affected photons
        percent.BPAD_lt[i,ifi] <- 100*length(unique(BPAD_lt_ids))/nsims # % fraction of affected photons
        percent.BPAD_ht[i,ifi] <- 100*length(unique(BPAD_ht_ids))/nsims # % fraction of affected photons
        # percentage of photons missed by AD/A1 
        percent.A1_lt[i,ifi] <- 100*length(unique(missedA1_lt_ids))/nsims # % fraction of affected photons
        percent.A1_ht[i,ifi] <- 100*length(unique(missedA1_ht_ids))/nsims # % fraction of affected photons
        percent.AD_lt[i,ifi] <- 100*length(unique(missedAD_lt_ids))/nsims # % fraction of affected photons
        percent.AD_ht[i,ifi] <- 100*length(unique(missedAD_ht_ids))/nsims # % fraction of affected photons
    } #each flux
}#each filter
par(mfrow=c(1,2))

# PLOTTING (UN)DETECTION
##########################
drawLogPlotBox(xlimits=c(1E-4,1E3),ylimits=c(0,2.5), logxy="x",xlabel="Intensity (mCrab)", 
               ylabel=expression("Fraction of photons (%)"), naxes=c(T,T,F,F))
title(main="Undetected Fraction of Photons")
for (ifi in 1:nfilters){
    points(fluxes.mcrab,percent.A1_lt[,ifi],col="red",type = "b",lty=2,pch=legendPch[ifi],cex=0.8) 
    points(fluxes.mcrab,percent.A1_ht[,ifi],col="red",type = "b",lty=1,pch=legendPch[ifi],cex=0.8) 
    points(fluxes.mcrab,percent.AD_lt[,ifi],col="blue",type = "b",lty=2,pch=legendPch[ifi],cex=0.8)
    points(fluxes.mcrab,percent.AD_ht[,ifi],col="blue",type = "b",lty=1,pch=legendPch[ifi],cex=0.8) 

    legend(legendPos[ifi],legend=c("A1, high threshold","A1, low threshold","AD, high threshold","AD, low threshold"), 
           col=c("red","red","blue","blue"), lty=c(1,2,1,2), 
           pch=c(legendPch[ifi],legendPch[ifi],legendPch[ifi],legendPch[ifi]),cex=0.7,bty="n", title=legendTitles[ifi])
    
}
# PLOTTING (UN)DETECTION * BAGPLTS flagging
##############################################
drawLogPlotBox(xlimits=c(1E-4,1E3),ylimits=c(0,0.5), logxy="x",xlabel="Intensity (mCrab)", 
               ylabel=expression("Fraction of photons (%)"), naxes=c(T,T,F,F))
title(main="Undetected Fraction of Photons \n(not spotted by derivative)",cex=0.8)
for (ifi in 1:nfilters){
    points(fluxes.mcrab,percent.BPA1_ht[,ifi],col="red",lty=1,type = "b",pch=legendPch[ifi],cex=0.8)
    points(fluxes.mcrab,percent.BPA1_ht_cnsrv[,ifi],col="grey",lty=1,type = "b",pch=legendPch[ifi],cex=0.8)
    points(fluxes.mcrab,percent.BPA1_lt[,ifi],col="red",lty=2,type = "b",pch=legendPch[ifi],cex=0.8)
    points(fluxes.mcrab,percent.BPAD_ht[,ifi],col="blue",lty=1,type = "b",pch=legendPch[ifi],cex=0.8)
    points(fluxes.mcrab,percent.BPAD_lt[,ifi],col="blue",lty=2,type = "b",pch=legendPch[ifi],cex=0.8)
    legend(legendPos[ifi],legend=c("A1, high threshold","A1, high threshold\n(200-300samples)",
                               "A1, low threshold","AD, high threshold","AD, low threshold"), 
           col=c("red","grey","red","blue","blue"), lty=c(1,1,2,1,2), 
           pch=c(legendPch[ifi],legendPch[ifi],legendPch[ifi],legendPch[ifi],legendPch[ifi]),cex=0.7,bty="n",
           title=legendTitles[ifi])
}
dev.off()
