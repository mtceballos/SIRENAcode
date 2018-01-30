#
# Calculate Frequency of pairs given a spectrum
#
rm(list=ls())
library(FITSio)
library(Hmisc)
source("~/R/Rfunctions/drawLogPlotBox.r")
setwd("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/baselineLPA2")
#------------------------
# THINGS THAT CAN CHANGE
#------------------------
timeRes <- 1E-6 #s
libTemplates <- "" # "" for full library or "_SHORT" for reduced library
jitter <- "_jitter"
#jitter=""
# Energies where detection matrices are calculated
EprimsMatrices <- c("0.2","1","2","2.5","4","4.5","6","6.5","8") 
EsecsMatrices  <- c("0.2","0.5",
                    "1","1.3","1.6","2","2.3","2.6",
                    "3","3.3","3.6","4","4.3","4.6",
                    "5","5.3","5.6","6","6.3","6.6",
                    "7","7.3","7.6","8") 
# Energies where baglots are calculated
EprimsBPs <- c("0.2","0.5","1","2","3","4","5","6","7","8") 
EsecsBPs  <- c("0.2","0.5","1","2","3","4","5","6","7","8") 
filterLengths <- c("4096", "512", "256") # filters with which events should have been reconstructed (by GRADE)
samprate <- 156250 # Hz-1 - sampling rate
XT <- "XT" # "XT" crosstalk for filenames
# separations for pairs in Detection Matrices
separationsMatrices <- c(4,5,7,10,14,20,28,39,54,75,105,146,202,281,389,540,749,1039,1442,2000) # log scale
# Bagplots Failure areas (as a funtion of secondary energy)
BP.failRange <- array(data=NA, dim=c(length(EprimsBPs),length(EsecsBPs),length(filterLengths),2), 
                      dimnames = list(EprimsBPs,EsecsBPs,filterLengths, c("min","max")))
# Esec=0.2 keV
BP.failRange[    ,1,1,1] <- "250"
BP.failRange[1:6 ,1,1,2] <- "350"
BP.failRange[7:10,1,1,2] <- "330"

BP.failRange[    ,1,2,1] <- "510"
BP.failRange[    ,1,2,2] <- "512"

BP.failRange[    ,1,3,1] <- "225"
BP.failRange[    ,1,3,2] <- "256"

# Esec=0.5 keV
BP.failRange[    ,2,1,1] <- "265"
BP.failRange[1:6 ,2,1,2] <- "305"
BP.failRange[7:10,2,1,2] <- "300"

BP.failRange[    ,2,2,1] <- "510"
BP.failRange[    ,2,2,2] <- "512"

BP.failRange[    ,2,3,1] <- "225"   # but no problem with 00240 (?)
BP.failRange[    ,2,3,2] <- "256"

# Esec=1 keV
BP.failRange[    ,3,1,1] <- "270"
BP.failRange[1:6 ,3,1,2] <- "300"
BP.failRange[7:10,3,1,2] <- "290"

BP.failRange[    ,3,2,1] <- "510"
BP.failRange[    ,3,2,2] <- "512"

BP.failRange[    ,3,3,1] <- "225"   # but no problem with 00240 or 00250 (?)
BP.failRange[    ,3,3,2] <- "256"

# Esec=2 keV
BP.failRange[    ,4,1,1] <- "275"
BP.failRange[    ,4,1,2] <- "285"

BP.failRange[    ,4,2,1] <- "510"
BP.failRange[    ,4,2,2] <- "512"

BP.failRange[1:2 ,4,3,1] <- "225"   # but no problem with 00240 or 00250 (?)
BP.failRange[3:10,4,3,1] <- "255" 
BP.failRange[    ,4,3,2] <- "256"

# Esec=3 keV
BP.failRange[    ,5,1,1] <- "275"
BP.failRange[1:8 ,5,1,2] <- "285"
BP.failRange[9:10,5,1,2] <- "280"

BP.failRange[1:2 ,5,2,1] <- "510"
BP.failRange[3:10,5,2,1] <- "512"
BP.failRange[    ,5,2,2] <- "512"

BP.failRange[    ,5,3,1] <- "256"
BP.failRange[    ,5,3,2] <- "256"
# Esec=4 keV
BP.failRange[    ,6,1,1] <- "272"
BP.failRange[1   ,6,1,2] <- "285"
BP.failRange[2:10,6,1,2] <- "280"

BP.failRange[1:2 ,6,2,1] <- "510"
BP.failRange[3:10,6,2,1] <- "512"
BP.failRange[    ,6,2,2] <- "512"

BP.failRange[    ,6,3,1] <- "256"
BP.failRange[    ,6,3,2] <- "256"
# Esec=5 keV
BP.failRange[    ,7,1,1] <- "272"
BP.failRange[    ,7,1,2] <- "280"

BP.failRange[1:2 ,7,2,1] <- "510"
BP.failRange[3:10,7,2,1] <- "512"
BP.failRange[    ,7,2,2] <- "512"

BP.failRange[    ,7,3,1] <- "256"
BP.failRange[    ,7,3,2] <- "256"
# Esec=6 keV
BP.failRange[    ,8,1,1] <- "270"
BP.failRange[    ,8,1,2] <- "275"

BP.failRange[1:2 ,8,2,1] <- "510"
BP.failRange[3:10,8,2,1] <- "512"
BP.failRange[    ,8,2,2] <- "512"

BP.failRange[    ,8,3,1] <- "256"
BP.failRange[    ,8,3,2] <- "256"

# Esec=7 keV
BP.failRange[1:7 ,9,1,1] <- "270"
BP.failRange[8:10,9,1,1] <- "265"
BP.failRange[    ,9,1,2] <- "275"

BP.failRange[1:2 ,9,2,1] <- "510"
BP.failRange[3:10,9,2,1] <- "512"
BP.failRange[    ,9,2,2] <- "512"

BP.failRange[    ,9,3,1] <- "256"
BP.failRange[    ,9,3,2] <- "256"

# Esec=8 keV
BP.failRange[    ,10,1,1] <- "265"
BP.failRange[1:9 ,10,1,2] <- "275"
BP.failRange[  10,10,1,2] <- "270"

BP.failRange[1:2 ,10,2,1] <- "510"
BP.failRange[3:10,10,2,1] <- "512"
BP.failRange[    ,10,2,2] <- "512"

BP.failRange[    ,10,3,1] <- "256"
BP.failRange[    ,10,3,2] <- "256"

# BP.failRange <- list("0.2"=c(250,360),"0.5"=c(200,300),
#                       "1"=c(240,260),"2"=c(240,260),"3"=c(240,260),
#                       "4"=c(235,245),"5"=c(235,245),"6"=c(235,245),
#                       "7"=c(235,245),"8"=c(235,245))

#BP.failRange <- list("0.2"=c(200,300), "1"=c(240,260),"4"=c(235,245))
#
pdf(paste("e2e/FrequencyPairsLogScale",jitter,XT,libTemplates,".pdf",sep=""),width=10, height=7,version="1.4")
#pdf(paste("e2e/FrequencyPairsLogScaleNOConserv02",XT,libTemplates,".pdf",sep=""),width=10, height=7,version="1.4")
#------------------------
# END THINGS THAT CAN CHANGE
#------------------------

EkeVrecons <- list()
evt.times  <- numeric()
evt.phids  <- numeric()
evt.pixids <- numeric()
phs.in.pix <- numeric()
pix.phids  <- numeric()

nSec <- length(EsecsMatrices)
nPrim <- length(EprimsMatrices)
#separations <- c(4,10,20,40,44,60,90,100,120,200,250,300,400,500,600,800,1000,1600,2000)
nseps <- length(separationsMatrices)
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
percentMiss   <- array(data=NA, dim=c(nFluxes,nfilters,nmodes),dimnames = list(rep("",nFluxes), c("BeFilter", "NoFilter"), c("AD","A1")))
percentBPfail <- array(data=NA, dim=c(nFluxes,nfilters,nmodes),dimnames = list(rep("",nFluxes), c("BeFilter", "NoFilter"), c("AD","A1")))
probDetMatrix <- array(data=NA, dim=c(nSec,nseps,nmodes,nPrim), dimnames = list(EsecsMatrices, separationsMatrices, 
                                                                                detectionModes, EprimsMatrices))
probDetMatrixOLD <- array(data=NA, dim=c(nSec,nseps,nmodes,nPrim), dimnames = list(EsecsMatrices, separationsMatrices, 
                                                                                detectionModes, EprimsMatrices))
# Read probability matrices

for (i1 in 1:nPrim){
    for(idet in 1:nmodes){
        matrixFile <- paste("imageMatrix_",detectionModes[idet],"_",EprimsMatrices[i1],"keV",libTemplates,".mat",sep="")
        load(matrixFile)
        for(i2 in 1:nSec){
            for(is in 1:nseps){
                probDetMatrixOLD[i2,is,idet,i1] <- mat_detectedPulses[i2,is]                
            }
        }
    }
}
for (i1 in 1:nPrim){
    for(idet in 1:nmodes){
        #matrixFile <- paste("imageMatrix_",detectionModes[idet],"_",EprimsMatrices[i1],"keV",libTemplates,".mat",sep="")
        matrixFile <- paste("D_NewimageMatrix_",detectionModes[idet],"_",EprimsMatrices[i1],"keV",libTemplates,jitter,".mat",sep="")
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
        # READ PIXIMPACT file
        pp <- file(description = pixFile, open = "rb")
        header0 <- readFITSheader(pp, fixHdr = 'remove') # read primary header
        header <- readFITSheader(pp, fixHdr = 'remove') # read extension header
        pp.hdr <- parseHdr(header)
        evtTablePIX <- readFITSbintable(pp, header)
        close(pp)
        nsimsPIX <-as.numeric(pp.hdr[which(pp.hdr=="NAXIS2")+1]) # source photons simulated (in piximpact)
        # list of piximpact photons ids
        idcol.ph_id.pix <- which(evtTablePIX$colNames == "PH_ID")
        pix.phids <- evtTablePIX$col[[idcol.ph_id.pix]][1:nsimsPIX]
        nimps <- max(pix.phids)
        # piximpact energies
        idcol.energy.pix <- which(evtTablePIX$colNames == "ENERGY")
        #          setNames(as.list(values),                                           names)
        EkeVsim <- setNames(as.list(evtTablePIX$col[[idcol.energy.pix]][1:nsimsPIX]), pix.phids)
            
        # READ EVT file
        cat("Working with ", evtFile,"\n")
        zz <- file(description = evtFile, open = "rb")
        header0  <- readFITSheader(zz, fixHdr = 'remove') # read primary header
        header   <- readFITSheader(zz, fixHdr = 'remove') # read extension header
        zz.hdr   <- parseHdr(header)
        evtTable <- readFITSbintable(zz, header)
        close(zz)
        nsims <-as.numeric(zz.hdr[which(zz.hdr=="NAXIS2")+1]) # simulated photons in EVT
        #stopifnot(nsimsPIX==nsims) 
        
        # list of EVT photons ids
        idcol.ph_id  <- which(evtTable$colNames == "PH_ID")
        evt.phids    <- evtTable$col[[idcol.ph_id]][1:nsims]
        # EVT energies (as list)
        idcol.signal <- which(evtTable$colNames == "SIGNAL")
        EkeVrecons   <- setNames(as.list(evtTable$col[[idcol.signal]][1:nsims]), evt.phids)
        idcol.time   <- which(evtTable$colNames == "TIME")
        evt.times    <- setNames(as.list(evtTable$col[[idcol.time]][1:nsims]), evt.phids)
        idcol.pixid  <- which(evtTable$colNames == "PIXID")
        evt.pixids   <- setNames(as.list(evtTable$col[[idcol.pixid]][1:nsims]), evt.phids)
        
        pixels <- unique(unlist(evt.pixids, use.names = F))
        
        for(idet in 1:nmodes){
            detMod <- detectionModes[idet]
            cat(" Detection mode: ",detectionModes[idet],"\n")    
            cat("   Photons loop (nsims=",nsims,"- nimps (max PH_ID in piximpact)=",nimps,")\n")
            probDetPh  <- setNames(as.list(rep(1,nsims)), evt.phids)
            probBPfail <- setNames(as.list(rep(0,nsims)), evt.phids)
            
            for (ipix in pixels){
                # which photons are in a given pixel? Time ordered list
                phs.in.pix <-as.character(evtTable$col[[idcol.ph_id]][evtTable$col[[idcol.pixid]]==ipix])
                #phs.in.pix <-as.character(evtTable$col[[idcol.ph_id]][evt.pixids==ipix]) # is the same?
                
                # how many?
                nphs <- length(phs.in.pix)
                # cat("     Pixel=", ipix, "  N.Photons=",nphs, "\n")
                # their initial probs
                if(nphs == 1) next
                for (ip in 2:nphs){
                    phName <- phs.in.pix[ip]
                    phNamePrev <- phs.in.pix[ip-1]

                    # save simulated energy instead of 0.00 energy
                    if(EkeVrecons[[phName]] == 0.0) EkeVrecons[[phName]] <- EkeVsim[[phName]]
                    BPfailsDet <- FALSE
                    BPfailsNoDet <- FALSE
                    
                    diffSamplesPrev <- (evt.times[[phName]]-evt.times[[phNamePrev]])*samprate
                    
                    # look for the closest energy of photons Prim & Sec and closest separation 
                    # in matrices 
                    idclosPrimMat <- whichClosest(EprimsMatrices,EkeVrecons[[phNamePrev]])
                    closestEnergyPrim <- EprimsMatrices[idclosPrimMat]
                    idclosSecMat  <- whichClosest(EsecsMatrices,EkeVrecons[[phName]])
                    closestEnergySec <- EsecsMatrices[idclosSecMat]
                    idclosSepMat  <- whichClosest(separationsMatrices,diffSamplesPrev)
                    closestSep <- separationsMatrices[idclosSepMat]
                    
                    
                    probDet    <- probDetMatrix[idclosSecMat,idclosSepMat,idet,idclosPrimMat] # 0-100% ip-1 <-> ip 
                    probDetOLD <- probDetMatrixOLD[idclosSecMat,idclosSepMat,idet,idclosPrimMat] # 0-100% ip-1 <-> ip 
                    
                    # if(probDetOLD != probDet){
                    #     cat("ProbDetOLD=",probDetOLD," ProbDet=",probDet,
                    #         "diffsamplesPrev=",diffSamplesPrev,"[",idclosSepMat,"]",
                    #         "CloseEprimMat=",closestEnergyPrim,"[",idclosPrimMat,"]",
                    #         "CloseEsecMat=",closestEnergySec,"[",idclosSecMat,"]","\n")
                    # }
                    probDet2 <- 0.
                    probDet2OLD <- 0.
                    
                    # same but with ip-2 (just in case it has to be taken into account)
                    
                    if(ip>2){
                        phNamePrevPrev <- phs.in.pix[ip-2]
                        diffSamplesPrevPrev <- (evt.times[[phName]]-evt.times[[phNamePrevPrev]])*samprate
                        diffSamplesPrevPrevPrev <- (evt.times[[phNamePrev]]-evt.times[[phNamePrevPrev]])*samprate
                        idclosPrim2Mat <- whichClosest(EprimsMatrices,EkeVrecons[ip-2])
                        closestEnergyPrim2 <- EprimsMatrices[idclosPrim2Mat]
                        idclosSep2Mat <- whichClosest(separationsMatrices,diffSamplesPrevPrev)
                        closestSep2 <- separationsMatrices[idclosSep2Mat]
                        probDet2    <- probDetMatrix[idclosSecMat,idclosSep2Mat,idet,idclosPrim2Mat] # 0-100% ip-2 <->ip
                        probDet2OLD <- probDetMatrixOLD[idclosSecMat,idclosSep2Mat,idet,idclosPrim2Mat] # 0-100% ip-2 <->ip
                        
                    }
                    
                    # just to compare with OLD matrices
                    #probDet <-probDetOLD
                    #probDet2 <-probDet2OLD
                    
                    if(diffSamplesPrev <= timeRes) probDet <- 0.
                    # ip = secondary
                    # ip -1 = primary
                    #======================================================================
                    # prob=140%     
                    #     Prim (50pulses)           Secondaries (50 pulses)
                    #         (50)ok                         (10) ok -> 10/50=20%  (well)         
                    #                                   (40+40 very close) -> 40/50=80% (bad)
                    # prob=70%     
                    #     Prim (50pulses)           Secondaries (50 pulses)
                    #         (50)ok                    (20) ok -> 20/50=40%  (well)         
                    #                                   (30) missing-> 30/50=60% (bad)
                    # prob=30% (if Eprim>0.3 keV) # missing secondaries
                    #     Prim (50pulses)           Secondaries (50 pulses)
                    #         (30)ok                    (0) ok -> 0/50=0%  (well)         
                    #         ->30/50=60%               (50) missing-> 50/50=100% (bad)    
                    # prob=30% (if Eprim <= 0.3keV) # missing primaries
                    #     Prim (50pulses)           Secondaries (50 pulses)
                    #     (0) ok -> 0/50=0%(well)       (30) ok -> 30/50=60% (well)
                    #     (50) missing->            
                    #          50/50=100% (bad)    
                
                    
                    # Computed Probs in 0-1 scale (matrix probs in 0-100%)
                    
                    # ip-1 <--> ip
                    if(probDet>=50 && probDet <=102){
                        probDetPrim <- 1.
                        probMissPrim <- 0.                       # primary
                        probMissSec <- abs(probDet-100)/50.   # secondary
                        probDetSec <- 1-probMissSec 
                    }else if (probDet < 50){
                        if(closestEnergyPrim > 0.3){
                            probDetPrim <- probDet/50.
                            probMissPrim <- 1.-probDetPrim      # primary
                            probDetSec <- 0.
                            probMissSec <- 1. #secondary
                        }else{
                            probDetPrim <- 0.
                            probMissPrim <- 1.      # primary
                            probDetSec <- probDet/50.
                            probMissSec <- 1.-probDetPrim   #secondary
                        }
                    }else if (probDet >102){
                        probDetPrim <- 0.
                        probDetSec <- 0.
                    }
                    # ip-2 <--> ip
                    if(probDet2>=50 && probDet2 <=102){
                        probMissSec2 <- abs(probDet2-100)/50.   # secondary
                        probDetSec2 <- 1-probMissSec2 
                    }else if (probDet2 < 50){
                        probDetSec2 <- 0.
                        probMissSec2 <- 1. #secondary
                    }else if (probDet2 >102){
                        probDetSec2 <- 0.
                    }
                    probDetPh[[phName]] <- probDetPh[[phNamePrev]]*probDetPrim*probDetSec 
                    #probDetPh[[phNamePrev]] <- probDetPh[[phNamePrev]]*probDetPrim 
                    probDetPh[[phName]] <- probDetPh[[phName]] + (1-probDetPh[[phNamePrev]])*probDetSec2
                    probDetPh[[phNamePrev]] <- probDetPh[[phNamePrev]]*probDetPrim 
                    
                    # Once known phName -> look for the filter used for phNamePrev
                    #=============================================================
                    # => calculate bagplots (a posteriori)
                    if(ip>2){
                        # 1) look for the closest energy of photons Prim & Sec  in bagplots
                        idclosPrimBP <- whichClosest(EprimsBPs,EkeVrecons[[phNamePrevPrev]])
                        idclosSecBP  <- whichClosest(EsecsBPs,EkeVrecons[[phNamePrev]])
                        
                        # 2) For probDetPh[phName] times -> filter with length=diffSamplesPrevPrev will be used for 
                        # (phNamePrevPrev,phNamePrev) pair
                        if(diffSamplesPrevPrev < as.numeric(min(filterLengths))) {
                             closestFilDet <- "LRes"
                        }else{
                             closestFilDet <- max(filterLengths[as.numeric(filterLengths)<=diffSamplesPrevPrev])
                        }
                        idclosFilDet  <- which(filterLengths == closestFilDet)
                        #    For 1-probDetPh[phName] times -> filter with length=HighRes will be used for phNamePrev
                        closestFilNoDet <- max(filterLengths)
                        idclosFilNoDet  <- which(filterLengths == closestFilNoDet)
                        
                        # # 3) BPA1AD: Pulses not rejected by bagplots and not detected by A1/AD
                        # #===============================================================
                        # # if in same pixel than previous && in required temporal range (bagplots) for given
                        # #    energy prim/sec. 
                        # 
                        if(closestFilDet == "LRes"){
                             BPfailsDet <- FALSE  # flagged as Lres and possibly not reconstructed
                        }else if((diffSamplesPrevPrev>=BP.failRange[idclosPrimBP,idclosSecBP,idclosFilDet,1] && 
                             diffSamplesPrevPrev<=BP.failRange[idclosPrimBP,idclosSecBP,idclosFilDet,2])){
                             BPfailsDet <- TRUE
                        }
                        if(diffSamplesPrevPrevPrev>=closestFilNoDet) { #no problem: phName does not interfere
                            if((diffSamplesPrevPrevPrev>=BP.failRange[idclosPrimBP,idclosSecBP,idclosFilNoDet,1] && 
                                diffSamplesPrevPrevPrev<=BP.failRange[idclosPrimBP,idclosSecBP,idclosFilNoDet,2]))
                                BPfailsNoDet <- TRUE
                        }else{ # pulse phName is not detected and enters in phNamePrevPrev-phNamePrev pair
                                # bagplots do not consider triple pulses...be conservative!
                            BPfailsNoDet <- TRUE
                        }     
                        if(BPfailsDet){ # if not in 'bad area' probBPfail=0
                             probBPfail[[phNamePrev]] <- probDetPh[[phName]]*(1-probDetPh[[phNamePrev]])
                        }
                        if(BPfailsNoDet){ # if not in 'bad area' probBPfail=0
                            probBPfail[[phNamePrev]] <- probBPfail[[phNamePrev]] + 
                                (1-probDetPh[[phName]])*(1-probDetPh[[phNamePrev]]) # just in case it's changed
                        }
                    }
                    # 
                    # # Discover with which filter this 'undetected' pair would be reconstructed
                    # if (ip<nphs) {
                    #     phNameNext <- phs.in.pix[ip+1]
                    #     diffSamples3 <- (evt.times[[phNameNext]] - evt.times[[phName]])*samprate
                    #     if(diffSamples3 < as.numeric(min(filterLengths))) {
                    #         closestFil <- "LRes"
                    #     }else{
                    #         closestFil <- max(filterLengths[as.numeric(filterLengths)<=diffSamples3])
                    #     }
                    # }else{
                    #     closestFil <- filterLengths[1]
                    # }
                    # idclosFil  <- which(filterLengths == closestFil)
                    # 
                    # 
                    
                } #photons in pixel 
            } #pixels (given a ctrate and detection mode)
            percentMiss[ifl,ifi,idet] <- sum(1-unlist(probDetPh, use.names = F))/nsims # 0-1
            percentBPfail[ifl,ifi,idet] <- sum(unlist(probBPfail, use.names = F))/nsims # 0-1
        } #detMod
    } #each flux
}#each filter


# PLOTTING (UN)DETECTION
##########################
par(mfrow=c(1,2))
subtit=paste("//",jitter,"//",XT)
y1max <- 6
y2max <- 0.06
if (XT == "XT") y2max <- 0.1
if(libTemplates == "_SHORT") {
    subtit <- "(SHORT lib templates)"
    y1max <- 80
    y2max <- 1
}
drawLogPlotBox(xlimits=c(1E-4,1E3),ylimits=c(0,y1max), logxy="x",xlabel="Intensity (mCrab)",
               ylabel=expression("Fraction of photons (%)"), naxes=c(T,T,F,F))
title(main=paste("Undetected Fraction of Photons\n Crab spectrum",subtit,"\n"), cex.main=0.8,
      sub="1 mCrab=2.11E-11 erg/cm2/s 2-10kev",cex.sub=0.6)
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

text(1E-3,4,"PSF=athena_psf_onaxis_20150602.fits",cex=0.4)
text(30,4,"PSF=\nathena_ladapt_defocus_35mm_kev_psf_20161013.fits",cex=0.4)
# PLOTTING (UN)DETECTION * BAGPLTS flagging
##############################################
drawLogPlotBox(xlimits=c(1E-4,1E3),ylimits=c(0,y2max), logxy="x",xlabel="Intensity (mCrab)", 
               ylabel=expression("Fraction of photons (%)"), naxes=c(T,T,F,F))
title(main=paste("Undetected Fraction of Photons \n(not spotted by derivative)",subtit,sep="\n"),
      cex.main=0.8)
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
text(1E-3,0.04,"PSF=athena_psf_onaxis_20150602.fits",cex=0.4)
text(30,0.04,"PSF=\nathena_ladapt_defocus_35mm_kev_psf_20161013.fits",cex=0.4)

dev.off()

