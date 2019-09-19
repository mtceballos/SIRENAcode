rm(list=ls())
library(FITSio)
library(Hmisc)
source("~/R/Rfunctions/drawLogPlotBox.r")
setwd("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/baselineLPA2")


is.closer.than <- function(vector,value,distance){
    # return TRUE if 'value' in closer than 'distance' 
    # (abs value) to any element in 'vector'
    res <- FALSE
    closest <- whichClosest(vector,value)
    res <- abs(vector[closest] - value) <= distance
    
}

expos=500 # 100, 1000 s
timeRes <- 1E-6 #s
libTemplates <- "" # "" for full library or "_SHORT" for reduced library
jitter <- "_jitter" # or ''
samprateStr <- ""      # or '_samprate2'
samprateStr <- "_samprate2"

simDir=paste("e2e/expos",expos,samprateStr,sep="")
pdffile <- paste(simDir,"/FrequencyPairsLogScale",jitter,samprateStr,libTemplates,".pdf",sep="")
# Energies where detection (AD/A1) matrices are calculated
EprimsMatrices <- c("0.2","1","2","2.5","4","4.5","6","6.5","8") 
EsecsMatrices  <- c("0.2","0.5","1","1.3","1.6","2","2.3","2.6",
                    "3","3.3","3.6","4","4.3","4.6","5","5.3",
                    "5.6","6","6.3","6.6","7","7.3","7.6","8") 
nSec <- length(EsecsMatrices)
nPrim <- length(EprimsMatrices)
# separations for pairs in Detection Matrices
if(samprateStr == ""){
    separationsMatrices <- c(4,5,7,10,14,20,28,39,54,75,105,146,202,281,
                         389,540,749,1039,1442,2000) # log scale
    filterLengths <- c("8192", "512", "256") # filters 
    samprate <- 156250 # Hz-1 - sampling rate
}else if(samprateStr == "_samprate2"){
    separationsMatrices <- c(2,3,5,7,10,14,19,27,37,52,73,101,140,194,
                         270,374,519,721,1000) # log scale
    filterLengths <- c("4096", "256", "128") # filters 
    samprate <- 78125 # Hz-1 - sampling rate
}
nseps <- length(separationsMatrices)

# Energies where baglots are calculated
EprimsBPs <- c("0.2","0.5","1","2","3","4","5","6","7","8") 
EsecsBPs  <- c("0.2","0.5","1","2","3","4","5","6","7","8") 


EkeVrecons <- list()
evt.times  <- numeric()
evt.phids  <- numeric()
evt.pixids <- numeric()
phs.in.pix <- numeric()
pix.phids  <- numeric()

detectionModes <- c("AD","A1")
nmodes <- length(detectionModes)
#flux210_mCrab=0.000000000021147 = 2.1147E-11 erg/cm2/s 2-10kev   ---> 1 mCrab 
fluxes.mcrab <- c("0.0001", "0.0005", "0.001", "0.005", "0.01", "0.036", "0.13", "0.46",
                  "0.60", "0.80", "1.67", "6.", "21.54", "77.43", "278.26", "1000.") # mCrab
nFluxes <- length(fluxes.mcrab)
filters <- c("_filter_", "_") # Be Filter & w/o filter
nfilters <- length(filters)
# info for the plots
legendTitles<-c("Be Filter","NO Be Filter")
legendPos <- c("topleft","top")
legendPch <- c(4,1)


#Create array 
#BP.failAreasStr <- array(data=NA, dim=c(length(EprimsBPs),length(EsecsBPs),length(filterLengths)), 
#                      dimnames = list(EprimsBPs,EsecsBPs,filterLengths))
#nsepsBP <- numeric(length(EsecsBPs))

BP.failCrossMaxStr <- array(data=NA, dim=c(length(EprimsBPs),length(EsecsBPs),length(filterLengths)), 
                         dimnames = list(EprimsBPs,EsecsBPs,filterLengths))
BP.failCrossMinStr <- array(data=NA, dim=c(length(EprimsBPs),length(EsecsBPs),length(filterLengths)), 
                         dimnames = list(EprimsBPs,EsecsBPs,filterLengths))

# One file per secondary energy
for (ieSec in 1:length(EsecsBPs)){
    ## for each secondary energy: read list of problematic bagplots
    Esec <- EsecsBPs[ieSec]
    fileBPfail <- paste("BPfail_",Esec,"keV",samprateStr,".dat",sep="")
    load(fileBPfail)
    for (ie in 1:length(EprimsBPs)){
     for(ig in 1:length(filterLengths)){
        BP.failCrossMinStr[ie,ieSec,ig] <- rootsMin[ie,ig] # string with crosses with bandMin
        BP.failCrossMaxStr[ie,ieSec,ig] <- rootsMax[ie,ig] # string with crosses with bandMax
     }
    }
    
    # # load separationsToPlot from derivative.R
    # load(fileBPfail)
    # for (ie in 1:length(EprimsBPs)){
    #     for(ifl in 1:length(filterLengths)){
    #         BP.failAreasStr[ie,ieSec,ifl] <- separationsToPlot[ie,ifl]
    #         #cat("BP fail areas:",separationsToPlot[ie,ifl])
    #     }
    # }
    # # for each secondary energy: read list of full separations in bagplots
    # fileBPseps <- paste("BPseps_",Esec,"keV",samprateStr,".dat",sep="")
    # load(fileBPseps)
    # nsepsBP[ieSec] <- length(separationsBGplots)
    # BP.seps <- array(data=NA, dim=c(length(EsecsBPs),nsepsBP), 
    #                  dimnames = list(EsecsBPs,separationsBGplots))
    # for(is in 1:nsepsBP[ieSec]){
    #     BP.seps[ieSec,is] <- separationsBGplots[is]
    # }
}


# Initialize probability arrays (percentage of photons)
namesPercent  <- list(rep("",nFluxes), c("BeFilter", "NoFilter"), c("AD","A1"))
namesProb     <- list(EsecsMatrices, separationsMatrices,detectionModes, EprimsMatrices,
                      c("Prim","Sec"))
photonsMiss   <- array(data=NA, dim=c(nFluxes,nfilters,nmodes), dimnames = namesPercent)
photonsBPfail <- array(data=NA, dim=c(nFluxes,nfilters,nmodes), dimnames = namesPercent)
enerBPfail    <- array(data=999, dim=c(nFluxes,nfilters,nmodes), dimnames = namesPercent)
probDetMatrix <- array(data=NA, dim=c(nSec,nseps,nmodes,nPrim,2), dimnames = namesProb)
 

# Read probability matrices (0-100% --> 0-1)
for (i1 in 1:nPrim){
    for(idet in 1:nmodes){
        #primary probabilities
        matrixFilePrim <- paste("SakaiMaps/D_imageMatrix_",detectionModes[idet],"_",EprimsMatrices[i1],"keV",
                                libTemplates,jitter,"P.mat",sep="")
        load(matrixFilePrim)
        for(i2 in 1:nSec){
            for(is in 1:nseps){
                probDetMatrix[i2,is,idet,i1,1] <- mat_detectedPulsesP[i2,is]/100.                 
            }
        }
        #secondary probabilities
        matrixFileSec <- paste("SakaiMaps/D_imageMatrix_",detectionModes[idet],"_",EprimsMatrices[i1],"keV",
                               libTemplates,jitter,"S.mat",sep="")
        load(matrixFileSec)
        for(i2 in 1:nSec){
            for(is in 1:nseps){
                probDetMatrix[i2,is,idet,i1,2] <- mat_detectedPulsesS[i2,is]/100.
            }
        }
    }
}


# Start processing
#==================
# allocate space for number of simulated photons
nsimsMat <- matrix(data=NA,nrow=nfilters, ncol=nFluxes)

for (ifi in 1:nfilters){
    cat("Working with filter ", filters[ifi],"\n")
    filter <- filters[ifi]
    for (ifl in 1:nFluxes){
        if(fluxes.mcrab[ifl] < 0.25) {
            # read evt file
            evtFile <- paste(simDir,"/crabSpec",fluxes.mcrab[ifl],"mCrab.fits",sep="")
            pixFile <- paste(simDir,"/crabSpec",fluxes.mcrab[ifl],"mCrab.piximpact",sep="")
        }else{
            
            evtFile <- paste(simDir,"/crabSpec",fluxes.mcrab[ifl],"mCrab",filter,"35mm.fits",sep="")
            pixFile <- paste(simDir,"/crabSpec",fluxes.mcrab[ifl],"mCrab",filter,"35mm.piximpact",sep="")
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
        nsimsMat[ifi,ifl] <-as.numeric(zz.hdr[which(zz.hdr=="NAXIS2")+1]) # simulated photons in EVT
        nsims <- nsimsMat[ifi,ifl] 
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
            #cat(" Detection mode: ",detectionModes[idet],"\n")    
            #cat("   Photons loop (nsims=",nsims,"- nimps (max PH_ID in piximpact)=",nimps,")\n")
            probDetPh  <- setNames(as.list(rep(1,nsims)), evt.phids)
            probBPfail <- setNames(as.list(rep(0,nsims)), evt.phids)
            
            enerBPfailVec <- numeric()
            for (ipix in pixels){
                # which photons are in a given pixel? Time ordered list
                phs.in.pix <-as.character(evtTable$col[[idcol.ph_id]][evtTable$col[[idcol.pixid]]==ipix])
                #phs.in.pix <-as.character(evtTable$col[[idcol.ph_id]][evt.pixids==ipix]) # is the same?
                
                # how many?
                nphs <- length(phs.in.pix)
                # their initial probs
                if(nphs <= 1) next
                #cat("detMode=",detectionModes[idet],"  Flux=",fluxes.mcrab[ifl]," Pixel=", ipix, "  N.Photons=",nphs, "\n")
                for (ip in 2:nphs){
                    #cat("Working with ph=",ip,"\n")
                    
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
                    
                    if(diffSamplesPrev > max(as.numeric(filterLengths))){
                        probDetPrim <- 1.
                        probDetSec  <- 1.
                        if(detectionModes[idet] == "AD" && closestEnergySec == 0.2) {
                            if(samprateStr == ""){
                                probDetPrim <- 0.8075
                                probDetSec  <- 0.8075
                            }else if (samprateStr == "_samprate2"){
                                probDetPrim <- 0.9995
                                probDetSec  <- 0.9995
                            }
                        }
                    }else{
                        probDetPrim <- probDetMatrix[idclosSecMat,idclosSepMat,idet,idclosPrimMat,1] # 0-1 ip-1 <-> ip
                        probDetSec  <- probDetMatrix[idclosSecMat,idclosSepMat,idet,idclosPrimMat,2] # 0-1 ip-1 <-> ip
                    }
                    # probDetPrim <- probDetMatrix[idclosSecMat,idclosSepMat,idet,idclosPrimMat,1] # 0-1 ip-1 <-> ip 
                    # probDetSec  <- probDetMatrix[idclosSecMat,idclosSepMat,idet,idclosPrimMat,2] # 0-1 ip-1 <-> ip 
                    
                    probDetSec2 <- 1. # 0-1 ip-2 <-> ip 
                    
                    # same but with ip-2 (just in case it has to be taken into account)
                    
                    if(ip>2){
                        phNamePrevPrev <- phs.in.pix[ip-2]
                        diffSamplesPrevPrev <- (evt.times[[phName]]-evt.times[[phNamePrevPrev]])*samprate
                        diffSamplesPrevPrevPrev <- (evt.times[[phNamePrev]]-evt.times[[phNamePrevPrev]])*samprate
                        idclosPrim2Mat <- whichClosest(EprimsMatrices,EkeVrecons[ip-2])
                        closestEnergyPrim2 <- EprimsMatrices[idclosPrim2Mat]
                        idclosSep2Mat <- whichClosest(separationsMatrices,diffSamplesPrevPrev)
                        closestSep2 <- separationsMatrices[idclosSep2Mat]
                        if(diffSamplesPrevPrev > max(as.numeric(filterLengths))){
                            probDetSec2 <- 1.
                        }else{
                            probDetSec2 <- probDetMatrix[idclosSecMat,idclosSep2Mat,idet,idclosPrim2Mat,2] # 0-100% ip-2 <->ip
                        }

                    }
                    
                    if(diffSamplesPrev <= timeRes) probDetSec <- 0.
                    # ip = secondary
                    # ip -1 = primary
                    
                    # Computed Probs in 0-1 scale 
                    probDetPh[[phName]] <- probDetPh[[phNamePrev]]*probDetSec 
                    #probDetPh[[phNamePrev]] <- probDetPh[[phNamePrev]]*probDetPrim 
                    probDetPh[[phName]] <- probDetPh[[phName]] + (1-probDetPh[[phNamePrev]])*probDetSec2
                    probDetPh[[phNamePrev]] <- probDetPh[[phNamePrev]]*probDetPrim 
                        
                    # if(fluxes.mcrab[ifl] == "1000." && ipix==2101 && (phName=="649" || phName=="126")){    
                    # #if(fluxes.mcrab[ifl] == "1000." && 
                    # #   (detectionModes[idet]=="AD" && (probDetPh[[phName]]<1 ||probDetPh[[phNamePrev]]<1 ))){    
                    #     cat("detMode=",detectionModes[idet]," ipix=",ipix,"\n")
                    #     cat("             ip=",phName," E=",EkeVrecons[[phName]]," Prob=", probDetPh[[phName]],"\n")
                    #     cat("             ip-1=",phNamePrev," E=",EkeVrecons[[phNamePrev]]," Prob=", probDetPh[[phNamePrev]],
                    #         " diffSamples",diffSamplesPrev," closestEnergySec=",closestEnergySec,
                    #         " closestEnergyPrim=",closestEnergyPrim,
                    #         " closestSep=",closestSep," probDetPrim=",probDetPrim," probDetSec=",probDetSec,"\n")
                    # }
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
 
                        BPfailsDet <- FALSE
                        BPfailsNoDet <- FALSE
                        
                        # Read/create intervals whith bagplots conflict
                        BP.failCrossMin <- unlist(strsplit(BP.failCrossMinStr[idclosPrimBP,idclosSecBP,idclosFilDet],","))
                        BP.failCrossMin <- as.numeric(BP.failCrossMin[BP.failCrossMin != ""])
                        BP.failCrossMax <- unlist(strsplit(BP.failCrossMaxStr[idclosPrimBP,idclosSecBP,idclosFilDet],","))
                        BP.failCrossMax <- as.numeric(BP.failCrossMax[BP.failCrossMax != ""])
                        nintervals <- length(BP.failCrossMin)
                        if(nintervals > 0){ # if there is at least one interval...
                            intervalMin <- numeric(nintervals)
                            intervalMax <- numeric(nintervals)
                            for (ir in 1:nintervals){
                                if(ir%%2 != 0){
                                    intervalMin[ir] <- BP.failCrossMax[ir]
                                    intervalMax[ir] <- BP.failCrossMin[ir]
                                }else{
                                    intervalMin[ir] <- BP.failCrossMin[ir]
                                    intervalMax[ir] <- BP.failCrossMax[ir]
                                }
                            }
                            irs <- 1:nintervals
                        }
                        
                        # if phName is detected...
                        if(length(BP.failCrossMin)>0){ # if there is at least one interval...
                            if(closestFilDet == "LRes"){
                                BPfailsDet <- FALSE  # flagged as Lres and possibly not reconstructed
                            }else{ 
                                # if distance PrevPrev<-->Prev is in one of the intervals
                                if(any(vapply(irs, function(irs) 
                                    diffSamplesPrevPrevPrev>intervalMin[irs] && 
                                    diffSamplesPrevPrevPrev<intervalMax[irs],logical(1))))
                                    BPfailsDet <- TRUE
                            }
                        }

                        # if phName is not detected...
                        # Read/create intervals whith bagplots conflict
                        BP.failCrossMin <- unlist(strsplit(BP.failCrossMinStr[idclosPrimBP,idclosSecBP,idclosFilNoDet],","))
                        BP.failCrossMin <- as.numeric(BP.failCrossMin[BP.failCrossMin != ""])
                        BP.failCrossMax <- unlist(strsplit(BP.failCrossMaxStr[idclosPrimBP,idclosSecBP,idclosFilNoDet],","))
                        BP.failCrossMax <- as.numeric(BP.failCrossMax[BP.failCrossMax != ""])
                        nintervals <- length(BP.failCrossMin)
                        if(nintervals > 0){ # if there is at least one interval...
                            intervalMin <- numeric(nintervals)
                            intervalMax <- numeric(nintervals)
                            for (ir in 1:nintervals){
                                if(ir%%2 != 0){
                                    intervalMin[ir] <- BP.failCrossMax[ir]
                                    intervalMax[ir] <- BP.failCrossMin[ir]
                                }else{
                                    intervalMin[ir] <- BP.failCrossMin[ir]
                                    intervalMax[ir] <- BP.failCrossMax[ir]
                                }
                            }
                            irs <- 1:nintervals
                        }
                        
                        if(length(BP.failCrossMin)>0){ # if there is at least one interval...   
                            # check if phName would be in the reconstructon distance
                            if((diffSamplesPrevPrevPrev+diffSamplesPrev)>=closestFilNoDet) { #no problem: phName does not interfere
                                # if distance PrevPrev<-->Prev is in one of the intervals
                                if(any(vapply(irs, function(irs) 
                                    diffSamplesPrevPrevPrev>intervalMin[irs] && 
                                    diffSamplesPrevPrevPrev<intervalMax[irs],logical(1))))
                                    BPfailsNoDet <- TRUE
                            }else{ # pulse phName is not detected and enters in phNamePrevPrev-phNamePrev pair
                                # bagplots do not consider triple pulses...be conservative!
                                BPfailsNoDet <- TRUE
                            }     
                        }
                            
                        if(BPfailsDet){ # if not in 'bad area' probBPfail=0
                             probBPfail[[phNamePrev]] <- probDetPh[[phName]]*(1-probDetPh[[phNamePrev]])
                        }
                        if(BPfailsNoDet){ # if not in 'bad area' probBPfail=0
                            probBPfail[[phNamePrev]] <- probBPfail[[phNamePrev]] + 
                                (1-probDetPh[[phName]])*(1-probDetPh[[phNamePrev]]) 
                        }
                        # Store energy of phNamePrev if BPfail
                        if(probBPfail[[phNamePrev]] > 0)
                            enerBPfailVec <- c(enerBPfailVec, EkeVrecons[[phNamePrev]]) 
                    }
                    
                    
                } #photons in pixel 
            } #pixels (given a ctrate and detection mode)
            photonsMiss[ifl,ifi,idet] <- sum(1-unlist(probDetPh, use.names = F)) 
            photonsBPfail[ifl,ifi,idet] <- sum(unlist(probBPfail, use.names = F))
            if(length(enerBPfailVec)>0) enerBPfail[ifl,ifi,idet] <- paste(enerBPfailVec,collapse=" ")
        } #detMod
    } #each flux
}#each filter

save(photonsMiss,file=paste(simDir,"/photonsMiss",samprateStr,jitter,".dat",sep=""))
save(photonsBPfail,file=paste(simDir,"/photonsBPfail",samprateStr,jitter,".dat",sep=""))
save(nsimsMat,file=paste(simDir,"/nsimsMat",samprateStr,jitter,".dat",sep=""))

pdf(pdffile)
par(mfrow=c(2,1))
subtit=paste("//",jitter,"//",samprateStr)

if(libTemplates == "_SHORT") {
    subtit <- "(SHORT lib templates)"
}
yminus <- as.numeric(nFluxes)
yplus  <- as.numeric(nFluxes)

# PLOTTING (UN)DETECTION
##########################
# draw log axes
y1maxAD <- photonsMiss[nFluxes,2,1]/nsimsMat[2,nFluxes]*100 + 
    sqrt(photonsMiss[nFluxes,2,1])/nsimsMat[2,nFluxes]*100
y1maxA1 <- photonsMiss[nFluxes,2,2]/nsimsMat[2,nFluxes]*100 + 
    sqrt(photonsMiss[nFluxes,2,2])/nsimsMat[2,nFluxes]*100
y1max <- max(y1maxAD, y1maxA1)
drawLogPlotBox(xlimits=c(1E-4,1E3),ylimits=c(0,y1max), logxy="x",xlabel="Intensity (mCrab)",
               ylabel=expression("Fraction of photons (%)"), naxes=c(T,T,F,F))
title(main=paste("Undetected Fraction of Photons\n Crab spectrum",subtit,"\n"), cex.main=0.8,
      sub="1 mCrab=2.11E-11 erg/cm2/s 2-10kev",cex.sub=0.6)
# draw limit of difusion
abline(v=0.25,col="grey",lty=2)
for (ifi in 1:nfilters){
    for(idet in 1:nmodes){
        coldet<- "blue"
        if(detectionModes[idet] == "A1") coldet <- "red"
        points(fluxes.mcrab,photonsMiss[,ifi,idet]/nsimsMat[ifi,]*100, 
               col=coldet,type = "b",pch=legendPch[ifi],cex=0.8) 
        yminus <- photonsMiss[,ifi,idet]/nsimsMat[ifi,]*100 - sqrt(photonsMiss[,ifi,idet])/nsimsMat[ifi,]*100
        yplus  <- photonsMiss[,ifi,idet]/nsimsMat[ifi,]*100 + sqrt(photonsMiss[,ifi,idet])/nsimsMat[ifi,]*100
        suppressWarnings(arrows(as.numeric(fluxes.mcrab), yminus,as.numeric(fluxes.mcrab),yplus,
               col=coldet,lty=1,length=0.025, angle=90, code=3))
    }
}
# plot legend
legend("topleft",legend=c(paste(detectionModes[1],legendTitles[1]),
                          paste(detectionModes[2],legendTitles[1]),
                          paste(detectionModes[1],legendTitles[2]),
                          paste(detectionModes[2],legendTitles[2])),        
                col=c("blue","red","blue","red"),
                pch=c(legendPch[1],legendPch[1],legendPch[2],legendPch[2]),
                cex=0.7,bty="n")

text(1E-3,y1max/3.,"PSF=athena_psf_onaxis_20150602.fits",cex=0.4)
text(30,y1max/2.,"PSF=\nathena_ladapt_defocus_35mm_kev_psf_20161013.fits",cex=0.4)


# PLOTTING (UN)DETECTION * BAGPLTS flagging
##############################################
y2maxAD <- max(photonsBPfail[nFluxes,,1]/nsimsMat[,nFluxes]*100 + sqrt(photonsBPfail[nFluxes,,1])/nsimsMat[,nFluxes]*100)
y2maxA1 <- max(photonsBPfail[nFluxes,,2]/nsimsMat[,nFluxes]*100 + sqrt(photonsBPfail[nFluxes,,2])/nsimsMat[,nFluxes]*100)
y2max <- max(y2maxAD,y2maxA1)
drawLogPlotBox(xlimits=c(1E-4,1E3),ylimits=c(0,y2max), logxy="x",xlabel="Intensity (mCrab)", 
               ylabel=expression("Fraction of photons (%)"), naxes=c(T,T,F,F))
title(main=paste("Undetected Fraction of Photons \n(not spotted by derivative)",subtit,sep="\n"),
      cex.main=0.8)
for (ifi in 1:nfilters){
    for(idet in 1:nmodes){
        coldet<- "blue"
        if(detectionModes[idet] == "A1") coldet <- "red"
        points(fluxes.mcrab,photonsBPfail[,ifi,idet]/nsimsMat[ifi,]*100,col=coldet,lty=1,type = "b",pch=legendPch[ifi],cex=0.8)
        yminus <- photonsBPfail[,ifi,idet]/nsimsMat[ifi,]*100 - sqrt(photonsBPfail[,ifi,idet])/nsimsMat[ifi,]*100
        yplus  <- photonsBPfail[,ifi,idet]/nsimsMat[ifi,]*100 + sqrt(photonsBPfail[,ifi,idet])/nsimsMat[ifi,]*100
        suppressWarnings(arrows(as.numeric(fluxes.mcrab), yminus,as.numeric(fluxes.mcrab),yplus,
               col=coldet,lty=1,length=0.025, angle=90, code=3))
        for (ifl in 1:nFluxes){
            photonsInvolved <- as.numeric(strsplit(enerBPfail[ifl,ifi,idet],split=" ")[[1]])
            photonsInvolved <- photonsInvolved[photonsInvolved<999]
            lessThan05 <- sum(photonsInvolved<=0.5)
            if(lessThan05 > 0) {
                cat("Low energy (<=0.5keV) photons involved for (Flux,A1/AD,Filter)=(",
                    fluxes.mcrab[ifl],",",detectionModes[idet],",",legendTitles[ifi],"): ",
                    lessThan05, "/",length(photonsInvolved),"\n")
                # text(fluxes.mcrab[ifl],photonsBPfail[ifl,ifi,idet]/nsimsMat[ifi,ifl]*100+0.05,
                #     col=coldet, cex=0.8)
            }
        }
    }
}
#plot legend
legend("topleft",legend=c(paste(detectionModes[1],legendTitles[1]),
                          paste(detectionModes[2],legendTitles[1]),
                          paste(detectionModes[1],legendTitles[2]),
                          paste(detectionModes[2],legendTitles[2])),        
                          col=c("blue","red","blue","red"),
                          pch=c(legendPch[1],legendPch[1],legendPch[2],legendPch[2]),
                          cex=0.7,bty="n")
abline(v=0.25,col="grey",lty=2)
text(1E-3,y2max/3.,"PSF=athena_psf_onaxis_20150602.fits",cex=0.4)
text(30,y2max/3.,"PSF=\nathena_ladapt_defocus_35mm_kev_psf_20161013.fits",cex=0.4)
dev.off()
