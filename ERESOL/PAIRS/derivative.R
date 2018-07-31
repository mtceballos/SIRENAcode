#
# BAGPLOTS: PACKAGE Documented on derivative.ipynb (also GitHub)
#
rm(list=ls())

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
npulses <- 1000      # Number of pulses at each energy
nPairs <- 500    
samprate <- 156250
samprateStr="" # to name files with full samprate
samprateStr="_samprate2" # to name files with 1/2 samprate
if (samprateStr == "_samprate2"){
    samprate <- samprate/2.    
    filterLengths <- c(4096,256,128)
    pulseLength<- 4096   # pulse length
    xmax <- c(4200, 260, 150 )
    separationsBGplots <- sprintf("%05d",
                           sort(c(15,seq(19,25,2),seq(16,34,2),35,seq(40,256,2),300,
                            seq(320,340,10),seq(360,400,20),450,500,750,800,seq(1000,4000,250),
                            seq(4010,4196,10),4200)))
}else{
    filterLengths <- c(8192,512,256)
    pulseLength<- 8192   # pulse length
    xmax <- c(8200, 550, 260 )
    separationsBGplots <- sprintf("%05d",
        sort(c(seq(30,300,5),42,47,52,57,62,seq(310,500,5),510,seq(520,800,20),900,1000,
               seq(1500,7500,500),seq(7600,8000,100),8050,8100,8150,8192,8250)))
}
fEnergy="6" #keV
energies <- c("0.2", "0.5", "1", "2", "3", "4", "5", "6", "7", "8") # pulses energies
nenergies <- length(energies)
Esec <- args[1] # keV : energy of secondary pulses
nseps <- length(separationsBGplots)
seps.ms <- as.numeric(separationsBGplots)/samprate*1E3 #separations in ms

findCrosses <- function(x1,x2,x3,y1,y2,y3,level){
    # parabola ax^2+bx+c
    A1 <- -x1^2 + x2^2
    B1 <- -x1 + x2
    D1 <- -y1 + y2
    A2 <- -x2^2 + x3^2
    B2 <- -x2 + x3
    D2 <- -y2 + y3
    A3 <- -(B2/B1)*A1 + A2
    D3 <- -(B2/B1)*D1 + D2
    apar <- D3/A3
    bpar <- (D1-A1*apar)/B1
    cpar <- y1 - apar*x1^2 - bpar*x1 - level
    root1 <- (-bpar + sqrt(bpar^2 - 4*apar*cpar))/(2*apar)
    root2 <- (-bpar - sqrt(bpar^2 - 4*apar*cpar))/(2*apar)
    return(c(root1,root2))
}
# --------------------------------------------------------------------------------------

derivArrayMean4samples <- array(data=NA,dim=c(npulses,length(energies),length(filterLengths)))
derivArrayMean4samplesPair <- array(data=NA, dim=c(nPairs, length(energies), nseps,length(filterLengths)))

#=====================================================================
#                                                                    #
#                              SINGLE PULSES                         #
#                                                                    #
#=====================================================================


# Get reconstructed energies AND <4-samples> derivative of SINGLE pulses
#=======================================================================
EkeVrecons <- array(data=NA,dim=c(npulses,length(energies),length(filterLengths)))
#BiasEkeVrecons <- array(data=NA,dim=c(npulses,length(energies),length(filterLengths)))
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
        #BiasEkeVrecons[,ie,ifl] <- EkeVrecons[,ie,ifl] - as.numeric(energies[ie])
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
EkeVreconsPairs <- array(NA, dim=c(nPairs, length(energies), nseps,length(filterLengths)))
#BiasEkeVreconsPairs <- array(NA, dim=c(nPairs, length(energies), nseps,length(filterLengths)))
# Get reconstructed energies of PAIRS of pulses
for (ie in 1:length(energies)){
    cat("Working with PAIRS: EnergyPrim=",energies[ie],"EnergySec=",Esec,"\n")
    for (is in 1:nseps){
        for (ifl in 1:length(filterLengths)){
            fl <- filterLengths[ifl]
            if(separationsBGplots[is] == "00002"){ # too short filter: cannot calculate energy
                EkeVreconsPairs[,ie,is,ifl] <- rep(NaN,nPairs)
                #BiasEkeVreconsPairs[,ie,is,ifl] <- rep(NaN,nPairs)
                next
            }
            # get reconstructed energies
            eventsFile <- paste("eresolLPA2shunt/nodetSP/events_sep",separationsBGplots[is],"sam_2000p_SIRENA",pulseLength,"_pL",
                                pulseLength,"_",energies[ie],"keV_",Esec,"keV_NTRIG_F0F_fixedlib",fEnergy,"OF_OPTFILT",fl,
                                samprateStr,".fits",sep="")
            #cat("  Reading file:", eventsFile,"\n")
            zz <- file(description = eventsFile, open = "rb")
            header0 <- readFITSheader(zz, fixHdr = 'none') # read primary header
            header <- readFITSheader(zz, fixHdr = 'none') # read extension header
            evtTable <- readFITSbintable(zz, header)
            close(zz)
            idcol <- which(evtTable$colNames == "SIGNAL")
            EkeVreconsPairs[,ie,is,ifl] <- evtTable$col[[idcol]][1:nPairs]
            #BiasEkeVreconsPairs[,ie,is,ifl] <- EkeVreconsPairs[,ie,is,ifl] - as.numeric(energies[ie])
            # get derivative
            idcol <- which(evtTable$colNames == "AVG4SD")
            derivArrayMean4samplesPair[,ie,is,ifl] <- evtTable$col[[idcol]][1:nPairs]
            
        } # foreach filter length    
    } # foreach separation    
} # foreach calib energy
cat("PAIR pulses (yellow bagplot) finished\n")

par(mfrow=c(2,3))
#
# DEFINE array of separations to Plot for Eprim, Filter
separationsToPlot <- array(data="", dim=c(nenergies,length(filterLengths)), 
                           dimnames = list(energies,filterLengths))
#
# DEFINE array of crosses with bandMax (rootsMax) and with bandMin (rootsMin) for Eprim, Filter
rootsMax <- array(data="", dim=c(nenergies,length(filterLengths)), 
                  dimnames = list(energies,filterLengths))
rootsMin <- array(data="", dim=c(nenergies,length(filterLengths)), 
                  dimnames = list(energies,filterLengths))
#
# PLOTTING
#
pdf(paste("baselineLPA2/derivativePairsStudy_Esec",Esec,"keV",samprateStr,".pdf",sep=""),width=10, height=7,version="1.4")

indexEsec <- which(energies == Esec)
lcols <- rainbow(length(energies))
cat("antes del bucle\n")
for (ie in 1:length(energies)){
    #
    # Plot variation of Reconstructed energy with separation for every filter length
    #
    bagplots.separations <- list()
    
    for (ifl in 1:length(filterLengths)){
        fLength=filterLengths[ifl]
        cat("Looking for roots for EnergyPrim=",energies[ie],"EnergySec=",Esec,"Filter=",fLength,"\n")
        nrootsMin <- 0
        nrootsMax <- 0
        meanErecon <- numeric(nseps)
        maxErecon <- numeric(nseps)
        minErecon <- numeric(nseps)
        # define orange band
        bandMin <- min(EkeVrecons[,ie,ifl])
        bandMax <- max(EkeVrecons[,ie,ifl])
        
        # look for separations to plot
        for(is in 1:nseps){
            meanErecon[is] <- mean(EkeVreconsPairs[,ie,is,ifl])
            maxErecon[is] <- max(EkeVreconsPairs[,ie,is,ifl])
            minErecon[is] <- min(EkeVreconsPairs[,ie,is,ifl])
            if((meanErecon[is]<=bandMax && meanErecon[is]>=bandMin ) || 
               (maxErecon[is] <=bandMax && maxErecon[is] >=bandMin ) ||
               (minErecon[is] <=bandMax && minErecon[is] >=bandMin ) ||
               (minErecon[is] <=bandMin && maxErecon[is] >=bandMax )){
                if(as.numeric(separationsBGplots[is]) > fLength) next
                if (separationsToPlot[ie,ifl] == ""){
                    separationsToPlot[ie,ifl] <- separationsBGplots[is]
                }else{
                    separationsToPlot[ie,ifl] <- paste(separationsToPlot[ie,ifl],",",separationsBGplots[is],sep="")
                }
                #cat("Adding",separationsBGplots[is]," ")
            }
        }
        roots<-numeric(2)
        # look for Crosses of bandMax
        is <- 3
        cat("Max crosses\n")
        while(is <= nseps && as.numeric(separationsBGplots[is])<fLength){
            cat("is=",is)
            crossSignMax <- (meanErecon[is]-bandMax)*(meanErecon[is-2]-bandMax)
            if(crossSignMax<0){
                x1 <- as.numeric(separationsBGplots[is-2])
                x2 <- as.numeric(separationsBGplots[is-1])
                x3 <- as.numeric(separationsBGplots[is])
                y1 <- meanErecon[is-2]
                y2 <- meanErecon[is-1]
                y3 <- meanErecon[is]
                # get cross with bandMax
                roots <- findCrosses(x1,x2,x3,y1,y2,y3,bandMax)
                okroot <- roots[whichClosest(roots,separationsBGplots[is])]
                # cat("okroot=",okroot," for is=",is,"")
                # get prev separation if minerr bar below bandMax
                altroot <- 0
                if(meanErecon[is-1] > meanErecon[is] && minErecon[is-1] < bandMax){ #en bajada
                    cat("en bajada\n")
                    altroot <- as.numeric(separationsBGplots[is-1])
                    isnew <- is-1
                    while (isnew > 1){
                        if(meanErecon[isnew] > meanErecon[is] && minErecon[isnew] < bandMax){
                            altroot <- as.numeric(separationsBGplots[isnew])
                            isnew <- isnew-1
                        }else{
                            break
                        }
                    }
                }else if(meanErecon[is] > meanErecon[is-1] && minErecon[is] < bandMax){ #en subida
                    #cat("en subida\n")
                    altroot <- as.numeric(separationsBGplots[is])
                    isnew <- is+1
                    while (isnew > 1){
                        if(meanErecon[isnew] > meanErecon[is] && minErecon[isnew] < bandMax){
                            altroot <- as.numeric(separationsBGplots[isnew])
                            isnew <- isnew+1
                        }else{
                            break
                        }
                    }#while isnew
                    is <- isnew # start again after new root
                } # if error bars in band
                newroot <- okroot
                if(altroot > 0) newroot<-altroot
                
                #cat("find root ",newroot," for is=", is,"\n")
                if(rootsMax[ie,ifl] == ""){
                    rootsMax[ie,ifl] <- newroot
                }else{
                    rootsMax[ie,ifl] <- paste(rootsMax[ie,ifl],",",newroot,sep="")
                }
                nrootsMax <- nrootsMax + 1
                is <- is + 1
            } # if cross
            is <- is+1
        }#while is
        
        # look for Crosses of bandMin
        roots<-numeric(2)
        is <- 3
        cat("Min crosses\n")
        while(is <= nseps && as.numeric(separationsBGplots[is])<fLength){
            crossSignMin <- (meanErecon[is]-bandMin)*(meanErecon[is-2]-bandMin)
            if(crossSignMin<0){
                cat("is in min=",is)
                x1 <- as.numeric(separationsBGplots[is-2])
                x2 <- as.numeric(separationsBGplots[is-1])
                x3 <- as.numeric(separationsBGplots[is])
                y1 <- meanErecon[is-2]
                y2 <- meanErecon[is-1]
                y3 <- meanErecon[is]
                # get cross with bandMin
                roots <- findCrosses(x1,x2,x3,y1,y2,y3,bandMin)
                okroot <- roots[whichClosest(roots,separationsBGplots[is])]
                # get prev separation if minerr bar below bandMax
                altroot <- 0
                if(meanErecon[is-1] > meanErecon[is] && maxErecon[is] > bandMin){#en bajada
                    cat("en bajada\n")
                    altroot <- as.numeric(separationsBGplots[is])
                    isnew <- is+1
                    while (isnew<nseps){
                        isnew <- isnew+1
                        if(meanErecon[is] > meanErecon[isnew] && maxErecon[isnew] > bandMin){
                            altroot <- as.numeric(separationsBGplots[isnew])
                            isnew <-isnew+1
                        }else{
                            break
                        }
                    }#while isnew 
                    is <- isnew # start again after new root
                }else if(meanErecon[is] > meanErecon[is-1] && maxErecon[is-1] > bandMin){#en subida
                    cat("en subida\n")
                    altroot <- as.numeric(separationsBGplots[is-1])
                    isnew <-is-1
                    while (isnew<nseps){
                        isnew <- isnew-1
                        if(meanErecon[is] > meanErecon[isnew] && maxErecon[isnew] > bandMin){
                            altroot <- as.numeric(separationsBGplots[isnew])
                        }else{
                            break
                        }
                    } #while isnew        
                }#if errbars in band
                newroot <- okroot
                if(altroot > 0) newroot<-altroot
                
                if(rootsMin[ie,ifl] == ""){
                    rootsMin[ie,ifl] <- newroot
                }else{
                    rootsMin[ie,ifl] <- paste(rootsMin[ie,ifl],",",newroot,sep="")
                }
                nrootsMin <- nrootsMin + 1
                is <- is + 1
            } # if cross
            is <- is+1
        } #while is
        
        stopifnot(abs(nrootsMax-nrootsMin)<2)
        if(nrootsMax>nrootsMin){
            rootsMin[ie,ifl] <- paste(rootsMin[ie,ifl],",",fLength,sep="")
        }else if(nrootsMax<nrootsMin){
            rootsMax[ie,ifl] <- paste(rootsMax[ie,ifl],",",fLength,sep="")
        }
        
        # save interesting quantities for detectionMaps y frequencyPairs
        # fileBPseps <- paste("baselineLPA2/BPseps_",Esec,"keV",samprateStr,".dat",sep="")
        # save(separationsBGplots,file=fileBPseps)
        
        fileBPfail <- paste("baselineLPA2/BPfail_",Esec,"keV",samprateStr,".dat",sep="")
        #save(separationsToPlot,file=fileBPfail)
        save(rootsMax,rootsMin,file=fileBPfail)
        
        bgpseps <-unlist(strsplit(separationsToPlot[ie,ifl],","))
        bgpseps <- bgpseps[bgpseps != ""]
        bagplots.separations[[ifl]] <-bgpseps
        
        #
        # Draw reconstruction curves
        #
        errbar(x=as.numeric(separationsBGplots),y=meanErecon, yplus=maxErecon, yminus=minErecon,
                col="darkmagenta", pch=1, typ="b",cex=0.7, errbar.col="cornflowerblue",
                xlab="Pair separation (samples)", ylab="<Erecon>", 
                xlim=c(min(as.numeric(separationsBGplots)),xmax[ifl]))
        title(main=(bquote(paste(E[prim],"=",.(energies[ie])," keV  ",
                                 E[sec],"=",.(Esec)," keV  ", "FLength=",.(fLength),
                                 "  FEnergy=",.(fEnergy)," keV",sep=""))),cex.main=0.8)
        abline(h=bandMax, col="orange",lty=2)
        abline(h=bandMin, col="orange",lty=2)
        
        # draw vertical lines in band crosses
        rootsMaxArr <-unlist(strsplit(rootsMax[ie,ifl],","))
        rootsMaxArr <- as.numeric(rootsMaxArr[rootsMaxArr != ""])
        rootsMinArr <-unlist(strsplit(rootsMin[ie,ifl],","))
        rootsMinArr <- as.numeric(rootsMinArr[rootsMinArr != ""])
        
        for (ir in 1:length(rootsMaxArr)){
            abline(v=rootsMaxArr[ir], col="gray", lty=2)
            abline(v=rootsMinArr[ir], col="gray", lty=2)
            #text(rootsMax[ir],ir*(max(maxErecon)-bandMax)/10+bandMax, cex=0.6,
            #     paste("(",sprintf("%.1f",sortedCrosses[2*ir-1]),",",
            #           sprintf("%.1f",sortedCrosses[2*ir]),")",sep=""))
            if(ir%%2 != 0){
                text(rootsMaxArr[ir],ir*(max(maxErecon)-bandMax)/10+bandMax, cex=0.6,
                     paste("(",sprintf("%.1f",rootsMaxArr[ir]),",",
                                sprintf("%.1f",rootsMinArr[ir]),")",sep=""))
            }else{
                text(rootsMaxArr[ir],ir*(max(maxErecon)-bandMax)/10+bandMax, cex=0.6,
                     paste("(",sprintf("%.1f",rootsMinArr[ir]),",",
                           sprintf("%.1f",rootsMaxArr[ir]),")",sep=""))
            }
        }
        abline(v=fLength, col="black",lty=2)
        
    } # foreach ifl (curves)
    #
    # Draw Bagplots 
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
                     EkeVreconsPairs[,ie,which(separationsBGplots %in% bagplots.separations[[ifl]]),ifl])
         ymax2 <- max(EkeVrecons[,ie,ifl],
                     EkeVreconsPairs[,ie,which(separationsBGplots %in% bagplots.separations[[ifl]]),ifl])
    
         # Yellow-ish bagplot (single pulses)
         bagplot(derivArrayMean4samples[,ie,ifl],EkeVrecons[,ie,ifl],xlab="<4 derivative samples>",
                 ylab="Reconstructed Energy (keV)",
                 main=(bquote(paste(E[prim],"=",.(energies[ie])," keV  ",
                                    E[sec],"=",.(Esec)," keV  ", "FLength=",.(fLength),
                                    "  FEnergy=",.(fEnergy)," keV",sep=""))), cex.main=0.8,
                 xlim=c(floor(xmin2),ceiling(xmax2*1.005)), ylim=c(ymin2,ymax2),
                 col.loophull="cornsilk",col.looppoints="peachpuff",col.baghull="orange",show.outlier = FALSE)
                 # blue-ish bagplots (pairs of pulses)
         
         text(xmax2,ymax2,"Separation(sam)", cex=0.8 )
         for(is in 1:nseps){
             if(!separationsBGplots[is] %in% bagplots.separations[[ifl]]) next
             if(as.numeric(separationsBGplots[is]) > fLength) next
             cat("............Plotting bagplot for sep=",separationsBGplots[is],"samples\n")
             bag <- try(compute.bagplot(derivArrayMean4samplesPair[,ie,is,ifl],EkeVreconsPairs[,ie,is,ifl]))
    
             if(class(bag) == "try-error") next
             plot.bagplot(bag,add=TRUE, transparency=TRUE)
    
             #text(xmax*1.001,bag$center[2],paste(separationsBGplots[is],"sam/",seps.ms[is],"ms",sep=""), cex=0.7)
             text(xmax2*1.001,bag$center[2],separationsBGplots[is], cex=0.7)
         } #separationsBGplots
         abline(h=max(EkeVrecons[,ie,ifl]), col="orange",lty=2)
         abline(h=min(EkeVrecons[,ie,ifl]), col="orange",lty=2)

     } # filter lengths
}
dev.off()

