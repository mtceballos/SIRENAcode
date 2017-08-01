#           SECONDARY DETECTION
# ShiftVSsecE:
# Plot EnergyReconsCorrected vs EnergyInputSecondary for different separations of Secondary pulse  (no Sec Detection)
# Plot EnergyReconsCorrected vs EnergyInputSecondary for different separations of Secondary pulse  (Sec Detection)
# ShiftVSsep:
# Plot EnergyReconsCorrected vs separation for different Primary Energies

PLOT <- "ShiftVSsep"

setwd("/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresolLPA2shunt")
pdf(paste("../baselineLPA2/Secondaries",PLOT,".pdf",sep=""))
sampPerSec <- 156250 # samprate
sepSingles <- "40000"
detectTHnoAD <- c()

library(rjson)
library(Hmisc)
nSamples <- "4096"
pulseLength <- "4096"
nSimPulses <- "2000"
nSimPulsesSingles <- "20000"
sepsStr <-c("00005", "00010", "00020", "00045", "00060", "00100", "00200", "00250","00300","00400", "00800")
#sepsStr <-c("00005", "00010", "00020", "00045", "00060", "00100", "00200","00400", "00800")
detectSPdirs <- c("nodetSP",".") # blank is for detection with automatic location of secondaries
detectSPdirs <- c("nodetSP") # blank is for detection with automatic location of secondaries
pulses       <- c("all", "primaries") # "all" if no sec detection // "primaries" if sec detection
#detectSPdirs <- c("nodetSP")
seps <- as.numeric(sepsStr)
colores <- rainbow(length(seps))

if(PLOT == "ShiftVSsecE"){
    primEnergy <- 6 # keV
    singlePrimEnergy <- 5.1914 # energy of 6 keV single pulses recons with OPTFILT 1 keV
    secEnergies <- c(0.2, 0.5, 1, 2, 3, 4, 5, 6, 7)
    energyCorrPrim    <- matrix(NA,nrow=length(seps), ncol=length(secEnergies))
    energyUnCorrPrim    <- matrix(NA,nrow=length(seps), ncol=length(secEnergies))
    ebiasCorr         <- matrix(NA,nrow=length(seps), ncol=length(secEnergies))
    ebiasUnCorr       <- matrix(NA,nrow=length(seps), ncol=length(secEnergies))
    diffPairSingle    <- matrix(NA,nrow=length(seps), ncol=length(secEnergies))
    fwhmPrimGAINCORRE <- matrix(NA,nrow=length(seps), ncol=length(secEnergies))
    fwhmPrimGAINCORREErr <- matrix(NA,nrow=length(seps), ncol=length(secEnergies))
    TRIGG <- "" # or "_NTRIGG"
    
    labels <- character(length(seps))
    
    # READ BIAS 
    # =========
    
    for (id in 1:length(detectSPdirs)){
        SPdir <- detectSPdirs[id]
        plot(secEnergies,seq(0.,0.8,length.out = length(secEnergies)), xlab="Secondaries Energies (keV)",
             ylab="(Eprim(pair)-Eprim(single))/EprimSingle", type="n", 
             main=paste("Reconstruction of 6 keV pulses \n Secondaries (",SPdir,")",sep=""))
        
        abline(v=secEnergies, lty=2, col="gray")
        
        for (is in 1:length(seps)){
            for (ie in 1:length(secEnergies)){
                eresolFile <- paste(SPdir,"/eresol_",nSimPulses,"p_SIRENA",nSamples,"_pL",pulseLength,
                                    "_",primEnergy,"keV_",secEnergies[ie],"keV_","F0F_", "fixedlib1OF_OPTFILT.json",sep="")
                if(file.exists(eresolFile)){
                    cat("Reading file ",eresolFile,"\n")
                    jsondata <- fromJSON(file=eresolFile)
                    idxSep <- which(sapply(jsondata,function(x) x$separation)==sepsStr[is])
                    if(length(idxSep)>0){
                        ebiasCorr[is,ie] <- as.numeric(jsondata[[idxSep]]$biasEreal$all) #eV
                        ebiasUnCorr[is,ie] <- as.numeric(jsondata[[idxSep]]$biasErecons$all) #eV
                        energyCorrPrim[is,ie] <- ebiasCorr[is,ie]/1000. + primEnergy
                        energyUnCorrPrim[is,ie] <- ebiasUnCorr[is,ie]/1000. + primEnergy
                        diffPairSingle[is,ie] <- (energyUnCorrPrim[is,ie]-singlePrimEnergy)/singlePrimEnergy
                        fwhmPrimGAINCORRE[is,ie] <- as.numeric(jsondata[[idxSep]]$fwhmEreal$all)
                        fwhmPrimGAINCORREErr[is,ie] <- as.numeric(jsondata[[idxSep]]$fwhmEreal_err$all)
                        
                    }else{
                        energyCorrPrim[is,ie] <- NaN
                        diffPairSingle[is,ie] <- NaN
                        fwhmPrimGAINCORRE[is,ie] <- NaN
                    }
                    
                }else{
                    warning("Not-existing file:", eresolFile)
                    energyCorrPrim[is,ie] <- NaN
                    fwhmPrimGAINCORRE[is,ie] <- NaN
                } # file exists
            }# secEnergies [ie]
            #lines(secEnergies, energyCorrPrim[is,], col=colores[is],lty=id, typ="b",cex=0.5)
            #
            lines(secEnergies, diffPairSingle[is,], col=colores[is],lty=id, typ="b",cex=0.5)
            #lines(secEnergies, fwhmPrimGAINCORRE[is,], col=colores[is],lty=id, typ="b",cex=0.5)
            #errbar(secEnergies,fwhmPrimGAINCORRE[is,], yplus=fwhmPrimGAINCORRE[is,]+fwhmPrimGAINCORREErr[is,],
            #       type="n",cap=0, yminus=fwhmPrimGAINCORRE[is,]-fwhmPrimGAINCORREErr[is,],add=TRUE,errbar.col=colores[is])
            #abline(h=2.05,)
            labels[is] <- paste("Separation ",seps[is], " samples", sep="")
        }# separations [is]
    }# detection of Secondaries [id]
    legend("bottomright",legend=labels, col=colores, cex=0.6, bty="n", lty=1)
    par(fig = c(0.07,0.5, 0.45, 0.95), new = T)
    plot(secEnergies,seq(0.,0.2,length.out = length(secEnergies)), type="n", xlab="", ylab="",cex.axis=0.5,
         xlim=c(0.2,1), ylim=c(0,0.15))
    for (is in 1:length(seps)){
        lines(secEnergies, diffPairSingle[is,], col=colores[is],lty=id, typ="b",cex=0.2)
    }
}else if (PLOT=="ShiftVSsep"){
    # Plot Shift vs separation for Eprim + 0.2
    secEnergy <- 0.2 # keV
    seps <- as.numeric(sepsStr)
    primEnergies <- c(0.2, 0.5, 1, 2, 3, 4, 5, 6, 7,8)
    
    energyCorrPrim    <- matrix(NA,nrow=length(seps), ncol=length(primEnergies))
    energyUnCorrPrim  <- matrix(NA,nrow=length(seps), ncol=length(primEnergies))
    ebiasCorr         <- matrix(NA,nrow=length(seps), ncol=length(primEnergies))
    ebiasUnCorr       <- matrix(NA,nrow=length(seps), ncol=length(primEnergies))
    diffPairSingle    <- matrix(NA,nrow=length(seps), ncol=length(primEnergies))
    fwhmPrimGAINCORRE <- matrix(NA,nrow=length(seps), ncol=length(primEnergies))
    fwhmPrimGAINCORREErr <- matrix(NA,nrow=length(seps), ncol=length(primEnergies))
    ebiasCorrSingles  <- numeric(length(primEnergies))
    ebiasUnCorrSingles <- numeric(length(primEnergies))
    fwhmPrimGAINCORRESingles <- numeric(length(primEnergies))
    fwhmPrimGAINCORREErrSingles <- numeric(length(primEnergies))
    
    TRIGG <- "_NTRIG"
    
    # READ BIAS 
    # =========
    #drawLogPlotBox(xlimits=c(3E-2,8.),ylimits=c(-0.1,1), logxy="x",
    #xlabel = "Separations(ms)", ylabel=expression(("<"~E[recons]~">"~-E[prim])/E[prim]),
    #naxes=c(T,T,F,F))
    split.screen(rbind(c(0.01,0.95,0.3, 0.98), c(0.01, 0.95, 0.01, 0.5)))
    screen(1)
    #par(mar=c(5,5,4,2)) # give some margin to Y axis label
    par(mgp=c(2,1,0))
    drawLogPlotBox(xlimits=c(3E-2,8.),ylimits=c(-10,200), logxy="x",xlabel="", 
                   ylabel=expression("(<"~E[recons]~">"~-E[prim]~") eV"), naxes=c(T,T,F,F))
    title(main=paste("Reconstruction of Pairs (Eprim + 0.2 keV)\n (gain-scale corrected)",sep=""))
    rect(xleft=1E3*200/sampPerSec, ybottom = -20, xright = 1E3*300/sampPerSec,ytop = 500,
         density = NULL, col="gray90",border = NA)
    grid(nx=NA,ny=NULL,col="gray70")
    
    pchs <- numeric()
    labels <- character()
    ltys <- numeric()
    for (id in 1:length(detectSPdirs)){
        SPdir <- detectSPdirs[id]
        for (ie in 1:length(primEnergies)){
            eresolFile <- paste(SPdir,"/eresol_",nSimPulses,"p_SIRENA",nSamples,"_pL",pulseLength,
                                "_",primEnergies[ie],"keV_",secEnergy,"keV_","F0F_", 
                                "fixedlib1OF_OPTFILT",TRIGG,".json",sep="")
            if(file.exists(eresolFile)){
                cat("Reading file ",eresolFile,"\n")
                jsondata <- fromJSON(file=eresolFile)
            }else{
                warning("Not-existing file:", eresolFile)
                energyCorrPrim[,ie] <- NaN
                energyUnCorrPrim[,ie] <- NaN
                fwhmPrimGAINCORRE[,ie] <- NaN
                next #next energy
            } # file exists
            for (is in 1:length(seps)){
                idxSep <- which(sapply(jsondata,function(x) x$separation)==sepsStr[is])
                if(length(idxSep)>0){
                    ebiasCorr[is,ie] <- as.numeric(jsondata[[idxSep]]$biasEreal[[pulses[id]]]) #eV
                    ebiasUnCorr[is,ie] <- as.numeric(jsondata[[idxSep]]$biasErecons[[pulses[id]]]) #eV
                    energyCorrPrim[is,ie] <- ebiasCorr[is,ie]/1000. + primEnergies[ie]
                    energyUnCorrPrim[is,ie] <- ebiasUnCorr[is,ie]/1000. + primEnergies[ie]
                    fwhmPrimGAINCORRE[is,ie] <- as.numeric(jsondata[[idxSep]]$fwhmEreal[[pulses[id]]])
                    fwhmPrimGAINCORREErr[is,ie] <- as.numeric(jsondata[[idxSep]]$fwhmEreal_err[[pulses[id]]])
                }else{
                    cat("Non-existing separation(",sepsStr[is]," in ",eresolFile,"\n")
                    energyCorrPrim[is,ie] <- NaN
                    energyUnCorrPrim[is,ie] <- NaN
                    fwhmPrimGAINCORRE[is,ie] <- NaN
                }
            }# separations [is]
            bty <-"b"
            pt <- id-1
            if(id==1){
                bty="b"
                pt <-3
            }
            #lines(seps/sampPerSec*1E3, ebiasCorr[,ie]/(primEnergies[ie]*1000.), col=colores[ie],lty=id, typ=bty,cex=0.5, pch=pt)
            lines(seps/sampPerSec*1E3, ebiasCorr[,ie], col=colores[ie],lty=id, typ=bty,cex=0.5, pch=pt)
            pchs <- append(pchs,pt)
            ltys <- append(ltys,id)
            lab <- paste("E_prim=",primEnergies[ie], " keV", sep="")
            if(id == 2) lab <- paste("PrimEner",primEnergies[ie], " keV - no Sec detect", sep="")
            labels <-append(labels, lab)
            colores <- append(colores, colores[ie])
            
            #Plot bias for singles
            eresolFile <- paste("./eresol_",nSimPulsesSingles,"p_SIRENA",nSamples,"_pL",pulseLength,
                                "_",primEnergies[ie],"keV_","F0F_", 
                                "fixedlib1OF_OPTFILT",TRIGG,".json",sep="")
            if(file.exists(eresolFile)){
                cat("Reading file ",eresolFile,"\n")
                jsondata <- fromJSON(file=eresolFile)
            }else{
                warning("Not-existing file:", eresolFile)
                ebiasCorrSingles[ie] <- NA
                ebiasUnCorrSingles[ie] <- NA
                fwhmPrimGAINCORRESingles[ie] <- NA
                fwhmPrimGAINCORREErrSingles[ie] <- NA
                next #next energy
            } # file exists
            idxSep <- which(sapply(jsondata,function(x) x$separation)==sepSingles)
            ebiasCorrSingles[ie] <- as.numeric(jsondata[[idxSep]]$biasEreal[["all"]]) #eV
            ebiasUnCorrSingles[ie] <- as.numeric(jsondata[[idxSep]]$biasErecons[["all"]]) #eV
            fwhmPrimGAINCORRESingles[ie] <- as.numeric(jsondata[[idxSep]]$fwhmEreal[["all"]])
            fwhmPrimGAINCORREErrSingles[ie] <- as.numeric(jsondata[[idxSep]]$fwhmEreal_err[["all"]])
            biasLevel <-ebiasCorrSingles[ie]/(primEnergies[ie]*1000.)
            #abline(h=biasLevel,col="gray67",lty=2)
        }# primEnergies [ie]
    }# detection of Secondaries [id]
    legend("topright",legend=labels,pch=pchs,lty=ltys, col=colores, cex=0.7,bty="n")
    #text(0.05,-0.05,"OF reconstruction of singles", cex=0.7, col="grey68")
    screen(2)
    drawLogPlotBox(xlimits=c(3E-2,8.),ylimits=c(-15,15), logxy="x",xlabel = "Separations(ms)", 
                   ylabel=expression("(<"~E[recons]~">"~-E[prim]~") eV"),naxes=c(T,T,F,F))
    rect(xleft=1E3*200/sampPerSec, ybottom = -20, xright = 1E3*300/sampPerSec,ytop = 500,
         density = NULL, col="gray90",border = NA)
    for (ie in 1:length(primEnergies)){
        lines(seps/sampPerSec*1E3, ebiasCorr[,ie], col=colores[ie],lty=id, typ=bty,cex=0.5, pch=pt)
        abline(h=abs(ebiasCorrSingles[ie]),col=colores[ie],lty=2)
    }
    grid(nx=NA,ny=NULL,col="gray70")
    text(0.1,2,"Single pulses @E_prim", cex=0.8, col="grey67")
    #par(mgp=c(2,1,0),mar=c(5,4,3,2)+0.1)
    #plot(c(1:10),seq(0.01,0.1,length.out=10),xlim=c(0,1),ylim=c(1E-2,0.1),xaxt="n",
    #     type = "n",log="y",xlab="Separation=2.56ms",ylab=expression("(<"~E[recons]~">"~-E[prim]~") eV"))
    #title(main="Reconstruction of Single Pulses")
    #abline(h=abs(ebiasCorrSingles),col=colores,lty=2)
    
    close.screen(all.screens = TRUE)
}


    
dev.off()
        
