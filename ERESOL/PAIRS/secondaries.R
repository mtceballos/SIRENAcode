#           SECONDARY DETECTION
# 
# Plot EnergyReconsCorrected vs EnergyInputSecondary for different separations of Secondary pulse  (no Sec Detection)
# Plot EnergyReconsCorrected vs EnergyInputSecondary for different separations of Secondary pulse  (Sec Detection)

setwd("/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresolLPA2shunt")
library(rjson)
library(Hmisc)
nSamples <- "4096"
pulseLength <- "4096"
nSimPulses <- "2000"
primEnergy <- 6 # keV
sepsStr <- c("00005","00010", "00020", "00045")
seps <- as.numeric(sepsStr)
secEnergies <- c(0.2, 0.5, 1, 2, 3, 4, 5, 6, 7)
detectSPdirs <- c("nodetSP","detSP")
detectSPdirs <- c("nodetSP")
    
energyCorrPrim    <- matrix(NA,nrow=length(seps), ncol=length(secEnergies))
ebiasCorr         <- matrix(NA,nrow=length(seps), ncol=length(secEnergies))
accuracy          <- matrix(NA,nrow=length(seps), ncol=length(secEnergies))
fwhmPrimGAINCORRE <- matrix(NA,nrow=length(seps), ncol=length(secEnergies))
fwhmPrimGAINCORREErr <- matrix(NA,nrow=length(seps), ncol=length(secEnergies))
TRIGG <- "" # or "_NTRIGG"

colores <- rainbow(length(seps))
labels <- character(length(seps))

# READ BIAS 
# =========
#plot(secEnergies,seq(6.0,6.5,length.out = length(secEnergies)), xlab="Secondaries Energies (keV)",
#                     ylab="(Mean) Primary Energies (keV; gain Corr)", type="n", main="Reconstruction of 6 keV pulses")
plot(secEnergies,seq(0.1,3000.,length.out = length(secEnergies)), xlab="Secondaries Energies (keV)",
     ylab="Bias for Primary (eV; gain Corr)", type="n", main="Reconstruction of 6 keV pulses")
#plot(secEnergies,seq(1.8,2.5,length.out = length(secEnergies)), xlab="Secondaries Energies (keV)",
#     ylab="FWHM Primary (eV; gain Corr)", type="n", main="Reconstruction of 6 keV pulses")

abline(v=secEnergies, lty=2, col="gray")

for (id in 1:length(detectSPdirs)){
    SPdir <- detectSPdirs[id]
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
                    energyCorrPrim[is,ie] <- ebiasCorr[is,ie]/1000. + primEnergy
                    accuracy[is,ie] <- (ebiasCorr[is,ie]/1000. )/primEnergy
                    fwhmPrimGAINCORRE[is,ie] <- as.numeric(jsondata[[idxSep]]$fwhmEreal$all)
                    fwhmPrimGAINCORREErr[is,ie] <- as.numeric(jsondata[[idxSep]]$fwhmEreal_err$all)
                    
                }else{
                    energyCorrPrim[is,ie] <- NaN
                    accuracy[is,ie] <- NaN
                    fwhmPrimGAINCORRE[is,ie] <- NaN
                }
                
            }else{
                warning("Not-existing file:", eresolFile)
                energyCorrPrim[is,ie] <- NaN
                fwhmPrimGAINCORRE[is,ie] <- NaN
            } # file exists
        }# secEnergies [ie]
        #lines(secEnergies, energyCorrPrim[is,], col=colores[is],lty=id, typ="b",cex=0.5)
        lines(secEnergies, ebiasCorr[is,], col=colores[is],lty=id, typ="b",cex=0.5)
        abline(h=0.014)
        #lines(secEnergies, fwhmPrimGAINCORRE[is,], col=colores[is],lty=id, typ="b",cex=0.5)
        #errbar(secEnergies,fwhmPrimGAINCORRE[is,], yplus=fwhmPrimGAINCORRE[is,]+fwhmPrimGAINCORREErr[is,],
        #       type="n",cap=0, yminus=fwhmPrimGAINCORRE[is,]-fwhmPrimGAINCORREErr[is,],add=TRUE,errbar.col=colores[is])
        #abline(h=2.05,)
        labels[is] <- paste("Separation ",seps[is], " samples", sep="")
    }# separations [is]
}# detection of Secondaries [id]
legend("topleft",legend=labels, col=colores, cex=0.6, bty="n", lty=1)


#################################



secEnergies <- c(0.2, 3, 7)
primEnergies <- c(0.2, 3, 6)
colsPrim <- rainbow(length(primEnergies))

MinimDist <- matrix(data=c(10,60,60, 10,45,45, 10,45,45), 
                    nrow=length(secEnergies), ncol=length(primEnergies), byrow = T)
dimnames(MinimDist) <- list(c("sec02", "sec3", "sec7"), c("prim02", "prim3", "prim6"))

plot(secEnergies, seq(5,60, length.out=length(secEnergies)),xlab="Energy of Secondaries (keV)",
     ylab="Closest pulses for detection (samples)", type="n")
for (ie in 1:length(primEnergies)){
    lines(secEnergies, MinimDist[,ie],col=colsPrim[ie], lty=ie, type="b")
    text(primEnergies[ie],MinimDist[,ie]+5,paste("Eprim=",primEnergies[ie]," keV",sep=""), col=colsPrim[ie], cex=0.8)
}


