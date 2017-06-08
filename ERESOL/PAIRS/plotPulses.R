#Plot pulses shape (when close)

require(FITSio)
setwd("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/")
pulseLength<- 4096
derivLength <- pulseLength - 1 
plotLength <- 800
baseline <- 706
sep <- 136
sep <- 5
sep <- 20
sepM1 <- sep - 1
sepP1 <- sep + 1 
sepstr <- sprintf("%05d",sep)
zeros <- rep(0,sep)
libEnergies <- c(0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8)
EkeV1 <- 6
EkeV2 <- 0.2
idx1 <- which(libEnergies == EkeV1)
idx2 <- which(libEnergies == EkeV2)
initsamp <- 10
#EkeV <- "02"

pdf(paste("PairsOfPulsesShapes_",EkeV1,"keV_",EkeV2,"keV_sep",sep,".pdf",sep=""),width=7, height=7)

if (0){
# Read additional data (always present) of 0.2 keV pulses for comparison
Data02keV<-read.table("../../testHarness/simulations/SIXTE/LIBRARIES/tessimLPA2shunt/GLOBAL/ADC/pulse02keV_noBSL.txt",
                      header=F)
pulse02keV <- Data02keV[,1]
Data02keV_pair1<-read.table("tessimLPA2shunt/sep00001_02keV_pulsesRow1.txt", header=F)
Data02keV_pair10<-read.table("tessimLPA2shunt/sep00010_02keV_pulsesRow1.txt", header=F)
Data02keV_pair56<-read.table("tessimLPA2shunt/sep00056_02keV_pulsesRow1.txt", header=F)
Data02keV_pair101<-read.table("tessimLPA2shunt/sep00101_02keV_pulsesRow1.txt", header=F)

pulse02keV_pair1 <- Data02keV_pair1[-(1:999),1] - baseline
pulse02keV_pair10 <- Data02keV_pair10[-(1:999),1] - baseline
pulse02keV_pair56 <- Data02keV_pair56[-(1:999),1] - baseline
pulse02keV_pair101 <- Data02keV_pair101[-(1:999),1] - baseline

derivPulse02 <- diff(pulse02keV)
derivPair02_sep1 <- diff(pulse02keV_pair1)
derivPair02_sep10 <- diff(pulse02keV_pair10)
derivPair02_sep56 <- diff(pulse02keV_pair56)
derivPair02_sep101 <- diff(pulse02keV_pair101)
}
# Take single pulses from the library of 200.000 pulses
filenameFITS <- "../../testHarness/simulations/SIXTE/LIBRARIES/tessimLPA2shunt/GLOBAL/ADC/libraryMultiE_GLOBAL_PL4096_200000p.fits"
zz <- file(description = filenameFITS, open = "rb")
header0 <- readFITSheader(zz, fixHdr = 'remove') # read primary header
header <- readFITSheader(zz, fixHdr = 'remove') # read extension header
Dpulses <- readFITSbintable(zz, header)
close(zz)
idcol <- which(Dpulses$colNames == "PULSEB0")
DatakeV1 <- Dpulses$col[[idcol]][idx1,]
pulsekeV1 <- c(rep(0,(initsamp-1)),DatakeV1)
DatakeV2 <- Dpulses$col[[idcol]][idx2,]
pulsekeV2 <- c(rep(0,(initsamp-1)),DatakeV2)


# if FITS read does not work...
#filenameTXT <- paste("../../testHarness/simulations/SIXTE/LIBRARIES/tessimLPA2shunt/GLOBAL/ADC/pulse",EkeV1,"keV_noBSL.txt",sep="")
#DatakeV1<-read.table(filenameTXT,header=F) 
#pulsekeV1 <- c(rep(0,(initsamp-1)),DatakeV1[,1])
#DatakeV2<-read.table(paste("../../testHarness/simulations/SIXTE/LIBRARIES/tessimLPA2shunt/GLOBAL/ADC/pulse",EkeV2,"keV_noBSL.txt",sep=""),
#                     header=F) 
#pulsekeV2 <- c(rep(0,(initsamp-1)),DatakeV2[,1])


# Take pair data from PAIRS @ EkeV keV && sep=SEP
#DatakeV_pairSEP<-read.table(paste("tessimLPA2shunt/sep",sepstr,"_",EkeV1,"keV_",EkeV2,"keV_pulsesRow1.txt",sep=""), header=F)
#pulsekeV_pairSEP <- DatakeV_pairSEP[-(1:999),1] - baseline
#pulsekeV_pairSEP <- DatakeV_pairSEP[-(1:(1000-initsamp)),1] - baseline # add initial sample in baseline
filenamePAIR <- paste("tessimLPA2shunt/sep",sepstr,"sam_2000p_",EkeV1,"keV_",EkeV2,"keV.fits",sep="")
kk <- file(description = filenamePAIR, open = "rb")
header0 <- readFITSheader(kk, fixHdr = 'remove') # read primary header
header <- readFITSheader(kk, fixHdr = 'remove') # read extension header
Dpair <- readFITSbintable(kk, header)
close(kk)
idcol <- which(Dpair$colNames == "ADC")
DatakeV_pairSEP <- Dpair$col[[idcol]][1,]
istart <- 1000 # zero samples before pulse Start
if (EkeV1 == 0.2) istart = 999
pulsekeV_pairSEP <- DatakeV_pairSEP[-(1:(istart-initsamp)),1] - baseline

offset <- append(zeros,pulsekeV2)
DataSum <- pulsekeV1[1:pulseLength] + offset[1:pulseLength]
derivDataSum <-diff(DataSum)
derivPulse1 <- diff(pulsekeV1)
derivPulse2 <- diff(pulsekeV2)
derivPairs <- diff(pulsekeV_pairSEP)   
derivOffset<- diff(offset)
derivSumMinusDerivPulse1 <- derivDataSum[1:derivLength]-derivPulse1[1:derivLength]
derivSumMinusDerivPulse12 <- derivDataSum[1:derivLength]-derivPulse1[1:derivLength]-derivOffset[1:derivLength]
derivPairsMinusDerivPulse1 <- derivPairs[1:derivLength]-derivPulse1[1:derivLength]
derivOffsetScaled <- max(derivPairsMinusDerivPulse1[1:derivLength])/max(derivOffset[1:derivLength])*derivOffset[1:derivLength]
derivPairsMinusDerivPulse12 <- derivPairsMinusDerivPulse1[1:derivLength] - derivOffsetScaled[1:derivLength]

colSumm <- "blue"
colPairs <- "red"

#--------------------------------------------------
# PLOT 2 pulses: SUMM and PAIR simulation (tessim)
#--------------------------------------------------
plotLength <- 600
plot(seq(1,pulseLength),DataSum, type="l", col=colSumm, 
      xlim=c(1,plotLength),ylab = "ADC", xlab="sample", main=paste("Pulses ",EkeV1," & ",EkeV2, " keV",sep="") )
lines(seq(1,pulseLength),pulsekeV_pairSEP[1:pulseLength], type="l", col=colPairs)
lines(seq(1,pulseLength),pulsekeV1[1:pulseLength], lty=2,cex=0.1)
lines(seq(1,pulseLength),offset[1:pulseLength],lty=2,cex=0.1)


abline(v=45,col="grey")
abline(v=sep+1,col="grey")
abline(v=45+sep+1,col="grey")
legend("topright",c("Single Pulses (lib)","Summed pulses",paste("Simulated Pair (sep",sepstr,")",sep="")),
             col=c("black",colSumm,colPairs), lty=c(2,1,1),cex=0.8)


#-----------------------------------------------------
# PLOT 2 pulses DERIVATIVE: deriv_SUMM and deriv_PAIR 
#-----------------------------------------------------
plotLength <- 100
plot(seq(1,derivLength),derivDataSum[1:derivLength], type="l",lty=3, col=colSumm,
     xlim=c(1,plotLength),ylab = "Derivate", main=paste("Pulses ",EkeV1," & ",EkeV2," keV",sep=""), xlab="sample")
#points(seq(1,derivLength),derivDataSum[1:derivLength], pch=1, cex=0.3)
lines(seq(1,derivLength),derivPairs[1:derivLength], lty=2, col=colPairs)
#points(seq(1,derivLength),derivPairs[1:derivLength], pch=2, col=colPairs,cex=0.4)


legend("topright",
       c("Deriv of Summed pulses",
         paste("Deriv of Simulated Pair (sep",sepstr,")",sep="")),
       col=c(colSumm,colPairs), 
       lty=c(3,2),
       pch=c(NA,NA),cex=0.8, bty="n")



#-----------------------------------------------------
# PLOT 2 pulses DERIVATIVE - derivPulse1: deriv_SUMM and deriv_PAIR 
#-----------------------------------------------------
plotLength <- 300
plot(seq(1,derivLength),derivSumMinusDerivPulse1[1:derivLength], pch=4, col=colSumm, 
     xlim=c(1,plotLength),ylab = "Deriv. all - Deriv. Pulse1", main=paste("Pulses ",EkeV1," & ",EkeV2," keV",sep=""), xlab="sample", cex=0.4)
lines(seq(1,derivLength),derivPairsMinusDerivPulse1[1:derivLength],col=colPairs)
points(seq(1,derivLength),derivPairsMinusDerivPulse1[1:derivLength],pch=4, cex=0.4,col=colPairs)
lines(seq(sepP1,pulseLength),derivPulse2[1:(derivLength-sepM1)], lty=2, col="black")

legend("topright", 
       c("Deriv of Summed -Deriv of Single1",
         "Deriv of Pairs -Deriv of Single1",
         "Deriv of Single2"),
       col=c(colSumm,colPairs,"black"), pch=c(4,4,NA), lty=c(NA,NA,2),cex=0.8, bty="n")


#-----------------------------------------------------
# PLOT 2 pulses DERIVATIVE - derivPulse1 - derivPulse2
#-----------------------------------------------------
plotLength <- 300
plot(seq(1,derivLength),derivSumMinusDerivPulse12[1:derivLength], pch=2, col="blue", 
     xlim=c(0,plotLength),ylab ="Deriv. all - Deriv. Pulse12", main=paste("Pulses ",EkeV1," & ", EkeV2," keV",paste=""), 
     xlab="sample", cex=0.4, ylim=c(-50,200))
points(seq(1,derivLength),derivPairsMinusDerivPulse12[1:derivLength], pch=2, cex=0.4,col="red")

mycols <- colorRampPalette(c("green","magenta", "cyan"))(5)
legend("topright", 
       c("Deriv of Summed -Deriv of Single12",
         "Deriv of Pairs -Deriv of Single12(scaled)",
         "Deriv of Pulse02keV","Deriv of Pair02keV_sep1",
         "Deriv of Pair02keV_sep10","Deriv of Pair02keV_sep56",
         "Deriv of Pair02keV_sep101"),
       col=c("blue","red",mycols), pch=c(2,2,4,4,4,4,4),cex=0.8, bty="n")
if(0){
loc02 <- 150
loc02 <- 50
loc02_2 <- loc02 - 2
points(seq(loc02,pulseLength),derivPulse02[1:(derivLength-loc02_2)], pch=4, cex=0.4, col=mycols[1])
lines(seq(loc02,pulseLength),derivPair02_sep1[1:(derivLength-loc02_2)], type="b",pch=4, cex=0.4, col=mycols[2])
lines(seq(loc02,pulseLength),derivPair02_sep10[1:(derivLength-loc02_2)], type="b", pch=4, cex=0.4, col=mycols[3])
lines(seq(loc02,pulseLength),derivPair02_sep56[1:(derivLength-loc02_2)], type="b", pch=4, cex=0.4, col=mycols[4])
lines(seq(loc02,pulseLength),derivPair02_sep101[1:(derivLength-loc02_2)], type="b", pch=4, cex=0.4, col=mycols[5])
abline(h=35,col="grey")
abline(h=18,col="grey")
}
dev.off()







