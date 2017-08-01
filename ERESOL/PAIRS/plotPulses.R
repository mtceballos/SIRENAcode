#Plot pulses shape (when close)

require(FITSio)
library(JBTools)
setwd("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/")
pulseLength<- 4096
derivLength <- pulseLength - 1 
plotLength <- 800
baseline <- 706
sigma11 <- 60
sigma6 <- 32
initsamp <- 10
libEnergies <- c("0.2", "0.5", "1", "2", "3", "4", "5", "6", "7", "8")
Der1     <- c(36.7518, 91.3329, 182.556, 364.723, 546.437, 727.802, 908.748, 1089.3, 1269.46, 1449.22)    

for (sep in c(5, 10, 20, 45)){
    sepM1 <- sep - 1
    sepP1 <- sep + 1 
    sepstr <- sprintf("%05d",sep)
    zeros <- rep(0,sep)

    for (EkeV1 in c("0.2", "3", "6")){
        for(EkeV2 in c("0.2","3", "7")){
            
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
            pulse <- list()
            sam1Der <- list()
            for (ie in 1:length(libEnergies)){
                DatakeV <- Dpulses$col[[idcol]][ie,]
                pulse[[libEnergies[ie]]] <-c(rep(0,(initsamp-1)),DatakeV)
                sam1Der[[libEnergies[ie]]] <- Der1[ie]
            }
            
            # if FITS read does not work...
            #filenameTXT <- paste("../../testHarness/simulations/SIXTE/LIBRARIES/tessimLPA2shunt/GLOBAL/ADC/pulse",EkeV1,"keV_noBSL.txt",sep="")
            #DatakeV1<-read.table(filenameTXT,header=F) 
            #pulsekeV1 <- c(rep(0,(initsamp-1)),DatakeV1[,1])
            #DatakeV2<-read.table(paste("../../testHarness/simulations/SIXTE/LIBRARIES/tessimLPA2shunt/GLOBAL/ADC/pulse",EkeV2,"keV_noBSL.txt",sep=""),
            #                     header=F) 
            #pulsekeV2 <- c(rep(0,(initsamp-1)),DatakeV2[,1])
            
            
            # Take pair data from PAIRS @ EkeV keV && sep=SEP
            filenamePAIR <- paste("tessimLPA2shunt/sep",sepstr,"_",EkeV1,"keV_",EkeV2,"keV_pulsesRow1.txt",sep="")
            if(!file.exists(filenamePAIR)){
                cat("Creating TXT file with pairs information")
                headas <- "export HEADAS=/home/ceballos/sw/heasoft/x86_64-unknown-linux-gnu-libc2.24"
                src <- "source /home/ceballos/sw/heasoft/x86_64-unknown-linux-gnu-libc2.24/headas-init.sh"
                # create file version without var len arrays
                fitsFile <- paste("tessimLPA2shunt/sep",sepstr,"sam_2000p_",EkeV1,"keV_",EkeV2,"keV.fits",sep="")
                if(!file.exists(fitsFile)){
                    # obtain file from larger one
                    fitsFileLarge <- paste("tessimLPA2shunt/sep",sepstr,"sam_20000p_",EkeV1,"keV_",EkeV2,"keV.fits",sep="")
                    stopifnot((file.exists(fitsFileLarge) || file.exists(fitsFile)))
                    fitsFile <- fitsFileLarge
                }
                toStilt <- paste("/home/ceballos/sw/topcat/stilts tpipe cmd='addcol TIME2 \"(TIME*2)\"' cmd='delcols TIME2' in=", 
                                 fitsFile, " out=pp_novar.fits", sep="")
                system(toStilt)  
                #comm <- paste(headas," && ", src," && cphead ", fitsFile, "+1 ", filenamePAIR, "+1",sep="")
                #system(comm)
                comm <- paste(headas," && ", src," && fdump infile=pp_novar.fits", "+1 outfile=", filenamePAIR, 
                              " columns='ADC' rows=1 prhead=no showcol=no showrow=no showunit=no",sep="")
                system(comm)
            }
            # does not WORK!
            # kk <- file(description = filenamePAIR, open = "rb")
            # header0 <- readFITSheader(kk, fixHdr = 'remove') # read primary header
            # header <- readFITSheader(kk, fixHdr = 'remove') # read extension header
            # Dpair <- readFITSbintable(kk, header)
            # close(kk)
            # idcol <- which(Dpair$colNames == "ADC")
            # DatakeV_pairSEP <- Dpair$col[[idcol]][1,]
            
            DatakeV_pairSEP<-read.table(filenamePAIR, header=F)
            istart <- 1000 # zero samples before pulse Start
            if (EkeV1 == "0.2") istart = 999
            pulsekeV_pairSEP <- DatakeV_pairSEP[-(1:(istart-initsamp)),1] - baseline
            
            offset <- append(zeros,pulse[[EkeV2]])
            DataSum <- pulse[[EkeV1]][1:pulseLength] + offset[1:pulseLength]
            derivDataSum <-diff(DataSum)
            derivPulse1 <- diff(pulse[[EkeV1]])
            derivPulse2 <- diff(pulse[[EkeV2]])
            derivPairs <- diff(pulsekeV_pairSEP)   
            derivOffset<- diff(offset)
            derivSumMinusDerivPulse1 <- derivDataSum[1:derivLength]-derivPulse1[1:derivLength]
            derivPairsMinusDerivPulse1 <- derivPairs[1:derivLength]-derivPulse1[1:derivLength]
            derivSumMinusDerivPulse12 <- derivDataSum[1:derivLength]-derivPulse1[1:derivLength]-derivOffset[1:derivLength]
            derivOffsetScaled <- max(derivPairsMinusDerivPulse1[1:derivLength])/max(derivOffset[1:derivLength])*derivOffset[1:derivLength]
            derivPairsMinusDerivPulse12 <- derivPairsMinusDerivPulse1[1:derivLength] - derivOffsetScaled[1:derivLength]
            
            # find new "closest"energy model for secondary pulse based on the 1st sample of dreivative
            newDer1 <- derivPairsMinusDerivPulse1[initsamp+sep]
            newEkeV2idx <- whichClosest(newDer1, Der1)
            newEkeV2 <- libEnergies[newEkeV2idx]
            newPulse2 <- pulse[[newEkeV2]] + (newDer1-Der1[newEkeV2idx])/(Der1[(newEkeV2idx+1)]-Der1[newEkeV2idx])*
                (pulse[[libEnergies[newEkeV2idx+1]]]-pulse[[newEkeV2]])
            newOffset <- append(zeros,newPulse2)
            derivNewOffset <- diff(newOffset)
            derivPairsMinusDerivPulse1new2 <- derivPairsMinusDerivPulse1[1:derivLength] - derivNewOffset[1:derivLength]
            
            
            colSumm <- "blue"
            colPairs <- "red"
            
            #--------------------------------------------------
            # PLOT 2 pulses: SUMM and PAIR simulation (tessim)
            #--------------------------------------------------
            plotLength <- 600
            plot(seq(1,pulseLength),DataSum, type="l", col=colSumm, 
                  xlim=c(1,plotLength),ylab = "ADC", xlab="sample", main=paste("Pulses ",EkeV1," & ",EkeV2, " keV",sep="") )
            lines(seq(1,pulseLength),pulsekeV_pairSEP[1:pulseLength], type="l", col=colPairs)
            lines(seq(1,pulseLength),pulse[[EkeV1]][1:pulseLength], lty=2,cex=0.1)
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
            abline(h=sigma6, lty=2, col="grey")
            text(80,sigma6,expression(6*sigma),cex=0.6)
            abline(h=sigma11, lty=2, col="grey")
            text(90,sigma11,expression(11*sigma),cex=0.6)
            
            
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
            points((initsamp+sep),newDer1,pch=1, cex=2)
            arrows((initsamp+sep+2),newDer1,100,newDer1,length=0.1)
            text(150,newDer1,paste("1st sam of pulses of ",newEkeV2," keV",sep=""),cex=0.6)
            points((initsamp+sep),Der1[libEnergies==EkeV2],pch=1, cex=2)
            arrows((initsamp+sep),Der1[libEnergies==EkeV2],100,Der1[libEnergies==EkeV2],length=0.1)
            text(150,Der1[libEnergies==EkeV2],paste("1st sam of pulses of ",EkeV2," keV",sep=""),cex=0.6)
            abline(h=sigma6, lty=2, col="grey")
            text(250,sigma6,expression(6*sigma),cex=0.6)
            abline(h=sigma11, lty=2, col="grey")
            text(260,sigma11,expression(11*sigma),cex=0.6)
            
            #-----------------------------------------------------
            # PLOT 2 pulses DERIVATIVE - derivPulse1 - derivPulse2
            #-----------------------------------------------------
            plotLength <- 300
            plot(seq(1,derivLength),derivSumMinusDerivPulse12[1:derivLength], pch=2, col="blue", 
                 xlim=c(0,plotLength),ylab ="Deriv. all - Deriv. Pulse12", main=paste("Pulses ",EkeV1," & ", EkeV2," keV",paste=""), 
                 xlab="sample", cex=0.4, ylim=c(-50,200))
            points(seq(1,derivLength),derivPairsMinusDerivPulse12[1:derivLength], pch=2, cex=0.4,col="orange")
            points(seq(1,derivLength),derivPairsMinusDerivPulse1new2[1:derivLength], pch=2, cex=0.4,col="red")
            mycols <- colorRampPalette(c("green","magenta", "cyan"))(5)
            loc02 <- 50
            loc02_2 <- loc02 -2 
            derivPulse02 <- diff(pulse[["0.2"]])
            points(seq(loc02,pulseLength),derivPulse02[1:(derivLength-loc02_2)], pch=4, cex=0.4, col=mycols[1])
            legend("topright", 
                   c(paste("Deriv of Summed -Deriv of Single1 (",EkeV1," keV) & Single2 (",EkeV2," keV)",sep=""),
                     paste("Deriv of Pairs -Deriv of Single1  (",EkeV1," keV) & Single2 (",newEkeV2," keV)",sep=""),
                     paste("Deriv of Pairs -Deriv of Single1  (",EkeV1," keV) & Single2 (",EkeV2," keV, scaled)",sep=""),
                     "Deriv of Pulse02keV"),
                   col=c("blue","red",mycols), pch=c(2,2,2,4),cex=0.6, bty="n")
            # legend("topright", 
            #        c("Deriv of Summed -Deriv of Single1 & single2",
            #          "Deriv of Pairs -Deriv of Single1 & single2(scaled)",
            #          "Deriv of Pulse02keV","Deriv of Pair02keV_sep1",
            #          "Deriv of Pair02keV_sep10","Deriv of Pair02keV_sep56",
            #          "Deriv of Pair02keV_sep101"),
            #        col=c("blue","red",mycols), pch=c(2,2,4,4,4,4,4),cex=0.8, bty="n")
            
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
        } # EkeV2
    } #EkeV1
}#sep






