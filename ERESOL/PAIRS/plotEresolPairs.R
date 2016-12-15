#
#  Plot Energy Resolution Curves based on the analysis 
#     (done by getEresolCurves.csh) of simulated pairs of pulses with
#     tesconstpileup + tessim (defaults)
#
#  SIRENA: no base 2 filter calculation; selection of model by max of derivative
#        POSSIBILITY OF USING MONOE LIBRARY (no energy selection)
#
#
#-------------------- VARIABLE PARAMS ---------------------
#labelLib <- "multilib" # "monolib" or "multilib" or nonoiselib
labelLib <- "monolib"
#labelLib <- "nonoiselib"
#dataDir <- "tessimv0"
#dataDir <- "tessim20150505"
#dataDir <- "tessim20150528"
#dataDir <- "tessim20150505m0"
#dataDir <- "tessim20150505Nonoise"
dataDir <- "tessimLPA1"
#fmeths <- c("F0","B0")
fmeths <- c("F0")
#----------------------------------------------------------

samprate <- 156250 #Hz
exptime<-100 #s
nsamples<-1024
if(labelLib == "multilib"){
    libLegend<-" Multi E library"
}else{
    libLegend<-" Mono E library"
}
#if(dataDir=="tessimv0"){
if(dataDir=="tessimLPA"){
    array <- "LPA"
#}else if(dataDir=="tessim20150505"){
}else if(dataDir=="tessimSPA"){
    array <- "SPA"
}else if(dataDir=="tessim20150505m0"){
    array <- "SPA (m0)"
}else if(dataDir=="tessim20150528"){
    array <- "SPA (tessim excess noise corrected)"
}
FWHMtitle<-paste("Energy Resolution FWHM (",array,")",sep="")
EBIAStitle <-paste("Energy BIAS (",array,")",sep="")
EBIASlabel <-expression(paste(E[BIAS]," [eV] =",frac(sum(E[i]-E[input]),N[pulses]),sep=""))
energies <- c("1keV","3keV","6keV") # simulated monochromatic energies

pdfName <- paste("./",dataDir,"/plotEresolPairs_",dataDir,"_",labelLib,".pdf",sep="")
#pdf(pdfName,width=10.,height=7.1) 
#pdf(pdfName,paper="a4r")
pdf(pdfName,width=10,height=7) 

for (fmeth in fmeths){
    if (fmeth == "B0"){
        PrimCol <- "blue"
        SecCol <- "red"
    }else if(fmeth == "F0"){
        PrimCol <- "cyan"
        SecCol <- "green4"
    }
    # READ data
    #--------------
    dataSIRENA_1keV <- read.table(paste("./",dataDir,"/eresol_",exptime,"s_SIRENA1024_1keV_",fmeth,"_",labelLib,".dat",sep=""),
                                         header=T,stringsAsFactors = FALSE)
    dataSIRENA_3keV <- read.table(paste("./",dataDir,"/eresol_",exptime,"s_SIRENA1024_3keV_",fmeth,"_",labelLib,".dat",sep=""),
                                         header=T,stringsAsFactors = FALSE)
    dataSIRENA_6keV <- read.table(paste("./",dataDir,"/eresol_",exptime,"s_SIRENA1024_6keV_",fmeth,"_",labelLib,".dat",sep=""),
                                         header=T,stringsAsFactors = FALSE)

    nsamplesUp <- c(2,2,2)
    nSgms <- c(10,10,20)
    if(ncol(dataSIRENA_1keV)>9){
        samplesUp <- c(dataSIRENA_1keV[1,10],dataSIRENA_3keV[1,10],dataSIRENA_6keV[1,10])
        nSgms     <- c(dataSIRENA_1keV[1,11],dataSIRENA_3keV[1,11],dataSIRENA_6keV[1,11])
    }
    labelPrim1 <- paste("Primary Pulses ",fmeth, "E=1keV, sU=",samplesUp[1]," nSg=",nSgms[1],sep="")
    labelPrim3 <- paste("Primary Pulses ",fmeth, "E=3keV, sU=",samplesUp[2]," nSg=",nSgms[2],sep="")
    labelPrim6 <- paste("Primary Pulses ",fmeth, "E=6keV, sU=",samplesUp[3]," nSg=",nSgms[3],sep="")
    labelSec1 <- paste("Secondary Pulses ",fmeth, "E=1keV, sU=",samplesUp[1]," nSg=",nSgms[1],sep="")
    labelSec3 <- paste("Secondary Pulses ",fmeth, "E=3keV, sU=",samplesUp[2]," nSg=",nSgms[2],sep="")
    labelSec6 <- paste("Secondary Pulses ",fmeth, "E=6keV, sU=",samplesUp[3]," nSg=",nSgms[3],sep="")
    
    samplesSep <- dataSIRENA_3keV[,1]
    pulseDistance <- (samplesSep/samprate)*1E6  #mus
    #taus <- pulseDistance/tf 

    #---------------------------------------------------------------------------------------
    # PLOT for different computation methods in SIRENA (BIAS corrected & UNcorrected values)
    #---------------------------------------------------------------------------------------
    for (biasCorr in c(T,F)){
        if(biasCorr){
            # read BIAS corrected FWHM
            FWHM1SIRENA_1keV <- dataSIRENA_1keV[,2]
            FWHM2SIRENA_1keV <- dataSIRENA_1keV[,3]
        
            FWHM1SIRENA_3keV <- dataSIRENA_3keV[,2]
            FWHM2SIRENA_3keV <- dataSIRENA_3keV[,3]
            
            FWHM1SIRENA_6keV <- dataSIRENA_6keV[,2]
            FWHM2SIRENA_6keV <- dataSIRENA_6keV[,3]
            
            ylimitsPri <- c(min(FWHM1SIRENA_1keV,
                                FWHM1SIRENA_3keV,
                                FWHM1SIRENA_6keV),
                            max(FWHM1SIRENA_1keV,
                                FWHM1SIRENA_3keV,
                                FWHM1SIRENA_6keV))
            ylimitsSec <- c(min(0,FWHM2SIRENA_1keV,
                                FWHM2SIRENA_3keV,
                                FWHM2SIRENA_6keV),
                            max(FWHM2SIRENA_1keV,
                                FWHM2SIRENA_3keV,
                                FWHM2SIRENA_6keV,6))*1.1
            ytextPri <- ylimitsPri[2]
            ytextSec <- ylimitsSec[2]
            biasLegend <- "(BIAS corrected)"
            FWHMlabel <-expression(paste("FWHM [eV] = ",
                        2.35 %*% sqrt(frac(sum(((E[i]-E[BIAS])-E[input])^2),N[pulses])), sep=""))
            
        }else{
            # read BIAS UNcorrected FWHM
            FWHM1SIRENA_1keV <- dataSIRENA_1keV[,6]
            FWHM2SIRENA_1keV <- dataSIRENA_1keV[,7]
            
            FWHM1SIRENA_3keV <- dataSIRENA_3keV[,6]
            FWHM2SIRENA_3keV <- dataSIRENA_3keV[,7]
            
            FWHM1SIRENA_6keV <- dataSIRENA_6keV[,6]
            FWHM2SIRENA_6keV <- dataSIRENA_6keV[,7]
            
            ylimitsPri <- c(min(FWHM1SIRENA_1keV,
                                FWHM1SIRENA_3keV,
                                FWHM1SIRENA_6keV),
                            max(FWHM1SIRENA_1keV,
                                FWHM1SIRENA_3keV,
                                FWHM1SIRENA_6keV))
            ylimitsSec <- c(min(0,FWHM2SIRENA_1keV,
                                FWHM2SIRENA_3keV,
                                FWHM2SIRENA_6keV),
                            max(FWHM2SIRENA_1keV,
                                FWHM2SIRENA_3keV,
                                FWHM2SIRENA_6keV,6))
            ytextPri <- ylimitsPri[2]
            ytextSec <- ylimitsSec[2]
            biasLegend <- "(BIAS UNcorrected)"
            FWHMlabel <-expression(paste("FWHM [eV] = ",
                 2.35 %*% sqrt(frac(sum((E[i]-E[input])^2),N[pulses]), sep="")))
            
            
        } # if biasCorr
    # Read FWHM err
    FWHM1ERRSIRENA_1keV <- dataSIRENA_1keV[,8]
    FWHM2ERRSIRENA_1keV <- dataSIRENA_1keV[,9]    
    FWHM1ERRSIRENA_3keV <- dataSIRENA_3keV[,8]
    FWHM2ERRSIRENA_3keV <- dataSIRENA_3keV[,9]    
    FWHM1ERRSIRENA_6keV <- dataSIRENA_6keV[,8]
    FWHM2ERRSIRENA_6keV <- dataSIRENA_6keV[,9]    
    
    #
    ############ PRIMARIES #######################
    #
    par(mar=c(5,6,4,2))
    drawLogPlotBox(xlimits=c(3,3000), ylimits=ylimitsPri,
                   logxy="x", 
                   xlabel="Pulse distance [samples]", 
                   ylabel=FWHMlabel,
                   naxes=c(T,T,F,F))
    title(main=paste(FWHMtitle,libLegend))
    
    points(samplesSep,FWHM1SIRENA_1keV,pch=1, col=PrimCol,lty=3,type="b",cex=1)
    points(samplesSep,FWHM1SIRENA_3keV,pch=10,col=PrimCol,lty=2,type="b",cex=1)
    points(samplesSep,FWHM1SIRENA_6keV,pch=19,col=PrimCol,lty=1,type="b",cex=1)
    arrows(samplesSep, FWHM1SIRENA_1keV-FWHM1ERRSIRENA_1keV, 
           samplesSep, FWHM1SIRENA_1keV+FWHM1ERRSIRENA_1keV, 
           length=0.05, angle=90, code=3, col=PrimCol)
    arrows(samplesSep, FWHM1SIRENA_3keV-FWHM1ERRSIRENA_3keV, 
           samplesSep, FWHM1SIRENA_3keV+FWHM1ERRSIRENA_3keV, 
           length=0.05, angle=90, code=3, col=PrimCol)
    arrows(samplesSep, FWHM1SIRENA_6keV-FWHM1ERRSIRENA_6keV, 
           samplesSep, FWHM1SIRENA_6keV+FWHM1ERRSIRENA_6keV, 
           length=0.05, angle=90, code=3, col=PrimCol)
    
    text(x=40,y=ytextPri,labels=paste("SIRENA+tessim PRIMARIES",biasLegend),col="black")
    legend("topright",legend=c(labelPrim1,labelPrim3,labelPrim6),
                               pch=c(1,10,19),col=c(PrimCol,PrimCol,PrimCol),
                               lty=c(3,2,1),cex=0.8)
    #
    ############ SECONDARIES #######################
    #
    par(mar=c(5,6,4,2))
    drawLogPlotBox(xlimits=c(3,3000), ylimits=ylimitsSec,
                   logxy="x", 
                   xlabel="Pulse Distance [samples]", 
                   ylabel=FWHMlabel,
                   naxes=c(T,T,F,F))
    title(main=paste(FWHMtitle,libLegend))
    
    points(samplesSep,FWHM2SIRENA_1keV,pch=1, col=SecCol,lty=3,type="b",cex=1)
    points(samplesSep,FWHM2SIRENA_3keV,pch=10,col=SecCol,lty=2,type="b",cex=1)
    points(samplesSep,FWHM2SIRENA_6keV,pch=19,col=SecCol,lty=1,type="b",cex=1)
    
    arrows(samplesSep, FWHM2SIRENA_1keV-FWHM2ERRSIRENA_1keV, 
           samplesSep, FWHM2SIRENA_1keV+FWHM2ERRSIRENA_1keV, 
           length=0.05, angle=90, code=3, col=SecCol)
    arrows(samplesSep, FWHM2SIRENA_3keV-FWHM2ERRSIRENA_3keV, 
           samplesSep, FWHM2SIRENA_3keV+FWHM2ERRSIRENA_3keV, 
           length=0.05, angle=90, code=3, col=SecCol)
    arrows(samplesSep, FWHM2SIRENA_6keV-FWHM2ERRSIRENA_6keV, 
           samplesSep, FWHM2SIRENA_6keV+FWHM2ERRSIRENA_6keV, 
           length=0.05, angle=90, code=3, col=SecCol)
    
    text(x=40,y=ytextSec,labels=paste("SIRENA+tessim SECONDARIES",biasLegend),col="black")
    legend("topright",legend=c(labelSec1,labelSec3,labelSec6),
           pch=c(1,10,19),col=c(SecCol,SecCol,SecCol),
           lty=c(3,2,1),cex=0.8)
    
           
           
    } # for  biasCorr
    
    #
    # PLOT BIAS 
    #
    BIAS1SIRENA_1keV <- dataSIRENA_1keV[,4]
    BIAS2SIRENA_1keV <- dataSIRENA_1keV[,5]
    BIAS1SIRENA_3keV <- dataSIRENA_3keV[,4]
    BIAS2SIRENA_3keV <- dataSIRENA_3keV[,5]
    BIAS1SIRENA_6keV <- dataSIRENA_6keV[,4]
    BIAS2SIRENA_6keV <- dataSIRENA_6keV[,5]
    maxy <- max(BIAS1SIRENA_1keV,BIAS2SIRENA_1keV,
                BIAS1SIRENA_3keV,BIAS2SIRENA_3keV,
                BIAS1SIRENA_6keV,BIAS2SIRENA_6keV)
    miny <- min(BIAS1SIRENA_1keV,BIAS2SIRENA_1keV,
                BIAS1SIRENA_3keV,BIAS2SIRENA_3keV,
                BIAS1SIRENA_6keV,BIAS2SIRENA_6keV)
    
    drawLogPlotBox(xlimits=c(3,3000), ylimits=c(miny,maxy),
                   logxy="x", 
                   xlabel="Pulse Distance [samples]",
                   ylabel=EBIASlabel,
                   naxes=c(T,T,F,F))
    title(main=EBIAStitle)
    # Primaries
    points(samplesSep,BIAS1SIRENA_1keV,pch=1,col=PrimCol,lty=1,type="p",cex=1)
    points(samplesSep,BIAS1SIRENA_3keV,pch=10,col=PrimCol,lty=1,type="p",cex=1)
    points(samplesSep,BIAS1SIRENA_6keV,pch=19,col=PrimCol,lty=1,type="p",cex=1)
    # Secondaries
    points(samplesSep,BIAS2SIRENA_1keV,pch=1,col=SecCol,lty=1,type="p",cex=1)
    points(samplesSep,BIAS2SIRENA_3keV,pch=10,col=SecCol,lty=1,type="p",cex=1)
    points(samplesSep,BIAS2SIRENA_6keV,pch=19,col=SecCol,lty=1,type="p",cex=1)
}
legend("topright",legend=c(paste("Primary pulses ",fmeth, " (E = 1 keV)",sep=""),
                           paste("Primary pulses ",fmeth, " (E = 3 keV)",sep=""),
                           paste("Primary pulses ",fmeth, " (E = 6 keV)",sep=""),
                           paste("Secondary pulses ",fmeth, " (E = 1 keV)",sep=""),
                           paste("Secondary pulses ",fmeth, " (E = 3 keV)",sep=""),
                           paste("Secondary pulses ",fmeth, " (E = 6 keV)",sep="")),
       pch=c(1,10,19,1,10,19),col=c(PrimCol,PrimCol,PrimCol,SecCol,SecCol,SecCol),
       lty=c(NA,NA,NA,2,2,2),cex=0.8)
dev.off()

