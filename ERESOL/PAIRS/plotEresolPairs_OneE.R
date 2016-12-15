#
#  Plot Energy Resolution Curves based on the analysis 
#     (done by getEresolCurves.csh) of simulated pairs of pulses with
#     tesconstpileup + tessim + SIRENA
#
#  
#
#
#-------------------- VARIABLE PARAMS ---------------------
labelLib <- "multilib" # "monolib" or "multilib" or nonoiselib
#labelLib <- "monolib"
array <- "SPA"
dataDir <- paste("eresol",array,sep="")
#fmeths <- c("F0","B0")
fmeths <- c("F0")
energy <- "1keV" # simulated monochromatic energy
nsamples<-1024
fdomain<-"F"
energyMethods <-c("LAGS","NOLAGS")
energyMethods<-c("NOLAGS_NTRIG")
PAIRSdir <- "/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS"
sepMin<-4
sepMax<-10000
FWHMmin<-1.6
FWHMmax<-20
BIASmin<--500
BIASmax<-50

#----------------------------------------------------------

samprate <- 156250 #Hz
exptime<-100 #s

if(array == "SPA"){
    FWHMmax <- 15
    lowRes <-128
}else if(array == "LPA1"){
    FWHMmax <- 15
    lowRes <-256
        
}else if(array == "LPA2"){
    FWHMmax <- 40
    lowRes <- 512
}else if(array == "LPA3"){
    #FWHMmax <- 50
    lowRes <- 1024
}
if(labelLib == "multilib"){
    libLegend<-" Multi E library"
    baseFWHM <- 4
}else{
    libLegend<-" Mono E library"
    baseFWHM <- 2
}
FWHMtitle<-paste("Energy Resolution FWHM (",array,")",sep="")
EBIAStitle <-paste("Energy BIAS (",array,")",sep="")
EBIASlabel <-expression(paste(E[BIAS]," [eV] =",frac(sum(E[i]-E[input]),N[pulses]),sep=""))
FWHMlabel <-expression(paste("FWHM [eV] = ",
                             2.35 %*% sqrt(frac(sum(((E[i]-E[BIAS])-E[input])^2),N[pulses])), sep=""))

pdfName <- paste(PAIRSdir,"/",dataDir,"/plotEresolPairs_",array,"_SIRENA",nsamples,"_",energy,"_B0F0",
                 fdomain,"_",labelLib,"_",energyMethod,".pdf",sep="")
#pdf(pdfName,width=10.,height=7.1) 
#pdf(pdfName,paper="a4r")
pdf(pdfName,width=10,height=7) 


###########################
# PLOT FWHM
###########################

par(mar=c(5,6,4,2))
drawLogPlotBox(xlimits=c(sepMin,sepMax), ylimits=c(FWHMmin,FWHMmax),
               logxy="x", xlabel="Pulse distance [samples]", 
               ylabel=FWHMlabel, naxes=c(T,T,F,F))
title(main=paste(FWHMtitle,libLegend))
axis(2,at=c(seq(1:10)),labels=F,lty=1)
#abline(v=lowRes,lty=2)
#rect(1E-3,-100,lowRes,FWHMmax+10,col="grey96",border=NA)
#text(lowRes-100,FWHMmin+(FWHMmax-FWHMmin)/2,"Low res events",cex=1.3, col="grey",srt=90)

#abline(h=baseFWHM,lty=5,col="grey")

legendCol <- c()
legendStr <- c()
legendLtype <- c()
legendPch <-c()

for (fmeth in fmeths){

    if (fmeth == "F0"){
        PrimCol <- "blue"
        SecCol <- "red"
        ltype<-0
    }else if(fmeth == "B0"){
        PrimCol <- "cyan"
        SecCol <- "green4"
        ltype<-3
    }
    for (energyMethod in energyMethods){
        if (energyMethod == "LAGS"){
            points <- 1
        }else if (energyMethod == "NOLAGS"){
            points <- 16
        }
            
            
        # READ data
        #--------------
        data <- paste(PAIRSdir,"/",dataDir,"/eresol_",exptime,"s_SIRENA",nsamples,"_",energy,"_",
                      fmeth,fdomain,"_",labelLib,"_",energyMethod,".dat",sep="")
        cat("Reading data for FWHM:",data,"\n")
        dataSIRENA <- read.table(data,header=T,stringsAsFactors = FALSE)
    
#         nsamplesUp <- c(2)
#         nSgms <- c(10)
#         if(ncol(dataSIRENA)>9){
#             samplesUp <- c(dataSIRENA[1,10])
#             nSgms     <- c(dataSIRENA[1,11])
#         }
        #labelPrim <- paste("Primary Pulses (",fmeth, ") E=",energy,", sU=",samplesUp[1]," nSg=",nSgms[1],sep="")
        #labelSec  <- paste("Secondary Pulses (",fmeth, ") E=",energy,", sU=",samplesUp[1]," nSg=",nSgms[1],sep="")
        labelPrim <- paste("Primary Pulses (",fmeth,fdomain,") E=",energy,", ",energyMethod,sep="")
        labelSec  <- paste("Secondary Pulses (",fmeth,fdomain,") E=",energy,", ",energyMethod,sep="")
    
        samplesSep <- dataSIRENA[,1]
        pulseDistance <- (samplesSep/samprate)*1E6  #mus
        #taus <- pulseDistance/tf 
    
        
        # read BIAS corrected FWHM
        FWHM1SIRENA <- dataSIRENA[,2]
        FWHM2SIRENA <- dataSIRENA[,3]
            
        ylimitsPri <- c(min(FWHM1SIRENA),max(FWHM1SIRENA))
        ylimitsSec <- c(min(0,FWHM2SIRENA),max(FWHM2SIRENA,6))*1.1
        ytextPri <- ylimitsPri[2]
        ytextSec <- ylimitsSec[2]

        # Read FWHM err
        FWHM1ERRSIRENA <- dataSIRENA[,8]
        FWHM2ERRSIRENA <- dataSIRENA[,9]    


        ############ PRIMARIES #######################    
        points(samplesSep,FWHM1SIRENA,pch=points, col=PrimCol,lty=ltype,type="b",cex=1)
        arrows(samplesSep, FWHM1SIRENA-FWHM1ERRSIRENA, 
               samplesSep, FWHM1SIRENA+FWHM1ERRSIRENA, 
               length=0.05, angle=90, code=3, col=PrimCol)
        
        ############ SECONDARIES #######################
        points(samplesSep,FWHM2SIRENA,pch=points, col=SecCol,lty=ltype,type="b",cex=1)
        arrows(samplesSep, FWHM2SIRENA-FWHM2ERRSIRENA, 
               samplesSep, FWHM2SIRENA+FWHM2ERRSIRENA, 
               length=0.05, angle=90, code=3, col=SecCol)
    
        legendStr <- c(legendStr,labelPrim,labelSec)
        legendCol <- c(legendCol,PrimCol,SecCol)
        legendLtype <- c(legendLtype,ltype,ltype)
        legendPch <- c(legendPch,points,points)
        
    } # energyMethods
} # fmeths for FWHM   
legend("topright",legend=legendStr, pch=legendPch,col=legendCol,
                   lty=legendLtype,cex=0.8)          

# subplot(
#     points(samplesSep,FWHM1SIRENA,pch=points, col=PrimCol,lty=ltype,type="b",cex=1), 
#     x=grconvertX(c(0.75,1), from='npc'),
#     y=grconvertY(c(0,0.25), from='npc'),
#     pars=list( mar=c(1.5,1.5,0,0)+0.1) )
#subplot(plot(samplesSep,FWHM1SIRENA,pch=points, col=PrimCol,lty=ltype,type="b",cex=1),
#        x=500,y=10)
# subplot( hist(rnorm(100)), 
#          x=grconvertX(c(0.75,1), from='npc'),
#          y=grconvertY(c(0,0.25), from='npc'),
#          pars=list( mar=c(1.5,1.5,0,0)+0.1) )
# 



##################################
# PLOT BIAS
##################################
drawLogPlotBox(xlimits=c(sepMin,sepMax), ylimits=c(BIASmin,BIASmax),
               logxy="x", 
               xlabel="Pulse Distance [samples]",
               ylabel=EBIASlabel,
               naxes=c(T,T,F,F))
title(main=EBIAStitle)
#abline(v=lowRes,lty=2)
#rect(1E-3,-100,lowRes,1500,col="grey96",border=NA)
#text(lowRes-100,BIASmin+(BIASmax-BIASmin)/2,"Low res events",cex=1.3, col="grey",srt=90)

legendCol <- c()
legendStr <- c()
legendLtype <- c()
legendPch <-c()

for (fmeth in fmeths){
    
    if (fmeth == "F0"){
        PrimCol <- "blue"
        SecCol <- "red"
        ltype<-0
    }else if(fmeth == "B0"){
        PrimCol <- "cyan"
        SecCol <- "green4"
        ltype<-3
    }
    for (energyMethod in energyMethods){
        if (energyMethod == "LAGS"){
            points <- 1
        }else if (energyMethod == "NOLAGS"){
            points <- 16
        }
        
        # READ data
        #--------------
        data <- paste(PAIRSdir,"/",dataDir,"/eresol_",exptime,"s_SIRENA",nsamples,"_",energy,"_",
                      fmeth,fdomain,"_",labelLib,"_",energyMethod,".dat",sep="")
        cat("Reading data for BIAS:",data,"\n")
        dataSIRENA <- read.table(data,header=T,stringsAsFactors = FALSE)
        
#         nsamplesUp <- c(2)
#         nSgms <- c(10)
#         if(ncol(dataSIRENA)>9){
#             samplesUp <- c(dataSIRENA[1,10])
#             nSgms     <- c(dataSIRENA[1,11])
#         }
        

        #
        # PLOT BIAS 
        #
        BIAS1SIRENA <- dataSIRENA[,4]
        BIAS2SIRENA <- dataSIRENA[,5]
        maxy <- max(BIAS1SIRENA,BIAS2SIRENA)
        miny <- min(BIAS1SIRENA,BIAS2SIRENA)
        
        # Primaries
        points(samplesSep,BIAS1SIRENA,pch=points,col=PrimCol,lty=ltype,type="b",cex=1)
        # Secondaries
        points(samplesSep,BIAS2SIRENA,pch=points,col=SecCol,lty=ltype,type="b",cex=1)
        
        legendStr <- c(legendStr,
                       paste("Primary pulses ",  fmeth,fdomain, " (E = ",energy,", ",energyMethod,")",sep=""),
                       paste("Secondary pulses ",fmeth,fdomain, " (E = ",energy,", ",energyMethod,")",sep=""))
        legendCol <- c(legendCol,PrimCol,SecCol)
        legendLtype <- c(legendLtype,ltype,ltype)
        legendPch <- c(legendPch,points,points)
        
    }# energyMethods
}#fmeths

legend("topright",legend=legendStr, pch=legendPch,col=legendCol,
       lty=legendLtype,cex=0.8)

dev.off()

