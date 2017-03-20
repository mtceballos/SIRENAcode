# Read eresol values for middle pulses in trios and plots image 
# for separations from the primary and to the secondary
library("fields")
library(rjson)
library("colorRamps")

array  <- "LPA2shunt"
eresolDir <- paste("/home/ceballos/INSTRUMEN/EURECA/ERESOL/TRIOS/eresol",array,sep="")
fmeth  <- "F0"
lib    <- "fixedlib1OF"
EkeVforTrios <- "7keV"
nSimPulses <- "30000"
nSamples <- "4096"
pulseLength <-"4096"
plotBias <- "N"
TRIGG <- "_NTRIG"
rootname <- paste("eresol_",nSimPulses,"p_SIRENA",nSamples,"_pL",pulseLength,"_",EkeVforTrios,
                  "_F0F_",lib,"_OPTFILT",TRIGG,sep="")

# GRADING TABLE from SPIE 2016
gradingTable <- data.frame(ARRAY=c("LPA1shunt","LPA2shunt","SPA"),
                breakBIAS=c(400,700,190), breakHIGHRES=c(1024,16384,512),
                breakMIDRES=c(256,512,128),stringsAsFactors=FALSE)
biasBreak<-gradingTable[gradingTable[,1]==array,"breakBIAS"]
HRbreak<-gradingTable[gradingTable[,1]==array,"breakHIGHRES"]
MRbreak<-gradingTable[gradingTable[,1]==array,"breakMIDRES"]



if(array == "SPA"){
    #nsamples <- 1024
    #biasBreak <- 128.
    #HRbreak <- 512
    #MRbreak <- 128
    biasLabelCol<-"white"
    biasLabelPos<-c(20,100)
    if(EkeVforTrios == "6keV")biasLabelPos<-c(50,200)
    HRlabelCol <- "white"
    HRlabelPos <- c(500,1000)
    MRlabelCol <- "white"
    MRlabelPos <- c(500,200)
    #LRlabelCol <- "white"
    LRlabelCol <- "black"
    #LRlabelPos <- c(500,30)
    LRlabelPos <- c(500,50)
    FWHMmin   <- 2
    FWHMmax   <- 8
    if(EkeVforTrios == "6keV") FWHMmax <- 6 #6keV
}else if(length(grep("LPA1",array))>0){    
    biasLabelCol<-"white"
    if(EkeVforTrios=="1keV"){
        biasLabelPos<-c(50,200)
        HRlabelCol <- "white"
        HRlabelPos <- c(2000,2000)
        MRlabelCol <- "white"
        MRlabelPos <- c(2000,500)
        LRlabelCol <- "white"
        LRlabelPos <- c(2000,50)
    }else if(EkeVforTrios=="7keV"){
        biasLabelPos<-c(50,500)
        HRlabelCol <- "white"
        HRlabelPos <- c(2000,2000)
        MRlabelCol <- "white"
        MRlabelPos <- c(2000,500)
        LRlabelCol <- "white"
        LRlabelPos <- c(2000,100)
    }
    FWHMmin   <- 2
    FWHMmax   <- 16
    sepsStr = c('00004', '00005', '00007', '00010', '00013', '00017', '00023', '00031', '00042', 
                '00056', '00075', '00101', '00136', '00182', '00244', '00328', '00439', '00589', 
                '00791', '01061', '01423', '01908', '02560', '03433', '04605', '06178', '08287')
    seps12 <- as.numeric(sepsStr)
    seps23 <- as.numeric(sepsStr)
    
}else if(array == "LPA2shunt"){
    #nsamples <- 1024
    #biasBreak <- 800.
    #HRbreak <- 16384
    #MRbreak <- 512
    biasLabelCol<-"white"
    if(EkeVforTrios=="1keV"){
        biasLabelPos<-c(50,500)
        HRlabelCol <- "white"
        HRlabelPos <- c(5000,18000)
        MRlabelCol <- "white"
        MRlabelPos <- c(5000,4000)
        LRlabelCol <- "white"
        LRlabelPos <- c(5000,30)
    }else if(EkeVforTrios=="7keV"){
        biasLabelPos<-c(50,500)
        HRlabelCol <- "white"
        HRlabelPos <- c(5000,18000)
        MRlabelCol <- "white"
        MRlabelPos <- c(5000,4000)
        LRlabelCol <- "white"
        LRlabelPos <- c(5000,30)
    }
    FWHMmin   <- 2
    FWHMmax   <- 16
    sepsStr = c('00050', '00061', '00076', '00093', '00114', '00140', '00173', '00212', '00261', 
               '00321', '00395', '00485', '00597', '00733', '00902', '01109', '01363', '01676', 
               '02061', '02534', '03115', '03830', '04709', '05790', '07119', '08752', '10761', 
               '13231', '16267', '20000')
    seps12 <- as.numeric(sepsStr)
    seps23 <- as.numeric(sepsStr)
}else if(array == "LPA3"){
    #nsamples <- 2048
    #biasBreak <- 1400.
    #HRbreak <- 16384
    #MRbreak <- 1024
    biasLabelCol<-"white"
    if(energy=="1keV"){
        biasLabelPos<-c(100,500)
        HRlabelCol <- "white"
        HRlabelPos <- c(5000,18000)
        MRlabelCol <- "white"
        MRlabelPos <- c(5000,5000)
        LRlabelCol <- "white"
        LRlabelPos <- c(5000,100)    
    }else if(energy=="6keV"){
        biasLabelPos<-c(100,500)
        HRlabelCol <- "white"
        HRlabelPos <- c(5000,18000)
        MRlabelCol <- "white"
        MRlabelPos <- c(5000,5000)
        LRlabelCol <- "white"
        LRlabelPos <- c(5000,100)    
    }
    FWHMmin   <- 2
    FWHMmax   <- 60 #1keV
    FWHMmax   <- 90 # 6keV
}


fwhmSecGAINCORRE  <- matrix(NA,nrow=length(seps12), ncol=length(seps23)) # FWHM of Ecorr
biasErealSec      <- matrix(NA,nrow=length(seps12), ncol=length(seps23)) # FWHM of Ecorr
eresolFile <- paste(eresolDir,"/",rootname,".json",sep="")
#eresolFile <- paste("eresol_",nSimPulses,"p_SIRENA",nSamples,"_", EkeVforTrios,"_F0F_",lib,"_OPTFILT",TRIGG,".json",sep="")
jsondata <- fromJSON(file=eresolFile)
cat("Reading file ",eresolFile,"\n")

for (s12 in 1:length(seps12)){
    for (s23 in 1:length(seps23)){
        idxSep12 <- which(sapply(jsondata,function(x) x$separation12)==as.numeric(sepsStr[s12]))
        idxSep23 <- which(sapply(jsondata,function(x) x$separation23)==as.numeric(sepsStr[s23]))
        idxSep <- intersect(idxSep12,idxSep23)
        stopifnot(idxSep>0)
        fwhmSecGAINCORRE[s12,s23] <- as.numeric(jsondata[[idxSep]]$fwhmEreal$secondaries)
        biasErealSec[s12,s23]     <- as.numeric(jsondata[[idxSep]]$biasEreal$secondaries)
    }
}

#data <- read.table(paste("/home/ceballos/INSTRUMEN/EURECA/ERESOL/TRIOS/",eresol,
#        "/eresol_100s_SIRENA",nsamples,"_",energy,"_",fmeth,"_",lib,".dat",sep=""),header=T)
#data <- read.table(paste("/home/ceballos/INSTRUMEN/EURECA/ERESOL/TRIOS/eresolDAT/",eresol,
#                         "/eresol_100s_SIRENA",nsamples,"_",energy,"_",fmeth,"_",lib,".dat",sep=""),header=T)

pdfName <- paste(eresolDir,"/PDFs/plotTrios_",rootname,".pdf",sep="")
nseps12 <- length(seps12)
nseps23 <- length(seps23)
minFWHM<-min(fwhmSecGAINCORRE)  
maxFWHM<-max(fwhmSecGAINCORRE)
##############
minzFWHM<-minFWHM 
maxzFWHM<-maxFWHM
#minzFWHM<-FWHMmin #pre-fixed limit
#maxzFWHM<-FWHMmax #pre-fixed limit
##############
# labels and scaleing for the color bar
ticks<-10**(seq(from=log10(minzFWHM),to=log10(maxzFWHM),length.out=10))
ticksLabels<-sprintf("%2.1f",ticks)

mat_data_eresol <- matrix(data=fwhmSecGAINCORRE, nrow=nseps12, ncol=nseps23,byrow=T)
mat_data_ebias  <- matrix(data=biasErealSec, nrow=nseps12, ncol=nseps23,byrow=T)

pdf(pdfName,width=7.,height=7.)
colorRamp = blue2green2red(900)
ColorLevels <- seq(minzFWHM,maxzFWHM,length=length(colorRamp))
#
# PLOT IMAGE FOR FWHM
#
par(mar=c(5,5,5,7))
# create template for image plot
image(x=seps12,y=seps23,z=log10(mat_data_eresol), axes=FALSE, log="xy",
      zlim=c(log10(minzFWHM),log10(maxzFWHM)),col=tim.colors(400),
      xlab="Time since previous pulse [samples]",
      ylab="Time until next pulse [samples]",
      main=paste("Energy Resolution FWHM [eV] - ",array," - (",EkeVforTrios,")",sep="")) 
box()
# add new log axes and new log labels (ugly in image.plot)
# X axis
expseq<-seq(from=round(log10(min(seps12))),to=round(log10(max(seps12))))
xmjt<-10**expseq
xmjtlabs<-c()
xmnt<-c()
for (expo in expseq){
    xmjtlabs<-append(xmjtlabs,as.expression(bquote(10^.(expo))))
    xmnt<-append(xmnt,seq(2,9)*10**expo)
}
axis(1,at=xmjt,labels=xmjtlabs,las=1,col.axis="black")
axis(1,at=xmnt,tcl=-0.2,labels=FALSE,col.axis="white")
# Y axis
expseq<-seq(round(log10(min(seps23))),round(log10(max(seps23))))
ymjt<-10**expseq
ymjtlabs<-c()
ymnt<-c()
for (expo in expseq){
    ymjtlabs<-append(ymjtlabs,as.expression(bquote(10^.(expo))))
    ymnt<-append(ymnt,seq(2,9)*10**expo)
}
axis(2,at=ymjt,labels=ymjtlabs,las=1)
axis(2,at=ymnt,tcl=-0.2,labels=FALSE,col.axis="white")

# add values color-coded, axes labels and legend
image.plot(x=seps12,y=seps23,z=log10(mat_data_eresol), axes=FALSE,
           axis.args=list( at=log10(ticks),  labels=ticksLabels),
           zlim=c(log10(minzFWHM),log10(maxzFWHM)),col=tim.colors(400),
           xlab="", ylab="", log="xy", add=TRUE)
abline(v=biasBreak,col="white",lwd=3)
segments(x0=biasBreak,y0=MRbreak,x1=30000,y1=MRbreak,col="white",lwd=3)
segments(x0=biasBreak,y0=HRbreak,x1=30000,y1=HRbreak,col="white",lwd=3)
text(biasLabelPos[1],biasLabelPos[2],"Rejected Secondaries",col="white",cex=1.6,srt=45)
text(LRlabelPos[1],LRlabelPos[2],"Low Res",col=LRlabelCol,cex=1.2)
text(MRlabelPos[1],MRlabelPos[2],"Mid Res",col=MRlabelCol,cex=1.2)
text(HRlabelPos[1],HRlabelPos[2],"High Res",col=HRlabelCol,cex=1.2)



if(plotBias == "Y"){
    #EBIAS
    minz <- min(biasErealSec)
    maxz <- max(biasErealSec)
    maxz<-100
    minz<--100
    #image.plot(x=seps,y=seps,z=mat_data_ebias, zlim=c(minz,maxz),col=colorRamp,
    image.plot(x=seps12,y=seps23,z=mat_data_ebias, zlim=c(minz,maxz),col=heat.colors(400),
               xlab="Time from previous pulse [samples]",
               ylab="Time to next pulse [samples]",
               main=paste("Energy BIAS [eV] - ",array," - (",EkeVforTrios,")",sep=""),log="xy")
    #abline(v=biasBreak,col="white",lwd=3)
    #segments(x0=biasBreak,y0=MRbreak,x1=30000,y1=MRbreak,col="white",lwd=3)
    #segments(x0=biasBreak,y0=HRbreak,x1=30000,y1=HRbreak,col="white",lwd=3)
    #text(biasLabelPos[1],biasLabelPos[2],"Rejected (high BIAS)",col="white",cex=1.6,srt=45)
    #text(LRlabelPos[1],LRlabelPos[2],"Low Res",col=LRlabelCol,cex=1.2)
    #text(MRlabelPos[1],MRlabelPos[2],"Mid Res",col=MRlabelCol,cex=1.2)
    #text(HRlabelPos[1],HRlabelPos[2],"High Res",col=HRlabelCol,cex=1.2)
    
} #plotBias?
dev.off()
