# Read eresol values for middle pulses in trios and plots image 
# for separations from the primary and to the secondary
library("fields")

array  <- "LPA3"
tessim <- paste("tessim",array,sep="")
fmeth  <- "F0"
lib    <- "monolib"
energy <- "6keV"
plotBias <- "N"
gradingTable <- read.table("~/INSTRUMEN/EURECA/ERESOL/gradingTable.dat",header=T)
nsamples<-gradingTable[gradingTable[,1]==array,"nsamples"]
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
    if(energy == "6keV")biasLabelPos<-c(50,200)
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
    if(energy == "6keV") FWHMmax <- 6 #6keV
}else if(array == "LPA1"){    
    #nsamples <- 1024
    #biasBreak <- 400.
    #HRbreak <- 1024
    #MRbreak <- 256
    biasLabelCol<-"white"
    if(energy=="1keV"){
        biasLabelPos<-c(50,200)
        HRlabelCol <- "white"
        HRlabelPos <- c(2000,2000)
        MRlabelCol <- "white"
        MRlabelPos <- c(2000,500)
        LRlabelCol <- "white"
        LRlabelPos <- c(2000,50)
    }else if(energy=="6keV"){
        biasLabelPos<-c(120,500)
        HRlabelCol <- "white"
        HRlabelPos <- c(2000,2000)
        MRlabelCol <- "white"
        MRlabelPos <- c(2000,500)
        LRlabelCol <- "white"
        LRlabelPos <- c(2000,100)
    }
    FWHMmin   <- 2
    FWHMmax   <- 16
}else if(array == "LPA2"){
    #nsamples <- 1024
    #biasBreak <- 800.
    #HRbreak <- 16384
    #MRbreak <- 512
    biasLabelCol<-"white"
    if(energy=="1keV"){
        biasLabelPos<-c(50,500)
        HRlabelCol <- "white"
        HRlabelPos <- c(5000,18000)
        MRlabelCol <- "white"
        MRlabelPos <- c(5000,4000)
        LRlabelCol <- "white"
        LRlabelPos <- c(5000,30)
    }else if(energy=="6keV"){
        biasLabelPos<-c(50,500)
        HRlabelCol <- "white"
        HRlabelPos <- c(5000,18000)
        MRlabelCol <- "white"
        MRlabelPos <- c(5000,4000)
        LRlabelCol <- "white"
        LRlabelPos <- c(5000,30)
    }
    FWHMmin   <- 2
    FWHMmax   <- 50
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


#data <- read.table(paste("/home/ceballos/INSTRUMEN/EURECA/ERESOL/TRIOS/",tessim,
#        "/eresol_100s_SIRENA",nsamples,"_",energy,"_",fmeth,"_",lib,".dat",sep=""),header=T)
data <- read.table(paste("/home/ceballos/INSTRUMEN/EURECA/ERESOL/TRIOS/eresolDAT/",tessim,
                         "/eresol_100s_SIRENA",nsamples,"_",energy,"_",fmeth,"_",lib,".dat",sep=""),header=T)
if(plotBias == "Y") {
    pdfName <- paste("/home/ceballos/INSTRUMEN/EURECA/ERESOL/TRIOS/",
                     tessim,"/plotEresolTrios_",tessim,"_",lib,"_",fmeth,"_",energy,".pdf",sep="")
}else{
    pdfName <- paste("/home/ceballos/INSTRUMEN/EURECA/ERESOL/TRIOS/eresolDAT/",
                     tessim,"/Eresol_",tessim,"_",lib,"_",fmeth,"_",energy,".pdf",sep="")
}

sep12 <- data[,1]
sep23 <- data[,2]
FWHM  <- data[,3]
ebias <- data[,5]
nseps <- sqrt(length(sep12))
seps <- sep23[1:nseps]
minFWHM<-min(FWHM)  
maxFWHM<-max(FWHM)
##############
minzFWHM<-minFWHM 
maxzFWHM<-maxFWHM
#minzFWHM<-FWHMmin #pre-fixed limit
#maxzFWHM<-FWHMmax #pre-fixed limit
##############
ticks<-10**(seq(from=log10(minzFWHM),to=log10(maxzFWHM),length.out=10))
ticksLabels<-sprintf("%2.1f",ticks)
mat_data_eresol <- matrix(data=FWHM, nrow=nseps, ncol=nseps,byrow=T)
mat_data_ebias  <- matrix(data=ebias, nrow=nseps, ncol=nseps,byrow=T)


#pdf(pdfName,width=10.,height=7.1) #landscape
pdf(pdfName,width=7.,height=7.)
library("colorRamps")
colorRamp = blue2green2red(900)
#colorRamp = matlab.like(900)
ColorLevels <- seq(minzFWHM,maxzFWHM,length=length(colorRamp))
#
# PLOT IMAGE FOR FWHM
#

image.plot(x=seps,y=seps,z=log10(mat_data_eresol), 
           axis.args=list( at=log10(ticks),                                                                labels=ticksLabels),
           zlim=c(log10(minzFWHM),log10(maxzFWHM)),col=tim.colors(400),
           xlab="Time from previous pulse [samples]",
           ylab="Time to next pulse [samples]",
           main=paste("Energy Resolution FWHM [eV] - ",array," - (",energy,")",sep=""),
           log="xy") 

abline(v=biasBreak,col="white",lwd=3)
#segments(x0=128,y0=64,x1=3000,y1=64,col="white",lwd=3)
#segments(x0=128,y0=1024,x1=3000,y1=1024,,col="white",lwd=3)
segments(x0=biasBreak,y0=MRbreak,x1=30000,y1=MRbreak,col="white",lwd=3)
segments(x0=biasBreak,y0=HRbreak,x1=30000,y1=HRbreak,col="white",lwd=3)
text(biasLabelPos[1],biasLabelPos[2],"Rejected Secondaries",col="white",cex=1.6,srt=45)
text(LRlabelPos[1],LRlabelPos[2],"Low Res",col=LRlabelCol,cex=1.2)
text(MRlabelPos[1],MRlabelPos[2],"Mid Res",col=MRlabelCol,cex=1.2)
text(HRlabelPos[1],HRlabelPos[2],"High Res",col=HRlabelCol,cex=1.2)

if(plotBias == "Y"){
    #EBIAS
    minz <- min(ebias)
    maxz <- max(ebias)
    maxz<-100
    minz<--100
    #image.plot(x=seps,y=seps,z=mat_data_ebias, zlim=c(minz,maxz),col=colorRamp,
    image.plot(x=seps,y=seps,z=mat_data_ebias, zlim=c(minz,maxz),col=heat.colors(400),
               xlab="Time from previous pulse [samples]",
               ylab="Time to next pulse [samples]",
               main=paste("Energy BIAS [eV] - ",array," - (",energy,")",sep=""),log="xy")
    abline(v=biasBreak,col="white",lwd=3)
    segments(x0=biasBreak,y0=MRbreak,x1=30000,y1=MRbreak,col="white",lwd=3)
    segments(x0=biasBreak,y0=HRbreak,x1=30000,y1=HRbreak,col="white",lwd=3)
    text(biasLabelPos[1],biasLabelPos[2],"Rejected (high BIAS)",col="white",cex=1.6,srt=45)
    text(LRlabelPos[1],LRlabelPos[2],"Low Res",col=LRlabelCol,cex=1.2)
    text(MRlabelPos[1],MRlabelPos[2],"Mid Res",col=MRlabelCol,cex=1.2)
    text(HRlabelPos[1],HRlabelPos[2],"High Res",col=HRlabelCol,cex=1.2)
    
} #plotBias?
dev.off()
