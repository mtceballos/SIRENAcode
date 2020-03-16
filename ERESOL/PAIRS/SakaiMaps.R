source("/disco07/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/readKeyword.R")

library("fields")
library(FITSio)

libraryFile <- "ALL" # "ALL" or "SHORT"
detectionMode <- c("AD","A1")
primEnergies <- c(0.2,1,2,4,6,8)
primEnergies_str <- c("0.2","1","2","4","6","8")
secondEnergies <- c(0.2,0.5,1,2,3,4,5,6,7,8)
secondEnergies_str <- c("0.2","0.5","1","2","3","4","5","6","7","8")
#separations <- c(4,10,20,40,44,60,90,100,120,200,250,300,400,500,600,800,1000,1600,2000)
#separations_str <- c("00004","00010","00020","00040","00044","00060","00090","00100","00120","00200","00250","00300","00400","00500","00600","00800","01000","01600","02000")
separations <- c(4,5,7,10,14,20,28,39,54,75,105,146,202,281,389,540,749,1039,1442,2000)
separations_str <- c("00004","00005","00007","00010","00014","00020","00028","00039","00054","00075","00105","00146","00202","00281","00389","00540","00749","01039","01442","02000")

minimum <- 1e6
maximum <- -1e6

totalPulses <- c()

detectedPulses  <- matrix(NA,nrow=length(secondEnergies), ncol=length(separations))  #Percentage of detected pulses (>0 if false pulses and <0 if losing pulses)

# for (idetectionMode in 1:length(detectionMode)){
#  for (ie in 1:length(primEnergies)){
#    for (je in 1:length(secondEnergies)){
#      index <- length(separations)
#      for (ke in length(separations_str):1){
#        if (libraryFile == "ALL"){
#          filename<-paste("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresolLPA2shunt/events_sep",separations_str[ke],"sam_2000p_",primEnergies_str[ie],"keV_",secondEnergies_str[je],"keV_",detectionMode[idetectionMode],".fits",sep="")
#         }else{
#           if (detectionMode[idetectionMode] == "AD"){
#             filename<-paste("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresolLPA2shunt/events_sep",separations_str[ke],"sam_2000p_",primEnergies_str[ie],"keV_",secondEnergies_str[je],"keV_",detectionMode[idetectionMode],"_SHORT.fits",sep="")
#           }else{
#             filename<-paste("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresolLPA2shunt/events_sep",separations_str[ke],"sam_2000p_",primEnergies_str[ie],"keV_",secondEnergies_str[je],"keV_",detectionMode[idetectionMode],".fits",sep="")
#           }
#        }
#        #cat(filename,"\n")
#        NAXIS2 <- readKeyword(filename,"NAXIS2")
#        NETTOT <- readKeyword(filename,"NETTOT")
#        totalPulses[ke] <- NETTOT
#        #cat("NAXIS2: ",str(NAXIS2),"\n")
#        #cat("NETTOT: ",str(NETTOT),"\n")
#        detectedPulses[je,index] <- NAXIS2*100/NETTOT
#        #cat("Pulsos=",str(NETTOT)," Detectados=",str(NAXIS2)," ",detectedPulses[je,index],"\n")
#        index <- index-1
#      }
#    }
# 
#    mat_detectedPulses <- matrix(data=detectedPulses, nrow=length(secondEnergies), ncol=length(separations),byrow=F)
#    # cat("min(detectedPulses)=",min(mat_detectedPulses)," ",log10(min(mat_detectedPulses)),"\n")
#    # cat("max(detectedPulses)=",max(mat_detectedPulses)," ",log10(max(mat_detectedPulses)),"\n")
# 
#    if (min(mat_detectedPulses)<minimum)
#    {
#      minimum <- min(mat_detectedPulses)
#    }
#    if (max(mat_detectedPulses)>maximum)
#    {
#      maximum <- max(mat_detectedPulses)
#    }
#  }
# }
# 
# cat("min=",minimum,"\n")
# cat("max=",maximum,"\n")

minimum <- 50
maximum <- 100
# minimum <- 27.7
# maximum <- 140.1

ticks<-10**(seq(from=log10(minimum),to=log10(maximum),length.out=10))
ticksLabels<-sprintf("%2.1f",ticks)  
par(mar=c(5,5,5,9))
pal = colorRampPalette(c("red","green4"))

file= paste("Sakai_ADA1_libraryALL.pdf",sep="")
if (libraryFile == "SHORT"){
  file= paste("Sakai_ADA1_librarySHORT.pdf",sep="")
}
pdf(file,width=24,height=14.)

for (idetectionMode in 1:length(detectionMode)){
  mat=matrix(data=1:6,nrow=2,ncol=3,byrow = T)
  layout(mat)
  for (ie in 1:length(primEnergies)){
    matrixFile <- paste("imageMatrix_",detectionMode[idetectionMode],"_",primEnergies[ie],"keV.mat",sep="")
    if (libraryFile == "SHORT"){
      matrixFile <- paste("imageMatrix_",detectionMode[idetectionMode],"_",primEnergies[ie],"keV_SHORT.mat",sep="")
    }
    if(file.exists(matrixFile)) {
        load(matrixFile)
    }else{
        for (je in 1:length(secondEnergies)){
          index <- length(separations)
          for (ke in length(separations_str):1){  
            if (libraryFile == "ALL"){
              filename<-paste("/disco07/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresolLPA2shunt/events_sep",separations_str[ke],"sam_2000p_",primEnergies_str[ie],"keV_",secondEnergies_str[je],"keV_",detectionMode[idetectionMode],".fits",sep="")
            }else{
              if (detectionMode[idetectionMode] == "AD"){
                filename<-paste("/disco07/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresolLPA2shunt/events_sep",separations_str[ke],"sam_2000p_",primEnergies_str[ie],"keV_",secondEnergies_str[je],"keV_",detectionMode[idetectionMode],"_SHORT.fits",sep="")
              }else{
                filename<-paste("/disco07/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresolLPA2shunt/events_sep",separations_str[ke],"sam_2000p_",primEnergies_str[ie],"keV_",secondEnergies_str[je],"keV_",detectionMode[idetectionMode],".fits",sep="")
              }
            }
            #cat(filename,"\n")
            NAXIS2 <- readKeyword(filename,"NAXIS2")
            NETTOT <- readKeyword(filename,"NETTOT")
            totalPulses[ke] <- NETTOT   # NAXIS2 value from header
            #cat("NAXIS2: ",str(NAXIS2),"\n")
            #cat("NETTOT: ",str(NETTOT),"\n")
            detectedPulses[je,index] <- NAXIS2*100/NETTOT
            #cat("Pulsos=",str(NETTOT)," Detectados=",str(NAXIS2)," ",detectedPulses[je,index],"\n")
            index <- index-1
          }
        }
        
        mat_detectedPulses <- matrix(data=detectedPulses, nrow=length(secondEnergies), ncol=length(separations),byrow=F)
        save(mat_detectedPulses,file=matrixFile)
    }
    
    par(mar=c(5,5,5,9))
    # make image + axis + labels + title
    image(x=secondEnergies,y=separations,z=log10(mat_detectedPulses),log="y",ylim=c(3,2200),
                      zlim=c(log10(minimum),log10(maximum)),
                      xlab="Secondaries Energy [keV]",ylab="separations [samples]",
                      main=paste("Ep=",primEnergies_str[ie],"keV ",detectionMode[idetectionMode],sep=""),cex=1.5)
    # add to previous plot: image (in the background) + legend  w/ scale (axis.args)
    image.plot(x=secondEnergies,y=separations,z=log10(mat_detectedPulses),axes=FALSE,
                                axis.args=list(at=log10(ticks),labels=ticksLabels), log="y",
                                zlim=c(log10(minimum),log10(maximum)),col=pal(200),
                                xlab="", ylab="", add=TRUE,legend.width=5, legend.shrink=.8,cex=150)
    # Put labels in each pixel with the detected pulses
    for (je in 1:length(secondEnergies)){
        for (ke in length(separations_str):1){ 
            txtcol = "white"
            txtsize = 0.7
            if (mat_detectedPulses[je,ke] >maximum || mat_detectedPulses[je,ke] < minimum) txtcol="black"
            if(secondEnergies[je]<1) txtsize=0.5
            text(secondEnergies[je],separations[ke],sprintf("%2.1f",mat_detectedPulses[je,ke]), col=txtcol,cex=txtsize)
        }
    }
    box()
  }
}
dev.off()
