#
# Plot values of simulated stream in sample 1000 (to check start of pulses for different energies)
#

# .txt file extracted from fv of fits file --> colum ADC --> expand --> 
#  Export as text --> column=1000 (all rows) --->User defined separator (blank)
#
#

array <- "LPA2shunt"
BASELINE <- 706.77
NOISESTD <- 7.28
pulseLength <- "4096"
files <- "analysis" #"calib" # "analysis"
#files <- "calib"

separations <- c("40000", "04605", "03433", "02560", "01908", "01423", "01061", "00791", "00589", "00439",
                 "00328", "00244", "00182", "00136", "00101", "00075", "00042", "00031", "00023", "00017",
                 "00013", "00010", "00005", "00002", "00001")
separations <- c("40000", "04605", "03433", "02560", "01908", "01423")
separations <- c("40000")

#sampleStr <- "1000sec" # to plot secondaries
sampleStr <- "1000prim" # to plot primaries


EkeV <- c(0.2,0.5,1,2,3,4,5,6,7,8)
#EkeV <- c(0.5)
EkeVStr = "ALL" #0.5
colors <- rainbow(length(EkeV))
if(length(EkeV) == 1) colors=c("red")

nSimPulses <- "200000"
nSimPulses <- "20000"
if (files == "calib"){
    dir <- paste("/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/tessim",
             array,"/",sep="")
    pdf(file=paste(dir,"pulsesCalib",nSimPulses,"AtSample",sample,".pdf", sep=""),width=10,height=7.)
}else if(files == "analysis"){
    dir <- paste("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/tessim", array,"/",sep="")
    pdf(file=paste(dir,"pulses",EkeVStr,"Simul",nSimPulses,"AtSample",sampleStr,".pdf",sep=""),width=10,height=7.)    
}
mat=matrix(data=1:27,nrow=3,ncol=3,byrow = T)
layout(mat)
for (ie in 1:length(EkeV)){
    col <- colors[ie]
    for (separation in separations){
        idxprim <- regexpr("p",sampleStr)[1]
        idxsec <- regexpr("s",sampleStr)[1]
        if(idxprim > 0){
            sampleNum <- as.numeric(substr(sampleStr,0,idxprim-1))
            cols <- c("999", "1000")  # to plot primaries            
        } else {
            if(separation == "40000") next
            # to plot secondaries
            sampleNum<-as.numeric(substr(sampleStr,0,idxsec-1)) + as.numeric(separation)
            col1 <- 999 + as.numeric(separation)
            col2 <- 999 + as.numeric(separation) +1
            cols <- list(col1, col2)
        }
        
        # ASCII file created with fv: Expand ADC, Export as Text (select all rows, column=998-1003),
        # User defined separator (blank)
        if (files == "calib"){
            file <- paste("mono",EkeV[ie],"_sep",separation,"_pix1_",nSimPulses,"p_",pulseLength,
                      "_init.txt",sep="")
        }else{
            file <- paste("sep",separation,"sam_",nSimPulses,"p_",EkeV[ie],"keV_cols998-1003.txt",sep="")
            file <- paste("sep",separation,"sam_",nSimPulses,"p_",EkeV[ie],"keV_sam",
                          cols[[1]],"_",cols[[-1]],".txt",sep="")
            
        }
        fileInit <- paste(dir,file,sep="")
        valuesTab <- read.table(fileInit, header=FALSE)
        values <- valuesTab[,which(cols==sampleNum)]
        cat("For ",file,":","\n")
        cat("     Max val in sample ",sampleNum,"=",max(values),"row(",which.max(values),")\n") 
        cat("     Min val in sample ",sampleNum,"=",min(values),"row(",which.min(values),")\n") 
        plot(seq(1:length(values)),values, main=paste("Energy = ",EkeV[ie]," keV\n sep=",separation,sep=""),
             pch=1,xlab="Pulse index",ylab=paste("Value of stream @ sample=",sampleNum,sep=""), col=colors[ie],
             ylim=c(min(675,min(values)),max(values)),cex=0.5)
        abline(h=BASELINE,lty=2,col="gray")
        abline(h=BASELINE+NOISESTD,lty=5,col="gray")
        abline(h=BASELINE-NOISESTD,lty=5,col="gray")
    }
}
dev.off()

