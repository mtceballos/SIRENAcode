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
separation <- "40000"
files <- "analysis" #"calib" # "analysis"
#files <- "calib"
sample <- "1002" # sample values to plot (998-1003)
cols <- list("998"=1, "999"=2, "1000"=3, "1001"=4, "1002"=5, "1003"=6)


EkeV <- c(0.2,0.5,1,2,3,4,5,6,7,8,9)
EkeV <- c(0.5)
nSimPulses <- "200000"
nSimPulses <- "20000"
if (files == "calib"){
    dir <- paste("/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/tessim",
             array,"/",sep="")
    pdf(file=paste(dir,"pulsesCalib",nSimPulses,"AtSample",sample,".pdf", sep=""),width=10,height=7.)
}else if(files == "analysis"){
    dir <- paste("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/tessim", array,"/",sep="")
    pdf(file=paste(dir,"pulsesSimul",nSimPulses,"AtSample",sample,".pdf",sep=""),width=10,height=7.)    
}
mat=matrix(data=1:12,nrow=3,ncol=4,byrow = T)
layout(mat)
colors <- rainbow(length(EkeV))
for (ie in 1:length(EkeV)){
    col <- colors[ie]
    # ASCII file created with fv: Expand ADC, Export as Text (select all rows, column=998-1003),
    # User defined separator (blank)
    if (files == "calib"){
        file <- paste("mono",EkeV[ie],"_sep",separation,"_pix1_",nSimPulses,"p_",pulseLength,
                  "_init.txt",sep="")
    }else{
        file <- paste("sep",separation,"sam_",nSimPulses,"p_",EkeV[ie],"keV_cols998-1003.txt",sep="")
    }
    fileInit <- paste(dir,file,sep="")
    valuesTab <- read.table(fileInit, header=FALSE)
    values <- valuesTab[,cols[[sample]]]
    cat("For ",file,":","\n")
    cat("     Max val in sample ",sample,"=",max(values),"row(",which.max(values),")\n") 
    cat("     Min val in sample ",sample,"=",min(values),"row(",which.min(values),")\n") 
    plot(seq(1:length(values)),values, main=paste("Energy = ",EkeV[ie]," keV",sep=""),
         pch=1,xlab="Pulse index",ylab="Value of stream @ sample=1000", col=colors[ie],
         ylim=c(min(650,min(values)),max(values)))
    abline(h=BASELINE,lty=2,col="gray")
    abline(h=BASELINE+NOISESTD,lty=5,col="gray")
    abline(h=BASELINE-NOISESTD,lty=5,col="gray")
}
dev.off()
