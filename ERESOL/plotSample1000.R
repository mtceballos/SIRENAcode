#
# Plot values of simulated stream in sample 1000 (to check start of pulses for different energies)
#

# .txt file extracted from fv of fits file --> colum SIGNAL --> expand --> 
#  Export as text --> column=1000 (all rows) --->User defined separator (blank)
#
#

array <- "LPA2shunt"
BASELINE <- 706.77
NOISESTD <- 7.28
pulseLength <- "4096"
separation <- "40000"

###### CALIBRATION FILES ######################################

EkeVcalib <- c(0.2,0.5,1,2,3,4,5,6,7,8,9)
nSimPulses <- "200000"
nSimPulses <- "20000"
dir <- paste("/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/tessim",
             array,"/",sep="")
pdf(file=paste(dir,"pulsesCalib",nSimPulses,"AtSample1000.pdf", sep=""),width=10,height=7.)
mat=matrix(data=1:12,nrow=3,ncol=4,byrow = T)
layout(mat)
colors <- rainbow(length(EkeVcalib))
for (ie in 1:length(EkeVcalib)){
    col <- colors[ie]
    # ASCII file created with fv: Expand ADC, Export as Text (select all rows, column=1000),
    # User defined separator (blank)
    file <- paste("mono",EkeVcalib[ie],"_sep",separation,"_pix1_",nSimPulses,"p_",pulseLength,
                  "_col1000.txt",sep="")
    file1000 <- paste(dir,file,sep="")
    val1000 <- read.table(file1000, header=FALSE)
    cat("Max val for ",file,"=",max(val1000[,1]),"\n") # it should be the baseline
    plot(seq(1:length(val1000[,1])),val1000[,1], main=paste("Energy = ",EkeVcalib[ie]," keV",sep=""),
         pch=1,xlab="Pulse index",ylab="Value of stream @ sample=1000", col=colors[ie],
         ylim=c(min(650,min(val1000[,1])),max(val1000[,1])))
    abline(h=BASELINE,lty=2,col="gray")
    abline(h=BASELINE+NOISESTD,lty=5,col="gray")
    abline(h=BASELINE-NOISESTD,lty=5,col="gray")
    # level where the first (non baseline sample should be): if no present => it's well above baseline
    #abline(h=pulseSample1[ie],lty=2, col="red") 
}
dev.off()


###### SIMULATION FILES ########################################

EkeVsimul <- c(0.1,0.2,0.5,1,2,3,4,5,6,7,8)
nSimPulses <- "20000"
dir <- paste("/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/tessim",
             array,"sep=")
pdf(file=paste(dir,"pulsesSimul",nSimPulses,"AtSample1000.pdf",sep=""),width=10,height=7.)
mat=matrix(data=1:18,nrow=3,ncol=3,byrow = T)
layout(mat)
colors <- rainbow(length(EkeVsimul))
for (ie in 1:length(EkeVsimul)){
    col <- colors[ie]
    file <- paste("sep",separation,"sam_",nSimPulses,"p_",EkeVsimul[ie],"keV_col1000.txt",sep="")
    file1000 <- paste(dir,file,sep="")
    val1000 <- read.table(file1000, header=FALSE)
    plot(seq(1:length(val1000[,1])),val1000[,1], main=paste("Energy = ",EkeVsimul[ie]," keV",sep=""),
         pch=1,xlab="Pulse index",ylab="Value of stream @ sample=1000", col=colors[ie],
         ylim=c(min(650,min(val1000[,1])),max(val1000[,1])))
    abline(h=660,lty=2,col="gray")
    abline(h=pulseSample1[ie],lty=2, col="red")
}
dev.off()
