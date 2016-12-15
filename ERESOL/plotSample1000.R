#
# Plot values of simulated stream in sample 1000 (to check start of pulses for different energies)
#

# .txt file extracted form fv of fits file --> colum SIGNAL --> expand --> 
#  Save as text --> column=1000 (all rows)
#
#

###### CALIBRATION FILES ######################################

#EkeVcalib <- c(0.1,0.2,0.5,1,2,3,4,5,6,7,8,9,10,11,12)
#pulseSample1 <- c(695,739,856,1110,1520,1960,2400,2780,3250,3679,4100,4530,4950,5370,5780 )
EkeVcalib <- c(0.2,0.5,1,2,3,4,5,6,7,8,9)
pulseSample1 <- c(739,856,1110,1520,1960,2400,2780,3250,3679,4100,4530 )
nSimPulses <- "200000"
dir <- "/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/tessimLPA1shunt/"
pdf(file=paste(dir,"pulsesCalib",nSimPulses,"AtSample1000.pdf", sep=""),width=10,height=7.)
mat=matrix(data=1:12,nrow=3,ncol=4,byrow = T)
layout(mat)
colors <- rainbow(length(EkeVcalib))
for (ie in 1:length(EkeVcalib)){
    col <- colors[ie]
    # ASCII file created with fv: Expand ADC, Export as Text (select all rows, column=1000),
    # User defined separator (blank)
    file <- paste("mono",EkeVcalib[ie],"_sep20000_pix1_",nSimPulses,"p_2048_col1000.txt",sep="")
    file1000 <- paste(dir,file,sep="")
    val1000 <- read.table(file1000, header=FALSE)
    print(max(val1000[,1]))
    #plot(seq(1:length(val1000[,1])),val1000[,1], main=paste("Energy = ",EkeVcalib[ie]," keV",sep=""),
    #     pch=1,xlab="Pulse index",ylab="Value of stream @ sample=1000", col=colors[ie],
    #     ylim=c(min(650,min(val1000[,1])),max(val1000[,1])))
    #abline(h=660,lty=2,col="gray")
    # level where the first (non baseline sample should be): if no present => it's well above baseline
    #abline(h=pulseSample1[ie],lty=2, col="red") 
}
dev.off()


###### SIMULATION FILES ########################################

EkeVsimul <- c(0.1,0.2,0.5,1,2,3,4,5,6,7,8,9,10,11,12)
pulseSample1 <- c(695,739,856,1110,1520,1960,2400,2780,3250,3679,4100,4530,4950,5370,5780 )
dir <- "/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/tessimLPA1shunt/"
nSimPulses <- "20000"
pdf(file=paste(dir,"pulsesSimul",nSimPulses,"AtSample1000.pdf",sep=""),width=10,height=7.)
mat=matrix(data=1:18,nrow=3,ncol=3,byrow = T)
layout(mat)
colors <- rainbow(length(EkeVsimul))
for (ie in 1:length(EkeVsimul)){
    col <- colors[ie]
    file <- paste("sep20000sam_",nSimPulses,"p_",EkeVsimul[ie],"keV_col1000.txt",sep="")
    file1000 <- paste(dir,file,sep="")
    val1000 <- read.table(file1000, header=FALSE)
    plot(seq(1:length(val1000[,1])),val1000[,1], main=paste("Energy = ",EkeVsimul[ie]," keV",sep=""),
         pch=1,xlab="Pulse index",ylab="Value of stream @ sample=1000", col=colors[ie],
         ylim=c(min(650,min(val1000[,1])),max(val1000[,1])))
    abline(h=660,lty=2,col="gray")
    abline(h=pulseSample1[ie],lty=2, col="red")
}
dev.off()
