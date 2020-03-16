#
# Plot PH vs arrival times for 300  pulses
#
library(FITSio)

dataPH <- read.csv("xifusimLPA75um/pulses_7keV_jitter_nonoise_bbfb_sam998sam1500.txt", sep = "")
nPulses <- length(dataPH[,1])
eventsFile <- "eresolLPA75um/events_sep40000sam_1000p_SIRENA8192_pL8192_7keV_STC_F0F_fixedlib6OF_OPTFILT8192_jitter_nonoise_bbfb_FT3.fits"
zz <- file(description = eventsFile, open = "rb")
header0 <- readFITSheader(zz, fixHdr = 'remove') # read primary header
header <- readFITSheader(zz, fixHdr = 'remove') # read extension header
evtTable <- readFITSbintable(zz, header)
close(zz)
idcol <- which(evtTable$colNames == "PHI")
dataPHI <- numeric(nPulses)
dataPHI <- evtTable$col[[idcol]][1:nPulses]
colors<-rainbow(30)
#plot(dataPHI,dataPH[,1], xlab="Arrival Time - Trigger time", 
#     ylab="TES signal (arbitrary units)", ylim=c(1750, 4000), pch=19,col=colors[1])
#for (is in 2:3){
#    points(dataPHI,dataPH[,is],pch=19,col=colors[is])
#}
plot(998:1001, as.numeric(dataPH[1,1:4]),xlab="Sample", 
     ylab="TES signal (arbitrary units)", 
     ylim=c(1750, 6500), pch=1,col=colors[1],cex=0.4)
for (ip in 1:30){
    lines(998:1001, as.numeric(dataPH[ip,1:4]),pch=1,col=colors[ip], cex=0.4)
    points(998:1001, as.numeric(dataPH[ip,1:4]),pch=1,col=colors[ip], cex=0.4)
    #points((1000+dataPHI[ip]), 6000, pch=8, col=colors[ip])
}
#hist(dataPH[,3])
