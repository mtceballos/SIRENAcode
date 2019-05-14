# Reads output of recon_resol.py (json and evt files) to obtain reconstructed
# energies vs input energies (gain scale curve).
# Get polynomial fit to gain scale curve:
#        a given array (SPA, LPA1, LPA2, LPA3)  
#        @ different energies (0.5,1,2,3,4,6,9) keV
#        for different methods (OPTFILT, I2R, WEIGHT, WEIGHTN, ...)

library(gridExtra)
library(Hmisc)
library(rjson)
library(FITSio)
library("RColorBrewer")

dcmt<- 100
dcmt<- 1
if (dcmt > 1){
    jitterStr=paste("_jitter_dcmt",dcmt,sep="")   
} else {
    jitterStr="_jitter"
}
#TRIGG = "STC"

# FUNCTIONS
#============
polyCurve <- function(x,coeffs) {
    polyres <- 0
    for (icoeff in 1:length(coeffs)){
        index <- icoeff - 1
        polyres <- polyres + coeffs[icoeff] * x^(index)
    }
    #a3[im]*x^3 + a2[im]*x^2 + a1[im]*x + a0[im]
    return(polyres)
}

# DEFINE and SAVE methods characteristics
#    Read them back with:
#    > load("/home/ceballos/INSTRUMEN/EURECA/ERESOL/methodsForR.Rdat")
#=======================================================================================
# samprate = 156250 Hz
# OF
fixed6OFsmprtAD <-   
    list(name="AD_F0F_fixedlib6OF_OPTFILT8192_jitter", nSamples=8192,
        samprateStr="", jitterStr="_jitter", detMethod="AD", ofLength=8192,
        lib="fixedlib6OF_OPTFILT",color="blue", point=1, ltype=1,
        lab=paste("OF_ADC (6keV, AD, s1",  jitterStr,")",sep=""))
fixed6OFsmprtSTC <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT8192",jitterStr,sep=""), 
         nSamples=8192, samprateStr="", jitterStr=jitterStr, detMethod="STC",
         noiseStr="", bbfbStr="", lib="fixedlib6OF_OPTFILT",  ofLength=8192,
         color="blue", point=0, ltype=2,
         lab=paste("OF_ADC (6keV, STC, s1", jitterStr,")",sep=""))
# OF NN
fixed6OFsmprtSTCnn <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT8192",jitterStr,"_nonoise",sep=""),
         nSamples=8192,samprateStr="", jitterStr=jitterStr, noiseStr="_nonoise",
         bbfbStr="", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color="blue", point=4, ltype=2,
         lab=paste("OF_ADC (6keV, AD, s1",jitterStr,", nonoise)",sep=""))
# OF BBFB
adccols = rev(brewer.pal(n = 9, name = "YlGnBu"))
adccols = rev(brewer.pal(n = 9, name = "Blues"))
#adccols = rainbow(6)
adccols=c("blue4", "blue", "cornflowerblue","cadetblue","blueviolet","darkmagenta")
fixed6OF8192smprtSTCBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT8192_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adccols[1], point=1, ltype=1,
         lab="OF_ADC (8192,6keV, STC, s1, bbfb)")
fixed6OF8192smprtSTCBbfb0.7Lc <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT8192_jitter_bbfb_0.7Lc",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="", LcStr="_0.7Lc",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adccols[1], point=1, ltype=1,
         lab="OF_ADC (8192,6keV, STC, s1, bbfb, 0.7Lc)")
fixed6OF8192smprtSTCBbfb0.5Lc <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT8192_jitter_bbfb_0.5Lc",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="", LcStr="_0.5Lc",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adccols[1], point=1, ltype=1,
         lab="OF_ADC (8192,6keV, STC, s1, bbfb, 0.5Lc)")
fixed6OF8192smprtSTCBbfb0.35Lc <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT8192_jitter_bbfb_0.35Lc",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="", LcStr="_0.35Lc",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adccols[1], point=1, ltype=1,
         lab="OF_ADC (8192,6keV, STC, s1, bbfb, 0.35Lc)")

fixed6OF4096smprtSTCBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT4096_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=4096,
         color=adccols[2], point=2, ltype=1,
         lab="OF_ADC (4096,6keV, STC, s1, bbfb)")
fixed6OF1024smprtSTCBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT1024_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=1024,
         color=adccols[3], point=3, ltype=1,
         lab="OF_ADC (1024,6keV, STC, s1, bbfb)")
fixed6OF512smprtSTCBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT512_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=512,
         color=adccols[4], point=4, ltype=1,
         lab="OF_ADC (512,6keV, STC, s1, bbfb)")
fixed6OF256smprtSTCBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT256_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=256,
         color=adccols[5], point=5, ltype=1,
         lab="OF_ADC (256,6keV, STC, s1, bbfb)")
fixed6OF128smprtSTCBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT128_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=128,
         color=adccols[6], point=6, ltype=1,
         lab="OF_ADC (128,6keV, STC, s1, bbfb)")

# OFNM BBFB
fixed6OF8192NM50000smprtSTCBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT8192NM50000_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color="cornflowerblue", point=2, ltype=0,
         lab="OF_ADC_NM(8192,50000int,6keV,STC,s1,bbfb)")
fixed6OF4096NM50000smprtSTCBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT4096NM50000_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=4096,
         color=adccols[2], point=17, ltype=0,
         lab="OF_ADC_NM(4096,50000int,6keV,STC,s1,bbfb)")

fixed6I2R8192NM150000smprtSTCBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_I2R8192NM150000_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color="red", point=2, ltype=0,
         lab="OF_R_NM(8192,150000int,6keV,STC,s1,bbfb)")
fixed6I2R4096NM150000smprtSTCBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_I2R4096NM150000_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=4096,
         color="orange", point=17, ltype=0,
         lab="OF_R_NM(4096,150000int,6keV,STC,s1,bbfb)")
fixed6I2R2048NM150000smprtSTCBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_I2R2048NM150000_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=2048,
         color="orange", point=2, ltype=0,
         lab="OF_R_NM(2048,150000int,6keV,STC,s1,bbfb)")

# OF BBFB NN
fixed6OF8192smprtSTCnnBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT8192_jitter_nonoise_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="_nonoise", 
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color="blue", point=19, ltype=1,
         lab="OF_ADC (8192,6keV, STC, s1, nonoise, bbfb)")
# I2R BBFB
fixed6I2R8192smprtSTCBbfb <-
    list(name="STC_F0F_fixedlib6OF_I2R8192_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=8192,
         color="orange", point=0, ltype=1,
         lab="OF_R (8192,6keV, STC, s1, bbfb)")
# I2RNOL BBFB
fixed6I2RNOL8192smprtSTCBbfb <-
    list(name="STC_F0F_fixedlib6OF_I2RNOL8192_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2RNOL", ofLength=8192,
         color="brown1", point=0, ltype=1,
         lab="OF_RNOL (6keV, STC, s1, bbfb)")

# samprate2 = 78125 Hz
fixed6OF4096smprt2STCnnBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT4096_samprate2_jitter_nonoise_bbfb",sep=""),
         nSamples=4096, samprateStr="_samprate2", jitterStr="_jitter", noiseStr="_nonoise",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=4096,
         color="blue", point=1, ltype=2,
         lab="OF_ADC (4096,6keV, STC, s2, nonoise, bbfb)")
fixed6OF4096smprt2STCBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT4096_samprate2_jitter_bbfb",sep=""),
         nSamples=4096, samprateStr="_samprate2", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=4096,
         color="blue", point=1, ltype=2,
         lab="OF_ADC (4096,6keV, STC, s2, bbfb)")






fixed6OFsmprt2AD <- 
    list(name="AD_F0F_fixedlib6OF_OPTFILT4096_samprate2_jitter", nSamples=4096, 
         samprateStr="_samprate2", jitterStr="_jitter",detMethod="AD",
         lib="fixedlib6OF_OPTFILT",color="cyan", point=1, ltype=1, ofLength=4096,
         lab=paste("OF_ADC (6keV, AD, s2", jitterStr,")",sep=""))
fixed6OFsmprt2STC <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT4096_samprate2", jitterStr,sep=""),
         nSamples=4096, samprateStr="_samprate2", jitterStr=jitterStr, detMethod="STC",
         noiseStr="",bbfbStr="", lib="fixedlib6OF_OPTFILT", ofLength=4096,
         color="blue", point=0, ltype=2,
         lab=paste("OF_ADC (6keV, STC, s2",jitterStr,sep=""))
fixed6OFsmprt2STCnn <-
        list(name=paste("STC_F0F_fixedlib6OF_OPTFILT4096_samprate2",jitterStr,"_nonoise",sep=""),
         nSamples=4096, samprateStr="_samprate2",jitterStr=jitterStr, noiseStr="_nonoise",
         bbfbStr="",detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=4096,
         color="blue", point=4, ltype=2,
         lab=paste("OF_ADC (6keV, AD, s2", jitterStr,", nonoise)",sep=""))

# samprate4: 39062.5 Hz
fixed6OFsmprt4STCBbfb <-
    list(name=paste("STC_F0F_fixedlib6OF_OPTFILT2048_samprate4_jitter_bbfb",sep=""),
         nSamples=2048, samprateStr="_samprate4", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", pLength=2048,
         color="blue", point=1, ltype=1,
         lab="OF_ADC (6keV, STC, s4, bbfb")
fixed6OFsmprt4STCnnBbfb <-
    list(name="STC_F0F_fixedlib6OF_OPTFILT2048_samprate4_jitter_bbfb",
         nSamples=2048, samprateStr="_samprate4",jitterStr="_jitter", noiseStr="_nonoise",
         bbfbStr="_bbfb",detMethod="STC", lib="fixedlib6OF_OPTFILT", pLength=2048,
         color="blue", point=4, ltype=2,
         lab="OF (6keV, STC, s4, nonoise, bbfb")

# 
# weight          <- list(name="multilib_WEIGHT",     color="darkgreen", point=15, ltype=1,
#                            lab="COVARIANCE MATRICES (Fixen)")
# weightn         <- list(name="multilib_WEIGHTN", color="darkgreen", point=0, ltype=3,
#                            lab="COVARIANCE MATRICES (0(n))")
# weightnOF       <- list(name="multilibOF_WEIGHTN", color="darkgreen", point=1, ltype=2,
#                         lab="COVARIANCE MATRICES OF (0(n))")
# 
save(fixed6OFsmprtSTCnn,
     fixed6OF8192smprtSTCBbfb,fixed6OF8192smprtSTCnnBbfb,
     fixed6OF4096smprtSTCBbfb,fixed6OF1024smprtSTCBbfb,
     fixed6OF512smprtSTCBbfb,fixed6OF256smprtSTCBbfb,
     fixed6OF128smprtSTCBbfb,fixed6OF8192NM50000smprtSTCBbfb,fixed6OF4096NM50000smprtSTCBbfb,
     fixed6I2R8192smprtSTCBbfb,fixed6I2R8192NM150000smprtSTCBbfb,fixed6I2R4096NM150000smprtSTCBbfb,
     fixed6I2RNOL8192smprtSTCBbfb,fixed6I2R2048NM150000smprtSTCBbfb,
     fixed6OFsmprt2STCnn,fixed6OFsmprt2STC,
     fixed6OF4096smprt2STCBbfb,fixed6OF4096smprt2STCnnBbfb,
     fixed6OFsmprt4STCBbfb,fixed6OFsmprt4STCnnBbfb,
     fixed6OF8192smprtSTCBbfb0.7Lc, fixed6OF8192smprtSTCBbfb0.5Lc,
     fixed6OF8192smprtSTCBbfb0.35Lc,
     file="/home/ceballos/INSTRUMEN/EURECA/ERESOL/methodsForR.Rdat")
#===========================================================================================

methods <- list(fixed6OF8192smprtSTCBbfb, fixed6OF4096smprtSTCBbfb,
                fixed6OF1024smprtSTCBbfb, fixed6OF512smprtSTCBbfb,
                fixed6OF256smprtSTCBbfb, fixed6OF128smprtSTCBbfb,
                fixed6OF8192NM50000smprtSTCBbfb,
                fixed6I2R8192smprtSTCBbfb, fixed6I2R8192NM150000smprtSTCBbfb
                )
methods <- list(fixed6OF8192smprtSTCBbfb, fixed6OF4096smprtSTCBbfb,
                fixed6OF1024smprtSTCBbfb, fixed6OF512smprtSTCBbfb,
                fixed6OF256smprtSTCBbfb, fixed6OF128smprtSTCBbfb,
                fixed6I2R8192smprtSTCBbfb, fixed6I2RNOL8192smprtSTCBbfb,
                fixed6OF4096NM50000smprtSTCBbfb,fixed6I2R4096NM150000smprtSTCBbfb,
                fixed6I2R2048NM150000smprtSTCBbfb,
                fixed6OF8192smprtSTCnnBbfb,
                fixed6OF4096smprt2STCnnBbfb,fixed6OF4096smprt2STCBbfb,
                fixed6OF8192smprtSTCBbfb0.7Lc,fixed6OF8192smprtSTCBbfb0.5Lc,
                fixed6OF8192smprtSTCBbfb0.35Lc
)
#methods <- list(fixed6OF8192smprtSTCBbfb,
#                fixed6I2R8192smprtSTCBbfb,
#                fixed6I2RNOL8192smprtSTCBbfb,
#                fixed6I2R8192NM150000smprtSTCBbfb)

nmethods <- length(methods)

# FWHM vs Energy Gain scale plot
#--------------------------------
npolyInit <- 4 # degree of polynomial to be fitted
pulsesCateg <- c("all")
#
# Initialize variables
#-----------------------
#array <- "LPA2shunt"
array <- "LPA75um" #"LPA2shunt"
nSimPulses <- "5000" # "20000"
separation <- 40000 #for samprate (corrected below for samprate2)

EkeV <- c(0.2,0.5,1,2,3,4,5,6,7,8)
EkeV <- c(1,2,3,4,5,6,7,8)
#nIntervals <- 50000
nIntervals <- 0
noiseMat<-""
if(nIntervals >0) noiseMat <-paste("_noiseMat",nIntervals,"i",sep="")
setwd(paste("/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresol",array,sep=""))
pdf(paste("./PDFs/polyfit2Bias_",nSimPulses,"p",noiseMat,".pdf",sep=""),width=10, height=7)
coeffsFile <- paste("coeffs_polyfit",noiseMat,".dat",sep="")
fwhmUNCORR  <- array(data=NA,dim=c(length(EkeV),nmethods)) # FWHM of Erecons

# READ CALIBRATION DATA
# ======================
meanEkeVfilt <- array(data=NA,dim=c(length(EkeV),nmethods))
errmean      <- array(data=NA,dim=c(length(EkeV),nmethods))
for (ie in 1:length(EkeV)){
    for (im in 1:nmethods){
        nSamples <- methods[[im]]$nSamples
        pulseLength <- methods[[im]]$nSamples
        ofLength <- methods[[im]]$ofLength
        TRIGG <- methods[[im]]$detMethod
        samprateStr <-methods[[im]]$samprateStr
        jitterStr <- methods[[im]]$jitterStr
        noiseStr <- methods[[im]]$noiseStr
        bbfbStr <- methods[[im]]$bbfbStr
        lib <- methods[[im]]$lib
        separation <- 40000
        if(samprateStr == "_samprate2") separation <- 20000
        if(samprateStr == "_samprate4") separation <- 10000
        #eresolFile <- paste("gainScale/eresol_",nSimPulses,"p_SIRENA",nSamples,
        #                    "_pL",pulseLength,"_", EkeV[ie],"keV_",TRIGG,"_F0F_", 
        #                    lib,ofLength,samprateStr,jitterStr,noiseStr,bbfbStr,".json",sep="")
        eresolFile <- paste("gainScale/eresol_",nSimPulses,"p_SIRENA",nSamples,
                    "_pL",pulseLength,"_", EkeV[ie],"keV_",methods[[im]]$name,".json",sep="")
        #eventsFile <- paste("gainScale/events_sep",separation,"sam_",nSimPulses,"p_SIRENA",nSamples,
        #                    "_pL",pulseLength,"_", EkeV[ie],"keV_",TRIGG,"_F0F_", 
        #                    lib,ofLength,samprateStr,jitterStr,noiseStr,bbfbStr,"_HR.fits",sep="")
        eventsFile <- paste("gainScale/events_sep",separation,"sam_",nSimPulses,"p_SIRENA",
                            nSamples,"_pL",pulseLength,"_", EkeV[ie],"keV_",
                            methods[[im]]$name,"_HR.fits",sep="")
        if(file.exists(eresolFile)){
            #data <- read.table(eresolFile,header=TRUE)
            # use data for selected separation (see initial definitions)
            cat("Reading file ",eresolFile,"\n")
            #cat("Separation= ",separation,"\n")
            #cat("   ie,im=",ie,im,"\n")
            jsondata <- fromJSON(file=eresolFile)
            idxSep <- which(sapply(jsondata,function(x) x$separation)==separation)
            #cat("   idxSep=",idxSep,"\n")
            fwhmUNCORR[ie,im]  <- as.numeric(jsondata[[idxSep]]$fwhmErecons[[pulsesCateg]])
            
        }else{
            warning("Not-existing file:", eresolFile)
            fwhmUNCORR[ie,im]  <- NaN
        }
        cat("Reading file ",eventsFile,"\n")
        stopifnot(file.exists(eventsFile))
        zz <- file(description = eventsFile, open = "rb")
        header0 <- readFITSheader(zz, fixHdr = 'remove') # read primary header
        header <- readFITSheader(zz, fixHdr = 'remove') # read extension header
        evtTable <- readFITSbintable(zz, header)
        close(zz)
        idcol <- which(evtTable$colNames == "SIGNAL")
        nreconPulses <- length(evtTable$col[[idcol]])
        cat("npulses=",nreconPulses,"\n")
        EkeVrecons <- numeric(nreconPulses)
        EkeVrecons <- evtTable$col[[idcol]]
        meanEkeVfilt[ie,im] <- mean(EkeVrecons)
        errmean[ie,im] <- sd(EkeVrecons)/sqrt(nreconPulses)

        if(all(is.nan(fwhmUNCORR[ie,im]))){
            warning("Error in ",eresolFile,"\n","  Non numerical values in eresol files: check event files")
        }
    } # for each method
} # foreach energy
cat("Writing coeffs:",coeffsFile,"\n")
fwhm  <- fwhmUNCORR

# Define plot layout (colors, tables)
# ====================================
colors <- sapply(methods, function(x) x$color)
MethodsLab <-  sapply(methods, function(x) x$lab)
alias <- sapply(methods, function(x) x$name)
points <- sapply(methods, function(x) x$point)
ltypes <- sapply(methods, function(x) x$ltype)

tt3 <- ttheme_minimal(
    core=list(bg_params = list(col=NA, fill="ivory"),
              fg_params=list(fontface=3,cex=0.6, col=colors)), 
    colhead=list(fg_params=list(col="navyblue")),
    rowhead=list(fg_params=list(col="navyblue", fontface=2L)))

# PLOT ENERGY_RECONS vs E input
#=================================
plot(EkeV,EkeV,type="n",mgp=c(2,1,0),ylim=c(0.,max(EkeV)+2),
     xlab="Input Energy keV (calibration points marked)", ylab="Reconstructed energy (keV)",
     main="GAIN SCALE (jitter)")
minor.tick(nx=5,ny=5,tick.ratio=0.5)
grid(nx=NA,ny=NULL)
lines(c(0,15),c(0,15),lty=2,col="gray")
for (e in EkeV){abline(v=EkeV,lty=3,col="grey80")}
coeffs <- matrix(0,nrow=nmethods,ncol=(npolyInit+1))
iemin <- c()
iemax <- c()
for (im in 1:nmethods){
    cat("Fitting method ",methods[[im]]$lab,"\n")
    # plot points and fit
    points(EkeV,meanEkeVfilt[,im],pch=methods[[im]]$point,col=methods[[im]]$color, 
           type="p",lty=1)
    #errbar(EkeV,meanEkeVfilt[,im], errbar.col = 'red', 
    #       yplus=meanEkeVfilt[,im]+errmean[,im],
    #       yminus=meanEkeVfilt[,im]-errmean[,im],
    #       pch=methods[[im]]$point, col=methods[[im]]$color)

    # POLYNOMIAL FIT
    #=================
    # check significance of npoly fit
    # use all calibration points 
    iemin[im] <- 1
    iemax[im] <- length(EkeV)
    badCoeffs <- min(npolyInit,(iemax[im]-iemin[im]))
    npoly <- min(npolyInit,(iemax[im]-iemin[im]))
    
    # first fit with orthogonal polynomia (assumes uncorrelated error)
    while (badCoeffs > 0){
        stopifnot(npoly > 0)
        fit <- lm(meanEkeVfilt[iemin[im]:iemax[im],im]~
                      poly(EkeV[iemin[im]:iemax[im]],npoly,raw=FALSE))
        probNrelev <- rep(-1,npoly)
        probNrelev <- as.numeric(summary(fit)$coefficients[2:(npoly+1),4])
        if (TRUE %in% is.na(probNrelev)){
            badCoeffs <- 1
        }else{
            badCoeffs <- length(probNrelev[probNrelev>0.001])  # not relevant coefficients
        }
        if (npoly== 1 && as.numeric(summary(fit)$coefficients[1,4])>0.001){ # bad Intercept
            stop("Error: bad intercept!")
        }else if(badCoeffs > 0){
            npoly <- npoly - badCoeffs
        }
    }
    # Once polynomia degree is set, fit with raw=TRUE
    coeffs[im,] <- 0
    fit <- lm(meanEkeVfilt[iemin[im]:iemax[im],im]~
                  poly(EkeV[iemin[im]:iemax[im]],npoly,raw=TRUE))
    coeffs[im,1:(npoly+1)] <- fit$coefficients[1:(npoly+1)]
    
    curve(polyCurve(x,coeffs[im,]), lty=methods[[im]]$ltype, 
          add=TRUE, col=methods[[im]]$color)
}

# PLOT polynomial fits over gain scale curve(s)
legend("topleft",legend=MethodsLab, col=colors, pch=points, lty=ltypes,
       cex=0.8,text.col=colors, bty="n")
text(7,1,paste("Lines= poly fits to calib energies",sep=""), cex=0.8)

# SAVE coefficients to coeffs table
coeffsTab <- data.frame(METHODS=MethodsLab, ALIAS=alias, coeffs)
cnames <- c("METHODS","ALIAS")
for (i in 0:npolyInit){
    colname <- paste("a",i,sep="")
    cnames <- append(cnames,colname)
}
colnames(coeffsTab) <- cnames
# save coefficients to be read by compareMethods.R
write.table(coeffsTab, file=coeffsFile, row.names=FALSE)

# PLOT table with results
plot(0:1,0:1, xlab="",ylab="", pch=NA_integer_, axes=FALSE)
mtext(bquote(paste(E[filt], " (keV)= a0 + a1 * ",E[cal]," + a2 * ", 
                   E[cal]^2, " + a3 * ",E[cal]^3, "+...", sep="")), side=3,cex=1.3)
grid.table(coeffsTab, rows=NULL, theme=tt3)
#legend("topleft", legend=fitLabel, title="Polynomial fit to Erecons vs Einput relation")

setwd("/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/")
dev.off()





