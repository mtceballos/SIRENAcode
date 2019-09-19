# Reads output of recon_resol.py (json and evt files) to obtain reconstructed
# energies vs input energies (gain scale curve).
# Get polynomial fit to gain scale curve:
#        a given array (SPA, LPA1, LPA2, LPA3)  
#        @ different energies (0.5,1,2,3,4,6,9) keV
#        for different methods (OPTFILT, I2R, WEIGHT, WEIGHTN, ...)
rm(list=ls())
library(gridExtra)
library(Hmisc)
library(rjson)
library(FITSio)
library("RColorBrewer")
source("~/R/Rfunctions/drawLogPlotBox.r")

#
# Initialize variables
#-----------------------
#array <- "LPA2shunt"
array <- "LPA75um" #"LPA2shunt"
nSimPulses <- "5000" # "20000"
separation <- 40000 #for samprate (corrected below for samprate2)
gainScaleID <-"methods_shortFilters"   # !!!! CHECK METHODS BELOW !!!!!
gainScaleID <-"methods_longFilter_zeroPadding"   # !!!! CHECK METHODS BELOW !!!!!
gainScaleID <-"methods_longFilter_zeroPadding_ADC_I2R"   # !!!! CHECK METHODS BELOW !!!!!
gainScaleID <-"methods_allFilter_ADC_I2R"   # !!!! CHECK METHODS BELOW !!!!!
EkeV <- c(0.2,0.5,1,2,3,4,5,6,7,8)
EkeV <- c(0.2,0.5,1,2,3,4,5,6,7)
#nIntervals <- 50000
nIntervals <- 0
noiseMat<-""
if(nIntervals >0) noiseMat <-paste("_noiseMat",nIntervals,"i",sep="")
setwd(paste("/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresol",array,sep=""))
#pdf(paste("./PDFs/polyfit2Bias_",nSimPulses,"p",noiseMat,".pdf",sep=""),width=10, height=7)
#coeffsFile <- paste("coeffs_polyfit",noiseMat,".dat",sep="")
pdf(paste("./PDFs/polyfit2Bias_",gainScaleID,".pdf",sep=""),width=10, height=7)
coeffsFile <- paste("coeffs_polyfit_",gainScaleID,".dat",sep="")
splineFile <- paste("spline_forfit_",gainScaleID,".json",sep="")
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
# METHODS NAME Changed from STC_F0F to STC_T_  (from Freq domain to T domain)
pBstr=""
fixed6OFsmprtSTC <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192",jitterStr,sep=""), 
         nSamples=8192, samprateStr="", jitterStr=jitterStr, detMethod="STC",
         noiseStr="", bbfbStr="", lib="fixedlib6OF_OPTFILT",  ofLength=8192,
         color="blue", point=0, ltype=2,
         lab=paste("OF_ADC (6keV, STC, s1", jitterStr,")",sep=""))
# OF NN
fixed6OFsmprtSTCnn <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192",jitterStr,"_nonoise",sep=""),
         nSamples=8192,samprateStr="", jitterStr=jitterStr, noiseStr="_nonoise",
         bbfbStr="", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color="blue", point=4, ltype=2,
         lab=paste("OF_ADC (6keV, AD, s1",jitterStr,", nonoise)",sep=""))
# OF BBFB
adccols = rev(brewer.pal(n = 9, name = "YlGnBu"))
adccols = rev(brewer.pal(n = 9, name = "Blues"))
adccols = rep("blue",9)
adc2cols= rev(brewer.pal(n = 9, name = "Reds"))
adc2cols = rep("red",9)
adc3cols= rev(brewer.pal(n = 9, name = "Greens"))
adc3cols = rep("green",9)
i2rcols = rev(brewer.pal(n = 9, name = "Oranges"))
#adccols = rainbow(6)
#adccols=c("blue4", "blue", "cornflowerblue","cadetblue","blueviolet","darkmagenta")

adc0cols = rep("blue4",9)
adcSHcols = rep("cornflowerblue",9)
adcPBcols = rep("blueviolet",9)
adcNMcols = rep("turquoise",9)

i2r0cols = rep("darkorange4",9)
i2rSHcols = rep("darkorange",9)
i2rPB50cols = rep("brown3",9)
i2rPB100cols = rep("gold4",9)

###### OPTFILT  ############

#--Lc
pL8192fixed6OF8192smprtSTCBbfb0.7Lc <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192_jitter_bbfb_0.7Lc",sep=""),pLength=8192,
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="", LcStr="_0.7Lc",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adccols[1], point=1, ltype=1,
         lab="OF_ADC (pL8192,ofL8192,6keV, STC, s1, bbfb, 0.7Lc)")
pL4096fixed6OF4096smprt2STCBbfb0.7Lc <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT4096_samprate2_jitter_bbfb_0.7Lc",sep=""),
         nSamples=4096, samprateStr="_samprate2", jitterStr="_jitter", noiseStr="", LcStr="_0.5Lc",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=4096,
         color=adccols[1], point=1, ltype=1,pLength=4096,
         lab="OF_ADC (pL4096,ofL4096,6keV, STC, s2, bbfb, 0.7Lc)")
pL8192fixed6OF8192smprtSTCBbfb0.5Lc <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192_jitter_bbfb_0.5Lc",sep=""),pLength=8192,
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="", LcStr="_0.5Lc",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adccols[1], point=1, ltype=1,
         lab="OF_ADC (pL8192,ofL8192,6keV, STC, s1, bbfb, 0.5Lc)")
pL4096fixed6OF4096smprt2STCBbfb0.5Lc <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT4096_samprate2_jitter_bbfb_0.5Lc",sep=""),
         nSamples=4096, samprateStr="_samprate2", jitterStr="_jitter", noiseStr="", LcStr="_0.5Lc",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=4096,pLength=4096,
         color=adccols[1], point=1, ltype=1,
         lab="OF_ADC (pL4096,ofL4096,6keV, STC, s2, bbfb, 0.5Lc)")
pL8192fixed6OF8192smprtSTCBbfb0.35Lc <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192_jitter_bbfb_0.35Lc",sep=""),pLength=8192,
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="", LcStr="_0.35Lc",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adccols[1], point=1, ltype=1,
         lab="OF_ADC (pL8192,ofL8192,6keV, STC, s1, bbfb, 0.35Lc)")
#--
pL8192fixed6OF8192smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adc0cols[1], point=1, ltype=1,
         lab="OF_ADC (pL8192,ofL8192,6keV, STC, s1, bbfb)")

pL4096fixed6OF8192smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=4096,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adc0cols[2], point=2, ltype=1,
         lab="OF_ADC (pL4096,ofL8192,6keV, STC, s1, bbfb)")
pL8192fixed6OF4096smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT4096_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=4096,
         color=adcSHcols[2], point=2, ltype=1,
         lab="OF_ADC (pL8192,ofL4096,6keV, STC, s1, bbfb)")
pL8192fixed6OF4096smprtSTCBbfb_pB50 <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT4096_pB50_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=4096,
         color=adcPBcols[2], point=2, ltype=1, 
         lab="OF_ADC (pL8192,ofL4096,6keV, STC, s1, bbfb, pB50)")

pL2048fixed6OF8192smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=2048,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adc0cols[3], point=3, ltype=1,
         lab="OF_ADC (pL2048,ofL8192,6keV, STC, s1, bbfb)")
pL8192fixed6OF2048smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT2048_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=2048,
         color=adcSHcols[3], point=3, ltype=1,
         lab="OF_ADC (pL8192,ofL2048,6keV, STC, s1, bbfb)")
pL8192fixed6OF2048smprtSTCBbfb_pB50 <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT2048_pB50_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=2048,
         color=adcPBcols[3], point=3, ltype=1,
         lab="OF_ADC (pL8192,ofL2048,6keV, STC, s1, bbfb, pB50)")

pL1024fixed6OF8192smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=1024,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adc0cols[4], point=4, ltype=1,
         lab="OF_ADC (pL1024,ofL8192,6keV, STC, s1, bbfb)")
pL8192fixed6OF1024smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT1024_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=1024,
         color=adcSHcols[4], point=4, ltype=1,
         lab="OF_ADC (pL8192,ofL1024,6keV, STC, s1, bbfb)")
pL8192fixed6OF1024smprtSTCBbfb_pB50 <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT1024_pB50_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=1024,
         color=adcPBcols[4], point=4, ltype=1, 
         lab="OF_ADC (pL8192,ofL1024,6keV, STC, s1, bbfb, pB50)")

pL512fixed6OF8192smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=512,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adc0cols[5], point=5, ltype=1,
         lab="OF_ADC (pL512,ofL8192,6keV, STC, s1, bbfb)")
pL512fixed6OF8192smprtSTCBbfbSUM0 <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192_jitter_bbfb_Sum0Filt",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=512,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color="black", point=5, ltype=1, Sum0Filt=1,
         lab="OF_ADC (pL512,ofL8192,6keV, STC, s1, bbfb, SUM0)")
pL8192fixed6OF512smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT512_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=512,
         color=adcSHcols[5], point=5, ltype=1,
         lab="OF_ADC (pL8192,ofL512,6keV, STC, s1, bbfb)")
pL8192fixed6OF512smprtSTCBbfb_pB50 <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT512_pB50_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=512,
         color=adcPBcols[5], point=5, ltype=1, 
         lab="OF_ADC (pL8192,ofL512,6keV, STC, s1, bbfb, pB50)")

pL256fixed6OF8192smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=256,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adc0cols[6], point=6, ltype=1,
         lab="OF_ADC (pL256,ofL8192,6keV, STC, s1, bbfb)")
pL8192fixed6OF256smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT256_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=256,
         color=adcSHcols[6], point=6, ltype=1,
         lab="OF_ADC (pL8192,ofL256,6keV, STC, s1, bbfb)")
pL8192fixed6OF256smprtSTCBbfb_pB50 <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT256_pB50_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=256,
         color=adcPBcols[6], point=6, ltype=1, 
         lab="OF_ADC (pL8192,ofL256,6keV, STC, s1, bbfb, pB50)")

pL128fixed6OF8192smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=128,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adc0cols[7], point=7, ltype=1,
         lab="OF_ADC (pL128,ofL8192,6keV, STC, s1, bbfb)")
pL8192fixed6OF128smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT128_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=128,
         color=adcSHcols[7], point=7, ltype=1,
         lab="OF_ADC (pL8192,ofL128,6keV, STC, s1, bbfb)")
pL8192fixed6OF128smprtSTCBbfb_pB50 <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT128_pB50_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=128,
         color=adcPBcols[7], point=7, ltype=1, 
         lab="OF_ADC (pL8192,ofL128,6keV, STC, s1, bbfb, pB50)")

##### I2R BBFB ########
pL8192fixed6I2R8192smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2R8192_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=8192,
         color=i2r0cols[1], point=1, ltype=1,
         lab="OF_R (pL8192,ofL8192,6keV, STC, s1, bbfb)")

pL4096fixed6I2R8192smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2R8192_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=4096,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=8192,
         color=i2r0cols[2], point=2, ltype=1,
         lab="OF_R (pL4096,ofL8192,6keV, STC, s1, bbfb)")
pL8192fixed6I2R4096smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2R4096_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=4096,
         color=i2rSHcols[2], point=2, ltype=1,
         lab="OF_R (pL8192,ofL4096,6keV, STC, s1, bbfb)")
pL8192fixed6I2R4096smprtSTCBbfb_pB50 <-
    list(name="STC_T_fixedlib6OF_I2R4096_pB50_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=4096,
         color=i2rPB50cols[2], point=2, ltype=1,
         lab="OF_R (pL8192,ofL4096,6keV, STC, s1, bbfb, pB50)")
pL8192fixed6I2R4096smprtSTCBbfb_pB100 <-
    list(name="STC_T_fixedlib6OF_I2R4096_pB100_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=4096,
         color=i2rPB100cols[2], point=2, ltype=2,
         lab="OF_R (pL8192,ofL4096,6keV, STC, s1, bbfb, pB100)")

pL2048fixed6I2R8192smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2R8192_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=2048,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=8192,
         color=i2r0cols[3], point=3, ltype=1,
         lab="OF_R (pL2048,ofL8192,6keV, STC, s1, bbfb)")
pL8192fixed6I2R2048smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2R2048_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=2048,
         color=i2rSHcols[3], point=3, ltype=1,
         lab="OF_R (pL8192,ofL2048,6keV, STC, s1, bbfb)")
pL8192fixed6I2R2048smprtSTCBbfb_pB50 <-
    list(name="STC_T_fixedlib6OF_I2R2048_pB50_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=2048,
         color=i2rPB50cols[3], point=3, ltype=1,
         lab="OF_R (pL8192,ofL2048,6keV, STC, s1, bbfb, pB50)")
pL8192fixed6I2R2048smprtSTCBbfb_pB100 <-
    list(name="STC_T_fixedlib6OF_I2R2048_pB100_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=2048,
         color=i2rPB100cols[3], point=3, ltype=2,
         lab="OF_R (pL8192,ofL2048,6keV, STC, s1, bbfb, pB100)")

pL1024fixed6I2R8192smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2R8192_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=1024,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=8192,
         color=i2r0cols[4], point=4, ltype=1,
         lab="OF_R (pL1024,ofL8192,6keV, STC, s1, bbfb)")
pL8192fixed6I2R1024smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2R1024_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=1024,
         color=i2rSHcols[4], point=4, ltype=1,
         lab="OF_R (pL8192,ofL1024,6keV, STC, s1, bbfb)")
pL8192fixed6I2R1024smprtSTCBbfb_pB50 <-
    list(name="STC_T_fixedlib6OF_I2R1024_pB50_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=1024,
         color=i2rPB50cols[4], point=4, ltype=1,
         lab="OF_R (pL8192,ofL1024,6keV, STC, s1, bbfb, pB50)")
pL8192fixed6I2R1024smprtSTCBbfb_pB100 <-
    list(name="STC_T_fixedlib6OF_I2R1024_pB100_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=1024,
         color=i2rPB100cols[4], point=4, ltype=2,
         lab="OF_R (pL8192,ofL1024,6keV, STC, s1, bbfb, pB100)")

pL512fixed6I2R8192smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2R8192_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=512,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=8192,
         color=i2r0cols[5], point=5, ltype=1,
         lab="OF_R (pL512,ofL8192,6keV, STC, s1, bbfb)")
pL8192fixed6I2R512smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2R512_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=512,
         color=i2rSHcols[5], point=5, ltype=1,
         lab="OF_R (pL8192,ofL512,6keV, STC, s1, bbfb)")
pL8192fixed6I2R512smprtSTCBbfb_pB50 <-
    list(name="STC_T_fixedlib6OF_I2R512_pB50_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=512,
         color=i2rPB50cols[5], point=5, ltype=1,
         lab="OF_R (pL8192,ofL512,6keV, STC, s1, bbfb, pB50)")
pL8192fixed6I2R512smprtSTCBbfb_pB100 <-
    list(name="STC_T_fixedlib6OF_I2R512_pB100_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=512,
         color=i2rPB100cols[5], point=5, ltype=2,
         lab="OF_R (pL8192,ofL512,6keV, STC, s1, bbfb, pB100)")

pL256fixed6I2R8192smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2R8192_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=256,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=8192,
         color=i2r0cols[6], point=6, ltype=1,
         lab="OF_R (pL256,ofL8192,6keV, STC, s1, bbfb)")
pL8192fixed6I2R256smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2R256_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=256,
         color=i2rSHcols[6], point=6, ltype=1,
         lab="OF_R (pL8192,ofL256,6keV, STC, s1, bbfb)")
pL8192fixed6I2R256smprtSTCBbfb_pB50 <-
    list(name=paste("STC_T_fixedlib6OF_I2R256_pB50_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_I2R", ofLength=256,
         color=i2rPB50cols[6], point=6, ltype=1, 
         lab="OF_R (pL8192,ofL256,6keV, STC, s1, bbfb, pB50)")
pL8192fixed6I2R256smprtSTCBbfb_pB100 <-
    list(name=paste("STC_T_fixedlib6OF_I2R256_pB100_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_I2R", ofLength=256,
         color=i2rPB100cols[6], point=6, ltype=2, 
         lab="OF_R (pL8192,ofL256,6keV, STC, s1, bbfb, pB100)")

pL128fixed6I2R8192smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2R8192_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=128,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=8192,
         color=i2r0cols[7], point=7, ltype=1,
         lab="OF_R (pL128,ofL8192,6keV, STC, s1, bbfb)")
pL8192fixed6I2R128smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2R128_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2R", ofLength=128,
         color=i2rSHcols[7], point=7, ltype=1,
         lab="OF_R (pL8192,ofL128,6keV, STC, s1, bbfb)")
pL8192fixed6I2R128smprtSTCBbfb_pB50 <-
    list(name=paste("STC_T_fixedlib6OF_I2R128_pB50_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_I2R", ofLength=128,
         color=i2rPB50cols[7], point=7, ltype=1, 
         lab="OF_R (pL8192,ofL128,6keV, STC, s1, bbfb, pB50)")
pL8192fixed6I2R128smprtSTCBbfb_pB100 <-
    list(name=paste("STC_T_fixedlib6OF_I2R128_pB100_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_I2R", ofLength=128,
         color=i2rPB100cols[7], point=7, ltype=1, 
         lab="OF_R (pL8192,ofL128,6keV, STC, s1, bbfb, pB100)")


# OFNM BBFB
pL8192fixed6OF8192NM50000smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192NM50000_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adcNMcols[1], point=2, ltype=0,
         lab="OF_ADC_NM(pL8192,ofL8192,50000int,6keV,STC,s1,bbfb)")
pL4096fixed6OF8192NM50000smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192NM50000_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=4096,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=adcNMcols[2], point=17, ltype=0,
         lab="OF_ADC_NM(pL4096,ofL8192,50000int,6keV,STC,s1,bbfb)")
pL1024fixed6OF1024NM50000smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT1024NM50000_jitter_bbfb",sep=""),
         nSamples=1024, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=1024,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=1024,
         color=adcNMcols[3], point=4, ltype=0,
         lab="OF_ADC_NM(pL1024,ofL1024,50000int,6keV,STC,s1,bbfb)")

pL8192fixed6I2R8192NM150000smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_I2R8192NM150000_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=i2rcols[1], point=2, ltype=0,
         lab="OF_R_NM(pL8192,ofL8192,150000int,6keV,STC,s1,bbfb)")
pL4096fixed6I2R8192NM150000smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_I2R8192NM150000_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=4096,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=i2rcols[2], point=17, ltype=0,
         lab="OF_R_NM(pL4096,ofL8192,150000int,6keV,STC,s1,bbfb)")
pL2048fixed6I2R8192NM150000smprtSTCBbfb <-
    list(name=paste("STC_T_fixedlib6OF_I2R8192NM150000_jitter_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=2048,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color=i2rcols[3], point=6, ltype=0,
         lab="OF_R_NM(pL2048,ofL8192,150000int,6keV,STC,s1,bbfb)")

# OF BBFB NN
fixed6OF8192smprtSTCnnBbfb <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT8192_jitter_nonoise_bbfb",sep=""),
         nSamples=8192, samprateStr="", jitterStr="_jitter", noiseStr="_nonoise", pLength=8192,
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=8192,
         color="blue", point=19, ltype=1,
         lab="OF_ADC (8192,6keV, STC, s1, nonoise, bbfb)")
# I2RNOL BBFB
pL8192fixed6I2RNOL8192smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2RNOL8192_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2RNOL", ofLength=8192,
         color="red", point=1, ltype=1,
         lab="OF_RNOL (pL8192,ofL8192, 6keV, STC, s1, bbfb)")
# I2RFITTED BBFB
pL8192fixed6I2RFITTED8192smprtSTCBbfb <-
    list(name="STC_T_fixedlib6OF_I2RFITTED8192_jitter_bbfb", 
         nSamples=8192, samprateStr="", jitterStr="_jitter", detMethod="STC",pLength=8192,
         noiseStr="", bbfbStr="_bbfb", lib="fixedlib6OF_I2RFITTED", ofLength=8192,
         color="green", point=0, ltype=1,
         lab="OF_RFITTED (pL8192,ofL8192,6keV, STC, s1, bbfb)")

# WEIGHTNOF   
pL1024multilibWEIGHTN1024NM50000smprtSTCBbfb <-
    list(name=paste("STC_T_multilibOF_WEIGHTN1024_jitter_bbfb",sep=""),
         nSamples=1024, samprateStr="", jitterStr="_jitter", noiseStr="",pLength=1024,
         bbfbStr="_bbfb", detMethod="STC", lib="multilibOF_WEIGHTN", ofLength=1024,
         color="darkgreen", point=4, ltype=0,
         lab="OF_ADC_WEIGHTN(pL1024,ofL1024,50000int,STC,s1,bbfb)")

# samprate2 = 78125 Hz
pL4096fixed6OF4096smprt2STCnnBbfb <-
    list(name="STC_T_fixedlib6OF_OPTFILT4096_samprate2_jitter_nonoise_bbfb",pLength=4096,
         nSamples=4096, samprateStr="_samprate2", jitterStr="_jitter", noiseStr="_nonoise",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=4096,
         color="blue", point=1, ltype=2,
         lab="OF_ADC (pL4096,ofL4096,6keV, STC, s2, nonoise, bbfb)")
pL4096fixed6OF4096smprt2STCBbfb <-
    list(name="STC_T_fixedlib6OF_OPTFILT4096_samprate2_jitter_bbfb",pLength=4096,
         nSamples=4096, samprateStr="_samprate2", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=4096,
         color="blue", point=1, ltype=2,
         lab="OF_ADC (pL4096,ofL4096,6keV, STC, s2, bbfb)")


pL4096fixed6OF4096smprt2STC <-
    list(name=paste("STC_T_fixedlib6OF_OPTFILT4096_samprate2", jitterStr,sep=""),
         nSamples=4096, samprateStr="_samprate2", jitterStr=jitterStr, detMethod="STC",
         noiseStr="",bbfbStr="", lib="fixedlib6OF_OPTFILT", ofLength=4096,
         color="blue", point=0, ltype=2,pLength=4096,
         lab=paste("OF_ADC (pL4096,ofL4096,6keV, STC, s2",jitterStr,sep=""))
pL4096fixed6OF4096smprt2STCnn <-
        list(name=paste("STC_T_fixedlib6OF_OPTFILT4096_samprate2",jitterStr,"_nonoise",sep=""),
         nSamples=4096, samprateStr="_samprate2",jitterStr=jitterStr, noiseStr="_nonoise",
         bbfbStr="",detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=4096,pLength=4096,
         color="blue", point=4, ltype=2,
         lab=paste("OF_ADC (pL4096,ofL4096,6keV, AD, s2", jitterStr,", nonoise)",sep=""))

# samprate4: 39062.5 Hz
pL2048fixed6OF2048smprt4STCnnBbfb <-
    list(name="STC_T_fixedlib6OF_OPTFILT2048_samprate4_jitter_nonoise_bbfb",pLength=2048,
         nSamples=2048, samprateStr="_samprate4", jitterStr="_jitter", noiseStr="_nonoise",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=2048,
         color="blue", point=1, ltype=2,
         lab="OF_ADC (pL2048,ofL2048,6keV, STC, s4, nonoise, bbfb)")
pL2048fixed6OF2048smprt4STCBbfb <-
    list(name="STC_T_fixedlib6OF_OPTFILT2048_samprate4_jitter_bbfb",pLength=2048,
         nSamples=2048, samprateStr="_samprate4", jitterStr="_jitter", noiseStr="",
         bbfbStr="_bbfb", detMethod="STC", lib="fixedlib6OF_OPTFILT", ofLength=2048,
         color="blue", point=1, ltype=2,
         lab="OF_ADC (pL2048,ofL2048,6keV, STC, s4, bbfb)")


# 
# weight          <- list(name="multilib_WEIGHT",     color="darkgreen", point=15, ltype=1,
#                            lab="COVARIANCE MATRICES (Fixen)")
# weightn         <- list(name="multilib_WEIGHTN", color="darkgreen", point=0, ltype=3,
#                            lab="COVARIANCE MATRICES (0(n))")

# 
#save(fixed6OFsmprtSTCnn,
#     fixed6OF8192smprtSTCBbfb,fixed6OF8192smprtSTCnnBbfb,
#     fixed6OF4096smprtSTCBbfb,fixed6OF1024smprtSTCBbfb,
#     fixed6OF512smprtSTCBbfb,fixed6OF256smprtSTCBbfb,
#     fixed6OF128smprtSTCBbfb,fixed6OF8192NM50000smprtSTCBbfb,fixed6OF4096NM50000smprtSTCBbfb,
#     fixed6I2R8192smprtSTCBbfb,fixed6I2R8192NM150000smprtSTCBbfb,fixed6I2R4096NM150000smprtSTCBbfb,
#     fixed6I2RNOL8192smprtSTCBbfb,fixed6I2R2048NM150000smprtSTCBbfb,
#     fixed6I2RFITTED8192smprtSTCBbfb,
#     fixed6OFsmprt2STCnn,fixed6OFsmprt2STC,
#     fixed6OF4096smprt2STCBbfb,fixed6OF4096smprt2STCnnBbfb,
#     fixed6OF2048smprt4STCBbfb,fixed6OF2048smprt4STCnnBbfb,
#     fixed6OF8192smprtSTCBbfb0.7Lc, fixed6OF4096smprt2STCBbfb0.7Lc,
#     fixed6OF8192smprtSTCBbfb0.5Lc,fixed6OF4096smprt2STCBbfb0.5Lc,
#     fixed6OF8192smprtSTCBbfb0.35Lc,
#     file="/home/ceballos/INSTRUMEN/EURECA/ERESOL/methodsForR.Rdat")
save(
     pL8192fixed6OF8192smprtSTCBbfb, pL8192fixed6OF4096smprtSTCBbfb,
     pL8192fixed6OF2048smprtSTCBbfb, pL8192fixed6OF1024smprtSTCBbfb, 
     pL8192fixed6OF512smprtSTCBbfb,  pL8192fixed6OF256smprtSTCBbfb,  
     pL8192fixed6OF128smprtSTCBbfb,
     pL8192fixed6OF4096smprtSTCBbfb_pB50, pL8192fixed6OF2048smprtSTCBbfb_pB50,
     pL8192fixed6OF1024smprtSTCBbfb_pB50, pL8192fixed6OF512smprtSTCBbfb_pB50,
     pL8192fixed6OF256smprtSTCBbfb_pB50,  pL8192fixed6OF128smprtSTCBbfb_pB50,
                                     pL4096fixed6OF8192smprtSTCBbfb, 
     pL2048fixed6OF8192smprtSTCBbfb, pL1024fixed6OF8192smprtSTCBbfb, 
     pL512fixed6OF8192smprtSTCBbfb,  pL256fixed6OF8192smprtSTCBbfb,  
     pL128fixed6OF8192smprtSTCBbfb,
     pL8192fixed6I2R8192smprtSTCBbfb, pL4096fixed6I2R8192smprtSTCBbfb,
     pL2048fixed6I2R8192smprtSTCBbfb, pL1024fixed6I2R8192smprtSTCBbfb,
     pL512fixed6I2R8192smprtSTCBbfb,  pL256fixed6I2R8192smprtSTCBbfb,
     pL128fixed6I2R8192smprtSTCBbfb,
     pL8192fixed6I2R4096smprtSTCBbfb_pB50, pL8192fixed6I2R2048smprtSTCBbfb_pB50,
     pL8192fixed6I2R1024smprtSTCBbfb_pB50, pL8192fixed6I2R512smprtSTCBbfb_pB50,
     pL8192fixed6I2R256smprtSTCBbfb_pB50,  pL8192fixed6I2R128smprtSTCBbfb_pB50,
     pL8192fixed6I2R4096smprtSTCBbfb_pB100, pL8192fixed6I2R2048smprtSTCBbfb_pB100,
     pL8192fixed6I2R1024smprtSTCBbfb_pB100, pL8192fixed6I2R512smprtSTCBbfb_pB100,
     pL8192fixed6I2R256smprtSTCBbfb_pB100,  pL8192fixed6I2R128smprtSTCBbfb_pB100,
     pL8192fixed6I2R4096smprtSTCBbfb, pL8192fixed6I2R2048smprtSTCBbfb,
     pL8192fixed6I2R1024smprtSTCBbfb, pL8192fixed6I2R512smprtSTCBbfb,
     pL8192fixed6I2R256smprtSTCBbfb,  pL8192fixed6I2R128smprtSTCBbfb,
     pL8192fixed6I2RNOL8192smprtSTCBbfb,
     pL8192fixed6I2RFITTED8192smprtSTCBbfb,
     pL1024fixed6OF1024NM50000smprtSTCBbfb,
     pL1024multilibWEIGHTN1024NM50000smprtSTCBbfb,
     pL512fixed6OF8192smprtSTCBbfbSUM0,
     file="/home/ceballos/INSTRUMEN/EURECA/ERESOL/methodsForR.Rdat")

#===========================================================================================

# for methods comparison

if (length(grep("shortFilters",gainScaleID)) > 0){
    methods <- list(pL8192fixed6OF8192smprtSTCBbfb, pL8192fixed6OF4096smprtSTCBbfb,
                    pL8192fixed6OF2048smprtSTCBbfb, pL8192fixed6OF1024smprtSTCBbfb, 
                    pL8192fixed6OF512smprtSTCBbfb,  pL8192fixed6OF256smprtSTCBbfb, 
                    pL8192fixed6OF128smprtSTCBbfb,
                    pL8192fixed6I2R8192smprtSTCBbfb, 
                    pL8192fixed6I2RNOL8192smprtSTCBbfb,
                    pL8192fixed6I2RFITTED8192smprtSTCBbfb
    )
    cat("Using method for ShortFilters\n")
}else if(length(grep("longFilter",gainScaleID)) > 0){
    methods <- list(pL8192fixed6OF8192smprtSTCBbfb, pL4096fixed6OF8192smprtSTCBbfb,
                    pL2048fixed6OF8192smprtSTCBbfb, pL1024fixed6OF8192smprtSTCBbfb, 
                    pL512fixed6OF8192smprtSTCBbfb,  pL256fixed6OF8192smprtSTCBbfb, 
                    pL128fixed6OF8192smprtSTCBbfb,  
                    pL8192fixed6I2R8192smprtSTCBbfb, pL4096fixed6I2R8192smprtSTCBbfb,
                    pL2048fixed6I2R8192smprtSTCBbfb, pL1024fixed6I2R8192smprtSTCBbfb,
                    pL512fixed6I2R8192smprtSTCBbfb,  pL256fixed6I2R8192smprtSTCBbfb,
                    pL128fixed6I2R8192smprtSTCBbfb #,
     #               pL8192fixed6I2RNOL8192smprtSTCBbfb,
     #               pL8192fixed6I2RFITTED8192smprtSTCBbfb
    )
    #methods <- list(pL8192fixed6OF8192smprtSTCBbfb, pL4096fixed6OF8192smprtSTCBbfb,
    #                pL2048fixed6OF8192smprtSTCBbfb, pL1024fixed6OF8192smprtSTCBbfb, 
    #                pL512fixed6OF8192smprtSTCBbfb,  pL256fixed6OF8192smprtSTCBbfb, 
    #                pL128fixed6OF8192smprtSTCBbfb,  
    #                pL8192fixed6I2R8192smprtSTCBbfb,
    #                pL8192fixed6I2RNOL8192smprtSTCBbfb,
    #                pL8192fixed6I2RFITTED8192smprtSTCBbfb
    #)
    #methods <- list(pL8192fixed6OF8192smprtSTCBbfb, 
    #                pL8192fixed6I2R8192smprtSTCBbfb,
    #                pL8192fixed6I2RNOL8192smprtSTCBbfb,
    #                pL8192fixed6I2RFITTED8192smprtSTCBbfb
    #)
    
    cat("Using method for LongFilter\n")
}else if(length(grep("allFilter",gainScaleID)) > 0){
    methods <- list(pL8192fixed6OF8192smprtSTCBbfb, pL4096fixed6OF8192smprtSTCBbfb,
                    pL2048fixed6OF8192smprtSTCBbfb, pL1024fixed6OF8192smprtSTCBbfb, 
                    pL512fixed6OF8192smprtSTCBbfb,  pL256fixed6OF8192smprtSTCBbfb, 
                    pL128fixed6OF8192smprtSTCBbfb, 
                    pL8192fixed6OF4096smprtSTCBbfb, pL8192fixed6OF2048smprtSTCBbfb,
                    pL8192fixed6OF1024smprtSTCBbfb, pL8192fixed6OF512smprtSTCBbfb,
                    pL8192fixed6OF256smprtSTCBbfb, pL8192fixed6OF128smprtSTCBbfb,
                    pL8192fixed6OF4096smprtSTCBbfb_pB50,
                    pL8192fixed6OF2048smprtSTCBbfb_pB50, pL8192fixed6OF1024smprtSTCBbfb_pB50, 
                    pL8192fixed6OF512smprtSTCBbfb_pB50,  pL8192fixed6OF256smprtSTCBbfb_pB50,  
                    pL8192fixed6OF128smprtSTCBbfb_pB50,
                    pL8192fixed6I2R8192smprtSTCBbfb, pL4096fixed6I2R8192smprtSTCBbfb,
                    pL2048fixed6I2R8192smprtSTCBbfb, pL1024fixed6I2R8192smprtSTCBbfb,
                    pL512fixed6I2R8192smprtSTCBbfb,  pL256fixed6I2R8192smprtSTCBbfb,
                    pL128fixed6I2R8192smprtSTCBbfb,
                    pL8192fixed6I2R4096smprtSTCBbfb, pL8192fixed6I2R2048smprtSTCBbfb,
                    pL8192fixed6I2R1024smprtSTCBbfb, pL8192fixed6I2R512smprtSTCBbfb,
                    pL8192fixed6I2R256smprtSTCBbfb, pL8192fixed6I2R128smprtSTCBbfb,
                    pL8192fixed6I2R4096smprtSTCBbfb_pB50, pL8192fixed6I2R2048smprtSTCBbfb_pB50,
                    pL8192fixed6I2R1024smprtSTCBbfb_pB50, pL8192fixed6I2R512smprtSTCBbfb_pB50,
                    pL8192fixed6I2R256smprtSTCBbfb_pB50,  pL8192fixed6I2R128smprtSTCBbfb_pB50,
                    pL8192fixed6I2R4096smprtSTCBbfb_pB100, pL8192fixed6I2R2048smprtSTCBbfb_pB100,
                    pL8192fixed6I2R1024smprtSTCBbfb_pB100, pL8192fixed6I2R512smprtSTCBbfb_pB100,
                    pL8192fixed6I2R256smprtSTCBbfb_pB100,  pL8192fixed6I2R128smprtSTCBbfb_pB100,
                    pL1024fixed6OF1024NM50000smprtSTCBbfb,
                    pL1024multilibWEIGHTN1024NM50000smprtSTCBbfb,
                    pL512fixed6OF8192smprtSTCBbfbSUM0
    )
    #methods <- list(pL1024multilibWEIGHTN1024NM50000smprtSTCBbfb)
    cat("Using method for allFilter\n")
    
}
nmethods <- length(methods)

# FWHM vs Energy Gain scale plot
#--------------------------------
npolyInit <- 4 # degree of polynomial to be fitted
pulsesCateg <- c("all")
fwhmUNCORR  <- array(data=NA,dim=c(length(EkeV),nmethods)) # FWHM of Erecons

# READ CALIBRATION DATA
# ======================
meanEkeVfilt <- array(data=NA,dim=c(length(EkeV),nmethods))
errmean      <- array(data=NA,dim=c(length(EkeV),nmethods))
sigrobust    <- array(data=NA,dim=c(length(EkeV),nmethods))
for (ie in 1:length(EkeV)){
    for (im in 1:nmethods){
        nSamples <- methods[[im]]$nSamples
        pulseLength <- methods[[im]]$pLength
        ofLength <- methods[[im]]$ofLength
        TRIGG <- methods[[im]]$detMethod
        samprateStr <-methods[[im]]$samprateStr
        jitterStr <- methods[[im]]$jitterStr
        noiseStr <- methods[[im]]$noiseStr
        bbfbStr <- methods[[im]]$bbfbStr
        lib <- methods[[im]]$lib
        pBstr <- methods[[im]]$pBstr
        separation <- 40000
        if(samprateStr == "_samprate2") separation <- 20000
        if(samprateStr == "_samprate4") separation <- 10000
        #eresolFile <- paste("gainScale/eresol_",nSimPulses,"p_SIRENA",nSamples,
        #                    "_pL",pulseLength,"_", EkeV[ie],"keV_",TRIGG,"_T_", 
        #                    lib,ofLength,samprateStr,jitterStr,noiseStr,bbfbStr,".json",sep="")
        eresolFile <- paste("gainScale/eresol_",nSimPulses,"p_SIRENA",nSamples,
                    "_pL",pulseLength,"_", EkeV[ie],"keV_",methods[[im]]$name,".json",sep="")
        #eventsFile <- paste("gainScale/events_sep",separation,"sam_",nSimPulses,"p_SIRENA",nSamples,
        #                    "_pL",pulseLength,"_", EkeV[ie],"keV_",TRIGG,"_T_", 
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
        sigrobust[ie,im] <- 0.7413*(quantile(EkeVrecons,0.75)-quantile(EkeVrecons,0.25))

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
alias <- sapply(methods, function(x) paste("pL",x$pLength,"_",x$name,sep=""))
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
       cex=0.5,text.col=colors, bty="n")
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

#
# Plot <Erecons>/Ereal vs <Erecons>
# ===================================
#plot(EkeV,EkeV,type="n",mgp=c(2,1,0),xlim=c(0,30),ylim=c(-2,3),
#     ylab="Input Energy/Recons Energy", xlab="Reconstructed energy (keV)",
#     main="GAIN SCALE (II)",log="y")
library("splines")
drawLogPlotBox(xlimits=c(0.1,28),ylimits=c(0.3,150),
               logxy="y", xlabel="Reconstructed energy (keV)", 
               ylabel="Input Energy/Recons Energy",
               naxes=c(T,T,F,F))
minor.tick(nx=5,ny=5,tick.ratio=0.5)
title(main="GAIN SCALE(II) - linear spline")
grid(nx=NA,ny=NULL)
abline(h=1,lty=2,col="gray")
splineList <- list(METHODS=MethodsLab, ALIAS=alias, 
                  xdata=array(data=NA,dim=c(length(EkeV),nmethods)),
                  ydata=array(data=NA,dim=c(length(EkeV),nmethods)))

for (im in 1:nmethods){
    # plot points and fit
    xmeth <- meanEkeVfilt[,im]
    ymeth <- meanEkeVfilt[,im]/EkeV
    linfun <- approxfun(x=xmeth, y=log10(ymeth), method="linear")
    pow10linfun <- function(x){
        10**linfun(x)
    }
    
    #splfun <- splinefun(x=xmeth, y=ymeth, method="monoH.FC")
    #x.knots <- quantile(xmeth,probs=seq(0.2,0.8,0.2))
    #splcub <- lm(ymeth~bs(xmeth,knots=x.knots))
    points(xmeth,ymeth,pch=methods[[im]]$point,col=methods[[im]]$color, 
               type="p",lty=1)
    arrows(x0=xmeth-sigrobust[,im],x1=xmeth+sigrobust[,im],
           y0=ymeth,y1=ymeth,length=0)
    curve(pow10linfun,from=min(xmeth),to=max(xmeth),add=TRUE)
    #curve(splfun,from=min(xmeth),to=max(xmeth),add=TRUE)
    #xfit <- seq(min(xmeth),max(xmeth),length.out=101)
    #yfit <- predict(splcub,newdata=data.frame(xmeth=xfit))
    #lines(xfit,yfit)
    splineList$xdata[,im] <- xmeth
    splineList$ydata[,im] <- ymeth
}
#legend("topright",legend=MethodsLab, col=colors, pch=points, lty=ltypes,
#       cex=0.8,text.col=colors, bty="n")
spline.json <- toJSON(splineList,indent=1)
write(spline.json, splineFile)

#
# Save fitting points (in JSON) for fitting/correction in convertEnergies.py
# 


setwd("/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/")
dev.off()

