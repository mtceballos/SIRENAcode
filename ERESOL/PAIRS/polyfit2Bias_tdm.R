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
array <- "LPA2.5a" # "LPA75um", "LPA2shunt"
dre="FLL" # or "BBFB"
Ifit=-21294.27
nSimPulses <- "5000" # "20000"
#nSimPulses <- "1"
separation <- 40000 #for samprate (corrected below for samprate2)
gainScaleID <-"methods_fll"   # !!!! CHECK METHODS BELOW !!!!!

EkeV <- c(0.2,0.5,1,2,3,4,5,6,7,8,9,10,11,12)

setwd(paste("/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresol",array,sep=""))
pdf(paste("./PDFs/polyfit2Bias_",gainScaleID,".pdf",sep=""),width=10, height=7)
coeffsFile <- paste("coeffs_polyfit_",gainScaleID,".dat",sep="")
splineFile <- paste("spline_forfit_",gainScaleID,".json",sep="")

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

# Create models
if (dre == "FLL" || dre == "TDM"){
    reconMeths=c("OPTFILT","I2R","I2RFITTED")
    pBs=c(0)
    jitterStr=""
    noiseStr=""
    IfitStr=""
    pLengths=c(8192)
    ofLengths=c(8192)#, 4096, 2048, 1024, 512, 256, 128, 32, 16, 8)
    nSamples=8192
    dreStr="_fll"
    samprateStr=""
    nmods <- length(reconMeths)*length(pLengths)*length(ofLengths)*length(pBs)
    models <- vector("list", nmods)
    maincol=list("OPTFILT"="blueviolet", "I2R"="darkorange", "I2RFITTED"="darkgreen")

    i=0
    for (meth in reconMeths){
        if (meth == "I2RFITTED"){
            IfitStr=paste("Ifit_", abs(as.integer(Ifit)),sep="")
            if (Ifit < 0) IfitStr=paste("_Ifit_m", abs(as.integer(Ifit)),sep="")
        } 
        for (pB in pBs){
            pBstr <- paste("_pB",pB,sep="")
            if (pB == 0) pBstr=""
            ipt = ceil(pB/10)
            il=1
            for (pLength in pLengths){
                for (ofLength in ofLengths){
                    i = i+1
                    label<-paste("OF_",meth,"(pL",pLength,",ofL",ofLength,",6keV, STC, s1,",dreStr,
                                ",pB",pB,"Ifit",Ifit,")",sep="")
                    name<-paste("STC_T_fixedlib6OF_",meth,ofLength,pBstr,IfitStr,dreStr,sep="")
                    cat("Defining model", name,"\n")
                    lib<-paste("fixedlib6OF_",meth,sep="")
                    models[[i]] <- list(name=name, nSamples=nSamples, samprateStr=samprateStr, 
                                     jitterStr=jitterStr, noiseStr=noiseStr, pLength=pLength, 
                                     dreStr=dreStr, IfitStr=IfitStr, detMethod="STC", lib=lib, 
                                     ofLength=ofLength, color=maincol[[meth]], point=ipt, ltype=il, lab=label)
                    il = il + 1
                } #ofL
            }#pL
        } #pB
    }#recon
}

if (dre == "FLL" || dre == "TDM"){
    save(models, file="/home/ceballos/INSTRUMEN/EURECA/ERESOL/methodsForR_TDM.Rdat")
}
#===========================================================================================

# for methods comparison
methods <- models
nmethods <- length(methods)

########## # FWHM vs Energy Gain scale plot
########## #--------------------------------
npolyInit <- 4 # degree of polynomial to be fitted
pulsesCateg <- c("all")

##########fwhmUNCORR  <- array(data=NA,dim=c(length(EkeV),nmethods)) # FWHM of Erecons

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
        dreStr <- methods[[im]]$dreStr
        lib <- methods[[im]]$lib
        base <- sub("_", "", methods[[im]]$baseStr)
        
        separation <- 40000
        if(samprateStr == "_samprate2") separation <- 20000
        if(samprateStr == "_samprate4") separation <- 10000
        
        eventsFile <- paste("gainScale/events_sep",separation,"sam_",nSimPulses,"p_SIRENA",
                            nSamples,"_pL",pulseLength,"_", EkeV[ie],"keV_",
                            methods[[im]]$name,"_HR.fits",sep="")
        if (length(base) > 0){
            eventsFile <- paste("gainScale/",base,"/events_sep",separation,"sam_",nSimPulses,"p_SIRENA",
                                nSamples,"_pL",pulseLength,"_", EkeV[ie],"keV_",
                                methods[[im]]$name,"_HR.fits",sep="")
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

    } # for each method
} # foreach energy
cat("Writing coeffs:",coeffsFile,"\n")
##########fwhm  <- fwhmUNCORR

# Define plot layout (colors, tables)
# ====================================
colors <- sapply(methods, function(x) x$color)
MethodsLab <-  sapply(methods, function(x) x$lab)
alias <- sapply(methods, function(x) paste("pL",x$pLength,"_",x$name,x$baseStr,sep=""))
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

