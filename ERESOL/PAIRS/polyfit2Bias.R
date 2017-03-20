# Reads first pass of getEresolCurves to obtain Ebias = <Ecalculated> - Einput fit coefficients
# A second run of getEresolCurves.py will use them to get corrected FWHM
# Get polynomial fit to a curve Ecalculated vs Einput :
#        a given array (SPA, LPA1, LPA2, LPA3)  
#        @ different energies (0.5,1,2,3,4,6,9) keV
#        for different methods (OPTFILT_mono, OPTFILT_multi, WEIGHT, WEIGHTN, I2R)

library(gridExtra)
library(Hmisc)
library(rjson)

# DEFINE and SAVE methods characteristics
#    Read them back with:
#    > load("/home/ceballos/INSTRUMEN/EURECA/ERESOL/methodsForR.Rdat")
#=======================================================================================

#
# WARNING: no need to calculate calibration gain scale curves with OFLib=yes because 
# on-the-fly filters are created at exactly the same lenght than precalc filters
#
multiMF         <- list(name="multilibMF_OPTFILT",  color="cyan", point=2, ltype=1,
                        lab="OPTIMAL FILTERING (interpolating filters MF)")

multi           <- list(name="multilib_OPTFILT",    color="cyan", point=15, ltype=1,
                        lab="OPTIMAL FILTERING (First order exp.)")
fixed1          <- list(name="fixedlib1_OPTFILT",   color="blue", point=0, ltype=1,
                        lab="OPTIMAL FILTERING (fixed 1keV filter)")
fixed1OF        <- list(name="fixedlib1OF_OPTFILT", color="blue", point=1, ltype=2,
                        lab="OPTIMAL FILTERING (fixed OF 1keV filter)")
fixed1OFNM      <- list(name="fixedlib1OF_OPTFILTNM", color="blue", point=2, ltype=2,
                        lab="OPTIMAL FILTERING NOISE MAT(fixed OF 1keV filter)")

multi_I2R       <- list(name="multilib_I2R",        color="magenta", point=15, ltype=1,
                        lab="Resistance Space R (First order exp.)")
fixed1_I2R      <- list(name="fixedlib1_I2R",       color="magenta", point=0, ltype=1,
                        lab="Resistance Space R (fixed 1keV filter)")
fixed1OF_I2R    <- list(name="fixedlib1OF_I2R",     color="magenta", point=1, ltype=2,
                        lab="Resistance Space R (fixed OF 1keV filter)")

multi_I2RALL    <- list(name="multilib_I2RALL",     color="purple3", point=15, ltype=1,
                        lab="Resistance Space RALL (First order exp.)")
fixed1_I2RALL   <- list(name="fixedlib1_I2RALL",    color="purple3", point=0, ltype=1,
                        lab="Resistance Space RALL (fixed 1keV filter)")
fixed1OF_I2RALL <- list(name="fixedlib1OF_I2RALL",  color="purple3", point=1, ltype=2,
                        lab="Resistance Space RALL (fixed OF 1keV filter)")

multi_I2RNOL    <- list(name="multilib_I2RNOL",     color="red", point=15, ltype=1,
                        lab="Resistance Space RNOL (First order exp.)")
fixed1_I2RNOL   <- list(name="fixedlib1_I2RNOL",    color="red", point=0, ltype=1,
                        lab="Resistance Space RNOL (fixed 1keV filter)")
fixed1OF_I2RNOL <- list(name="fixedlib1OF_I2RNOL",  color="red", point=1, ltype=2,
                        lab="Resistance Space RNOL (fixed OF 1keV filter)")

multi_I2RFITTED <- list(name="multilib_I2RFITTED",  color="orange", point=15, ltype=1,
                           lab="Resistance Space RFITTED (First order exp.)")
fixed1_I2RFITTED<- list(name="fixedlib1_I2RFITTED", color="orange", point=0, ltype=1,
                           lab="Resistance Space RFITTED (fixed 1keV filter)")
fixed1OF_I2RFITTED<- list(name="fixedlib1OF_I2RFITTED", color="orange", point=1, ltype=2,
                        lab="Resistance Space RFITTED (fixed OF 1keV filter)")

weight          <- list(name="multilib_WEIGHT",     color="darkgreen", point=15, ltype=1,
                           lab="COVARIANCE MATRICES (Fixen)")
weightn         <- list(name="multilib_WEIGHTN", color="darkgreen", point=0, ltype=3,
                           lab="COVARIANCE MATRICES (0(n))")
weightnOF       <- list(name="multilibOF_WEIGHTN", color="darkgreen", point=1, ltype=2,
                        lab="COVARIANCE MATRICES OF (0(n))")

save(multi, multi_I2R, multi_I2RALL, multi_I2RNOL, multi_I2RFITTED, 
     fixed1, fixed1_I2R, fixed1_I2RALL, fixed1_I2RNOL, fixed1_I2RFITTED, 
     fixed1OF, fixed1OFNM,fixed1OF_I2R, fixed1OF_I2RALL, fixed1OF_I2RNOL, fixed1OF_I2RFITTED,
     weight, weightn, weightnOF,file="/home/ceballos/INSTRUMEN/EURECA/ERESOL/methodsForR.Rdat")
#===========================================================================================
methods<- list(multi, fixed1, multi_I2R, fixed1_I2R,  multi_I2RALL, fixed1_I2RALL, 
               multi_I2RNOL, fixed1_I2RNOL, multi_I2RFITTED, fixed1_I2RFITTED,
               weight, weightn)
methods <- list(fixed1OF, fixed1OFNM, fixed1OF_I2R,fixed1OF_I2RNOL, fixed1OF_I2RFITTED, weight,weightnOF)
#methods <- list(fixed1OF, fixed1OF_I2R,fixed1OF_I2RNOL, fixed1OF_I2RFITTED, weight,weightnOF)
nmethods <- length(methods)


# FWHM vs Energy Gain scale plot
#--------------------------------
npolyInit <- 4 # degree of polynomial to be fitted
pulsesCateg <- c("primaries", "secondaries","all")
pulsesCateg <- c("all")
#
# Initialize variables
#-----------------------
array <- "LPA2shunt"
triggEners <- c(0.2,0.5,1,2) # energies where manual triggering has to be done; not anymore hopefully
triggEners <- c()
nSimPulses <- "20000"
nSamples <- "4096"
pulseLength <- "4096"
separation <- 40000
setwd(paste("/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresol",array,sep=""))
pdf(paste("./PDFs/polyfit2Bias_",nSamples,"s_",nSimPulses,"p.pdf",sep=""),width=10, height=7)
EkeV <- c(0.2,0.5,1,2,3,4,5,6,7)

fwhmUNCORR  <- array(data=NA,dim=c(length(EkeV),nmethods,length(pulsesCateg))) # FWHM of Erecons
ebiasUNCORR <- array(data=NA,dim=c(length(EkeV),nmethods,length(pulsesCateg)))

# READ CALIBRATION DATA
# ======================
for (ie in 1:length(EkeV)){
    # due to a problem in tessim, initial sample in some pulse may not be always in the same position
    # TRIGGERING is required
    if (EkeV[ie] %in% triggEners) {
        TRIGG = ""
    }else{ # TRIGGERING is not required
        TRIGG = "_NTRIG"
    }
    for (im in 1:nmethods){
        eresolFile <- paste("gainScale/eresol_",nSimPulses,"p_SIRENA",nSamples,"_pL",pulseLength,"_", 
                            EkeV[ie],"keV_F0F_", methods[[im]]$name,TRIGG,".json",sep="")
        
        if(file.exists(eresolFile)){
            #data <- read.table(eresolFile,header=TRUE)
            # use data for selected separation (see initial definitions)
            cat("Reading file ",eresolFile,"\n")
            jsondata <- fromJSON(file=eresolFile)
            idxSep <- which(sapply(jsondata,function(x) x$separation)==separation)
            
            for(ip in 1:length(pulsesCateg)){
                fwhmUNCORR[ie,im,ip]  <- as.numeric(jsondata[[idxSep]]$fwhmErecons[[pulsesCateg[ip]]])
                ebiasUNCORR[ie,im,ip] <- as.numeric(jsondata[[idxSep]]$biasErecons[[pulsesCateg[ip]]])
            }
        }else{
            warning("Not-existing file:", eresolFile)
            for(ip in 1:length(pulsesCateg)){
                fwhmUNCORR[ie,im,ip]  <- NaN
                ebiasUNCORR[ie,im,ip] <- NaN
            }
        }
        if(all(is.nan(fwhmUNCORR[ie,im,])) || all(is.nan(ebiasUNCORR[ie,im,]))){
            warning("Error in ",eresolFile,"\n","  Non numerical values in eresol files: check event files")
        }
    } # for each method
} # foreach energy

for (ip in 1:length(pulsesCateg)){ # only for LPA1; LPA2 simulates only single pulses
    coeffsFile <- paste("coeffs_polyfit_",nSamples,"_",pulsesCateg[ip],".dat",sep="")
    cat("Writing coeffs:",coeffsFile,"\n")
    fwhm  <- fwhmUNCORR
    ebias <- ebiasUNCORR

    colors <- sapply(methods, function(x) x$color)
    MethodsLab <-  sapply(methods, function(x) x$lab)
    alias <- sapply(methods, function(x) x$name)
    points <- sapply(methods, function(x) x$point)
    ltypes <- sapply(methods, function(x) x$ltype)
        
    tt3 <- ttheme_minimal(
        core=list(bg_params = list(fill = colors,col=NA),
                  fg_params=list(fontface=3,cex=0.8)), 
        colhead=list(fg_params=list(col="navyblue")),
        rowhead=list(fg_params=list(col="navyblue", fontface=2L)))
    
    tt3 <- ttheme_minimal(
        core=list(bg_params = list(col=NA, fill="ivory"),
                  fg_params=list(fontface=3,cex=0.8, col=colors)), 
        colhead=list(fg_params=list(col="navyblue")),
        rowhead=list(fg_params=list(col="navyblue", fontface=2L)))
    
    # PLOT ENERGY_CALC vs E input
    #=================================
    
    ebiaskeV <- ebias/1000.
    meanEkeVfilt <- ebiaskeV + EkeV  # Ebias = <Erecons> - Einput
    plot(EkeV,EkeV,type="n",mgp=c(2,1,0),ylim=c(min(EkeV),max(EkeV)+2),
         xlab="Input Energy keV (calibration points marked)", ylab="Reconstructed energy (keV)",
         main=paste("GAIN SCALE: AC simulations - ",array," - ",nSamples,"samples"," - ",nSimPulses,
                    "pulses"," - ", pulsesCateg[ip],sep=""))
    
    #grid()
    minor.tick(nx=5,ny=5,tick.ratio=0.5)
    grid(nx=NA,ny=NULL)
    lines(c(0,15),c(0,15),lty=2,col="gray")
    for (e in EkeV){abline(v=EkeV,lty=3,col="grey80")}
        
    coeffs <- matrix(0,nrow=nmethods,ncol=(npolyInit+1))
    polyCurve <- function(x,im) {
        polyres <- 0
        for (icoeff in 1:(npolyInit)){
            index <- icoeff - 1
            polyres <- polyres + coeffs[im,icoeff] * x^(index)
        }
        #a3[im]*x^3 + a2[im]*x^2 + a1[im]*x + a0[im]
        return(polyres)
    }
    iemin <- c()
    iemax <- c()
    for (im in 1:nmethods){
        cat("Fitting method ",methods[[im]]$lab,"\n")
        # plot points and fit
        points(EkeV,meanEkeVfilt[,im,ip],pch=methods[[im]]$point,col=methods[[im]]$color, 
               type="p",lty=1)
        
        # check significance of npoly fit
        # use all calibration points available for multilib or fixedlib
        iemin[im] <- 1
        iemax[im] <- length(EkeV)
        badCoeffs <- min(npolyInit,(iemax[im]-iemin[im]))
        npoly <- min(npolyInit,(iemax[im]-iemin[im]))
        
        # first fit with orthogonal polynomia (uncorrelated error)
        while (badCoeffs > 0){
            stopifnot(npoly > 0)
            fit <- lm(meanEkeVfilt[iemin[im]:iemax[im],im,ip]~
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
        fit <- lm(meanEkeVfilt[iemin[im]:iemax[im],im,ip]~
                      poly(EkeV[iemin[im]:iemax[im]],npoly,raw=TRUE))
        coeffs[im,1:(npoly+1)] <- fit$coefficients[1:(npoly+1)]
        
        curve(polyCurve(x,im), lty=methods[[im]]$ltype,add=TRUE, col=methods[[im]]$color)
    }
    
    legend("topleft",legend=MethodsLab, col=colors, pch=points, lty=ltypes,
           cex=0.8,text.col=colors, bty="n")
    text(7,1,paste("Lines= poly fits to calib energies",sep=""), cex=0.8)
    
    coeffsTab <- data.frame(METHODS=MethodsLab, ALIAS=alias, coeffs)
    cnames <- c("METHODS","ALIAS")
    for (i in 0:npolyInit){
        colname <- paste("a",i,sep="")
        cnames <- append(cnames,colname)
    }
    colnames(coeffsTab) <- cnames
    # save coefficients to be read by compareMethods.R
    write.table(coeffsTab, file=coeffsFile, row.names=FALSE)
    
    # plot table with results
    plot(0:1,0:1, xlab="",ylab="", pch=NA_integer_, axes=FALSE)
    mtext(bquote(paste(E[filt], " (keV)= a0 + a1 * ",E[cal]," + a2 * ", 
                       E[cal]^2, " + a3 * ",E[cal]^3, "+...", sep="")), side=3,cex=1.3)
    grid.table(coeffsTab, rows=NULL, theme=tt3)
    #legend("topleft", legend=fitLabel, title="Polynomial fit to Erecons vs Einput relation")

} # PRIMARIES & SECONDARIES & ALL...
setwd("/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/")
dev.off()





