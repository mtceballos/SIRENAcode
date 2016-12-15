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
fixed1          <- list(name="fixedlib1_OPTFILT",   color="blue", point=0, ltype=2,
                        lab="OPTIMAL FILTERING (fixed 1keV filter)")
fixed1OF        <- list(name="fixedlib1OF_OPTFILT", color="blue", point=1, ltype=2,
                        lab="OPTIMAL FILTERING (fixed OF 1keV filter)")

multi_I2R       <- list(name="multilib_I2R",        color="magenta", point=15, ltype=1,
                        lab="Resistance Space R (First order exp.)")
fixed1_I2R      <- list(name="fixedlib1_I2R",       color="magenta", point=0, ltype=2,
                        lab="Resistance Space R (fixed 1keV filter)")
fixed1OF_I2R    <- list(name="fixedlib1OF_I2R",     color="magenta", point=1, ltype=2,
                        lab="Resistance Space R (fixed OF 1keV filter)")

multi_I2RALL    <- list(name="multilib_I2RALL",     color="purple3", point=15, ltype=1,
                        lab="Resistance Space RALL (First order exp.)")
fixed1_I2RALL   <- list(name="fixedlib1_I2RALL",    color="purple3", point=0, ltype=2,
                        lab="Resistance Space RALL (fixed 1keV filter)")
fixed1OF_I2RALL <- list(name="fixedlib1OF_I2RALL",  color="purple3", point=1, ltype=2,
                        lab="Resistance Space RALL (fixed OF 1keV filter)")

multi_I2RNOL    <- list(name="multilib_I2RNOL",     color="red", point=15, ltype=1,
                        lab="Resistance Space RNOL (First order exp.)")
fixed1_I2RNOL   <- list(name="fixedlib1_I2RNOL",    color="red", point=0, ltype=2,
                        lab="Resistance Space RNOL (fixed 1keV filter)")
fixed1OF_I2RNOL <- list(name="fixedlib1OF_I2RNOL",  color="red", point=1, ltype=2,
                        lab="Resistance Space RNOL (fixed OF 1keV filter)")

multi_I2RFITTED <- list(name="multilib_I2RFITTED",  color="orange", point=15, ltype=1,
                           lab="Resistance Space RFITTED (First order exp.)")
fixed1_I2RFITTED<- list(name="fixedlib1_I2RFITTED", color="orange", point=0, ltype=2,
                           lab="Resistance Space RFITTED (fixed 1keV filter)")
fixed1OF_I2RFITTED<- list(name="fixedlib1OF_I2RFITTED", color="orange", point=1, ltype=2,
                        lab="Resistance Space RFITTED (fixed OF 1keV filter)")

weight          <- list(name="multilib_WEIGHT",     color="darkgreen", point=15, ltype=1,
                           lab="COVARIANCE MATRICES (Fixen)")
weightn         <- list(name="multilib_WEIGHTN", color="darkgreen", point=0, ltype=2,
                           lab="COVARIANCE MATRICES (0(n))")
weightnOF       <- list(name="multilibOF_WEIGHTN", color="darkgreen", point=1, ltype=2,
                        lab="COVARIANCE MATRICES OF (0(n))")

save(multi, multi_I2R, multi_I2RALL, multi_I2RNOL, multi_I2RFITTED, 
     fixed1, fixed1_I2R, fixed1_I2RALL, fixed1_I2RNOL, fixed1_I2RFITTED, 
     fixed1OF, fixed1OF_I2R, fixed1OF_I2RALL, fixed1OF_I2RNOL, fixed1OF_I2RFITTED,
     weight, weightn, weightnOF,file="/home/ceballos/INSTRUMEN/EURECA/ERESOL/methodsForR.Rdat")
#===========================================================================================
methods<- list(multi, fixed1, multi_I2R, fixed1_I2R,  multi_I2RALL, fixed1_I2RALL, 
               multi_I2RNOL, fixed1_I2RNOL, multi_I2RFITTED, fixed1_I2RFITTED,
               weight, weightn)

#methods <-c("multilib_OPTFILT_NTRIG","fixedlib1_OPTFILT_NTRIG",
#            "multilib_I2R_NTRIG", "fixedlib1_I2R_NTRIG",
#            "multilib0.1keV0.5keVDAB_OPTFILT_NTRIG", "multilib0.5keV1keVDAB_OPTFILT_NTRIG",
#            "multilib1keV2keVDAB_OPTFILT_NTRIG", "multilib2keV3keVDAB_OPTFILT_NTRIG",
#            "multilib3keV4keVDAB_OPTFILT_NTRIG", "multilib4keV5keVDAB_OPTFILT_NTRIG",
#            "multilib5keV6keVDAB_OPTFILT_NTRIG", "multilib6keV7keVDAB_OPTFILT_NTRIG",
#            "multilib7keV8keVDAB_OPTFILT_NTRIG", "multilib8keV9keVDAB_OPTFILT_NTRIG",
#            "multilib9keV10keVDAB_OPTFILT_NTRIG", "multilib10keV11keVDAB_OPTFILT_NTRIG"
#            )
ifim <- 100 # model index for interval methods  (100 for no interval methods at all)
nmethods <- length(methods)


# FWHM vs Energy Gain scale plot
#--------------------------------
plotEEfit   <- 1  # FWHM vs. Energy
plotEEall   <- 0  # FWHM vs. Energy also for intermediate energies
npolyInit <- 4 # degree of polynomial to be fitted

#
# Initialize variables
#-----------------------
array <- "LPA1shunt"
nSimPulses <- "20000"
nSamples <- "2048"
pulseLength <- "2048"
separation <- 20000
setwd(paste("/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresol",array,sep=""))
pdf(paste("./PDFs/polyfit2Bias_",nSamples,"s_",nSimPulses,"p.pdf",sep=""),width=10, height=7)
EkeV <- c(0.2,0.5,1,2,3,4,5,6,7,8,9,10)
EkeV <- c(0.2,0.5,1,2,3,4,5,6,7,8)
EkeVall <- c(0.2,0.2,0.4,0.5,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,
             4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7,7.2,7.4,7.6,7.8,8,8.2,
             8.4,8.6,8.8,9,9.2,9.4,9.6,9.8,10,10.2,10.4,10.6,10.8)


fwhmPrimUNCORR <- matrix(NA,nrow=length(EkeV), ncol=length(methods)) # FWHM of Erecons
fwhmSecUNCORR  <- matrix(NA,nrow=length(EkeV), ncol=length(methods)) # FWHM of Erecons
ebiasPrim     <- matrix(NA,nrow=length(EkeV), ncol=length(methods))
ebiasSec      <- matrix(NA,nrow=length(EkeV), ncol=length(methods))

# READ CALIBRATION DATA
# ======================
for (ie in 1:length(EkeV)){
    # due to a problem in tessim, initial sample in 0.2 keV and 1 keV is not always fixed:
    # TRIGGERING is required
    if (EkeV[ie] == 0.2 || EkeV[ie] == 1){
        TRIGG = ""
    }else{ # TRIGGERING is not required
        TRIGG = "_NTRIG"
    }
    for (im in 1:length(methods)){
        # Use 3000s libraries for all methods
        
        eresolFile <- paste("eresol_",nSimPulses,"p_SIRENA",nSamples,"_pL",pulseLength,"_", EkeV[ie],"keV_F0F_", 
                            methods[[im]]$name,TRIGG,".json",sep="")
        
        if(file.exists(eresolFile)){
            #data <- read.table(eresolFile,header=TRUE)
            # use data for selected separation (see initial definitions)
            cat("Reading file ",eresolFile,"\n")
            jsondata <- fromJSON(file=eresolFile)
            idxSep <- which(sapply(jsondata,function(x) x$separation)==separation)
            fwhmPrimUNCORR[ie,im] <- as.numeric(jsondata[[idxSep]]$fwhmErecons$primaries)
            fwhmSecUNCORR[ie,im] <- as.numeric(jsondata[[idxSep]]$fwhmErecons$secondaries)
            ebiasPrim[ie,im] <- as.numeric(jsondata[[idxSep]]$biasErecons$primaries)
            ebiasSec[ie,im] <- as.numeric(jsondata[[idxSep]]$biasErecons$secondaries)
        }else{
            warning("Not-existing file:", eresolFile)
            fwhmPrimUNCORR[ie,im] <- NaN
            fwhmSecUNCORR[ie,im]  <- NaN
            ebiasPrim[ie,im]     <- NaN
            ebiasSec[ie,im]      <- NaN
        }
        if(is.nan(fwhmPrimUNCORR[ie,im]) || is.nan(fwhmSecUNCORR[ie,im]) || 
           is.nan(ebiasPrim[ie,im]) || is.nan(ebiasSec[ie,im])){
            warning("Error in ",eresolFile,"\n","  Non numerical values in eresol files: check event files")
        }
    } # for each method
} # foreach energy

for (use in c("PRIM","SEC")){
    coeffsFile <- paste("coeffs_polyfit_",nSamples,"_",use,".dat",sep="")
    cat("Using coeffs:",coeffsFile,"\n")
    if (use == "PRIM"){
        fwhm  <- fwhmPrimUNCORR
        ebias <- ebiasPrim
    } else if(use == "SEC"){
        fwhm  <- fwhmSecUNCORR
        ebias <- ebiasSec
    }

    if(plotEEall){
        # PLOT ENERGY_CALC vs E input for all energies (also intermediate energies)
        #===========================================================================
        if (EkeVall[ie] == 0.2 || EkeVall[ie] == 0.2){
            TRIGG = ""
        }else{
            TRIGG = "_NTRIG"
        }
        # Read ALL data BIAS
        ebiasAllPrim     <- matrix(NA,nrow=length(EkeVall), ncol=length(methods))
        ebiasAllSec      <- matrix(NA,nrow=length(EkeVall), ncol=length(methods))
        for (ie in 1:length(EkeVall)){
            for (im in 1:length(methods)){
            #for (im in 1:1){
                # Use 3000s libraries for all methods
                eresolFile <- paste("eresol_",nSimPulses,"p_SIRENA",nSamples,"_", EkeVall[ie],"keV_F0F_", 
                                    methods[[im]]$name,TRIGG,".dat",sep="")
                #cat("Reading file ",eresolFile,"\n")
                data <- read.table(eresolFile,header=TRUE)
                ebiasAllPrim[ie,im]     <- data[nrow(data),4]
                ebiasAllSec[ie,im]      <- data[nrow(data),5]
            }
        }
        if (use == "PRIM"){
            ebiasAll <- ebiasAllPrim
        } else if(use == "SEC"){
            ebiasAll <- ebiasAllSec
        }
        ebiasAllkeV <- ebiasAll/1000.
        meanEkeVallfilt <- ebiasAllkeV + EkeVall  # Ebias = <Erecons> - Einput
    }



    if(plotEEfit){
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
        plot(EkeV,EkeV,type="n",mgp=c(2,1,0),
              xlab="Input Energy keV (calibration points marked)", ylab="Reconstructed energy (keV)",
              main=paste("GAIN SCALE: AC simulations - ",array," - ",nSamples,"samples"," - ",nSimPulses,
                         "pulses"," - ", use,sep=""))
        
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
            # plot points and fit for non-interval methods (these are plotted afterwards)
            if(im < ifim){
                points(EkeV,meanEkeVfilt[,im],pch=methods[[im]]$point,col=methods[[im]]$color,
                       type="p",lty=1)
                if(plotEEall){
                    points(EkeVall,meanEkeVallfilt[,im],pch=methods[[im]]$point,
                           col=methods[[im]]$color,type="p",lty=1)}
            }
            
            # check significance of npoly fit
            # use all calibration points available for multilib or fixedlib
            iemin[im] <- 1
            iemax[im] <- length(EkeV)
            if(im > (ifim-1)){ # use only few points around interval model
                iemin[im] <- max(im-3,1)
                iemin[im] <- min(iemin[im],(length(EkeV)-3))
                iemax[im] <- min(im,length(EkeV))
                iemax[im] <- max(4,iemax[im]) #at least 4 points to fit
            }
            badCoeffs <- min(npolyInit,(iemax[im]-iemin[im]))
            npoly <- min(npolyInit,(iemax[im]-iemin[im]))
            
            # first fit with orthogonal polynomia (uncorrelated error)
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
            
            
            #fitLabel <- bquote(paste(.(MethodsLab[im]),":  ", E[filt], "(KeV) = ", .(round(a0[im],4)), " + ",
            #                           .(round(a1[im],4)), E[keV], " + ", .(round(a2[im],4)), E[keV]^2, " + ",
            #                           .(round(a3[im],5)), E[keV]^3))
            
            
            if(im<ifim){
                curve(polyCurve(x,im), lty=methods[[im]]$ltype,add=TRUE, col=methods[[im]]$color)
            
            }else{
                curve(polyCurve(x,im),from=EkeV[iemin[im]], to=EkeV[iemax[im]], 
                      lty=methods[[im]]$ltype,add=TRUE, col=methods[[im]]$color)
            }
        }
        
        legend("topleft",legend=MethodsLab, col=colors, pch=points, lty=ltypes,
               cex=0.8,text.col=colors, bty="n")
        text(8,1,paste("Lines= poly fits to calib energies",sep=""), cex=0.8)
        
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
    } #if plotEEfit

    # plot MULTILIB (im=1) vs. other methods
    if(ifim<=nmethods){
        for (im in ifim:nmethods){
            plot(EkeV,meanEkeVfilt[,im],#xlim=c(EkeV[iemin[im]],EkeV[iemax[im]]), 
                 #ylim=c(EkeV[iemin[im]],EkeV[iemax[im]]),
                pch=points[im],col=colors[10],type="p",lty=1,
                xlab="Input Energy keV (calibration points marked)", ylab="Reconstructed energy (keV)")
            for (e in EkeV){abline(v=EkeV,lty=3,col="grey80")}
            curve(polyCurve(x,im),from=EkeV[iemin[im]], to=EkeV[iemax[im]],lty=ltype[im],add=TRUE, col="blue")
            # multilib
            points(EkeV,meanEkeVfilt[,1],pch=points[1],col=colors[1],type="p",lty=1)
            curve(polyCurve(x,1),lty=ltype[1],add=TRUE, col=colors[1])
            if(plotEEall){points(EkeVall,meanEkeVallfilt[,1],pch=points[1],col=colors[1],type="p")}
        
            legend("topleft",legend=c(MethodsLab[1],MethodsLab[im]), col=c(colors[1],"blue"), 
               pch=c(points[1],points[im]),lty=c(1,3),
               cex=0.8,text.col=c(colors[1],"blue"), bty="n")
        }
    }


} # PRIMARIES & SECONDARIES
setwd("/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/")
dev.off()





