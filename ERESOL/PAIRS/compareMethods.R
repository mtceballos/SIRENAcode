#
# Plot energy resolution for:
#        a given array (SPA, LPA1, LPA2, LPA3, LPA75)  
#        @ different energies (0.2,0.5,1,2,3,4,5,6,7,8) keV
#        for different methods  (fixed and multilib): OPTFILT, WEIGHT, WEIGHTN, I2R, I2RALL, I2RNOL, I2RFITTED
# calculated with pulses separated by 40000 samples (samprate 156250Hz)

rm(list=ls())
library(rjson)
library(Hmisc)
source("~/R/Rfunctions/drawLogPlotBox.r")

# LOAD methods characteristics
load("/home/ceballos/INSTRUMEN/EURECA/ERESOL/methodsForR.Rdat") #created @ polyfit2Bias.R

plotFWHM_GAINCORRS <- 0   # FWHM vs. Energy
plotFWHM_GAINCORRE <- 0   # FWHM vs. Energy
plotBiasCorrFit    <- 0   # Bias-gainScaleCorrected vs. separation
plotFWHM_rlength   <- 1   # FWHM vs. Record Length     
# if(plotFWHM_GAINCORRS || plotFWHM_GAINCORRE){
#     pdfName <- "fwhmVSenergy_methods"
# }else if (plotBiasCorrFit){
#     pdfName <- "biasVSsep_methods"
# }else if (plotFWHM_rlength){    
#     pdfName <- "fwhmVSrlenght_methods"
# }
subtitle <- "ALL"
#subtitle <- "PERF"
subtitle <- "OPTFILT_long_basePadding"
subtitle <- "OPTFILT_long_zeroPadding_ADC_I2R"
#subtitle <- "OPTFILT_short"
subtitle <- "allFilter_ADC_I2R"
#subtitle <- "OPTFILT"
#subtitle <- "SPIE2016" # SPIE2016/SPIE2016PP (to plot also PP points)
plttype="b"
useGainCorr <- "all" # or "SEC" (pulses used for Gain Curve)
pulsesCateg <- c("primaries", "secondaries","all")
pulsesCateg <- c("all")
par(pty="s")
separation <- "40000"
calib = "surface" # or "curve"
#calib = "curve"
array <- "LPA2shunt"
array <- "LPA75um"
nSimPulses <- "2000"
nIntervals <- "150000"
noiseMat <-paste("_noiseMat",nIntervals,sep="")
noiseMat=""
invalids <- 494
bbfb <-"_bbfb"

setwd(paste("/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresol",array,sep=""))

if (!subtitle %in% c("PERF" )){ # NOT PERF
    
    if(subtitle == "ALL"){
        methods <- list(multi, fixed1, fixed1OF, 
                        multi_I2R, fixed1_I2R, fixed1OF_I2R,
                        multi_I2RALL, fixed1_I2RALL, fixed1OF_I2RALL, 
                        multi_I2RNOL, fixed1_I2RNOL, fixed1OF_I2RNOL,
                        multi_I2RFITTED, fixed1_I2RFITTED, fixed1OF_I2RFITTED,
                        weight, weightn)
    }else if(length(i<-grep("SPIE2016",subtitle))){
        methods <- list(fixed1OF, fixed1OF_I2RNOL, fixed1OF_I2R, fixed1OF_I2RFITTED, 
                        weightnOF,weight)
    }else if(subtitle == "MULTILIB"){
        methods <- list(multi, multi_I2R, multi_I2RALL, multi_I2RNOL,
                        weight, weightn)
    }else if(subtitle == "FIXEDLIBOF"){
        methods <- list(fixed1OF, fixed1OFNM, fixed1OF_I2R, fixed1OF_I2RNOL, fixed1OF_I2RFITTED, 
                        weightnOF)
    }else if(subtitle == "FIXEDLIB6OF"){
        methods <- list(fixed6OFsmprtAD, fixed6OFsmprt2AD,fixed6OFsmprtA1, fixed6OFsmprt2A1)
    }else if(subtitle == "FIXEDLIB"){
        methods <- list(fixed1OF, fixed1OF_I2R, fixed1OF_I2RNOL, fixed1OF_I2RFITTED,fixed1OFNM,
                        fixed1, fixed1_I2R, fixed1_I2RNOL, fixed1_I2RFITTED, weightnOF,  weight)        
        #methods <- list(fixed1, fixed1OFNM, fixed1_I2R, fixed1_I2RNOL, fixed1_I2RFITTED)        
    }else if(subtitle == "OPTFILT_short"){
        methods <- list(pL8192fixed6OF8192smprtSTCBbfb, pL8192fixed6OF4096smprtSTCBbfb,
                        pL8192fixed6OF1024smprtSTCBbfb, pL8192fixed6OF512smprtSTCBbfb,
                        pL8192fixed6OF256smprtSTCBbfb, pL8192fixed6OF128smprtSTCBbfb,
                        pL8192fixed6I2R8192smprtSTCBbfb,pL8192fixed6I2RFITTED8192smprtSTCBbfb,
                        pL8192fixed6I2RNOL8192smprtSTCBbfb
                        )
    }else if(subtitle == "OPTFILT_long_basePadding"){
        methods <- list(pL8192fixed6OF8192smprtSTCBbfb, pL4096fixed6OF8192smprtSTCBbfb,
                        pL1024fixed6OF8192smprtSTCBbfb, pL512fixed6OF8192smprtSTCBbfb,
                        pL256fixed6OF8192smprtSTCBbfb, pL128fixed6OF8192smprtSTCBbfb,
                        pL8192fixed6I2R8192smprtSTCBbfb,pL8192fixed6I2RFITTED8192smprtSTCBbfb,
                        pL8192fixed6I2RNOL8192smprtSTCBbfb
            )
    }else if(subtitle == "OPTFILT_long_zeroPadding_ADC_I2R"){
        methods <- list(pL8192fixed6OF8192smprtSTCBbfb, pL4096fixed6OF8192smprtSTCBbfb,
                        pL2048fixed6OF8192smprtSTCBbfb, pL1024fixed6OF8192smprtSTCBbfb, 
                        pL512fixed6OF8192smprtSTCBbfb,  pL256fixed6OF8192smprtSTCBbfb, 
                        pL128fixed6OF8192smprtSTCBbfb,
                        pL8192fixed6I2R8192smprtSTCBbfb,pL4096fixed6I2R8192smprtSTCBbfb,
                        pL2048fixed6I2R8192smprtSTCBbfb, pL1024fixed6I2R8192smprtSTCBbfb, 
                        pL512fixed6I2R8192smprtSTCBbfb,  pL256fixed6I2R8192smprtSTCBbfb, 
                        pL128fixed6I2R8192smprtSTCBbfb
        )
    }else if(subtitle == "allFilter_ADC_I2R"){
        cat("reading ok\n")
        methods <- list(pL8192fixed6OF8192smprtSTCBbfb, pL4096fixed6OF8192smprtSTCBbfb,
                        pL2048fixed6OF8192smprtSTCBbfb, pL1024fixed6OF8192smprtSTCBbfb, 
                        pL512fixed6OF8192smprtSTCBbfb,  pL256fixed6OF8192smprtSTCBbfb, 
                        pL128fixed6OF8192smprtSTCBbfb,
                        pL8192fixed6OF4096smprtSTCBbfb, pL8192fixed6OF2048smprtSTCBbfb,
                        pL8192fixed6OF1024smprtSTCBbfb, pL8192fixed6OF512smprtSTCBbfb,
                        pL8192fixed6OF256smprtSTCBbfb, pL8192fixed6OF128smprtSTCBbfb,
                        #pL8192fixed6OF4096smprtSTCBbfb_pB50,pL8192fixed6OF2048smprtSTCBbfb_pB50,
                        #pL8192fixed6OF1024smprtSTCBbfb_pB50, pL8192fixed6OF512smprtSTCBbfb_pB50,
                        #pL8192fixed6OF256smprtSTCBbfb_pB50, pL8192fixed6OF128smprtSTCBbfb_pB50,
                        pL8192fixed6OF4096smprtSTCBbfb_pB75, pL8192fixed6OF2048smprtSTCBbfb_pB75,
                        pL8192fixed6OF1024smprtSTCBbfb_pB75, pL8192fixed6OF512smprtSTCBbfb_pB75,
                        pL8192fixed6OF256smprtSTCBbfb_pB75,  pL8192fixed6OF128smprtSTCBbfb_pB75,
                        pL8192fixed6I2R8192smprtSTCBbfb, pL4096fixed6I2R8192smprtSTCBbfb,
                        pL2048fixed6I2R8192smprtSTCBbfb, pL1024fixed6I2R8192smprtSTCBbfb, 
                        pL512fixed6I2R8192smprtSTCBbfb,  pL256fixed6I2R8192smprtSTCBbfb, 
                        pL128fixed6I2R8192smprtSTCBbfb,
                        pL8192fixed6I2R4096smprtSTCBbfb, pL8192fixed6I2R2048smprtSTCBbfb,
                        pL8192fixed6I2R1024smprtSTCBbfb, pL8192fixed6I2R512smprtSTCBbfb,
                        pL8192fixed6I2R256smprtSTCBbfb, pL8192fixed6I2R128smprtSTCBbfb,
                        #pL8192fixed6I2R4096smprtSTCBbfb_pB25,pL8192fixed6I2R2048smprtSTCBbfb_pB25,
                        #pL8192fixed6I2R1024smprtSTCBbfb_pB25, pL8192fixed6I2R512smprtSTCBbfb_pB25,
                        #pL8192fixed6I2R256smprtSTCBbfb_pB25, pL8192fixed6I2R128smprtSTCBbfb_pB25,
                        #pL8192fixed6I2R4096smprtSTCBbfb_pB50,pL8192fixed6I2R2048smprtSTCBbfb_pB50,
                        #pL8192fixed6I2R1024smprtSTCBbfb_pB50, pL8192fixed6I2R512smprtSTCBbfb_pB50,
                        #pL8192fixed6I2R256smprtSTCBbfb_pB50, pL8192fixed6I2R128smprtSTCBbfb_pB50,
                        pL8192fixed6I2R4096smprtSTCBbfb_pB75, pL8192fixed6I2R2048smprtSTCBbfb_pB75,
                        pL8192fixed6I2R1024smprtSTCBbfb_pB75, pL8192fixed6I2R512smprtSTCBbfb_pB75,
                        pL8192fixed6I2R256smprtSTCBbfb_pB75,  pL8192fixed6I2R128smprtSTCBbfb_pB75,
                        #pL8192fixed6I2R4096smprtSTCBbfb_pB85, pL8192fixed6I2R2048smprtSTCBbfb_pB85,
                        #pL8192fixed6I2R1024smprtSTCBbfb_pB85, pL8192fixed6I2R512smprtSTCBbfb_pB85,
                        #pL8192fixed6I2R256smprtSTCBbfb_pB85,  pL8192fixed6I2R128smprtSTCBbfb_pB85,
                        pL1024fixed6OF1024NM150000smprtSTCBbfb,
                        pL1024multilibWEIGHTN1024NM150000smprtSTCBbfb
                        #pL512fixed6OF8192smprtSTCBbfbSUM0
        )
        preBuffer=75 # selected prebuffer for fwhmvsrlength plotting
    }else if(subtitle == "OPTFILT"){
        methods <- list(pL8192fixed6OF8192smprtSTCBbfb,
                        pL8192fixed6I2R8192smprtSTCBbfb,pL8192fixed6I2RFITTED8192smprtSTCBbfb,
                        pL8192fixed6I2RNOL8192smprtSTCBbfb
        )
    }else if(subtitle == "RSPACE"){
        methods <- list(fixed6OF8192smprtSTCBbfb, fixed6I2R8192smprtSTCBbfb, 
                        fixed6I2RNOL8192smprtSTCBbfb,fixed6OF8192NM50000smprtSTCBbfb)
    }else if(subtitle == "WEIGHTS"){
        methods <- list(weight,weightn)#,weightnOF)
    }else{
        stop("Undefined list of methods\n")
    }
    
    nmethods <- length(methods)
    colors <- sapply(methods, function(x) x$color)
    labs <- sapply(methods, function(x) x$lab)
    points <- sapply(methods, function(x) x$point)
    ltys <- sapply(methods, function(x) x$ltype)
    names <- sapply(methods, function(x) x$name)
    pLengths <- sapply(methods, function(x) x$pLength)
    aliases <- paste("pL",pLengths,"_",names,sep="")
    
    if(plotFWHM_GAINCORRE || plotFWHM_GAINCORRS){
        pdfName <- "fwhmVSenergy_methods"
        outPDF <- paste("./PDFs/",pdfName,"_",subtitle,noiseMat,".pdf",sep="")
        pdf(outPDF,width=7, height=7)
        
        Emin <- 0.1
        Emax <- 9.
        FWmin<-1.7
        FWmax<-10.
        coeffsFile <- "coeffs_polyfit.dat"
        EkeV <- c(0.2,0.5,1,2,3,4,5,6,7,8) # CALIBRATION
        #EkeV <- c(0.2,1,2,3,4,5,6,7,8) # CALIBRATION
        
        # INITIALIZE MATRICES #
        # =======================
        
        ## For CALIBRATION data
        fwhmUNCORR  <- matrix(NA,nrow=length(EkeV), ncol=nmethods) 
        ebiasUNCORR <- matrix(NA,nrow=length(EkeV), ncol=nmethods) 
        fwhmGAINCORRS <- matrix(NA,nrow=length(EkeV), ncol=nmethods) 
        fwhmGAINCORRSErr <- matrix(NA,nrow=length(EkeV), ncol=nmethods)
        fwhmGAINCORRE    <- matrix(NA,nrow=length(EkeV), ncol=nmethods)
        fwhmGAINCORREErr <- matrix(NA,nrow=length(EkeV), ncol=nmethods)
        ebiasGAINCORRE   <- matrix(NA,nrow=length(EkeV), ncol=nmethods)
    
        # READ energy and resolution from DATA  (EkeV)
        # ========================================================
        for (ie in 1:length(EkeV)){
            for (im in 1:nmethods){
                nSamples = methods[[im]]$nSamples
                pulseLength = methods[[im]]$pLength
                ofLength = methods[[im]]$ofLength
                TRIGG <- methods[[im]]$detMethod
                samprateStr <-methods[[im]]$samprateStr
                jitterStr <- methods[[im]]$jitterStr
                bbfbStr <- methods[[im]]$bbfbStr
                lib <- methods[[im]]$lib
                eresolFile1D <- paste("eresol_",nSimPulses,"p_SIRENA",nSamples,
                                      "_pL",pulseLength,"_", EkeV[ie],"keV_",
                                      methods[[im]]$name,".json1D",sep="")
                eresolFile2D <- paste("eresol_",nSimPulses,"p_SIRENA",nSamples,
                                      "_pL",pulseLength,"_", EkeV[ie],"keV_",
                                      methods[[im]]$name,".json2D",sep="")
                if (calib == "curve" || !file.exists(eresolFile2D)){
                    eresolFile <- eresolFile1D
                }else if (calib == "surface"){
                    eresolFile <- eresolFile2D
                }
                if(file.exists(eresolFile)){
                    # use data for selected separation (see initial definitions)
                    cat("Reading file ",eresolFile,"\n")
                    jsondata <- fromJSON(file=eresolFile)
                    idxSep <- which(sapply(jsondata,function(x) x$separation)==separation)
                    fwhmUNCORR[ie,im] <- as.numeric(jsondata[[idxSep]]$fwhmErecons$all)
                    fwhmGAINCORRE[ie,im] <- as.numeric(jsondata[[idxSep]]$fwhmEreal$all)
                    ebiasUNCORR[ie,im] <- as.numeric(jsondata[[idxSep]]$biasErecons$all)
                    ebiasGAINCORRE[ie,im] <- as.numeric(jsondata[[idxSep]]$biasEreal$all)
                    fwhmGAINCORREErr[ie,im] <- as.numeric(jsondata[[idxSep]]$fwhmEreal_err$all)
                        
                }else{
                    warning("Not-existing file:", eresolFile)
                    fwhmUNCORR[ie,im]  <- NaN
                    ebiasUNCORR[ie,im] <- NaN
                    fwhmGAINCORRE[ie,im] <- NaN
                    ebiasGAINCORRE[ie,im] <- NaN
                }
            } #method
        }#energy

        # Read coefficients of poly fit done by polyfit2Bias.R (Erecons vs. Ecalib)
        coeffsTable <- read.table(coeffsFile, header=T)
        
        if(plotFWHM_GAINCORRS){
            # PLOT ENERGY RESOLUTION CORRECTED HERE FROM POLYFIT coefficientes (APPROX??)
            #==============================================================================
            plot(seq(Emin,Emax,length.out=20),seq(FWmin,FWmax,length.out=20),type="n",cex=2,
                 xlab="Input Energy (keV)", ylab="Energy Resolution FWHM (eV)",
                 main=paste("ENERGY RESOLUTION (",array,
                            ")\n(Calculated from error derivation)",sep=""))
            axis(4,labels=FALSE)
            minor.tick(nx=5,ny=5,tick.ratio=0.5)
            #title(main="Comparison of reconstruction methods \n for monochromatic sources")     
            grid(nx=NA,ny=NULL)
                
            for (e in EkeV){abline(v=EkeV,lty=3,col="grey80") }
            for (im in 1:nmethods){
                # correct FWHM for bias # calculations done in keV
                # if I correct *each* Ecalc with the curve Ecalc vs Ereal, then
                #   Ecalc = f(Ereal) =  a0 + a1*(Ereal) + a2*(Ereal)² + a3*(Ereal)³  + a4*(Ereal)⁴
                #  var(Ecal) = [d(f)/d(Ereal)]^2 * var(Ereal)
                #  var(Ereal) = var(Ecal) / [d(f)/d(Ereal)]^2
                #  FWHM(Ereal) = 2.35 * sigma(Ereal)
                    
                # get coefficients for method
                alias <- paste("pL",methods[[im]]$pLength,"_",methods[[im]]$name,sep="")
                methodTable <- coeffsTable[coeffsTable$ALIAS==alias,]
                a0 <- methodTable$a0
                a1 <- methodTable$a1
                a2 <- methodTable$a2
                a3 <- methodTable$a3
                a4 <- methodTable$a4
                
                for (ie in 1:length(EkeV)){
                    EeV <- EkeV[ie] * 1000. 
                    #varEcalc <- (fwhm[ie,im] / 2.35)^2 #eV
                    #df_dEreal <- a1/1E3 + 2.*a2/1E6*EeV + 3.*a3/1E9*EeV^2
                    #varEreal <- varEcalc / (df_dEreal)^2
                    #fwhmSecGain[ie,im] <- 2.35 * sqrt(varEreal) # eV
                    
                    varEcalc <- ((fwhmUNCORR[ie,im]/1000.) / 2.35)^2 #keV
                    df_dEreal <- a1 + 2.*a2*EkeV[ie] + 3.*a3*EkeV[ie]^2 + 4.*a4*EkeV[ie]^3
                    varEreal <- varEcalc / (df_dEreal)^2
                    fwhmGAINCORRS[ie,im] <- 2.35 * sqrt(varEreal)*1000 # eV
                    fwhmGAINCORRSErr[ie,im] <- fwhmGAINCORRS[ie,im]/sqrt(as.numeric(nSimPulses))
                }
                points(EkeV,fwhmGAINCORRS[,im],pch=methods[[im]]$point,col=methods[[im]]$color,
                       type=plttype,lty=methods[[im]]$ltype,pty="s")
                errbar(EkeV,fwhmGAINCORRS[,im], yplus=fwhmGAINCORRS[,im]+fwhmGAINCORRSErr[,im],type="n",cap=0,
                       yminus=fwhmGAINCORRS[,im]-fwhmGAINCORRSErr[,im],add=TRUE,errbar.col=methods[[im]]$color,
                       pty="s")
            }
            
            legend("topleft", legend=labs, col=colors, pch=points, cex=0.6,
                   text.col=colors, bty="n",y.intersp=1., lty=ltys)
        }
            
        if(plotFWHM_GAINCORRE){ 
            # PLOT all ENERGY RESOLUTION CORRECTED at getEresolCurves.py with the polyfit coefficients
            #=====================================================================================
            menuOpts <- c("ADC 0-padding", "ADC short filters", "ADC preBuffer",
                          "I2R 0-padding", "I2R short filters", "I2R preBuffer")
            imeth <- menu(menuOpts, title="Select a method for base plot:")
            menuMethsIds <- list()
            menuMethsIds[[1]] <- grep("OPTFILT8192", aliases, perl=TRUE) # ADC long filter
            menuMethsIds[[2]] <- grep("OPTFILT[1-7]{1}[0-9]{2,}_[^pB]", aliases, perl=TRUE) #ADC short
            menuMethsIds[[3]] <- grep("OPTFILT.*_pB", aliases, perl=TRUE) # ADC preBuffer
            menuMethsIds[[4]] <- grep("I2R8192", aliases, perl=TRUE) # I2R long
            menuMethsIds[[5]] <- grep("I2R[1-7]{1}[0-9]{2,}_[^pB]", aliases, perl=TRUE) #I2R short
            menuMethsIds[[6]] <- grep("I2R.*_pB", aliases, perl=TRUE) #I2R preBuffer
            cat("Color for base plot is:", methods[[menuMethsIds[[imeth]][[1]]]]$color)
            alt <- readline("Alternative plots (y/n):")
            if (alt == "y" || alt == "Y"){
                colAlt <- readline("What color for altenative models?: ")
            }
            
            # Plot secondary plot(s)
            for (iblock in 1:length(menuMethsIds)){
                if (iblock == imeth) next  # this is base plot
                if (length(menuMethsIds[[iblock]]) == 0) next  # no methods to be plotted here
                # Plot dummy box
                plot(seq(Emin,Emax,length.out=20),seq(FWmin,FWmax,length.out=20),type="n",cex=2,
                     xlab="Input Energy (keV)", ylab="Energy Resolution FWHM (eV)",
                     main=paste("ENERGY RESOLUTION (",array,
                                ")\n(Calculated from gain-scale-", calib," calibrated energies)",sep=""))
                minor.tick(nx=5,ny=5,tick.ratio=0.5)
                axis(4,label=FALSE)
                grid(nx=NA,ny=NULL)
                legendCols<-character()
                for (e in EkeV){abline(v=EkeV,lty=3,col="grey80")}
                # tabulated resolutions
                abline(h=2.5, lty=2, col="green")
                text(9,2.5,"HR (>8192)",cex=0.4)
                abline(h=3, lty=2, col="green3")
                text(9,3,"MR (512-8192)",cex=0.4)
                abline(h=7, lty=2, col="darkgreen")
                text(9,7,"LR (256-512)",cex=0.4)
                
                labsLegPlot <- character() 
                colsLegPlot <- character() 
                ptsLegPlot  <- numeric()
                ltysLegPlot  <- numeric()
                
                # Plot 0-SUM method (512) - last methods in list (only if ADC 0-pad is involved)
                # if(iblock == 1 || imeth == 1){
                #     points(EkeV,fwhmGAINCORRE[,nmethods],pch=methods[[nmethods]]$point,
                #             col=methods[[nmethods]]$color,
                #             type=plttype,lty=methods[[nmethods]]$ltype)
                #     errbar(EkeV,fwhmGAINCORRE[,nmethods], type="n",cap=0, 
                #             yplus=fwhmGAINCORRE[,nmethods]+fwhmGAINCORREErr[,nmethods],
                #             yminus=fwhmGAINCORRE[,nmethods]-fwhmGAINCORREErr[,nmethods],
                #             add=TRUE,errbar.col=methods[[nmethods]]$color)
                #     labsLegPlot <- append(methods[[nmethods]]$lab, labsLegPlot)
                #     colsLegPlot <- append(methods[[nmethods]]$color, colsLegPlot)
                #     ptsLegPlot <- append(methods[[nmethods]]$point, ptsLegPlot)
                #     ltysLegPlot <- append(methods[[nmethods]]$ltype, ltysLegPlot)
                # }
                if (alt == "y" || alt == "Y"){
                    # locate 1024NM150000 and plot
                    if(length(inm <- grep("1024NM150000", aliases, perl=TRUE))>0){
                        points(EkeV,fwhmGAINCORRE[,inm],pch=methods[[inm]]$point,
                                col=methods[[inm]]$color,
                                type=plttype,lty=methods[[inm]]$ltype)
                        errbar(EkeV,fwhmGAINCORRE[,inm], type="n",cap=0, 
                                yplus=fwhmGAINCORRE[,inm]+fwhmGAINCORREErr[,inm],
                                yminus=fwhmGAINCORRE[,inm]-fwhmGAINCORREErr[,inm],
                                add=TRUE,errbar.col=methods[[inm]]$color)
                        labsLegPlot <- append(methods[[inm]]$lab, labsLegPlot)
                        colsLegPlot <- append(methods[[inm]]$color, colsLegPlot)
                        ptsLegPlot <- append(methods[[inm]]$point, ptsLegPlot)
                        ltysLegPlot <- append(methods[[inm]]$ltype, ltysLegPlot)
                    }
                    # # locate WEIGHTN1024 and plot
                    if(length(inm <- grep("WEIGHTN1024", aliases, perl=TRUE)) > 0){
                         points(EkeV,fwhmGAINCORRE[,inm],pch=methods[[inm]]$point,
                                col=methods[[inm]]$color,
                                type=plttype,lty=methods[[inm]]$ltype)
                         errbar(EkeV,fwhmGAINCORRE[,inm], type="n",cap=0, 
                                yplus=fwhmGAINCORRE[,inm]+fwhmGAINCORREErr[,inm],
                                yminus=fwhmGAINCORRE[,inm]-fwhmGAINCORREErr[,inm],
                                add=TRUE,errbar.col=methods[[inm]]$color)
                         labsLegPlot <- append(methods[[inm]]$lab, labsLegPlot)
                         colsLegPlot <- append(methods[[inm]]$color, colsLegPlot)
                         ptsLegPlot <- append(methods[[inm]]$point, ptsLegPlot)
                         ltysLegPlot <- append(methods[[inm]]$ltype, ltysLegPlot)
                    }
                
                    # plot methods in alternative blocks
                    for (im in menuMethsIds[[iblock]]){ 
                        if(length(grep("Sum0Filt", aliases[im], perl=TRUE))>0)next
                        points(EkeV,fwhmGAINCORRE[,im],pch=methods[[im]]$point,
                               col=colAlt, type=plttype, lty=methods[[im]]$ltype)
                        errbar(EkeV,fwhmGAINCORRE[,im], type="n",cap=0,
                               yplus=fwhmGAINCORRE[,im]+fwhmGAINCORREErr[,im],
                               yminus=fwhmGAINCORRE[,im]-fwhmGAINCORREErr[,im],
                               add=TRUE,errbar.col=colAlt)
                        labsLegPlot <- append(methods[[im]]$lab, labsLegPlot)
                        #colsLegPlot <- append(methods[[im]]$color, colsLegPlot)
                        colsLegPlot <- append(colAlt, colsLegPlot)
                        ptsLegPlot <- append(methods[[im]]$point, ptsLegPlot)
                        ltysLegPlot <- append(methods[[im]]$ltype, ltysLegPlot)
                        
                        
                        if(length(grep("R",methods[[im]]$name)>0)){
                            legendCols <- append(legendCols,methods[[im]]$color)
                        }
                    }
                }
                # Base plot are methods in menuMethsIds[[imeth]]
                for (im in menuMethsIds[[imeth]]){
                    if(length(grep("Sum0Filt", aliases[im], perl=TRUE))>0)next
                    points(EkeV,fwhmGAINCORRE[,im],pch=methods[[im]]$point,col=methods[[im]]$color,
                           type=plttype,lty=methods[[im]]$ltype)
                    errbar(EkeV,fwhmGAINCORRE[,im], yplus=fwhmGAINCORRE[,im]+fwhmGAINCORREErr[,im],
                           type="n",cap=0, yminus=fwhmGAINCORRE[,im]-fwhmGAINCORREErr[,im],
                           add=TRUE,errbar.col=methods[[im]]$color)
                    labsLegPlot <- append(methods[[im]]$lab, labsLegPlot)
                    colsLegPlot <- append(methods[[im]]$color, colsLegPlot)
                    ptsLegPlot <- append(methods[[im]]$point, ptsLegPlot)
                    ltysLegPlot <- append(methods[[im]]$ltype, ltysLegPlot)
                }
                
                # Print legend of base plot and alternative block
                #legend("topleft", legend=labs, col=colors, pch=points, cex=0.6,
                #       text.col=colors, bty="n",y.intersp=1., lty=ltys)
                legend("topright", legend=labsLegPlot, col=colsLegPlot, pch=ptsLegPlot, cex=0.4,
                              text.col=colsLegPlot, bty="n",y.intersp=1., lty=ltysLegPlot)
            
                RNOL <- expression(R[NOL]*"=" * over(V[0],I)- R[L])
                RFIT <- expression(R[FIT]*"=" * over(V[0],I[fit]+I))
                R    <- expression("R="*R[0]-R[0]*"(" * over(abs( I ),I[bias]+abs( I ))*")")
                Ibias <- expression(I[bias]*"=1.36E-5")
                Ifit <- expression(I[fit]*"=1.46E-5")
            
                #legend("bottomright",legend=c(R, RFIT, RNOL),bty="n",cex=0.8,
                #       text.col=legendCols)
                #legend("bottomleft",legend=c(Ibias,Ifit),bty="n",cex=0.8)
                if (alt == "n" || alt == "N") break
            }
        } # if plotFWHM_GAINCORRE
    } # if plotFWHM_GAINCORRS/plotFWHM_GAINCORRE
        
    if(plotFWHM_rlength){
        
        cat("#\n#Plotting FWHM vs. Record length \n#\n")
        rlens <- c(8192, 4096, 2048, 1024, 512, 256, 128)
        EkeV_rl <- 7
            
        # INITIALIZE MATRICES #
        # =======================
            
        fwhmGAINCORRE    <- matrix(NA,nrow=length(rlens), ncol=nmethods)
        fwhmGAINCORREErr <- matrix(NA,nrow=length(rlens), ncol=nmethods)
        
        # READ recordlength and resolution from DATA
        # ===========================================
        ims1D <- numeric()
        for (il in 1:length(rlens)){
            for (im in 1:nmethods){
                nSamples = methods[[im]]$nSamples
                pulseLength = methods[[im]]$pLength
                ofLength = methods[[im]]$ofLength
                TRIGG <- methods[[im]]$detMethod
                samprateStr <-methods[[im]]$samprateStr
                jitterStr <- methods[[im]]$jitterStr
                bbfbStr <- methods[[im]]$bbfbStr
                lib <- methods[[im]]$lib
                name <- methods[[im]]$name
                alias <- paste("pL",pulseLength,"_",name,sep="")
                
                # for 8192, pulseLength AND ofLength must be 8192
                if (rlens[il] == 8192 && 
                    (pulseLength != rlens[il] || ofLength != rlens[il])) next
                # for other record lengths, any of them must be == rlens[il]
                if (pulseLength != rlens[il] && ofLength != rlens[il]) next
                
                eresolFile1D <- paste("eresol_",nSimPulses,"p_SIRENA",nSamples,
                                     "_pL",pulseLength,"_", EkeV_rl,"keV_",
                                     name,".json1D",sep="")
                eresolFile2D <- paste("eresol_",nSimPulses,"p_SIRENA",nSamples,
                                      "_pL",pulseLength,"_", EkeV_rl,"keV_",
                                      name,".json2D",sep="")
                if (calib == "curve" || !file.exists(eresolFile2D)){
                    eresolFile <- eresolFile1D
                    ims1D <- append(ims1D,im)
                }else if (calib == "surface"){
                    eresolFile <- eresolFile2D
                }
                
                if(file.exists(eresolFile)){
                    # use data for selected separation (see initial definitions)
                    cat("Reading file ",eresolFile,"\n")
                    jsondata <- fromJSON(file=eresolFile)
                    idxSep <- which(sapply(jsondata,function(x) x$separation)==separation)
                    fwhmGAINCORRE[il,im] <- as.numeric(jsondata[[idxSep]]$fwhmEreal$all)
                    fwhmGAINCORREErr[il,im] <- as.numeric(jsondata[[idxSep]]$fwhmEreal_err$all)
                    
                }else{
                    warning("Not-existing file:", eresolFile)
                    fwhmGAINCORRE[ie,im] <- NaN
                    fwhmGAINCORREErr[ie,im] <- NaN
                }
            } #method
        }# record length
        
        # PLOT ENERGY RESOLUTION CORRECTED vs. record length
        #===================================================
        # common part to all plots (1/method to create a movie)
        drawLogPlotBox(xlimits=c(min(100),max(1E4)),ylimits=c(1,50),
                       x2limits=c(min(rlens),max(rlens)),y2limits=c(1.5,7),
                       logxy="xy", xlabel="Record length (samples)", 
                       ylabel="Energy Resolution FWHM (eV)",
                       naxes=c(T,T,T,T))
        title(main=paste("ENERGY RESOLUTION (",array,")\n Surface Correction for jitter",sep=""))
        grid(nx=NA,ny=NULL)
        abline(v=2^(5:13),lty=3,col="grey50") # where pre-calc filters are
        # tabulated resolutions
        abline(h=2.5, lty=3, col="gray4")
        text(0.9E4,2.5,"HR (>8192)",cex=0.4)
        abline(h=3, lty=3, col="gray4")
        text(0.9E4,3,"MR (512-8192)",cex=0.4)
        abline(h=7, lty=3, col="gray4")
        text(0.9E4,7,"LR (256-512)",cex=0.4)
        #abline(v=invalids,lty=3, col="gray4")# invalids
        
        uniqCols <- unique(colors)
        #put last color 'black' for SUM0
        # legend("topright", legend=c("ADC 0-padding","ADC short filter", "ADC preBuffer",
        #                             "I2R 0-padding","I2R short filter", "I2R preBuffer",
        #                             "ADC 1024NM50000", "ADC WEIGHTN1024NM50000", 
        #                             "ADC 0-padding-SUM0"),
        #        col=uniqCols, pch=c(1:(length(uniqCols)-1),0), cex=0.6,
        #        text.col=uniqCols, bty="n",y.intersp=1.)
        legend(x=3000,y=60, legend=c("ADC 0-padding","ADC short filter", 
                                     paste("ADC preBuffer",preBuffer,sep=""),
                                     "I2R 0-padding","I2R short filter", 
                                     paste("I2R preBuffer",preBuffer,sep=""),
                                     "ADC 1024NM150000", 
                                     "ADC WEIGHTN1024NM150000"),
               col=uniqCols, pch=c(1:length(uniqCols)), cex=0.6,
               text.col=uniqCols, bty="n",y.intersp=0.5)
        
        for (im in 1:nmethods){
            pdfName <- paste("fwhmVSrlenght_methods_",sprintf("%02d",im),sep="")
            outPDF <- paste("./PDFs/",pdfName,"_",subtitle,noiseMat,".pdf",sep="")
            #pdf(outPDF,width=7, height=7)
            
            alias <- paste("pL",methods[[im]]$pLength,"_",methods[[im]]$name,sep="")
            pt <- which(methods[[im]]$color == uniqCols)
            #if (im == nmethods) pt=0
            points(rlens,fwhmGAINCORRE[,im],col=methods[[im]]$color,type=plttype,
                   pch=pt)
            errbar(rlens,fwhmGAINCORRE[,im], yplus=fwhmGAINCORRE[,im]+fwhmGAINCORREErr[,im],
                type="n",cap=0, yminus=fwhmGAINCORRE[,im]-fwhmGAINCORREErr[,im],
                add=TRUE,errbar.col=methods[[im]]$color)
            #cat("im=",im," (",alias,")\n")
            if (im %in% ims1D){
                cat("Using 1D calibration for im=",im," (",alias,")\n")
                points(rlens,fwhmGAINCORRE[,im], pch=22,cex=2)
            }
            #if(readline(prompt = "Press <Enter> to continue...(q to quit):") == "q") {
            #    dev.off()
            #    break
            #}
            #dev.off()
            dev.print(pdf, outPDF, width=7, height=7)
        }
    } # plot rlength
    
    if(plotBiasCorrFit){ # pairs of pulses
        cat("#\n#Plotting Energy bias vs. Separation \n#\n")
        # PLOT ENERGY BIAS
        #====================
        #
        sepsStr <- c("00023", "00031", "00042", "00056", "00075", "00101", "00136", 
                     "00182", "00244", "00328", "00439", "00589", "00791", "01061", 
                     "01423", "01908", "02560", "03433", "04605", "40000")
        #sepsStr <- c("00023", "00031", "00075", "00136", "00244", "00589", "01061", 
        #             "01908", "03433", "04605", "40000")
        seps <- as.numeric(sepsStr)
        
        ebiasCorrSec      <- matrix(NA,nrow=length(seps), ncol=length(methods))
        fwhmSecGAINCORRE  <- matrix(NA,nrow=length(seps), ncol=nmethods) # FWHM of Ecorr 
        fwhmPrimGAINCORRE  <- matrix(NA,nrow=length(seps), ncol=nmethods) # FWHM of Ecorr 
        
        EkeVforBias <- 7
        TRIGG <- "_NTRIG"
        
        # READ BIAS & FWHM DATA
        # ======================
        for (is in 1:length(seps)){
            for (im in 1:length(methods)){
                eresolFile <- paste("eresol_",nSimPulses,"p_SIRENA",nSamples,"_pL",pulseLength,
                                    "_",EkeVforBias,"keV_F0F_", methods[[im]]$name,TRIGG,
                                    ".json",sep="")
                if(file.exists(eresolFile)){
                    cat("Reading file ",eresolFile,"\n")
                    jsondata <- fromJSON(file=eresolFile)
                    idxSep <- which(sapply(jsondata,function(x) x$separation)==sepsStr[is])
                    #stopifnot(idxSep>0)
                    if(length(idxSep)>0){
                        ebiasCorrSec[is,im] <- as.numeric(jsondata[[idxSep]]$biasEreal$secondaries)
                        fwhmSecGAINCORRE[is,im] <- as.numeric(jsondata[[idxSep]]$fwhmEreal$secondaries)
                        fwhmPrimGAINCORRE[is,im] <- as.numeric(jsondata[[idxSep]]$fwhmEreal$primaries)
                    }else{
                        ebiasCorrSec[is,im] <- NaN
                        fwhmSecGAINCORRE[is,im] <- NaN
                        fwhmPrimGAINCORRE[is,im] <- NaN
                    }
                    
                }else{
                    warning("Not-existing file:", eresolFile)
                    ebiasCorrSec[is,im] <- NaN
                    fwhmSecGAINCORRE[is,im] <- NaN
                    fwhmPrimGAINCORRE[is,im] <- NaN
                }
                #cat("FWHM=",fwhmSecGAINCORRE[is,im],"\n")
                # cat("FWHM=",fwhmPrimGAINCORRE[is,im],"\n")
            }
        }
        
        drawLogPlotBox(xlimits=c(min(seps),max(seps)),ylimits=c(1E-3,1E4),
                       x2limits=c(min(seps),max(seps)),y2limits=c(1E-3,1E4),
                       logxy="xy", xlabel="Previous Pulse separation (samples)", 
                       ylabel="Energy Bias Gain Corrected (eV)",
                       naxes=c(T,T,T,T))
        title(main=paste(array," AC simulations - 7 keV - secondaries (pL=",pulseLength,")",sep=""))
        grid(nx=NA,ny=NULL)
        abline(v=c(32,64,128,256,512,1024,2048,4096),lty=3,col="grey80")
        abline(v=invalids,lty=2, col="cyan")
        text(80,0.01,"Invalid Events",col="cyan")
        for (im in 1:nmethods){
            points(seps,abs(ebiasCorrSec[,im]),pch=methods[[im]]$point,col=methods[[im]]$color,
                   type=plttype,lty=methods[[im]]$ltype)
        }
        legend("topright",legend=labs, col=colors, pch=points,cex=1,
               text.col=colors, bty="n", lty=ltys)
        
        #Plot also FWHM as a function of Time separation
        # For PRIMARIES...
        drawLogPlotBox(xlimits=c(min(seps),max(seps)),ylimits=c(1.8,4.5),
                       x2limits=c(min(seps),max(seps)),y2limits=c(1.8,4.5),
                       logxy="x", xlabel="Previous Pulse separation (samples)", 
                       ylabel="Energy resolution (FWHM) Gain Corrected (eV)",
                       naxes=c(T,T,T,T))
        title(main=paste(array," AC simulations - 7 keV - PRIMARIES (pL=separation)",sep=""))
        grid(nx=NA,ny=NULL)
        abline(v=c(32,64,128,256,512,1024,2048),lty=3,col="grey80")
        abline(v=invalids,lty=2, col="cyan")
        text(80,1.8,"Invalid Events",col="cyan")
        
        for (im in 1:nmethods){
            points(seps,fwhmPrimGAINCORRE[,im],pch=methods[[im]]$point,col=methods[[im]]$color,
                   type=plttype,lty=methods[[im]]$ltype)
        }
        legend("topright",legend=labs, col=colors, pch=points,cex=1,
               text.col=colors, bty="n", lty=ltys)
        
        # For SECONDARIES...
        drawLogPlotBox(xlimits=c(min(seps),max(seps)),ylimits=c(1.8,4.5),
                       x2limits=c(min(seps),max(seps)),y2limits=c(1.8,4.5),
                       logxy="x", xlabel="Previous Pulse separation (samples)", 
                       ylabel="Energy resolution (FWHM) Gain Corrected (eV)",
                       naxes=c(T,T,T,T))
        title(main=paste(array," AC simulations - 7 keV - SECONDARIES (pL=",pulseLength,")",sep=""))
        grid(nx=NA,ny=NULL)
        abline(v=c(32,64,128,256,512,1024,2048),lty=3,col="grey80")
        abline(v=invalids,lty=2, col="cyan")
        text(80,1.8,"Invalid Events",col="cyan")
        
        for (im in 1:nmethods){
            points(seps,fwhmSecGAINCORRE[,im],pch=methods[[im]]$point,col=methods[[im]]$color,
                   type=plttype,lty=methods[[im]]$ltype)
        }
        legend("topright",legend=labs, col=colors, pch=points,cex=1,
               text.col=colors, bty="n", lty=ltys)
    } # if plotBias
    
    
} else { # plot PERF

    # TIME/PERF
    #==============
    # Results for a simulated file in Intel(R) Xeon(R) CPU   X5550  @ 2.67GHz (jupiter)
    # 100 rows * 40009samples/row /156250 Hz = 25.6 s (200 pulsos; ~ 8 ct/s) 
    simtime <- 25.6
    #             "OPTFILT--LIB-F","OPTFILT--F", "WEIGHT 0(n)", "WEIGHT", "OPTFILT--LIB-T", "OPTFILT--T"  ),
    MFLOP_r0110  = c( 389.8,              394.4,       1020,         1019.5,       385.3,            396.8)
    MFLOP_r0210  = c( 0,                      0,          0,              0,        0,                0)
    MFLOP_r0410  = c( 165.9,              181.4,       945.9,          947.3,       103,              196)
    MFLOP_r0810  = c( 15,                  14.9,       15.8,           16.6,       14.1,             14.4)
    MFLOP_r1010  = c( 0,                      0,          0,              0,          0,              0)
    MFLOP_r2010  = c( 165.8,              181.4,      945.9,          947.3,       103,              196.2)
    MFLOP_r4010  = c( 0,                      0,          0,              0,         0,               0)
    MFLOP_r8010  = c( 165.9,              181.3,      945.9,          947.1,       103,              196.1)
    MFLOP_Xeon   = MFLOP_r0110 + MFLOP_r0210 + MFLOP_r0410 + MFLOP_r0810 + 
                   MFLOP_r1010 + MFLOP_r2010 + MFLOP_r4010 + MFLOP_r8010
    maxMFLOP     = max(MFLOP_Xeon)    
    MFLOPS_Xeon  = MFLOP_Xeon/simtime
    maxMFLOPS     = max(MFLOPS_Xeon)   
    # (5000 rows) = 10000 pulses for CPU and RAM
    CPUTIME      = c( 3.7,                 3.65,     35.44,             60.37,      3.79,            3.71)
    maxCPU       = max(CPUTIME)
    MaxRSS       = c( 5217.636,         5216.388,     5734,        5619.340,    5217.016,         5217.016) 
    maxRAM       = max(MaxRSS)
    

    # Results for a simulated file in Inter(R) Core(TM) i7-2620M CPU @ 2.70 GHz (rhea)
    # 100 rows * 40009samples/row /156250 Hz = 25.6 s (200 pulsos; ~ 8 ct/s) 
    #                 "OPTFILT--LIB", "OPTFILT-", "WEIGHT 0(n)", "WEIGHT"),
    #MFLOP_r0110_rhea = c( 402.9,               407,       1286.7,           840)
    #MFLOP_r0210_rhea = c( 0,                     0,            0,             0)
    #MFLOP_r0410_rhea = c( 165.9,               181,       2968.3,        1189.5)
    #MFLOP_r0810_rhea = c( 99,                107.5,           73,          73.8)
    #MFLOP_r1010_rhea = c( 0,                     0,            0,             0)
    #MFLOP_r2010_rhea = c( 0,                     0,            0,             0)
    #MFLOP_r4010_rhea = c( 0,                     0,            0,             0)
    #MFLOP_r8010_rhea = c( 160.1,              175.1,      2540.5,         971.4)
    #MFLOP_rhea   = MFLOP_r0110_rhea + MFLOP_r0210_rhea + MFLOP_r0410_rhea + MFLOP_r0810_rhea +
    #               MFLOP_r1010_rhea + MFLOP_r2010_rhea + MFLOP_r4010_rhea + MFLOP_r8010_rhea
    
    
    
    timePerf <- data.frame(
        "METHODS"      = c("OPTFILT--LIB-T", "OPTFILT--T", "WEIGHT 0(n)", "WEIGHT"),
        "CPUTIME"      = c( CPUTIME[5]/maxCPU, CPUTIME[6]/maxCPU, CPUTIME[3]/maxCPU, CPUTIME[4]/maxCPU ),
        "MaxRSS[Mb]"   = c( MaxRSS[5]/maxRAM,  MaxRSS[6]/maxRAM,  MaxRSS[3]/maxRAM,  MaxRSS[4]/maxRAM),
        "MFLOPS (8 ct/s)" = c( MFLOPS_Xeon[5]/maxMFLOPS, MFLOPS_Xeon[6]/maxMFLOPS, MFLOPS_Xeon[3]/maxMFLOPS, MFLOPS_Xeon[4]/maxMFLOPS),
                           stringsAsFactors=FALSE)
    mat<-as.matrix(timePerf[,2:4])
    for (i in 1:ncol(mat))   mat[,i]<-mat[,i]/max(mat[,i])
    barplot(height = mat,beside=TRUE,
         col=c("red","indianred","green","darkgreen"),
         legend=timePerf[,1], args.legend = list(x="topleft"))#rownames(timePerf))
    text(x=11.5, y=mat[1,3]+0.05, sprintf("%.2f",MFLOPS_Xeon[5]))
    text(x=12.5, y=mat[2,3]+0.05, sprintf("%.2f",MFLOPS_Xeon[6]))
    text(x=13.5, y=mat[3,3]-0.05, sprintf("%.2f",MFLOPS_Xeon[3]))
    text(x=14.5, y=mat[4,3]+0.05, sprintf("%.2f",MFLOPS_Xeon[4]))
    
} #if plotFWHM or plotPERF

dev.off()
setwd("/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/")
#setwd("/run/media/ceballos/Elements/backup/current/ERESOL/PAIRS/")






