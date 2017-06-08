#
# Plot energy resolution for:
#        a given array (SPA, LPA1, LPA2, LPA3)  
#        @ different energies (0.5,1,2,3,4,6,9) keV
#        for different methods  (fixed and multilib): OPTFILT, WEIGHT, WEIGHTN, I2R, I2RALL, I2RNOL, I2RFITTED
# calculated with pulses separated by 20000 samples      

library(rjson)
library(Hmisc)

# LOAD methods characteristics
load("/home/ceballos/INSTRUMEN/EURECA/ERESOL/methodsForR.Rdat")

plotFWHM_GAINCORRS <- 0   # FWHM vs. Energy
plotFWHM_GAINCORRE <- 1   # FWHM vs. Energy
plotBiasCorrFit    <- 0   # Bias-gainScaleCorrected vs. separation
plotFWHM_rlength   <- 0   # FWHM vs. Record Length     
if(plotFWHM_GAINCORRS || plotFWHM_GAINCORRE){
    pdfName <- "fwhmVSenergy_methods"
}else if (plotBiasCorrFit){
    pdfName <- "biasVSsep_methods"
}else if (plotFWHM_rlength){    
    pdfName <- "fwhmVSrlenght_methods"
}
subtitle <- "ALL"
subtitle <- "FIXEDLIB"
subtitle <- "FIXEDLIBOF"
#subtitle <- "WEIGHTS"
#subtitle <- "MULTILIB"
#subtitle <- "PERF"
#subtitle <- "OPTFILT"
#subtitle <- "RSPACE"
#subtitle <- "SPIE2016" # SPIE2016/SPIE2016PP (to plot also PP points)
plttype="b"
useGainCorr <- "all" # or "SEC" (pulses used for Gain Curve)
pulsesCateg <- c("primaries", "secondaries","all")
pulsesCateg <- c("all")
par(pty="s")
separation <- "40000"

array <- "LPA2shunt"
nSimPulses <- "20000"
nIntervals <- "150000"
nSamples <- "4096" # samples for the noise (also pulseLength in library)
pulseLength <- "4096" # pulse(record) length
invalids <- 700

outPDF <- paste("./PDFs/",pdfName,"_",subtitle,"_noiseMat",nIntervals,".pdf",sep="")
setwd(paste("/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresol",array,sep=""))
pdf(outPDF,width=7, height=7)

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
                        weight,weightnOF)
    }else if(subtitle == "FIXEDLIB"){
        methods <- list(fixed1OF, fixed1OF_I2R, fixed1OF_I2RNOL, fixed1OF_I2RFITTED,fixed1OFNM,
                        fixed1, fixed1_I2R, fixed1_I2RNOL, fixed1_I2RFITTED, weightnOF,  weight)        
        #methods <- list(fixed1, fixed1OFNM, fixed1_I2R, fixed1_I2RNOL, fixed1_I2RFITTED)        
    }else if(subtitle == "OPTFILT"){
        methods <- list(fixed1,fixed1OF,multi)
    }else if(subtitle == "RSPACE"){
        methods <- list(fixed1_I2R,fixed1_I2RALL, fixed1_I2RNOL, fixed1_I2RFITTED)
    }else if(subtitle == "WEIGHTS"){
        methods <- list(weight,weightn)#,weightnOF)
    }
    
    nmethods <- length(methods)
    colors <- sapply(methods, function(x) x$color)
    labs <- sapply(methods, function(x) x$lab)
    points <- sapply(methods, function(x) x$point)
    ltys <- sapply(methods, function(x) x$ltype)

    
    if(plotFWHM_GAINCORRE || plotFWHM_GAINCORRS){
        Emin <- 0.1
        Emax <- 9.
        FWmin<-1.7
        FWmax<-2.3
        coeffsFile <- paste("coeffs_polyfit_",nSamples,"_",useGainCorr,".dat",sep="")
        EkeV <- c(0.2,0.5,1,2,3,4,5,6,7,8) # CALIBRATION
        
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
            TRIGG = "_NTRIG"
            for (im in 1:nmethods){
                eresolFile <- paste("eresol_",nSimPulses,"p_SIRENA",nSamples,"_pL",pulseLength,"_", EkeV[ie],"keV_F0F_", 
                                    methods[[im]]$name,TRIGG,".json",sep="")
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
                 main=paste("ENERGY RESOLUTION: GAINCORRS (with ",useGainCorr," pulses) \n",
                            array," - ",nSamples,sep=""),
                 sub="(Calibration Points marked)",pty="s")
            axis(4,labels=FALSE)
            minor.tick(nx=5,ny=5,tick.ratio=0.5)
            #title(main="Comparison of reconstruction methods \n for monochromatic sources")     
            grid(nx=NA,ny=NULL)
            if(subtitle=="SPIE2016PP"){
                points(8,2.19,pch=8,col=methods[[1]]$color) # ADC fixed1
                points(8,2.18,pch=8,col=methods[[2]]$color) # RNOL fixed1
                points(8,2.15,pch=8,col=methods[[3]]$color) # RFITTED fixed1
                points(8,2.11,pch=8,col=methods[[4]]$color) # WEIGHT
                points(8,2.03,pch=8,col=methods[[5]]$color) # WEIGHTN
                
                points(7,2.122,pch=8,col=methods[[1]]$color) # ADC fixed1
                points(7,2.11,pch=8,col=methods[[2]]$color) # RNOL fixed1
                points(7,2.09,pch=8,col=methods[[3]]$color) # RFITTED fixed1
                points(7,2.07,pch=8,col=methods[[4]]$color) # WEIGHT
                points(7,1.99,pch=8,col=methods[[5]]$color) # WEIGHTN
                
                points(6,2.06,pch=8,col=methods[[1]]$color) # ADC fixed1
                points(6,2.02,pch=8,col=methods[[4]]$color) # WEIGHT
                points(6,1.95,pch=8,col=methods[[5]]$color) # WEIGHTN 
                
                points(4,1.97,pch=8,col=methods[[1]]$color) # ADC fixed1
                points(4,1.94,pch=8,col=methods[[4]]$color) # WEIGHT
                points(4,1.91,pch=8,col=methods[[5]]$color) # WEIGHTN 
                
                points(2,1.88,pch=8,col=methods[[1]]$color) # ADC fixed1
                points(2,1.86,pch=8,col=methods[[4]]$color) # WEIGHT
                points(2,1.82,pch=8,col=methods[[5]]$color) # WEIGHTN 
            }
                
            for (e in EkeV){abline(v=EkeV,lty=3,col="grey80") }
            for (im in 1:nmethods){
                # correct FWHM for bias # calculations done in keV
                # if I correct *each* Ecalc with the curve Ecalc vs Ereal, then
                #   Ecalc = f(Ereal) =  a0 + a1*(Ereal) + a2*(Ereal)² + a3*(Ereal)³  + a4*(Ereal)⁴
                #  var(Ecal) = [d(f)/d(Ereal)]^2 * var(Ereal)
                #  var(Ereal) = var(Ecal) / [d(f)/d(Ereal)]^2
                #  FWHM(Ereal) = 2.35 * sigma(Ereal)
                    
                # get coefficients for method
                alias <- methods[[im]]$name
                #alias <- gsub("OF","",alias) # coeffs are equal for OFLib=yes & OFLib=no methods
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
                   text.col=colors, bty="n",y.intersp=2, lty=ltys)
        }
            
        if(plotFWHM_GAINCORRE){ 
            # PLOT all ENERGY RESOLUTION CORRECTED at getEresolCurves.py with the polyfit coefficients
            #=====================================================================================
            
            plot(seq(Emin,Emax,length.out=20),seq(FWmin,FWmax,length.out=20),type="n",cex=2,
                 xlab="Input Energy (keV)", ylab="Energy Resolution FWHM (eV)",
                 main=paste("ENERGY RESOLUTION: GAINCORRE (with ",
                            useGainCorr," pulses) \n",array," - ",nSamples,sep=""),
                 sub="(Calibration Points marked)")
            minor.tick(nx=5,ny=5,tick.ratio=0.5)
            axis(4,label=FALSE)
            grid(nx=NA,ny=NULL)
            for (e in EkeV){abline(v=EkeV,lty=3,col="grey80")}
            for (im in 1:nmethods){
                points(EkeV,fwhmGAINCORRE[,im],pch=methods[[im]]$point,col=methods[[im]]$color,
                       type=plttype,lty=methods[[im]]$ltype)
                errbar(EkeV,fwhmGAINCORRE[,im], yplus=fwhmGAINCORRE[,im]+fwhmGAINCORREErr[,im],type="n",cap=0,
                       yminus=fwhmGAINCORRE[,im]-fwhmGAINCORREErr[,im],add=TRUE,errbar.col=methods[[im]]$color)
                
            }
            if(subtitle=="SPIE2016PP"){
                points(8,2.19,pch=8,col=methods[[1]]$color) # ADC fixed1
                points(8,2.18,pch=8,col=methods[[2]]$color) # RNOL fixed1
                points(8,2.15,pch=8,col=methods[[3]]$color) # RFITTED fixed1
                points(8,2.11,pch=8,col=methods[[4]]$color) # WEIGHT
                points(8,2.03,pch=8,col=methods[[5]]$color) # WEIGHTN
                
                points(7,2.122,pch=8,col=methods[[1]]$color) # ADC fixed1
                points(7,2.11,pch=8,col=methods[[2]]$color) # RNOL fixed1
                points(7,2.09,pch=8,col=methods[[3]]$color) # RFITTED fixed1
                points(7,2.07,pch=8,col=methods[[4]]$color) # WEIGHT
                points(7,1.99,pch=8,col=methods[[5]]$color) # WEIGHTN
                
                points(6,2.06,pch=8,col=methods[[1]]$color) # ADC fixed1
                points(6,2.02,pch=8,col=methods[[4]]$color) # WEIGHT
                points(6,1.95,pch=8,col=methods[[5]]$color) # WEIGHTN 
                
                points(4,1.97,pch=8,col=methods[[1]]$color) # ADC fixed1
                points(4,1.94,pch=8,col=methods[[4]]$color) # WEIGHT
                points(4,1.91,pch=8,col=methods[[5]]$color) # WEIGHTN 
                
                points(2,1.88,pch=8,col=methods[[1]]$color) # ADC fixed1
                points(2,1.86,pch=8,col=methods[[4]]$color) # WEIGHT
                points(2,1.82,pch=8,col=methods[[5]]$color) # WEIGHTN 
            }
            legend("topleft", legend=labs, col=colors, pch=points, cex=0.6,
                   text.col=colors, bty="n",y.intersp=2, lty=ltys)
        }
    } # if plotFWHM_GAINCORRS/plotFWHM_GAINCORRE
        
    if(plotFWHM_rlength){
        cat("#\n#Plotting FWHM vs. Record length \n#\n")
        rlens <- c(4096, 2048, 1024, 750, 512, 400, 256, 200, 128, 90, 64, 45, 32)
        EkeV_rl <- 7
        EkeV_rl <- 1
        TRIGG <- "_NTRIG"
            
        # INITIALIZE MATRICES #
        # =======================
            
        fwhmGAINCORRE    <- matrix(NA,nrow=length(rlens), ncol=nmethods) 
        fwhmGAINCORREErr <- matrix(NA,nrow=length(rlens), ncol=nmethods) 
        ebiasCorr        <- matrix(NA,nrow=length(rlens), ncol=nmethods)
        ebiasRecons      <- matrix(NA,nrow=length(rlens), ncol=nmethods)
        
        # READ recordlength and resolution from DATA
        # ===========================================
        for (il in 1:length(rlens)){
            for (im in 1:nmethods){
                eresolFile <- paste("eresol_",nSimPulses,"p_SIRENA",pulseLength,"_pL",rlens[il],
                                "_",EkeV_rl,"keV_F0F_", methods[[im]]$name,TRIGG,".json",sep="")
                isWEIGHTN <- length(i<-grep("WEIGHTN",methods[[im]]$name))
                isWEIGHT <- ("multilib_WEIGHT"==methods[[im]]$name)
                if(isWEIGHT) next
                if(isWEIGHTN){
                    eresolFile <- paste("eresol_","5000","p_SIRENA",pulseLength,"_pL",rlens[il],
                                        "_",EkeV_rl,"keV_F0F_", methods[[im]]$name,TRIGG,".json",sep="")
                }
                    
                if(file.exists(eresolFile)){
                    cat("Reading file ",eresolFile,"\n")
                    jsondata <- fromJSON(file=eresolFile)
                    idxSep <- which(sapply(jsondata,function(x) x$separation)==separation)
                    fwhmGAINCORRE[il,im] <- as.numeric(jsondata[[idxSep]]$fwhmEreal$all)
                    fwhmGAINCORREErr[il,im] <- as.numeric(jsondata[[idxSep]]$fwhmEreal_err$all)
                    ebiasCorr[il,im] <- as.numeric(jsondata[[idxSep]]$biasEreal$all)
                    ebiasRecons[il,im] <- as.numeric(jsondata[[idxSep]]$biasErecons$all)
                }else{
                    warning("Non-existing file:", eresolFile,"\n")
                    fwhmGAINCORRE[il,im] <- NaN
                    ebiasCorr[il,im]  <- NaN
                }
            }
        }
        
        # PLOT ENERGY RESOLUTION CORRECTED vs. record length
        #===================================================
        drawLogPlotBox(xlimits=c(min(rlens),max(rlens)),ylimits=c(2,7),
                       x2limits=c(min(rlens),max(rlens)),y2limits=c(2,7),
                       logxy="x", xlabel="Record length (samples)", 
                       ylabel="Energy Resolution FWHM (eV)",
                       naxes=c(T,T,T,T))
        title(main=paste("AC simulations - ",EkeV_rl," keV - ", separation," samples separation (All)\n",
                         array,sep=""))
        grid(nx=NA,ny=NULL)
        abline(v=2^(5:12),lty=3,col="grey50") # where pre-calc filters are
        abline(v=invalids,lty=2, col="cyan")
        text(80,2.5,"Invalid Events",col="cyan")
        
        for (im in 1:nmethods){
            isWEIGHT <- ("multilib_WEIGHT"==methods[[im]]$name)
            if(isWEIGHT) next
            points(rlens,fwhmGAINCORRE[,im],pch=methods[[im]]$point,col=methods[[im]]$color,
                   type=plttype,lty=methods[[im]]$ltype)
            errbar(rlens,fwhmGAINCORRE[,im], yplus=fwhmGAINCORRE[,im]+fwhmGAINCORREErr[,im],
                   type="n",cap=0, yminus=fwhmGAINCORRE[,im]-fwhmGAINCORREErr[,im],
                   add=TRUE,errbar.col=methods[[im]]$color)
            
        }
        legend("topright", legend=labs, col=colors, pch=points, cex=0.5,
               text.col=colors, bty="n",y.intersp=2, lty=ltys)
        
        # Plot also Erecons vs Rlength for OPTFILT and OPTFILTNM and WEIGHTNOF
        # ======================================================================
        cat("#\n#Plotting FWHM vs. Record length \n#\n")
        drawLogPlotBox(xlimits=c(min(rlens),max(rlens)),ylimits=c(EkeV_rl-0.5,EkeV_rl+0.5),
                       x2limits=c(min(rlens),max(rlens)),y2limits=c(EkeV_rl-0.5,EkeV_rl+0.5),
                       logxy="x", xlabel="Record length (samples)", 
                       ylabel="Gain scale-corrected Energy (keV)",
                       naxes=c(T,T,T,T))
        title(main=paste("AC simulations - 7 keV - ", separation," samples separation (All)\n",
                         array,sep=""))
        grid(nx=NA,ny=NULL)
        abline(v=2^(5:12),lty=3,col="grey50") # where pre-calc filters are
        abline(v=invalids,lty=2, col="cyan")
        text(80,2.5,"Invalid Events",col="cyan")
        ims <- c()
        for (im in 1:nmethods){
            isOPTFILT <- length(i<-grep("OPTFILT",methods[[im]]$name))
            isWEIGHTNOF <- length(i<-grep("OF_WEIGHTN",methods[[im]]$name))
            
            if(!isOPTFILT && !isWEIGHTNOF) next
            ims <- append(ims,im)
            points(rlens,ebiasCorr[,im]/1000.+EkeV_rl,pch=methods[[im]]$point,col=methods[[im]]$color,
                   type=plttype,lty=methods[[im]]$ltype)
            #points(rlens,ebiasRecons[,im]/1000.+EkeV_rl,pch=methods[[im]]$point,col=methods[[im]]$color,
            #       type=plttype,lty=methods[[im]]$ltype)
        }
        legend("topright", legend=labs[ims], col=colors[ims], pch=points[ims], cex=0.6,
               text.col=colors[ims], bty="n",y.intersp=2, lty=ltys[ims])
    }
        

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
        legend("topright",legend=labs, col=colors, pch=points,cex=0.6,
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
        legend("topright",legend=labs, col=colors, pch=points,cex=0.6,
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
        legend("topright",legend=labs, col=colors, pch=points,cex=0.6,
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






