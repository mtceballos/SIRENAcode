# plot histogram from reconstructed energies
#===========================================

# READ data and locate multiple-pulses records
#---------------------------------------------
rm(list=ls())
library("MASS")
library(stats)
library(Hmisc)
noise <- "2"
file <- "2"
hostname <- Sys.getenv("HOSTNAME")
if (length(grep("jupiter",hostname))>0){
    HEADAS="/dataj4/software/64/heasoft/x86_64-unknown-linux-gnu-libc2.17" #for jupiter
}else{
    HEADAS="/home/ceballos/sw/heasoft/x86_64-pc-linux-gnu-libc2.28" # for rhea
    library(RcppFaddeeva)
}
cat("Using HEADAS=", HEADAS,"\n")
GS <- function (x,C,mn,sg){
    g <- C * exp(-(x-mn)**2/(2 * sg**2)) 
    return(g)
}
ratio_th <- 5 # threshold for G2/G1 ratio to select Ka2 pulses

setwd("/dataj6/ceballos/INSTRUMEN/EURECA/realData/h5Files")
dataname <- paste("evt_file",file,"ph_lib",noise,"file",file,"phKas_OPTFILT8192.fits+1",sep="") #reconstructed with double-ka-lib @ 5.895
colname <- "'SIGNAL, PH_ID, GRADE1, GRADE2'" # PH_ID identifies record
command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                 "fdump wrap=yes infile=",dataname," columns=",colname, " rows='-' prhead=no ",
                 "showcol=yes showunit=no showrow=no outfile=pulse.txt clobber=yes", sep="")
system(command)
dataKas <- read.table(file='pulse.txt', header=TRUE )
n_occur <- data.frame(table(dataKas$PH_ID)) # table with PH_ID and times they occur
dataKas_single <-dataKas[!dataKas$PH_ID %in% dataKas[duplicated(dataKas$PH_ID),]$PH_ID,] # records with ONLY 1 pulse
file.remove("pulse.txt")

# READ data from HR file for FWHM analysis
#-------------------------------------------
dataname <- paste("evt_file",file,"ph_lib",noise,"file",file,"phKas_OPTFILT8192_HR.fits+1",sep="") #reconstructed with double-ka-lib @ 5.895
colname <- "'SIGNAL, PH_ID, GRADE1, GRADE2'" # PH_ID identifies record
command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                 "fdump wrap=yes infile=",dataname," columns=",colname, " rows='-' prhead=no ",
                 "showcol=yes showunit=no showrow=no outfile=pulse.txt clobber=yes", sep="")
system(command)
dataKas_HR <- read.table(file='pulse.txt', header=TRUE )
file.remove("pulse.txt")

# Plot data (ka1,ka2)   (ka1=5.89875 keV (16.2%); ka2=5.88765 (8.2%))
#-----------------------------------------------------------------------------
histo<-hist(dataKas_HR$SIGNAL, xlab='Reconstructed PH', breaks = 200,
            main="Reconstructed PHs")
dataKas_HR <- dataKas_HR[dataKas_HR$SIGNAL>5.80 & dataKas_HR$SIGNAL<5.95,]
histo<-hist(dataKas_HR$SIGNAL, xlab='Reconstructed PH', breaks = 100,
            main="Reconstructed PHs")
rug(dataKas_HR$SIGNAL)

# fit 2 Gaussians over *density* histogram
#-----------------------------------------
histo<-hist(dataKas_HR$SIGNAL, xlab='Reconstructed PH', breaks=100,freq = FALSE,
            main="Ka1+Ka2 PH", xlim=c(5.85,5.95))
x <- histo$mids
y <- histo$density
df <- data.frame(x, y)
fit2G <- nls(y ~ (C1 * exp(-(x-mean1)**2/(2 * sigma1**2)) +
                      C2 * exp(-(x-mean2)**2/(2 * sigma2**2))), data=df,
                      start=list(C1=56, mean1=5.89, sigma1=0.01,
                      C2=90, mean2=5.90, sigma2=0.01), algorithm="port")  
newx <- seq(5.80,5.95,length.out = 500)
fit2G_prediction <- predict(fit2G, newdata=list(x=newx))
lines(newx, fit2G_prediction, col="red")
CG2_1  <- coef(fit2G)[1]
meanG2_1 <- coef(fit2G)[2] # keV 
sigG2_1  <- coef(fit2G)[3]
fwhmG2_1 <- sigG2_1*2.35*1000 #eV
CG2_2  <- coef(fit2G)[4]
meanG2_2 <- coef(fit2G)[5] # keV
sigG2_2  <- coef(fit2G)[6]
fwhmG2_2 <- sigG2_2*2.35*1000 #eV

GS1 <- GS(newx,C=CG2_1,mn=meanG2_1, sg=sigG2_1) # Gaussian1
GS2 <- GS(newx,C=CG2_2,mn=meanG2_2, sg=sigG2_2) # Gaussian2

text(5.86,90, "Double Gaussian", col="red")
text(5.86,80, paste("Mean(Ka1)=",sprintf("%.3f",meanG2_1),"a.u"), col="blue")
text(5.86,70, paste("FWHM(Ka1)=",sprintf("%.3f",fwhmG2_1), "ma.u"), col="blue")
text(5.86,60, paste("Mean(Ka2)=",sprintf("%.3f",meanG2_2),"a.u"), col="green")
text(5.86,50, paste("FWHM(Ka2)=",sprintf("%.3f",fwhmG2_2), "ma.u"), col="green")

lines(newx,GS1, col="blue")
lines(newx,GS2, col="green")
abline(v=meanG2_1, col="blue")
abline(v=meanG2_2, col="green")
ratio <- GS2/GS1
lines(newx,ratio, col="orange")
abline(h=ratio_th, col="orange",lty=2) # level 5 in ratio GS2/GS1
text(5.86,20, "ratio = GS2/GS1", col="orange")
PH1 <- newx[which.max(ratio >= ratio_th)]
PH2 <- newx[max(which(ratio >= ratio_th))]
cat("Ka2 PHs in [",PH1, ",",PH2,"] a.u. \n",sep="")

PHID_Ka2        <- dataKas_HR[dataKas_HR$SIGNAL>=PH1 & dataKas_HR$SIGNAL <=PH2, ]$PH_ID #correct energy
PHID_Ka2_single <- PHID_Ka2[PHID_Ka2 %in% dataKas_single$PH_ID] # not in multiple-pulse records
PHID_noKa2      <- dataKas$PH_ID[!dataKas$PH_ID %in% PHID_Ka2_single]



#=======================================
# FINE TUNNING
#=======================================
fileKa2 <- paste("file",file,"ph_Ka2_noise",noise,"_ratio",ratio_th,".fits", sep="")
finetun <- "Y"
if (finetun == "Y" | finetun == "y"){
    # select frome initial FITS file those pulses in meanG2_2+/-sigG2_2 and create a new library
    nphs <- length(PHID_noKa2)
    expr <- paste("'PH_ID != ", PHID_noKa2[1],"'",sep="")
    command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                     "fselect infile=file",file,"ph_Kas.fits+8 outfile=",
                     fileKa2," clobber=yes expr=",expr, sep="")
    cat(command,"\n")
    system(command)
    iph <- 2
    while (iph <= nphs){
        expr <- "'"
        iiph <- 1
        while (iiph < 20 && iph < nphs){
            expr <- paste(expr,"PH_ID != ", PHID_noKa2[iph]," && ", sep="")
            #cat("Adding:", iiph, "\n")
            iiph <- iiph + 1
            iph <- iph + 1
        } 
        #cat("Adding iph=", iph,"\n")
        expr <- paste(expr,"PH_ID != ", PHID_noKa2[iph],"'", sep="")
        iph <- iph + 1
        command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                         "fselect infile=",fileKa2,"+8 outfile=pp.fits clobber=yes expr=",expr, sep="")
        cat(command,"\n")
        cat("iph=", iph, "/",nphs, "\n")
        system(command)
        file.copy("pp.fits", fileKa2, overwrite = TRUE)
    }
    cat("Finshed selection of Ka2 events\n")
}

#
# ========================== NEW LIBRARY ==================================================
#
# =========================================================================================
# Reconstruct again fileX_Kas.fits with new Ka2 library and get new histograms
# =========================================================================================
newlib <- "N"
if (newlib == "Y" | newlib == "y"){
    # READ data from HR file for FWHM analysis
    #-------------------------------------------
    dataname <- paste("evt_file",file,"ph_lib",noise,"file",file,"phKa2_ratio",ratio_th,"_OPTFILT8192_HR.fits+1",sep="") #reconstructed with Ka2-lib @ 5.89875
    #dataname <- paste("evt_file",file,"phid_lib",noise,"file",file,"phKa2_ratio",ratio_th,"_OPTFILT8192_HR.fits+1",sep="")
    dataname <- paste("evt_file",file,"ph_lib",noise,"file2phKa2_ratio",ratio_th,"_OPTFILT8192_HR.fits+1",sep="") #reconstructed with Ka2-lib @ 5.89875    
    dataname <- paste("evt_file",file,"ph_lib",noise,"file2phKa2_ratio",ratio_th,"_OPTFILT4096_HR.fits+1",sep="") #reconstructed with Ka2-lib @ 5.89875    
    dataname <- paste("evt_file",file,"ph_lib",noise,"file2phKa2_ratio",ratio_th,"_OPTFILT2048_HR.fits+1",sep="") #reconstructed with Ka2-lib @ 5.89875    
    dataname <- paste("evt_file",file,"ph_lib",noise,"file2phKa2_ratio",ratio_th,"_OPTFILT1024_HR.fits+1",sep="") #reconstructed with Ka2-lib @ 5.89875    
    dataname <- paste("evt_file",file,"ph_lib",noise,"file2phKa2_ratio",ratio_th,"_pL4096_HR.fits+1",sep="") #reconstructed with Ka2-lib @ 5.89875    
    dataname <- paste("evt_file",file,"ph_lib",noise,"file2phKa2_ratio",ratio_th,"_pL2048_HR.fits+1",sep="") #reconstructed with Ka2-lib @ 5.89875    
    dataname <- paste("evt_file",file,"ph_lib",noise,"file2phKa2_ratio",ratio_th,"_pL1024_HR.fits+1",sep="") #reconstructed with Ka2-lib @ 5.89875    
    dataname <- paste("evt_file",file,"ph_lib",noise,"file2phKa2_ratio",ratio_th,"_OPTFILT8192pB75_HR.fits+1",sep="") #reconstructed with Ka2-lib @ 5.89875    
    dataname <- paste("evt_file",file,"ph_lib",noise,"file2phKa2_ratio",ratio_th,"_OPTFILT4096pB75_HR.fits+1",sep="") #reconstructed with Ka2-lib @ 5.89875    
    dataname <- paste("evt_file",file,"ph_lib",noise,"file2phKa2_ratio",ratio_th,"_OPTFILT2048pB75_HR.fits+1",sep="") #reconstructed with Ka2-lib @ 5.89875    
    #dataname <- paste("evt_file",file,"_lib",noise,"file2phKa2_ratio",ratio_th,"_OPTFILT1024pB75_HR.fits+1",sep="") #reconstructed with Ka2-lib @ 5.89875    
    colname <- "'SIGNAL, PH_ID, GRADE1, GRADE2'" # PH_ID identifies record
    command <- paste("export HEADAS=",HEADAS,";. $HEADAS/headas-init.sh;",
                     "fdump wrap=yes infile=",dataname," columns=",colname, 
                     " rows='-' prhead=no ",
                     "showcol=yes showunit=no showrow=no outfile=pulse.txt clobber=yes", sep="")
    system(command)
    data_HR <- read.table(file='pulse.txt', header=TRUE )
    dataKas_HR <- data_HR[data_HR$SIGNAL>5.80 & data_HR$SIGNAL<5.95,]
    dataKb_HR <- data_HR[data_HR$SIGNAL>6.1 & data_HR$SIGNAL<6.7,]
    
    # fit 1 GAUSS for Kb
    #--------------------
    Kb="N"
    if (Kb == "Y" | Kb == "y"){
        histo0<-hist(dataKb_HR$SIGNAL, xlab='Reconstructed PH', breaks=500, freq = FALSE,
                    main="Kb PH (Lib Ka2)", xlim=c(6.1,6.4))
        x <- histo0$mids
        y <- histo0$density
        df <- data.frame(x, y)
        fit1G_Kb <- nls(y ~ (C1 * exp(-(x-mean1)**2/(2 * sigma1**2))), data=df,
                         start=list(C1=25, mean1=6.25, sigma1=0.01), algorithm="port")  
        newx <- seq(6.1,6.40,length.out = 500)
        fit1G_Kb_prediction <- predict(fit1G_Kb, newdata=list(x=newx))
        lines(newx, fit1G_Kb_prediction, col="red")
        CG1  <- coef(fit1G_Kb)[1]
        meanG1 <- coef(fit1G_Kb)[2] # keV 
        sigG1  <- coef(fit1G_Kb)[3]
        fwhmG1 <- sigG1*2.35*1000 #eV
        
        dataKbCalib <- dataKb_HR$SIGNAL / meanG1*6.49045
        histo1<-hist(dataKbCalib, xlab='Calibrated energy', breaks=500, freq = FALSE,
                    main="Kb (keV) (Lib Ka2)", xlim=c(6.4,6.6))
        x <- histo1$mids
        y <- histo1$density
        df <- data.frame(x, y)
        fit1G_Kb <- nls(y ~ (C1 * exp(-(x-mean1)**2/(2 * sigma1**2))), data=df,
                        start=list(C1=25, mean1=6.49, sigma1=0.01), algorithm="port")  
        newx <- seq(6.4,6.6,length.out = 500)
        fit1G_Kb_prediction <- predict(fit1G_Kb, newdata=list(x=newx))
        lines(newx, fit1G_Kb_prediction, col="red")
        CG1  <- coef(fit1G_Kb)[1]
        meanG1 <- coef(fit1G_Kb)[2] # keV 
        sigG1  <- coef(fit1G_Kb)[3]
        fwhmG1 <- sigG1*2.35*1000 #eV
        xtext <- 6.44
        text(xtext,60, "Single Gaussian", col="red")
        text(xtext,50, paste("Mean(Kb)=",sprintf("%.3f",meanG1),"keV"), col="red")
        text(xtext,40, paste("FWHM(Kb)=",sprintf("%.3f",fwhmG1), "eV"), col="red")
    }    

    
    # fit 2 Gaussians to calibrate
    #------------------------------
    histo2<-hist(dataKas_HR$SIGNAL, xlab='Reconstructed PH', breaks=100,freq = FALSE,
                 main="Ka1+Ka2 PH (Lib Ka2)", xlim=c(5.87,5.93))
    x <- histo2$mids
    y <- histo2$density
    df <- data.frame(x, y)
    fit2G_Ka2 <- nls(y ~ (C1 * exp(-(x-mean1)**2/(2 * sigma1**2)) +
                              C2 * exp(-(x-mean2)**2/(2 * sigma2**2))), data=df,
                     start=list(C1=25, mean1=5.895, sigma1=0.01,
                                C2=30, mean2=5.905, sigma2=0.01), algorithm="port")  
    newx <- seq(5.80,5.95,length.out = 500)
    fit2G_Ka2_prediction <- predict(fit2G_Ka2, newdata=list(x=newx))
    lines(newx, fit2G_Ka2_prediction, col="red")
    CG2_1  <- coef(fit2G_Ka2)[1]
    meanG2_1 <- coef(fit2G_Ka2)[2] # keV 
    sigG2_1  <- coef(fit2G_Ka2)[3]
    err_sigG2_1 <- coef(summary(fit2G_Ka2))[3,2]
    fwhmG2_1 <- sigG2_1*2.35*1000 #eV
    CG2_2  <- coef(fit2G_Ka2)[4]
    meanG2_2 <- coef(fit2G_Ka2)[5] # keV
    sigG2_2  <- coef(fit2G_Ka2)[6]
    err_sigG2_2 <- coef(summary(fit2G_Ka2))[6,2]
    fwhmG2_2 <- sigG2_2*2.35*1000 #eV
    GS1 <- GS(newx,C=CG2_1,mn=meanG2_1, sg=sigG2_1) # Gaussian1
    GS2 <- GS(newx,C=CG2_2,mn=meanG2_2, sg=sigG2_2) # Gaussian2
    
    xtext <- 5.88
    text(xtext,90, "Double Gaussian", col="red")
    text(xtext,80, paste("Mean(Ka1)=",sprintf("%.3f",meanG2_1),"a.u."), col="blue")
    text(xtext,70, paste("FWHM(Ka1)=",sprintf("%.3f",fwhmG2_1), "ma.u."), col="blue")
    text(xtext,60, paste("Mean(Ka2)=",sprintf("%.3f",meanG2_2),"a.u."), col="green")
    text(xtext,50, paste("FWHM(Ka2)=",sprintf("%.3f",fwhmG2_2), "ma.u."), col="green")
    lines(newx,GS1, col="blue")
    lines(newx,GS2, col="green")
    
    y <- c(5.88765,5.89875)
    x <- c(min(meanG2_1,meanG2_2), max(meanG2_1,meanG2_2))
    calib_fit <- lm(y ~ x)
    a0 <- calib_fit$coefficients[1]
    b  <- calib_fit$coefficients[2]
    Energies_Kas = a0 + b * dataKas_HR$SIGNAL
    
    # Plot calibrated histogram and fit again
    # -----------------------------------------
    fG <- FALSE
    if (fG){
        # fit 2 Gaussians
        #----------------
        histo3<-hist(Energies_Kas, xlab='Calibrated energies', breaks=100,freq = FALSE,
                    main=paste("Ka1+Ka2 energies (File",file,", Noise",noise,", Lib file2 Ka2)",sep=""), 
                    xlim=c(5.85,5.95))
        x <- histo3$mids
        y <- histo3$density
        df <- data.frame(x, y)
        fit2G_Ka2 <- nls(y ~ (C1 * exp(-(x-mean1)**2/(2 * sigma1**2)) +
                                  C2 * exp(-(x-mean2)**2/(2 * sigma2**2))), data=df,
                         start=list(C1=25, mean1=5.88, sigma1=0.01,
                                    C2=30, mean2=5.89, sigma2=0.01), algorithm="port")  
        newx <- seq(5.80,5.95,length.out = 500)
        fit2G_Ka2_prediction <- predict(fit2G_Ka2, newdata=list(x=newx))
        lines(newx, fit2G_Ka2_prediction, col="red")
        CG2_1  <- coef(fit2G_Ka2)[1]
        meanG2_1 <- coef(fit2G_Ka2)[2] # keV 
        sigG2_1  <- coef(fit2G_Ka2)[3]
        err_sigG2_1 <- coef(summary(fit2G_Ka2))[3,2]
        fwhmG2_1 <- sigG2_1*2.35*1000 #eV
        err_fwhmG2_1 <- err_sigG2_1 * 1000* 2.35
        CG2_2  <- coef(fit2G_Ka2)[4]
        meanG2_2 <- coef(fit2G_Ka2)[5] # keV
        sigG2_2  <- coef(fit2G_Ka2)[6]
        err_sigG2_2 <- coef(summary(fit2G_Ka2))[6,2]
        fwhmG2_2 <- sigG2_2*2.35*1000 #eV
        err_fwhmG2_2 <- err_sigG2_2 * 1000 * 2.35
        GS1 <- GS(newx,C=CG2_1,mn=meanG2_1, sg=sigG2_1) # Gaussian1
        GS2 <- GS(newx,C=CG2_2,mn=meanG2_2, sg=sigG2_2) # Gaussian2
        
        xtext <- 5.86
        text(xtext,50, "Double Gaussian", col="red")
        text(xtext,45, paste("Mean(Ka1)=",sprintf("%.3f",meanG2_1),"keV"), col="blue", cex=0.8)
        text(xtext,40, paste("FWHM(Ka1)=",sprintf("%.1f",fwhmG2_1), "eV"), col="blue", cex=0.8)
        text(xtext,35, paste("Mean(Ka2)=",sprintf("%.3f",meanG2_2),"KeV"), col="green", cex=0.8)
        text(xtext,30, paste("FWHM(Ka2)=",sprintf("%.1f",fwhmG2_2), "eV"), col="green", cex=0.8)
        abline(v=5.88765, col="blue", lty=2)   # Ka1
        abline(v=5.89875, col="green", lty=2)  # Ka2
        lines(newx,GS1, col="blue")
        lines(newx,GS2, col="green")
        cat("Err FWHM(1,2)=", err_fwhmG2_1, err_fwhmG2_2)
    }
    # fit Voigt profile
    #----------------
    histo3<-hist(Energies_Kas, xlab='Calibrated energies', breaks=200,freq = FALSE,
                 main=paste("Ka1+Ka2 energies (File",file,", Noise",noise,", Lib file2 Ka2)",sep=""), 
                 xlim=c(5.85,5.95))
    x <- histo3$mids
    y <- histo3$density
    df <- data.frame(x, y)
    fit2V_Ka2 <- nls(y ~ (C1*Voigt(x,mean1,gamma1,sigma1) + C2*Voigt(x,mean2,gamma2,sigma2)), 
                     data=df, start=list(C1=0.5,mean1=5.895, gamma1=0.004, sigma1=0.002, 
                                         C2=3,mean2=5.905, gamma2=0.004, sigma2=0.002), 
                     algorithm = "port") 
    fit2V_Ka2_prediction <- predict(fit2V_Ka2, newdata=list(x=newx))
    lines(newx, fit2V_Ka2_prediction, col="red")
    
    CV2_1  <- coef(fit2V_Ka2)[1]
    meanV2_1 <- coef(fit2V_Ka2)[2] # keV 
    gamV2_1  <- coef(fit2V_Ka2)[3]
    sigV2_1  <- coef(fit2V_Ka2)[4]
    err_sigV2_1 <- coef(summary(fit2V_Ka2))[4,2]
    fwhmV2_1 <- sigV2_1*2.35*1000 #eV
    err_fwhmV2_1 <- err_sigV2_1 * 1000* 2.35
    CV2_2  <- coef(fit2V_Ka2)[5]
    meanV2_2 <- coef(fit2V_Ka2)[6] # keV
    gamV2_2  <- coef(fit2V_Ka2)[7]
    sigV2_2  <- coef(fit2V_Ka2)[8]
    err_sigV2_2 <- coef(summary(fit2V_Ka2))[8,2]
    fwhmV2_2 <- sigV2_2*2.35*1000 #eV
    err_fwhmV2_2 <- err_sigV2_2 * 1000 * 2.35
    V1 <- CV2_1*Voigt(newx,meanV2_1, gamV2_1, sigV2_1) # Voigt1
    V2 <- CV2_2*Voigt(newx,meanV2_2, gamV2_2, sigV2_2) # Voigt2
    
    xtext <- 5.86
    text(xtext,50, "Double Voigt", col="red")
    text(xtext,45, paste("Mean(Ka1)=",sprintf("%.3f",meanV2_1),"keV"), col="blue", cex=0.8)
    text(xtext,40, paste("FWHM(Ka1)=",sprintf("%.1f",fwhmV2_1),"+/-",
                         sprintf("%.1f",err_fwhmV2_1), "eV"), col="blue", cex=0.8)
    text(xtext,35, paste("Gamma(Ka1)=",sprintf("%.1f",gamV2_1*1000), "eV"), col="blue", cex=0.8)
    text(xtext,30, paste("Mean(Ka2)=",sprintf("%.3f",meanV2_2),"KeV"), col="green", cex=0.8)
    text(xtext,25, paste("FWHM(Ka2)=",sprintf("%.1f",fwhmV2_2),"+/-",
                         sprintf("%.1f",err_fwhmV2_2), "eV"), col="green", cex=0.8)
    text(xtext,20, paste("Gamma(Ka2)=",sprintf("%.1f",gamV2_2*1000), "eV"), col="green", cex=0.8)
    abline(v=5.88765, col="blue", lty=2)   # Ka1
    abline(v=5.89875, col="green", lty=2)  # Ka2
    lines(newx,V1, col="blue")
    lines(newx,V2, col="green")
    cat("Err FWHM(1,2)=", err_fwhmV2_1, err_fwhmV2_2)
    
    
    
    #fun.to.minimize <- function(pars,x,y){
        C1     <- pars[1]
        mean1  <- pars[2]
        gamma1 <- pars[3]
        sigma1 <- pars[4]
        C2     <- pars[5]
        mean2  <- pars[6]
        gamma2 <- pars[7]
        #sigma2 <- pars[8]
        result <- sum(y-C1*Voigt(x,mean1,gamma1,sigma1)-C2*Voigt(x,mean2,gamma2,sigma1))**2
        return(result)
    }
    #fit2V_Ka2 <-optim(c(0.45,5.88765,0.004,0.002,0.5,5.89875,0.002), fun.to.minimize, x=x, y=y)
    #C1     <- fit2V_Ka2$par[1]
    #mean1  <- fit2V_Ka2$par[2]
    #gamma1 <- fit2V_Ka2$par[3]
    #sigma1 <- fit2V_Ka2$par[4]
    #C2     <- fit2V_Ka2$par[5]
    #mean2  <- fit2V_Ka2$par[6]
    ##gamma2 <- gamma1
    #sigma2 <- sigma1
    #gamma2 <- fit2V_Ka2$par[7]
    ##sigma2 <- fit2V_Ka2$par[8]
    #newx <- seq(5.80,5.95,length.out = 500)
    #fit2V_Ka2_prediction <- C1*Voigt(newx,mean1,gamma1,sigma1)+C2*Voigt(newx,mean2,gamma2,sigma2)
    #lines(newx, fit2V_Ka2_prediction, col="red")
}


#################################################
plot=TRUE
if (plot){
    len <- c(1024,2048,4096,8192)
    fwhmOPT <- data.frame(L1024=c(14.086, 0.4, 9.597, 0.16),
                          L2048=c(14.086, 0.5, 9.597, 0.16),
                          L4096=c(13.933, 0.4, 9.016, 0.15),
                          L8192=c(12.544, 0.4, 8.810, 0.15),
                          row.names = c("fwhmKa1", "err1", "fwhmKa2", "err2"))
    fwhm0PAD <- data.frame(L1024=c(13.718, 0.4, 8.750, 0.14), 
                           L2048=c(13.777, 0.4, 8.668, 0.16),
                           L4096=c(13.370, 0.5, 8.610, 0.16), 
                           L8192=c(12.544, 0.4, 8.810, 0.15),
                           row.names = c("fwhmKa1", "err1", "fwhmKa2", "err2"))
    fwhmpB75 <- data.frame(L1024=c(NA,NA,NA,NA), 
                           L2048=c(19.772, 0.7, 12.918, 0.3),
                           L4096=c(19.072, 0.8, 11.827, 0.3),
                           L8192=c(18.828, 0.7, 11.537, 0.2),
                       row.names = c("fwhmKa1", "err1", "fwhmKa2", "err2"))
    
    plot(len, len, xlab = "record length (samples)", ylab="FWHM(eV)", type="n",
         ylim=c(8,13))
    errbar(len, fwhmOPT[3,], yplus=as.numeric(fwhmOPT[3,]+fwhmOPT[4,]), 
           yminus=as.numeric(fwhmOPT[3,]-fwhmOPT[4,]), col="blue", 
           errbar.col="blue",pch=1, add=TRUE)
    errbar(len, fwhm0PAD[3,], yplus=as.numeric(fwhm0PAD[3,]+fwhm0PAD[4,]), 
           yminus=as.numeric(fwhm0PAD[3,]-fwhm0PAD[4,]), col="green", 
           errbar.col="green", pch=2, add=TRUE)
    errbar(len, fwhmpB75[3,], yplus=as.numeric(fwhmpB75[3,]+fwhmpB75[4,]), 
           yminus=as.numeric(fwhmpB75[3,]-fwhmpB75[4,]), col="red",
           errbar.col="red", pch=4, add=TRUE)
    legend("topright", legend=c("OPTFILT","0PAD", "pB75"), pch=c(1,2,4), 
           col=c("blue", "green", "red"), bty="n")
}
