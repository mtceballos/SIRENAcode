# By NCL, MTC 2014
#

instrument<-"XIFU" #or ASTROH

if(instrument == "ASTROH"){
    #
    #   For ASTRO-H
    #
    # Poissonian Distribution of Photons
    #   tp = Tprevious
    #   tn = Tnext
    #   L : event with "Low" Quality
    #   M : event with "Medium" Quality
    #   H : event with "High" Quality
    #===================================================
    #             |  tp<t1     t1<tp<t2   tp>t2        =
    # -------------------------------------------------=
    #   tn<t1    |     L          L          L        =
    # t1<tn<t2   |     L          M          M        =
    #   tn>t2    |     L          M          H        =
    #===================================================
    # set plot arrangement
    mat = matrix(data=1:2, nrow=1, ncol=2,byrow=T)
    layout(mat)
    # Initial (example) valuesto generate photons(events)
    tmin <- 0.0         # initial time (seconds)
    tmax <- 1000         # final time (seconds)
    
    # A) CHECK that histogram of differences between consecutive 
    # photons follow an exponential law
    # ============================================================
    lambda <- 10        # mean number of photons/second
    # expected mean number of photons
    N.mean.photons <- (tmax-tmin)*lambda
    # actual number of photons received
    N.real.photons <- rpois(1,N.mean.photons)
    #
    cat(">>> N.mean.photons: ",N.mean.photons,"\n")
    cat(">>> N.real.photons: ",N.real.photons,"\n")
    # generate random arrival of photons as a function of time
    points <- runif(N.real.photons,tmin,tmax)
    # sort arrival times
    points <- sort(points)
    # differences (in seconds) between consecutive photons
    differences <- points[-1] - points[-N.real.photons]
    #
    cat(">>> Expected mean differences (seconds): ",1/lambda,"\n")
    cat(">>> Measured mean differences (seconds): ",mean(differences),"\n")
    # histogram of time between events (consecutive arrival of photons) 
    hist(differences, freq=FALSE, main="Differences (lambda=10 ct/s)",
         xlab="time difference between consecutive photons (s)",
         ylab="Probability density")
    rug(differences, col="blue")
    # overplot expected negative exponential
    curve(dexp(x, rate=lambda), from=0, to=10*mean(differences), 
          col="red", add=TRUE)
    
    #
    # B) CHECK that the distribution of the quality of the photons
    # follow the theoretical law based on Poissonian statistics
    #===============================================================
    # Example values to compare with Boyce et al 1999 (Maite's stuff)
    t1 <- 35.5E-3 # s
    t2 <- 142.0E-3 # s
    
    # FRACTION OF PHOTONS for each quality at different count rates
    
    # initialize count rates values
    ctrates <- seq(0,25,1)
    
    # initalize vector to store fractions of events
    FracHQ <- c()
    FracMQ <- c()
    FracLQ <- c()
    pL <- c()
    pH <- c()
    pM <- c()
    # run over list of count rates
    for (lambda in ctrates){
        N.mean.photons <- (tmax-tmin)*lambda
        N.real.photons <- rpois(1,N.mean.photons)
        points <- runif(N.real.photons,tmin,tmax)
        points <- sort(points)
        N.HQ.photons <- 0.
        N.MQ.photons <- 0.
        N.LQ.photons <- 0.
        if(N.real.photons < 2){  # if less than 2 photons in total
            Frac.HQ.photons <- 1.
            Frac.MQ.photons <- 0.
            Frac.LQ.photons <- 0.
        }else{
            for (i in 1:N.real.photons){
                if(i==1){            # first photon case
                    if((points[i+1]-points[i])>t2){  # HQ
                        N.HQ.photons <- N.HQ.photons +1
                    }else if((points[i+1]-points[i])>t1 && (points[i+1]-points[i])<t2){  # MQ
                        N.MQ.photons <- N.MQ.photons +1  
                    }else{                           # LQ
                        N.LQ.photons <- N.LQ.photons +1
                    }
                }else if(i==N.real.photons){ # last photon case
                    if((points[i]-points[i-1])>t2){   # HQ
                        N.HQ.photons <- N.HQ.photons +1
                    }else if((points[i]-points[i-1])>t1 && (points[i]-points[i-1])<t2){  # MQ
                        N.MQ.photons <- N.MQ.photons +1  
                    }else{                           # LQ
                        N.LQ.photons <- N.LQ.photons +1
                    }
                }else{
                    if((points[i]-points[i-1])>t2 && (points[i+1]-points[i])>t2){  # HQ
                        N.HQ.photons <- N.HQ.photons +1
                    }else if( ((points[i]-points[i-1])>t1 && 
                               (points[i]-points[i-1])<t2 && 
                               (points[i+1]-points[i])>t1)  
                              ||((points[i+1]-points[i])>t1 && 
                                 (points[i+1]-points[i])<t2 &&
                                 (points[i]-points[i-1])>t2) ){   # MQ
                        N.MQ.photons <- N.MQ.photons +1
                    }else{                          # LQ
                        N.LQ.photons <- N.LQ.photons +1
                    }
                } # if for i values
            }# foreach ph
            Frac.HQ.photons <- N.HQ.photons/N.real.photons
            Frac.MQ.photons <- N.MQ.photons/N.real.photons
            Frac.LQ.photons <- N.LQ.photons/N.real.photons
        }# nof ph
        FracHQ <- append(FracHQ,Frac.HQ.photons)
        FracMQ <- append(FracMQ,Frac.MQ.photons)
        FracLQ <- append(FracLQ,Frac.LQ.photons)
        
        # Expression from Poissonian theory
        # HQ
        p1 <- ppois(0,lambda*t2)*ppois(0,lambda*t2)
        pH <- append(pH,p1)
        # MQ
        p2 <- ppois(0,lambda*t1)*ppois(0,lambda*t1)-p1
        pM <- append(pM,p2)
        # LQ
        pL <- 1 - pH -pM
    }
    plot(ctrates,FracHQ,ylim=c(0,1),col="red",log="x",
         xlab="Count rate (ct/s)", ylab="Fraction of photons",
         main="Quality Grading of photons")
    points(ctrates,FracMQ,col="blue")
    points(ctrates,FracLQ,col="green")
    abline(h=0.8)
    legend("top",c("H.Quality","M.Quality","L.Quality","(Lines from theory)"),cex=0.8,
           pch=c(21,21,21,NA_integer_),
           col=c("red","blue","green","black"),
           lty=c(1,1,1,0))
    lines(ctrates,pH,col="red")
    lines(ctrates,pM,col="blue")
    lines(ctrates,pL,col="green")
    
} else if(instrument == "XIFU"){

    #######################################################################################
    #
    #     FOR X-IFU
    #
    #######################################################################################
    
    pixelFracs <- function(ctrate, deltat, t0, t1, t2, tolsep){
        #
        # Calculate (Quality) Fractions of photons for XIFU pixels, using photon simulation AND
        # Poisson statistics
        #
        # Poissonian Distribution of Photons
        #   tp = Tprevious
        #   tn = Tnext
        #   L : event with "Low" Quality
        #   M : event with "Medium" Quality
        #   H : event with "High" Quality
        #   I : event "Invalid"
        #   P : event piled-up (same arrival time)
        #===================================================
        #            | tolsep>tp  tolsep<tp<t0   tp>t0     =
        # -------------------------------------------------=
        #   tn<t1    |      P            I        L        =
        # t1<tn<t2   |      P            I        M        =
        #   tn>t2    |      P            I        H        =
        #  tn<tolsep |      P            P        P        =
        #===================================================
        #
        # ctrate (in) : count rate (phs/s)
        # deltat (in) : time interval for simulations
        # t0, t1, t2, tolsep as in table above
        # fractions (out): list with H(igh)Q, M(edium)Q, L(ow)Q, I(nvalid)Q, P(i)L(eup) quality
        
        N.mean.photons <- deltat*ctrate
        stopifnot(N.mean.photons>0)
        N.real.photons <- rpois(1,N.mean.photons)
        #cat("psfrac=",psf, "N.mean.photons=", N.mean.photons,"\n")
        points <- runif(N.real.photons,0.,deltat)
        points <- sort(points)
        
        # Calculate fractions by photon counting (points in plots)
        # ----------------------------------------------------------
        Frac.HQ <- 0.
        Frac.MQ <- 0.
        Frac.LQ <- 0.
        Frac.IQ <- 0.
        Frac.PL <- 0.
        
        N.HQ.photons <- 0.
        N.MQ.photons <- 0.
        N.LQ.photons <- 0.
        N.IQ.photons <- 0.
        N.PL.photons <- 0.
        
        if(N.real.photons < 2){  # if less than 2 photons in total
            Frac.HQ <- 1.
            Frac.MQ <- 0.
            Frac.LQ <- 0.
            Frac.IQ <- 0.
            Frac.PL <- 0.
        }else{
            for (i in 1:N.real.photons){
                if(i==1){            # first photon case
                    if((points[i+1]-points[i])>t2){  # HQ
                        N.HQ.photons <- N.HQ.photons +1
                    }else if((points[i+1]-points[i])>t1 && (points[i+1]-points[i])<t2){  # MQ
                        N.MQ.photons <- N.MQ.photons +1  
                    }else if((points[i+1]-points[i])<tolsep){
                        N.PL.photons <- N.PL.photons +1  
                    }else{                           # LQ
                        N.LQ.photons <- N.LQ.photons +1
                    }
                }else if(i==N.real.photons){ # last photon case
                    if((points[i]-points[i-1])>t0){   # HQ
                        N.HQ.photons <- N.HQ.photons +1
                    }else if((points[i]-points[i-1])<tolsep){
                        N.PL.photons <- N.PL.photons +1  
                    }else{                           # LQ
                        N.I.photons <- N.IQ.photons +1
                    }
                }else{
                    if((points[i]-points[i-1])>t0 && (points[i+1]-points[i])>t2){  # HQ
                        N.HQ.photons <- N.HQ.photons +1 
                    }else if((points[i]-points[i-1])>t0 && (points[i+1]-points[i])>t1 && 
                             (points[i+1]-points[i])<t2){  # MQ
                        N.MQ.photons <- N.MQ.photons +1 
                    }else if((points[i]-points[i-1])>t0 && ((points[i+1]-points[i])<t1) &&
                             (points[i+1]-points[i])>tolsep){ # LQ
                        N.LQ.photons <- N.LQ.photons +1
                    }else if((points[i+1]-points[i])<tolsep || (points[i]-points[i-1])<tolsep){
                        N.PL.photons <- N.PL.photons +1  
                    }else{                          # IQ
                        N.IQ.photons <- N.IQ.photons +1
                    }
                } # if for first/last/interm. photons
            }# foreach photon
            Frac.HQ <- N.HQ.photons/N.real.photons
            Frac.MQ <- N.MQ.photons/N.real.photons
            Frac.LQ <- N.LQ.photons/N.real.photons
            Frac.IQ <- N.IQ.photons/N.real.photons
            Frac.PL <- N.PL.photons/N.real.photons
        }# if >2 photons
        
        # Calculate fractions from Expression by Poissonian theory: given a photon...
        #----------------------------------------------------------------------------
        # HQ => pH: exp(-ctrate*t0) * exp(-ctrate*t2)
        # Prob of having 0 prev. photons in t0  and     0 next photons in t2
        pH <- dpois(0,lambda*t0)*dpois(0,lambda*t2)
        
        # MQ => pM: exp(-ctrate*t1) *  exp(-ctrate*t1)- pH
        # Prob of having 0 prev. photons in t1  and     0 next photons in t1    - pH
        pM <- dpois(0,lambda*t0)*dpois(0,lambda*t1)-pH
        
        # LQ => pL: exp(-ctrate*t0) -pH - pM
        pL <- dpois(0,lambda*t0)-pH-pM
        
        # PL => pP: two pulses arriving at t<ti
        # 2*exp(-ctrate*ti)  * (1- exp(-ctrate*ti))  : one pulse < ti, the other >ti
        # + (1- exp(-ctrate*ti))*(1- exp(-ctrate*ti)) : tp and tn closer than ti
        pP <- 1-dpois(0,2*lambda*tolsep)  
        # also = ppois(0,2*lambda*tolsep, lower.tail=FALSE) -> given 1 photon, which is the prob of finding >0 ph in the double interval?
        
        # IQ => pI
        pI <- 1-pH-pM-pL-pP
        
        return(list("HQ_ph"=Frac.HQ,"MQ_ph"=Frac.MQ,"LQ_ph"=Frac.LQ,"IQ_ph"=Frac.IQ,"PL_ph"=Frac.PL,
                    "HQ_Ps"=pH,     "MQ_Ps"=pM,     "LQ_Ps"=pL,     "IQ_Ps"=pI,     "PL_Ps"=pP))
    }
    
    deltat <- 1000.         # interval time (seconds)
    samprate <- 156250.
    pixel <- "LPA1"
    mCrab <- 94. # ct/s
    timeConst <- list("LPA1" = list("t0"=2.56E-3,"t1"=1.64e-3,"t2"=6.55e-3,"tolsep"=1./samprate),
                      "LPA2" = list("t0"=2.56E-3,"t1"=1.64e-3,"t2"=6.55e-3,"tolsep"=1./samprate))
    t0 <- timeConst[[pixel]][["t0"]]
    t1 <- timeConst[[pixel]][["t1"]]
    t2 <- timeConst[[pixel]][["t2"]]
    tolsep <- timeConst[[pixel]][["tolsep"]]
    psffracs <- list("LPA1"=sort(c(1e-3, 3e-3, 6e-3, 3e-3, 3e-3,0.032, 0.091, 0.033, 3e-3, 6e-3, 0.091, 0.441, 
                                   0.091, 6e-3, 3e-3, 0.032, 0.090, 0.032, 3e-3, 3e-3, 6e-3, 3e-3), decreasing = TRUE),
                     "LPA2"=sort(c(0.002, 0.005, 0.003, 0.002, 0.030, 0.089, 0.031, 0.002, 0.005, 0.089, 0.463, 
                                   0.090, 0.005, 0.002, 0.030, 0.089, 0.030, 0.003, 0.003, 0.005, 0.003 ), decreasing = TRUE))
    sumPSF <- sum(psffracs[[pixel]])
    npsfs <- length(psffracs[[pixel]])
    
    # initialize count rates values
    ctrates <- 10**(seq(0, log10(1000*mCrab),length.out=25))
    
    # initalize vector to store fractions of events
    FracHQ <- c()
    FracMQ <- c()
    FracLQ <- c()
    FracIQ <- c()
    FracPL <- c()
    pL <- c()
    pH <- c()
    pM <- c()
    pI <- c()
    pP <- c()
    # run over list of count rates
    for (ctr in ctrates){
        cat("Calculating for ctr=",ctr/mCrab," mCrab\n")
        Frac.HQ.photons <- 0.
        Frac.MQ.photons <- 0.
        Frac.LQ.photons <- 0.
        Frac.IQ.photons <- 0.
        Frac.PL.photons <- 0.
        
        pH.pois <- 0.
        pM.pois <- 0.
        pL.pois <- 0.
        pI.pois <- 0.
        pP.pois <- 0.
        
        for (psf in psffracs[[pixel]]){ # foreach pixel with incident photons
            # Correct count rate for PSF distribution
            lambda <- ctr * psf # PSF (if source is centred on a pixel, it receives psf*100% of incident flux)
            
            # calculate fractions for each pixel in the array (different psf)
            listPixelFracs <- pixelFracs(lambda,deltat, t0,t1,t2,tolsep)
            Frac.HQ.photons.psf <- listPixelFracs[["HQ_ph"]]
            Frac.MQ.photons.psf <- listPixelFracs[["MQ_ph"]]
            Frac.LQ.photons.psf <- listPixelFracs[["LQ_ph"]]
            Frac.IQ.photons.psf <- listPixelFracs[["IQ_ph"]]
            Frac.PL.photons.psf <- listPixelFracs[["PL_ph"]]
            
            pH.psf <- listPixelFracs[["HQ_Ps"]]
            pM.psf <- listPixelFracs[["MQ_Ps"]]
            pL.psf <- listPixelFracs[["LQ_Ps"]]
            pI.psf <- listPixelFracs[["IQ_Ps"]]
            pP.psf <- listPixelFracs[["PL_Ps"]]
            
            # add fractions for given pixel taking into account contribution to total gathering of photons
            Frac.HQ.photons <- Frac.HQ.photons + Frac.HQ.photons.psf*psf/sumPSF
            Frac.MQ.photons <- Frac.MQ.photons + Frac.MQ.photons.psf*psf/sumPSF
            Frac.LQ.photons <- Frac.LQ.photons + Frac.LQ.photons.psf*psf/sumPSF
            Frac.IQ.photons <- Frac.IQ.photons + Frac.IQ.photons.psf*psf/sumPSF
            Frac.PL.photons <- Frac.PL.photons + Frac.PL.photons.psf*psf/sumPSF
            
            pH.pois <- pH.pois + pH.psf*psf/sumPSF
            pM.pois <- pM.pois + pM.psf*psf/sumPSF
            pL.pois <- pL.pois + pL.psf*psf/sumPSF
            pI.pois <- pI.pois + pI.psf*psf/sumPSF
            pP.pois <- pP.pois + pP.psf*psf/sumPSF
            
        } # foreach psffrac
        
        FracHQ <- append(FracHQ,Frac.HQ.photons)
        FracMQ <- append(FracMQ,Frac.MQ.photons)
        FracLQ <- append(FracLQ,Frac.LQ.photons)
        FracIQ <- append(FracIQ,Frac.IQ.photons)
        FracPL <- append(FracPL,Frac.PL.photons)
        
        pH <- append(pH,pH.pois)
        pM <- append(pM,pM.pois)
        pL <- append(pL,pL.pois)
        pI <- append(pI,pI.pois)
        pP <- append(pP,pP.pois)
    
    } #foreach ctrate
    
    #
    # PLOT photons simulations and Poisson statistics
    #
    drawLogPlotBox(xlimits=c(0.01,1000.),ylimits=c(0,1), x2limits=c(0.01,1000.),y2limits=c(0,1),
                   logxy="x", naxes=c(T,T,F,T), xlabel="Intensity (mCrab)", ylabel="Fraction of photons")
    title(main="Quality Grading of photons")
    points(ctrates/mCrab,FracHQ,col="red")
    points(ctrates/mCrab,FracMQ,col="blue")
    points(ctrates/mCrab,FracLQ,col="green")
    points(ctrates/mCrab,FracIQ,col="gray")
    points(ctrates/mCrab,FracPL,col="violet")
    abline(h=0.8, lty=2, col="gray")
    abline(v=1, lty=2, col="gray") # 1mCrab
    legend("left",c("High Quality","Medium Quality","Low Quality","Invalid","Pile-Up"),cex=1., bty="n",
           pch=c(21,21,21,21,21,NA_integer_),
           col=c("red","blue","green","gray","violet"),
           lty=c(1,1,1,1,1,0))
    lines(ctrates/mCrab,pH,col="red")
    lines(ctrates/mCrab,pM,col="blue")
    lines(ctrates/mCrab,pL,col="green")
    lines(ctrates/mCrab,pI,col="gray")
    lines(ctrates/mCrab,pP,col="violet")
}
