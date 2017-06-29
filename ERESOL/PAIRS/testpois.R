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
    
    #
    # Poissonian Distribution of Photons
    #   tp = Tprevious
    #   tn = Tnext
    #   L : event with "Low" Quality
    #   M : event with "Medium" Quality
    #   H : event with "High" Quality
    #   I : event "Invalid"
    #===================================================
    #             |  tp<t0   tp>t0        =
    # -------------------------------------------------=
    #   tn<t1    |      I        L        =
    # t1<tn<t2   |      I        M        =
    #   tn>t2    |      I        H        =
    #===================================================
    # Initial (example) valuesto generate photons(events)
    tmin <- 0.0         # initial time (seconds)
    tmax <- 1000         # final time (seconds)
    
    #
    # B) CHECK that the distribution of the quality of the photons
    # follow the theoretical law based on Poissonian statistics
    #===============================================================
    # Example values to compare with Boyce et al 1999 (Maite's stuff)
    # XIFU-LPA1
    t0 <- 2.56E-3 #s
    t1 <- 1.64e-3 #s 
    t2 <- 6.55e-3 #s 
    mCrab <- 94. # ct/s
    psffracs <- sort(c(1e-3, 3e-3, 6e-3, 3e-3, 3e-3,0.032, 0.091, 0.033, 3e-3, 6e-3, 
                       0.091, 0.441, 0.091, 6e-3, 3e-3, 0.032, 0.090, 0.032, 3e-3, 
                       3e-3, 6e-3, 3e-3), decreasing = TRUE)
    sumPSF <- sum(psffracs)
    #psffracs <- c(0.441, 0.091, 0.091)
    npsfs <- length(psffracs)
    # FRACTION OF PHOTONS for each quality at different count rates
    
    # initialize count rates values
    ctrates <- 10**(seq(0, log10(400*95),length.out=25))
    
    # initalize vector to store fractions of events
    FracHQ <- c()
    FracMQ <- c()
    FracLQ <- c()
    FracIQ <- c()
    pL <- c()
    pH <- c()
    pM <- c()
    pI <- c()
    # run over list of count rates
    for (ctr in ctrates){
        
        Frac.HQ.photons <- 0.
        Frac.MQ.photons <- 0.
        Frac.LQ.photons <- 0.
        Frac.IQ.photons <- 0.
        
        p1 <- 0.
        p2 <- 0.
        p3 <- 0.
        
        for (psf in psffracs){ # foreach pixel with incident photons
            lambda <- ctr * psf # PSF (if source is centred on a pixel, it receives psf*100% of incident flux)
            N.mean.photons <- (tmax-tmin)*lambda
            stopifnot(N.mean.photons>0)
            N.real.photons <- rpois(1,N.mean.photons)
            #cat("psfrac=",psf, "N.mean.photons=", N.mean.photons,"\n")
            points <- runif(N.real.photons,tmin,tmax)
            points <- sort(points)
            N.HQ.photons <- 0.
            N.MQ.photons <- 0.
            N.LQ.photons <- 0.
            N.IQ.photons <- 0.
            
            Frac.HQ.photons.psf <- 0.
            Frac.MQ.photons.psf <- 0.
            Frac.LQ.photons.psf <- 0.
            Frac.IQ.photons.psf <- 0.
            
            #cat("N.real.photons=",N.real.photons)
            if(N.real.photons < 2){  # if less than 2 photons in total
                # Frac.HQ.photons <- Frac.HQ.photons + 1.*psf/sumPSF
                # Frac.MQ.photons <- Frac.MQ.photons + 0.
                # Frac.LQ.photons <- Frac.LQ.photons + 0.
                # Frac.IQ.photons <- Frac.IQ.photons + 0.
                Frac.HQ.photons.psf <- 1.
                Frac.MQ.photons.psf <- 0.
                Frac.LQ.photons.psf <- 0.
                Frac.IQ.photons.psf <- 0.
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
                        if((points[i]-points[i-1])>t0){   # HQ
                            N.HQ.photons <- N.HQ.photons +1
                        }else{                           # LQ
                            N.I.photons <- N.IQ.photons +1
                        }
                    }else{
                        if((points[i]-points[i-1])>t0 && (points[i+1]-points[i])>t2){  # HQ
                            N.HQ.photons <- N.HQ.photons +1 
                        }else if((points[i]-points[i-1])>t0 && (points[i+1]-points[i])>t1 && 
                                 (points[i+1]-points[i])<t2){  # MQ
                            N.MQ.photons <- N.MQ.photons +1 
                        }else if((points[i]-points[i-1])>t0 && (points[i+1]-points[i])<t1){ # LQ
                            N.LQ.photons <- N.LQ.photons +1
                        }else{                          # IQ
                            N.IQ.photons <- N.IQ.photons +1
                        }
                    } # if for first/last/interm. photons
                }# foreach photon
                # Frac.HQ.photons <- Frac.HQ.photons + (N.HQ.photons/N.real.photons)*psf/sumPSF
                # Frac.MQ.photons <- Frac.MQ.photons + (N.MQ.photons/N.real.photons)*psf/sumPSF
                # Frac.LQ.photons <- Frac.LQ.photons + (N.LQ.photons/N.real.photons)*psf/sumPSF
                # Frac.IQ.photons <- Frac.IQ.photons + (N.IQ.photons/N.real.photons)*psf/sumPSF
                Frac.HQ.photons.psf <- N.HQ.photons/N.real.photons
                Frac.MQ.photons.psf <- N.MQ.photons/N.real.photons
                Frac.LQ.photons.psf <- N.LQ.photons/N.real.photons
                Frac.IQ.photons.psf <- N.IQ.photons/N.real.photons
            }# if >2 photons
            #cat("Frac.HQ=",Frac.HQ.photons," Frac.MQ=",Frac.MQ.photons," Frac.LQ=",Frac.LQ.photons,"\n")
            Frac.HQ.photons <- Frac.HQ.photons + Frac.HQ.photons.psf*psf/sumPSF
            Frac.MQ.photons <- Frac.MQ.photons + Frac.MQ.photons.psf*psf/sumPSF
            Frac.LQ.photons <- Frac.LQ.photons + Frac.LQ.photons.psf*psf/sumPSF
            Frac.IQ.photons <- Frac.IQ.photons + Frac.IQ.photons.psf*psf/sumPSF
            
            # Expression from Poissonian theory
            # HQ
            p1f <- ppois(0,lambda*t0)*ppois(0,lambda*t2)
            p1 <- p1 + p1f*psf/sumPSF
            # MQ
            p2f <- ppois(0,lambda*t0)*ppois(0,lambda*t1)-p1f
            p2 <- p2 + p2f*psf/sumPSF
            # LQ
            p3f <- ppois(0,lambda*t0)-p1f-p2f
            p3 <- p3 + p3f*psf/sumPSF
            #cat("    p1=", p1f, " p2=",p2f, "p3=",p3f,"\n")
            
        } # foreach psffrac
        
        FracHQ <- append(FracHQ,Frac.HQ.photons)
        FracMQ <- append(FracMQ,Frac.MQ.photons)
        FracLQ <- append(FracLQ,Frac.LQ.photons)
        FracIQ <- append(FracIQ,Frac.IQ.photons)
        
        pH <- append(pH,p1)
        pM <- append(pM,p2)
        pL <- append(pL,p3)
        pI <- 1 - pH -pM - pL
    
    } #foreach ctrate
    
    drawLogPlotBox(xlimits=c(0.01,1000.),ylimits=c(0,1), x2limits=c(0.01,1000.),y2limits=c(0,1),
                   logxy="x", naxes=c(T,T,F,F), xlabel="Intensity (mCrab)", ylabel="Fraction of photons")
    title(main="Quality Grading of photons")
    points(ctrates/mCrab,FracHQ,col="red")
    points(ctrates/mCrab,FracMQ,col="blue")
    points(ctrates/mCrab,FracLQ,col="green")
    points(ctrates/mCrab,FracIQ,col="gray")
    abline(h=0.8, lty=2, col="gray")
    abline(v=1, lty=2, col="gray") # 1mCrab
    legend("top",c("H.Quality","M.Quality","L.Quality","I.Quality","(Lines from theory)"),cex=0.8,
           pch=c(21,21,21,21,NA_integer_),
           col=c("red","blue","green","gray","black"),
           lty=c(1,1,1,1,0))
    lines(ctrates/mCrab,pH,col="red")
    lines(ctrates/mCrab,pM,col="blue")
    lines(ctrates/mCrab,pL,col="green")
    lines(ctrates/mCrab,pI,col="gray")
}
