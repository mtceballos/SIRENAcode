#
#  Script to evaluate how often pairs (or multipulses) ocurr in a XIFU array at different separations and ctrates
#  ==============================================================================================================
#  
#  Using simulations and poisson statistics (based on testpois.R script)
#
setwd("/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/baselineLPA2")
pixel<- "LPA2 35mmBe"
pixelShort <- gsub("^([A-Z]*[1-9]*) .*", "\\1", pixel)
instrument <- "XIFU4Q"
pdf(paste("MultiPulses_",instrument,"_",gsub(' ','_',pixel),".pdf",sep=""))
separations <- round(10**(seq(0, log10(400), length.out = 5))) #samples
nseps <- length(separations)
TimeInterval <- 1000. # interval time (seconds)
samprate <- 156250.
mCrab <- 94. # ct/s
ctrmax <-1000.*mCrab
# initialize count rates values
ctrates <- 10**(seq(0, log10(ctrmax),length.out=10))
ncrts <- length(ctrates)
colors <- rainbow(nseps)

timeConst <- list("LPA2" = list("t0"=1234./samprate, "t1"=146./samprate, "t2"=354./samprate,"t3"=7080./samprate,"tolsep"=4./samprate))
t0 <- timeConst[[pixelShort]][["t0"]]
t1 <- timeConst[[pixelShort]][["t1"]]
t2 <- timeConst[[pixelShort]][["t2"]]
t3 <- timeConst[[pixelShort]][["t3"]]
tolsep <- timeConst[[pixelShort]][["tolsep"]]
psffracs <- list("LPA1"=sort(c(1e-3, 3e-3, 6e-3, 3e-3, 3e-3,0.032, 0.091, 0.033, 3e-3, 6e-3, 0.091, 0.441, 
                               0.091, 6e-3, 3e-3, 0.032, 0.090, 0.032, 3e-3, 3e-3, 6e-3, 3e-3), decreasing = TRUE),
                 "LPA2 in focus"=sort(c(0.002, 0.005, 0.003, 0.002, 0.030, 0.089, 0.031, 0.002, 0.005, 0.089, 0.463, 
                                        0.090, 0.005, 0.002, 0.030, 0.089, 0.030, 0.003, 0.003, 0.005, 0.003 ), decreasing = TRUE))
LPA2Ppsf <- read.csv(file="/home/ceballos/INSTRUMEN/ATHENA/XIFU/psf_35mm_Be_ratios.txt", header=TRUE,sep = ",")
psffracs[["LPA2 35mmBe"]] <- LPA2Ppsf$Ratio
sumPSF <- sum(psffracs[[pixel]])
npsfs <- length(psffracs[[pixel]])

Prob.photons <- matrix(nrow=nseps, ncol=ncrts)   # piled-up photons (whichever ph with a <sep prev/next)
Prob1.photons <- matrix(nrow=nseps, ncol=ncrts)  # piled-up photons (PRIMARY: ph with a <sep next but prev>sep)
Prob2.photons <- matrix(nrow=nseps, ncol=ncrts)  # piled-up photons (SECONDARY: ph with a <sep prev)
Prob.Poisson <- matrix(nrow=nseps, ncol=ncrts)
Prob1.Poisson <- matrix(nrow=nseps, ncol=ncrts)
Prob2.Poisson <- matrix(nrow=nseps, ncol=ncrts)

legend <- character(nseps)
drawLogPlotBox(xlimits=c(0.01,1000.),ylimits=c(0,1), x2limits=c(0.01,1000.),y2limits=c(0,1),
logxy="x", naxes=c(T,T,F,T), xlabel="Intensity (mCrab)", ylabel="Fraction of photons")
title(main=paste("Pairs/multipulses Ocurrence \n (",pixel,"; pixel w/ max ctrate)",sep=""))

listPixelFracs <- list()
for (is in 1:nseps){ # time interval for "pile-up": once we have a photon, prob of having another one in the given separation
    cat("For separation=", separations[is], " and...\n")
    tolsep <- separations[is]/ samprate  # s
    tolsepStr <- sprintf("%5.3f",tolsep*1000.) #ms
    
    for (ic in 1:ncrts){ # different count rates
        cat("Ctrate=", ctrates[ic]/mCrab, " mCrab\n")
        lambda <- ctrates[ic] * max(psffracs[[pixel]])  # in the central pixe (max ctrate)
        listPixelFracs <- XIFUpixelFracsHMTL(lambda, TimeInterval, t0,t1,t2,t3,tolsep)
        Prob.photons[is,ic] <- listPixelFracs[["PL_ph"]]
        Prob1.photons[is,ic] <- listPixelFracs[["P1L_ph"]]
        Prob2.photons[is,ic] <- listPixelFracs[["P2L_ph"]]
        Prob.Poisson[is,ic] <- listPixelFracs[["PL_Ps"]]
        Prob1.Poisson[is,ic] <- listPixelFracs[["P1L_Ps"]]
        Prob2.Poisson[is,ic] <- listPixelFracs[["P2L_Ps"]]
    }
    lines(ctrates/mCrab, Prob.Poisson[is,], col=colors[is] )
    lines(ctrates/mCrab, Prob1.Poisson[is,], col=colors[is], lty=2 )
    lines(ctrates/mCrab, Prob2.Poisson[is,], col=colors[is], lty=2 )
    points(ctrates/mCrab, Prob.photons[is,], col=colors[is], pch=1 )
    points(ctrates/mCrab, Prob1.photons[is,], col=colors[is], pch=8 )
    points(ctrates/mCrab, Prob2.photons[is,], col=colors[is], pch=20 )
    legend[is] <- paste(separations[is], "samps =",tolsepStr, "ms")
    #text(ctrates[30],1.5*Prob.Poisson[is,30],srt=30, lab=paste(separations[is], "samples"))
    
}
legend("topleft", legend=c(legend, "Multiple Pulse", "Primary in multiple", "Secondary in multiple"), cex=0.8, bty="n",
                           col=c(colors,"black", "black", "black"), pch=c(rep(1,nseps),1,8,20),lty=c(rep(1,nseps),1,2,2))
abline(h=seq(0,1,by=0.1),col = "gray", lty = "dotted",)

dev.off()
