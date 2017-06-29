#
# Given an average number of events (lambda) in an interval of time (T), 
# the probability of finding k events in that interval is:
#
#                  lambda**k
#  P(k) = exp(-k) ------------            ctrate=lambda/T
#                     k!
#
samprate <- 156250.
separations <- round(10**(seq(0, log10(500), length.out = 10)))
nseps <- length(separations)
colors <- rainbow(nseps)
ctrates <- seq(1,1000,length.out = 50)
ncrts <- length(ctrates)
P_k <- matrix(nrow=nseps, ncol=ncrts)
legend <- character(nseps)
drawLogPlotBox(xlimits=c(min(ctrates),max(ctrates)),ylimits=c(1E-8,1),
               x2limits=c(min(ctrates),max(ctrates)),y2limits=c(1E-8,1),
               logxy="xy", xlabel="Count rate ct/s", 
               ylabel="Prob(>=2 events)",
               naxes=c(T,T,F,T))
title(main="Probability of having >=2 photons in interval given by separation")
for (is in 1:nseps){
    Time <- separations[is]/ samprate
    Timesrt <- sprintf("%5.3f",Time*1000.) #ms
    for (ic in 1:ncrts){
        lambda <- ctrates[ic] * Time
        P_k[is,ic] <- ppois(2,lambda,lower.tail=FALSE) # prob of having >=2 events in a given interval
    }
    lines(ctrates, P_k[is,], col=colors[is], pch=is, type="b" )
    legend[is] <- paste(separations[is], "samps =",Timesrt, "ms")
    text(ctrates[30],1.5*P_k[is,30],srt=30, lab=paste(separations[is], "samples"))
    text(100, 1E-8, lab="mCrab")
    text(1000, 1E-8, lab="Crab")
}
legend("topleft", legend=legend, col=colors, pch=(1:nseps),lty=rep(1,nseps))
grid()

