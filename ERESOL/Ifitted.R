#
#
# Ifit fitting
#
library(FITSio)
library(Hmisc)
library(RColorBrewer)


polyCurve <- function(x,coeffs) {
    polyres <- 0
    for (icoeff in 1:length(coeffs)){
        index <- icoeff - 1
        polyres <- polyres + coeffs[icoeff] * x^(index)
    }
    #a3[im]*x^3 + a2[im]*x^2 + a1[im]*x + a0[im]
    return(polyres)
}
npolyInit <- 3 # degree of polynomial to be fitted

#gainDir <- "/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresolLPA75um/gainScale"
gainDir <- "/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresol/gainScale"
Ndata<- 5016
EkeVs <- c(0.2,0.5,1,2,3,4,5,6,7,8)
#Ifits<-c(13,14,15,16,17)
Ifits<-c(-19000, -20000, -21000)
#rownames<-paste(Ifits, "muA",sep="")
rownames<-paste(Ifits, "ADU",sep="")
colnames<-paste(EkeVs,"keV",sep="")
Erecons <-matrix(NA,nrow=length(Ifits), ncol=length(EkeVs),dimnames=list(rownames,colnames))
Erecons_SE <-matrix(NA,nrow=length(Ifits), ncol=length(EkeVs),dimnames=list(rownames,colnames))

# read reconstructed data at different Ifits and Energies
for (ii in 1:length(Ifits)){
    Ifit <- Ifits[ii]
    IfitStr <- paste("Ifit_",Ifit,sep="")
    if (Ifit < 0) IfitStr <- paste("Ifit_m",abs(Ifit),sep="")
    cat("IfitStr=", IfitStr)
    for (ie in 1:length(EkeVs)){
        EkeV <- EkeVs[ie]
        #evtfile <- paste(gainDir,"/Ifit",Ifit,"muA/events_sep40000sam_5000p_SIRENA8192_pL8192_",
        #                 EkeV,"keV_STC_F0F_fixedlib6OF_I2RFITTED8192_jitter_bbfb_HR.fits",sep="")
        evtfile <- paste(gainDir,"/events_sep40000sam_5000p_SIRENA8192_pL8192_",
                         EkeV,"keV_STC_T_fixedlib6OF_I2RFITTED8192_",IfitStr,"_fll_HR.fits",sep="")
        #evtfile <- paste(gainDir,"/events_pp",EkeV,".fits",sep="")
        zz <- file(description = evtfile, open = "rb")
        header0 <- readFITSheader(zz, fixHdr = 'remove') # read primary header
        header <- readFITSheader(zz, fixHdr = 'remove') # read extension header
        evtTable <- readFITSbintable(zz, header)
        zz.hdr   <- parseHdr(header)
        close(zz)
        idcol <- which(evtTable$colNames == "SIGNAL")
        npulses <- as.numeric(zz.hdr[which(zz.hdr=="NAXIS2")+1])
        Erecons[ii,ie] <- mean(evtTable$col[[idcol]])
        Erecons_SE [ii,ie] <- sd(evtTable$col[[idcol]])/sqrt(npulses)
    }
}

pdf(paste(gainDir,"/Ifits.pdf",sep=""), width=10, height=7)
layout(matrix(1:length(EkeVs),2,2,byrow = TRUE))

# foreach checking energy...
coeffs <- matrix(0,nrow=ncol(Erecons),ncol=(npolyInit+1))
colors=rainbow(length(EkeVs))

for (i in 1:length(EkeVs)){
    # plot data points
    errbar(Ifits,as.numeric(Erecons[,i]), col=colors[i],pch=16, 
           yplus=as.numeric(Erecons[,i])+as.numeric(Erecons_SE[,i]), 
           yminus=as.numeric(Erecons[,i])-as.numeric(Erecons_SE[,i]),
           xlab="Ifit", ylab="E+/-SE (keV)", 
           errbar.col=colors[i])
           title(main=paste("Reconstructed energy (",colnames(Erecons)[i],")",sep=""))
    #lfit <- lm(Erecons[,i]~Ifits, weights = Erecons_sigma[,i]**2)
    #abline(lfit)

    #npoly <- min(npolyInit, (length(Ifits)-1))
    #badCoeffs <- npoly 
    
    # first fit with orthogonal polynomia (assumes uncorrelated error)
#    while (badCoeffs > 0){
#        stopifnot(npoly > 0)
#        fit <- lm(Erecons[,i]~ poly(Ifits,npoly,raw=FALSE))
#        probNrelev <- rep(-1,npoly)
#        probNrelev <- as.numeric(summary(fit)$coefficients[2:(npoly+1),4])
#        if (TRUE %in% is.na(probNrelev)){
#            badCoeffs <- 1
#        }else{
#            badCoeffs <- length(probNrelev[probNrelev>0.001])  # not relevant coefficients
#        }
#        if (npoly== 1 && as.numeric(summary(fit)$coefficients[1,4])>0.001){ # bad Intercept
#            stop("Error: bad intercept!")
#        }else if(badCoeffs > 0){
#            npoly <- npoly - badCoeffs
#        }
#    }
    # Once polynomia degree is set, fit with raw=TRUE
    npoly=2
    coeffs[i,] <- 0
    fit <- lm(Erecons[,i]~ poly(Ifits,npoly,raw=TRUE))
    coeffs[i,1:(npoly+1)] <- fit$coefficients[1:(npoly+1)]
    curve(polyCurve(x,coeffs[i,]), add=TRUE)
    abline(h=EkeVs[i], lty=2, col=colors[i], lw=2)
}
dev.off()

############## LOOK FOR MINIMIZATION  #####################
par(mfrow=c(1,2))
# Plot surface
nEs<-length(EkeVs)
colors=brewer.pal(n = length(Ifits), name = "Set1")
plot(EkeVs[1:nEs/2], EkeVs[1:nEs/2],type="n",log = "xy", 
     xlab="Input energy (keV)", ylab="Reconstructed PH")
for (ii in 1:length(Ifits)){
 lines(EkeVs[1:nEs/2], Erecons[ii,1:nEs/2], col=colors[ii])   
}
lines(EkeVs[1:nEs/2], EkeVs[1:nEs/2], lty=2, col="black")
legend("topleft",legend=c(Ifits,"linear"), col=c(colors,"black"),
       lty=c(rep(1,length(Ifits)),2), cex=0.8, bty="n")

plot(EkeVs[(nEs/2+1):nEs], EkeVs[(nEs/2+1):nEs],type="n",log = "xy",
     xlab="Input energy (keV)", ylab="Reconstructed PH")
for (ii in 1:length(Ifits)){
    lines(EkeVs[(nEs/2+1):nEs], Erecons[ii,(nEs/2+1):nEs], col=colors[ii])   
}
lines(EkeVs[(nEs/2+1):nEs], EkeVs[(nEs/2+1):nEs], lty=2, col="black")
legend("topleft",legend=c(Ifits,"linear"), col=c(colors,"black"),
       lty=c(rep(1,length(Ifits)),2), cex=0.8, bty="n")

# Plot kind of chi.square
chi2<-numeric(length(Ifits))
for (ii in 1:length(Ifits)){
    chi2[ii]<-sum(((Erecons[ii,]-EkeVs)/Erecons_SE[ii,])**2)/1E8
}
plot(Ifits,chi2)

# Try fitting a parabola
npoly=2
coeffs<-numeric(npoly+1)
fit <- lm(chi2~ poly(Ifits,npoly,raw=TRUE))
fstat2 <- summary(fit)$fstatistic
coeffs <- fit$coefficients[1:(npoly+1)]
curve(polyCurve(x,coeffs), add=TRUE,lty=1)
fpol <- function(x){
    return(polyCurve(x,coeffs))
}
opt2<-optimize(fpol, interval = c(12,18))
mini2<-opt2$minimum
chi22<-as.numeric(opt2$objective)

# Try fitting a 3-order polynomia
npoly=3
coeffs<-numeric(npoly+1)
fit <- lm(chi2~ poly(Ifits,npoly,raw=TRUE))
coeffs <- fit$coefficients[1:(npoly+1)]
curve(polyCurve(x,coeffs), add=TRUE, lty=2)
fpol <- function(x){
    return(polyCurve(x,coeffs))
}
opt3<-optimize(fpol, interval = c(12,18))
mini3<-opt3$minimum
chi23<-as.numeric(opt3$objective)

F23 <-((chi22 - chi23)/1) / (chi23/1)
pF23 <- pf(F23,1,1,lower.tail=FALSE)
cat("pF23=",pF23)
if((1-pF23)>0.9){
    mini=mini3
}else{
    mini=mini2
}
miniStr2<-sprintf("%.2f",mini2)
miniStr3<-sprintf("%.2f",mini3)
cat("Best fit Ifit(pol2)=",miniStr2,"muA\n")
cat("Best fit Ifit(pol3)=",miniStr3,"muA\n")
text(x=min(Ifits)+(max(Ifits)-min(Ifits))/2,y=0.9*max(chi2),
     paste("BF Ifit(p2)=",miniStr2,"muA\n"),cex=0.8)
text(x=min(Ifits)+(max(Ifits)-min(Ifits))/2,y=0.8*max(chi2),
     paste("BF Ifit(p3)=",miniStr3,"muA\n"),cex=0.8)
legend("bottomright",c("pol2","pol3"), lty=c(1,2),bty="n", cex=0.5)
