 #
# SVD decomposition for pulses
#

# .txt file extracted from fv of fits file --> colum SIGNAL --> expand --> 
#  Export as text --> rows:1-1000; columns=1000-5096; User defined separator (blank)
#
#

array <- "LPA2shunt"
separation <- "40000"
nSimPulses <- "20000"
pulseLength <- 4096
size <- 4096
npulsesPerEnergy <- 1000 # 1000
ncomps <- 4 #9
EkeVcalib <- c(0.2,0.5,1,2,3,4,5,6,7,8)
ncols <- npulsesPerEnergy*length(EkeVcalib)
dir <- paste("/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/tessim",
             array,"/",sep="")
for (ie in 1:length(EkeVcalib)){
    # file with npulsesPerEnergy for each energy; each pulse is pulselength long
    file <- paste("mono",EkeVcalib[ie],"_sep",separation,"_pix1_",nSimPulses,"p_",pulseLength,
                      "_pulse.txt",sep="")
    fileInit <- paste(dir,file,sep="")
    valuesTab <- read.table(fileInit, header=FALSE)
    valuesMat <- as.matrix(valuesTab[1:npulsesPerEnergy,1:size])
    valuesMat <- cbind(valuesMat, c(rep(EkeVcalib[ie],npulsesPerEnergy)))
    if(ie==1){
        svdMatrixT <- valuesMat
    }else{
        svdMatrixT <- rbind(svdMatrixT,valuesMat)    
    }
}
# From svdMatrixT create svdMatrix for SVD:
svdMatrix <- t(svdMatrixT)
# pulse1(1).....pulse1000(1)    pulse1(1).......pulse1000(1) ...... pulse1(1)........pulse1000(1) 
# pulse1(2).....pulse1000(2)    pulse1(2).......pulse1000(2) ...... pulse1(2)........pulse1000(2)     
# ..................................................................................................
# pulse1(4096)..pulse1000(4096) pulse1(4096)....pulse1000(4096) ... pulse1(4096).....pulse1000(4096)
#    e0 <..........> e0           e1 <............> e1                  en <..........> en
#
# plot pulses linearity
#-----------------------
#plot(svdMatrix[90,],svdMatrix[130,],cex=0.5, 
#     xlab="Sample #90", ylab="Sample #130", typ="b")
#drawLogPlotBox(xlimits=c(min(svdMatrix[90,]),max(svdMatrix[90,])),
#               ylimits=c(min(svdMatrix[130,]),max(svdMatrix[130,])),
#               x2limits=c(min(svdMatrix[90,]),max(svdMatrix[90,])),
#               y2limits=c(min(svdMatrix[130,]),max(svdMatrix[130,])),
#               logxy="xy", xlabel="Sample #90", 
#               ylabel="Sample #130",
#               naxes=c(T,T,T,T))
#points(svdMatrix[90,],svdMatrix[130,],typ="b")


# # take SVD and get "ncomps" right-singular vectors and plot SVD coefficients
#-----------------------------------------------------------------------------
#    S = U * SIGMA * V'
#    S ~ U * SIGMA_gus * V'    # https://en.wikipedia.org/wiki/Singular_value_decomposition
#
SVD <- svd(svdMatrix[1:nrow(svdMatrix)-1,], nu=ncomps, nv=ncomps) 
s_vals <- SVD$d # vector with singular values
SIGMA <- diag(s_vals)
SIGMA_gus <- diag(s_vals[1:ncomps])

# Right singular vectors
#
#      RSV1            RSV2         .......        RSVi             ...... 
#  comp1_pulse1    comp2_pulse1                 compi_pulse1        ......
#  comp1_pulse2    comp2_pulse2                 compi_pulse2        ......
#   ........         .........                   ..........         ......
#  comp1_pulseN    comp2_pulseN                 compi_pulseN        ......
#
#
rsv    <- SVD$v # matrix with right singular vectors (columns)
trsv <- t(rsv)
lsv    <- SVD$u
MatTrunc <- lsv %*% SIGMA_gus %*% t(rsv)
rsv1 <- rsv[,1]  
rsv2 <- rsv[,2]
plot(rsv1,rsv2, main="", pch=1,xlab="SVD component 1",ylab="SVD component 2", 
     col="red",typ="b")

# plot coefficientes versus energy
mat=matrix(data=1:ncomps,nrow=2,ncol=2,byrow = T)
layout(mat)
pred1 <- numeric(ncomps)
for (icom in 1:ncomps){
    RSV <-rsv[,icom]
    x <- svdMatrix[nrow(svdMatrix),]
    plot(x,RSV,xlab="Energy(keV)")
    eners <- svdMatrix[nrow(svdMatrix),]
    smsp <- smooth.spline(eners,RSV,nknots=10)
    lines(smsp,col="blue")
    
    # derivate of spline in Eo:
    Eo <- 1
    pred0 <- predict(smsp, x=Eo, deriv=0) # prediction of spline in Eo
    pred1[icom] <- predict(smsp, x=Eo, deriv=1) # derivative of spline in Eo
    yint <- pred0$y - (pred1[icom]$y*Eo)
    points(pred0, col=2, pch=19) 
    lines(x, yint + pred1[icom]$y*x, col=3)
}

# Calculate approx derivative of the pulses signal (to be used in column 1 of M matrix in 
# tangential filtering method)
deriv <- numeric(ncomps)
derivate <- numeric(size)
for (is in 1:size)
    for (ic in 1:ncomps){
        deriv[ic] = lsv[is,ic]*s_vals[ic]*pred1[ic]
    }
    derivate[is] = sum(deriv)
