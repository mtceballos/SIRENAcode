#
# Plots in 3D weight and covariance matrices for pulses library
# Crear el fichero de datos con:
#  fdump ../testHarness/simulations/SIXTE/LIBRARIES/tessimLPA1shunt/GLOBAL/ADC/libraryMultiE_GLOBAL_PL2048_tessimLPA1shunt.fits clobber=yes prhead=no showcol=no showrow=no
#   Name of optional output file[LPA1shunt_COVAR_1keV.txt] 
#   Names of columns[COVARM] 
#   Lists of rows[3] 

#
setwd("/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/tessimLPA1shunt/GLOBAL/RFITTED")
pulseLength <- 1024
#energy<-"1KeV08072015"
energy<-"7keV"
array <- "LPA1shunt"
#plot<-"dataW"
plot<-"dataC"
#plot<-"noiseW"
#plot<-"noiseC"
#noise<-"nonoise"
noise<-""
nx<-pulseLength
#nx<-100 # to test plotting with smaller numbers
nx2<-nx*nx
x<-seq(1:nx)
y<-seq(1:nx)
# persp3d (to be able to use rotation, rgl)
library(rgl)
library(evd)
xlabel <- "samples"
ylabel <- "samples"
nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))

if(plot == "dataW"){
    title<-paste("Weight Matrix ",energy)
    dataWEIGHT<-read.table(paste("./SPA_WEIGHT_",energy,".txt",sep=""),header=F)    
    matWEIGHT <- matrix(data=dataWEIGHT[,1], nrow=pulseLength, ncol=pulseLength, byrow=T)
    # (subsample) data
    dataWEIGHTSub <- dataWEIGHT[1:nx2,1]/10E18
    matWEIGHTSub  <- matrix(data=dataWEIGHTSub,nrow=nx,ncol=nx,byrow=T)
    
    zcol  = cut(matWEIGHTSub, nbcol)
    persp3d(x,y,z=matWEIGHTSub,main=title,
            shade=0.3,image=T,col=color[zcol],
            xlab=xlabel, ylab=ylabel,zlab="weights", specular="black")
}else if(plot == "dataC"){
    #title<-paste("Covariance Matrix (",noise,") ",energy)
    title<-paste("Covariance Matrix ",energy)
    file<-paste("./LPA1shunt_COVAR_",energy,"_",pulseLength,".txt",sep="")
    cat("Reading:",file)
    dataCOV   <-read.table(file,header=F)
    matCOV    <- matrix(data=dataCOV[,1],    nrow=pulseLength, ncol=pulseLength, byrow=T)
    dataCOVSub <- dataCOV[1:nx2,1]
    matCOVSub  <- matrix(data=dataCOVSub,nrow=nx,ncol=nx,byrow=T)
    
    zcol  = cut(matCOVSub, nbcol)
    persp3d(x,y,z=matCOVSub,main=title,
            shade=0.3,image=T,col=color[zcol],
            xlab=xlabel, ylab=ylabel,zlab="covariance", specular="black")
    
}else if(plot == "noiseW"){
    title<-paste("Weight Noise Matrix ",energy)
    dataWEIGHTnoise<-read.table(paste("./SPAnoise1024_WEIGHT_1KeV.txt",sep=""),header=F)
    matWEIGHTnoise <- matrix(data=dataWEIGHTnoise[,1], nrow=pulseLength, 
                             ncol=pulseLength, byrow=T)   
    dataWEIGHTnoiseSub <- dataWEIGHTnoise[1:nx2,1]/10E18
    matWEIGHTnoiseSub  <- matrix(data=dataWEIGHTnoiseSub,nrow=nx,ncol=nx,byrow=T)
    
    zcol  = cut(matWEIGHTnoiseSub, nbcol)
    persp3d(x,y,z=matWEIGHTnoiseSub,main=title,
            shade=0.3,image=T,col=color[zcol],
            xlab=xlabel, ylab=ylabel,zlab="weights", specular="black")
    
}else if(plot == "noiseC"){
    #title<-paste("Covariance Noise Matrix ",energy)
    title<-paste("Covariance Noise Matrix ")
    dataCOVnoise   <-read.table(paste("./SPAnoise1024_COVAR_1KeV.txt",sep=""),header=F)
    matCOVnoise    <- matrix(data=dataCOVnoise[,1],    nrow=pulseLength, ncol=pulseLength, byrow=T)   
    dataCOVnoiseSub <- dataCOVnoise[1:nx2,1]
    matCOVnoiseSub  <- matrix(data=dataCOVnoiseSub,nrow=nx,ncol=nx,byrow=T)
    
    zcol  = cut(matCOVnoiseSub, nbcol)
    persp3d(x,y,z=matCOVnoiseSub,main=title,
            shade=0.3,image=T,col=color[zcol],
            xlab=xlabel, ylab=ylabel,zlab="covariance", specular="black")
}



# Diffrent tests...

# Image ??
#=========
#minCov<--20
#maxCov<-+20
#ticks<-10**(seq(from=log10(minCov),to=log10(maxCov),length.out=10))
#ticksLabels<-sprintf("%2.1f",ticks)
#image.plot(x,y,log10(matSmall),axis.args=list( at=log10(ticks), labels=ticksLabels),
#           zlim=c(log10(minCov),log10(maxCov)))
#image.plot(x,y,matSmall, zlim=c(minCov,maxCov))

# persp & persp3D ???
#====================
#persp(seq(1:nx),seq(1:nx),matSmall)
#par(mfrow=c(2,2))
# persp3D(z=matSmall,main="Covariance 1 keV")
# persp3D(z=matSmall,main="Covariance 1 keV",shade=0.3,theta=0,image=T)
# persp3D(z=matSmall,main="Covariance 1 keV",shade=0.3,theta=0,phi=10,image=T)
# persp3D(z=matSmall,main="Covariance 1 keV",shade=0.3,theta=30,image=T)

