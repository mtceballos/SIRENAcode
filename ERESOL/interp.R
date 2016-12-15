

data<-read.table("interp_data.dat")
xdata <- data[,1]
ydata <- data[,2]
alldata <- c(2,2.5,xdata)
allexpo <-exp(-alldata/2.)

dataInterp <- read.table("interp_spline.dat")
xinterp <- dataInterp[,1]
yinterp <- dataInterp[,2]

plot(xdata,ydata,xlim=c(2,8),ylim=c(0,0.4))
lines(xinterp,yinterp,col="blue")        # spline interpolation
lines(alldata,allexpo,lty=2,col="red")   # real function
points(c(2,2.5),c(exp(-1),exp(-2.5/2.)),pch=4,col="red")  # extrapolation true points

points(c(2,2.5),c(0.363777,0.285688),pch=5,col="cyan")

legend("topright", c("Data values","Cubic GSL interpolation",
                     "True exponential function (exp(-x/2)",
                     "True additional values following exp",
                     "Extrapolation GSL"),
       pch=c(1,NA,NA,4,5),col=c("black","blue","red","red","cyan"),
       lty=c(NA,1,2,NA,NA))
