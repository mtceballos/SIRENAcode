data<-read.csv("philippe300.csv")
mus <- data[,1]
fwhm<- data[,2]
samprate<-156250

samples<-mus*1E-6*samprate
