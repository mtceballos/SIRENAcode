Data3000<-read.table("sep3000sam_100s_3keV_t1.txt",header=F)
Data0003<-read.table("sep0003sam_100s_3keV_t1.txt",header=F)
plot(seq(1:10000),Data[1,], type="l")

base <- as.numeric(Data3000[1,5001:10000])
baseline <- c(base,base)

pulse3 <- c(rep(0,499),as.numeric(Data3000[1,500:2500]),rep(0,7500))

#pulse3nb1 <- c(rep(0,1002),as.numeric(Data3000[1,1003:2500])-130,rep(0,7500))
#pulse3nb2 <- c(rep(0,3002),as.numeric(Data3000[1,1003:2500])-130,rep(0,5500))

#                                 1500 samples 
pulse3nb1 <- c(rep(0,1002),as.numeric(Data3000[1,1003:2500])-130,rep(0,7500))
#pulse3nb2 <- c(rep(0,1005),as.numeric(Data3000[1,1003:2500])-130,rep(0,7497))
pulse3nb2 <- c(rep(0,1005),as.numeric(Data3000[1,1003:2500])-130,rep(0,7497))
pulse3nb1Deriv <- diff(pulse3nb1)

#plot(seq(500:2500),pulse3,type="l", col="blue")
pulses3pair <- baseline +pulse3nb1 + pulse3nb2
plot(seq(1:10000),pulses3pair,type="l", col="blue",xlim=c(1000,3100))
pulses3pairDeriv <- diff(pulses3pair)/diff(seq(1:10000))
plot(seq(1:9999),pulses3pairDeriv,col="red",type="l",xlim=c(1000,3000), ylim=c(-500,1800))

resta <-pulses3pairDeriv - pulse3nb1Deriv

Data0003Deriv <- diff(as.numeric(Data0003))

#plot(seq(1:10000),as.numeric(Data0003[1,]), type="l",col="blue",xlim=c(1001,1006), ylim=c(0,7000))
#lines(seq(1:10000),pulses3pair,col="red")
plot(seq(1:9999),Data0003Deriv, type="l",col="blue",xlim=c(1001,1006))
lines(seq(1:9999),pulses3pairDeriv,col="red")
