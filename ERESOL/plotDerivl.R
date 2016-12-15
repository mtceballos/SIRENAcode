#
#  Plot derivatives of pulses & threshold for detection
#
deriv1<-read.table("~/INSTRUMEN/EURECA/ERESOL/record4_1keV_sep0001.txt")
deriv2<-read.table("~/INSTRUMEN/EURECA/ERESOL/record4_1keV_sep0002.txt")
deriv3<-read.table("~/INSTRUMEN/EURECA/ERESOL/record4_1keV_sep0003.txt")
deriv5<-read.table("~/INSTRUMEN/EURECA/ERESOL/record4_1keV_sep0005.txt")
threshold=198.05
indices<-deriv1[,1]+1

#pdf("~/INSTRUMEN/EURECA/ERESOL/PAIRS/tessim20150505/derivPlot.pdf",width=10.,height=7.1) 
pdf("pp.pdf")
plot(indices,deriv1[,2], xlab="Samples", ylab="Pulse Stream Derivative",xlim=c(999,1010),type="b",col="red",
     main="Two pulses Derivative",pch=19)
points(indices,deriv3[,2],col="magenta",typ="b",pch=19)
points(indices,deriv5[,2],col="blue",typ="b",pch=19)
points(indices,deriv2[,2],col="green",typ="b",pch=19)

# pulses at 1002 & 1003
points(x=c(1002,1003),y=c(deriv1[1002,2],deriv1[1003,2]),pch=1,col="red",cex=2)
points(x=c(1002,1004),y=c(deriv2[1002,2],deriv2[1004,2]),pch=1,col="green",cex=2)
points(x=c(1002,1005),y=c(deriv3[1002,2],deriv3[1005,2]),pch=1,col="magenta",cex=2)
points(x=c(1002,1007),y=c(deriv5[1002,2],deriv5[1007,2]),pch=1,col="blue",cex=2)
abline(h=threshold,lty=2)

legend("bottomright",c("sep0001, 1keV, SPA",
                       "sep0002, 1keV, SPA",
                       "sep0003, 1keV, SPA",
                       "sep0005, 1keV, SPA"),
       col=c("red","green","magenta","blue"),pch=c(19,19),cex=0.8)

slope1_0001<-deriv1[1003,2]-deriv1[1002,2]
slope1_0002<-deriv2[1003,2]-deriv2[1002,2]
slope1_0003<-deriv3[1003,2]-deriv3[1002,2]
slope1_0005<-deriv5[1003,2]-deriv5[1002,2]

slope2_0001<-deriv1[1004,2]-deriv1[1003,2]
slope2_0002<-deriv2[1004,2]-deriv2[1003,2]
slope2_0003<-deriv3[1004,2]-deriv3[1003,2]
slope2_0005<-deriv5[1004,2]-deriv5[1003,2]


text(x=1003.5,y=deriv1[1003,2],srt=39,labels="slope2")
text(x=1002.5,y=2500,srt=70,labels="slope1")

percent1<-round(slope2_0001/slope1_0001*100)
percent2<-round(slope2_0002/slope1_0002*100)
percent3<-round(slope2_0003/slope1_0003*100)
percent5<-round(slope2_0005/slope1_0005*100)
legend("topleft",c(paste("slope2=",percent1,"% of slope1",sep=""),
                   paste("slope2=",percent2,"% of slope1",sep=""),
                   paste("slope2=",percent3,"% of slope1",sep=""),
                   paste("slope2=",percent5,"% of slope1",sep="")),
                   text.col=c("red","green","magenta","blue"))

dev.off()
