# PLOT limits of grading table for XIFU and pre-XIFU

#par(mfrow=c(2,2))
#pdf(file="gradinTable.pdf",width=10,height=7.)

mat=matrix(data=c(1:4),nrow=2,ncol=2,byrow = T)
layout(mat)
arrays<-c("SPA","LPA1","LPA2","LPA3")
gradingTable <- read.table("~/INSTRUMEN/EURECA/ERESOL/gradingTable.dat",header=T)
maxsep<-50000
minsep<-10

for (array in arrays){
    # read resolution breaks
    nsamples<-gradingTable[gradingTable[,1]==array,"nsamples"]
    biasBreak<-gradingTable[gradingTable[,1]==array,"breakBIAS"]
    HRbreak<-gradingTable[gradingTable[,1]==array,"breakHIGHRES"]
    MRbreak<-gradingTable[gradingTable[,1]==array,"breakMIDRES"]

    # plot grading as Astro-H/XRS classification but with X-IFU values
    x<-exp(log(seq(minsep,maxsep,length.out = 200)))
    y<-exp(log(seq(minsep,maxsep,length.out = 200)))
    labelx="Time to previous pulse (samples)"
    labely="Time to next pulse (samples)"
    #plot(x,y,log="xy",type="n",xaxp=c(min(x),max(x),3),main="SPA",xlab=labelx, ylab=labely)
    plot(x,y,log="xy",type="n",main=array,xlab=labelx, ylab=labely,cex=0.8,asp=1)

    # drawLogPlotBox(xlimits=c(min(x),max(x)),ylimits=c(min(y),max(y)),
    #                logxy="xy",naxes = c(T,T,F,F))
    rect(HRbreak,HRbreak,maxsep,maxsep,border = "transparent",col = "blue") #HR
    rect(HRbreak,MRbreak,maxsep,HRbreak,border = "transparent",col = "green") # Mp
    rect(MRbreak,MRbreak,HRbreak,HRbreak,border = "transparent",col = "green3") # Ms
    rect(HRbreak,minsep,maxsep,MRbreak,border = "transparent",col = "pink") # Lp
    rect(minsep,minsep,HRbreak,MRbreak,border = "transparent",col = "pink3") # Ls
    rect(minsep,MRbreak,MRbreak,maxsep,border = "transparent",col = "pink3") # Ls
    rect(MRbreak,HRbreak,HRbreak,maxsep,border = "transparent",col = "green3") # Ms
    abline(v=MRbreak,lty=2)
    abline(v=HRbreak,lty=2)
    abline(h=MRbreak,lty=2)
    abline(h=HRbreak,lty=2)
    mtext(text=MRbreak,side=1,at=MRbreak,cex=0.7)
    mtext(text=HRbreak,side=1,at=HRbreak,cex=0.7)
    mtext(text=MRbreak,side=2,at=MRbreak,cex=0.7,las=1)
    mtext(text=HRbreak,side=2,at=HRbreak,cex=0.7,las=1)
    # HR labels
    text(x=10**((log10(max(x))-log10(HRbreak))/2+log10(HRbreak)), 
         y=10**((log10(max(y))-log10(HRbreak))/2+log10(HRbreak)),
         labels = "HR",cex=1.2)
    text(x=10**((log10(max(x))-log10(HRbreak))/3+log10(HRbreak)), 
         y=10**((log10(max(y))-log10(HRbreak))/3+log10(HRbreak)),
         labels = "Full",col="black")
    text(x=10**(2*(log10(max(x))-log10(HRbreak))/3+log10(HRbreak)), 
         y=10**((log10(max(y))-log10(HRbreak))/3+log10(HRbreak)),
         labels = "Full",col="white")
    # MR labels
    text(x=10**((log10(max(x))-log10(HRbreak))/2+log10(HRbreak)), 
         y=10**((log10(HRbreak)-log10(MRbreak))/2+log10(MRbreak)),
         labels = "Mp")
    text(x=10**((log10(max(x))-log10(HRbreak))/4+log10(HRbreak)), 
         y=10**((log10(HRbreak)-log10(MRbreak))/3+log10(MRbreak)),
         labels = "1/4 Full")
    text(x=10**(3*(log10(max(x))-log10(HRbreak))/4+log10(HRbreak)), 
         y=10**((log10(HRbreak)-log10(MRbreak))/3+log10(MRbreak)),
         labels = "VarLength",col="white")
    
    text(x=10**((log10(HRbreak)-log10(MRbreak))/2+log10(MRbreak)), 
         y=10**((log10(HRbreak)-log10(MRbreak))/2+log10(MRbreak)),
         labels = "Ms")
    text(x=10**((log10(HRbreak)-log10(MRbreak))/4+log10(MRbreak)), 
         y=10**((log10(HRbreak)-log10(MRbreak))/3+log10(MRbreak)),
         labels = "1/4",cex=0.8)
    text(x=10**(3*(log10(HRbreak)-log10(MRbreak))/4+log10(MRbreak)), 
         y=10**((log10(HRbreak)-log10(MRbreak))/3+log10(MRbreak)),
         labels = "VL",cex=0.8,col="white")
    
    text(x=10**((log10(HRbreak)-log10(MRbreak))/2+log10(MRbreak)), 
         y=10**((log10(max(y))-log10(HRbreak))/2+log10(HRbreak)),
         labels = "Ms")
    text(x=10**((log10(HRbreak)-log10(MRbreak))/4+log10(MRbreak)), 
         y=10**((log10(max(y))-log10(HRbreak))/3+log10(HRbreak)),
         labels = "1/4",cex=0.8)
    text(x=10**(3*(log10(HRbreak)-log10(MRbreak))/4+log10(MRbreak)), 
         y=10**((log10(max(y))-log10(HRbreak))/3+log10(HRbreak)),
         labels = "Full",col="white",cex=0.8)
    
    # LR labels
    text(x=10**((log10(max(x))-log10(HRbreak))/2+log10(HRbreak)), 
         y=10**((log10(MRbreak)-log10(min(y)))/2+log10(min(y))),
         labels = "Lp")
    
    text(x=10**((log10(HRbreak)-log10(MRbreak))/2+log10(MRbreak)), 
         y=10**((log10(MRbreak)-log10(min(y)))/2+log10(min(y))),
         labels = "Ls")
    text(x=10**((log10(MRbreak)-log10(min(x)))/2+log10(min(x))), 
         y=10**((log10(HRbreak)-log10(MRbreak))/2+log10(MRbreak)),
         labels = "Ls")
    text(x=10**((log10(MRbreak)-log10(min(x)))/2+log10(min(x))), 
         y=10**((log10(MRbreak)-log10(min(y)))/2+log10(min(y))),
         labels = "Ls")
    text(x=10**((log10(MRbreak)-log10(min(x)))/2+log10(min(x))), 
         y=10**((log10(max(y))-log10(HRbreak))/2+log10(HRbreak)),
         labels = "Ls")
    
    # plot lines for X-IFU classification in grading table
    abline(v=biasBreak,col="white",lw=2)
    segments(biasBreak,HRbreak,x1=max(x),col="white",lw=2)
    segments(biasBreak,MRbreak,x1=max(x),col="white",lw=2)
}

#dev.off()

