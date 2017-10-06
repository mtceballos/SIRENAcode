# ##############     elmNN     #################################################
# Training and predict functions for SLFN ( Single Hidden-layer Feedforward Neural Networks ) using the ELM
# algorithm. ELM algorithm differs from the traditional gradient-based algorithms for very short training times 
# ( it doesn't need any iterative tuning, this makes learning time very fast ) and there is no need to set any 
# other parameters like learning rate, momentum, epochs, etc.

setwd("/home/ceballos/NN")
pdf(file="ELMNN.pdf",width=10,height=7.)    
rm(list = ls())
options(stringsAsFactors = FALSE)

library(elmNN)
library(ELMR)
library(JBTools)
EkeVs <- c(0.2, 0.5, 1., 2., 3., 4., 5., 6., 7., 8.)
EkeVs <- list("1"=0.2, "2"=0.5, "3"=1, "4"=2, "5"=3, "6"=4, "7"=5, "8"=6, "9"=7, "10"=8)
nEkeVs <- length(EkeVs)
#set.seed(100)

readDataYXs <- function(nsamples,ntrain,ntest){
    # reads files of data according to input number of rows(npulses) and columns(nsamples)
    # full data is 19997 pulses & 4096 samples
    # nsamples: initial samples for each pulse
    # ntrain: number of pulses per Energy for training purposes
    # ntest: number of pulses per Energy for testing purposes
    pulsesTotal <- 19997 # rows in files (number of pulses per file)
    samplesTotal <- 4096 # total samples  in file 
    stopifnot((ntrain+ntest) <= pulsesTotal)
    
    nosamples <- samplesTotal -nsamples # samples not to be read

    dataTrain <- data.frame()
    dataTest  <- data.frame()
    #for (EkeV in EkeVs){
    for (EkeVCatg in names(EkeVs)){
        Eval <- EkeVs[[EkeVCatg]]
        file <- paste("mono",Eval,"_sep40000_1000-5095.txt",sep="")
        if(Eval < 1)file <- paste("mono",Eval,"_sep40000_999-5094.txt",sep="")
        colclasses <- c(rep("numeric",nsamples),c(rep("NULL",nosamples)))
        
        dtr <- read.csv(file, header=FALSE, colClasses=colclasses, nrows=ntrain)
        dtr <- cbind(EkeVCatg=EkeVCatg,dtr)
        dtr <- cbind(EkeV=Eval,dtr)
        dataTrain <- rbind(dataTrain, dtr)
        
        dts <- read.csv(file, header=FALSE, colClasses=colclasses, skip=ntrain, nrows=ntest)
        dts <- cbind(EkeVCatg=EkeVCatg,dts)
        dts <- cbind(EkeV=Eval,dts)
        dataTest <- rbind(dataTest, dts)
    }
    return(list(train=dataTrain,test=dataTest))
}

nsamples <- 300 #300
npulsesToTrain <- 300 #100
npulsesToTest <- 100 #300
npulses <- npulsesToTrain + npulsesToTest #input pulses per energy value #20
FWHMtest.ELMNN <- numeric()
RMSEtest.ELMNN <- numeric()
FWHMtest.ELM <- numeric()

# randomly separate train and test data 
#index <- sample(1:nrow(data),round(0.75*nrow(data)))
#train <- data[index,]
#test <- data[-index,]
train <- readDataYXs(nsamples,npulsesToTrain,npulsesToTest)$train
test <- readDataYXs(nsamples,npulsesToTrain,npulsesToTest)$test
# preProcess data for ELMR methods
train_ <- data.frame()
train_ <-data.frame(preProcess(train[,3:(nsamples+2)]))
train_ <- cbind(EkeVCatg=train$EkeVCatg,train_)
train_ <- cbind(EkeV=train$EkeV,train_)
test_ <- data.frame()
test_ <-data.frame(preProcess(test[,3:(nsamples+2)]))
test_ <- cbind(EkeVCatg=test$EkeVCatg,test_)
test_ <- cbind(EkeV=test$EkeV,test_)

n <- names(train)
freg <- as.formula(paste("EkeV ~", paste(n[!n %in% c("EkeV","EkeVCatg")], collapse = " + ")))
fclas <- as.formula(paste("EkeVCatg ~", paste(n[!n %in% c("EkeV","EkeVCatg")], collapse = " + ")))

cat("Doing ELMNN/regression\n")
# ELMNN - regression
modelELMNN = elmtrain(freg, train, nhid=300, actfun="poslin")
predELMNN <- predict.elmNN(modelELMNN,newdata=test) # same as predict(modelELMNN,newdata=test)

cat("Doing ELMN/regression\n")
# ELM - regression
modelELM_r = OSelm_train.formula(freg, train_, "regression", 300, "sig", 10, 10)
predELM = predict_elm(modelELM_r, test_)
predELM_r = predELM$predicted
accuracy_r <-format(predELM$testAccuracy, digits=3)

cat("Doing ELMN/classification\n")
# ELM - classification
modelELM_c = OSelm_train.formula(fclas, train_, "classification", 200, "sig", 10, 10)
predELM = predict_elm(modelELM_c, test_)
predELM_c = predELM$prediction
accuracy_c <- format(predELM$accuracy.reg, digits=3)

for (EkeV in EkeVs){
    fwhm <- 1E3*2.35*sd(predELMNN[test[,1] == EkeV]) #eV
    rmse <- RMSE(predELMNN[test[,1] == EkeV], test[test[,1] == EkeV,1])
    FWHMtest.ELMNN <- append(FWHMtest.ELMNN, fwhm) 
    RMSEtest.ELMNN <- append(RMSEtest.ELMNN, rmse) 
    
    fwhm <- 1E3*2.35*sd(predELM_r[test[,1] == EkeV]) #eV
    FWHMtest.ELM <- append(FWHMtest.ELM, fwhm)
    
}
par(mfrow=c(2,2))
par(pty="s")
plot(test[,1],predELMNN, xlab="Test Energies (keV)", ylab="Predicted energies (keV)")
title(main=paste("ELMNN-regression"))
abline(0,1)
plot(EkeVs,FWHMtest.ELMNN,xlab="Test Energies (keV)", ylab="FWHM(eV)")
title(main=bquote("FWHM.ELMNN (eV)=2.35*" * sqrt(frac(sum((E[i]-bar(E))^2),.(npulsesToTest)))), cex.main=0.7)
#plot(EkeVs,RMSEtest.ELMNN,xlab="Test Energies (keV)", ylab="RMSE")
#title(main="Residual mean Square Error", cex.main=0.7)

plot(test[,1],predELM_r, xlab="Test Energies (keV)", ylab="Predicted Energy (keV)")
title(main=paste("EMLR-regression\n Accuracy=",accuracy_r,sep=""))
abline(0,1)
#plot(EkeVs,FWHMtest.ELM,xlab="Test Energies (keV)", ylab="FWHM(eV)")
#title(main=bquote("FWHM.ELM (eV)=2.35*" * sqrt(frac(sum((E[i]-bar(E))^2),.(npulsesToTest)))), cex.main=0.7)

plot(as.character(test[["EkeVCatg"]]),predELM_c, xlab="Test categories", ylab="Predicted categories")
title(main=paste("EMLR-classification\n  Accuracy=",accuracy_c,sep=""))
abline(0,1)

dev.off()
