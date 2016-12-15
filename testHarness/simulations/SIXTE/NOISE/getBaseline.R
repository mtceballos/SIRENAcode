#!/usr/bin/env Rscript
#
# fit a line to a stream (noise) file to get constant baseline
# Script return baseline value
# Run as getBaseline.R inputTXTfile valuesColNum

args <- commandArgs(trailingOnly = TRUE)
#args<-c("./tessimSPA/forBaseline1024samples_tessimSPA_10s_R1.txt",2)
#cat("Running getBAseline.R with ",args[1]," and ", args[2],"\n")
streamData <- read.table(args[1],header=T)
npoints <- nrow(streamData)
time <- streamData[500:npoints,1]
ADC  <- streamData[500:npoints,as.numeric(args[2])]

ajuste <- lm(ADC~time)
baseline <- summary(ajuste)$coefficients[1]
baseline2 <- ajuste$coefficients[1]

#cat("Baseline=",baseline)
cat(baseline)
