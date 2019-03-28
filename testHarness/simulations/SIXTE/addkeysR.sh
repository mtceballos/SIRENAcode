#!/bin/sh

#
# Shell script to add keywords for XIFUSIM
#

fits0=${1}+0
fits1=${1}+1

R0=7.12e-3 
fparkey value=${R0} fitsfile=${fits1}  keyword='R0' add=y
I0_START=13.56e-6 
fparkey value=${I0_START} fitsfile=${fits1}  keyword='I0_START' add=y
ADUCNV=1.9142e-10 
fparkey value=${ADUCNV} fitsfile=${fits1}  keyword='ADUCNV' add=y
RPARA=1e-3
fparkey value=${RPARA} fitsfile=${fits1}  keyword='RPARA' add=y
TTR=0.9144 
fparkey value=${TTR} fitsfile=${fits1}  keyword='TTR' add=y
LFILTER=1e-6
fparkey value=${LFILTER} fitsfile=${fits1}  keyword='LFILTER' add=y


