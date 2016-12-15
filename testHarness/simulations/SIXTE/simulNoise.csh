#!/bin/csh
#
# NOISE spectrum simulation: baseline can be given as input parameter of fitted
#
# source simulNoise.csh [baseline]
#
#  Input parameters: 
#          tessim (see possibilities/examples below)
#          [baseline]

#
# WARNING: Do not initialize HEASOFT!!!!!!!!  Only SIXTE
#
alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'

if($#argv < 3) then 
    echo "Please complete input command line with all the parameters"
    echo "   source simulNoise.csh array pulseLength ADC/R [baseline]"
    exit
endif 
set array=$1
set tessim="tessim$1"
set pulseLength=$2  # also noise samples
set space=$3
set baseline=0

if($#argv >3) then
    set baseline=$4
endif

sixteinit
set XMLrect=$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/CALIBRATION/newFormat/xifu_detector_300us_156kHz_rect.xml
set pixel=1
set pixelStr=`printf "%05d" $pixel`
set dataCol="PXL$pixelStr"
if($space == "R") set dataCol="RESIST"
set samprate=156250 #Hz
set simTimeB=10     # s of simulation for baseline calculations
set simTimeN=100    # s of simulation for noise spectra calculations
set RSHUNT=91E-6 #Ohm

set PixType="file:/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/tespixels.fits[${array}]"

set scaleFactor=0.005
set samplesUp=2
set nSgms=20
set cpsN=pairs
set pwd=`pwd`
#
# Create NOISE file
# 
cd NOISE/$tessim
set root1="forBaseline${pulseLength}samples_${tessim}_${simTimeB}s"
set root2="forNoise${pulseLength}samples_${tessim}_${simTimeN}s_${cpsN}cps"
#set root="forNoiseBaseline_${tessim}_${simTimeN}s"
set noiseFile="noise${pulseLength}samples_${tessim}_B0_${simTimeN}s_${cpsN}cps_${space}.fits"

# get baseline
#--------------
if($4 == "") then
	echo "\nPIXDETILLUM: Generating no events file for BASELINE \n"
	pixdetillum PixImpList=${root1}.piximpact XMLFile=$XMLrect tstart=0. Tstop=$simTimeB pixels=$pixel rate=0 energy=0 clobber=yes seed=-1

	echo "\nTESSIM: Generating stream file for BASELINE \n"      
	tessim PixID=$pixel PixImpList=${root1}.piximpact Streamfile=${root1}.stream tstart=0. tstop=$simTimeB sample_rate=$samprate clobber=yes PixType="$PixType"

	if ($space == "R")then
	    echo "\nFCALC: Adding Resistance column for BASELINE \n"
	    (heainit;fcalc infile=${root1}.stream outfile=${root1}.stream clobber=yes clname=RESIST expr="${RSHUNT}*(#I0_START/(${dataCol}*#ADUCNV)-1)")
	    set dataCol="RESIST"
	endif
	    
	echo "\nFDUMP: Dumping stream file for BASELINE\n"
	(heainit;fdump infile=${root1}.stream outfile=${root1}.txt prhead=no showrow=no clobber=yes showcol=no columns="TIME,${dataCol}" rows="-")

	set baseline=`../getBaseline.R ${root1}.txt 2`
	rm ${root1}.stream
	rm ${root1}.txt
endif
echo "############################\n"
echo "    Baseline=$baseline"
echo "\n############################"

#process empty trigger file
#echo "\nPIXDETILLUM: Generating low events file for NOISE\n"
##pixdetillum PixImpList=${root2}.piximpact XMLFile=$XMLrect tstart=0. Tstop=$simTimeN pixels=$pixel rate=$cpsN energy=3 clobber=yes seed=-1
#tesconstpileup PixImpList=${root2}.piximpact XMLFile=$XMLrect tstop=$simTimeN energy=0 pulseDistance=1 TriggerSize=10000 clobber=yes
# 
echo "\nTESSIM: Generating stream file for NOISE\n"            
tessim PixID=$pixel PixImpList=${root2}.piximpact Streamfile=${root2}.stream tstart=0. tstop=$simTimeN sample_rate=$samprate clobber=yes PixType="$PixType"

echo "\nSTREAMTOTRIGGERS: Generating trigger file for NOISE\n"      
streamtotriggers PixImpList=${root2}.piximpact XMLFile=$XMLrect tstart=0. tstop=$simTimeN Streamfile=${root2}.stream TesTriggerFile=${root2}.fits TriggerSize=10000 PreBufferSize=1000 pixels=$pixel clobber=yes
(heainit;cphead ${root2}.stream+1 ${root2}.fits+1)
rm ${root2}.stream
#rm first record (noisy)
(heainit;fdelrow infile=${root2}.fits+1 firstrow=1 nrows=1 confirm="no" proceed=yes)
if ($space == "R")then
    echo "\nFCALC: Adding Resistance column for NOISE \n"
    (heainit;fcalc infile=${root2}.fits outfile=${root2}.fits clobber=yes clname=ADC expr="${RSHUNT}*(#I0_START/(ADC*#ADUCNV)-1)")
endif
echo "\nGENNOISESPEC: Generating NOISE spectrum file \n"
gennoisespec --inFile=${root2}.fits --outFile=$noiseFile --intervalMinSamples=$pulseLength --nintervals=1000 --scaleFactor=$scaleFactor --samplesUp=$samplesUp --nSgms=$nSgms --pulse_length=$pulseLength --baseline=$baseline --ntausPF=0 --clobber=yes verbosity=0 
(heainit;fparkey fitsfile=$noiseFile value=$baseline keyword=BASELINE add=yes)

cd $pwd



  
  
