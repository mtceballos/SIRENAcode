#!/bin/csh
#
# NOISE spectrum simulation in Resistance Space: baseline can be given as input parameter of fitted
#
# source simulNoiseResist.csh [baseline]
#
#  Input parameters: 
#          tessim (see possibilities/examples below)
#          [baseline]

#
# WARNING: Do not initialize HEASOFT!!!!!!!!  Only SIXTE
#
if($#argv < 2) then 
    echo "Please complete input command line with all the parameters"
    echo "   source simulNoiseResist.csh array pulseLength [baseline]"
    exit
endif 
set array=$1
set tessim="tessim$1"
set pulseLength=$2  # also noise samples
set baseline=0

if($#argv >2) then
    set baseline=$3
endif

sixteinit
set XMLrect=$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/CALIBRATION/newFormat/xifu_detector_300us_156kHz_rect.xml
set pixel=1
set samprate=156250 #Hz
set simTimeB=10
set simTimeN=100

set RSHUNT=91E-6 #Ohm
set ADUCNV=7.6311532944731E-10 #A/ADU
set IBIAS=7.35E-05

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
set noiseFile="noise${pulseLength}samples_${tessim}_B0_${simTimeN}s_${cpsN}cps_R1.fits"

# get baseline
#--------------
if($3 == "") then
	echo "\nPIXDETILLUM: Generating no events file for BASELINE \n"
	#pixdetillum PixImpList=${root1}.piximpact XMLFile=$XMLrect tstart=0. Tstop=$simTimeB pixels=$pixel rate=0 energy=0 clobber=yes seed=-1

	echo "\nTESSIM: Generating stream file for BASELINE \n"      
	#tessim PixID=$pixel PixImpList=${root1}.piximpact Streamfile=${root1}.stream tstart=0. tstop=$simTimeB sample_rate=$samprate clobber=yes PixType="$PixType"
	echo "\nFCAL: transform stream to Resitance space\n"
	(heainit;fcalc infile=${root1}.stream outfile=${root1}_R1.stream clobber=yes clname=RESIST expr="${RSHUNT}*(#I0_START/(PXL00001*#ADUCNV)-1)")
	
	echo "\nFDUMP: Dumping stream file for BASELINE\n"
	(heainit;fdump infile=${root1}_R1.stream outfile=${root1}_R1.txt prhead=no showrow=no showunit=no clobber=yes showcol=no columns="-" rows="-")
	echo "Running ${pwd}/NOISE/getBaseline.R ${root1}_R1.txt 4"
	set baseline=`${pwd}/NOISE/getBaseline.R ${root1}_R1.txt 4`
	#rm ${root1}.stream
	#rm ${root1}.txt
endif
echo "############################\n"
echo "    Baseline=$baseline"
echo "\n############################"

#exit
#process empty trigger file
echo "\nTESCONSTPILEUP: Generating no events file for NOISE\n"
##pixdetillum PixImpList=${root2}.piximpact XMLFile=$XMLrect tstart=0. Tstop=$simTimeN pixels=$pixel rate=$cpsN energy=3 clobber=yes seed=-1
#tesconstpileup PixImpList=${root2}.piximpact XMLFile=$XMLrect tstop=$simTimeN energy=0 pulseDistance=1 TriggerSize=10000 clobber=yes
# 
echo "\nTESSIM: Generating stream file for NOISE\n"            
tessim PixID=$pixel PixImpList=${root2}.piximpact Streamfile=${root2}.stream tstart=0. tstop=$simTimeN sample_rate=$samprate clobber=yes PixType="$PixType"

echo "\nFCAL: transform stream to Resitance space\n"
(heainit;fcalc infile=${root2}.stream outfile=${root2}_R1.stream clobber=yes clname=PXL00001 expr="${RSHUNT}*(#I0_START/(PXL00001*#ADUCNV)-1)")

echo "\nSTREAMTOTRIGGERS: Generating trigger file for NOISE\n"      
streamtotriggers PixImpList=${root2}.piximpact XMLFile=$XMLrect tstart=0. tstop=$simTimeN Streamfile=${root2}_R1.stream TesTriggerFile=${root2}_R1.fits TriggerSize=10000 PreBufferSize=1000 pixels=$pixel clobber=yes
(heainit;cphead ${root2}_R1.stream+1 ${root2}_R1.fits+1)
#rm ${root2}.stream
#rm first record (noisy)
(heainit;fdelrow infile=${root2}_R1.fits+1 firstrow=1 nrows=1 confirm="no" proceed=yes)

echo "\nGENNOISESPEC: Generating NOISE spectrum file \n"
#echo "Command: gennoisespec --inFile=${root2}.fits --outFile=$noiseFile --intervalMinSamples=$pulseLength --nintervals=600 --scaleFactor=$scaleFactor --samplesUp=$samplesUp --nSgms=$nSgms --pulse_length=$pulseLength --baseline=$baseline --ntausPF=0 --clobber=yes verbosity=0 --tauFALL=3E-5"

gennoisespec --inFile=${root2}_R1.fits --outFile=$noiseFile --intervalMinSamples=$pulseLength --nintervals=1000 --scaleFactor=$scaleFactor --samplesUp=$samplesUp --nSgms=$nSgms --pulse_length=$pulseLength --baseline=$baseline --ntausPF=0 --clobber=yes verbosity=0 
(heainit;fparkey fitsfile=$noiseFile value=$baseline keyword=BASELINE add=yes)

cd $pwd



  
  
