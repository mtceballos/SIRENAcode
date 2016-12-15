#!/bin/csh
#
#  Input parameters: 
#          tessim (see possibilities/examples below)
#          
# Unset variables to be used in script
#
unset calibDir XMLrect pixel samprate tessim cps simTime XMLFile monoen root 

if($#argv < 1) then 
    echo "Please complete input command line with all the parameters"
    echo "   source simulPairs.csh tessim"
    exit
endif    

set pixel=1
set samprate=156250 #Hz
set pulseLength=1024
#set tessim="tessim20150505"
set tessim=$1
if($tessim =~ *m0) set m_unknown="m_unknown=0"
set cps=10 
set simTime=300
set calibDir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/CALIBRATION/pulsetempl_v16092014/156khz"
set PIXIMPACTdir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/CALIBRATION/PIXIMPACT"
set XMLrect="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/CALIBRATION/newFormat/xifu_detector_300us_156kHz_rect.xml";
set XMLfile=$XMLrect
set libDir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/${tessim}"
set libFile=${libDir}/libraryMultiE_PL${pulseLength}_${tessim}.fits

set pwd=`pwd`
cd $tessim

    
#
# SIMULATE files for calibration
#
foreach monoen (0.1 0.5 1 2 3 4 6 9 10) #keV
#foreach monoen (0.1) #keV
echo "\nSimulating monochromatic energy: $monoen keV for $simTime seconds";
      
    set root=mono${monoen}_${cps}cps_pix${pixel}_${simTime}s
    #pixdetillum PixImpList=${root}.piximpact XMLFile=$XMLrect tstart=0. Tstop=$simTime pixels=$pixel rate=$cps energy=$monoen clobber=yes seed=-1
      
    tessim PixelID=$pixel PixImpList=$PIXIMPACTdir/${root}.piximpact Streamfile=${root}.stream tstart=0. tstop=$simTime sample_rate=$samprate clobber=yes
      
    streamtotriggers PixImpList=${root}.piximpact XMLFile=$XMLrect tstart=0. tstop=$simTime Streamfile=${root}.stream TesTriggerFile=${root}.fits TriggerSize=10000 PreBufferSize=1000 pixels=$pixel clobber=yes  
    rm ${root}.stream
		  
end # energies for calib

#
# Calculate CLAIBRATION factors (all possible combinations)
#
foreach monoen1 (0.1 0.5 1 2 3 4 6 9 10) #keV
    foreach monoen2 (0.1 0.5 1 2 3 4 6 9 10) #keV
	  set root1=mono${monoen1}_${cps}cps_pix${pixel}_${simTime}s
	  set root2=mono${monoen2}_${cps}cps_pix${pixel}_${simTime}s
	  rm evtcal.fits
	  tesreconstruction Recordfile=${root1}.fits TesEventFile="evtcal1.fits" Rcmethod="SIRENA" PulseLength=$pulseLength LibraryFile=$libFile scaleFactor=$scaleFactor samplesUp=$samplesUp nSgms=$nSgms crtLib=1 mode=0 clobber=yes intermediate=0 monoenergy=${monoen}E3 baseline=$baseline EventListSize=1000 filterMethod="B0"


cd $pwd

  
  
  
