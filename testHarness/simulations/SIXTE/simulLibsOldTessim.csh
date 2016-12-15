#!/bin/csh
#
# Create libraries and files for libraries
# TO BE USED ONLY WITH tessim-pretrigger
#
#  Input parameters: 
#          array (see possibilities/examples below)
#          pulseLength 
#          reconMethod  (OPTFILT, WEIGHT, I2R)
#          [OFLib]  (lib inlcude optimal filters?: Y/N, def=N)
#          [nonoise]
#
# Unset variables to be used in script
#
unset calibDir XMLrect pixel samprate pulseLength baseline tessim cps simTime samplesUp nSgms scaleFactor XMLFile libFile m_unknown simnoise PixType PIXIMPACTdir
unset monoen root nsimpulses ndetpulses array PreBufferSize triggerSize 
alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'

#---------------
# PARAMETERS
#---------------
if($#argv < 3) then 
    echo "Please complete input command line with all the parameters"
    echo "   source simulLibs.csh array pulseLength reconMethod [OPTF] [nonoise]"
    exit
endif 
set array=$1
set tessim="tessim$array"
set pulseLength=$2 # also noise samples
MATH separation = (10*$pulseLength)
set PreBufferSize=1000
MATH tstartPulse1 = ($PreBufferSize + 1)
MATH tstartPulse2 = ($tstartPulse1 + $separation)
set tstartPulse3=0 # simulated files for calibration are pairs of pulses
MATH triggerSize = ($PreBufferSize+$separation+$separation+1000)
set reconMethod=$3
set OFLib="no"
if($#argv > 3) set OFLib=$4
set noise=""
set Fil=""
set simnoise="simnoise=y"
set simTime=1500 #120
set pixel=1
set samprate=156250 #Hz
set libEnergies=(0.1 0.5 1 2 3 4 5 6 7 8 9 10)
set lastEnergy=$libEnergies[$#libEnergies]

# if($#argv == 4) then
#     set noise="NONOISE"
#     set simnoise="simnoise=n"
#     set simTime=10
# endif

set tmpdir="simulLibs${array}_${reconMethod}"
mkdir -p /tmp/${tmpdir}/pfiles
setenv PFILES "/tmp/${tmpdir}/pfiles;${HEADAS}/syspfiles:${SIXTE}/share/simput/pfiles:$SIXTE/share/sixte/pfiles"
echo "PFILES=$PFILES"
setenv HEADASNOQUERY
setenv HEADASPROMPT /dev/null


if($array == "SPA" || $array == "LPA1") then
    set samplesUp=3
    set nSgms=20
    set scaleFactor=0.005
else if($array == "LPA2") then
    set samplesUp=2
    set scaleFactor=0.005
else if($array == "LPA3") then
      set samplesUp=2
      set scaleFactor=0.02
      set Fil="Fil"
endif
    
if($tstartPulse1 > 0) then
      set nSgms=0
      set samplesUp=0
      set scaleFactor=0
      set Fil=""
endif

# CALIB STUFF
set calibDir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/CALIBRATION/pulsetempl_v16092014/156khz"
set XMLrect="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/CALIBRATION/newFormat/xifu_detector_300us_156kHz_rect.xml";
set XMLfile=$XMLrect
set PIXIMPACTdir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/PIXIMPACT"
set SIMFILESdir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/${tessim}"
set PixType="file:/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/tespixels.fits[${array}${noise}]"


# NOISE FILES
set noiseDir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/NOISE/${tessim}"
set noiseFile=${noiseDir}/noise${pulseLength}samples_${tessim}_B0_100s_pairscps.fits
if($reconMethod == "I2R") then
#   (heainit;fkeypar fitsfile=${noiseFile}+1 keyword="BASELINR")
   set noiseFile=${noiseDir}/noise${pulseLength}samples_${tessim}_B0_100s_pairscps_R.fits
endif

# always use BASELINE=BASELINE(I) because input data is always in current
(heainit;fkeypar fitsfile=${noiseFile}+1 keyword="BASELINE")
set baseline=`heainit;pget fkeypar value`
#echo "baseline=$baseline and $simnoise"

#
# Create LIBRARIES and FILES for libraries
# 
set pwd=`pwd`
set libDir=${tessim}/${reconMethod}
if($OFLib =~ y* || $OFLib =~ Y*) set libDir=${tessim}/${reconMethod}LIB
cd LIBRARIES/$libDir

set libFile = libraryMultiE_PL${pulseLength}_${tessim}${noise}${Fil}.fits
if (-e $libFile) rm $libFile

foreach monoen ($libEnergies) #keV

    set lastElibrary=0
    if($monoen == $lastEnergy) set lastElibrary=1
    echo "\nSimulating monochromatic energy: $monoen keV for $PixType";
    
    if($array == "LPA2" && $monoen == 0.5) then
      set nSgms=11.733
    else if($array == "LPA2") then
      set nSgms=20
    endif
    
    if($array == "LPA3" && $monoen == 0.5) then
      set nSgms=9.8788
    else if($array == "LPA3" && $monoen == 1) then
      set nSgms=18.69
    else if($array == "LPA3" && $monoen == 2)then
      set nSgms=37.03
    else if($array == "LPA3" && $monoen == 3)then
      set nSgms=45
    endif
    
    if($tstartPulse1 > 0) then
	set nSgms=0
	set samplesUp=0
	set scaleFactor=0
	set Fil=""
    endif
    
    #
    # Simulate SIXTE file
    #
   
    set root0=mono${monoen}_sep${separation}_pix${pixel}_${simTime}s
    set root=mono${monoen}_sep${separation}_pix${pixel}_${simTime}s_${pulseLength}$noise
    set pixFile=${PIXIMPACTdir}/${root0}_trSz${triggerSize}.piximpact
    set streamFile=${SIMFILESdir}/${root0}.stream
    set simFile=${SIMFILESdir}/${root}.fits
    
    if(! -e $pixFile) then
	tesconstpileup PixImpList=$pixFile XMLFile=$XMLrect tstop=$simTime energy=$monoen pulseDistance=$separation TriggerSize=$triggerSize clobber=yes
    endif	
    
    if(! -e ${simFile}) then
	tessim PixID=$pixel PixImpList=$pixFile Streamfile=$streamFile tstart=0. tstop=$simTime sample_rate=$samprate clobber=yes $simnoise PixType="$PixType" 
     
	streamtotriggers PixImpList=$pixFile XMLFile=$XMLrect tstart=0. tstop=$simTime Streamfile=$streamFile TesTriggerFile=$simFile TriggerSize=$triggerSize PreBufferSize=$PreBufferSize pixels=$pixel clobber=yes 
	(heainit;cphead ${streamFile}+1 $simFile+1)
	rm $streamFile
	
	#rm first (and LAST) record and update NETTOT (first starts high and last can be cut)
	fkeypar fitsfile=$simFile key="NAXIS2"
 	set nrows=`pget fkeypar value`
 	fdelrow infile=${simFile}+1 firstrow=$nrows nrows=1 confirm="no" proceed=yes
 	fdelrow infile=${simFile}+1 firstrow=1 nrows=1 confirm="no" proceed=yes
 	fkeypar fitsfile=${simFile} key="NAXIS2"
 	set nrows=`pget fkeypar value`
 	MATH nettot = ($nrows*2)
 	fparkey fitsfile=${simFile} value=$nettot keyword="NETTOT"
    endif	  
    #continue
    #
    # Create/update LIBRARY
    #
    if(-s "evtcal.fits") rm evtcal.fits
    
#     set command="tesreconstruction Recordfile=${simFile} TesEventFile="evtcal.fits" Rcmethod="SIRENA" PulseLength=$pulseLength LibraryFile=$libFile scaleFactor=$scaleFactor samplesUp=$samplesUp nSgms=$nSgms crtLib=1 mode=0 clobber=yes intermediate=0 monoenergy=${monoen}E3 baseline=$baseline EventListSize=1000 filterMethod="F0" EnergyMethod=$reconMethod lastElibrary=$lastElibrary tstartPulse1=$tstartPulse1 tstartPulse2=$tstartPulse2 tstartPulse3=$tstartPulse3"
    set command="tesreconstruction Recordfile=${simFile} TesEventFile="evtcal.fits" Rcmethod="SIRENA" PulseLength=$pulseLength LibraryFile=$libFile scaleFactor=$scaleFactor samplesUp=$samplesUp nSgms=$nSgms mode=0 OFLib=$OFLib clobber=yes intermediate=0 monoenergy=${monoen}E3 baseline=$baseline EventListSize=1000 filterMethod="F0" EnergyMethod=$reconMethod lastElibrary=$lastElibrary tstartPulse1=$tstartPulse1 tstartPulse2=$tstartPulse2 tstartPulse3=$tstartPulse3 NoiseFile=$noiseFile"

    echo "Running library creation in LIBRARIES/${libDir} as:\n$command"
    $command
    
    (heainit;fkeypar fitsfile=${simFile} keyword="NETTOT")
    unset nsimpulses
    set nsimpulses=`heainit;pget fkeypar value`
    if(-s "evtcal.fits") then
	(heainit;fkeypar fitsfile=evtcal.fits keyword="NAXIS2")
	unset ndetpulses
	set ndetpulses=`heainit;pget fkeypar value`
    else
	echo "ERROR: event file has not been created"
        set ndetpulses="XXXX"
    endif
    echo "For sim file ${simFile} : simpulses=$nsimpulses, detpulses=$ndetpulses"
    rm "evtcal.fits"
    if($ndetpulses > $nsimpulses) then
	 cd $pwd
	 echo "Library creation stopped: invented pulses"
	 exit
    endif
    if($ndetpulses == 0) then
	 cd $pwd
	 exit
    endif
end # energies for lib

cd $pwd
rm -rf /tmp/$tmpdir

  
  
  
