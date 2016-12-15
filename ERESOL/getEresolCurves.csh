#!/usr/bin/csh -fe

#
# Get energy resolution curves as a function of pulses pair separation

#
# IT USES INTEGRATION OF SIRENA INTO SIXTE, so SIXTE must be initialized
#
# RUN as:
#        source getEresolCurves.csh monoenergy tessim lib samplesUp nSgms
#
unset monoEkeV monoEeV calibDir pulseLength baseline tessim samplesUp nSgms scaleFactor libDir simDir noiseDir noiseFile filterMeth array energyMethod eresolFile evtFile root rootEvt sepsStr
alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'  #!!!! BE CAREFUL - bc does not support Scientific Notation

# Read input options
#---------------------
echo "Run as:"
echo "source getEresolCurves.csh  monoenergy array lib samplesUp nSgms filterMeth nsamples fdomain energyMethod OFLib ACDC "
echo ""
echo "       monoenergy: monochromatic energy (keV) of the pulses"
echo "       array: array configuration (LPA1,LPA2,LPA3,SPA)"
echo "       lib: 'monolib' or 'multilib' or 'multilibOF'"
echo "       samplesUp: number of samples above threshold"
echo "       nSgms: number of sigmas to define threshold"
echo "       filterMeth: filter method F0 (freq. bin=0) or B0 (subtract Baseline)"
echo "       nsamples: pulse length & noise samples"
echo "       fdomain: filter domain (freq:F, time:T)"
echo "       energyMethod: Energy reconstruction method (OPTFILT, WEIGHT, WEIGHTN, I2R)"
echo "       ACDC: AC or DC (if empty -> old tessim -> DC"
echo ""   

if($#argv < 9) then 
    echo "Please complete input command line with all the parameters"
    echo "   source getEresolCurves.csh monoenergy array lib samplesUp nSgms F0/B0 nsamples fdomain energyMethod ACDC"
    exit
endif    

# Read input parameters
#######################
set monoEkeV=$1 # file monochromatic energy (keV)
MATH monoEeV = ($monoEkeV*1000.)
set array=$2 # LPA1, LPA2, LPA3, SPA
set tessim="tessim${array}"   # ex. tessimLPA1
set labelLib=$3 # monolib or multilib or nonoiselib
set OFLib="no"
if($labelLib == "multilibOF") set OFLib="yes"

set samplesUp=$4  # 2
set nSgms=$5      # 10
set filterMeth=$6 # F0/B0
set nsamples=$7    # 1024
set fdomain=$8     # F(requency) or T(ime)
set energyMethod=$9 # LAGS, OPTFILT, WEIGHT, WEIGHTN, I2R
set ACDC=""
if ($10 != "") set ACDC=`echo $10 | tr '[:upper:]' '[:lower:]'` # AC or DC

set PreBufferSize=1000
#MATH tstartPulse1 = ($PreBufferSize + 1)
MATH tstartPulse1 = ($PreBufferSize - 1)
set Fil=""
if($tstartPulse1 > 0) then
  set nSgms=0
  set samplesUp=0
  set scaleFactor=0
endif

set pwd=`pwd`

#set exptime="100" #s 
set nSimPulses=5000 # minimum number of simulated pulses in simulPairs.csh (approx)
set simDir="${pwd}/PAIRS/tessim${array}"
set resultsDir="${pwd}/PAIRS/eresol${array}"
set PixType="file:/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/newpixels.fits[${array}${ACDC}]"

# Define tmp dir
set tmpdir="eresolPairs${array}_${monoEkeV}_${energyMethod}_${labelLib}"
mkdir -p /tmp/${tmpdir}/pfiles
setenv PFILES "/tmp/${tmpdir}/pfiles;${HEADAS}/syspfiles:${SIXTE}/share/simput/pfiles:$SIXTE/share/sixte/pfiles"
echo "PFILES=$PFILES"
setenv HEADASNOQUERY
setenv HEADASPROMPT /dev/null


# SIRENA definitions
#######################

set scaleFactor=0.005

# LIBDIR && LIBFILE && NOISEDIR && NOISEFILE
#############################################

set noiseDir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/NOISE/${tessim}"
set noiseFile=${noiseDir}/noise${nsamples}samples_${tessim}_B0_100s_pairscps_ADC.fits
if($energyMethod == "I2R") set noiseFile=${noiseDir}/noise${nsamples}samples_${tessim}_B0_100s_pairscps_R.fits  

set libDirRoot="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/${tessim}"
#set libDirRoot="$HOME/INSTRUMEN/EURECA/ERESOL/PAIRS/eresolSPA/OKforToulouse/LIBRARIES"
set libDir=${libDirRoot}/${energyMethod}
if($OFLib == "yes") set libDir=${libDir}LIB

unset libFile
# OPTFILT (mono & multi & fixed lib) WEIGHT (multilib) I2R (multi, mono & fixed lib)
if($labelLib =~ multi*) set libFile=${libDir}/libraryMultiE_GLOBAL_PL${nsamples}_${tessim}${Fil}.fits
if($labelLib == "monolib" ) set libFile=${libDir}/libraryMultiE_PL${nsamples}_${monoEkeV}keV_${tessim}${Fil}.fits
if($labelLib == "fixedlib" ) set libFile=${libDir}/libraryMultiE_PL${nsamples}_1keV_${tessim}${Fil}.fits
if($labelLib == "nonoiselib") set libFile=${libDir}/libraryMultiE_PL${nsamples}_${monoEkeV}keV_${tessim}Nonoise.fits

#
# IN case OLD tessim (no trigger must be used)
#
if ($ACDC == "") then
    set nSimPulses=3000 # minimum number of simulated pulses in simulPairs.csh (approx) -> old_tessim
    set PixType="file:/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/tespixels.fits[${array}]"
    set noiseFile=${noiseDir}/noise${nsamples}samples_${tessim}_B0_100s_pairscps.fits
endif


# get baseline from noise file (in current)
fkeypar fitsfile=${noiseFile}+1 keyword="BASELINE"
set baseline=`pget fkeypar value`

cd $resultsDir

if($array == "SPA") then
    
    set sepsStr=(00004 00006 00009 00014 00021 00031 00047 00071 00108 00163 00246 00371 00560 00845 01276 01926 02907 04389)
    set sepsStr=(10000)
    
else if($array == "LPA1") then
    
    set sepsStr=(00004 00006 00009 00014 00021 00031 00047 00071 00108 00163 00246 00371 00560 00845 01276 01926 02907 04389 06625 10000)
        
else if(${array} == "LPA2") then
    
    set sepsStr=(00004 00006 00009 00014 00021 00031 00047 00071 00108 00163 00246 00371 00560 00845 01276 01926 02907 04389 06625 10000)
    
else if(${array} == "LPA3") then
    set scaleFactor=0.02  
    # if triggering is required, do filtering
    if($tstartPulse1 == 0) set Fil="Fil"
    set sepsStr=(00004 00006 00009 00014 00021 00031 00047 00071 00108 00163 00246 00371 00560 00845 01276 01926 02907 04389 06625 10000)
    #set sepsStr=(20000)
endif


#
# Create output ERESOL file
#---------------------------
set root="${nSimPulses}p_SIRENA${nsamples}_${monoEkeV}keV_${filterMeth}${fdomain}_${labelLib}_${energyMethod}"
#set root="100s_SIRENA${nsamples}_${monoEkeV}keV_${filterMeth}${fdomain}_${labelLib}_${energyMethod}"
if($tstartPulse1 > 0) set root=${root}_NTRIG
set eresolFile="eresol_${root}.dat"

# write header in eresol file
echo "samSep  FWHMP(E-Eb)  FWHMS(E-Eb)  EBIASP      EBIASS     FWHMP(E)  FWHMS(E)   FWHMP_ERR   FWHMS_ERR">$eresolFile

#
# Process input data files
#--------------------------
foreach sep12 ($sepsStr)
	MATH sep = ($sep12 * 1.)
	set inFile="${simDir}/sep${sep12}sam_${nSimPulses}p_${monoEkeV}keV.fits"
	#set inFile="${simDir}/sep${sep12}sam_100s_${monoEkeV}keV.fits"
	set evtFile="events_sep${sep12}sam_${root}.fits"
	echo "============================================="
	echo "Using file: $inFile"
	echo "Using library: $libFile"
	echo "Using noisefile: $noiseFile"
	echo "Using Baseline: $baseline"
	echo "Setting evtFile: $evtFile"
	echo "Setting eresolFile: $eresolFile"
	echo "============================================="
        
        # when run also in detection mode, run it iteratively to get the best possible combination
        # of sigmas/samples able to detect all pulses
        set sigmasMin=0
	if($monoEkeV == 1)then 
	    set sigmasMax=50
	else
	    set sigmasMax=100
	endif
	set iter=0
	
	#MATH tstartPulse1 = ($PreBufferSize + 1)
	MATH tstartPulse1 = ($PreBufferSize - 1) #999
	#MATH tstartPulse2 = ($tstartPulse1 + $sep)
	MATH tstartPulse2 = ($tstartPulse1 + 1 + $sep) #1000 + sep
	set tstartPulse3=0
	
	if($tstartPulse1 > 0) then
	  set nSgms=0
	  set samplesUp=0
	  set scaleFactor=0
	endif
            
	# SIRENA processing
	#-------------------     
	if (-e $evtFile) then
	
	    # if events already exist get number nSgms used to get events
	    echo "Event file $evtFile already DOES exist"
	    set nSgmsStr=`fkeyprint $evtFile+0 HISTORY|grep P19`
	    set nSgms=`expr "$nSgmsStr" : '.* \(.*\)'`
	    
	else
	    if ($ACDC == "") then # OLD tessim
		set command="tesreconstruction Recordfile=$inFile TesEventFile=$evtFile Rcmethod="SIRENA" PulseLength=$nsamples LibraryFile=$libFile scaleFactor=$scaleFactor samplesUp=$samplesUp nSgms=$nSgms crtLib=0 mode=1 NoiseFile=$noiseFile FilterDomain=$fdomain FilterMethod=$filterMeth calibLQ=1 b_cF=$b_cF c_cF=$c_cF clobber=yes intermediate=0 baseline=$baseline  EnergyMethod=$energyMethod PixelType=$array tstartPulse1=$tstartPulse1 tstartPulse2=$tstartPulse2 tstartPulse3=$tstartPulse3 OFInterp=DAB"
	    else
	    
		set command="tesreconstruction Recordfile=$inFile TesEventFile=$evtFile Rcmethod="SIRENA" PulseLength=$nsamples LibraryFile=$libFile scaleFactor=$scaleFactor samplesUp=$samplesUp nSgms=$nSgms OFLib=$OFLib mode=1 NoiseFile=$noiseFile FilterDomain=$fdomain FilterMethod=$filterMeth clobber=yes intermediate=0 baseline=$baseline  EnergyMethod=$energyMethod PixelType=$array tstartPulse1=$tstartPulse1 tstartPulse2=$tstartPulse2 tstartPulse3=$tstartPulse3 OFInterp=DAB"
	    endif
	    echo "Running SIRENA as: $command"
	    $command
	    
	    if ($status) then
		  echo "Error running tesreconstruction"
		  cd $pwd
		  exit
	    endif
	    echo "SIRENA proc with $filterMeth$fdomain for $array with $energyMethod method finished OK"
	    #cd $pwd; exit
		      
	    # check that detection was ok
	    fkeypar fitsfile=$inFile keyword="NETTOT"
	    set nsimpulses=`pget fkeypar value`
	    fkeypar fitsfile=$evtFile keyword="NAXIS2"
	    set nrows=`pget fkeypar value`
	    set ndetpulses=0
	    if ($nrows > 0) then
		fstatistic infile=$evtFile colname="GRADE1" rows="-" minval=0 >/dev/null
		set ndetpulses=`pget fstatistic numb`
	    endif
	      
	    while ($nsimpulses != $ndetpulses)
		if($iter>50) break
		MATH iter = ($iter+1)
		echo "With nSgms=$nSgms => sim/det: $nsimpulses/$ndetpulses" 
		if($nsimpulses<$ndetpulses) then
		    #increase sigmas
		    set sigmasMin=$nSgms
		else if($nsimpulses>$ndetpulses) then
		    #decrease sigmas
		    set sigmasMax=$nSgms
		endif
		MATH nSgms = ($sigmasMin + ($sigmasMax-$sigmasMin)/2.)
		echo "   Trying with nSgms=$nSgms (sigmasMin/Max:$sigmasMin/$sigmasMax)"
		
		if ($ACDC == "") then # OLD tessim
		    set command="tesreconstruction Recordfile=$inFile TesEventFile=$evtFile Rcmethod="SIRENA" PulseLength=$nsamples LibraryFile=$libFile scaleFactor=$scaleFactor samplesUp=$samplesUp nSgms=$nSgms crtLib=0 mode=1 NoiseFile=$noiseFile FilterDomain=$fdomain FilterMethod=$filterMeth calibLQ=1 b_cF=$b_cF c_cF=$c_cF clobber=yes intermediate=0 baseline=$baseline  EnergyMethod=$energyMethod PixelType=$array tstartPulse1=$tstartPulse1 tstartPulse2=$tstartPulse2 tstartPulse3=$tstartPulse3 OFInterp=DAB"
		else
	    
		    set command="tesreconstruction Recordfile=$inFile TesEventFile=$evtFile Rcmethod="SIRENA" PulseLength=$nsamples LibraryFile=$libFile scaleFactor=$scaleFactor samplesUp=$samplesUp nSgms=$nSgms OFLib=$OFLib mode=1 NoiseFile=$noiseFile FilterDomain=$fdomain FilterMethod=$filterMeth clobber=yes intermediate=0 baseline=$baseline  EnergyMethod=$energyMethod PixelType=$array tstartPulse1=$tstartPulse1 tstartPulse2=$tstartPulse2 tstartPulse3=$tstartPulse3 OFInterp=DAB"
		endif
		echo "Running SIRENA as: $command"
		$command
		
		fkeypar fitsfile=$evtFile keyword="NAXIS2"
		set nrows=`pget fkeypar value`
		set ndetpulses=0
		if ($nrows > 0) then
		    fstatistic infile=$evtFile colname="GRADE1" rows="-" minval=0 >/dev/null
		    set ndetpulses=`pget fstatistic numb`      
		endif
		  
	    end #while
	    if($iter>50)then
		rm $evtFile
		set iter=0
		echo "Reconstruction NOT finished with nSgms=$nSgms"
		continue
	    endif
	      
	endif #if evtFile exists...
	set nSgms=`printf "%.1f" $nSgms`
	echo "Reconstruction finished with nSgms=$nSgms"
        
            
	#-----------------------------------------
	# EVENT file processing to calculate FWHM
	#-----------------------------------------
        #echo "EVT file: $evtFile"
	set rootEvt=$evtFile:r
	set evt1="${rootEvt}_primaries.fits"
	set evt2="${rootEvt}_secondaries.fits"
	if(-e $evt1) rm $evt1
	if(-e $evt2) rm $evt2
	
	#use only rows with GRADE1>0 (non initially truncated)
 	fselect infile=$evtFile outfile=\!$evt1 expr="#ROW%2!=0 && GRADE1>0"
 	fkeypar fitsfile=$evt1 keyword="NAXIS2"
 	set nrows=`pget fkeypar value`

	# PRIMARIES BIAS and FWHM 
	#-------------------------
	# create column Ei-Eo (in eV)
	fcalc infile=$evt1 outfile=\!$evt1 clname=EeVBIAS_i expr="(SIGNAL*1000.-$monoEeV)"
	fstatistic infile=$evt1 colname=EeVBIAS_i rows="-" >/dev/null 
	set ebias=`pget fstatistic mean`
	set ebias1=`printf "%.3f" $ebias`
 	# create column with BIAS corrected energies (in eV)
 	fcalc infile=$evt1 outfile=\!$evt1 clname=EeVbiasCorr expr="(SIGNAL*1000.-$ebias1)"
# 	# create column SQR(Ei(biasCorr)-Eo)
# 	fcalc infile=$evt1 outfile=\!$evt1 clname=EbcDIFFSQR expr="(EeVbiasCorr-$monoEeV)^2"
# 	# calculate FWHM for primaries (BIAS corrected)
# 	fstatistic infile=$evt1 colname=EbcDIFFSQR rows="-" >/dev/null
# 	set sigma1sqr=`pget fstatistic mean`
# 	MATH fwhm = (sqrt($sigma1sqr) * 2.35)   # FWHM for PRIMARIES in eV
# 	set fwhm1=`printf "%.3f" $fwhm`
	fstatistic infile=$evt1 colname=EeVbiasCorr rows="-" >/dev/null
 	set sigma=`pget fstatistic sigma`
	MATH fwhm = ($sigma * 2.35)   # FWHM for PRIMARIES in eV
	set fwhm1=`printf "%.3f" $fwhm`
	
	# calculate FWHM for primaries (NO BIAS corrected) : NEW METHOD AFTER GAIN-CORRECTION
# 	#fcalc infile=$evt1 outfile=\!$evt1 clname=EDIFFSQR expr="(SIGNAL*1000.-$monoEeV)^2"
# 	fstatistic infile=$evt1 colname=EDIFFSQR rows="-" >/dev/null
# 	set sigma1sqr=`pget fstatistic mean`
# 	set sum1sqr=`pget fstatistic sum`
# 	MATH fwhm_nocorr = (sqrt($sigma1sqr) * 2.35)   # FWHM for PRIMARIES in eV
# 	set fwhm1_nocorr=`printf "%.3f" $fwhm_nocorr`
	fstatistic infile=$evt1 colname=SIGNAL rows="-" >/dev/null
	set sigma=`pget fstatistic sigma`
	MATH fwhm_nocorr = ($sigma * 2.35 *1000.)   # FWHM for PRIMARIES in eV
	set fwhm1_nocorr=`printf "%.3f" $fwhm_nocorr`
	
	#calculate error for fwhm1
# 	fstatistic infile=$evt1 colname=SIGNAL rows="-" >/dev/null
# 	set sigma_keV=`pget fstatistic sigma`
# 	#MATH fwhm_err = (2.35 * $sigma_keV * 1000. / sqrt($nrows))
# 	MATH fwhm_err = (2.35 * $sigma_keV * 1000. / sqrt(2*$nrows-2))
	MATH fwhm_err = ($fwhm_nocorr / sqrt(2*$nrows-2))
	set fwhm1_err=`printf "%.3f" $fwhm_err`
	
	
	# SECONDARIES BIAS and FWHM 
	#-------------------------
	fselect infile=$evtFile outfile=\!$evt2 expr="#ROW%2==0 && GRADE1>0"
 	fkeypar fitsfile=$evt2 keyword="NAXIS2"
 	set nrows=`pget fkeypar value`
 	
	# create column Ei-Eo (in eV)
	fcalc infile=$evt2 outfile=\!$evt2 clname=EeVBIAS_i expr="(SIGNAL*1000.-$monoEeV)"
	fstatistic infile=$evt2 colname=EeVBIAS_i rows="-" >/dev/null 
	set ebias=`pget fstatistic mean`
	set ebias2=`printf "%.3f" $ebias`
	# create column with BIAS corrected energies (in eV)
	fcalc infile=$evt2 outfile=\!$evt2 clname=EeVbiasCorr expr="(SIGNAL*1000.-$ebias2)"
# 	# create column SQR(Ei(biasCorr)-Eo)
# 	fcalc infile=$evt2 outfile=\!$evt2 clname=EbcDIFFSQR expr="(EeVbiasCorr-$monoEeV)^2"
# 	# calculate FWHM for secondaries (BIAS corrected)
# 	fstatistic infile=$evt2 colname=EbcDIFFSQR rows="-" >/dev/null
# 	set sigma2sqr=`pget fstatistic mean`
# 	MATH fwhm = (sqrt($sigma2sqr) * 2.35)   # FWHM for SECONDARIES in eV
# 	set fwhm2=`printf "%.3f" $fwhm`
	fstatistic infile=$evt2 colname=EeVbiasCorr rows="-" >/dev/null
 	set sigma=`pget fstatistic sigma`
	MATH fwhm = ($sigma * 2.35)   # FWHM for PRIMARIES in eV
	set fwhm2=`printf "%.3f" $fwhm`

	# calculate FWHM for secondaries (NO BIAS corrected)
# 	fcalc infile=$evt2 outfile=\!$evt2 clname=EDIFFSQR expr="(SIGNAL*1000.-$monoEeV)^2"
# 	fstatistic infile=$evt2 colname=EDIFFSQR rows="-" >/dev/null
# 	set sigma2sqr=`pget fstatistic mean`
# 	MATH fwhm_nocorr = (sqrt($sigma2sqr) * 2.35)   # FWHM for SECONDARIES in eV
# 	set fwhm2_nocorr=`printf "%.3f" $fwhm_nocorr`
	fstatistic infile=$evt2 colname=SIGNAL rows="-" >/dev/null
	set sigma=`pget fstatistic sigma`
	MATH fwhm_nocorr = ($sigma * 2.35 *1000.)   # FWHM for PRIMARIES in eV
	set fwhm2_nocorr=`printf "%.3f" $fwhm_nocorr`
	
	#calculate error for fwhm2
# 	fstatistic infile=$evt2 colname=SIGNAL rows="-" >/dev/null
# 	set sigma_keV=`pget fstatistic sigma`
# 	#MATH fwhm_err = (2.35 * $sigma_keV * 1000. / sqrt($nrows))
# 	MATH fwhm_err = (2.35 * $sigma_keV * 1000. / sqrt(2*$nrows-2))
	MATH fwhm_err = ($fwhm_nocorr / sqrt(2*$nrows-2))
	set fwhm2_err=`printf "%.3f" $fwhm_err`

	echo "$sep12   $fwhm1        $fwhm2         $ebias1      $ebias2      $fwhm1_nocorr      $fwhm2_nocorr      $fwhm1_err        $fwhm2_err" >>$eresolFile
      
end # infiles sep

cd $pwd
rm -rf /tmp/${tmpdir}

