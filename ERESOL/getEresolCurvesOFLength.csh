#!/usr/bin/csh

#
# Get energy resolution curves as a function of pulses pair separation

#
# IT USES INTEGRATION OF SIRENA INTO SIXTE, so SIXTE must be initialized
#
# RUN as:
#        source getEresolCurves.csh monoenergy tessim lib samplesUp nSgms
#
unset monoEkeV monoEeV calibDir pulseLength baseline tessim samplesUp nSgms scaleFactor libDir simDir noiseDir noiseFile filterMeth array energyMethod eresolFile evtFile root rootEvt sepsStr
alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'
set dateNow=`date +%s`
mkdir -p /tmp/${dateNow}/pfiles
setenv PFILES "/tmp/${dateNow}/pfiles;${HEADAS}/syspfiles:${SIXTE}/share/simput/pfiles:$SIXTE/share/sixte/pfiles"
echo "PFILES=$PFILES"
setenv HEADASNOQUERY
setenv HEADASPROMPT /dev/null


# Read input options
#---------------------
echo "Run as:"
echo "source getEresolCurvesOFLength.csh  monoenergy array lib samplesUp nSgms filterMeth nsamples fdomain energyMethod OFLength"
echo ""
echo "       monoenergy: monochromatic energy (keV) of the pulses"
echo "       array: array configuration (LPA1,LPA2,LPA3,SPA)"
echo "       lib: 'monolib' o 'multilib'"
echo "       samplesUp: number of samples above threshold"
echo "       nSgms: number of sigmas to define threshold"
echo "       filterMeth: filter method F0 (freq. bin=0) or B0 (subtract Baseline)"
echo "       nsamples: pulse length & noise samples"
echo "       fdomain: filter domain (freq:F, time:T)"
echo "       energyMethod: Energy reconstruction method (LAGS, NOLAGS, WEIGHT)"
echo "       OFLength: optimal filter length for OFStrategy=FIXED"
echo ""   

if($#argv < 9) then 
    echo "Please complete input command line with all the parameters"
    echo "   source getEresolCurvesOFLength.csh monoenergy array lib samplesUp nSgms F0/B0 nsamples fdomain energyMethod OFLength"
    exit
endif    

# Read input parameters
#######################
set monoEkeV=$1 # file monochromatic energy (keV)
MATH monoEeV = ($monoEkeV*1000.)
set array=$2 # LPA1, LPA2, LPA3, SPA
set tessim="tessim${array}"   # ex. tessimLPA1
set labelLib=$3 # monolib or multilib or nonoiselib
set samplesUp=$4  # 2
set nSgms=$5      # 10
set filterMeth=$6 # F0/B0
set nsamples=$7    # 1024
set fdomain=$8     # F(requency) or T(ime)
set energyMethod=$9 # LAGS, NOLAGS, WEIGHT
set oflength=$10

# SIRENA definitions
#######################

set scaleFactor=0.005
#set b_cF=1.6823832496783
#set c_cF=-0.216608489549083
set b_cF=1.
set c_cF=0.

#echo "Using baseline=$baseline"

set libDir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/${tessim}"
set exptime="100" #s 
set simDir="./PAIRS/${tessim}"
set noiseDir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/NOISE/${tessim}"
set noiseFile=${noiseDir}/noise${nsamples}samples_${tessim}_B0_100s_pairscps.fits
fkeypar fitsfile=${noiseFile}+1 keyword="BASELINE"
set baseline=`pget fkeypar value`
#set noiseDir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/NOISE/tessim20150505PP"
#set noiseFile=${noiseDir}/noise${nsamples}samples_tessim20150505PP_B0.fits
#if($tessim =~ *Nonoise) then
#      set noiseDir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/NOISE/tessim20150505"
#      set noiseFile="${noiseDir}/noise${nsamples}samples_tessim20150505_B0.fits"
#endif
set PixType="file:/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/tespixels.fits[${array}]"
#
unset libFile
if($labelLib == "multilib") set libFile="${libDir}/libraryMultiE_PL${nsamples}_${tessim}.fits"
if($labelLib == "monolib" ) then
    #
    if($monoEkeV == 1) set libFile=${libDir}/libraryMultiE_PL${nsamples}_1keV_${tessim}.fits
    #if($monoEkeV == 1) set libFile=${libDir}/libraryMonoE_PL${nsamples}_1keV_${tessim}.fits
    if($monoEkeV == 3) set libFile=${libDir}/libraryMultiE_PL${nsamples}_3keV_${tessim}.fits
    if($monoEkeV == 6) set libFile=${libDir}/libraryMultiE_PL${nsamples}_6keV_${tessim}.fits
    #
endif
if($labelLib == "nonoiselib") then
    if($monoEkeV == 1) set libFile=${libDir}/libraryMultiE_PL${nsamples}_1keV_${tessim}Nonoise.fits
    if($monoEkeV == 3) set libFile=${libDir}/libraryMultiE_PL${nsamples}_3keV_${tessim}Nonoise.fits
    if($monoEkeV == 6) set libFile=${libDir}/libraryMultiE_PL${nsamples}_6keV_${tessim}Nonoise.fits
    #
endif

set pwd=`pwd`
cd $simDir

if($array == "LPA1") then
    #set sepsStr=(0004 0008 0016 0032 0064 0128 0256 0512 1024 2048 4096 8192)
    set sepsStr=(00004 00005 00007 00010 00013 00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287)
    set sepsStr=(0010 3000)
else if(${array} == "LPA2") then
    set sepsStr=(00004 00005 00007 00010 00013 00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287 11115 14910 20000)
else if(${array} == "LPA3") then
    set scaleFactor=0.02    
    if($labelLib == "monolib") then
	if($monoEkeV == 1) set libFile=${libDir}/libraryMultiE_PL${nsamples}_1keV_${tessim}Fil.fits
	if($monoEkeV == 6) set libFile=${libDir}/libraryMultiE_PL${nsamples}_6keV_${tessim}Fil.fits
    endif
    if($labelLib == "multilib") set libFile="${libDir}/libraryMultiE_PL${nsamples}_${tessim}Fil.fits"
    set sepsStr=(00004 00005 00007 00010 00013 00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287 11115 14910 20000)
    #set sepsStr=(20000)
endif

#
# Create output ERESOL file
#---------------------------
set root="${exptime}s_SIRENA${nsamples}_${monoEkeV}keV_${filterMeth}${fdomain}_${labelLib}_${energyMethod}_OFL${oflength}"
set eresolFile="eresol_${root}.dat"

# write header in eresol file
echo "samSep  FWHM(prim)  FWHM(secon)  EBIAS(prim)    EBIAS(sec)   FWHM(prim,nocorr)  FWHM(secon,nocorr) FWHM(prim)ERR   FWHM(secon)ERR samplesUp nSgms">$eresolFile

#
# Process input data files
#--------------------------
#foreach sep ($sepsStr)
foreach sep (3000)
	set inFile="sep${sep}sam_${exptime}s_${monoEkeV}keV.fits"
	#set inFile="sep3000sam_${exptime}s_${monoEkeV}keV.fits"
	set evtFile="events_sep${sep}sam_${root}.fits"
	echo "============================================="
	echo "\nUsing file: $inFile"
	echo "Using library: $libFile"
	echo "Using noisefile: $noiseFile"
	echo "Using Baseline: $baseline"
	echo "Setting evtFile: $evtFile"
	echo "============================================="
            
	# SIRENA processing
	#-------------------     
	#set filtersFile="filters_sep${sep}sam_${exptime}s_${rcmethod}_${monoEkeV}keV_${filterMeth}_${labelLib}.fits"
	#if(! -e $evtFile) then
	    set command="tesreconstruction Recordfile=$inFile TesEventFile=$evtFile Rcmethod="SIRENA" PulseLength=$nsamples LibraryFile=$libFile scaleFactor=$scaleFactor samplesUp=$samplesUp nSgms=$nSgms crtLib=0 mode=1 NoiseFile=$noiseFile FilterDomain=$fdomain FilterMethod=$filterMeth calibLQ=1 b_cF=$b_cF c_cF=$c_cF clobber=yes intermediate=0 baseline=$baseline  EnergyMethod=$energyMethod OFStrategy=FIXED OFLength=$oflength PixelType=$array" 
	    echo "command=$command"
	    tesreconstruction Recordfile=$inFile TesEventFile=$evtFile Rcmethod="SIRENA" PulseLength=$nsamples LibraryFile=$libFile scaleFactor=$scaleFactor samplesUp=$samplesUp nSgms=$nSgms crtLib=0 mode=1 NoiseFile=$noiseFile FilterDomain=$fdomain FilterMethod=$filterMeth calibLQ=1 b_cF=$b_cF c_cF=$c_cF clobber=yes intermediate=0 baseline=$baseline  EnergyMethod=$energyMethod OFStrategy=FIXED OFLength=$oflength PixelType=$array 
	    echo "SIRENA proc with $filterMeth$fdomain for $array with $energyMethod method and OFLength $oflength finished OK"
	#endif
	#cd $pwd
	#exit
	# check that detection was ok
	fkeypar fitsfile=$inFile keyword="NETTOT"
	set nsimpulses=`pget fkeypar value`
	fstatistic infile=$evtFile colname="GRADE1" rows="-" minval=0 >/dev/null
	set ndetpulses=`pget fstatistic numb`
	if($nsimpulses != $ndetpulses) then
	    echo "ERROR in detection for sim file $inFile (simpulses=$nsimpulses, detpulses=$ndetpulses)"
	    cd $pwd
	    exit
	endif      
            
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
	# create column SQR(Ei(biasCorr)-Eo)
	fcalc infile=$evt1 outfile=\!$evt1 clname=EbcDIFFSQR expr="(EeVbiasCorr-$monoEeV)^2"
	# calculate FWHM for primaries (BIAS corrected)
	fstatistic infile=$evt1 colname=EbcDIFFSQR rows="-" >/dev/null
	set sigma1sqr=`pget fstatistic mean`
	MATH fwhm = (sqrt($sigma1sqr) * 2.35)   # FWHM for PRIMARIES in eV
	set fwhm1=`printf "%.3f" $fwhm`
	
	# calculate FWHM for primaries (NO BIAS corrected)
	fcalc infile=$evt1 outfile=\!$evt1 clname=EDIFFSQR expr="(SIGNAL*1000.-$monoEeV)^2"
	fstatistic infile=$evt1 colname=EDIFFSQR rows="-" >/dev/null
	set sigma1sqr=`pget fstatistic mean`
	set sum1sqr=`pget fstatistic sum`
	MATH fwhm_nocorr = (sqrt($sigma1sqr) * 2.35)   # FWHM for PRIMARIES in eV
	set fwhm1_nocorr=`printf "%.3f" $fwhm_nocorr`
	
	#calculate error for fwhm1
	fstatistic infile=$evt1 colname=SIGNAL rows="-" >/dev/null
	set sigma_keV=`pget fstatistic sigma`
	#MATH fwhm_err = (2.35 * $sigma_keV * 1000. / sqrt($nrows))
	MATH fwhm_err = (2.35 * $sigma_keV * 1000. / sqrt(2*$nrows-2))
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
	# create column SQR(Ei(biasCorr)-Eo)
	fcalc infile=$evt2 outfile=\!$evt2 clname=EbcDIFFSQR expr="(EeVbiasCorr-$monoEeV)^2"
	# calculate FWHM for secondaries (BIAS corrected)
	fstatistic infile=$evt2 colname=EbcDIFFSQR rows="-" >/dev/null
	set sigma2sqr=`pget fstatistic mean`
	MATH fwhm = (sqrt($sigma2sqr) * 2.35)   # FWHM for SECONDARIES in eV
	set fwhm2=`printf "%.3f" $fwhm`
	# calculate FWHM for secondaries (NO BIAS corrected)
	fcalc infile=$evt2 outfile=\!$evt2 clname=EDIFFSQR expr="(SIGNAL*1000.-$monoEeV)^2"
	fstatistic infile=$evt2 colname=EDIFFSQR rows="-" >/dev/null
	set sigma2sqr=`pget fstatistic mean`
	MATH fwhm_nocorr = (sqrt($sigma2sqr) * 2.35)   # FWHM for SECONDARIES in eV
	set fwhm2_nocorr=`printf "%.3f" $fwhm_nocorr`

	#calculate error for fwhm2
	fstatistic infile=$evt2 colname=SIGNAL rows="-" >/dev/null
	set sigma_keV=`pget fstatistic sigma`
	#MATH fwhm_err = (2.35 * $sigma_keV * 1000. / sqrt($nrows))
	MATH fwhm_err = (2.35 * $sigma_keV * 1000. / sqrt(2*$nrows-2))
	set fwhm2_err=`printf "%.3f" $fwhm_err`

	echo "$sep   $fwhm1    $fwhm2   $ebias1  $ebias2  $fwhm1_nocorr  $fwhm2_nocorr $fwhm1_err  $fwhm2_err $samplesUp $nSgms">>$eresolFile
      
  end # infiles 
  cd $pwd
endif


