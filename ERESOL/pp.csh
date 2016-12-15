#!/usr/bin/csh

#
# Get energy resolution curves as a function of pulses pair separation

#
# IT USES INTEGRATION OF SIRENA INTO SIXTE, so SIXTE must be initialized
#
# RUN as:
#        source getEresolCurves.csh monoenergy tessim lib samplesUp nSgms
#
unset monoEkeV monoEeV calibDir pulseLength baseline tessim samplesUp nSgms scaleFactor libDir simDir noiseDir noiseFile SIRENAmeths
alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'

# Read input options
#---------------------
echo "Run as:"
echo "source getEresolCurves.csh  monoenergy tessim lib samplesUp nSgms F0/B0 nsamples"
echo ""
echo "       monoenergy: monochromatic energy (keV) of the pulses"
echo "       tessim: tessim version"
echo "       lib: 'monolib' o 'multilib'"
echo "       samplesUp: number of samples above threshold"
echo "       nSgms: number of sigmas to define threshold"
echo "       F0/B0: filter method (freq. bin 0 or subtract Baseline)"
echo "       nsamples: pulse length & noise samples"
echo ""   

if($#argv < 7) then 
    echo "Please complete input command line with all the parameters"
    echo "   source getEresolCurves.csh monoenergy tessim lib samplesUp nSgms F0/B0 nsamples"
    exit
endif    

# Read input parameters
#######################
set monoEkeV=$1 # file monochromatic energy (keV)
MATH monoEeV = ($monoEkeV*1000.)
set tessim=$2
set labelLib=$3 # monolib or multilib or nonoiselib
set samplesUp=$4  # 2
set nSgms=$5      # 10
set SIRENAmeths=$6
set nsamples=$7

# SIRENA definitions
#######################

set scaleFactor=0.005
#set rcmethod=$2 # PP or SIRENA or ALL
set rcmethod="SIRENA"
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

#
unset libFile
if($labelLib == "multilib") set libFile="${libDir}/libraryMultiE_PL${nsamples}_${tessim}.fits"
if($labelLib == "monolib" ) then
    #!!!!!!! TAKE CARE HERE : MONO LIBS
    #
    if($monoEkeV == 1) set libFile=${libDir}/libraryMultiE_PL${nsamples}_1keV_${tessim}.fits
    if($monoEkeV == 3) set libFile=${libDir}/libraryMultiE_PL${nsamples}_3keV_${tessim}.fits
    if($monoEkeV == 6) set libFile=${libDir}/libraryMultiE_PL${nsamples}_6keV_${tessim}.fits
    #
    #!!!!!!! TAKE CARE HERE : MONO LIBS
endif
if($labelLib == "nonoiselib") then
    if($monoEkeV == 1) set libFile=${libDir}/libraryMultiE_PL${nsamples}_1keV_${tessim}Nonoise.fits
    if($monoEkeV == 3) set libFile=${libDir}/libraryMultiE_PL${nsamples}_3keV_${tessim}Nonoise.fits
    if($monoEkeV == 6) set libFile=${libDir}/libraryMultiE_PL${nsamples}_6keV_${tessim}Nonoise.fits
    #
    #!!!!!!! TAKE CARE HERE : MONO LIBS
endif

set pwd=`pwd`
cd $simDir

set sepSamples=(0003 0005 0007 0010 0014 0019 0026 0037 0051 0070 0097 0135 0187 0260 0361 0420 0500 0700 1000 3000)
#set sepSamples=(3000)

#
# Create output ERESOL files
#----------------------------
set eresolFiles=""
#if("$2" == "ALL") set rcmethod="SIRENA"
if($rcmethod == "SIRENA") then
      foreach fmeth ($SIRENAmeths) 
	  set erS="eresol_${exptime}s_SIRENA${nsamples}_${monoEkeV}keV_${fmeth}_${labelLib}.dat"
	  if(-e $erS) rm $erS
	  set eresolFiles=($eresolFiles $erS)
      end
endif

# write header in eresol files
foreach erf ($eresolFiles)
      echo "samSep  FWHM(prim)  FWHM(secon)  EBIAS(prim)    EBIAS(sec)   FWHM(prim,nocorr)  FWHM(secon,nocorr) FWHM(prim)ERR   FWHM(secon)ERR samplesUp nSgms">$erf
end
  
#
# Process input data files
#--------------------------
foreach sep ($sepSamples)
      echo "============================================="
      set inFile="sep${sep}sam_${exptime}s_${monoEkeV}keV.fits"
      #set inFile="sep3000sam_${exptime}s_${monoEkeV}keV.fits"
      echo "\nUsing file: $inFile"
      echo "Using library: $libFile"
      echo "Using noisefile: $noiseFile"
      echo "Using Baseline: $baseline"
      echo "============================================="
      
      set evtFiles=""
      
      # SIRENA processing
      #-------------------
      if("$2" == "ALL")set rcmethod="SIRENA"
      if($rcmethod == "SIRENA") then
	      
	    foreach fmeth ($SIRENAmeths)
		#MATH sep = ($sep*1.)
		set evtFile="events_sep${sep}sam_${exptime}s_${rcmethod}_${monoEkeV}keV_${fmeth}_${labelLib}.fits"
		set filtersFile="filters_sep${sep}sam_${exptime}s_${rcmethod}_${monoEkeV}keV_${fmeth}_${labelLib}.fits"
		set command="tesreconstruction Recordfile=$inFile TesEventFile=$evtFile Rcmethod="SIRENA" PulseLength=$nsamples LibraryFile=$libFile scaleFactor=$scaleFactor samplesUp=$samplesUp nSgms=$nSgms crtLib=0 mode=1 NoiseFile=$noiseFile FilterDomain=F FilterMethod=$fmeth calibLQ=1 b_cF=$b_cF c_cF=$c_cF clobber=yes intermediate=0"
		echo "command=$command"
		tesreconstruction Recordfile=$inFile TesEventFile=$evtFile Rcmethod="SIRENA" PulseLength=$nsamples LibraryFile=$libFile scaleFactor=$scaleFactor samplesUp=$samplesUp nSgms=$nSgms crtLib=0 mode=1 NoiseFile=$noiseFile FilterDomain="F" FilterMethod=$fmeth calibLQ=1 b_cF=$b_cF c_cF=$c_cF clobber=yes intermediate=0 baseline=$baseline intermediate=1 filterFile=$filtersFile #resize=$sep
		echo "SIRENA proc with $fmeth finished"

		#cd $pwd
		#exit
		# check that detection was ok
		fkeypar fitsfile=$inFile keyword="NETTOT"
		set nsimpulses=`pget fkeypar value`
		#fkeypar fitsfile=$evtFile keyword="NAXIS2"
		#set ndetpulses=`pget fkeypar value`
		fstatistic infile=$evtFile colname="GRADE1" rows="-" minval=0 >/dev/null
		set ndetpulses=`pget fstatistic numb`
		if($nsimpulses != $ndetpulses) then
		    echo "ERROR in detection for sim file $inFile (simpulses=$nsimpulses, detpulses=$ndetpulses)"
		    cd $pwd
		    exit
		endif
		set evtFiles = ($evtFiles $evtFile)
		#exit
	    end
	    
      endif
            
      #-----------------------------------------
      # EVENT file processing to calculate FWHM
      #-----------------------------------------
      set i=1
      #echo "EVT files: $evtFiles"
      set nfiles=$#evtFiles
      while ($i <= $nfiles)
	set root=$evtFiles[$i]:r
	set evt=$evtFiles[$i]
	set evt1="${root}_primaries.fits"
	set evt2="${root}_secondaries.fits"
	#echo "Processing EVT file: $evt"
	if(-e $evt1) rm $evt1
	if(-e $evt2) rm $evt2
	
	#use only rows with GRADE1>0 (non initially truncated)
 	fselect infile=$evt outfile=\!$evt1 expr="#ROW%2!=0 && GRADE1>0"
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
	fselect infile=$evt outfile=\!$evt2 expr="#ROW%2==0 && GRADE1>0"
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

	echo "$sep   $fwhm1    $fwhm2   $ebias1  $ebias2  $fwhm1_nocorr  $fwhm2_nocorr $fwhm1_err  $fwhm2_err $samplesUp $nSgms">>$eresolFiles[$i]
	
	@ i++
      end #evtfiles
      
  end # infiles 
  cd $pwd
endif


