#!/usr/bin/csh

#
# Get energy resolution curves as a function of pulses pair separation

#
# IT USES INTEGRATION OF SIRENA INTO SIXTE, so SIXTE must be initialized
#
# RUN as:
#        source getEresolCurves.csh monoenergy tessim lib samplesUp nSgms
#
unset monoEkeV monoEeV calibDir pulseLength baseline tessim samplesUp nSgms scaleFactor libDir simDir noiseDir noiseFile SIRENAmeths nsimpulses ndetpulses sigmasMin sigmasMax

alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'


# Read input options
#---------------------
echo "Run as:"
echo "source getEresolTrios.csh  monoenergy array lib samplesUp nSgms [F0/B0] nsamples"
echo ""
echo "       monoenergy: monochromatic energy (keV) of the pulses"
echo "       array: array configuration"
echo "       lib: 'monolib' o 'multilib'"
echo "       samplesUp: number of samples above threshold"
echo "       nSgms: number of sigmas to define threshold"
echo "       [F0/B0]: filter method (freq. bin 0 or subtract Baseline)"
echo "       nsamples: number of samples for High resolution record and noise "
echo ""

if($#argv < 7) then 
    echo "Please complete input command line with all the parameters"
    echo "   source getEresolTrios.csh monoenergy array lib samplesUp nSgms F0|B0 nsamples"
    exit
endif    

# Read input parameters
#######################
set monoEkeV=$1 # file monochromatic energy (keV)
MATH monoEeV = ($monoEkeV*1000.)
set array=$2
set tessim="tessim$2"
set labelLib=$3 # monolib or multilib or nonoiselib
set samplesUp=$4  # 2
set nSgms=$5      # 10
set SIRENAmeth=$6
set nsamples=$7 # pulse_length & noise samples
set PreBufferSize=1000
set Fil=""

# Define tmp dir
set tmpdir="eresolTrios${array}_${monoEkeV}"
mkdir -p /tmp/${tmpdir}/pfiles
setenv PFILES "/tmp/${tmpdir}/pfiles;${HEADAS}/syspfiles:${SIXTE}/share/simput/pfiles:$SIXTE/share/sixte/pfiles"
echo "PFILES=$PFILES"
setenv HEADASNOQUERY
setenv HEADASPROMPT /dev/null

# SIRENA definitions
#######################
#set scaleFactor=0.02
set b_cF=1.
set c_cF=0.

#echo "Using baseline=$baseline"


set libDir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/LIBRARIES/${tessim}/OPTFILT"
set exptime="100" #s 
set simDir="./TRIOS/${tessim}"
set noiseDir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/NOISE/${tessim}"
set noiseFile=${noiseDir}/noise${nsamples}samples_${tessim}_B0_100s_pairscps.fits
fkeypar fitsfile=${noiseFile}+1 keyword="BASELINE"
set baseline=`pget fkeypar value`
unset libFile
if($labelLib == "multilib") set libFile="${libDir}/libraryMultiE_PL${nsamples}_${tessim}.fits"
if($labelLib == "monolib") then
    if($monoEkeV == 1) set libFile=${libDir}/libraryMultiE_PL${nsamples}_1keV_${tessim}.fits
    if($monoEkeV == 3) set libFile=${libDir}/libraryMultiE_PL${nsamples}_3keV_${tessim}.fits
    if($monoEkeV == 6) set libFile=${libDir}/libraryMultiE_PL${nsamples}_6keV_${tessim}.fits
endif

if($labelLib == "nonoiselib") then
    if($monoEkeV == 1) set libFile=${libDir}/libraryMultiE_PL${nsamples}_1keV_${tessim}Nonoise.fits
    if($monoEkeV == 3) set libFile=${libDir}/libraryMultiE_PL${nsamples}_3keV_${tessim}Nonoise.fits
    if($monoEkeV == 6) set libFile=${libDir}/libraryMultiE_PL${nsamples}_6keV_${tessim}Nonoise.fits
endif
set pwd=`pwd`
cd $simDir

set PixType="file:/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/tespixels.fits[${array}]"
if($array == "SPA") then

    #set sepsStr=(0004 0008 0016 0032 0064 0128 0256 0512 1024 2048 ) # 4096 )# 8192)
    #set sepsStr=(00002 00003 00004 00005 00007 00008 00011 00013 00017 00022 00028 00035 00045 00057 00072 00092 00116 00148  00188 00238 00303 00384 00488 00620 00787 01000 01270 01613 02048)
    set sepsStr=(00003 00004 00005 00007 00008 00011 00013 00017 00022 00028 00035 00045 00057 00072 00092 00116 00148  00188 00238 00303 00384 00488 00620 00787 01000 01270 01613 02048)

    if($monoEkeV ==  6) set sepsStr=(00022 00028 00035 00045 00057 00072 00092 00116 00148  00188 00238 00303 00384 00488 00620 00787 01000 01270 01613 02048) #saturation!
    set scaleFactor=0.005
    
else if($array == "LPA1") then
    set sepsStr=(00004 00005 00007 00010 00013 00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287)
    #set sepsStr=(00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178)
    if($monoEkeV ==  6) set sepsStr=(00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287) #to avoid saturation
    set scaleFactor=0.005
    
else if($array == "LPA2") then
   
    #set sepsStr=(00004 00005 00007 00010 00013 00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287 11115 14910 20000)
    set sepsStr=(00004 00005 00007 00010 00013 00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287 11115 14910 20000)
    set scaleFactor=0.005
    
else if($array == "LPA3") then
    
    #set sepsStr=(00004 00005 00007 00010 00013 00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287 11115 14910 20000)
    set sepsStr=(00010 00013 00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287 11115 14910 20000) # very complx detection # high DE values for <00010
    

    set scaleFactor=0.02    
    if($labelLib == "monolib") then
	if($monoEkeV == 1) set libFile=${libDir}/libraryMultiE_PL${nsamples}_1keV_${tessim}${Fil}.fits
	if($monoEkeV == 6) set libFile=${libDir}/libraryMultiE_PL${nsamples}_6keV_${tessim}${Fil}.fits
    endif
    if($labelLib == "multilib") set libFile="${libDir}/libraryMultiE_PL${nsamples}_${tessim}${Fil}.fits"

endif
#
# Create output ERESOL files
#----------------------------
set eresolFile="${pwd}/TRIOS/${tessim}/eresol_${exptime}s_SIRENA${nsamples}_${monoEkeV}keV_${SIRENAmeth}_${labelLib}.dat"
set eresolFileFinal="${pwd}/TRIOS/eresolDAT/${tessim}/eresol_${exptime}s_SIRENA${nsamples}_${monoEkeV}keV_${SIRENAmeth}_${labelLib}.dat"

#set eresolFile="pp.dat"
if(-e $eresolFile) rm $eresolFile

# write header in ereso file
echo "samSep12   samSep23  FWHM(middle)  FWHM(middle)ERR  EBIAS(middle)  samplesUp nSgms nsamples">$eresolFile
  
#
# Process input data files
#--------------------------
set sigmas=$nSgms
foreach sepA ($sepsStr)
#foreach sepA (02048)
    MATH sep12 = ($sepA * 1.)
    
    foreach sepB ($sepsStr)
    #foreach sepB (02048)
	MATH sep23 = ($sepB * 1.)

	MATH tstartPulse1 = ($PreBufferSize + 1)
	MATH tstartPulse2 = ($tstartPulse1 + $sep12)
	MATH tstartPulse3 = ($tstartPulse2 + $sep23)
	
	if($tstartPulse1 > 0) then
	  set nSgms=0
	  set samplesUp=0
	  set scaleFactor=0
	  set Fil=""
	endif

	set sigmasMin=0
	if($monoEkeV == 1)then 
	    set sigmasMax=50
	else
	    set sigmasMax=100
	endif
	set iter=0
	
	echo "============================================="
	set inFile="sep${sepA}sep${sepB}sam_${exptime}s_${monoEkeV}keV.fits"
	echo "\nUsing file: $inFile"
	echo "Using library: $libFile"
	echo "Using noisefile: $noiseFile"
	echo "============================================="
	
	
	# SIRENA processing
	#-------------------    
	set evtFile="events_sep${sepA}sep${sepB}sam_${exptime}s_SIRENA_${monoEkeV}keV_${SIRENAmeth}_${labelLib}.fits"
	if(! -e $evtFile) then
	    #set evtFile="pp.evt"
	    set command="tesreconstruction Recordfile=$inFile TesEventFile=$evtFile Rcmethod="SIRENA" PulseLength=$nsamples LibraryFile=$libFile scaleFactor=$scaleFactor samplesUp=$samplesUp nSgms=$sigmas crtLib=0 mode=1 NoiseFile=$noiseFile FilterDomain="F" FilterMethod=$SIRENAmeth calibLQ=1 b_cF=$b_cF c_cF=$c_cF clobber=yes intermediate=0 baseline=$baseline EventListSize=1000 EnergyMethod=NOLAGS tstartPulse1=$tstartPulse1 tstartPulse2=$tstartPulse2 tstartPulse3=$tstartPulse3"
	    echo "command=$command"
	    $command
	    
	endif
	
	#cd $pwd
	#exit
	# check that detection was ok
	fkeypar fitsfile=$inFile keyword="NETTOT"
	set nsimpulses=`pget fkeypar value`
	#set nsimpulses = 4638
	#fparkey value=$nsimpulses fitsfile=$inFile keyword=NETTOT
	#continue
	fkeypar fitsfile=$evtFile keyword="NAXIS2"
	set nrows=`pget fkeypar value`
	set ndetpulses=0
	if ($nrows > 0) then
	    fstatistic infile=$evtFile colname="GRADE1" rows="-" minval=0 >/dev/null
	    set ndetpulses=`pget fstatistic numb`
	endif
# 	if($nsimpulses != $ndetpulses) then
# 	    echo "ERROR in detection for sim file $inFile (simpulses=$nsimpulses, detpulses=$ndetpulses)"
# 	    cd $pwd
# 	    exit
# 	endif
        while ($nsimpulses != $ndetpulses)
	    if($iter>50) break
	    MATH iter = ($iter+1)
	    echo "With nSgms=$sigmas => sim/det: $nsimpulses/$ndetpulses" 
	    if($nsimpulses<$ndetpulses) then
		#increase sigmas
		set sigmasMin=$sigmas
	    else if($nsimpulses>$ndetpulses) then
		#decrease sigmas
		set sigmasMax=$sigmas
	    endif
	    MATH sigmas = ($sigmasMin + ($sigmasMax-$sigmasMin)/2.)
	    echo "   Trying with nSgms=$sigmas (sigmasMin/Max:$sigmasMin/$sigmasMax)"
	    tesreconstruction Recordfile=$inFile TesEventFile=$evtFile Rcmethod="SIRENA" PulseLength=$nsamples LibraryFile=$libFile scaleFactor=$scaleFactor samplesUp=$samplesUp nSgms=$sigmas crtLib=0 mode=1 NoiseFile=$noiseFile FilterDomain="F" FilterMethod=$SIRENAmeth calibLQ=1 b_cF=$b_cF c_cF=$c_cF clobber=yes intermediate=0 baseline=$baseline EventListSize=1000 EnergyMethod=NOLAGS tstartPulse1=$tstartPulse1 tstartPulse2=$tstartPulse2 tstartPulse3=$tstartPulse3
	    
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
	    echo "Reconstruction NOT finished with nSgms=$sigmas"
	    continue
	endif
        echo "Reconstruction finished with nSgms=$sigmas"
        
	#-----------------------------------------
	# EVENT file processing to calculate FWHM
	#-----------------------------------------
	set root=$evtFile:r
	set evt1="${root}_primaries.fits"
	set evt2="${root}_secondaries.fits"
	set evt3="${root}_tertiaries.fits"
	#echo "Processing EVT file: $evtFile"
	if(-e $evt1) rm $evt1
	if(-e $evt2) rm $evt2
	if(-e $evt3) rm $evt3
	
	
	# print key with sigmas/nsamplesUp params
	fparkey fitsfile=${evtFile} value=$nSgms keyword="NSGMS" add=y
	fparkey fitsfile=${evtFile} value=$samplesUp keyword="SMPLSUP" add=y
	
	# SECONDARIES BIAS and FWHM 
	#-------------------------
	#use only rows with GRADE1>0 (non initially truncated) and storing middle pulses
 	fselect infile=$evtFile outfile=\!$evt2 expr="#ROW%3==2 && GRADE1>0"
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

	echo "$sepA  $sepB  $fwhm2  $fwhm2_err $ebias2  $samplesUp $nSgms $nsamples">>$eresolFile
      
  end # sep12
end #sep23
#make a backup copy
cp $eresolFile $eresolFileFinal
cd $pwd


