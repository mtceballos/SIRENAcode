#! /bin/bash
#
# run different reconstruction methods for all the energies/all record lengths
#
#  Option 1: calculate reconstructed energies and set coefficients table for gain scale curves (filter always 8192 + zero padding)
#  Option 2: with the Erec vs. Ecal coeffs, calculate FWHM corrected curves (filter always 8192 + zero padding)
#                      - correction performed with Gain Scale curves
#                      - correction performed with Gain Scale 2D surface
# Option 3: calculate reconstructed energies and set coefficients table for gain scale curves (reduced filter)
#  Option 4: with the Erec vs. Ecal coeffs, calculate FWHM corrected curves (reduced filter)
#                      - correction performed with Gain Scale curves
#                      - correction performed with Gain Scale 2D surface
#
clear
option=$1
instrument="LPA75um" #"LPA2shunt"
domain="T"
# ----------------------------------------------------------------------------------------------------------------------
# BBFB and NOISE and Jitter and samprate
bbfbStr="_bbfb" # "" or "_bbfb"
Lc="1" # inductance with respect to lcrit "1" or '0.35' or '0.5' or '0.7' 
jitterStr="_jitter"
noiseStr="" # "" or"_nonoise"
smpStr="_samprate2" # "" for samprate (156.25kHz)  or "_samprate2" for 78125Hz or "_samprate4" for 39062.5Hz
smpStr=""
coeffsFile="coeffs_polyfit.dat"  # or coeffs_poly2Dfit_pL8192 _OPTFILT8192.dat, coeffs_poly2Dfit_pL512 _OPTFILT8192.dat, coeffs_poly2Dfit_pL8192 _I2R8192.dat
# ----------------------------------------------------------------------------------------------------------------------

jitterParam=""
bbfbParam=""
dcmt=100
LcParam=""
if [ "$bbfbStr"  ==  "_bbfb" ];  then
    bbfbParam="--bbfb bbfb"
    jitterParam="--jitter jitter"
    dcmt=1
fi
if [ "$jitterStr"  ==  "_jitter" ] && [ $dcmt -gt 1 ];  then 
    jitterParam="--jitter jitter --decimation $dcmt"
    jitterStr="_jitter_dcmt$dcmt"
fi
noiseParam=""
if [ "$noiseStr"  ==  "_nonoise" ];  then
    noiseParam="--noise nonoise"
fi
LcStr="" #  '_0.35Lc' or '_0.5Lc' or '_0.7Lc'
if [ "$(echo ${Lc} '==' 0.7 | bc -l)" -eq 1 ] || [  "$(echo ${Lc} '==' 0.5 | bc -l)" -eq 1 ]   || [ "$(echo ${Lc} '==' 0.35 | bc -l)" -eq 1 ]; then
    LcParam="--Lc ${Lc}"
    LcStr="_${Lc}Lc"
elif [ ! "$(echo ${Lc} '==' 1 | bc -l)" -eq 1 ]; then
    echo "Error in inductance (L) value selection nLc=${Lc}Lc (nor 0.35Lc, 0.5Lc, 0.7Lc or 1Lc values given)"
    exit
fi
#########################################################
detMethod="AD"
detMethod="STC"
samplesUp=0
samplesDown=0
nSigmas=0

if [ "$smpStr"  ==  "_samprate2" ];  then
    separation=20000
    smpParam="--samprate samprate2"
    pulseLength=4096
    nSamples=4096
    largeFilter=4096
    if [ ${detMethod}  ==  "AD" ];  then
        nSgms=12 # BE CAREFUL!!! => Check the hard-coded slope criteria 
    elif [ ${detMethod}  ==  "STC" ]; then
        nSgms=4.5 #6
        samplesUp=2 #2 
        samplesDown=3 #3
    fi
    
elif [ "$smpStr"  ==  "_samprate4" ];  then
    separation=10000
    smpParam="--samprate samprate4"
    pulseLength=2048
    nSamples=2048
    largeFilter=2048
    if [ ${detMethod}  == "AD" ];  then
        nSgms=23.5 # BE CAREFUL!!! => Check the hard-coded slope criteria 
    elif [ ${detMethod} == "STC" ]; then
        nSgms=4 #6
        samplesUp=2 #2 
        samplesDown=3 #3
    fi
    
elif [ -z "$smpStr" ]; then  # samprate 156250Hz
    separation=40000
    smpParam=""
    pulseLength=8192
    nSamples=8192
    largeFilter=8192
    if [ ${detMethod}  ==  "AD" ];  then
        nSgms=4.4 # BE CAREFUL!!! => Check the hard-coded slope criteria 
    elif [ ${detMethod}  ==  "STC" ];  then
        nSgms=3.5 #5
        samplesUp=3 #3
        samplesDown=4 #3
    fi
fi
flengths=($pulseLength)
echo "Using det=$detMethod nSgms=$nSgms smplsUp=$samplesUp smplsDown=$samplesDown "
echo "(smpParam: $smpParam) (jitterParam=$jitterParam) (noiseParam=$noiseParam) (bbfb=$bbfbParam) (Lc=$LcParam) "

# USE one BY ONE
fixedlib6OF_OPTFILT=1
fixedlib6OF_OPTFILTNM=0
fixedlib6OF_I2RNM=0
fixedlib6OF_I2R=0
fixedlib6OF_I2RNOL=0
fixedlib6OF_I2RFITTED=0
multilibOF_WEIGHTN=0
multilib_WEIGHTN=0
multilib_WEIGHT=0

if  [ $fixedlib6OF_OPTFILT -eq 1 ]; then
    # Global OPTFILT AC fixedlib6 OFLib=yes
    # -------------------------------------------
    lib="fixedlib6OF"
    meth="OPTFILT"
    nSimPulsesLib=20000 
    flengths=(8192 2048 4096 1024 512 256 128)
    flengths=(8192)
    
elif [ $fixedlib6OF_OPTFILTNM -eq 1 ]; then
    
    # Global OPTFILT AC fixedlib6 OFLib=yes NoiseMAT
    # -------------------------------------------------------
    lib="fixedlib6OF"
    meth="OPTFILTNM50000"   
    nSimPulsesLib=20000
    flengths=(8192 4096)     # 8192 gives bad results in reconstruction
    
elif [ $fixedlib6OF_I2RNM -eq 1 ]; then
    
    # Global I2R AC fixedlib6 OFLib=yes NoiseMAT
    # -------------------------------------------------------
    lib="fixedlib6OF"
    meth="I2RNM150000"
    nSimPulsesLib=20000
    flengths=(8192 4096 2048)   # 8192 give bad results in reconstruction
    
elif [ $fixedlib6OF_I2R -eq 1 ]; then
    
    # Global OPTFILT I2R fixedlib6  --OFLib yes 
    # ---------------------------------------------
    lib="fixedlib6OF"
    meth="I2R"
    nSimPulsesLib=20000
    flengths=(8192 2048 4096 1024 512 256 128)

elif [ $fixedlib6OF_I2RNOL -eq 1 ]; then

# Global OPTFILT I2RNOL fixedlib6  --OFLib yes 
    # ---------------------------------------------------
    lib="fixedlib6OF"
    meth="I2RNOL"
    nSimPulsesLib=20000

elif [ $fixedlib6OF_I2RFITTED -eq 1 ]; then
    
    # Global OPTFILT I2RFITTED fixedlib6  --OFLib yes 
    # -------------------------------------------------------
    lib="fixedlib6OF"
    meth="I2RFITTED"
    nSimPulsesLib=20000

elif [ $multilib_WEIGHT -eq 1 ]; then
    
    # WEIGHT AC multilib
    #---------------------
    lib="multilib"
    meth="WEIGHT"
    nSimPulsesLib=200000

elif [ $multilib_WEIGHTN -eq 1 ]; then
    
    # WEIGHTN AC multilib
    #---------------------
    lib="multilib"
    meth="WEIGHTN"
    nSimPulsesLib=200000

elif [ $multilibOF_WEIGHTN -eq 1 ]; then
    
    # WEIGHTN AC OF
    #---------------------
    lib="multilibOF"
    meth="WEIGHTN"
    nSimPulsesLib=200000
fi

lags=(1 1 1 1 1 1 1 1 1 1) 
nSimPulses=20000
energies=(0.2 0.5 1 2 3 4 5 6 7 8)
#energies=( 1 2 3 4 5 6 7 8)
nenergies=${#energies[@]}
nfls=${#flengths[@]}

if  [ $option -eq 1 ]; then
	# =================================================
	# To create gain scale curves & coefficients table
	# =================================================
	# All methods can be run OFLib=yes because length to be used is maximum (pulseLength=largeFilter)
	#set tstartPulse1All=(999 999 1000 1000 1000 1000 1000 1000 1000 1000) # default 0 in getEresolCurves
        #set tstartPulse1All=(0 0 0 0 0 0 0 0 0 0) # default 0 in getEresolCurves
	nSimPulses=5000
	#energies=(0.2 0.5)
        
	for (( ie=0; ie<$nenergies; ie++ )); do
		mono1EkeV=${energies[$ie]}
		mono2EkeV=0
		echo "Runing gain scale creation (Erecons vs Einput)  for energy=$mono1EkeV"
		#tstartPulse1=${tstartPulse1All[$ie]}
		tstartPulse1=0
		tstartPulse2=0
		lgs=${lags[$ie]}
                
		for (( ifl=0; ifl<$nfls; ifl++ )); do
                    filterLength=${flengths[$ifl]}
                    echo "Launching $meth $lib for $mono1EkeV "
                    logf="${instrument}_${detMethod}_${meth}_${filterLength}_${mono1EkeV}${smpStr}${jitterStr}${noiseStr}${bbfbStr}.log"
                    command="recon_resol.py --pixName ${instrument} --labelLib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $filterLength --fdomain $domain --tstartPulse1 $tstartPulse1  --libTmpl LONG --detMethod $detMethod --nSgms $nSgms --samplesUp ${samplesUp} --samplesDown ${samplesDown} --filterLength $pulseLength $smpParam --resultsDir gainScale ${jitterParam} ${noiseParam} ${bbfbParam} ${LcParam}"
                    nohup python $command >& $logf &
                    echo "Command=python $command >> $logf" 
                done
               # sleep 80
	done
fi

if [ $option -eq 2 ]; then
	# To create resolution curves with the coefficients from gain scale curves or surfaces
	# ====================================================
	nSimPulses=2000
        for (( ie=0; ie<$nenergies; ie++ )); do
		mono1EkeV=${energies[$ie]}
		mono2EkeV=0
		echo "Runing method comparison (FWHM vs Energy) for energy=$mono1EkeV"
		#tstartPulse1=${tstartPulse1All[$ie]}
		tstartPulse1=0
		tstartPulse2=0
		lgs=${lags[$ie]}
                
		for (( ifl=0; ifl<$nfls; ifl++ )); do
                    filterLength=${flengths[$ifl]}
                    echo "Launching $meth $lib for $mono1EkeV  for flength= ${filterLength}"
                    logf="${instrument}_$meth${lib}_${detMethod}${filterLength}_${mono1EkeV}${smpStr}${jitterStr}${noiseStr}${bbfbStr}.log"
                    command="recon_resol.py --pixName ${instrument} --labelLib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $filterLength --fdomain $domain --tstartPulse1 $tstartPulse1  --libTmpl LONG --detMethod $detMethod --nSgms $nSgms --samplesUp ${samplesUp} --samplesDown ${samplesDown} --filterLength $pulseLength ${smpParam} ${jitterParam} ${noiseParam}  ${bbfbParam} ${LcParam} --coeffsFile coeffs_polyfit.dat"
                    nohup python $command >&  $logf &
                    echo "Command=python $command >> $logf" 
                done    
                #sleep 80
	done
fi

if  [ $option -eq 3 ]; then
	# =================================================
	# To create gain scale curves & coefficients table
	# =================================================
	# All methods can be run OFLib=yes because length to be used is maximum (pulseLength=largeFilter)
	#set tstartPulse1All=(999 999 1000 1000 1000 1000 1000 1000 1000 1000) # default 0 in getEresolCurves
        #set tstartPulse1All=(0 0 0 0 0 0 0 0 0 0) # default 0 in getEresolCurves
	nSimPulses=5000
	#energies=(0.2 0.5)
        
	for (( ie=0; ie<$nenergies; ie++ )); do
		mono1EkeV=${energies[$ie]}
		mono2EkeV=0
		echo "Runing gain scale creation (Erecons vs Einput)  for energy=$mono1EkeV"
		#tstartPulse1=${tstartPulse1All[$ie]}
		tstartPulse1=0
		tstartPulse2=0
		lgs=${lags[$ie]}
                
		for (( ifl=0; ifl<$nfls; ifl++ )); do
                    filterLength=${flengths[$ifl]}
                    echo "Launching $meth $lib for $mono1EkeV "
                    logf="${instrument}_${detMethod}_${meth}_${filterLength}_${mono1EkeV}${smpStr}${jitterStr}${noiseStr}${bbfbStr}.log"
                    command="recon_resol.py --pixName ${instrument} --labelLib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain $domain --tstartPulse1 $tstartPulse1  --libTmpl LONG --detMethod $detMethod --nSgms $nSgms --samplesUp ${samplesUp} --samplesDown ${samplesDown} --filterLength $filterLength $smpParam --resultsDir gainScale ${jitterParam} ${noiseParam} ${bbfbParam} ${LcParam}"
                    nohup python $command >& $logf &
                    echo "Command=python $command >> $logf" 
                done
                #sleep 80
	done
fi

if [ $option -eq 4 ]; then
	# To create resolution curves with the coefficients from gain scale curves
	# =============================================
	nSimPulses=2000
        for (( ie=0; ie<$nenergies; ie++ )); do
		mono1EkeV=${energies[$ie]}
		mono2EkeV=0
		echo "Runing method comparison (FWHM vs Energy) for energy=$mono1EkeV"
		#tstartPulse1=${tstartPulse1All[$ie]}
		tstartPulse1=0
		tstartPulse2=0
		lgs=${lags[$ie]}
                
		for (( ifl=0; ifl<$nfls; ifl++ )); do
                    filterLength=${flengths[$ifl]}
                    echo "Launching $meth $lib for $mono1EkeV  for flength= ${filterLength}"
                    logf="${instrument}_$meth${lib}_${detMethod}${filterLength}_${mono1EkeV}${smpStr}${jitterStr}${noiseStr}${bbfbStr}.log"
                    command="recon_resol.py --pixName ${instrument} --labelLib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain $domain --tstartPulse1 $tstartPulse1  --libTmpl LONG --detMethod $detMethod --nSgms $nSgms --samplesUp ${samplesUp} --samplesDown ${samplesDown} --filterLength $filterLength ${smpParam} ${jitterParam} ${noiseParam}  ${bbfbParam} ${LcParam} --coeffsFile coeffs_polyfit.dat"
                    nohup python $command >&  $logf &
                    echo "Command=python $command >> $logf" 
                done    
                sleep 80
	done
fi

