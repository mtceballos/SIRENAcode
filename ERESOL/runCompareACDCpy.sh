#! /bin/bash
#
# run different reconstruction methods for all the energies/all record lengths
#
#  Option 1: calculate reconstructed energies and set coefficients table for gain scale curves (filter always 8192 + zero padding)
#  Option 2: with the Erec vs. Ecal coeffs, calculate FWHM corrected curves (filter always 8192 + zero padding)
#                      - correction performed with Gain Scale curves
#                      - correction performed with Gain Scale 2D surface
#
# Option 3: calculate reconstructed energies and set coefficients table for gain scale curves (reduced filter + optionally preBuffer)
# Option 4: with the Erec vs. Ecal coeffs, calculate FWHM corrected curves (reduced filter + optionally preBuffer)
#                      - correction performed with Gain Scale curves
#                      - correction performed with Gain Scale 2D surface
#
clear
option=$1
instrument="LPA75um" #"LPA2shunt"
domain="T"
M82Str="_M82_040"
NewParStr="" #"_NewPar" (for Tbath sensitivity)
#NewParStr="_040" # (for bias voltage sensitivity) & baseADU sensitivity
#ctStr="_ct" #  for centre-replaced-by-constant filters  or "" for normal/original filters or 
#ctStr="_fit" # for filters derived from large filter (after a fitting process)
ctStr=""

# ----------------------------------------------------------------------------------------------------------------------
# BBFB and NOISE and Jitter and samprate and preBuffer and LAGS
bbfbStr="_bbfb" # "" or "_bbfb"
Lc="1" # inductance with respect to lcrit "1" or '0.35' or '0.5' or '0.7' 
jitterStr="_jitter" # "" or "_jitter"
noiseStr="" # "" or"_nonoise"
smpStr="_samprate2" # "" for samprate (156.25kHz)  or "_samprate2" for 78125Hz or "_samprate4" for 39062.5Hz
smpStr=""
preBuffer=0 # 0 or 75
pBstr=""
lags=1
lagsStr=""

#LbT=0.64E-3 # 100 samples to average baseline or "0" if baseline is not to be subtracted
#LbT=3.2E-3 # 500 samples to average baseline or "0" if baseline is not to be subtracted
#LbT=5.76E-3 # 900 samples
LbT=0
LbTstr=""
B0=1000 # if >0 run with B0 (baseline subtraction) instead of F0
B0=0

calib="1D" # "1D" for gainScale curve or "2D" for "gain scale surface"
 # energies to be calibrated
if [ "${calib}" == "1D" ]; then
    energies=(0.2 0.5 1 2 3 4 5 6 7 8)
elif [ "${calib}" == "2D" ]; then
    energies=(7)
fi

# ----------------------------------------------------------------------------------------------------------------------

jitterParam=""
bbfbParam=""
dcmt=100
LcParam=""
pBparam=""
B0param=""
LbTparam=""
sum0Param="--Sum0Filt 0"  # No need to 0-sum filter
lagsParam=""

if [ "$lags" == 0 ]; then
        lagsStr="_nolags"
        lagsParam="--lags 0"
fi  

if [ ${preBuffer} -gt 0 ]; then
    pBstr="_pB${preBuffer}"
fi
if [ ${B0} -gt 0 ]; then
    B0str="_B0-${B0}"
fi
if [ "$bbfbStr"  ==  "_bbfb" ];  then
    bbfbParam="--bbfb bbfb"
    #jitterParam="--jitter jitter"
    dcmt=1
fi
if [ "$jitterStr"  ==  "_jitter" ]; then
    jitterParam="--jitter jitter"
fi    
if [[ "$NewParStr"  ==  "_NewPar" ||  "$NewParStr"  ==  "_040" ||  "$NewParStr"  ==  "_040_ct" ]];  then
    bbfbParam="--bbfb bbfb${NewParStr}"
    jitterStr=""
fi
if [ "$M82Str"  ==  "_M82_040" ] ; then
        jitterStr="_jitter${M82Str}"
        jitterParam="--jitter jitter${M82Str}"
        bbfbParam=""
        bbfbStr=""
fi

LbTparam=""
if [ ${LbT} != 0 ]; then
    LbTstr="_LbT${LbT}"
    LbTparam="--LbT=${LbT}"
fi
ctParam=""
if [ "$ctStr"  ==  "_ct" ]; then
    ctParam="--ct ct"
fi    
if [ "$ctStr"  ==  "_fit" ]; then
    ctParam="--ct fit"
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
    return 1
fi

if [ "${calib}" == "1D" ]; then
    coeffsFile="coeffs_polyfit.dat"  
    coeffsFile="coeffs_polyfit_methods${NewParStr}${LbTstr}.dat" 
    coeffsFile="coeffs_polyfit_methods${M82Str}.dat" 
    #coeffsFile="coeffs_polyfit_methods_040_nojitter_nonoise_ct.dat"
    #coeffsFile="coeffs_polyfit_methods${ctStr}.dat"
    #coeffsFile="coeffs_polyfit_methods_040_Rs.dat"
    #energies=(7)
fi
#coeffsFile="spline_forfit_methods_longFilter_zeroPadding_ADC_I2R.json"

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
        nSgms=4.5 #6    #HARDCODED in auxpy.py
        samplesUp=2 #2   #HARDCODED in auxpy.py 
        samplesDown=3 #3  #HARDCODED in auxpy.py
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
        nSgms=4 #6  #HARDCODED in auxpy.py
        samplesUp=2 #2   #HARDCODED in auxpy.py
        samplesDown=3 #3  #HARDCODED in auxpy.py
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
        nSgms=3.5 #5  #HARDCODED in auxpy.py
        samplesUp=3 #3  #HARDCODED in auxpy.py
        samplesDown=4 #3  #HARDCODED in auxpy.py
    fi
fi

flengths=($pulseLength)

# USE one BY ONE
fixedlib6OF_OPTFILT=1    # includes 0-padding with SUM/=0 (no action over SUM(filter))
fixedlib6OF_OPTFILT_B0=0     # includes 0-padding with baseline subtraction (read from noise file)
fixedlib6OF_OPTFILT_0SUM=0     # includes 0-padding with SUM(filter)=0 
fixedlib6OF_OPTFILT_pB=0  # filters done with preBuffer 
fixedlib6OF_OPTFILTNM=0
fixedlib6OF_I2RNM=0
fixedlib6OF_I2R=1
fixedlib6OF_I2R_pB=0
fixedlib6OF_I2R_pB0p=0 # 0 padding + preBuffer
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
    flengths=(8192 4096 2048 1024 512 256 128 4)
    flengths=(512 256)
    preBuffer=0
    pBparam=""
     (( ${B0} > 0 )) && { echo "Error: if B0 > 0 you should run fixedlib6OF_OPTFILT_B0"; return 1;}
    

elif  [ $fixedlib6OF_OPTFILT_B0 -eq 1 ]; then
    # Global OPTFILT AC fixedlib6 OFLib=yes B0
    # -------------------------------------------
    lib="fixedlib6OF"
    meth="OPTFILT"
    nSimPulsesLib=20000 
    flengths=(8192)
    preBuffer=0
    pBparam=""
    B0param="--B0=${B0}"

elif  [ $fixedlib6OF_OPTFILT_0SUM -eq 1 ]; then
    # Global OPTFILT AC fixedlib6 OFLib=yes
    # -------------------------------------------
    lib="fixedlib6OF"
    meth="OPTFILT"
    nSimPulsesLib=20000 
    flengths=(512)
    preBuffer=0
    pBparam=""
    sum0Param="--Sum0Filt 1"

elif  [ $fixedlib6OF_OPTFILT_pB -eq 1 ]; then
    # Global OPTFILT AC fixedlib6 OFLib=yes preBuffer
    # --------------------------------------------------
    lib="fixedlib6OF"
    meth="OPTFILT"
    nSimPulsesLib=20000 
    flengths=(4096 2048 1024 512 256 128)
    pBparam="--preBuffer ${preBuffer}"

elif [ $fixedlib6OF_OPTFILTNM -eq 1 ]; then
    
    # Global OPTFILT AC fixedlib6 OFLib=yes NoiseMAT
    # -------------------------------------------------------
    lib="fixedlib6OF"
    meth="OPTFILTNM150000"   
    nSimPulsesLib=200000
    pulseLength=1024
    nSamples=1024
    flengths=(1024)    
    
    
elif [ $fixedlib6OF_I2RNM -eq 1 ]; then
    
    # Global I2R AC fixedlib6 OFLib=yes NoiseMAT
    # -------------------------------------------------------
    lib="fixedlib6OF"
    meth="I2RNM150000"
    nSimPulsesLib=200000
    pulseLength=1024
    nSamples=1024
    flengths=(1024)    
    
elif [ $fixedlib6OF_I2R -eq 1 ]; then
    
    # Global OPTFILT I2R fixedlib6  --OFLib yes 
    # ---------------------------------------------
    lib="fixedlib6OF"
    meth="I2R"
    nSimPulsesLib=20000
    flengths=(8192 4096 2048 1024 512 256 128)
    #flengths=(4)
    preBuffer=0
    pBparam=""

elif  [ $fixedlib6OF_I2R_pB -eq 1 ]; then
    # Global I2R AC fixedlib6 OFLib=yes preBuffer
    # --------------------------------------------------
    lib="fixedlib6OF"
    meth="I2R"
    nSimPulsesLib=20000 
    flengths=(4096 2048 1024 512 256 128)
    pBparam="--preBuffer ${preBuffer}"

elif  [ $fixedlib6OF_I2R_pB0p -eq 1 ]; then
    # Global I2R AC fixedlib6 OFLib=yes preBuffer + 0-padding
    # ---------------------------------------------------------
    lib="fixedlib6OF"
    meth="I2R"
    nSimPulsesLib=20000 
    flengths=(4096 2048 1024 512 256 128)
    pBparam="--preBuffer ${preBuffer}"
    pBstr="_pB${preBuffer}0p"
    
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
    flengths=(8192 4096 2048 1024 512 256 128)
    flengths=(8192)

elif [ $multilib_WEIGHT -eq 1 ]; then
    
    # WEIGHT AC multilib
    #---------------------
    lib="multilib"
    meth="WEIGHT"
    nSimPulsesLib=200000
    #nSimPulsesLib=20000
    pulseLength=1024
    nSamples=1024
    flengths=(1024) 
    
elif [ $multilib_WEIGHTN -eq 1 ]; then
    
    # WEIGHTN AC multilib
    #---------------------
    lib="multilib"
    meth="WEIGHTN"
    nSimPulsesLib=200000
    #nSimPulsesLib=20000
    pulseLength=1024
    nSamples=1024
    flengths=(1024) 
    
elif [ $multilibOF_WEIGHTN -eq 1 ]; then
    
    # WEIGHTN AC OF
    #---------------------
    lib="multilibOF"
    meth="WEIGHTN"
    nSimPulsesLib=200000
    #nSimPulsesLib=20000
    pulseLength=1024
    nSamples=1024
    flengths=(1024)
fi
echo "Using det=$detMethod nSgms=$nSgms smplsUp=$samplesUp smplsDown=$samplesDown "
echo "(smpParam: $smpParam) (jitterParam=$jitterParam) (noiseParam=$noiseParam) (bbfbParam=$bbfbParam) (LcParam=$LcParam) (pBparam=$pBparam) (lags=$lags) (LbTparam=$LbTparam) (ctParam=$ctParam) (B0=$B0param)"
nenergies=${#energies[@]}
nfls=${#flengths[@]}

###########    LONG FILTER #####################################################
if  [ $option -eq 1 ]; then
	# =================================================
	# To create gain scale curves & coefficients table
	# =================================================
	# All methods can be run OFLib=yes because length to be used is maximum (pulseLength=largeFilter)
	nSimPulses=5000
	#nSimPulses=1
        (( ${preBuffer} > 0 )) && { echo "Warning: You are running options '1' and '2' for 0-pad +  preBuffer"; }
	for (( ie=0; ie<$nenergies; ie++ )); do
		mono1EkeV=${energies[$ie]}
		mono2EkeV=0
		echo "Runing gain scale creation (Erecons vs Einput)  for energy=$mono1EkeV"
		#tstartPulse1=${tstartPulse1All[$ie]}
		tstartPulse1=0
		tstartPulse2=0
                
		for (( ifl=0; ifl<$nfls; ifl++ )); do
                    filterLength=${flengths[$ifl]}
                    echo "Launching $meth $lib for $mono1EkeV "
                    logf="${instrument}_${detMethod}_${meth}_${filterLength}_${mono1EkeV}${smpStr}${jitterStr}${noiseStr}${bbfbStr}${B0str}.log"
                    command="recon_resol.py --pixName ${instrument} --labelLib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth  --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $filterLength --fdomain $domain --tstartPulse1 $tstartPulse1  --libTmpl LONG --detMethod $detMethod  --filterLength $pulseLength $smpParam --resultsDir gainScale ${jitterParam} ${noiseParam} ${bbfbParam} ${LcParam} ${pBparam} ${sum0Param} ${lagsParam} ${LbTparam} ${ctParam} ${B0param}"
                    nohup python $command >& $logf &
                    echo "Command=python $command "                     
                    echo "Command=python $command >> $logf" 
                done
                sleep 20
	done
fi

if [ $option -eq 2 ]; then
	# To create resolution curves with the coefficients from gain scale curves or surfaces
	# ====================================================
        (( ${preBuffer} > 0 )) && { echo "Warning: You are running options '1' and '2' for 0-pad +  preBuffer"; }
	 
	nSimPulses=2000
        nSimPulses=5000
        #nSimPulses=10000
        for (( ie=0; ie<$nenergies; ie++ )); do
		mono1EkeV=${energies[$ie]}
		mono2EkeV=0
		echo "Runing method comparison (FWHM vs Energy) for energy=$mono1EkeV"
		#tstartPulse1=${tstartPulse1All[$ie]}
		tstartPulse1=0
		tstartPulse2=0
                
		for (( ifl=0; ifl<$nfls; ifl++ )); do
                    filterLength=${flengths[$ifl]}

                    if [ "${calib}" == "2D" ]; then
                        coeffsFile="coeffs_poly2Dfit_pL${filterLength}_${meth}${pulseLength}${ctStr}.dat"
                    fi

                    echo "Using calibration file ${coeffsFile}"
                    echo "Launching $meth $lib for $mono1EkeV  for flength= ${filterLength}"

                    logf="${instrument}_$meth${lib}_${detMethod}${filterLength}_${mono1EkeV}${smpStr}${jitterStr}${noiseStr}${bbfbStr}.log"
                    command="recon_resol.py --pixName ${instrument} --labelLib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $filterLength --fdomain $domain --tstartPulse1 $tstartPulse1  --libTmpl LONG --detMethod $detMethod --filterLength $pulseLength ${smpParam} ${jitterParam} ${noiseParam}  ${bbfbParam} ${LcParam} ${pBparam} ${sum0Param} ${lagsParam}  --coeffsFile ${coeffsFile}"
                    echo "Command=python $command" 
                    #nohup python $command >&  $logf &
                    echo "Command=python $command >> $logf" 
                done    
                #sleep 10
	done
fi

###########   SHORT  FILTER  (AND optionally preBuffer)  ##################################
if  [ $option -eq 3 ]; then
	# =================================================
	# To create gain scale curves & coefficients table
	# =================================================
	# All methods can be run OFLib=yes because length to be used is maximum (pulseLength=largeFilter)
    
	nSimPulses=5000
	#nSimPulses=1
	#energies=(0.2 0.5)
        
	for (( ie=0; ie<$nenergies; ie++ )); do
		mono1EkeV=${energies[$ie]}
		mono2EkeV=0
		echo "Runing gain scale creation (Erecons vs Einput)  for energy=$mono1EkeV"
		#tstartPulse1=${tstartPulse1All[$ie]}
		tstartPulse1=0
		tstartPulse2=0
                
		for (( ifl=0; ifl<$nfls; ifl++ )); do
                    filterLength=${flengths[$ifl]}
                    echo "Launching $meth $lib for $mono1EkeV  and $pBparam"

                    logf="${instrument}_${detMethod}_${meth}_${filterLength}_${mono1EkeV}${smpStr}${jitterStr}${noiseStr}${bbfbStr}${pBstr}${B0str}.log"
                    command="recon_resol.py --pixName ${instrument} --labelLib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth  --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain $domain --tstartPulse1 $tstartPulse1  --libTmpl LONG --detMethod $detMethod --filterLength $filterLength $smpParam --resultsDir gainScale ${jitterParam} ${noiseParam} ${bbfbParam} ${LcParam} ${pBparam}${lagsParam} ${LbTparam} ${ctParam} ${B0param}"
                    nohup python $command >& $logf &
                    echo "Command=python $command >> $logf" 
                done
                #sleep 30
	done
fi

if [ $option -eq 4 ]; then
	# To create resolution curves with the coefficients from gain scale curves
	# =============================================
	nSimPulses=2000
	nSimPulses=5000
        for (( ie=0; ie<$nenergies; ie++ )); do
		mono1EkeV=${energies[$ie]}
		mono2EkeV=0
		echo "Runing method comparison (FWHM vs Energy) for energy=$mono1EkeV"
		#tstartPulse1=${tstartPulse1All[$ie]}
		tstartPulse1=0
		tstartPulse2=0
                
		for (( ifl=0; ifl<$nfls; ifl++ )); do
                    filterLength=${flengths[$ifl]}

                     if [ "${calib}" == "2D" ]; then
                        coeffsFile="coeffs_poly2Dfit_pL${pulseLength}_${meth}${filterLength}${ctStr}.dat"
                    fi
                    echo "Using calibration file ${coeffsFile}"
                    echo "Launching $meth $lib for $mono1EkeV  for flength= ${filterLength}"
                    logf="${instrument}_$meth${lib}_${detMethod}${filterLength}_${mono1EkeV}${smpStr}${jitterStr}${noiseStr}${bbfbStr}${pBstr}${B0str}.log"
                    command="recon_resol.py --pixName ${instrument} --labelLib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain $domain --tstartPulse1 $tstartPulse1  --libTmpl LONG --detMethod $detMethod  --filterLength $filterLength ${smpParam} ${jitterParam} ${noiseParam}  ${bbfbParam} ${LcParam} ${pBparam}${lagsParam} ${LbTparam} ${ctParam} ${B0param} --coeffsFile ${coeffsFile}"
                    nohup python $command >&  $logf &
                    echo "Command=python $command >> $logf" 
                done    
               #sleep 15
	done
fi

