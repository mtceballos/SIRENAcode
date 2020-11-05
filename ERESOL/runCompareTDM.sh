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
instrument="LPA2.5a" #"LPA2shunt"
domain="T"

# ----------------------------------------------------------------------------------------------------------------------
bbfbStr="_fll" # "" or "_bbfb" or "_fll"
noiseStr="" # "" or"_nonoise"
smpStr="" # "" for samprate (156.25kHz)  or "_samprate2" for 78125Hz or "_samprate4" for 39062.5Hz
smpStr=""
pBstr=""
preBuffer=0 # 

Ifit=-23000
IfitStr="_Ifit_m23000"


LbT=0.64E-2 # 1000 samples to average baseline or "0" if baseline is not to be subtracted
LbTparam="--LbT ${LbT}"

bbfParam=""
if [ "$bbfbStr"  ==  "_fll" ];  then
    bbfbParam="--bbfb fll"
elif [ "$bbfbStr"  ==  "_bbfb" ];  then
    bbfbParam="--bbfb bbfb"
fi

calib="1D" # "1D" for gainScale curve or "2D" for "gain scale surface"
 # energies to be calibrated
if [ "${calib}" == "1D" ]; then
    energies=(0.2 0.5 1 2 3 4 5 6 7 8 9 10 11 12)
elif [ "${calib}" == "2D" ]; then
    energies=(7)
fi

# ----------------------------------------------------------------------------------------------------------------------

pBparam=""
noiseParam=""
IfitParam=""

if [ ${preBuffer} -gt 0 ]; then
    pBstr="_pB${preBuffer}"
    pBparam="--preBuffer ${preBuffer}" 
fi
if [ "$noiseStr"  ==  "_nonoise" ];  then
    noiseParam="--noise nonoise"
fi

if [ "${calib}" == "1D" ]; then
    coeffsFile="coeffs_polyfit_methods"$bbfbStr".dat"  
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
fixedlib6OF_OPTFILT=0    # includes 0-padding with SUM/=0 (no action over SUM(filter))
fixedlib6OF_I2R=0
fixedlib6OF_I2RFITTED=1

if  [ $fixedlib6OF_OPTFILT -eq 1 ]; then
    # Global OPTFILT AC fixedlib6 OFLib=yes
    # -------------------------------------------
    lib="fixedlib6OF"
    meth="OPTFILT"
    nSimPulsesLib=20000 
    flengths=(8192 4096 2048 1024 512 256 128 32 16 8)
    #flengths=(8192)    
elif [ $fixedlib6OF_I2R -eq 1 ]; then
    
    # Global OPTFILT I2R fixedlib6  --OFLib yes 
    # ---------------------------------------------
    lib="fixedlib6OF"
    meth="I2R"
    nSimPulsesLib=20000
    flengths=(8192 4096 2048 1024 512 256 128 32 16 8)
    #flengths=(8192)
    
elif [ $fixedlib6OF_I2RFITTED -eq 1 ]; then
    
    # Global OPTFILT I2RFITTED fixedlib6  --OFLib yes 
    # -------------------------------------------------------
    lib="fixedlib6OF"
    meth="I2RFITTED"
    nSimPulsesLib=20000
    flengths=(8192 4096 2048 1024 512 256 128 32 16 8)
    flengths=(8192)
    IfitParam="--Ifit ${Ifit}"
fi

echo "Using det=$detMethod nSgms=$nSgms smplsUp=$samplesUp smplsDown=$samplesDown "
echo "(smpParam: $smpParam) (noiseParam=$noiseParam)  (pBparam=$pBparam)  (IfitParam=$IfitParam)"
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
                    logf="${instrument}_${detMethod}_${meth}_${filterLength}_${mono1EkeV}${smpStr}${noiseStr}${bbfbStr}${pBStr}.log"
                    command="recon_resol.py --pixName ${instrument} --labelLib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth  --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $filterLength --fdomain $domain --tstartPulse1 $tstartPulse1  --libTmpl LONG --detMethod $detMethod  --filterLength $pulseLength $smpParam --resultsDir gainScale  ${noiseParam}  ${bbfbParam} ${pBparam} ${LbTparam} ${IfitParam}"
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
                        coeffsFile="coeffs_poly2Dfit_pL${filterLength}_${meth}${pulseLength}.dat"
                    fi

                    echo "Using calibration file ${coeffsFile}"
                    echo "Launching $meth $lib for $mono1EkeV  for flength= ${filterLength}"
                    
                    logf="${instrument}_${detMethod}_${meth}_${filterLength}_${mono1EkeV}${smpStr}${noiseStr}${bbfbStr}${pBStr}.log"
                    command="recon_resol.py --pixName ${instrument} --labelLib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $filterLength --fdomain $domain --tstartPulse1 $tstartPulse1 --detMethod $detMethod --filterLength $pulseLength ${smpParam} ${noiseParam} ${bbfbParam} ${pBparam} ${IfitParam} --coeffsFile ${coeffsFile}"
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
                    echo "Launching $meth $lib for $mono1EkeV  and $pBparam and $IfitParam"

                    logf="${instrument}_${detMethod}_${meth}_${filterLength}_${mono1EkeV}${smpStr}${noiseStr}${bbfbStr}${pBStr}${IfitStr}.log"
                    command="recon_resol.py --pixName ${instrument} --labelLib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth  --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain $domain --tstartPulse1 $tstartPulse1 --detMethod $detMethod --filterLength $filterLength $smpParam --resultsDir gainScale ${noiseParam} ${bbfbParam} ${pBparam} ${LbTparam} ${IfitParam}"
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
                        coeffsFile="coeffs_poly2Dfit_pL${pulseLength}_${meth}${filterLength}.dat"
                    fi
                    echo "Using calibration file ${coeffsFile}"
                    echo "Launching $meth $lib for $mono1EkeV  for flength= ${filterLength}"
                    logf="${instrument}_${detMethod}_${meth}_${filterLength}_${mono1EkeV}${smpStr}${noiseStr}${bbfbStr}${pBStr}.log"
                    command="recon_resol.py --pixName ${instrument} --labelLib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain $domain --tstartPulse1 $tstartPulse1  --detMethod $detMethod  --filterLength $filterLength ${smpParam} ${noiseParam}  ${bbfbParam}  ${pBparam} ${LbTparam} ${IfitParam}  --coeffsFile ${coeffsFile}"
                    nohup python $command >&  $logf &
                    echo "Command=python $command >> $logf" 
                done    
               #sleep 15
	done
fi

