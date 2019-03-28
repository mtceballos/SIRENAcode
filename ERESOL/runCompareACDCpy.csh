#! /bin/csh -f
#
# run different reconstruction methods for all the energies/all record lengths
#
#  Option 1: calculate reconstructed energies and set coefficients table for gain scale curves
#  Option 2: with the Erec vs. Ecal coeffs, calculate FWHM corrected curves

set option=$1
set instrument="LPA75um" #"LPA2shunt"

# ----------------------------------------------------------------------------------------------------------------------
# BBFB and NOISE and Jitter and samprate
set bbfbStr="_bbfb" # "" or "_bbfb"
set jitterStr="_jitter"
set noiseStr="" # "" or"_nonoise"
set smpStr = "" # "" for samprate (156.25kHz)  or "_samprate2" for 78125Hz or "_samprate4" for 39062.5Hz
#set smpStr = ""
# ----------------------------------------------------------------------------------------------------------------------

set jitterParam=""
set bbfbParam=""
set dcmt=100
if($bbfbStr == "_bbfb")  then
    set bbfbParam="--bbfb bbfb"
    set jitterParam="--jitter jitter"
    set dcmt=1
endif
if($jitterStr == "_jitter" && $dcmt > 1)  then 
    set jitterParam="--jitter jitter --decimation $dcmt"
    set jitterStr="_jitter_dcmt$dcmt"
endif
set noiseParam=""
if($noiseStr == "_nonoise")  set noiseParam="--noise nonoise"

#########################################################
set detMethod="AD"
set detMethod="STC"
set samplesUp=0
set samplesDown=0
set nSigmas=0

if($smpStr == "_samprate2")  then
    set separation=20000
    set smpParam="--samprate samprate2"
    set  pulseLength = 4096
    set nSamples=4096
    set largeFilter=4096
    if(${detMethod} == "AD") then
        set nSgms=12 # BE CAREFUL!!! => Check the hard-coded slope criteria 
    else if (${detMethod} == "STC") then
        set nSgms=4 #6
        set samplesUp=2 #2 
        set samplesDown=3 #3
    endif
else # samprate 156250Hz
    set separation=40000
    set smpParam=""
    set  pulseLength = 8192
    set nSamples=8192
    set largeFilter=8192
    if(${detMethod} == "AD") then
        set nSgms=4.4 # BE CAREFUL!!! => Check the hard-coded slope criteria 
    else if (${detMethod} == "STC") then
        set nSgms=3 #5
        set samplesUp=3 #3
        set samplesDown=4 #3
    endif
endif
if($smpStr == "_samprate4")  then
    set separation=10000
    set smpParam="--samprate samprate4"
    set  pulseLength = 2048
    set nSamples=2048
    set largeFilter=2048
    if(${detMethod} == "AD") then
        set nSgms=23.5 # BE CAREFUL!!! => Check the hard-coded slope criteria 
    else if (${detMethod} == "STC") then
        set nSgms=4 #6
        set samplesUp=2 #2 
        set samplesDown=3 #3
    endif
endif

echo "Using det=$detMethod nSgms=$nSgms smplsUp=$samplesUp smplsDown=$samplesDown "
echo "smpParam: $smpParam \njitterParam=$jitterParam \nnoiseParam=$noiseParam \nbbfb=$bbfbParam "


set fixedlib6OF_OPTFILT=1
set fixedlib6OF_OPTFILTNM=0
set fixedlib6OF_I2R=1
set fixedlib6OF_I2RNOL=1
set fixedlib6OF_I2RFITTED=0
set multilibOF_WEIGHTN=0
set multilib_WEIGHTN=0
set multilib_WEIGHT=0

if ($fixedlib6OF_OPTFILT == 1) then
    # Global OPTFILT AC fixedlib6 OFLib=yes
    # -------------------------------------------
    set lib="fixedlib6OF"
    set meth="OPTFILT"
    set nSimPulsesLib=20000 
endif
                
if ($fixedlib6OF_OPTFILTNM == 1) then
    # Global OPTFILT AC fixedlib6 OFLib=yes NoiseMAT
    # -------------------------------------------------------
    set lib="fixedlib6OF"
    set meth="OPTFILTNM"
    set nSimPulsesLib=20000
endif
                
if ($fixedlib6OF_I2R == 1) then
    # Global OPTFILT I2R fixedlib6  --OFLib yes 
    # ---------------------------------------------
    set lib="fixedlib6OF"
    set meth="I2R"
    set nSimPulsesLib=20000
endif
                
if ($fixedlib6OF_I2RNOL == 1) then
    # Global OPTFILT I2RNOL fixedlib6  --OFLib yes 
    # ---------------------------------------------------
    set lib="fixedlib6OF"
    set meth="I2RNOL"
    set nSimPulsesLib=20000
endif

if ($fixedlib6OF_I2RFITTED == 1) then
    # Global OPTFILT I2RFITTED fixedlib6  --OFLib yes 
    # -------------------------------------------------------
    set lib="fixedlib6OF"
    set meth="I2RFITTED"
    set nSimPulsesLib=20000
endif

if ($multilib_WEIGHT == 1) then
    # WEIGHT AC multilib
    #---------------------
    set lib="multilib"
    set meth="WEIGHT"
    set nSimPulsesLib=200000
endif

if ($multilib_WEIGHTN == 1) then
    # WEIGHTN AC multilib
    #---------------------
    set lib="multilib"
    set meth="WEIGHTN"
    set nSimPulsesLib=200000
endif

if ($multilibOF_WEIGHTN == 1) then
    # WEIGHTN AC OF
    #---------------------
    set lib="multilibOF"
    set meth="WEIGHTN"
    set nSimPulsesLib=200000
endif

set energies=(0.2 0.5 1 2 3 4 5 6 7 8)
set energies=(1 2 3 4 5 6 7 8)
set nenergies=$#energies
set lags=(1 1 1 1 1 1 1 1 1 1) # 0.2 and 0.5 keV still difficult to locate
set nSimPulses=20000
set nSimPulses=5000

if ($option == 1) then
	# =================================================
	# To create gain scale curves & coefficients table
	# =================================================
	# All methods can be run OFLib=yes because length to be used is maximum (pulseLength=largeFilter)
	#set tstartPulse1All=(999 999 1000 1000 1000 1000 1000 1000 1000 1000) # default 0 in getEresolCurves
        #set tstartPulse1All=(0 0 0 0 0 0 0 0 0 0) # default 0 in getEresolCurves
	
	foreach ie (`seq 1 $nenergies`)
		set mono1EkeV=${energies[$ie]}
		set mono2EkeV=0
		echo "Runing gain scale creation (Erecons vs Einput)  for energy=$mono1EkeV"
		#set tstartPulse1=${tstartPulse1All[$ie]}
		set tstartPulse1=0
		set tstartPulse2=0
		set lgs=${lags[$ie]}
                
		                
                echo "Launching $meth $lib for $mono1EkeV "
                set logf="${instrument}_$meth${lib}_${detMethod}_${mono1EkeV}${smpStr}${jitterStr}${noiseStr}${bbfbStr}.log"
                set command="recon_resol.py --pixName ${instrument} --labelLib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --tstartPulse1 $tstartPulse1  --libTmpl LONG --detMethod $detMethod --nSgms $nSgms --samplesUp ${samplesUp} --samplesDown ${samplesDown} --filterLength $pulseLength $smpParam --resultsDir gainScale ${jitterParam} ${noiseParam} ${bbfbParam} "
                echo "Command=python $command >& $logf" 
                nohup python $command >>& $logf &
	end
endif

if ($option == 2) then
	# To create resolution curves with the coefficients from gain scale curves
	# =============================================
        foreach ie (`seq 1 $nenergies`)
		set mono1EkeV=${energies[$ie]}
		set mono2EkeV=0
		echo "Runing method comparison (FWHM vs Energy) for energy=$mono1EkeV"
		#set tstartPulse1=${tstartPulse1All[$ie]}
		set tstartPulse1=0
		set tstartPulse2=0
		set lgs=${lags[$ie]}
                
		                
                echo "Launching $meth $lib for $mono1EkeV "
                set logf="${instrument}_$meth${lib}_${detMethod}_${mono1EkeV}${smpStr}${jitterStr}${noiseStr}${bbfbStr}.log"
                set command="recon_resol.py --pixName ${instrument} --labelLib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --tstartPulse1 $tstartPulse1  --libTmpl LONG --detMethod $detMethod --nSgms $nSgms --samplesUp ${samplesUp} --samplesDown ${samplesDown} --filterLength $pulseLength ${smpParam} ${jitterParam} ${noiseParam}  ${bbfbParam} --coeffsFile coeffs_polyfit.dat"
                echo "Command=python $command >& $logf" 
                nohup python $command >>& $logf &
	end
endif


