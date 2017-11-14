#
# run different reconstruction methods for all the energies/all record lengths
#
#  Option 1: calculate reconstructed energies and set coefficients table for gain scale curves
#  Option 2: with the Erec vs. Ecal coeffs, calculate FWHM corrected curves
#  Option 3: With the Erecons vs. Ecalc coeffs, calculate FWHM corrected curves for different record lengths
#  Option 4: With the Erecons vs. Ecalc coeffs, calculate FWHM/Ebias corrected curves for different separations @ 7keV

# #### LPA2 hay que detectar para E=0.2, 0.5 keV

set option=$1
set instrument="LPA2shunt"
set separation=40000
if ($instrument == "LPA1shunt") then
   set  pulseLength=2048
   set nSamples=2048
else if ($instrument == "LPA2shunt") then
    set  pulseLength = 4096
    set nSamples=4096
    set largeFilter=32768
endif
set smpStr = "" # samprate = 156.25kHz 
#set smpStr = "samprate2" # samprate = 156.25/2kHz 

set fixedlib1_OPTFILT=0
set fixedlib1OF_OPTFILT=1
set fixedlib1OF_OPTFILTNM=0
set fixedlib1_I2R=0
set fixedlib1OF_I2R=0
set fixedlib1_I2RNOL=0
set fixedlib1OF_I2RNOL=0
set fixedlib1_I2RFITTED=0
set fixedlib1OF_I2RFITTED=0
set multilibOF_WEIGHTN=0
set multilib_WEIGHTN=0
set multilib_WEIGHT=0

if ($option == 1) then
	# =================================================
	# To create gain scale curves & coefficients table
	# =================================================
	# All methods can be run OFLib=yes because length to be used is maximum (pulseLength=largeFilter)
	
	set energies=(0.2 0.5 1 2 3 4 5 6 7 8)
        set energies=(0.2)
	set nenergies=$#energies
	#set nSamples=$largeFilter
	#set pulseLength=$largeFilter
	set scaleFactor=0.
	set nSimPulses=20000
	# For LPA2: do detection for 0.2, 0.5, 1 and 2 (problem with threshold for these energies in simulated files; 
	# only in libraries; corrected in singles/pairs simulations)
	set lags=(1 1 1 1 1 1 1 1 1 1) # 0.2 and 0.5 keV still difficult to locate
	set tstartPulse1All=(999 999 1000 1000 1000 1000 1000 1000 1000 1000) # default 0 in getEresolCurves
        set tstartPulse1All=(0 0 0 0 0 0 0 0 0 0) # default 0 in getEresolCurves
	# if LPA1, calibration files are pairs (in LPA2 are single pulses: no tstartPulse2All nor tstartPulse3All)
	#set nSgmsAll=(4.6 4.6 12 12 0 0 0 0 0 0)
        #set nSgmsAll=(0 0 0 0 0 0 0 0 0 0)
        #set nSgmsAll=(11.3 11.3 11.3 11.3 11.3 11.3 11.3 11.3 11.3 11.3)
        
	foreach ie (`seq 1 $nenergies`)
		set mono1EkeV=${energies[$ie]}
		set mono2EkeV=0
		echo "Runing method comparison (FWHM vs Energy) for energy=$mono1EkeV"
		set tstartPulse1=${tstartPulse1All[$ie]}
		set tstartPulse2=0
		set lgs=${lags[$ie]}
		
		################ ON-THE-FLY ##########################################
		# Global OPTFILT AC multilib DAB
		# --------------------------------
	    	#echo "Launching OPTFILT multilib DAB for $monoEkeV"
                #set nSimPulsesLib=20000
		#nohup python getEresolCurves.py --pixName ${instrument} --lib multilib --monoEnergy $monoEkeV --reconMethod OPTFILT --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& OPTmultiDAB$monoEkeV.log &
 		if ($fixedlib1OF_OPTFILT == 1) then
                    # Global OPTFILT AC fixedlib1 OFLib=yes
                    # -------------------------------------------
                    set lib="fixedlib1OF"
                    set meth="OPTFILT"
                    echo "Launching $meth $lib for $mono1EkeV"
                    set nSimPulsesLib=20000
                    foreach detMethod ("AD" "A1")
                        if($detMethod == "AD") then
                            set nSgms=11.3 # samprate=156.25 kHz
                            set nSamplesUp=0
                            set nSamplesDown=0
                            if($smpStr == "samprate2")   set nSgms=23.5 # samprate=78.125 kHz
                        endif
                        if($detMethod == "A1") then
                            set nSgms=5 # samprate=156.25 kHz
                            set nSamplesUp=3
                            set nSamplesDown=3
                            if($smpStr == "samprate2") then
                                set nSgms=6 # samprate=78.125 kHz
                                set nSamplesUp=2
                                set nSamplesDown=3
                            endif       
                        endif
                        nohup python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy1 $mono1EkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --detMethod=$detMethod --samplesUp=$nSamplesUp --samplesDown=$nSamplesDown  --libTmpl="LONG" --samprate=$smpStr --resultsDir="gainScale">& $logf &
                endif
		if ($fixedlib1OF_OPTFILTNM == 1) then
                    # Global OPTFILT AC fixedlib1 OFLib=yes NoiseMAT
                    # -------------------------------------------------------
                    set lib="fixedlib1OF"
                    set meth="OPTFILTNM"
                    echo "Launching $meth $lib for $monoEkeV"
                    set nSimPulsesLib=20000
                    set logf="${instrument}_$meth${lib}_${monoEkeV}.log"
                    nohup python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& $logf &
                endif
	 	# Global OPTFILT I2R multilib DAB
	 	# --------------------------------
    	 	#echo "Launching OPTFILT I2R multilib DAB for $monoEkeV"
    	 	#set nSimPulsesLib=20000
		#nohup python getEresolCurves.py --pixName ${instrument} --lib multilib --monoEnergy $monoEkeV --reconMethod I2R --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RmultiDAB$monoEkeV.log &
 		if ($fixedlib1OF_I2R == 1) then
                    # Global OPTFILT I2R fixedlib1  --OFLib yes 
                    # ---------------------------------------------
                    set lib="fixedlib1OF"
                    set meth="I2R"
                    echo "Launching $meth $lib for $monoEkeV"
                    set nSimPulsesLib=20000
                    set logf="${instrument}_$meth${lib}_${monoEkeV}.log"
                    nohup python getEresolCurves.py --pixName ${instrument} --lib $lib  --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& $logf &
                endif
                # AVOID 0.2 keV FITS FOR I2RBALL (NOT RELIABLE - BAD PULSE INITIAL SAMPLE -> BAD ENERGY -> BAD GAIN SCALE
		#if ( $monoEkeV != 0.2) then
                    # Global OPTFILT I2RALL multilib DAB
                    # ------------------------------------
                #    echo "Launching OPTFILT I2RALL multilib DAB for $monoEkeV"
                #    set nSimPulsesLib=20000
                #    nohup python getEresolCurves.py --pixName ${instrument} --lib multilib --monoEnergy $monoEkeV --reconMethod I2RALL --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RALLmultiDAB$monoEkeV.log &
                    # Global OPTFILT I2RALL fixedlib1 
                    # --------------------------------
                #    echo "Launching OPTFILT I2RALL fixedlib1 for $monoEkeV"
                #    set nSimPulsesLib=20000
                #    nohup python getEresolCurves.py --pixName ${instrument} --lib fixedlib1 --monoEnergy $monoEkeV --reconMethod I2RALL --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RALLfixedlib1$monoEkeV.log &
                #endif
 		# Global OPTFILT I2RNOL multilib DAB
 		# ------------------------------------
    	 	#echo "Launching OPTFILT I2RNOL multilib DAB for $monoEkeV"
    	 	#set nSimPulsesLib=20000
		#nohup python getEresolCurves.py --pixName ${instrument} --lib multilib --monoEnergy $monoEkeV --reconMethod I2RNOL --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --pulseLength $pulseLength  --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RNOLmultiDAB$monoEkeV.log  &
		if ($fixedlib1OF_I2RNOL == 1) then
                    # Global OPTFILT I2RNOL fixedlib1  --OFLib yes 
                    # ---------------------------------------------------
                    set lib="fixedlib1OF"
                    set meth="I2RNOL"
                    echo "Launching $meth $lib for $monoEkeV"
                    set nSimPulsesLib=20000
                    set logf="${instrument}_$meth${lib}_${monoEkeV}.log"
                    nohup python getEresolCurves.py --pixName ${instrument} --lib $lib  --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& $logf &
                endif
		# Global OPTFILT I2RFITTED multilib DAB
 		# ------------------------------------
    	 	#echo "Launching OPTFILT I2RFITTED multilib DAB for $monoEkeV"
    	 	#set nSimPulsesLib=20000
		#nohup python getEresolCurves.py --pixName ${instrument} --lib multilib --monoEnergy $monoEkeV --reconMethod I2RFITTED --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RFITTEDmultiDAB$monoEkeV.log  &
		if ($fixedlib1OF_I2RFITTED == 1) then
                    # Global OPTFILT I2RFITTED fixedlib1  --OFLib yes 
                    # -------------------------------------------------------
                    set lib="fixedlib1OF"
                    set meth="I2RFITTED"
                    echo "Launching $meth $lib for $monoEkeV"
                    set nSimPulsesLib=20000
                    set logf="${instrument}_$meth${lib}_${monoEkeV}.log"
                    nohup python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& $logf
                endif
                if ($multilib_WEIGHT == 1) then
                    # WEIGHT AC multilib
                    #---------------------
                    set lib="multilib"
                    set meth="WEIGHT"
                    set nSimPulsesLib=200000
                    echo "Launching WEIGHT multilib for $monoEkeV"
                    set logf="${instrument}_$meth${lib}_${monoEkeV}.log"
                    nohup python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& $logf &
                endif
                if ($multilib_WEIGHTN == 1) then
                    # WEIGHTN AC multilib
                    #---------------------
                    set lib="multilib"
                    set meth="WEIGHTN"
                    set nSimPulsesLib=200000
                    echo "Launching WEIGHT multilib for $monoEkeV"
                    set logf="${instrument}_$meth${lib}_${monoEkeV}.log"
                    nohup python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& $logf &
                endif
                if ($multilibOF_WEIGHTN == 1) then
                    # WEIGHTN AC OF
                    #---------------------
                    set lib="multilibOF"
                    set meth="WEIGHTN"
                    set nSimPulsesLib=200000
                    echo "Launching $meth $lib for $monoEkeV"
                    set logf="${instrument}_$meth${lib}_${monoEkeV}.log"
                    nohup python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& $logf
                endif
	end
endif

if ($option == 2) then
	# To create resolution curves (no need to run onthefly & OFLib=yes because length is maximum for pulses (=OFLib)
	# ================================================================================================================
	
	set energies=(0.2 0.5 1 2 3 4 5 6 7 8)
	#set energies=(0.2 0.5)
	set lags=(0 0 0 0 0 0 0 0 0 0)
        set tstartPulse1All=(999 999 1000 1000 1000 1000 1000 1000 1000 1000) # default 0 in getEresolCurves
        #set tstartPulse1All=(0 0 1000 1000 1000 1000 1000 1000 1000 1000) # default 0 in getEresolCurves
	set nSgmsAll=(9 9 0 0 0 0 0 0 0 0)
        set samplesUp=0 # does not mind in production mode (it will be removed afterwards)
        set scaleFactor=0.
        #set energies=(2)
        #set tstartPulse1All=(1000)
	set nenergies=$#energies
        set use="all" # all pulses; not PRIM or SEC
        	foreach ie (`seq 1 $#energies`)
                set nSimPulses=20000
                set nSimPulsesLib=20000
                set monoEkeV=${energies[$ie]}
                set lgs=${lags[$ie]}
                
                #if($monoEkeV == 7) then
                #    continue                 # so as not to rewrite 7 keV files with several separations
                #endif   
		
		set tstartPulse1=${tstartPulse1All[$ie]}
		set tstartPulse2=0
		set nSgms=${nSgmsAll[$ie]}
		
		# Global OPTFILT AC multilib DAB
		# --------------------------------
		#set lib="multilib"
		#set meth="OPTFILT"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &
		if ($fixedlib1_OPTFILT == 1) then
                    # Global OPTFILT AC fixedlib1 
                    # --------------------------------
                    set lib="fixedlib1"
                    set meth="OPTFILT"
                    set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf &
                endif
                if ($fixedlib1OF_OPTFILT == 1) then
                    # Global OPTFILT AC fixedlib1OF
                    # --------------------------------
                    set lib="fixedlib1OF"
                    set meth="OPTFILT"
                    set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf &
                endif
                if ($fixedlib1OF_OPTFILTNM == 1) then
                    # Global OPTFILTNM AC fixedlib1OF
                    # --------------------------------
                    set lib="fixedlib1OF"
                    set meth="OPTFILTNM"
                    set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf &
                endif
                # Global OPTFILT I2R multilib DAB
	 	# --------------------------------
	 	#set lib="multilib"
		#set meth="I2R"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F  --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &
		if ($fixedlib1_I2R == 1) then
                    # Global OPTFILT I2R fixedlib1 
                    # --------------------------------
                    set lib="fixedlib1"
                    set meth="I2R"
                    set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F  --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf &
                endif
                if ($fixedlib1OF_I2R == 1) then
                    # Global OPTFILT I2R fixedlib1OF
                    # --------------------------------
                    set lib="fixedlib1OF"
                    set meth="I2R"
                    set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf &
                endif
		#if ( $monoEkeV != 0.2) then
                    # Global OPTFILT I2RALL multilib DAB
                    # ------------------------------------
                    #set lib="multilib"
                    #set meth="I2RALL"
                    #set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    #set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F  --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    #nohup $comma >&$logf &
                    # Global OPTFILT I2RALL fixedlib1 
                    # --------------------------------
                    #set lib="fixedlib1"
                    #set meth="I2RALL"
                    #set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    #set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F  --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    #nohup $comma >&$logf &
                    # Global OPTFILT I2RALL fixedlib1OF 
                    # --------------------------------
                    #set lib="fixedlib1OF"
                    #set meth="I2RALL"
                    #set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    #set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F  --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    #nohup $comma >&$logf &
                #endif
 		# Global OPTFILT I2RNOL multilib DAB
 		# ------------------------------------
 		#set lib="multilib"
		#set meth="I2RNOL"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F  --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &
		if ($fixedlib1_I2RNOL == 1) then
                    # Global OPTFILT I2RNOL fixedlib1
                    # ------------------------------------
                    set lib="fixedlib1"
                    set meth="I2RNOL"
                    set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F  --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf &
                endif
                if ($fixedlib1OF_I2RNOL == 1) then
                    # Global OPTFILT I2RNOL fixedlib1OF
                    # ------------------------------------
                    set lib="fixedlib1OF"
                    set meth="I2RNOL"
                    set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf &
                endif    
		# Global OPTFILT I2RFITTED multilib 
 		# ------------------------------------
 		#set lib="multilib"
		#set meth="I2RFITTED"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F  --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &
		if ($fixedlib1_I2RFITTED == 1) then
                    # Global OPTFILT I2RFITTED fixedlib1
                    # ------------------------------------
                    set lib="fixedlib1"
                    set meth="I2RFITTED"
                    set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F  --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf &
                endif
                if ($fixedlib1OF_I2RFITTED == 1) then
                    # Global OPTFILT I2RFITTED fixedlib1OF
                    # ------------------------------------
                    set lib="fixedlib1OF"
                    set meth="I2RFITTED"
                    set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F  --lags $lgs --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf &
                endif
                if ($multilib_WEIGHT == 1) then
                    # WEIGHT AC multilib
                    #---------------------
                    set lib="multilib"
                    set meth="WEIGHT"
                    set nSimPulsesLib=200000
                    set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 ${tstartPulse1All[$ie]} --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf 
                endif
                if ($multilib_WEIGHTN == 1) then
                    # WEIGHTN AC multilib
                    #---------------------
                    set lib="multilib"
                    set meth="WEIGHTN"
                    set nSimPulsesLib=200000
                    set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses #$nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 ${tstartPulse1All[$ie]}  --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf 
		endif
		if ($multilibOF_WEIGHTN == 1) then
                    # WEIGHTN AC multilib OF
                    #---------------------
                    set lib="multilibOF"
                    set meth="WEIGHTN"
                    set nSimPulsesLib=200000
                    set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 ${tstartPulse1All[$ie]}  --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf 
                endif
	end
endif

if ($option == 3) then
    # To calculate resolutions for different record lengths 
    echo "Calculate resolutions for different record lengths"
    set scaleFactor=0.
    set monoEkeV=7
    #set monoEkeV=1
    set tstartPulse1=1000
    set nSgms=0
    set samplesUp=0
    set recordLens=(4096 2048 1024 750 512 400 256 200 128 90 64 45 32)
    #set recordLens=(4096)
    set use="all"
    foreach rl ($recordLens)
        echo " ........For record length: $rl"
        set nSimPulses=20000
        set nSimPulsesLib=20000
        # Global OPTFILT AC multilib DAB
	# --------------------------------
	#set lib="multilib"
	#set meth="OPTFILT"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &
        if ($fixedlib1_OPTFILT == 1) then
            # Global OPTFILT AC fixedlib1 
            # --------------------------------
            set lib="fixedlib1"
            set meth="OPTFILT"
            set nSimPulsesLib=20000
            set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
            echo "$logf"
            set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1  --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
            nohup $comma >& $logf &
	endif
	if ($fixedlib1OF_OPTFILT == 1) then
            # Global OPTFILT AC fixedlib1OF
            # --------------------------------
            set lib="fixedlib1OF"
            set meth="OPTFILT"
            set nSimPulsesLib=20000
            set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
            set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
            nohup $comma >& $logf &
	endif
	if ($fixedlib1OF_OPTFILTNM == 1) then
            # Global OPTFILT AC fixedlib1OF NoiseMAT
 	# ----------------------------------------------
            set lib="fixedlib1OF"
            set meth="OPTFILTNM"
            set nSimPulses=20000
            set nSimPulsesLib=20000
            set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
            set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
            nohup $comma >& $logf &
	endif
	# Global I2R AC multilib DAB
	# --------------------------------
	#set lib="multilib"
	#set meth="I2R"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1   --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &
	if ($fixedlib1_I2R == 1) then	
            # Global I2R fixedlib1 
            # --------------------------------
            set lib="fixedlib1"
            set meth="I2R"
            set nSimPulsesLib=20000
            set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
            set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
            nohup $comma >& $logf &
	endif
	if ($fixedlib1OF_I2R == 1) then	
            # Global OPTFILT I2R fixedlib1OF
            # --------------------------------
            set lib="fixedlib1OF"
            set meth="I2R"
            set nSimPulsesLib=20000
            set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
            set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
            nohup $comma >& $logf &
        endif
        # Global I2RNOL multilib DAB
	# --------------------------------
	#set lib="multilib"
	#set meth="I2RNOL"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &
	if ($fixedlib1_I2RNOL == 1) then		
            # Global I2RNOL fixedlib1 
            # --------------------------------
            set lib="fixedlib1"
            set meth="I2RNOL"
            set nSimPulsesLib=20000
            set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
            set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
            nohup $comma >& $logf &
	endif
	if ($fixedlib1OF_I2RNOL == 1) then		
            # Global I2RNOL fixedlib1OF
            # --------------------------------
            set lib="fixedlib1OF"
            set meth="I2RNOL"
            set nSimPulsesLib=20000
            set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
            set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
            nohup $comma >& $logf &
        endif
        # Global I2RFITTED multilib DAB
	# --------------------------------
	#set lib="multilib"
	#set meth="I2RFITTED"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &
	if ($fixedlib1_I2RFITTED == 1) then			
            # Global I2RFITTED fixedlib1 
            # --------------------------------
            set lib="fixedlib1"
            set meth="I2RFITTED"
            set nSimPulsesLib=20000
            set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
            set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
            nohup $comma >& $logf &
	endif
	if ($fixedlib1OF_I2RFITTED == 1) then		
            # Global I2RFITTED fixedlib1OF
            # --------------------------------
            set lib="fixedlib1OF"
            set meth="I2RFITTED"
            set nSimPulsesLib=20000
            set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
            set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
            nohup $comma >& $logf 
        endif
        if ($multilibOF_WEIGHTN == 1) then		
            # Global WEIGHTN OF
            # --------------------------------
            set nSimPulsesLib=200000
            set nSimPulses=5000
            set lib="multilibOF"
            set meth="WEIGHTN"
            set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
            set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
            nohup $comma >& $logf &
        endif        
        if ($multilib_WEIGHTN == 1) then		
            # Global WEIGHTN multilib 
            # --------------------------------
            set nSimPulsesLib=200000
            set nSimPulses=5000
            set lib="multilib"
            set meth="WEIGHTN"
            set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
            set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
            nohup $comma >& $logf &
	endif
        if ($multilib_WEIGHT == 1) then		
            # Global WEIGHT multilib
            # --------------------------------
            set nSimPulsesLib=200000
            set nSimPulses=5000
            set lib="multilib"
            set meth="WEIGHT"
            set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
            set comma="python getEresolCurves.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
            nohup $comma >& $logf 
        endif
    end	
    
endif

if($option == 4) then
    #  Option 4: With the Erecons vs. Ecalc coeffs, calculate FWHM/Ebias corrected curves for different separations @ 7keV
    set samplesUp=0 # does not mind in production mode (it will be removed afterwards)
    set scaleFactor=0.
    set nSimPulses=20000
    set nSimPulsesLib=20000
    set use="all"
    set monoEkeV=7
    set tstartPulse1=1000
    set tstartPulse2=-1 
    set nSgms=0
    set separations=('00023'  '00031'  '00042'  '00056'  '00075'  '00101'  '00136'  '00182'  '00244'  '00328'  '00439'  '00589'  '00791'  '01061'  '01423'  '01908'  '02560'  '03433'  '04605'  '40000')
    
    if ($fixedlib1_OPTFILT == 1) then
            # Global OPTFILT AC fixedlib1 
            # --------------------------------
            set lib="fixedlib1"
            set meth="OPTFILT"
            set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --separations $separations"
            echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            nohup $comma >&$logf &
    endif
    if ($fixedlib1OF_OPTFILT == 1) then
            # Global OPTFILT AC fixedlib1OF
            # --------------------------------
            set lib="fixedlib1OF"
            set meth="OPTFILT"
            set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --separations $separations"
            echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            nohup $comma >&$logf &
    endif	
    if ($fixedlib1OF_OPTFILTNM == 1) then
            # Global OPTFILTNM AC fixedlib1OF
            # --------------------------------
            set lib="fixedlib1OF"
            set meth="OPTFILTNM"
            set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --separations $separations"
            echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            nohup $comma >&$logf &
    endif	
    if ($fixedlib1_I2R == 1) then
            # Global OPTFILT I2R fixedlib1 
            # --------------------------------
            set lib="fixedlib1"
            set meth="I2R"
            set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --separations $separations"
            echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            nohup $comma >&$logf &
    endif
    if ($fixedlib1OF_I2R == 1) then
            # Global OPTFILT I2R fixedlib1OF
            # --------------------------------
            set lib="fixedlib1OF"
            set meth="I2R"
            set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --separations $separations"
            echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            nohup $comma >&$logf &
    endif
    if ($fixedlib1_I2RNOL == 1) then
            # Global OPTFILT I2RNOL fixedlib1
            # ------------------------------------
            set lib="fixedlib1"
            set meth="I2RNOL"
            set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --separations $separations"
            echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            nohup $comma >&$logf &
    endif
    if ($fixedlib1OF_I2RNOL == 1) then
            # Global OPTFILT I2RNOL fixedlib1OF
            # ------------------------------------
            set lib="fixedlib1OF"
            set meth="I2RNOL"
            set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --separations $separations"
            echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            nohup $comma >&$logf &
    endif
    if ($fixedlib1_I2RFITTED == 1) then
            # Global OPTFILT I2RFITTED fixedlib1
            # ------------------------------------
            set lib="fixedlib1"
            set meth="I2RFITTED"
            set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --separations $separations"
            echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            nohup $comma >&$logf &
    endif
    if ($fixedlib1OF_I2RFITTED == 1) then
            # Global OPTFILT I2RFITTED fixedlib1OF
            # ------------------------------------
            set lib="fixedlib1OF"
            set meth="I2RFITTED"
            set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --separations $separations"
            echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            nohup $comma >&$logf &
    endif
    if ($multilib_WEIGHT == 1) then
            # WEIGHT AC multilib
            #---------------------
            set lib="multilib"
            set meth="WEIGHT"
            set nSimPulsesLib=200000
            set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --separations $separations"
            echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            nohup $comma >&$logf 
    endif
    if ($multilib_WEIGHTN == 1) then
            # WEIGHTN AC multilib
            #---------------------
            set lib="multilib"
            set meth="WEIGHTN"
            set nSimPulsesLib=200000
            set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --separations $separations"
            echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            nohup $comma >&$logf 
    endif
    if ($multilibOF_WEIGHTN == 1) then
            # WEIGHTN AC multilib OF
            #---------------------
            set lib="multilibOF"
            set meth="WEIGHTN"
            set nSimPulsesLib=200000
            set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --separations $separations"
            echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            nohup $comma >&$logf 
    endif
endif

if ($option == 5) then
    #  Option 5: With the Erecons vs. Ecalc coeffs, calculate FWHM corrected curves for different separations of two pulses @ 6keV & XkeV
        set samplesUp=0 # does not mind in production mode (it will be removed afterwards)
        set scaleFactor=0.
        set nSimPulses=2000
        set nSimPulsesLib=200000 
	set use="all"
        set monoEkeV=6
        set energies=(0.2 0.5 1 2 3 4 5 6 7)
        set tstartPulse1=0
        set tstartPulse2=0 
        set nSgms=11
        set detectSP=1   # choose whether secondary detection (adjusted derivative) is to be taken into account
	set separations=('00005' '00010'  '00020'  '00045')
	
	foreach ie (`seq 1 $#energies`)
                set mono2EkeV=${energies[$ie]}
                
                if ($fixedlib1OF_OPTFILT == 1) then
                    # Global OPTFILT AC fixedlib1OF
                    # --------------------------------
                    set lib="fixedlib1OF"
                    set meth="OPTFILT"
                    set logf="$meth${lib}_${monoEkeV}_${mono2EkeV}_coeffs.log"
                    set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --detectSP=${detectSP} --separations $separations"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf 
                endif
                if ($fixedlib1OF_OPTFILTNM == 1) then
                    # Global OPTFILTNM AC fixedlib1OF
                    # --------------------------------
                    set lib="fixedlib1OF"
                    set meth="OPTFILTNM"
                    set logf="$meth${lib}_${monoEkeV}_${mono2EkeV}_coeffs.log"
                    set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --detectSP=${detectSP} --separations $separations"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf &
                endif
                if ($fixedlib1OF_I2R == 1) then
                    # Global OPTFILT I2R fixedlib1OF
                    # --------------------------------
                    set lib="fixedlib1OF"
                    set meth="I2R"
                    set logf="$meth${lib}_${monoEkeV}_${mono2EkeV}_coeffs.log"
                    set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat  --detectSP=${detectSP} --separations $separations"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf &
                endif    
                if ($fixedlib1OF_I2RNOL == 1) then
                    # Global OPTFILT I2RNOL fixedlib1OF
                    # ------------------------------------
                    set lib="fixedlib1OF"
                    set meth="I2RNOL"
                    set logf="$meth${lib}_${monoEkeV}_${mono2EkeV}_coeffs.log"
                    set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --detectSP=${detectSP} --separations $separations"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf &
                endif
                if ($fixedlib1OF_I2RFITTED == 1) then
                    # Global OPTFILT I2RFITTED fixedlib1OF
                    # ------------------------------------
                    set lib="fixedlib1OF"
                    set meth="I2RFITTED"
                    set logf="$meth${lib}_${monoEkeV}_${mono2EkeV}_coeffs.log"
                    set comma="python getEresolCurves_manySeps.py --pixName ${instrument} --lib $lib --monoEnergy $monoEkeV --monoEnergy2 $mono2EkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat --detectSP=${detectSP} --separations $separations"
                    echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    nohup $comma >&$logf 
                endif
        end 
endif

