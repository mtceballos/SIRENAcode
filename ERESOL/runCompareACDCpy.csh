#
# run different reconstruction methods for all the energies/all record lengths
#
#  Option 1: calculate reconstructed energies and set coefficients table for gain scale curves
#  Option 2: with the Erec vs. Ecal coeffs, calculate FWHM corrected curves
#  Option 3: With the Erecons vs. Ecalc coeffs, calculate FWHM corrected curves for different record lengths
#  Option 4: With the Erecons vs. Ecalc coeffs, calculate FWHM/Ebias corrected curves for different separations @ 7keV

set option=$1

if ($option == 1) then
	# =================================================
	# To create gain scale curves & coefficients table
	# =================================================
	#set energies=(0.2 0.5 1 2 3 4 5 6 7 8 9 10 11)
	#set energiesLib=(0.2 0.5 1 2 3 4 5 6 7 8 9 10 11 12 13)
	set energies=(0.2 0.5 1 2 3 4 5 6 7 8)
	#set energiesLib=(0.2 0.5 1 2 3 4 5 6 7 8 9)
	#@ nL1 = $#energiesLib - 1
#	set energies=(0.2 0.5 1)
	set nenergies=$#energies
	set nSamples=2048
	set pulseLength=2048
	set scaleFactor=0.
	set nSimPulses=20000
	# only detect for 0.2 and 1 (problematic energies in simulated files)
	set tstartPulse1All=(0 1000 0 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000) # default 0 in getEresolCurves
	set tstartPulse2All=(0 21000 0 21000 21000 21000 21000 21000 21000 21000 21000 21000 21000) 
	set tstartPulse3All=(0 0 0 0 0 0 0 0 0 0 0 0 0) # default in getEresolCurves
        set nSgmsAll=(7 0 12 0 0 0 0 0 0 0 0 0 0)
        
	foreach ie (`seq 1 $nenergies`)
	#foreach monoEkeV ($energies)
	#       # Intervals multilib
	#	set i = 1
	#	while ($i < $#energiesLib)
	#		@ i2 = $i + 1
	#		# OPTFILT AC duallib DAB
	#   		echo "Launching OPTFILT multilib (${energiesLib[$i]}keV-${energiesLib[$i2]}keV) DAB for $monoEkeV"
	#		nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib${energiesLib[$i]}keV${energiesLib[$i2]}keV --monoEnergy $monoEkeV --reconMethod OPTFILT --filter F0 --nsamples $nSamples --fdomain F --interp DAB >& OPTdual${energiesLib[$i]}-${energiesLib[$i2]}DAB$monoEkeV.log  
	#		@ i++
	#    	end
		set monoEkeV=${energies[$ie]}
		echo "Runing method comparison (FWHM vs Energy) for energy=$monoEkeV"
		set tstartPulse1=${tstartPulse1All[$ie]}
		set tstartPulse2=${tstartPulse2All[$ie]}
		set nSgms=${nSgmsAll[$ie]}
		################ ON-THE-FLY ##########################################
		# Global OPTFILT AC multilib DAB
		# --------------------------------
	    	#echo "Launching OPTFILT multilib DAB for $monoEkeV"
                #set nSimPulsesLib=20000
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod OPTFILT --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& OPTmultiDAB$monoEkeV.log &
 		# Global OPTFILT AC fixedlib1 
 		# --------------------------------
		#echo "Launching OPTFILT fixedlib1 for $monoEkeV"
		#set nSimPulsesLib=20000
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib fixedlib1 --monoEnergy $monoEkeV --reconMethod OPTFILT --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& OPTfixedlib1$monoEkeV.log &

	 	# Global OPTFILT I2R multilib DAB
	 	# --------------------------------
    	 	echo "Launching OPTFILT I2R multilib DAB for $monoEkeV"
    	 	set nSimPulsesLib=20000
		nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod I2R --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RmultiDAB$monoEkeV.log &
 		# Global OPTFILT I2R fixedlib1 
 		# --------------------------------
    	 	echo "Launching OPTFILT I2R fixedlib1 for $monoEkeV"
    	 	set nSimPulsesLib=20000
		nohup python getEresolCurves.py --pixType LPA1shunt --lib fixedlib1 --monoEnergy $monoEkeV --reconMethod I2R --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2Rfixedlib1$monoEkeV.log &

                # AVOID 0.2 keV FITS FOR I2RBALL (NOT RELIABLE - BAD PULSE INITIAL SAMPLE -> BAD ENERGY -> BAD GAIN SCALE
		#if ( $monoEkeV != 0.2) then
                    # Global OPTFILT I2RALL multilib DAB
                    # ------------------------------------
                #    echo "Launching OPTFILT I2RALL multilib DAB for $monoEkeV"
                #    set nSimPulsesLib=20000
                #    nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod I2RALL --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RALLmultiDAB$monoEkeV.log &
                    # Global OPTFILT I2RALL fixedlib1 
                    # --------------------------------
                #    echo "Launching OPTFILT I2RALL fixedlib1 for $monoEkeV"
                #    set nSimPulsesLib=20000
                #    nohup python getEresolCurves.py --pixType LPA1shunt --lib fixedlib1 --monoEnergy $monoEkeV --reconMethod I2RALL --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RALLfixedlib1$monoEkeV.log &
                #endif

 		# Global OPTFILT I2RNOL multilib DAB
 		# ------------------------------------
    	 	#echo "Launching OPTFILT I2RNOL multilib DAB for $monoEkeV"
    	 	#set nSimPulsesLib=20000
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod I2RNOL --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --pulseLength $pulseLength  --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RNOLmultiDAB$monoEkeV.log  &
 		# Global OPTFILT I2RNOL fixedlib1 
 		# --------------------------------
    	 	#echo "Launching OPTFILT I2RNOL fixedlib1 for $monoEkeV"
    	 	#set nSimPulsesLib=20000
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib fixedlib1 --monoEnergy $monoEkeV --reconMethod I2RNOL --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RNOLfixedlib1$monoEkeV.log &

		# Global OPTFILT I2RFITTED multilib DAB
 		# ------------------------------------
    	 	#echo "Launching OPTFILT I2RFITTED multilib DAB for $monoEkeV"
    	 	#set nSimPulsesLib=20000
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod I2RFITTED --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RFITTEDmultiDAB$monoEkeV.log  &
 		# Global OPTFILT I2RFITTED fixedlib1 
 		# --------------------------------
    	 	#echo "Launching OPTFILT I2RFITTED fixedlib1 for $monoEkeV"
    	 	#set nSimPulsesLib=20000
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib fixedlib1 --monoEnergy $monoEkeV --reconMethod I2RFITTED --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RFITTEDfixedlib1$monoEkeV.log 

    		# WEIGHT AC multilib
		#---------------------
                #set nSimPulsesLib=200000
	    	#echo "Launching WEIGHT multilib for $monoEkeV"
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod WEIGHT --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& WEIGHT${monoEkeV}${nSimPulsesLib}p.log &

    		# WEIGHTN AC multilib
		#---------------------
		#set nSimPulsesLib=200000
    		#echo "Launching WEIGHTN multilib for $monoEkeV"
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod WEIGHTN --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& WEIGHTN${monoEkeV}${nSimPulsesLib}.log 

		
    		# WEIGHTN AC multilib  TEST TEST TEST TEST TEST
		#---------------------
		set nSimPulsesLib=200000
    		echo "Launching WEIGHTN multilibOF for $monoEkeV"
		nohup python getEresolCurves.py --pixType LPA1shunt --lib multilibOF --monoEnergy $monoEkeV --reconMethod WEIGHTN --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& WEIGHTNOF${monoEkeV}${nSimPulsesLib}.log &

	end
endif

if ($option == 2) then
	# To create resolution curves (no need to run onthefly & OFLib=yes because length is maximum for pulses (=OFLib)
	# ================================================================================================================
	#set energies=(0.2 0.4 0.5 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 4.4 4.6 4.8 5 5.2\
	#	5.4 5.6 5.8 6 6.2 6.4 6.6 6.8 7 7.2 7.4 7.6 7.8 8 8.2 8.4 8.6 8.8 9 9.2 9.4 9.6 9.8 10 10.2 10.4 10.6 10.8 11)
	set energies=(0.2 0.5 1 2 3 4 5 6 7 8)
	set nenergies=$#energies
        set tstartPulse1All=(0 1000 0 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000) # default 0 in getEresolCurves
	set tstartPulse2All=(0 -1 0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1) #calculate inside getEresolCurves using separation
	set tstartPulse3All=(0 0 0 0 0 0 0 0 0 0 0 0 0) # default in getEresolCurves
	set nSgmsAll=(7 0 12 0 0 0 0 0 0 0 0 0 0)
        set samplesUp=0 # does not mind in production mode (it will be removed afterwards)
        set scaleFactor=0.
        set nSamples=2048
	set pulseLength=2048
        set nSimPulses=20000
        set nSimPulsesLib=20000
	set use="PRIM"
        #set energies=(7)
        #set tstartPulse1All=(1000)
        #set tstartPulse2All=(-1) 
        #set nSgmsAll=(0)
	foreach ie (`seq 1 $#energies`)
                set monoEkeV=${energies[$ie]}
                if($monoEkeV == 7) then
                    continue                 # !!!WARNING
                endif   
		set tstartPulse1=${tstartPulse1All[$ie]}
		set tstartPulse2=${tstartPulse2All[$ie]}
		set nSgms=${nSgmsAll[$ie]}
		# Global OPTFILT AC multilib DAB
		# --------------------------------
		#set lib="multilib"
		#set meth="OPTFILT"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &
		
 		# Global OPTFILT AC fixedlib1 
 		# --------------------------------
 		#set lib="fixedlib1"
		#set meth="OPTFILT"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &

		# Global OPTFILT AC fixedlib1OF
 		# --------------------------------
 		#set lib="fixedlib1OF"
		#set meth="OPTFILT"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &
		
	 	# Global OPTFILT I2R multilib DAB
	 	# --------------------------------
	 	#set lib="multilib"
		#set meth="I2R"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &
 		# Global OPTFILT I2R fixedlib1 
 		# --------------------------------
 		#set lib="fixedlib1"
		#set meth="I2R"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &
                # Global OPTFILT I2R fixedlib1OF
 		# --------------------------------
 		#set lib="fixedlib1OF"
		#set meth="I2R"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &
		#if ( $monoEkeV != 0.2) then
                    # Global OPTFILT I2RALL multilib DAB
                    # ------------------------------------
                    #set lib="multilib"
                    #set meth="I2RALL"
                    #set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    #set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    #nohup $comma >&$logf &
                    # Global OPTFILT I2RALL fixedlib1 
                    # --------------------------------
                    #set lib="fixedlib1"
                    #set meth="I2RALL"
                    #set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    #set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    #nohup $comma >&$logf &
                    # Global OPTFILT I2RALL fixedlib1OF 
                    # --------------------------------
                    #set lib="fixedlib1OF"
                    #set meth="I2RALL"
                    #set logf="$meth${lib}_${monoEkeV}_coeffs.log"
                    #set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
                    #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
                    #nohup $comma >&$logf &
                #endif
 		# Global OPTFILT I2RNOL multilib DAB
 		# ------------------------------------
 		#set lib="multilib"
		#set meth="I2RNOL"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &
                # Global OPTFILT I2RNOL fixedlib1
 		# ------------------------------------
 		#set lib="fixedlib1"
		#set meth="I2RNOL"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &
                # Global OPTFILT I2RNOL fixedlib1OF
 		# ------------------------------------
 		#set lib="fixedlib1OF"
		#set meth="I2RNOL"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &
		
		# Global OPTFILT I2RFITTED multilib 
 		# ------------------------------------
 		#set lib="multilib"
		#set meth="I2RFITTED"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &
                # Global OPTFILT I2RFITTED fixedlib1
 		# ------------------------------------
 		#set lib="fixedlib1"
		#set meth="I2RFITTED"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &
                # Global OPTFILT I2RFITTED fixedlib1OF
 		# ------------------------------------
 		#set lib="fixedlib1OF"
		#set meth="I2RFITTED"
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf 
		
    		# WEIGHT AC multilib
		#---------------------
		#set lib="multilib"
		#set meth="WEIGHT"
		#set nSimPulsesLib=200000
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf &

    		# WEIGHTN AC multilib
		#---------------------
		#set lib="multilib"
		#set meth="WEIGHTN"
		#set nSimPulsesLib=200000
		#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		#nohup $comma >&$logf 
		
                # WEIGHTN AC multilib OF
		#---------------------
		set lib="multilibOF"
		set meth="WEIGHTN"
		set nSimPulsesLib=200000
		set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	    	set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	    	echo "Launching $meth $lib for $monoEkeV with:\n $comma"
		nohup $comma >&$logf 
	end
endif

if ($option == 3) then
    # To calculate resolutions for different record lengths (with highly separated pulses)
    echo "Calculate resolutions for different record lengths (with highly separated pulses)"
    set scaleFactor=0.
    set monoEkeV=7
    set separation=20000
    set tstartPulse1=1000
    set tstartPulse2=21000
    set nSgms=(0)
    set samplesUp=0
    set nSamples=2048
    set nSimPulses=20000
    set nSimPulsesLib=20000
    set recordLens=(1024 512 256 128 64 32)
    #set recordLens=(128)
    set use="PRIM"
    foreach rl ($recordLens)
    
        # Global OPTFILT AC multilib DAB
	# --------------------------------
	#set lib="multilib"
	#set meth="OPTFILT"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &
		
 	# Global OPTFILT AC fixedlib1 
 	# --------------------------------
 	#set lib="fixedlib1"
	#set meth="OPTFILT"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
	#echo "$logf"
    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &
	
	# Global OPTFILT AC fixedlib1OF
 	# --------------------------------
 	#set lib="fixedlib1OF"
	#set meth="OPTFILT"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &
	
	# Global I2R AC multilib DAB
	# --------------------------------
	#set lib="multilib"
	#set meth="I2R"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &
		
 	# Global I2R fixedlib1 
 	# --------------------------------
 	#set lib="fixedlib1"
	#set meth="I2R"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &
	
	# Global OPTFILT I2R fixedlib1OF
 	# --------------------------------
 	#set lib="fixedlib1OF"
	#set meth="I2R"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &

        # Global I2RNOL multilib DAB
	# --------------------------------
	#set lib="multilib"
	#set meth="I2RNOL"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &
		
 	# Global I2RNOL fixedlib1 
 	# --------------------------------
 	#set lib="fixedlib1"
	#set meth="I2RNOL"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &
	
	# Global I2RNOL fixedlib1OF
 	# --------------------------------
 	#set lib="fixedlib1OF"
	#set meth="I2RNOL"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &

        # Global I2RFITTED multilib DAB
	# --------------------------------
	#set lib="multilib"
	#set meth="I2RFITTED"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &
		
 	# Global I2RFITTED fixedlib1 
 	# --------------------------------
 	#set lib="fixedlib1"
	#set meth="I2RFITTED"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &
	
	# Global I2RFITTED fixedlib1OF
 	# --------------------------------
 	#set lib="fixedlib1OF"
	#set meth="I2RFITTED"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf 

	# Global WEIGHT multilib
	# --------------------------------
	set nSimPulsesLib=200000
	set lib="multilib"
	set meth="WEIGHT"
	set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	nohup $comma >& $logf &
		
 	# Global WEIGHTN multilib 
 	# --------------------------------
 	set nSimPulsesLib=200000
 	set lib="multilib"
	set meth="WEIGHTN"
	set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	nohup $comma >& $logf 
	
	# Global WEIGHTN OF
 	# --------------------------------
 	#set nSimPulsesLib=200000
 	#set lib="multilibOF"
	#set meth="WEIGHTN"
	#set logf="$meth${lib}_${monoEkeV}_rl${rl}_coeffs.log"
    	#set comma="python getEresolCurves.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength ${rl} --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with ${rl}:\n $comma"
	#nohup $comma >& $logf &

    end	
    
endif

if ($option == 4) then
        set samplesUp=0 # does not mind in production mode (it will be removed afterwards)
        set scaleFactor=0.
        set nSamples=2048
	set pulseLength=2048
        set nSimPulses=20000
        set nSimPulsesLib=20000
	set use="PRIM"
        set monoEkeV=7
        set tstartPulse1=1000
        set tstartPulse2=-1 
        set nSgms=0
	# Global OPTFILT AC multilib DAB
	# --------------------------------
	#set lib="multilib"
	#set meth="OPTFILT"
	#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
    	#set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
    	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	#nohup $comma >&$logf &
		
        # Global OPTFILT AC fixedlib1 
 	# --------------------------------
 	#set lib="fixedlib1"
	#set meth="OPTFILT"
	#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	#set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	#nohup $comma >&$logf &

	# Global OPTFILT AC fixedlib1OF
 	# --------------------------------
 	#set lib="fixedlib1OF"
	#set meth="OPTFILT"
	#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	#set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	#nohup $comma >&$logf &
		
        # Global OPTFILT I2R multilib DAB
	# --------------------------------
	#set lib="multilib"
	#set meth="I2R"
	#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	#set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	#nohup $comma >&$logf &
 	
 	# Global OPTFILT I2R fixedlib1 
 	# --------------------------------
 	#set lib="fixedlib1"
	#set meth="I2R"
	#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	#set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	#nohup $comma >&$logf &
        
        # Global OPTFILT I2R fixedlib1OF
        # --------------------------------
 	#set lib="fixedlib1OF"
	#set meth="I2R"
	#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	#set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	#nohup $comma >&$logf &
	
	#if ( $monoEkeV != 0.2) then
            # Global OPTFILT I2RALL multilib DAB
            # ------------------------------------
            #set lib="multilib"
            #set meth="I2RALL"
            #set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            #set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            #nohup $comma >&$logf &
            # Global OPTFILT I2RALL fixedlib1 
            # --------------------------------
            #set lib="fixedlib1"
            #set meth="I2RALL"
            #set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            #set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            #nohup $comma >&$logf &
            # Global OPTFILT I2RALL fixedlib1OF 
            # --------------------------------
            #set lib="fixedlib1OF"
            #set meth="I2RALL"
            #set logf="$meth${lib}_${monoEkeV}_coeffs.log"
            #set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
            #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
            #nohup $comma >&$logf &
        #endif
 	# Global OPTFILT I2RNOL multilib DAB
 	# ------------------------------------
 	#set lib="multilib"
	#set meth="I2RNOL"
	#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	#set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	#nohup $comma >&$logf &
        # Global OPTFILT I2RNOL fixedlib1
 	# ------------------------------------
 	#set lib="fixedlib1"
	#set meth="I2RNOL"
	#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	#set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	#nohup $comma >&$logf &
        # Global OPTFILT I2RNOL fixedlib1OF
 	# ------------------------------------
 	#set lib="fixedlib1OF"
	#set meth="I2RNOL"
	#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	#set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	#nohup $comma >&$logf &
		
	# Global OPTFILT I2RFITTED multilib 
 	# ------------------------------------
 	#set lib="multilib"
	#set meth="I2RFITTED"
	#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	#set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	#nohup $comma >&$logf &
        # Global OPTFILT I2RFITTED fixedlib1
 	# ------------------------------------
 	#set lib="fixedlib1"
	#set meth="I2RFITTED"
	#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	#set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	#nohup $comma >&$logf &
        # Global OPTFILT I2RFITTED fixedlib1OF
        # ------------------------------------
 	#set lib="fixedlib1OF"
	#set meth="I2RFITTED"
	#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	#set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	#nohup $comma >&$logf &
		
    	# WEIGHT AC multilib
	#---------------------
	#set lib="multilib"
	#set meth="WEIGHT"
	#set nSimPulsesLib=200000
	#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	#set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
        #echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	#nohup $comma >&$logf &

    	# WEIGHTN AC multilib
	#---------------------
	#set lib="multilib"
	#set meth="WEIGHTN"
	#set nSimPulsesLib=200000
	#set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	#set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	#echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	#nohup $comma >&$logf &
		
        # WEIGHTN AC multilib OF
	#---------------------
	set lib="multilibOF"
	set meth="WEIGHTN"
	set nSimPulsesLib=200000
	set logf="$meth${lib}_${monoEkeV}_coeffs.log"
	set comma="python getEresolCurves_manySeps.py --pixType LPA1shunt --lib $lib --monoEnergy $monoEkeV --reconMethod $meth --filter F0 --nsamples $nSamples --nSimPulses $nSimPulses --nSimPulsesLib $nSimPulsesLib --pulseLength $pulseLength --fdomain F --scaleFactor $scaleFactor --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit_${nSamples}_${use}.dat"
	echo "Launching $meth $lib for $monoEkeV with:\n $comma"
	nohup $comma >&$logf &
endif
