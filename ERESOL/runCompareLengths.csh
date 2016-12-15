#
# run different reconstruction methods for a given energy  and several record Length
#
#  Option 1: calculate reconstructed energies and set coefficients table for gain scale curves
#  Option 2: with the Erec vs. Ecal coeffs, calculate FWHM corrected curves

set option=$1

if ($option == 1) then
	# To create gain scale curves & coefficients table
	# =================================================
	set energies=(0.2 0.5 1 2 3 4 5 6 7 8 9 10 11)
	set energiesLib=(0.2 0.5 1 2 3 4 5 6 7 8 9 10 11 12 13)
	set nenergies=$#energies
	set scaleFactor=0.
	# only detect for 0.2 and 1
	set tstartPulse1All=(0 1000 0 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000) # default 0 in getEresolCurves
	set tstartPulse2All=(0 21000 0 21000 21000 21000 21000 21000 21000 21000 21000 21000 21000) 
	set tstartPulse3All=(0 0 0 0 0 0 0 0 0 0 0 0 0) # default in getEresolCurves
        set nSgmsAll=(7 0 12 0 0 0 0 0 0 0 0 0 0)
        set samplesUp=2 # does not mind in production mode (it will be removed afterwards)
	@ nL1 = $#energiesLib - 1
	foreach ie (`seq 1 $nenergies`)
	#foreach monoEkeV ($energies)
	#       # Intervals multilib
	#	set i = 1
	#	while ($i < $#energiesLib)
	#		@ i2 = $i + 1
	#		# OPTFILT AC duallib DAB
	#   		echo "Launching OPTFILT multilib (${energiesLib[$i]}keV-${energiesLib[$i2]}keV) DAB for $monoEkeV"
	#		nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib${energiesLib[$i]}keV${energiesLib[$i2]}keV --monoEnergy $monoEkeV --reconMethod OPTFILT --filter F0 --nsamples 2048 --fdomain F --interp DAB >& OPTdual${energiesLib[$i]}-${energiesLib[$i2]}DAB$monoEkeV.log  
	#		@ i++
	#    	end
		set monoEkeV=${energies[$ie]}
		set tstartPulse1=${tstartPulse1All[$ie]}
		set tstartPulse2=${tstartPulse2All[$ie]}
		set nSgms=${nSgmsAll[$ie]}
		# Global OPTFILT AC multilib DAB
		# --------------------------------
	    	#echo "Launching OPTFILT multilib DAB for $monoEkeV"
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod OPTFILT --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& OPTmultiDAB$monoEkeV.log &
 		# Global OPTFILT AC fixedlib1 
 		# --------------------------------
		#echo "Launching OPTFILT fixedlib1 for $monoEkeV"
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib fixedlib1 --monoEnergy $monoEkeV --reconMethod OPTFILT --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& OPTfixedlib1$monoEkeV.log &

	 	# Global OPTFILT I2R multilib DAB
	 	# --------------------------------
    	 	#echo "Launching OPTFILT I2R multilib DAB for $monoEkeV"
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod I2R --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RmultiDAB$monoEkeV.log &
 		# Global OPTFILT I2R fixedlib1 
 		# --------------------------------
    	 	#echo "Launching OPTFILT I2R fixedlib1 for $monoEkeV"
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib fixedlib1 --monoEnergy $monoEkeV --reconMethod I2R --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2Rfixedlib1$monoEkeV.log 

	 	# Global OPTFILT I2RBALL multilib DAB
	 	# ------------------------------------
    	 	#echo "Launching OPTFILT I2RBALL multilib DAB for $monoEkeV"
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod I2RBALL --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RBALLmultiDAB$monoEkeV.log &
 		# Global OPTFILT I2RBALL fixedlib1 
 		# --------------------------------
    	 	#echo "Launching OPTFILT I2RBALL fixedlib1 for $monoEkeV"
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib fixedlib1 --monoEnergy $monoEkeV --reconMethod I2RBALL --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RBALLfixedlib1$monoEkeV.log &

 		# Global OPTFILT I2RBNOL multilib DAB
 		# ------------------------------------
    	 	#echo "Launching OPTFILT I2RBNOL multilib DAB for $monoEkeV"
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod I2RBNOL --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RBNOLmultiDAB$monoEkeV.log  &
 		# Global OPTFILT I2RBNOL fixedlib1 
 		# --------------------------------
    	 	#echo "Launching OPTFILT I2RBNOL fixedlib1 for $monoEkeV"
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib fixedlib1 --monoEnergy $monoEkeV --reconMethod I2RBNOL --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& I2RBNOLfixedlib1$monoEkeV.log 

    		# WEIGHT AC multilib
		#---------------------
	    	#echo "Launching WEIGHT multilib for $monoEkeV"
		#nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod WEIGHT --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& WEIGHT$monoEkeV.log 

    		# WEIGHTN AC multilib
		#---------------------
    		echo "Launching WEIGHTN multilib for $monoEkeV"
		nohup python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod WEIGHTN --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 >& WEIGHTN$monoEkeV.log 

	end
endif

if ($option == 2) then
	# To create resolution curves
	# =================================
	#set energies=(0.2 0.4 0.5 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4 4.2 4.4 4.6 4.8 5 5.2\
	#	5.4 5.6 5.8 6 6.2 6.4 6.6 6.8 7 7.2 7.4 7.6 7.8 8 8.2 8.4 8.6 8.8 9 9.2 9.4 9.6 9.8 10 10.2 10.4 10.6 10.8 11)
	set energies=(0.2 0.5 1 2 3 4 5 6 7 8 9 10 11)
	set nenergies=$#energies
        set tstartPulse1All=(0 1000 0 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000) # default 0 in getEresolCurves
	set tstartPulse2All=(0 21000 0 21000 21000 21000 21000 21000 21000 21000 21000 21000 21000) 
	set tstartPulse3All=(0 0 0 0 0 0 0 0 0 0 0 0 0) # default in getEresolCurves
	set nSgmsAll=(7 0 12 0 0 0 0 0 0 0 0 0 0)
        set samplesUp=0 # does not mind in production mode (it will be removed afterwards)
        set scaleFactor=0.
        
        set energies=(7)
        set tstartPulse1All=(1000)
        set tstartPulse2All=(-1) #calculate inside getEresolCurves using separation
        set tstartPulse3All=(0)
        set nSgms=(0)
	#set energies=(7)
	#foreach monoEkeV ($energies)
	foreach ie (`seq 1 $#energies`)
                set monoEkeV=${energies[$ie]}
		set tstartPulse1=${tstartPulse1All[$ie]}
		set tstartPulse2=${tstartPulse2All[$ie]}
		set nSgms=${nSgmsAll[$ie]}
		
		# Global OPTFILT AC multilib DAB
		# --------------------------------
	    	set comm="python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod OPTFILT --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit.dat"
	    	echo "Launching OPTFILT multilib DAB for $monoEkeV with:\n $comm"
		nohup $comm >& OPTmultiDAB${monoEkeV}_coeffs.log &
		
 		# Global OPTFILT AC fixedlib1 
 		# --------------------------------
 		set comm="python getEresolCurves.py --pixType LPA1shunt --lib fixedlib1 --monoEnergy $monoEkeV --reconMethod OPTFILT --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit.dat"
		echo "Launching OPTFILT fixedlib1 for $monoEkeV with:\n $comm"
		nohup $comm >& OPTfixedlib1${monoEkeV}_coeffs.log &

	 	# Global OPTFILT I2R multilib DAB
	 	# --------------------------------
	 	set comm="python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod I2R --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit.dat"
    	 	echo "Launching OPTFILT I2R multilib DAB for $monoEkeV with:\n $comm"
		nohup $comm >& I2RmultiDAB${monoEkeV}_coeffs.log &
 		# Global OPTFILT I2R fixedlib1 
 		# --------------------------------
 		set comm="python getEresolCurves.py --pixType LPA1shunt --lib fixedlib1 --monoEnergy $monoEkeV --reconMethod I2R --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit.dat"
    	 	echo "Launching OPTFILT I2R fixedlib1 for $monoEkeV with:\n $comm"
		nohup $comm >& I2Rfixedlib1${monoEkeV}_coeffs.log & 

	 	# Global OPTFILT I2RBALL multilib DAB
	 	# ------------------------------------
	 	#set comm="python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod I2RBALL --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit.dat"
    	 	#echo "Launching OPTFILT I2RBALL multilib DAB for $monoEkeV with:\n $comm"
		#nohup $comm >& I2RBALLmultiDAB${monoEkeV}_coeffs.log &
 		# Global OPTFILT I2RBALL fixedlib1 
 		# --------------------------------
 		#set comm="python getEresolCurves.py --pixType LPA1shunt --lib fixedlib1 --monoEnergy $monoEkeV --reconMethod I2RBALL --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit.dat"
    	 	#echo "Launching OPTFILT I2RBALL fixedlib1 for $monoEkeV with:\n $comm"
		#nohup $comm >& I2RBALLfixedlib1${monoEkeV}_coeffs.log &

 		# Global OPTFILT I2RBNOL multilib DAB
 		# ------------------------------------
 		set comm="python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod I2RBNOL --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit.dat"
    	 	echo "Launching OPTFILT I2RBNOL multilib DAB for $monoEkeV with:\n $comm"
		nohup $comm >& I2RBNOLmultiDAB${monoEkeV}_coeffs.log  &
 		# Global OPTFILT I2RBNOL fixedlib1 
 		# --------------------------------
 		set comm="python getEresolCurves.py --pixType LPA1shunt --lib fixedlib1 --monoEnergy $monoEkeV --reconMethod I2RBNOL --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit.dat"
    	 	echo "Launching OPTFILT I2RBNOL fixedlib1 for $monoEkeV with:\n $comm"
		nohup $comm >& I2RBNOLfixedlib1${monoEkeV}_coeffs.log &

    		# WEIGHT AC multilib
		#---------------------
		set comm="python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod WEIGHT --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit.dat"
	    	echo "Launching WEIGHT multilib for $monoEkeV with:\n $comm"
		nohup $comm >& WEIGHT${monoEkeV}_coeffs.log &

    		# WEIGHTN AC multilib
		#---------------------
		set comm="python getEresolCurves.py --pixType LPA1shunt --lib multilib --monoEnergy $monoEkeV --reconMethod WEIGHTN --filter F0 --nsamples 2048 --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile coeffs_polyfit.dat"
    		echo "Launching WEIGHTN multilib for $monoEkeV with:\n $comm"
		nohup $comm >& WEIGHTN${monoEkeV}_coeffs.log 

	end	
endif

