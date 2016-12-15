#
# Run getEresolTrios.py to create (JSON) data files to plot TRIOS maps
#
#
# OPTIMAL FILTERING, FIXED FILTER 1 keV @ 7 keV (SPIE2016)
#

set array="LPA1shunt"
set lib="fixedlib1"
set monoEkeV=7
set recon="OPTFILT"
set scaleFactor=0.
set samplesUp=0
set nSgms=0
set tstartPulse1=1000
set tstartPulse2=-1
set tstartPulse3=-1
set nsamples=2048
set coeffs="/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresol${array}/coeffs_polyfit.dat"

nohup python getEresolTrios.py --pixType $array --lib $lib --monoEnergy $monoEkeV --reconMethod $recon --filter F0 --nsamples $nsamples --fdomain F --interp DAB --scaleFactor $scaleFactor --samplesUp $samplesUp --nSgms $nSgms --tstartPulse1 $tstartPulse1 --tstartPulse2 $tstartPulse2 --coeffsFile $coeffs >& OPTmultiDABtrios${monoEkeV}_coeffs.log &  
