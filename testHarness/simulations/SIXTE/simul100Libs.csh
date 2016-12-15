#!/bin/csh
#
#  Simulate 100 libs (1 keV) and reconsruct 1 input fits file 100 times (once for each library)
#  to get intermediate filters and check "energy of template"
#
#  Do not initialize HEASOFT!!!!!!!
#  To use it:
#     (heainit;fdump infile=${root1}.stream outfile=${root1}.txt prhead=no showrow=no clobber=yes showcol=no columns="-" rows="-")

unset libFile inFile evtFile noiseDir noiseFile baseline intermFile pwd iter energy
sixteinit

set libFile = "libraryMonoE_PL1024_1keV_tessimLPA1.fits"
set inFile="data10_onerecord.fits" #only one record (only one filter)
#set inFile="mono1_sep3072_pix1_120s_1024.fits"
set evtFile = "evt10_onerecord.fits"
set noiseDir="$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/NOISE/tessimLPA1"
set noiseFile="${noiseDir}/noise1024samples_tessimLPA1_B0_100s_pairscps.fits"
(heainit;fkeypar fitsfile=${noiseFile}+1 keyword="BASELINE")
set baseline=`heainit;pget fkeypar value`
set intermFile="filters.fits"
set pwd=`pwd`
set outfile="energies.dat"
if(-e $outfile) rm $outfile

foreach iter (`seq 1`)
    echo "Running iteration $iter...\n"
    #`source simulLibs.csh LPA1 1024`
    #echo "  Created $libFile"
    cd $pwd/LIBRARIES/tessimLPA1
    
    # Reconstruct only one data file to get filter and normfactor
    if($iter == 1) then 
	echo "tesreconstruction Recordfile=$inFile TesEventFile=$evtFile Rcmethod="SIRENA" PulseLength=1024 LibraryFile=$libFile scaleFactor=0.005 samplesUp=2 nSgms=20 crtLib=0 mode=1 NoiseFile=$noiseFile FilterDomain="F" FilterMethod="F0" calibLQ=1 b_cF=1 c_cF=0 clobber=yes intermediate=1 filterFile=$intermFile baseline=$baseline EnergyMethod=NOLAGS OFStrategy=FREE PixelType=LPA1"
	(heainit;tesreconstruction Recordfile=$inFile TesEventFile=$evtFile Rcmethod="SIRENA" PulseLength=1024 LibraryFile=$libFile scaleFactor=0.005 samplesUp=2 nSgms=20 crtLib=0 mode=1 NoiseFile=$noiseFile FilterDomain="F" FilterMethod="F0" calibLQ=1 b_cF=1 c_cF=0 clobber=yes intermediate=1 detectFile=$intermFile baseline=$baseline EnergyMethod=NOLAGS OFStrategy=FREE PixelType=LPA1)
	echo "   First file reconstructed"
	
	# Dump first filter and first normalization scaleFactor
	#set intermFile="${intermFile}[FILTER]"
	(heainit;fdump ${intermFile}+3 outfile="nrmfctr.txt" columns="NRMFCTR" rows=1 prhead=no showcol=no showunit=no showrow=no clobber=yes)
	(heainit;fdump ${intermFile}+3 outfile="filter1024.txt" columns="OPTIMALF" rows=1 prhead=no showcol=no showunit=no showrow=no clobber=yes)
	echo "   Filter extracted"
    endif
    #(heainit;fdump $libFile outfile="template.txt" columns="PULSEB0" rows=1 prhead=no showcol=no showunit=no showrow=no clobber=yes)
    #echo "   Template extracted"
    
    # Run R script to calculate energy of template
    #set energy=`../../check.R`
    #echo "ENERGY ($iter): $energy ">>$outfile
    #rm $libFile
    #rm $evtFile
    #rm $intermFile
    cd $pwd
    #exit
end
#rm $inFile
cd $pwd
# 