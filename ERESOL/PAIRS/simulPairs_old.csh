#!/usr/bin/csh
#
#  Simulate pairs of pulses with Philippe's tesconstpileup tool
#
#  Input parameter: tessim (see possibilities/examples below)

alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'

if($#argv < 1) then 
    echo "Please complete input command line with all the parameters"
    echo "   source simulPairs.csh tessim"
    exit
endif    
set tessim=$1
set PixType=""
if($tessim =~ *SPA) then
    set PixType="file:/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/tespixels.fits[SPA]"
    set sepsStr=(0004 0008 0016 0032 0064 0128 0256 0512 1024 2048 ) # 4096 )# 8192)
endif
if($tessim =~ *LPA1) then
    set PixType="file:/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/tespixels.fits[LPA1]"
    #set sepsStr=(0004 0008 0016 0032 0064 0128 0256 0512 1024 2048 4096 8192)
    set sepsStr=(00004 00005 00007 00010 00013 00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287)
endif
if($tessim =~ *LPA2) then
    set PixType="file:/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/tespixels.fits[LPA2]"
    set sepsStr=(00004 00005 00007 00010 00013 00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287 11115 14910 20000)
endif
if($tessim =~ *LPA3) then
    set PixType="file:/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/tespixels.fits[LPA3]"
    set sepsStr=(00004 00005 00007 00010 00013 00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287 11115 14910 20000)
    #set sepsStr=(20000)
endif

set simTime=100
set pixel=1
set samprate=156250.
set pulseLength=2048 # just to define triggerSize
set PreBufferSize=1000 #   "

set pwd=`pwd`
cd $tessim

foreach monoEkeV (1)

  set XMLrect="/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/CALIBRATION/newFormat/xifu_detector_300us_156kHz_rect.xml";
  
    foreach sep12 ($sepsStr)
	MATH sep = ($sep12 * 1.)
	MATH triggerSize = ($PreBufferSize+$sep+$pulseLength+10)
	
	set root=sep${sep12}sam_${simTime}s_${monoEkeV}keV
	echo "-------------------------------------------\n"
	echo "Simulating ${root}.fits\n"
	echo "-------------------------------------------\n"
	if(-e ${root}.fits) continue
	if(! -e ${pwd}/PIXIMPACT/${root}.piximpact) then 
	    tesconstpileup PixImpList=${pwd}/PIXIMPACT/${root}.piximpact XMLFile=$XMLrect tstop=$simTime energy=$monoEkeV pulseDistance=$sep TriggerSize=$triggerSize clobber=yes
	endif
	#runtes PixImpList=${root}.piximpact XMLFile=$XMLfile pixels=1 tstart=0.0 tstop=$simTime WriteStreamFile=no Streamfile=pp.stream Reconstruct=no WriteRecordFile=yes TesTriggerfile=${root}.fits clobber=yes seed=-1
	echo "\n##### Runing tessim #########"

	tessim PixID=$pixel PixImpList=${pwd}/PIXIMPACT/${root}.piximpact Streamfile=${root}.stream tstart=0. tstop=$simTime sample_rate=$samprate clobber=yes PixType="$PixType"

	echo "\n##### Runing streamtotriggers #########"
	streamtotriggers PixImpList=${pwd}/PIXIMPACT/${root}.piximpact XMLFile=$XMLrect tstart=0. tstop=$simTime Streamfile=${root}.stream TesTriggerFile=${root}.fits TriggerSize=$triggerSize PreBufferSize=$PreBufferSize pixels=$pixel clobber=yes
	cphead ${root}.stream ${root}.fits
	rm ${root}.stream
	#rm first (and LAST) record and update NETTOT (first starts high and last can be cut)
	fkeypar fitsfile=${root}.fits key="NAXIS2"
 	set nrows=`pget fkeypar value`
 	fdelrow infile=${root}.fits+1 firstrow=$nrows nrows=1 confirm="no" proceed=yes
 	fdelrow infile=${root}.fits+1 firstrow=1 nrows=1 confirm="no" proceed=yes
 	fkeypar fitsfile=${root}.fits key="NAXIS2"
 	set nrows=`pget fkeypar value`
 	MATH nettot = ($nrows*2)
 	fparkey fitsfile=${root}.fits value=$nettot keyword="NETTOT"
		
    end #sep
end #monoen
cd $pwd
#end
