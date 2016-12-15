#!/usr/bin/csh
#
#  Simulate thresome grid of pulses with Philippe's tesconstpileup tool
#
#  Input parameter: tessim (see possibilities/examples below)

if($#argv < 1) then 
    echo "Please complete input command line with all the parameters"
    echo "   source simulTrios.csh tessim"
    exit
endif    
set tessim=$1
set PixType=""
if($tessim =~ *SPA) set PixType="SPA"
if($tessim =~ *LPA) set PixType="LPA"
set XMLrect="/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/CALIBRATION/newFormat/xifu_detector_300us_156kHz_rect.xml"


set sepsStr=(0002 0004 0008 0016 0032 0064 0128 0256 0512 1024 2048 4096 8192)

#set sepsStr=(0001)
#set seps=(1)
set simTime=100
set pixel=1
set samprate=156250.

set pwd=`pwd`
cd $tessim

set monoEkeV=1

@ i = 1
#foreach sep12 ($seps)
foreach sep12 (2 4 8 16 32 64 128)    
    @ j = 8
    #foreach sep23 ($seps)
    foreach sep23 (0256 0512 1024 2048 4096 8192)
        
	set root=sep${sepsStr[$i]}sep${sepsStr[$j]}sam_${simTime}s_${monoEkeV}keV
	echo "-------------------------------------------\n"
	echo "Simulating ${root}.fits\n"
	echo "-------------------------------------------\n"
	tesconstpileup PixImpList=${root}.piximpact XMLFile=$XMLrect tstop=$simTime energy=$monoEkeV pulseDistance=$sep12 pulseDistance2=$sep23 energy2=$monoEkeV energy3=$monoEkeV TriggerSize=10000 clobber=yes
	
	echo "\n##### Runing tessim #########"
	#tessim PixID=$pixel PixImpList=${pwd}/PIXIMPACT/${root}.piximpact Streamfile=${root}.stream tstart=0. tstop=$simTime sample_rate=$samprate clobber=yes PixType=$PixType $m_unknown
	tessim PixID=$pixel PixImpList=${root}.piximpact Streamfile=${root}.stream tstart=0. tstop=$simTime sample_rate=$samprate clobber=yes PixType=$PixType 
	
	echo "\n##### Runing streamtotriggers #########"
	streamtotriggers PixImpList=${root}.piximpact XMLFile=$XMLrect tstart=0. tstop=$simTime Streamfile=${root}.stream TesTriggerFile=${root}.fits TriggerSize=10000 PreBufferSize=1000 pixels=$pixel clobber=yes
	cphead ${root}.stream ${root}.fits
	rm ${root}.stream
	#rm first record and update NETTOT
	fdelrow infile=${root}.fits+1 firstrow=1 nrows=1 confirm="no" proceed=yes
	fkeypar fitsfile=${root}.fits key="NETTOT"
	set nettot=`pget fkeypar value`
	@ nettot = $nettot - 2
	fparkey fitsfile=${root}.fits value=$nettot keyword="NETTOT"
	
	@ j = $j + 1
	
    end #sep23
    @ i = $i + 1
end #sep12

cd $pwd
