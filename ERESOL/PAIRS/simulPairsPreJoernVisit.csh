#!/usr/bin/csh
#
#  Simulate pairs of pulses with Philippe's tesconstpileup tool
#
#  Input parameters: array (SPA, LPA1, LPA2, LPA3)
#                    monoEkeV (keV)

# Global definitions and system settings
alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'
#set dateNow=`date +%s`

if($#argv < 1) then 
    echo "Please complete input command line with all the parameters"
    echo "   source simulPairs.csh array monoEkeV ACDC"
    exit
endif    
set array=$1
set tessim="tessim"$1
set monoEkeV=$2
set ACDC=`echo $3 | tr '[:upper:]' '[:lower:]'` # AC or DC

set XMLrect="/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/CALIBRATION/newFormat/xifu_detector_300us_156kHz_rect.xml"
set nSimPulses=20000
set pixel=1
set samprate=156250. # Hz
set PreBufferSize=1000 # 128?
set pwd=`pwd`

# Define tmp dir
###################
set tmpdir="simulPairs${array}_${monoEkeV}"
mkdir -p /tmp/${tmpdir}/pfiles
setenv PFILES "/tmp/${tmpdir}/pfiles;${HEADAS}/syspfiles:${SIXTE}/share/simput/pfiles:$SIXTE/share/sixte/pfiles"
echo "PFILES=$PFILES"
setenv HEADASNOQUERY
setenv HEADASPROMPT /dev/null

unset sepsStr
unset PixType
#set PixType="file:/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/tespixels.fits[$array]"
set PixType="file:/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/newpixels.fits[${array}${ACDC}]"

if($array == "SPA") then
    
    #set sepsStr=(0004 0008 0016 0032 0064 0128 0256 0512 1024 2048 ) # 4096 )# 8192)
    #set sepsStr=(00002 00003 00004 00005 00007 00008 00011 00013 00017 00022 00028 00035 00045 00057 00072 00092 00116 00148  00188 00238 00303 00384 00488 00620 00787 01000 01270 01613 02048)
    set sepsStr=(00004 00006 00009 00014 00021 00031 00047 00071 00108 00163 00246 00371 00560 00845 01276 01926 02907 04389)
    set sepsStr=(20000)

else if($array == "LPA1") then
   
    #set sepsStr=(0004 0008 0016 0032 0064 0128 0256 0512 1024 2048 4096 8192)
    #set sepsStr=(00004 00005 00007 00010 00013 00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287)
    set sepsStr=(00004 00006 00009 00014 00021 00031 00047 00071 00108 00163 00246 00371 00560 00845 01276 01926 02907 04389 06625 10000)    

else if($array == "LPA2") then
    
    #set sepsStr=(00004 00005 00007 00010 00013 00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287 11115 14910 20000)
    set sepsStr=(00004 00006 00009 00014 00021 00031 00047 00071 00108 00163 00246 00371 00560 00845 01276 01926 02907 04389 06625 10000)    
    
else if($array == "LPA3") then
    
    #set sepsStr=(00004 00005 00007 00010 00013 00017 00023 00031 00042 00056 00075 00101 00136 00182 00244 00328 00439 00589 00791 01061 01423 01908 02560 03433 04605 06178 08287 11115 14910 20000)
    set sepsStr=(00004 00006 00009 00014 00021 00031 00047 00071 00108 00163 00246 00371 00560 00845 01276 01926 02907 04389 06625 10000)    
    #set sepsStr=(20000)
    
endif


cd $tessim

foreach sepA ($sepsStr)
	MATH sep12 = ($sepA * 1.)
	#MATH triggerSize = ($PreBufferSize+$sep12+$sep12+1000) -> for pretrigg tessim
	# for new trigger-tessim (no streamtotriggers)
	MATH triggerSizeTC = ($PreBufferSize+$sep12+$sep12+$PreBufferSize+1000)
	MATH triggerSizeTS = ($PreBufferSize+$sep12)
	
	# calculate sim time to have at least nSimPulses pulses:          
	#   simTime= (nSimPulses/2)recs * (triggerSize sam/rec)/(samprate sam/s)
	MATH simTime =  (0.5*$nSimPulses*$triggerSizeTC/$samprate)
	set simTime=`printf "%.0f" $simTime`
	set root0=sep${sepA}sam_${simTime}s_${monoEkeV}keV
	set root=sep${sepA}sam_${nSimPulses}p_${monoEkeV}keV
	set pixFile=${pwd}/PIXIMPACT/${root0}_trSz${triggerSizeTC}.piximpact
	set streamFile=${root0}.stream
	set fitsFile=${root}.fits
	echo "-------------------------------------------\n"
	echo "Simulating $fitsFile\n"
	echo "-------------------------------------------\n"
	
	if(! -e $pixFile) then
	    echo "\n##### Runing tesconstpileup #########"
	    tesconstpileup PixImpList=$pixFile XMLFile=$XMLrect tstop=$simTime energy=$monoEkeV pulseDistance=$sep12 TriggerSize=$triggerSizeTC clobber=yes
	endif
	
	#continue  # to simulate only piximpact files
	
	#echo "\n##### Runing tessim #########"
	#echo "tessim PixType=$PixType PixID=$pixel PixImpList=$pixFile Streamfile=${root}.stream tstart=0. tstop=$simTime triggertype=stream triggersize=$triggerSize prebuffer=$PreBufferSize clobber=yes" 
	    
# 	tessim PixType="$PixType" PixID=$pixel PixImpList=$pixFile Streamfile=$streamFile  
#	tstart=0. tstop=$simTime clobber=yes
	echo "\n##### Runing tessim #########"
	echo "tessim PixID=$pixel PixImpList=$pixFile Streamfile=$fitsFile tstart=0. tstop=$simTime sample_rate=$samprate triggerSize=$triggerSizeTS preBuffer=$PreBufferSize triggertype='movavg:5:1.1:0' PixType='$PixType'"
	tessim PixID=$pixel PixImpList=$pixFile Streamfile=$fitsFile tstart=0. tstop=$simTime sample_rate=$samprate triggerSize=$triggerSizeTS preBuffer=$PreBufferSize triggertype='movavg:5:1.1:0' PixType="$PixType"
	#   continue
	#echo "\n##### Runing streamtotriggers #########"
	#streamtotriggers PixImpList=$pixFile XMLFile=$XMLrect tstart=0. tstop=$simTime Streamfile=$streamFile TesTriggerFile=$fitsFile TriggerSize=$triggerSize PreBufferSize=$PreBufferSize pixels=$pixel clobber=yes
	#cphead $streamFile+1 $fitsFile+1
	#rm $streamFile
	
	#rm first (and LAST) record and update NETTOT
	fkeypar fitsfile=$fitsFile key="NAXIS2"
	set nrows=`pget fkeypar value`
	fdelrow infile=$fitsFile+1 firstrow=$nrows nrows=1 confirm="no" proceed=yes
	fdelrow infile=$fitsFile+1 firstrow=1 nrows=1 confirm="no" proceed=yes
	fkeypar fitsfile=$fitsFile key="NAXIS2"
	set nrows=`pget fkeypar value`
	MATH nettot = ($nrows*2)
	fparkey fitsfile=$fitsFile value=$nettot keyword="NETTOT"
		    
end #sep12

cd $pwd
rm -rf /tmp/${tmpdir}

