#!/bin/sh

#
# Shell script to add keywords for XIFUSIM
# param1: file
# param2: decimation factor

fits0=${1}+0
fits1=${1}[TESRECORDS]
dcmt=${2}

fparkey value=" " fitsfile=${fits0}  keyword='TELESCOP' add=y
fparkey value=" " fitsfile=${fits0}  keyword='INSTRUME'  add=y
fparkey value=" " fitsfile=${fits0}  keyword='FILTER' add=y
fparkey value=" " fitsfile=${fits0}  keyword='ANCRFILE' add=y
fparkey value=" " fitsfile=${fits0}  keyword='RESPFILE' add=y
fparkey value=0.0 fitsfile=${fits0}  keyword='MJDREF' add=y
fparkey value=0.0 fitsfile=${fits0}  keyword='TSTART' add=y
fparkey value=0.0 fitsfile=${fits0}  keyword='TSTOP' add=y
fparkey value=0.0 fitsfile=${fits0}  keyword='TIMEZERO' add=y

fparkey value=" " fitsfile=${fits1}  keyword='TELESCOP' add=y
fparkey value=" " fitsfile=${fits1}  keyword='INSTRUME'  add=y
fparkey value=" " fitsfile=${fits1}  keyword='FILTER' add=y
fparkey value=" " fitsfile=${fits1}  keyword='ANCRFILE' add=y
fparkey value=" " fitsfile=${fits1}  keyword='RESPFILE' add=y
fparkey value=0.0 fitsfile=${fits1}  keyword='MJDREF' add=y
fparkey value=0.0 fitsfile=${fits1}  keyword='TSTART' add=y
fparkey value=0.0 fitsfile=${fits1}  keyword='TSTOP' add=y
fparkey value=0.0 fitsfile=${fits1}  keyword='TIMEZERO' add=y
fparkey value=0 fitsfile=${fits1}  keyword='FIRSTPIX' add=y
fparkey value=0 fitsfile=${fits1}  keyword='LASTPIX' add=y
fparkey value=0 fitsfile=${fits1}  keyword='NPIX' add=y
fparkey value=-1 fitsfile=${fits1}  keyword='MONOEN' add=y
fparkey value=0 fitsfile=${fits1}  keyword='NESTOT' add=y

fkeypar fitsfile=${fits1} keyword='NAXIS2'
NETTOT=`pget fkeypar value`
fparkey value=${NETTOT} fitsfile=${fits1}  keyword='NETTOT' add=y

fkeypar fitsfile=${fits1} keyword='DELTA_T'
DELTA_T=`pget fkeypar value`
DELTAT=`python -c "print ($DELTA_T * $dcmt)"`
#echo "DELTAT=$DELTAT"

# Be careful here with fixed/variable-length columns...
fkeypar fitsfile=${fits1} keyword='TFORM2'
colsinfo=`pget fkeypar value`
#colsinfo=`flcol ${fits1}`
#echo "colsinfo=$colsinfo"
if echo "$colsinfo" | grep -q "("; then
    TRIGGSZ=`echo $colsinfo| awk -F"[()]" '{print $2}'`
else
    TRIGGSZ=`echo $colsinfo| awk -F"[D]" '{print $1}'`
    TRIGGSZ=`echo $TRIGGSZ|tr "'" " "`
    TRIGGSZ=`echo $TRIGGSZ|tr ' 1' "1"`
fi
#echo "TRIGGSZ=${TRIGGSZ}"

fparkey value=${TRIGGSZ} fitsfile=${fits0}  keyword='TRIGGSZ' add=y
fparkey value=${DELTAT}  fitsfile=${fits0}  keyword='DELTAT'  add=y
fparkey value=${TRIGGSZ} fitsfile=${fits1}  keyword='TRIGGSZ' add=y
fparkey value=${DELTAT}  fitsfile=${fits1}  keyword='DELTAT'  add=y

#fparkey value="RECORDS" fitsfile=${fits1}  keyword='EXTNAME' add=y

