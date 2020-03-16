set energyStr=${1} # ex. "0.2keV" or "0.2keV_0.2keV"
set filter=${2}
echo "Creating png files for $energyStr"
#foreach i ( `seq 10 35` )
set i=0
foreach infile  ( `ls scProd_${energyStr}_filter${filter}_*pdf` )
    set ii=`printf "%02d" $i`
    echo "Proccessing $infile for i=${ii}"
    #set infile="scProd_${energyStr}_filter${filter}_${i}.pdf"
    set outfile="scProd${energyStr}_filter${filter}_${ii}"
    pdftoppm $infile $outfile -png
    @ i++
end
set moviefile="movie${energyStr}_filter${filter}.mp4"
echo "Creating movie for $energyStr and filter $filter"
convert -delay 0.5x3 scProd${energyStr}_filter${filter}_*.png $moviefile
