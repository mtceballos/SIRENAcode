set energies=(0.2 0.5 1 2 3 4 5 6 7 8 9)
foreach monoEkeV ($energies)
    set outfile="LIBRARIES/tessimLPA2shunt/mono${monoEkeV}_sep40000_pix1_200000p_4096.fits"
    set firstfile="sim${monoEkeV}.fits.1"
    # rm first/last row of first simulated file
    cp ${firstfile} pp.tmp
    fkeypar fitsfile=pp.tmp key=NAXIS2
    set nrows=`pget fkeypar value`
    #       rm last row
    fdelrow infile=pp.tmp+1 firstrow=${nrows} nrows=1 confirm=no proceed=yes 
    #       rm first row
    fdelrow infile=pp.tmp+1 firstrow=1 nrows=1 confirm=no proceed=yes
    # cp first file to output file
    cp pp.tmp $outfile

    foreach if (`seq 2 21`)
	set infile="sim${monoEkeV}.fits.${if}"
	cp ${infile} pp.tmp
	# rm first/last row of simulated file
	fkeypar fitsfile=pp.tmp key=NAXIS2
	set nrows=`pget fkeypar value`
	# rm last row
	fdelrow infile=pp.tmp+1 firstrow=${nrows} nrows=1 confirm=no proceed=yes 
	# rm first row
	fdelrow infile=pp.tmp+1 firstrow=1 nrows=1 confirm=no proceed=yes

	echo "Merging $infile to energy=$monoEkeV"
	./tabmergeADC pp.tmp+1 ${outfile}+1
    end
    fkeypar fitsfile=${outfile} key=NAXIS2
    set ntotalrows=`pget fkeypar value`
    fparkey fitsfile=${outfile} key=NAXIS2 value=${ntotalrows} 
    rm pp.tmp
end
