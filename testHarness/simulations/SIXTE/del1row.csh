set energies=(0.5 1 2 3 4 5 6 7)
foreach monoEkeV ($energies)
    foreach if (`seq 2 21`)
	fdelrow infile=sim${monoEkeV}.fits.${if}+1 firstrow=1 nrows=1 confirm=no proceed=yes
    end
end
