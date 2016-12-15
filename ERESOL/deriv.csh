fcalc infile=pp.fits outfile=pp.fits clobber=yes clname=DER  expr="seqdiff(ADC)"
fcalc infile=pp.fits outfile=pp.fits clobber=yes clname=UNOS expr=1
fcollen pp.fits UNOS 22000

fcalc infile=pp.fits outfile=pp.fits clobber=yes clname=DER1COL expr=DER
fcollen pp.fits DER1COL 1
fcollen pp.fits DER1COL 22000
fcalc infile=pp.fits outfile=pp.fits clobber=yes clname=DER0COL1 expr="DER-DER1COL"

fcalc infile=pp.fits outfile=pp.fits clobber=yes clname=DER2 expr="UNOS*DER[2]"
fcalc infile=pp.fits outfile=pp.fits clobber=yes clname=DERFINAL expr="DER2+DER0COL1"

fdelcol infile=pp.fits+1 colname=UNOS confirm=no proceed=yes
fdelcol infile=pp.fits+1 colname=DER confirm=no proceed=yes
fdelcol infile=pp.fits+1 colname=DER1COL confirm=no proceed=yes
fdelcol infile=pp.fits+1 colname=DER0COL1 confirm=no proceed=yes
fdelcol infile=pp.fits+1 colname=DER2 confirm=no proceed=yes

