"""
APPEND FITS files extensions having variable array columns

Extensions to be appended should have the same columns in both files
Header of first file[0] will be copied to the output file if input parameter "head=True"
Header of first file[extension] will be copied to the output file

python appendVarFits.py --fitsList "file1.fits[1]" "file2.fits[3]" --outfits <output.fits>

"""

from __future__ import print_function
from astropy.io.fits import getheader, update, getdata, Header
from astropy.io import fits
import numpy as np

def appendVarFits(fitsList, head, outfits):

    """
        :type fitsList: list
        :param fitsList: (separated by blanks) list of FITS files with HDU number to be merged
        :type outfits: str
        :param outfits: final fits file with both input extensions appended
        :type head: bool
        :param head: copy the HDU=0 header  (if not one of the selected) of the initial file into the final one?
    """

    print("Appending ", fitsList)
    fits1 = fitsList[0][0:fitsList[0].index("[")]
    hdu1 = int(fitsList[0][fitsList[0].index("[")+1: fitsList[0].index("]")])

    # Read headers from first file: HDU1 (and HDU0) will be (optinally) added to final file
    print("Reading initial values from file:", fits1)
    f1 = fits.open(fits1, memmap=True)
    if head:
        header0 = getheader(fits1, 0)
    header1 = getheader(fits1, hdu1)
    ncols = len(f1[hdu1].columns)
    f1.close()

    COLs = list()
    # run along list of columns in first file
    for icol in range(ncols):
        COLdataList = []
        # for each column add data in every input file
        # one row => numpy array
        # all rows => list of numpy arrays
        for infile in fitsList:
            filename = infile[0:infile.index("[")]
            hdu = int(infile[infile.index("[")+1: infile.index("]")])
            f = fits.open(filename, memmap=True)
            nrows = f[hdu].data.shape[0]
            # consider possibility of column being scaled
            # see https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/f_user/node50.html
            # empirically found that:
            # if valfits > TZERO => astropy sets (valfits-2*TZERO)
            tscalKey = "TSCAL" + str(icol+1)
            tzeroKey = "TZERO" + str(icol+1)
            try:
                TSCAL = f[hdu].header[tscalKey]
                TZERO = f[hdu].header[tzeroKey]
            except:
                TSCAL = 1
                TZERO = 0
            #print(tscalKey,"/",tzeroKey,"=",TSCAL,",",TZERO)
            for irow in range(nrows):
                dataForRow = f[hdu].data.field(icol)[irow]  #int16
                # negInds = np.where(dataForRow < 0.)
                # npTzero = np.zeros(dataForRow.size)
                # npTzero[negInds] = 2*TZERO
                #dataForRowScaled = dataForRow + npTzero  # float64
                dataForRowScaled = dataForRow
                COLdataList.append(dataForRowScaled)
                # if icol == 1 and irow == 0:  #ADC
                    # print("data read by astropy around pulse=", dataForRow[1010:1030])
                    # print("data after scale=",dataForRowScaled[1010:1030])

            colname = f[hdu].columns[icol].name
            colformat = f[hdu].columns[icol].format
            # print("colformat antes=", colformat)
            # colformat = colformat.replace("PI", "PJ")
            # print("colformat despues=", colformat)
            f.close()
        # convert list of numpy arrays into a numpy array of numpy arrays
        # as required by  fits.Column
        COLarray = np.array(COLdataList)
        # Add new column to list of columns for output HDU
        COLs.append(fits.Column(name=colname, format=colformat, array=COLarray))

    # write new fits file
    hduTotal = fits.BinTableHDU.from_columns(COLs)
    hduTotal.writeto(outfits)
    # update header of first extension of output file
    outdata1 = getdata(outfits, ext=1)
    update(outfits, data=outdata1, header=header1, ext=1)

    # update TSCAL, TZERO values for columns
    out = fits.open(outfits, 'update')
    for icol in range(ncols):
        tscalKey = "TSCAL" + str(icol + 1)
        tzeroKey = "TZERO" + str(icol + 1)
        try:
            out[1].header[tscalKey] = header1[tscalKey]
            out[1].header[tzeroKey] = header1[tzeroKey]
        except:
            continue
    print("data written by astropy around pulse=", out[1].data.field(1)[0][1010:1030])
    if head:
        # if required, add header of HDU0 from first file to utput file
        # some weird things made impossible the use of same "update" call before
        header0new = getheader(outfits, ext=0)
        # check for keywords in HDU0 of first file and add them to
        # output file if not present
        for key in header0.keys():
            try:
                k0n = header0new[key]
            except:
                out[0].header[key] = header0[key]
    out.close()

        # print(repr(getheader(outfits, ext=0)))

#
# --------- Read input parameters and append FITS files -----------------
#
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Append 2 FITS files with variable array columns', prog='appendVarFits')

    parser.add_argument('--fitsList', nargs='*', help='Python list of FITS files and extension numbers')
    parser.add_argument('--head', help='copy the HDU=0 header of the initial file into the final one?',
                        action="store_true")
    parser.add_argument('--outfits', help='Output FITS file', default="merged.fits")

    args = parser.parse_args()

    appendVarFits(fitsList=args.fitsList, outfits=args.outfits, head=args.head)

