"""
RESOLUTION CURVES for single/pairs of pulses

python getEresolCurves_manySeps.py

!!!!!!! Things to review evry time it is run (they can change):
        =======================================================
        ** Check resultsDir (nodetSP in case secondaries are not detected)
                            (detSP in case secondaries are also detected
                        gainScale in case results are for gain scale curves)
        ** Check library (LONG, SHORT) & filenames
        ** Check jitter

"""

# ----IMPORT MODULES --------------
from __future__ import print_function
import shutil
import numpy as np
import auxpy
from astropy.io import fits


def convertEnergies(inFile, outFile, coeffsFile, alias):
    """
    :param inFile: absolute path to file with reconstructed
                    (non-calibrated) energies
    :param coeffsFile: file with coefficients of polynomial
                    fit to gain curves from polyfit2bias.R
    :param alias: string to select reconstruction type in
                    the coefficients table
    :param outFile: file with reconstructed/calibrated energies
                    for the input pulses
    """

    # ------------------------------------
    # --- Process input data file  -------
    # ------------------------------------

    f = fits.open(inFile, memmap=True)
    nrows = f[1].header["NAXIS2"]
    assert nrows > 0, "Empty evt file (%s): nrows=0 " % inFile

    # read Erecons (SIGNAL) column in numpy array (in keV)
    ftab = f[1].data
    EreconKeV = np.array(ftab['SIGNAL'])

    # calculate corrected energies with polyfit coeffs
    print("...Calculating corrected energies for pulses in ", inFile)
    EcorrKeV = auxpy.enerToCalEner(EreconKeV, coeffsFile, alias)

    # close input file
    f.close()
    del f[1].data

    # ------------------------------------
    # --- Create and populate output file
    # ------------------------------------
    shutil.copy(inFile, outFile)
    f = fits.open(outFile, memmap=True, mode='update')
    hdr = f[0].header
    hdr['HISTORY'] = ("Updated by convertEnergies.py to correct "
                      "energies by gainscale")
    ftab = f[1].data
    ftab['SIGNAL'] = EcorrKeV
    f.close()
    del f[1].data


#
# --------- Read input parameters and PROCESS simulations -----------------
#


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Convert input energies to \
                                     calibrated according to gainScale',
                                     prog='convertEnergies.py')

    parser.add_argument('--inFile', required=True,
                        help='Input FITS file (EVENTS) with \
                        reconstructed energies')
    parser.add_argument('--outFile', required=True,
                        help='Output FITS file (EVENTS) with \
                        calibrated energies')
    parser.add_argument('--coeffsFile', help='file with polynomial fit \
                        coefficients from R script polyfit2bias.R', default='')
    parser.add_argument('--alias', help='String for reconstruction method \
                        in coefficient file (i.e. fixedlib6OF_OPTFILT4096)',
                        required=True)

    inargs = parser.parse_args()

    # print("array=",inargs.array)

    convertEnergies(inFile=inargs.inFile, outFile=inargs.outFile,
                    coeffsFile=inargs.coeffsFile, alias=inargs.alias)
