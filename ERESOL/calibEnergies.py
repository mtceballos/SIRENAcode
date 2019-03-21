"""
Calibrate reconstructed energies file based on coefficients of gain scale

python convertEnergies.py

"""

# ----IMPORT MODULES --------------
from __future__ import print_function
import auxpy

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

    auxpy.convertEnergies(inFile=inargs.inFile, outFile=inargs.outFile,
                          coeffsFile=inargs.coeffsFile, alias=inargs.alias)
