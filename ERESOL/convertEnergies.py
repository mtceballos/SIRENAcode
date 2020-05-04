#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 12:36:32 2020

@author: ceballos
"""

from os import copy
import numpy as np
from astropy.io import fits
from enerToCalEner import enerToCalEner

def convertEnergies(inFile, outFile, coeffsFile, alias):
    """
    :param inFile: absolute path to file with reconstructed
                    (non-calibrated) energies
    :param coeffsFile: file with coefficients of polynomial
                    fit to gain curves from polyfit2bias.R
                    or JSON file with data points for spline fitting
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
    # reconPhase = np.array(ftab['PHI'])
    reconPhase = np.array(ftab['PHI']) + np.array(ftab['LAGS'])

    # calculate corrected energies with polyfit coeffs
    # 1) check if pulses are of different length (different gain scale)
    print("...Calculating corrected energies for pulses in ", inFile)
    EcorrKeV = enerToCalEner(EreconKeV, reconPhase, coeffsFile, alias)

    # close input file
    f.close()
    del f[1].data

    # ------------------------------------
    # --- Create and populate output file
    # ------------------------------------
    copy(inFile, outFile)
    f = fits.open(outFile, memmap=True, mode='update')
    hdr = f[0].header
    hdr['HISTORY'] = ("Updated by convertEnergies.py to correct "
                      "energies by gainscale using", coeffsFile)
    ftab = f[1].data
    ftab['SIGNAL'] = EcorrKeV
    f.close()
    del f[1].data
