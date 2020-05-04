#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 13:19:32 2020

@author: ceballos
"""

import numpy as np
import json
from scipy.interpolate import interp1d


def enerToCalEner(inEner, inPhase, coeffsFile, alias):
    """
    :param inEner: numpy array with input uncorrected energies
    :param inPhase: numpy array with input phases (jitter)
    :param coeffsFile: file with coefficients of polynomial fit to gain curves
                        from polyfit2bias.R or surface 2D polynomial or
                        JSON file with data points for spline fit
    :param alias: string to select reconstruction type in the
                        coefficients table (if gainScale curve)
    :return calEner: numpy vector with calibrated energies
    """

    # locate coefficients in calibration table
    # ----------------------------------------
    coeffsDict = dict()
    with open(coeffsFile, "rt") as f:
        fileCont = f.read()   # JSON file
        if 'surface' in fileCont:
            ftype = "surface"
        elif fileCont[0] == '{':
            ftype = 'json'
        else:
            ftype = 'poly'

    if ftype == 'poly':
        codata = ascii.read(coeffsFile, guess=False, format='basic')
        # codata[1] : row 2
        # codata[1][1]: row 2, col 2
        print("Reading curve coefficients from", coeffsFile, "\n")
        for i in range(0, len(codata)):
            #  METHOD   ALIAS  a0  a1  a2  a3  a4
            coeffsDict[codata[i][1]] = (codata[i][2], codata[i][3],
                                        codata[i][4], codata[i][5],
                                        codata[i][6])
        npCoeffs = np.array(coeffsDict[alias])
    elif ftype == 'surface':
        print("Reading surface coefficients from", coeffsFile, "\n")
        npCoeffs = np.loadtxt(coeffsFile, comments="#")
    elif ftype == 'json':
        with open(coeffsFile, 'r') as f:
            jsonDict = json.load(f)
        nalias = len(jsonDict["ALIAS"])
        aliasidx = jsonDict["ALIAS"].index(alias)
        if len(jsonDict["xdata"]) % nalias == 0:
            nEnerCal = len(jsonDict["xdata"])//nalias
        else:
            raise ValueError("Length of xdata is not a multiple of number"
                             " of calibration energies")
        stridx = nEnerCal*aliasidx
        endidx = stridx + nEnerCal
        xdata = jsonDict["xdata"][stridx:endidx]
        ydata = jsonDict["ydata"][stridx:endidx]
        funinterp = interp1d(xdata, ydata, kind="linear",
                             fill_value="extrapolate")
    else:
        raise ValueError("Incorrect Coeffs file type")

    calEner = np.zeros(inEner.size, dtype=float)
    ie = 0

    # print("npCoeffs=",npCoeffs)
    #
    # convert energies
    #
    if ftype == 'surface':
        # print("Using surface to calibrate reconstructed energies\n")
        calEner = np.polynomial.polynomial.polyval2d(inEner, inPhase, npCoeffs)
    elif ftype == 'poly':
        for ie in range(0, inEner.size):
            # read fitting coeffs taken from polyfit2Bias.R (a0, a1, a2, a3)
            #  as in y = a0 + a1*x + a2*x^2 + a3*x^3
            # where y=E_reconstructed and x=Ecalibration (keV)

            # print("Using curve to calibrate reconstructed energies\n")
            # subtract "y" (0 = a0 + a1*x + a2*x^2 + ... - y) :
            npCoeffs = np.array(coeffsDict[alias])
            npCoeffs[0] -= inEner[ie]
            # reversed to say fit with poly1d definition:
            npCoeffsRev = npCoeffs[::-1]
            polyfit = np.poly1d(npCoeffsRev)
            # get real root (value of Ereal for a given Erecons )
            r = np.roots(polyfit)
            # real && >0 roots
            # print("r=", r)
            rreal = r.real[abs(r.imag) < 1e-5]
            # print("rreal=", rreal)
            rrealpos = rreal[rreal > 0]
            # print("rrealpos=", rrealpos)
            # closest root
            # print(inEner[ie])
            rclosest = min(enumerate(rrealpos),
                           key=lambda x: abs(x[1]-inEner[ie]))[1] # (idx,value)
            # print("For:", alias, " Recon energy=", rclosest)

            calEner[ie] = rclosest

    elif ftype == 'json':
        calEner = inEner/funinterp(inEner)

    # return calibrated energies
    # ------------------------------------
    return calEner
