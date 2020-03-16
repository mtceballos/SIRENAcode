#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 13:58:49 2019

@author: ceballos
"""
from __future__ import print_function
from astropy.io import fits
import numpy as np


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
            description='Recalculate GRADE2 according to arrival time in pixel',
            prog='recalcGrade2')

    parser.add_argument('--inFile', required=True,
                        help='File with photons in same pixel where\
                        GRADE2 columns is to be corrected')
    parser.add_argument('--maxGrade2', type=int, default=9999,
                        help='Max grade2 for photons without previous')
    parser.add_argument('--samprate', type=float, default=156250.,
                        help='Sampling rate')
    inargs = parser.parse_args()
    inFile = inargs.inFile
    max2 = inargs.maxGrade2
    samprate = inargs. samprate

    f = fits.open(inFile, mode='update')
    nrows = f['EVENTS'].header["NAXIS2"]
    ftab = f[1].data
    times = np.array(ftab['TIME'])
    grades2 = np.array(ftab['GRADE2'])

    newgrades2 = np.array(np.zeros(nrows), dtype=np.int32)
    newgrades2[0] = max2
    for i in range(1, nrows):
        newgrades2[i] = np.int32((times[i] - times[i-1])*samprate)
        # if not newgrades2[i] == grades2[i]:
        #    print("GRADE2 for photon", i+1, "is", grades2[i],
        #          " and it should be", newgrades2[i])

    ftab['GRADE2'] = newgrades2
    f.close()
