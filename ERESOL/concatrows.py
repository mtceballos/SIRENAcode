#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 10:47:51 2020

@author: ceballos
"""

from astropy.io import fits
import numpy as np

def concatrows(inFile, ext, outFile, ncon):
    """
    Concatenate nrows contiguous rows in a fitsfile to get larger records

    :type inFile: str
    :param inFile: 1-col input FITS file name (with short records)

    :type outFile: str
    :param outFile: output FITS file name (with larger records)

    :type ext: int
    :param ext: extension number

    :type ncon: int
    :param ncon: number of rows to be concatenated

    """

    f = fits.open(inFile)
    data = f[ext].data
    colname = f[ext].header["TTYPE1"]
    coldim0 = f[ext].header["TFORM1"]
    coldim0 = coldim0.rstrip('E')
    nrows0 = f[ext].header["NAXIS2"]
    coldata = data[colname]
    f.close()

    nrows1 = nrows0//ncon
    coldim1 = int(coldim0) * ncon
    data2 = np.zeros((nrows1, coldim1))
    coldim1 = str(coldim1) + 'E'

    j = 0
    for i in range(0, nrows0, ncon):
        if j > nrows1-1:
            break
        # data2[j] = np.concatenate((coldata[i], coldata[i+1], coldata[i+2]))
        pp = coldata[i]
        for k in range(1, ncon):
            pp = np.concatenate((pp, coldata[i+k]))
        data2[j] = pp
        j = j + 1

    newcol1 = fits.Column(name='ADC', format=coldim1, array=data2)
    t = fits.BinTableHDU.from_columns([newcol1, ])
    t.writeto(outFile)
