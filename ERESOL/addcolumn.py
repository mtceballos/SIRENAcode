#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 10:55:07 2020

@author: ceballos
"""

from astropy.io import fits

def addcolumn(inFile, ext, colname, datacol):
    """
    Add column to existing file
    :type inFile: str
    :param inFile: input FITS file name

    :type colname: str
    :param colname: column name

    :type outFile: str
    :param outFile: output FITS file name (with larger records)

    :type ext: int
    :param ext: extension number

    :type datacol: np array
    :param datacol: data for the column (same length as existing columns)


    """

    f = fits.open(inFile)
    nrows0 = f[ext].header["NAXIS2"]
    f.close()

    coldim1 = str(nrows0) + 'E'
    newcol1 = fits.Column(name=colname, format=coldim1, array=datacol)
    t = fits.BinTableHDU.from_columns([newcol1, ])
    t.writeto(inFile)
