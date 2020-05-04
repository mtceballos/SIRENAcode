#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 10:59:03 2020

@author: ceballos
"""

from astropy.io import fits


def addkeys(fitsfile, ext, keynames, keyvals):
    """
    Add key to FITS file name in given extension
    :type fitsfile: str
    :param fitsfile: FITS file name

    :type ext: int
    :param ext: extension number

    :type keynames: list
    :param keynames: keyword name

    :type keyvals: list
    :param keyvals: keyword value

    """

    f = fits.open(fitsfile, mode='update')
    for i in range(len(keynames)):
        print("Open file ", fitsfile, "to add", keynames[i])
        f[ext].header[keynames[i]] = keyvals[i]
    f.close()
    print("Closed")
