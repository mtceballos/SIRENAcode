"""
Substitute central part of optimal filter by a constant

python set_ct_filter.py

"""

# ----IMPORT MODULES --------------
from __future__ import print_function
import os
import shlex
from shutil import copy, rmtree
import tempfile
from time import gmtime, strftime
import numpy as np
import math
from subprocess import check_call, STDOUT
from astropy.io import fits, ascii
import json
from scipy.interpolate import interp1d

def VLtoFLnew(inputFile, outputFile):
    """
    Transform variable length column to fixed length using stilts
    Astropy & ftools & stilts does not work with variable-length arrays,
    so a previous conversion to fixed-format is required
    """

    print("\n Converting variable length column to fixed length column")
    extnum = 2
    print("Using:", inputFile, outputFile)
    try:
        # save inputFile into outFile
        copy(inputFile, outputFile)
            
        comm = ("/dataj4/software/32/starjava/bin/stilts tpipe cmd='addcol " +
                "ENERGY2 \"(ENERGY*2)\"' cmd='delcols ENERGY2' ofmt=fits in=" +
                outputFile + "#" + str(extnum) + " out=" + outputFile)
        print(comm)
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        
        comm = ("fappend " + inputFile + "+1 " + outputFile)
        print(comm)
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)

        comm = ("fappend " + inputFile + "+3 " + outputFile)
        print(comm)
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
            
        # copy header (modified by stilts) from one FITS header to another
        # existing FITS file header
        comm = ("cphead " + inputFile + "+" + str(extnum) + " " + outputFile +
                    "+1")
        print(comm)
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
            
    except RuntimeError:
        print("Error Coverting columns: ")
        print(comm)
        raise
        
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
            description='Set central part of optimal filter as a ct',
            prog='set_ct_filter')
    parser.add_argument('--infile', required=True, help='Inpult filter file+ext')
    parser.add_argument('--outfile', required=True, help='Inpult filter file')
    parser.add_argument('--ct', default="0", help="Constant value to replace central part of filter")
    parser.add_argument('--init', type=int, default="0", help="initial sample to replace central part of filter")
    parser.add_argument('--fin', type=int, default="0", help="final sample to replace central part of filter")
    inargs = parser.parse_args()
    infile = inargs.infile
    outfile = inargs.outfile
    ct = inargs.ct
    init = inargs.init
    fin = inargs.fin
    
    # convert to fixed length columns
    VLtoFLnew(inputFile=infile,outputFile="pp.fits")
        
    # add new column to outfile
    f = fits.open("pp.fits", mode='update')
    data = f[1].data
    coldata = data['T8192']
    coldata[0,init:fin] = ct
    f.flush()
    f.close()

    #reorder extensions
    try:
        comm = ("fextract pp.fits+2 " + outfile)
        print(comm)
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)

        comm = ("fappend pp.fits+1 " + outfile)
        print(comm)
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)

        comm = ("fappend pp.fits+3 " + outfile)
        print(comm)
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error sorting extensions: ")
        print(comm)
        raise
