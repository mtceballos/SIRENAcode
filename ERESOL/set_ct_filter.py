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

        
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
            description='Set central part of optimal filter as a ct',
            prog='set_ct_filter')
    parser.add_argument('--infile', required=True, help='File with filter to be replaced')
    parser.add_argument('--ext', required=True, type=int, help='FITS Extension of filter ')
    parser.add_argument('--colname', required=True, help='Filter colname')
    parser.add_argument('--ct', default="0", help="Constant value to replace central part of filter")
    parser.add_argument('--init', type=int, default="0", help="initial sample to replace central part of filter")
    parser.add_argument('--fin', type=int, default="0", help="final sample to replace central part of filter")
    inargs = parser.parse_args()
    infile = inargs.infile
    ext = inargs.ext
    colname = inargs.colname
    ct = inargs.ct
    init = inargs.init
    fin = inargs.fin
    
    # add new column to outfile
    f = fits.open(infile, mode='update')
    data = f[ext].data
    coldata = data[colname]
    coldata[0,(init-1):fin] = ct
    f.flush()
    f.close()

