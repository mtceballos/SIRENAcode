"""
Substitute optimal filter by a new array of values

python set_fit_filter.py

"""

# ----IMPORT MODULES --------------
from astropy.io import fits, ascii
import numpy as np
        
#if __name__ == "__main__":
#    import argparse

#    parser = argparse.ArgumentParser(
#            description='Set central part of optimal filter as a ct',
#            prog='set_ct_filter')
#    parser.add_argument('--infile', required=True, help='File with filter to be replaced')
#    parser.add_argument('--ext', required=True, type=int, help='FITS Extension of filter ')
#    parser.add_argument('--colname', required=True, help='Filter colname')
#    parser.add_argument('--arrayfile', help="File with values for the replacement")
#    parser.add_argument('--init', type=int, default="0", help="initial sample to replace central part of filter")
#    parser.add_argument('--fin', type=int, default="0", help="final sample to replace central part of filter")
#    inargs = parser.parse_args()
#    infile = inargs.infile
#    ext = inargs.ext
#    colname = inargs.colname
#    arrayfile = inargs.arrayfile
#    init = inargs.init
#    fin = inargs.fin
    

def replaceCol(infile, ext, colname, arrayfile, init, fin):
    f = fits.open(infile, mode='update')
    data = f[int(ext)].data
    coldata = data[colname]
    newdata = np.loadtxt(arrayfile)
    coldata[0,(int(init)-1):int(fin)] = newdata
    f.flush()
    f.close()

