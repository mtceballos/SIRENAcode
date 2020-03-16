import numpy as np
from astropy.io import fits
import shlex
from subprocess import check_call, STDOUT
import os
from shutil import copy


def toXIFUSIMfits(input_nonstd_file, model_fits_file, output_fits_file):

    """Reads in an non-std file and converts it into a FITS file
    (similar to model_fits_file, i.e., with TESRECORDS)

    Arguments:

        - input_nonstd_file: original non-std file to read
        - model_fits_file: Model FITS file with TESRECORDS extension
        - output_fits_file: std-XIFUSIM-like output FITS file

    """

    copy(model_fits_file, output_fits_file)

    f = fits.open(input_nonstd_file)
    data = f[1].data
    colname = f[1].header["TTYPE1"]
    coldim0 = f[1].header["TFORM1"]
    coldim0 = coldim0.rstrip('E')
    # print("coldim0=", coldim0)
    record_length = int(coldim0)
    nb_records = f[1].header["NAXIS2"]
    coldata = data[colname]
    f.close()

    # create temp file with ADC column
    coldimADC = coldim0 + "E"
    newcolADC = fits.Column(name='ADC', format=coldimADC, array=coldata)
    ADC = fits.BinTableHDU.from_columns([newcolADC, ])
    ADC.writeto('fileADC.fits')
    print("ADC column created")

    # create temp file with TIME column
    coldimTIME = 1
    coldimTIME = str(coldimTIME) + 'D'
    time = []
    samprate = 156250
    for i in range(0, nb_records):
        time.insert(i, i*record_length/samprate)
    newcolTIME = fits.Column(name='TIME', format=coldimTIME, array=time)
    TIME = fits.BinTableHDU.from_columns([newcolTIME, ])
    TIME.writeto('fileTIME.fits')
    print("TIME column created")

    # create temp file with PIXID column
    coldimPIXID = 1
    coldimPIXID = str(coldimPIXID) + 'J'
    pixid = np.ones(nb_records)
    newcolPIXID = fits.Column(name='PIXID', format=coldimPIXID, array=pixid)
    PIXID = fits.BinTableHDU.from_columns([newcolPIXID, ])
    PIXID.writeto('filePIXID.fits')
    print("PIXID column created")

    # create temp file with PH_ID column
    coldimPH_ID = 1
    coldimPH_ID = str(coldimPH_ID) + 'J'
    ph_id = np.zeros(nb_records)
    newcolPH_ID = fits.Column(name='PH_ID', format=coldimPH_ID, array=ph_id)
    PH_ID = fits.BinTableHDU.from_columns([newcolPH_ID, ])
    PH_ID.writeto('filePH_ID.fits')
    print("PH_ID column created")

    # Add new columns to first TIME file
    comm = "faddcol infile=fileTIME.fits colfile=fileADC.fits colname=ADC"
    # print(comm)
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm,
              " to add ADC column from fileADC.fits to fileTIME.fits")
        raise
    os.remove("fileADC.fits")
    print("TIME+ADC => File fileADC.fits removed!")

    comm = "faddcol fileTIME.fits filePIXID.fits PIXID"
    # print(comm)
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm,
              " to add PIXID column from filePIXID.fits to fileTIME.fits")
        raise
    os.remove("filePIXID.fits")
    print("TIME+ADC+PIXID => File filePIXID.fits removed!")

    comm = "faddcol fileTIME.fits filePH_ID.fits PH_ID"
    # print(comm)
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm,
              " to add PH_ID column from filePH_ID.fits to fileTIME.fits")
        raise
    os.remove("filePH_ID.fits")
    print("TIME+ADC+PIXID+PH_ID => File filePH_ID.fits removed!")

    # append new fits file (new columns) to copy of template STD-XIFUSIM file
    comm = "fappend infile=fileTIME.fits[1] outfile=" + str(output_fits_file)
    # print(comm)
    print("Copying new extension in output FITS file...")
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm,
              " to append extension from fileTIME.fits to output_fits_file")
        raise
    os.remove("fileTIME.fits")
    print("File fileTIME.fits removed!")

    print("Copying TESRECORDS header in the new extension...")
    comm = ("cphead infile=" + model_fits_file + "[1] outfile=" +
            output_fits_file + "[9]")
    # print(comm)
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " to cphead TESRECORDS")
        raise

    comm = "fdelhdu infile=" + output_fits_file + "[1] confirm=N proceed=Y"
    # print(comm)
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " cphead")
        raise
    print("Old TESRECORDS extension in the FITS file removed!")

    # modify RECLEN in TRIGGERPARAM extension
    try:
        comm = ("fparkey value=" + str(record_length) + " fitsfile=" +
                output_fits_file + "[TRIGGERPARAM] keyword=RECLEN")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " to change RECLEN in TRIGGERPARAM")
        raise

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
            description='Reads in a non-std-XIFUSIM FITS file and converts'
            ' it into a std-XIFUSIM FITS file', prog='toXIFUSIMfits')
    parser.add_argument('--inFile', help='non-std-XIFUSIM FITS file',
                        required=True)
    parser.add_argument('--modelFile', help='Model FITS file',
                        required=True)
    parser.add_argument('--outFile', help='output std-XIFUSIM file',
                        required=True)

    inargs = parser.parse_args()

    toXIFUSIMfits(input_nonstd_file=inargs.inFile,
                  model_fits_file=inargs.modelFile,
                  output_fits_file=inargs.outFile)
