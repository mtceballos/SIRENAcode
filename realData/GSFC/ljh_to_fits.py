#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri May 20 10:47:47 2020

@author: GSFC + Maite Ceballos
"""

import xcalpy.pulse as xcp
import xcalpy.pulse.Units as xcu
import numpy as np
import shlex
import os
import argparse
from astropy.io import fits
from datetime import datetime
from subprocess import check_call, STDOUT


def ljh_to_fits(input_ljh_file, output_fits_file, model_fits_file, verbosity=0, clobber='no'):

    """Reads in an GSFC LJH file and converts it in a FITS file
    (similar to XIFUSIM model_fits_file but with TESRECORDS as in LHJ)

    Arguments:

        - input_ljh_file: LJH file to read in

        - output_fits_file: Resulting FITS file

        - model_fits_file: Model FITS file (from XIFUSIM simulations) whose TESRECORDS extension will be changed

        - verbosity: verbosity level for output messages (0 for silence)

        - clobber: "yes" or "no". If "yes", it overwrites output file.
    """

    if verbosity > 0:
        print("Reading from ljh file...")


    if (clobber.lower()=="yes" or clobber.lower()=="y") and os.path.isfile(output_fits_file):
        print("Output file", output_fits_file, "exists and it will be overwritten")
        os.remove(output_fits_file)

    # read information form LJH file
    units = xcp.UnitPars(xcu.RAW)
    ljh = xcp.LJHFile(input_ljh_file)
    numrecs = ljh.numRecords
    channel = ljh.channel
    rec = ljh.readRecord(0,units=units)
    record_length = ljh.numSamples
    deltat = rec.source.sampleTime

    # check consistency of file name with channel info
    channel_from_name = input_ljh_file[10:-5]
    if not int(channel_from_name) == channel:
        print("Filename does not seems to be consistent with channel info stored in LJH file")
        proceed = input("Do you want to continue?[Y]/N:")
        if proceed == "N" or proceed=="n":
            print("File ", input_ljh_file, "has not been converted to FITS")
            return

    # inititalize storage of ADC and TIME columns
    records = np.zeros((numrecs,record_length))
    times = np.zeros((numrecs))
    ph_id = np.zeros(numrecs)
    pixid = np.zeros(numrecs)
    pixid[:] = channel

    for ir in range(numrecs):
        rec = ljh.readRecord(ir,units=units)
        if not rec.source.numSamples == record_length:
            print("Warning: need to check variable length records")
        # flip record about horizontal line of ~baseline and then subtract baseline
        constant = 40000
        flip = 47000
        newrec = (2.*flip - rec.record) - constant
        records[ir,:] = newrec
        times[ir] = rec.time
        ph_id[ir] = ir


    coldimADC = record_length
    coldimADC = str(coldimADC) + 'E'
    newcolADC = fits.Column(name='ADC', format=coldimADC, array=records)
    ADC = fits.BinTableHDU.from_columns([newcolADC, ])
    ADC.writeto('fileADC.fits')
    if verbosity > 0:
        print("ADC column created")

    coldimTIME = 1
    coldimTIME = str(coldimTIME) + 'D'
    newcolTIME = fits.Column(name='TIME', format=coldimTIME, array=times)
    TIME = fits.BinTableHDU.from_columns([newcolTIME, ])
    TIME.writeto(output_fits_file)
    if verbosity > 0:
        print("TIME column created")

    coldimPIXID = 1
    coldimPIXID = str(coldimPIXID) + 'J'
    newcolPIXID = fits.Column(name='PIXID', format=coldimPIXID, array=pixid)
    PIXID = fits.BinTableHDU.from_columns([newcolPIXID, ])
    PIXID.writeto('filePIXID.fits')
    if verbosity > 0:
        print("PIXID column created")

    coldimPH_ID = 1
    coldimPH_ID = str(coldimPH_ID) + 'J'
    newcolPH_ID = fits.Column(name='PH_ID', format=coldimPH_ID, array=ph_id)
    PH_ID = fits.BinTableHDU.from_columns([newcolPH_ID, ])
    PH_ID.writeto('filePH_ID.fits')
    if verbosity > 0:
        print("PH_ID column created")

    comm = "faddcol " + output_fits_file + " fileADC.fits ADC"
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except:
        print("Error running ", comm, " to add ADC column from fileADC.fits to ", output_fits_file)
        os.remove("fileADC.fits")
        os.remove("filePH_ID.fits")
        os.remove("filePIXID.fits")
        raise
    os.remove("fileADC.fits")
    if verbosity > 0:
        print("ADC => File fileADC.fits removed!")

    comm = "faddcol " + output_fits_file + " filePIXID.fits PIXID"
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except:
        print("Error running ", comm, " to add PIXID column from filePIXID.fits to ", output_fits_file)
        os.remove("filePH_ID.fits")
        os.remove("filePIXID.fits")
        raise
    os.remove("filePIXID.fits")
    if verbosity > 0:
        print("PIXID => File filePIXID.fits removed!")

    comm = "faddcol " + output_fits_file + " filePH_ID.fits PH_ID"
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except:
        print("Error running ", comm, " to add PH_ID column from filePH_ID.fits to ", output_fits_file)
        os.remove("filePH_ID.fits")
        raise
    os.remove("filePH_ID.fits")
    if verbosity > 0:
        print("PH_ID => File filePH_ID.fits removed!")
        print("Copying TESRECORDS header in the new extension...")

    comm = "cphead infile=" + model_fits_file + "[1] outfile=" + output_fits_file + "[1]"
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError as err:
        print("Error running ", comm, " cphead")
        print(err)
        raise

    # Update header with real values
    f = fits.open(output_fits_file, mode='update')
    f["TESRECORDS"].header["DELTAT"] = (deltat, "Sample time (s); inverse of sampling rate")
    f["TESRECORDS"].header["RECLEN"] = (record_length, "Record length (samples)")
    now = datetime.now()
    date_time = now.strftime("%Y-%m-%dT%H:%M:%S")
    f["TESRECORDS"].header["DATE-HDU"] = date_time.rstrip()
    comment = "GSFC data filename"
    f[0].header["LJHFILE"] = (os.path.basename(input_ljh_file), comment)
    f.close()


desc="Convert LJH file to a xifusim-like FITS file to be processed with SIRENA."
parser=argparse.ArgumentParser(description=desc)
parser.add_argument("--infile",metavar="FILE",type=str,help="name of LJH file", required=True)
parser.add_argument("--outfile",metavar="FILE",type=str,help="name of FITS file", required=True)
parser.add_argument("--modelfile",metavar="FILE",type=str,help="name of xifusim-like FITS template file", required=True)
parser.add_argument("--clobber",type=str,choices=['yes','no'],help="Overwrite output if clobber=yes", default="no")
parser.add_argument("--verbosity",type=int,help="verbosity level", default=0)
inargs=parser.parse_args()

ljh_to_fits(input_ljh_file=inargs.infile,
            output_fits_file=inargs.outfile,
            model_fits_file=inargs.modelfile,
            clobber=inargs.clobber,
            verbosity=inargs.verbosity)
