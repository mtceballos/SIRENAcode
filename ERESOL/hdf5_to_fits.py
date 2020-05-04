#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 10:47:47 2020

@author: Philippe Peille + Bea Cobo
"""


import h5py
import numpy as np
from astropy.io import fits
import shlex
from subprocess import check_call, STDOUT
import os
from shutil import copy

def hdf5_to_fits(input_hdf5_file, output_fits_file, model_fits_file, use_phase = True,
                 channel_number = 0, freq_number = 6, sampling_rate=156250):

    """Reads in an HDF5 file and converts it in a FITS file
    (similar to XIFUSIM model_fits_file but with TESRECORDS as in HDF5)

    Arguments:

        - input_hdf5_file: HDF5 file to read in

        - output_fits_file: Resulting FITS file

        - model_fits_file: Model FITS file (from XIFUSIM simulations) whose TESRECORDS extension will be changed

        - use_phase: option to return A*cos(phase) instead of just the amplitude column

        - channel_number: channel number identifier

        - freq_number: frequency/pixel number identifier

        - sampling_rate: sampling frequency of data acquisition

    """

    output_fits_fileAUX = input_hdf5_file[0:len(input_hdf5_file)-3] + 'AUX.fits'
    copy(model_fits_file,output_fits_fileAUX)

    print("Reading from HDF5 file...")

    f = h5py.File(input_hdf5_file, 'r')

    data = f["mux"]["channel_00" + str(channel_number)]["freq_00" + str(freq_number)]

    nb_records = len(data.keys())-1 # For file2.h5 because it has a table more called 'bea'
    print("nb_records: ",nb_records)

    record_length = len(data["000000000"].value[:,0])
    records = np.zeros((nb_records,record_length))


    for record_number in range(len(data.keys())):

        record = data["%09d" % record_number]
        amplitude = record.value[:,0]

        if use_phase:

            phase = record.value[:,0]*record.attrs["scale_phase_rad"]

            amplitude = amplitude*np.cos(phase)

        #records[record_number] = -amplitude
        records[record_number] = -amplitude+2*amplitude[0] # To change the pulse polarity

    coldimADC = record_length
    coldimADC = str(coldimADC) + 'E'
    newcolADC = fits.Column(name='ADC', format=coldimADC, array=records)
    ADC = fits.BinTableHDU.from_columns([newcolADC, ])
    ADC.writeto('fileADC.fits')
    print("ADC column created")

    coldimTIME = 1
    coldimTIME = str(coldimTIME) + 'D'
    time = []
    samprate = sampling_rate
    for i in range(0,nb_records):
        time.insert(i,i*record_length/samprate)
    newcolTIME = fits.Column(name='TIME', format=coldimTIME, array=time)
    TIME = fits.BinTableHDU.from_columns([newcolTIME, ])
    TIME.writeto('fileTIME.fits')
    print("TIME column created")

    coldimPIXID = 1
    coldimPIXID = str(coldimPIXID) + 'J'
    pixid = np.ones(nb_records)
    newcolPIXID = fits.Column(name='PIXID', format=coldimPIXID, array=pixid)
    PIXID = fits.BinTableHDU.from_columns([newcolPIXID, ])
    PIXID.writeto('filePIXID.fits')
    print("PIXID column created")

    coldimPH_ID = 1
    coldimPH_ID = str(coldimPH_ID) + 'J'
    ph_id = np.zeros(nb_records)
    newcolPH_ID = fits.Column(name='PH_ID', format=coldimPH_ID, array=ph_id)
    PH_ID = fits.BinTableHDU.from_columns([newcolPH_ID, ])
    PH_ID.writeto('filePH_ID.fits')
    print("PH_ID column created")

    comm = "faddcol fileTIME.fits fileADC.fits ADC"
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " to add ADC column from fileADC.fits to fileTIME.fits")
        raise
    os.remove("fileADC.fits")
    print("TIME+ADC => File fileADC.fits removed!")

    comm = "faddcol fileTIME.fits filePIXID.fits PIXID"
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " to add PIXID column from filePIXID.fits to fileTIME.fits")
        raise
    os.remove("filePIXID.fits")
    print("TIME+ADC+PIXID => File filePIXID.fits removed!")

    comm = "faddcol fileTIME.fits filePH_ID.fits PH_ID"
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " to add PH_ID column from filePH_ID.fits to fileTIME.fits")
        raise
    os.remove("filePH_ID.fits")
    print("TIME+ADC+PIXID+PH_ID => File filePH_ID.fits removed!")

    comm = "fappend infile=fileTIME.fits[1] outfile=" + str(output_fits_fileAUX)

    print("Copying new extension in output FITS file...")
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " to append extension from fileTIME.fits to output_fits_file")
        raise
    os.remove("fileTIME.fits")
    print("File fileTIME.fits removed!")

    print("Copying TESRECORDS header in the new extension...")
    copy(output_fits_fileAUX,output_fits_file)

    comm = "cphead infile=" + output_fits_fileAUX + "[1] outfile=" + output_fits_file + "[9]"
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " cphead")
        raise

    os.remove(output_fits_fileAUX)
    print("File output_fits_fileAUX removed!")

    comm = "fdelhdu infile=" + output_fits_file + "[1] confirm=N proceed=Y"
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " cphead")
        raise
    print("Old TESRECORDS extension in the FITS file removed!")

    return np.array(records)
