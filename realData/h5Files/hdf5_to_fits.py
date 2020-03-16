import h5py
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import shlex
from subprocess import check_call, STDOUT
import os
from shutil import copy

def hdf5_to_fits(input_hdf5_file, output_fits_file, model_fits_file, use_phase = True, channel_number = 0, freq_number = 6):

    """Reads in an HDF5 file and converts it in a FITS file (similar to model_fits_file but with TESRECORDS as in HDF5)

    Arguments:

        - input_hdf5_file: HDF5 file to read in
        
        -output_fits_file: Resulting FITS file
        
        - model_fits_file: Model FITS file whose TESRECORDS extension is going to be changed

        - use_phase: option to return A*cos(phase) instead of just the amplitude column

        - channel_number: channel number identifier

        - freq_number: frequency/pixel number identifier

    """
    
    #output_fits_file = input_hdf5_file[0:len(input_hdf5_file)-3]  # Delete the '.h5'
    output_fits_fileAUX = input_hdf5_file[0:len(input_hdf5_file)-3] + 'AUX.fits'
    copy(model_fits_file,output_fits_fileAUX)
    #output_fits_file = output_fits_file + '.fits'
    
    print("Reading from HDF5 file...")

    f = h5py.File(input_hdf5_file, 'r')

    #data = f["mux"]["channel_%03d" % channel_number]["freq_%03d" % freq_number]
    #data = f["mux"]["channel_000"]["freq_006"]
    data = f["mux"]["channel_00" + str(channel_number)]["freq_00" + str(freq_number)]

    #nb_records = len(data.keys())
    nb_records = len(data.keys())-1 # For file2.h5 because it has a table more called 'bea'
    print("nb_records: ",nb_records)

    record_length = len(data["000000000"].value[:,0])
    #print("record_length: ",record_length)

    records = np.zeros((nb_records,record_length))
    
    #print("len(data.keys()): ",len(data.keys()))

    #for record_number in range(len(data.keys())):
    for record_number in range(len(data.keys())-1): # For file2.h5 because it has a table more called 'bea'

        #record = data[str(record_number)]
        record = data["%09d" % record_number]
        #print("record: ", record)

        amplitude = record.value[:,0]
        #print("amplitude: ", amplitude)

        if use_phase:

            phase = record.value[:,0]*record.attrs["scale_phase_rad"]

            amplitude = amplitude*np.cos(phase)

        #records[record_number] = -amplitude
        records[record_number] = -amplitude+2*amplitude[0] # To change the pulse polarity
        
    #print("records: ",records)
    #print(np.array(records[0]))
    #print(np.array(records[0,0]))
       
    """index = 8001
    samples = list(range(record_length))
    amplitudeRecordi = []
    for i in range(0,len(samples)):
        amplitudeRecordi.insert(i,np.array(records[index][i]))
    plt.plot(samples,amplitudeRecordi)
    plt.xlabel('(samples)')
    plt.ylabel('Record_amplitude (adu?) (i=' + str(index) + ')')
    plt.show()"""
    
    coldimADC = record_length 
    coldimADC = str(coldimADC) + 'E'
    newcolADC = fits.Column(name='ADC', format=coldimADC, array=records)
    ADC = fits.BinTableHDU.from_columns([newcolADC, ])
    ADC.writeto('fileADC.fits')
    print("ADC column created")
    
    coldimTIME = 1 
    coldimTIME = str(coldimTIME) + 'D'
    time = []
    samprate = 156250
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
    #print(comm)
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " to add ADC column from fileADC.fits to fileTIME.fits")
        raise
    os.remove("fileADC.fits")
    print("TIME+ADC => File fileADC.fits removed!")
        
    comm = "faddcol fileTIME.fits filePIXID.fits PIXID"
    #print(comm)
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " to add PIXID column from filePIXID.fits to fileTIME.fits")
        raise
    os.remove("filePIXID.fits")
    print("TIME+ADC+PIXID => File filePIXID.fits removed!")
        
    comm = "faddcol fileTIME.fits filePH_ID.fits PH_ID"
    #print(comm)
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " to add PH_ID column from filePH_ID.fits to fileTIME.fits")
        raise
    os.remove("filePH_ID.fits")
    print("TIME+ADC+PIXID+PH_ID => File filePH_ID.fits removed!")
    
    comm = "fappend infile=fileTIME.fits[1] outfile=" + str(output_fits_fileAUX)
    #print(comm)
    print("Copying new extension in output FITS file...")
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " to append extension from fileTIME.fits to output_fits_file")
        raise
    os.remove("fileTIME.fits")
    print("File fileTIME.fits removed!")
    
    #copy("prueba.fits","pruebaBIG.fits")
    print("Copying TESRECORDS header in the new extension...")
    copy(output_fits_fileAUX,output_fits_file)
    
    #comm = "cphead infile=prueba.fits[1] outfile=pruebaBIG.fits[9]"
    comm = "cphead infile=" + output_fits_fileAUX + "[1] outfile=" + output_fits_file + "[9]"
    #print(comm)
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " cphead")
        raise
    #os.remove("prueba.fits")
    os.remove(output_fits_fileAUX)
    print("File output_fits_fileAUX removed!")

    #comm = "fdelhdu infile=pruebaBIG.fits[1] confirm=N proceed=Y"
    comm = "fdelhdu infile=" + output_fits_file + "[1] confirm=N proceed=Y"
    #print(comm)
    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running ", comm, " cphead")
        raise
    print("Old TESRECORDS extension in the FITS file removed!")

    return np.array(records)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Reads in an HDF5 file and converts it in a FITS file (similar to model_fits_file but with TESRECORDS as in HDF5)',
                                     prog='hdf5_to_fits')
    parser.add_argument('--input_hdf5_file', help='HDF5 file to read in', required=True)
    parser.add_argument('--output_fits_file', help='Resulting FITS file', required=True)
    parser.add_argument('--model_fits_file', help='Model FITS file whose TESRECORDS extension is going to be changed', required=True)
    parser.add_argument('--use_phase', help='option to return A*cos(phase) instead of just the amplitude column', required=False)
    parser.add_argument('--channel_number', help='channel number identifier', required=True)
    parser.add_argument('--freq_number', help='frequency/pixel number identifier', required=True)

inargs = parser.parse_args()

hdf5_to_fits(input_hdf5_file=inargs.input_hdf5_file,output_fits_file=inargs.output_fits_file, model_fits_file=inargs.model_fits_file, use_phase=inargs.use_phase,channel_number=inargs.channel_number,freq_number=inargs.freq_number)

