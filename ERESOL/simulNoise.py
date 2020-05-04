#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 11:19:26 2020

@author: ceballos
"""

from os import getcwd, path, chdir
from math import ceil
import shlex
from astropy.io import fits
from subprocess import check_call, STDOUT
from shutil import rmtree
from sixtevars import XIFUSIMinst, XMLsixte, sampfreqs, sampids, sampStrs
import tempfile

tmpDir = tempfile.mkdtemp()

def simulNoise(pixName, samprate, jitter, bbfb, pulseLength,
               space, acbias, scaleFactor, samplesUp, nSgms,
               nintervals, dcmt):
    """ simulate data in input parameter space, calculate the data baseline
            and create Noise file to be ingested in SIRENA processing tasks

          :param pixName : name of extension in pixel definition file
                          SPA*, LPA1*, LPA2*, LPA3*
          :param samprate: Samprate value with respect to baseline
              of 156250 Hz: "" (baseline), "samprate2" (half_baseline)
          :param jitter: jitter option ("" for no_jitter and
                                        "jitter" for jitter)
          :param bbfb: ("") for dobbfb=n or ("bbfb") for dobbfb=y
          :param pulseLength : pulse length in samples to select sampling
                  in noise spectrum
          :param space :  Data space: ADC for current, R or RALL or RNOL
                  or RFITTED for resistance
          :param acbias : AC (acbias=yes) or DC (acbias=no)
          :param scaleFactor : Param scaleFactor for gennoise tool
          :param samplesUp : Param samplesUp for gennoise tool
          :param nSgms : Param nSgms for gennoise tool
          :param nintervals : Number of intervals for gennoisespec
                  noise calculation
          :param dcmt: decimation factor for xifusim jitter simulation
          :return fits file with Noise spectra in specified space and sampling

    """
    # ----GLOBAL VARIABLES -------------
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA"
    simSIXTEdir = EURECAdir + "/testHarness/simulations/SIXTE"
    cwd = getcwd()

    preBuffer = 0
    diffth = 0

    # ---- param dependent variables
    triggerSize = pulseLength + 10

    # samprate
    smprtStr = ""
    idxsmp = sampids.index(samprate)
    smprtStr = sampStrs[idxsmp]
    samplingrate = sampfreqs[idxsmp]

    # set simulation time based on number of intervals required
    simTimeN = ceil(nintervals * triggerSize / samplingrate)

    # jitter
    jitterStr = ""
    offset = ""
    if jitter == "jitter":
        jitterStr = "_jitter"
        jitterStrPix = "_jitter"
        offset = " offset=-1"
    if dcmt > 1:
        jitterStr = "_jitter_dcmt" + str(dcmt)

    # bbfb
    bbfbStr = ""
    XMLxifusim = XIFUSIMinst + "1pix_noise_nobbfb.xml"
    assert not(jitter == "" and bbfb == "bbfb"),\
        "Error: BBFB will produce jitter"

    if bbfb == "bbfb":
        XMLxifusim = XIFUSIMinst + "1pix_noise_bbfb.xml"
        bbfbStr = "_bbfb"

    # energyMethod
    spacePar = "NONE"
    if space == "R":
        spacePar = "I2R=I2R"
    elif space == "RALL":
        spacePar = "I2R=RALL"
    elif space == "RNOL":
        spacePar = "I2R=RNOL"
    elif space == "RFITTED":
        spacePar = "I2R=RFITTED"
    elif space == "ADC":
        spacePar = "I2R=I"

    xifusim = "xifusim" + pixName
    wdir = simSIXTEdir + "/NOISE/" + xifusim
    chdir(wdir)

    # define files
    rootN = ("forNoise" + str(pulseLength) + "samples_" + xifusim +
             "_" + str(simTimeN) + "s")
    noiseFile = ("noise" + str(pulseLength) + "samples_" + xifusim + "_B0_" +
                 space + smprtStr + jitterStr + bbfbStr + ".fits")
    pixFileN = rootN + smprtStr + jitterStrPix + ".piximpact"
    fitsFileN = rootN + smprtStr + jitterStr + bbfbStr + ".fits"

    # -------------------------------------------------------------------------
    # get NOISE file: process empty trigger file
    # -------------------------------------------------------------------------
    if not path.isfile(pixFileN):
        print("\nTESCONSPILEUP: Generating no-events file for NOISE...")
        comm = ("tesconstpileup PixImpList=" + pixFileN +
                " sample_freq=" + str(samplingrate) +
                " XMLFile=" + XMLsixte +
                " tstop=" + str(simTimeN) + offset +
                " energy=0 pulseDistance=1 TriggerSize=" +
                str(triggerSize) + " clobber=yes")
        args = shlex.split(comm)
        try:
            check_call(args, stderr=STDOUT)
        except RuntimeError:
            print("Error running tool tesconstpileup (NOISE):")
            print(comm)
            chdir(cwd)
            rmtree(tmpDir)
            raise

        print("TESCONSPILEUP: ...................................END")

    if not path.isfile(fitsFileN):
        print("\nxifusim: Generating file for NOISE:", fitsFileN)
        tstart = 0.
        tstop = float(simTimeN)

        comm = ("xifusim PixImpList=" + pixFileN +
                " Streamfile=" + fitsFileN +
                " sample_rate=" + str(samplingrate) +
                " tstart=" + str(tstart) + " tstop=" + str(tstop) +
                " trig_reclength=" + str(triggerSize) +
                " trig_n_pre=" + str(preBuffer) +
                " trig_thresh=" + str(diffth) +
                " trig_n_suppress=" + str(pulseLength) +
                " acbias=" + acbias +
                " decimate_factor=" + str(dcmt) +
                " clobber=yes XMLfilename=" + XMLxifusim)

        args = shlex.split(comm)
        try:
            print("Running tool xifusim (NOISE):", comm)
            check_call(args, stderr=STDOUT)
        except RuntimeError:
            print("Error running tool xifusim (NOISE):")
            print(comm)
            chdir(cwd)
            rmtree(tmpDir)
            raise
        print("xifusim: ................................END")
        # sys.exit()

        # it seems (27/03/2019) that there are no problems now
#        try:
#            # rm first record in fits file (noisy)
#            print("rm first record in fits file (noisy)")
#            rmLastAndFirst(fitsFileN, 'TESRECORDS', 0)
#            updateHISTORY(fitsFileN, comm)
#        except RuntimeError:
#            print("Error making up trigger file  (NOISE):")
#            os.chdir(cwd)
#            rmtree(tmpDir)
#            raise

        # and Add TRIGGSZ to xifusim simulated file
        fN = fits.open(fitsFileN, mode='update')
        # fNhdr = fN[1].header
        fN[0].header["TRIGGSZ"] = triggerSize
        fN[0].header["DELTAT"] = dcmt * fN['TESRECORDS'].header["DELTA_T"]
        fN['TESRECORDS'].header["TRIGGSZ"] = triggerSize
        fN['TESRECORDS'].header["DELTAT"] = \
            dcmt * fN['TESRECORDS'].header["DELTA_T"]
        fN.close()

    print("\nGENNOISESPEC: Generating NOISE spectrum file "
          "in (", space, " space)")
    comm = ("gennoisespec --inFile=" + fitsFileN +
            " --" + spacePar +
            " --outFile=" + noiseFile +
            " --intervalMinSamples=" + str(pulseLength) +
            " --nintervals=" + str(nintervals) +
            " --scaleFactor=" + str(scaleFactor) +
            " --samplesUp=" + str(samplesUp) +
            " --nSgms=" + str(nSgms) +
            " --pulse_length=" + str(pulseLength) +
            " --clobber=yes verbosity=0")
    print("             ", comm)
    args = shlex.split(comm)
    try:
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running tool gennoise (NOISE):")
        print(comm)
        chdir(cwd)
        rmtree(tmpDir)
        raise
    print("GENNOISESPEC: .......................................END")
    print("Noise file ", noiseFile, " has been successfully created\n")
    chdir(cwd)
