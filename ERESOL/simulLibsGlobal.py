#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 12:26:42 2020

@author: ceballos
"""

from os import getcwd, path, chdir
import shlex
from astropy.io import fits
from subprocess import check_call, STDOUT
from shutil import rmtree
from time import gmtime, strftime
from sixtevars import XIFUSIMinst, XMLsixte, sampfreqs, sampids, sampStrs, \
    preBufferPulses, separations, samplesUps, samplesDowns, nSigmss
import tempfile
from rmLastAndFirst import rmLastAndFirst
from updateHISTORY import updateHISTORY

tmpDir = tempfile.mkdtemp()

def simulLibsGlobal(pixName, space, samprate, jitter, noise, bbfb,
                    pulseLength, libEnergies, largeFilter, nsamples,
                    nSimPulses, acbias, tstartPulse1All,
                    createLib, noiseMat, nintervals, weightMat, dcmt):
    """
    :type pixName: str
    :param pixName: Extension name in FITS file pixel definition
                    (SPA*, LPA1*, LPA2*, LPA3*)
    :type space: str
    :param space: Input Data Space (ADC, R, RALL, RNOL, RFITTED)
    :type samprate: str
    :param samprate: Samprate value with respect to baseline of 156250 Hz:
                     "" (baseline), "samprate2" (half_baseline)
    :type jitter: str
    :param jitter: jitter option ("" for no_jitter and "jitter" for jitter)
    :type noise: str
    :param noise: simulations done with ("") or without ("nonoise") noise
    :type bbfb: str
    :param bbfb: ("") for dobbfb=n or ("bbfb") for dobbfb=y
    :type pulseLength: int
    :param pulseLength: Pulse length
    :type libEnergies: str
    :param libEnergies: list of calibration energies (keV)
    :type largeFilter: int
    :param largeFilter: size of the extra-large filter to avoid record-length
                        effects (samples)
    :type nsamples: int
    :param nsamples: noise length samples
    :type nSimPulses: int
    :param nSimPulses: number of pulses in simulated files (uses to
                       create/choose files and to name library)
    :type acbias: str
    :param acbias: Operating Current: AC (acbias=yes) or DC (acbias=no)
    :type tstartPulse1All: int
    :param tstartPulse1All: list of initial sample for 1st pulses in each
                            record. If empty, just simulate calibration files
                            **  0 if detection is to be performed
    :type createLib: int
    :param createLib: if library should be created (1) or script should
                      just create simulated files (0)
    :type noiseMat: str
    :param noiseMat: should the Noise matrices HDU be created? (yes/no)
    :type nintervals: int
    :param nintervals: number of intervals if WEIGHTM noise used (inst.of NSD)
    :type weigthMat: str
    :param weightMat: should the Weight matrices HDU be created? (yes/no)
    :param dcmt: decimation factor for jitter in xifusim
    :return: simulated calibration pulses pulses && Global library from them

    """

    cwd = getcwd()
    preBufferSize = preBufferPulses
    xifusim = "xifusim" + pixName

    # samprate
    smprtStr = ""
    idxsmp = sampids.index(samprate)
    smprtStr = sampStrs[idxsmp]
    samplingrate = sampfreqs[idxsmp]

    # deal with separations from definitions above
    separation = int(separations[idxsmp])

    # Sigmas and Samples and scaleFactor for Detection
    maxFilterLength = pulseLength
    if largeFilter > 0:
        maxFilterLength = largeFilter

    samplesUp = samplesUps[idxsmp]
    samplesDown = samplesDowns[idxsmp]
    nSgms = nSigmss[idxsmp]

    XMLxifusim = XIFUSIMinst + "1pix_nobbfb.xml"

    # jitter
    jitterStr = ""
    jitterStrPix = ""
    offset = ""
    if jitter == "jitter":
        jitterStr = "_jitter"
        jitterStrPix = "_jitter"
        offset = " offset=-1"
    if dcmt > 1:
        jitterStr = "_jitter_dcmt" + str(dcmt)

    # noise
    noiseStr = ""
    simnoise = " simnoise=y"
    if noise == "nonoise":
        noiseStr = "_nonoise"
        simnoise = " simnoise=n"

    # noiseMat
    if noiseMat == "yes" and nintervals > 0:
        noiseMatStr = "_" + str(nintervals) + "int"
    # bbfb
    bbfbStr = ""
    assert not(jitter == "" and bbfb == "bbfb"),\
        "Error: BBFB will produce jitter"

    if bbfb == "bbfb":
        bbfbStr = "_bbfb"
        XMLxifusim = XIFUSIMinst + "1pix_bbfb.xml"

    # -- CALIB/SIM/NOISE/LIB/XML dirs & files --
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA"
    simSIXTEdir = EURECAdir + "/testHarness/simulations/SIXTE"
    PIXIMPACTdir = simSIXTEdir + "/LIBRARIES/PIXIMPACT"
    SIMFILESdir = simSIXTEdir + "/LIBRARIES/" + xifusim
    NOISEdir = simSIXTEdir + "/NOISE/" + xifusim
    libDir = SIMFILESdir + "/GLOBAL/" + space

    noiseFile = (NOISEdir + "/noise" + str(nsamples) + "samples_" + xifusim +
                 "_B0_" + space + smprtStr + jitterStr +
                 bbfbStr + noiseMatStr + ".fits")
    libFile = (libDir + "/libraryMultiE_GLOBAL_PL" + str(pulseLength) +
               "_" + str(nSimPulses) + "p" + smprtStr + jitterStr +
               noiseStr + bbfbStr + noiseMatStr + ".fits")
    evttmpFile = tempfile.NamedTemporaryFile()

    # added to solve floating point inaccuracies due to sampling rate
    # (Christian's mail 31/03/2017):
    tstart = 0.5/samplingrate

    # Calibration energies and Tstarts of pulses
    tstartPulse1 = dict(zip(libEnergies, tstartPulse1All))

    # calculate simTime so that nSimPulses are simulated (1 pulses per record)
    simTime = nSimPulses * (preBufferSize + separation) / float(samplingrate)
    tstop = tstart + simTime
    simTime = '{0:0.0f}'.format(simTime)
    tstop = '{0:0.0f}'.format(tstop)

    # Trigger sizes in tesconstpileup
    triggerSizeTC = preBufferSize + separation + separation + 1000
    # and xifusim ___|\_______________ooo
    triggerSizeXF = preBufferSize + maxFilterLength + 1000
    triggerSuppXF = triggerSizeXF - preBufferSize
    diffth = 60
    #
    # Select processing space
    #
    energyMethod = "NONE"
    if space == "R":
        energyMethod = " EnergyMethod=I2R"
    elif space == "RALL":
        energyMethod = " EnergyMethod=I2RALL"
    elif space == "RNOL":
        energyMethod = " EnergyMethod=I2RNOL"
    elif space == "RFITTED":
        energyMethod = " EnergyMethod=I2RFITTED"
    elif space == "ADC":
        energyMethod = " EnergyMethod=OPTFILT"

    # --- Simulate and Process calibration data files -------
    for monoEkeV in libEnergies:  # keV
        monoEeV = float(monoEkeV) * 1.E3  # eV

        print("=============================================")
        print("Adding monochromatic energy", monoEkeV, "keV")
        print("=============================================")
        print("simTime=", simTime, "\n")
        # print("Temporary event file in: ", evttmpFile.name)
        # simulate SIXTE file with tesconstpileup + xifusim
        root0 = ("mono" + str(monoEkeV) + "_sep" + str(separation) + "_pix1" +
                 "_" + str(nSimPulses) + "p")
        root = root0 + "_" + str(pulseLength)
        # pixFile = PIXIMPACTdir +
        #          "/" + root0 + "_trSz" + str(triggerSizeTC) + ".piximpact"
        pixFile = (PIXIMPACTdir + "/" + root0 + smprtStr + jitterStrPix +
                   ".piximpact")
        simFile = (SIMFILESdir + "/" + root + smprtStr + jitterStr + noiseStr +
                   bbfbStr + ".fits")

        # -- TESCONSTPILEUP: generate impacts for well separated single pulses
        # -- xifusim: simulate well separated pulses --
        if not path.isfile(simFile):
            print("Simulating & triggering pulses with xifusim to ", simFile)
            comm = ("tesconstpileup PixImpList=" + pixFile +
                    " XMLFile=" + XMLsixte + " timezero=" + str(tstart) +
                    " tstop=" + str(tstop) + " energy=" + str(monoEkeV) +
                    " pulseDistance=" + str(separation) + offset +
                    " TriggerSize=" + str(triggerSizeTC) + " clobber=yes" +
                    " sample_freq= " + str(samplingrate))
            print("Generating PIXIMPACT ", pixFile, " running:\n", comm)
            try:
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except RuntimeError:
                print("Error running task for piximpact list generation")
                chdir(cwd)
                rmtree(tmpDir)
                raise

            # continue # to only generate piximpact
            commxifusim = ("xifusim PixImpList=" + pixFile +
                           " Streamfile=" + simFile +
                           " tstart=0 tstop=" + str(simTime) +
                           " trig_reclength=" + triggerSizeXF +
                           " trig_n_pre=" + preBufferSize +
                           " trig_thresh=" + str(diffth) +
                           " trig_n_suppress=" + triggerSuppXF +
                           " acbias=" + acbias +
                           " sample_rate=" + str(samplingrate) +
                           " decimate_factor=" + str(dcmt) + simnoise +
                           " XMLfilename=" + XMLxifusim)

            print("Running xifusim for simulation\n", commxifusim)
            try:
                args = shlex.split(commxifusim)
                check_call(args, stderr=STDOUT)
            except RuntimeError:
                print("Error running xifusim for simulation\n")
                chdir(cwd)
                rmtree(tmpDir)
                raise

            rmLastAndFirst(simFile, 1)

            # update HISTORY in header[0]
            dateTime = strftime("%Y-%m-%d %H:%M:%S", gmtime())
            history = ("Created & Updated by simulLibsGlobal.py on " +
                       dateTime + "with command: " + commxifusim)
            updateHISTORY(simFile, history)

        # print("Antes de evaluar createLib\n")
        if not createLib:
            print("NO Library creation\n")
            continue  # (to just simulate pulses files)

        print("Creating library\n")
        # -- Create/update LIBRARY --
        if tstartPulse1[monoEkeV] == 0:
            # Detecting
            detectStr = (" samplesUp=" + str(samplesUp) +
                         " samplesDown=" + str(samplesDown) +
                         " nSgms=" + str(nSgms) + " detectionMode=STC")

        else:
            # Not detecting
            detectStr = " tstartPulse1=" + str(tstartPulse1[monoEkeV])

        comm = ("tesreconstruction Recordfile=" + simFile +
                " TesEventFile=" + evttmpFile.name +
                " Rcmethod=SIRENA" + detectStr +
                " PulseLength=" + str(pulseLength) +
                " LibraryFile=" + libFile + energyMethod +
                " NoiseFile=" + noiseFile +
                " opmode=0 clobber=yes intermediate=0" +
                " monoenergy=" + str(monoEeV) + " EventListSize=1000" +
                " XMLFile=" + XMLsixte +
                " hduPRECALWN=" + weightMat + " hduPRCLOFWM=" + noiseMat)

        args = shlex.split(comm)
        print("SIRENA reconstruction to add a line to the library, "
              "running command:\n", comm)
        try:
            check_call(args, stderr=STDOUT)
        except RuntimeError:
            print("Error running SIRENA to add new line to library")
            chdir(cwd)
            rmtree(tmpDir)
            raise
        print("SIRENA reconstruction to add a line to the library: done\n",
              comm)
        # -- check that number of reconstructed pulses is the same as the
        # number of simulated pulses
        fsim = fits.open(simFile)
        # number of simulated/triggered pulses:
        nettot = fsim[1].header['NETTOT']
        fsim.close()
        fevt = fits.open(evttmpFile.name)
        ndetpulses = fevt[1].header['NAXIS2']  # number of detected pulses
        fevt.close()
        # assert nettot == ndetpulses, \
        #    "Detection failure: all pulses (%d) should be detected (%d) in %s:
        #    " % (nettot, ndetpulses, simFile)
        print("Simpulses=", nettot, "and detected pulses=",
              ndetpulses, "in ", simFile)

    chdir(cwd)
