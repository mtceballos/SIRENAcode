#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 12:10:40 2020

@author: ceballos
"""

from os import getcwd, path, chdir, remove
import shlex
from astropy.io import fits
from subprocess import check_call, STDOUT
from shutil import rmtree
import tempfile
from sixtevars import XIFUSIMinst, XMLsixte, sampfreqs, sampids, sampStrs, \
    preBufferPulses, separations
from fitsVerify import fitsVerify
from updateHISTORY import updateHISTORY
from rmLastAndFirst import rmLastAndFirst


tmpDir = tempfile.mkdtemp()

def simulPairs(pixName, monoEkeV1, monoEkeV2, acbias, samprate, jitter,
               noise, bbfb, nSimPulses, pulseLength, dcmt, sepsStr):
    """
    :param pixName: Extension name in the FITS pixel definition file
                     (SPA*, LPA1*, LPA2*, LPA3*)
    :param monoEkeV1: Monochromatic energy (keV) of input first
                    simulated pulses
    :param monoEkeV2: Monochromatic energy (keV) of input second
                    simulated pulses
    :param acbias: Operating Current (AC if acbias=yes or DC if acbias=no)
    :param samprate: Samprate value with respect to baseline of 156250 Hz:
                    "" (baseline), "samprate2" (half_baseline)
    :param jitter: jitter option ("" for no_jitter and "jitter" for jitter)
    :param noise: noise option ("" for noisy and "nonoise" for no_noise)
    :param bbfb: ("") for dobbfb=n or ("bbfb") for dobbfb=y
    :param nSimPulses: number of pulses to be simulated
    :param pulseLength: length of pulses to calculate trigger size
    :param sepsStr: separations between pulses
    :return: files with simulated PAIRS
    """
    cwd = getcwd()
    print("cwd=", cwd)
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA"
    ERESOLdir = EURECAdir + "/ERESOL"
    PAIRSdir = ERESOLdir + "/PAIRS"
    preBufferSize = preBufferPulses

    # samprate
    smprtStr = ""
    idxsmp = sampids.index(samprate)
    smprtStr = sampStrs[idxsmp]
    samplingrate = sampfreqs[idxsmp]

    XMLxifusim = XIFUSIMinst + "1pix_nobbfb.xml"

    # deal with record separations from definitions above
    recordSeparation = separations[idxsmp]

    # added to solve floating point inaccu. due to sampling rate
    # (Christian's mail 31/03/2017)
    tstart = 0.5 / samplingrate

    # jitter
    jitterStr = ""
    jitterStrPix = ""
    offset = ""
    if jitter == "jitter":
        jitterStr = "_jitter_dcmt" + str(dcmt)
        jitterStrPix = "_jitter"
        offset = " offset=-1"

    # noise
    noiseStr = ""
    simnoise = " simnoise=y"
    if noise == "nonoise":
        noiseStr = "_nonoise"
        simnoise = " simnoise=n"

    # bbfb
    bbfbStr = ""
    assert not(jitter == "" and bbfb == "bbfb"),\
        "Error: BBFB will produce jitter"

    if bbfb == "bbfb":
        bbfbStr = "_bbfb"
        XMLxifusim = XIFUSIMinst + "1pix_bbfb.xml"

    xifusim = "xifusim" + pixName
    SIMFILESdir = PAIRSdir + "/" + xifusim

    for sepA in sepsStr:
        sep12 = int(sepA)
        triggerSizeTC = preBufferSize + sep12 + recordSeparation + 1000
        triggerSizeXF = preBufferSize + sep12 + pulseLength + 1000
        triggerSuppXF = triggerSizeXF - preBufferSize
        triggerSizeTC = int(triggerSizeTC)
        diffth = 60

        # calculate sim time to have at least nSimPulses pulses:
        #  simTime= (nSimPulses/2)recs * (triggerSize sam/rec)/(samprate sam/s)
        simTime = nSimPulses / 2. * triggerSizeTC / samplingrate
        simTime = '{0:0.0f}'.format(simTime)

        # for piximpact:
        root0 = ("sep" + sepA + "sam_" + simTime + "s_" + monoEkeV1 +
                 "keV_" + monoEkeV2 + "keV" + smprtStr + jitterStrPix)
        # for fits:
        root = ("sep" + sepA + "sam_" + str(nSimPulses) + "p_" +
                monoEkeV1 + "keV_" + monoEkeV2 + "keV" +
                smprtStr + jitterStr + noiseStr + bbfbStr)
        pixFile = (PAIRSdir + "/PIXIMPACT/" + root0 + "_trSz" +
                   str(triggerSizeTC) + ".piximpact")
        fitsFile = SIMFILESdir + "/" + root + ".fits"
        print("-------------------------------------------\n")
        print("Simulating ", fitsFile, "\n")
        print("  with ", pixFile, "\n")
        print("-------------------------------------------\n")

        if not path.isfile(pixFile):
            comm = ("tesconstpileup PixImpList=" + pixFile +
                    " XMLFile=" + XMLsixte + " timezero=" + str(tstart) +
                    " tstop=" + str(simTime) + " energy=" + monoEkeV1 +
                    " energy2=" + monoEkeV2 + offset + " clobber=yes" +
                    " pulseDistance=" + str(sep12) + " sample_freq=" +
                    str(samplingrate) + " TriggerSize=" +
                    str(triggerSizeTC))
            print("\n##### Runing tesconstpileup #########")
            print(comm, "\n")
            try:
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except RuntimeError:
                print("Error running TESCONSTPILEUP for piximpact list",
                      " generation with:", comm)
                chdir(cwd)
                rmtree(tmpDir)
                raise
            # continue  # to simulate only piximpact files

        if path.isfile(fitsFile):
            # verify existing file
            numerrs, numwrns = fitsVerify(fitsFile)
            if numerrs > 0:
                print("numerrs = ", numerrs, " for ", fitsFile,
                      ": repeating simulation")
                remove(fitsFile)

        if not path.isfile(fitsFile):
            commxifusim = ("xifusim PixImpList=" + pixFile +
                           " Streamfile=" + fitsFile +
                           " tstart=0." + " tstop=" + simTime +
                           " trig_reclength=" + triggerSizeXF +
                           " trig_n_pre=" + preBufferSize +
                           " trig_thresh=" + str(diffth) +
                           " trig_n_suppress=" + triggerSuppXF +
                           " acbias=" + acbias +
                           " sample_rate=" + str(samplingrate) + simnoise +
                           " decimate_factor=" + str(dcmt) +
                           " XMLfilename=" + XMLxifusim)

            print("\n##### Runing xifusim #########")
            print(commxifusim, "\n")
            try:
                args = shlex.split(commxifusim)
                check_call(args, stderr=STDOUT)
            except RuntimeError:
                print("Error running xifusim for data simulation")
                chdir(cwd)
                rmtree(tmpDir)
                raise
            # continue
            # update HISTORY  for tesconstpileup run in header[0]
            updateHISTORY(fitsFile, commxifusim)

            # rm first (and LAST) record and update NETTOT
            fsim = fits.open(fitsFile)
            nrows = fsim[1].header["NAXIS2"]
            assert nrows > 1, "xifusim failed: just one huge row present!"
            fsim.close()
            # continue
        try:
            print("Removing first & last row, just in case, ",
                  "and updating NETTOT")
            rmLastAndFirst(fitsFile, 1)
        except RuntimeError:
            print("Error running FTOOLS to remove initial & last rows in ",
                  fitsFile)
            chdir(cwd)
            rmtree(tmpDir)
            raise
        # Modify EXTNAME if not EXTNAME=RECORDS in empty trigger file
        # and Add TRIGGSZ to xifusim simulated file
        fsim = fits.open(fitsFile, mode='update')
        fsim[0].header["TRIGGSZ"] = triggerSizeXF
        fsim[0].header["DELTAT"] = dcmt * fsim[1].header["DELTA_T"]
        fsim[1].header["TRIGGSZ"] = triggerSizeXF
        fsim[1].header["DELTAT"] = dcmt * fsim[1].header["DELTA_T"]
        oldExtName = fsim[1].header["EXTNAME"]
        if not oldExtName == "RECORDS":
            print("Changing name...")
            fsim[1].header['EXTNAME'] = 'RECORDS'
            fsim.close()

    # some cleaning before exiting the function
    chdir(cwd)
