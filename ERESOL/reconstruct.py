#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 12:30:32 2020

@author: ceballos

Reconstruct events ~energy (no calibration)
"""

from os import path, chdir, stat, remove
import shlex
from astropy.io import fits
from subprocess import check_call, STDOUT
from shutil import rmtree
from sixtevars import XMLfll, XMLsixte, sampids, sampStrs, separations, samplesUps,\
    samplesDowns, nSigmss, scaleFactor
import tempfile


tmpDir = tempfile.mkdtemp()

def reconstruct(pixName, labelLib, samprate, jitter, dcmt, noise, bbfb, Lc,
                mono1EkeV, mono2EkeV, reconMethod, filterLength,
                nsamples, pulseLength, nSimPulses, fdomain, detMethod,
                tstartPulse1, tstartPulse2, nSimPulsesLib, coeffsFile,
                libTmpl, simDir, outDir, detSP, pB, Ifit, LbT, s0, lags,
                filterct, B0, sepsStr):

    """
    :param pixName: Extension name for FITS pixel definition file
                    (SPA*, LPA1*, LPA2*, LPA3*)
    :param labelLib: Label identifying the library
                    ( multilib, multilibOF, fixedlib1 )
    :param samprate: Samprate value with respect to baseline of 156250 Hz:
                    "" (baseline), "samprate2" (half_baseline),
                    samprate4 (quarter baseline)
    :param jitter: jitter option ("" for no_jitter and "jitter" for jitter or
                                  "jitter_M82" for M82 special case)
    :param dcmt: decimation factor for xifusim jitter
    :param bbfb: ("") for dobbfb=n or ("bbfb") for dobbfb=y
    :param Lc: ("" for L=Lcrit; otherwise, L/Lcrit)
    :param noise: simulations done with ("") or without ("nonoise") noise
    :param mono1EkeV: Monochromatic energy (keV)
                        of input primary simulated pulse
    :param mono2EkeV: Monochromatic energy (keV)
                        of input secondary simulated pulse
    :param reconMethod: Energy reconstruction Method
                        (OPTFILT, WEIGHT, WEIGHTN, I2R, I2RALL, I2RNOL,...)
    :param filterLength: if set to 0 (default) filters are taken
                        from library in BASE2, otherwise filterStartegy=FIXED
    :param nsamples: noise length samples
    :param pulseLength: pulse length in samples
    :param nSimPulses: number of simulated pulses
                        (to locate filenames to be analysed)
    :param nSimPulsesLib: number of simulated pulses for the library
    :param fdomain: Optimal Filtering domain (F(req) or T(ime))
    :param detMethod: Detection method (AD or STC)
    :param tstartPulse1: Start sample for first pulse (PreBufferSize)
    :param tstartPulse2:
                    if /=0 Start sample for second pulse
                    if = -1 modified to (tstartPulse1 + separation)
                    if = 0 && tstartPulse1 /=0 => no secondary pulse
    :param nSimPulsesLib: number of simulated pulses for the library
    :param libTmpl: label to identify the library regarding
                    the templates used (LONG or SHORT)
    :param simDir: directory for simulated files
    :param outDir: directory for resulting evt and json files
    :param detSP: 1 secondary pulses will be detected (default), 0 otherwise
    :param pB: preBuffer value for optimal filters
    :param Ifit: fitted constant for I2RFITTED conversion
    :param LbT: Time (s) to average baseline to be subtracted
    :param s0: Optimal Filters' SUM should be '0'?: s0=0 (NO); s0=1 (YES)
    :param lags: Do parabola fit? lags=1 (YES), lags=0 (NO)
    :param filterct: Filters central part have been replaced by ct value?
                    or use derived filters from largest 8192 filter (fit)))
    :param B0: if B0>0 processing will be done by B0 (baseline subtraction)
                Otherwise, F0 (bin 0 freq.)
    :param sepsStr: blank spaces separated list of pulses separations
    :return: file with energy resolutions for the input pairs of pulses
    """

    # --- Define some initial values and conversions ----
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA/"
    simSIXTEdir = EURECAdir + "testHarness/simulations/SIXTE"

    # samprate
    smprtStr = ""
    idxsmp = sampids.index(samprate)
    smprtStr = sampStrs[idxsmp]

    # separations
    separation = separations[idxsmp]
    if sepsStr == "":
        sepsStr = [separation]

    # detection related params
    samplesUp = samplesUps[idxsmp]
    samplesDown = samplesDowns[idxsmp]
    nSgms = nSigmss[idxsmp]

    # jitter
    jitterStr = ""
    if jitter == "jitter":
        jitterStr = "_jitter"
    if jitter == "jitter" and dcmt > 1:
        jitterStr = "_jitter_dcmt" + str(dcmt)
    if jitter == "jitter_M82":
        jitterStr = "_jitter_M82"

    # noise
    noiseStr = ""
    if noise == "nonoise":
        noiseStr = "_nonoise"

    # bbfb
    bbfbStr = ""
    XML = XMLsixte
    if bbfb == "bbfb":
        bbfbStr = "_bbfb"
    elif bbfb == "fll":
        bbfbStr = "_fll"
        XML = XMLfll

    if "NewPar" in bbfb or "040" in bbfb:
        bbfbStr = "_" + bbfb
        bbfb = "bbfb"

    # optimal filters' preBuffer
    pBStr = ""
    if pB > 0:
        pBStr = "_pB" + str(pB)

    # I2RFITTED Ifit value
    IfitStr = ""
    if Ifit > 0:
        IfitStr = "_Ifit_" + str(int(Ifit))
    elif Ifit < 0:
        IfitStr = "_Ifit_m" + str(int(abs(Ifit)))

    # baseline average samples
    LbTStr = ""
    #if float(LbT) > 0:
    #    LbTStr = "_LbT" + str(LbT)
    # optimal filter's SUM
    s0Str = ""
    s0Param = ""
    if s0 == 1:
        s0Str = "_Sum0Filt"
        s0Param = " Sum0Filt=1"
    # lags
    lagsStr = ""
    if lags == 0:
        lagsStr = "_nolags"

    # optimal filters' (constant central part or derived-from-8192 filters)
    ctStr = ""
    if filterct:
        ctStr = "_" + filterct
    # optimal filters B0 (baseline subtraction) or F0 (bin f=0) processing
    B0str = ""
    FMparam = " FilterMethod=F0"
    if B0 > 0:
        B0str = "_B0-" + str(B0)
        FMparam = " FilterMethod=B0"

    # Lcrit
    LcStr = ""
    if not Lc == "":
        LcStr = "_" + str(Lc) + "Lc"

    # pulses start
    TRIGG = ""
    if tstartPulse1 > 0:
        TRIGG = "NTRIG"
    else:
        TRIGG = detMethod

    # pulses types
    xifusim = "xifusim" + pixName

    # data space and reconstruction method
    space = "ADC"
    if "I2R" in reconMethod:
        space = reconMethod.lstrip('I2')
        space = space.split('NM')[0]

    ofnoise = "NSD"
    reconMethod2 = ""
    if "NM" in reconMethod:
        ofnoise = "WEIGHTM"
        reconMethod2 = "NM" + reconMethod.split('NM')[1]
        reconMethod = reconMethod.split('NM')[0]

    # if "WEIGHT" in reconMethod:
    #    lags = 0
    # libraries
    OFLib = "no"
    OFstrategy = ""
    # noiseDir = simSIXTEdir + "/NOISE/" + xifusim
    # noiseFile = (noiseDir + "/noise" + str(nsamples) + "samples_" +
    #             xifusim + "_B0_" + space + smprtStr +
    #             jitterStr + bbfbStr + LcStr + ".fits")
    # noiseParam = " NoiseFile=" + noiseFile
    if "OF" in labelLib:
        OFLib = "yes"
        OFstrategy = " OFStrategy=FIXED OFLength=" + str(filterLength)
        noiseParam = ""
        # if "WEIGHT" in reconMethod:
        #    noiseParam = " NoiseFile=" + noiseFile

    # -- LIB & NOISE & SIMS & RESULTS dirs and files ----------
    # simDir = ERESOLdir + "PAIRS/" + xifusim
    # resultsDir = resultsDir.rstrip('\\')
    # resultsDir = resultsDir.lstrip('\\')
    # if ("gainScale" in resultsDir) or ("base" in resultsDir):
    #    simDir += "/" + resultsDir + "/"

    # outDir = ERESOLdir + "PAIRS/eresol" + pixName + "/" + resultsDir
    libDirRoot = simSIXTEdir + "/LIBRARIES/" + xifusim
    libDir = libDirRoot + "/GLOBAL/" + space + "/"
    if libTmpl == "SHORT":
        libTmpl = "_SHORT"
    else:
        libTmpl = ""

    if 'multilib' in labelLib:
        libFile = (libDir + "/libraryMultiE_GLOBAL_PL" + str(nsamples) +
                   "_" + str(nSimPulsesLib) + "p" + smprtStr +
                   jitterStr + noiseStr + bbfbStr + LcStr + ".fits")
        filtEeV = 6000.  # eV
    elif 'fixedlib' in labelLib:  # fixedlib1,...
        fixedEkeV = labelLib.replace("OF", "")[8:]
        filtEeV = float(fixedEkeV) * 1E3  # eV
        # detection to be performed -> require different models:
        if tstartPulse1 == 0 and detMethod == "AD":
            libFile = (libDir + "/libraryMultiE_GLOBAL_PL" + str(nsamples) +
                       "_" + str(nSimPulsesLib) + "p" + smprtStr + B0str +
                       jitterStr + noiseStr + bbfbStr + LcStr + pBStr + IfitStr +
                       ctStr + ".fits")
        else:
            libFile = (libDir + "/library" + fixedEkeV + "keV_PL" +
                       str(nsamples) + "_" + str(nSimPulsesLib) + "p" +
                       smprtStr + B0str + jitterStr + noiseStr + bbfbStr +
                       LcStr + pBStr + IfitStr + ctStr + ".fits")
    libFile = libFile.replace(".fits", libTmpl + ".fits")

    root = ''.join([str(nSimPulses), 'p_SIRENA', str(nsamples), '_pL',
                    str(pulseLength), '_', mono1EkeV, 'keV_', mono2EkeV,
                    'keV_', TRIGG, "_", str(fdomain), '_',
                    str(labelLib), '_', str(reconMethod), str(filterLength),
                    reconMethod2, B0str, pBStr, IfitStr, LbTStr, smprtStr, jitterStr,
                    noiseStr, bbfbStr, LcStr, s0Str, lagsStr, ctStr])
    if mono2EkeV == "0":
        root = ''.join([str(nSimPulses), 'p_SIRENA', str(nsamples), '_pL',
                        str(pulseLength), '_', mono1EkeV, 'keV_', TRIGG, "_",
                        str(fdomain), '_', str(labelLib),
                        '_', str(reconMethod), str(filterLength),
                        reconMethod2, B0str, pBStr, LbTStr, IfitStr, smprtStr,
                        jitterStr, noiseStr, bbfbStr, LcStr, s0Str,
                        lagsStr, ctStr])

    eresolFile = "eresol_" + root + ".json"
    eresolFile = eresolFile.replace(".json", libTmpl+".json")

    chdir(outDir)

    # ------------------------------------
    # --- Process input data files -------
    # ------------------------------------
    #
    for sep12 in sepsStr:
        sep = float(sep12)
        # calculate tstartPulse2 using separation:
        if tstartPulse1 > 0 and tstartPulse2 == -1:
            tstartPulse2 = tstartPulse1 + int(sep)

        inFile = (simDir + "/sep" + sep12 + "sam_" + str(nSimPulses) +
                  "p_" + str(mono1EkeV) + "keV_" + str(mono2EkeV) + "keV" +
                  smprtStr + jitterStr + noiseStr + bbfbStr + LcStr +
                  ".fits")
        if mono2EkeV == "0":  # for single Pulses
            inFile = (simDir + "/sep" + sep12 + "sam_" + str(nSimPulses) +
                      "p_" + str(mono1EkeV) + "keV" +
                      smprtStr + jitterStr + noiseStr + bbfbStr + LcStr +
                      ".fits")

        evtFile = "events_sep" + sep12 + "sam_" + root + ".fits"
        evtFile = evtFile.replace(
                jitterStr + noiseStr + bbfbStr + LcStr + lagsStr + ".fits",
                libTmpl + jitterStr + noiseStr + bbfbStr + LcStr + lagsStr +
                ".fits")

        print("=============================================")
        print("RECONSTRUCTING ENERGIES.....................")
        print("Working in:", outDir)
        print("Using file: ", inFile)
        print("Using library: ", libFile)
        # print("Using noisefile: ", noiseFile)
        print("Setting evtFile: ", evtFile)
        print("=============================================")

        # -- SIRENA processing -----

        if path.isfile(evtFile):
            if stat(evtFile).st_size == 0:
                remove(evtFile)
                print("Event file", evtFile,
                      " already DOES exist but size 0: removing...")
            else:
                print("Event file", evtFile,
                      " already DOES exist: recalculating FWHMs...")

        if not path.isfile(evtFile):
            ndetpulses = 0
            print("\nRunning SIRENA for detection & reconstruction")
            print("filtEev=", filtEeV)
            print("XMLFile=", XML)
            comm = ("tesreconstruction Recordfile=" + inFile +
                    " TesEventFile=" + evtFile +
                    " Rcmethod='SIRENA'" +
                    " PulseLength=" + str(pulseLength) +
                    " LibraryFile=" + libFile +
                    " scaleFactor=" + str(scaleFactor) +
                    " samplesUp=" + str(samplesUp) + " nSgms=" + str(nSgms) +
                    " samplesDown=" + str(samplesDown) +
                    " opmode=1 " + noiseParam +
                    " OFLib=" + OFLib + " FilterDomain=" + fdomain +
                    " detectionMode=" + detMethod + " detectSP=" + str(detSP) +
                    " clobber=yes" +
                    " EventListSize=1000" + " EnergyMethod=" + reconMethod +
                    " LagsOrNot=" + str(lags) +
                    " tstartPulse1=" + str(tstartPulse1) +
                    " tstartPulse2=" + str(tstartPulse2) +
                    " OFNoise=" + ofnoise + " XMLFile=" + XML +
                    " filtEeV=" + str(filtEeV) + OFstrategy +
                    " preBuffer=" + str(pB) +
                    " Ifit=" + str(Ifit) +
                    " LbT=" + str(LbT) + s0Param + FMparam)
            try:
                print(comm)
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except RuntimeError:
                print("Error running SIRENA for detection & reconstruction:")
                print(comm)
                rmtree(tmpDir)
                raise

            # Check that detection went well
            evtf = fits.open(evtFile)
            nTrigPulses = evtf[1].header["NETTOT"]
            ndetpulses = evtf[1].header["NAXIS2"]
            print("Checking DETECTION............................")
            assert ndetpulses > 0, "Empty evt file (%s): nrows=0 " % evtFile
            print("Reconstruction finished => sim/det: ",
                  nTrigPulses, "/", ndetpulses)
            evtf.close()

    return smprtStr, jitterStr, noiseStr, bbfbStr, LcStr, pBStr, IfitStr, LbTStr,\
        s0Str, lagsStr, ctStr, B0str, evtFile, eresolFile
