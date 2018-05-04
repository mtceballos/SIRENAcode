"""
RESOLUTION CURVES for single/pairs of pulses

python getEresolCurves_manySeps.py

!!!!!!! Things to review evry time it is run (they can change):
        =======================================================
        ** Check resultsDir (nodetSP in case secondaries are not detected)
                            (detSP in case secondaries are also detected
                           gainScale in case results are for gain scale curves)

        ** SINCE NEWSIMULATIONS (Feb2018): library is always the jiiter one
        although applied to non-jitter simulations(bagplots)
"""

# ----IMPORT MODULES --------------
from __future__ import print_function
import os
import math
import shlex
import shutil
import sys
import tempfile
import numpy as np
import json
import auxpy
from subprocess import check_call, STDOUT
from astropy.io import fits
import xml.etree.ElementTree as ET


# ----GLOBAL VARIABLES -------------
PreBufferSize = 1000
separation = '40000'  # if unique separation
cwd = os.getcwd()
tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"


def getEresolCurves(pixName, labelLib, samprate, jitter, stoch, noise,
                    mono1EkeV, mono2EkeV, reconMethod, filterMeth,
                    filterLength, nsamples, pulseLength, nSimPulses,
                    fdomain, scaleFactor, samplesUp, samplesDown,
                    nSgms, detMethod, tstartPulse1, tstartPulse2,
                    nSimPulsesLib, coeffsFile, libTmpl, resultsDir,
                    sepsStr):
    """
    :param pixName: Extension name for FITS pixel definition file
                    (SPA*, LPA1*, LPA2*, LPA3*)
    :param labelLib: Label identifying the library
                    ( multilib, multilibOF, fixedlib1 )
    :param samprate: Samprate value with respect to baseline of 156250 Hz:
                    "" (baseline), "samprate2" (half_baseline)
    :param jitter: jitter option ("" for no_jitter and "jitter" for jitter)
    :param stoch: stochastic option ("" for no_stoch and "stoch" for stochast)
    :param noise: simulations done with ("") or without ("nonoise") noise
    :param mono1EkeV: Monochromatic energy (keV)
                        of input primary simulated pulse
    :param mono2EkeV: Monochromatic energy (keV)
                        of input secondary simulated pulse
    :param reconMethod: Energy reconstruction Method
                        (OPTFILT, WEIGHT, WEIGHTN, I2R, I2RALL, I2RNOL)
    :param filterMeth: Optimal Filtering Method (F0 or B0)
    :param filterLength: if set to 0 (default) filters are taken
                        from library in BASE2, otherwise filterStartegy=FIXED
    :param nsamples: noise length samples
    :param pulseLength: pulse length in samples
    :param nSimPulses: number of simulated pulses
                        (to locate filenames to be analysed)
    :param nSimPulsesLib: number of simulated pulses for the library
    :param fdomain: Optimal Filtering domain (F(req) or T(ime))
    :param scaleFactor: Param scaleFactor for detection
    :param samplesUp: Initial Param samplesUp for detection
    :param samplesDown: Initial Param samplesDown for detection
    :param nSgms: Initial Param nSgms for detection
    :param detMethod: Detection method (AD or STC)
    :param tstartPulse1: Start sample for first pulse (PreBufferSize)
    :param tstartPulse2:
                    if /=0 Start sample for second pulse
                    if = -1 modified to (tstartPulse1 + separation)
                    if = 0 && tstartPulse1 /=0 => no secondary pulse
    :param nSimPulsesLib: number of simulated pulses for the library
    :param coeffsFile: file with coefficients of polynomial fit to
                    gain curves from polyfit2bias.R
    :param libTmpl: label to identify the library regarding
                    the templates used (LONG or SHORT)
    :param resultsDir: directory for resulting evt and json files
                    (from .../PAIRS/eresol+pixName ), tipically
                    'nodetSP', 'detSP' or 'gainScale' or ''
    :param sepsStr: blank spaces separated list of pulses separations
    :return: file with energy resolutions for the input pairs of pulses
    """

    TRIGG = ""
    global PreBufferSize, separation
    # --- Define some initial values and conversions ----
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA/"
    ERESOLdir = EURECAdir + "/ERESOL/"

    # samprate
    smprtStr = ""
    if samprate == 'samprate2':
        smprtStr = "_samprate2"

    XMLfile = (ERESOLdir + "xifu_detector_hex_baselineNEWgrades" +
               smprtStr + ".xml")
    XMLtree = ET.parse(XMLfile)
    XMLroot = XMLtree.getroot()
    for key in XMLroot.findall('grading'):
        num = key.get('num')
        if num == "1":
            invalid = key.get('pre')
            Hres = key.get('post')
    print("Hres,invalid=", Hres, invalid)

    # jitter
    jitterStr = ""
    if jitter == "jitter":
        jitterStr = "_jitter"

    # noise
    noiseStr = ""
    if noise == "nonoise":
        noiseStr = "_nonoise"

    # stochastic
    stochStr = ""
    if stoch == "stoch":
        stochStr = "_stoch"

    # pulses start
    if tstartPulse1 > 0:
        TRIGG = "NTRIG"
    else:
        TRIGG = detMethod

    # pulses energies
    mono1EeV = float(mono1EkeV) * 1000.
    mono2EeV = float(mono2EkeV) * 1000.

    # pulses types
    tessim = "tessim" + pixName
    classAries = ["primaries", "secondaries", "all"]
    if tstartPulse2 == 0:
        classAries = ["all"]

    # separations
    if sepsStr == "":
        sepsStr = [separation]

    # data space and reconstruction method
    space = "ADC"
    if reconMethod in ("I2R", "I2RALL", "I2RNOL", "I2RFITTED"):
        space = reconMethod.lstrip('I2')

    # libraries
    OFLib = "no"
    OFstrategy = ""
    if "OF" in labelLib:
        OFLib = "yes"
        OFstrategy = " OFStrategy=FIXED OFLength=" + str(filterLength)

    # -- LIB & NOISE & SIMS & RESULTS dirs and files ----------
    simDir = cwd + "/PAIRS/" + tessim
    if resultsDir == "gainScale":
        simDir += "/gainScale/"

    outDir = cwd + "/PAIRS/eresol" + pixName + "/" + resultsDir
    simSIXTEdir = EURECAdir + "/testHarness/simulations/SIXTE"
    noiseDir = simSIXTEdir + "/NOISE/" + tessim
    noiseFile = (noiseDir + "/noise" + str(nsamples) + "samples_" +
                 tessim + "_B0_" + space + smprtStr +
                 jitterStr + stochStr + ".fits")
    libDirRoot = simSIXTEdir + "/LIBRARIES/" + tessim
    libDir = libDirRoot + "/GLOBAL/" + space + "/"
    if libTmpl == "SHORT":
        libTmpl = "_SHORT"
    else:
        libTmpl = ""

    if 'multilib' in labelLib:
        libFile = (libDir + "/libraryMultiE_GLOBAL_PL" + str(nsamples) +
                   "_" + str(nSimPulsesLib) + "p" + smprtStr +
                   jitterStr + noiseStr + stochStr + ".fits")
        filtEeV = 1000.  # eV
    elif 'fixedlib' in labelLib:  # fixedlib1,...
        fixedEkeV = labelLib.replace("OF", "")[8:]
        filtEeV = float(fixedEkeV) * 1E3  # eV
        # detection to be performed -> require different models:
        if tstartPulse1 == 0:
            libFile = (libDir + "/libraryMultiE_GLOBAL_PL" + str(nsamples) +
                       "_" + str(nSimPulsesLib) + "p" + smprtStr +
                       jitterStr + noiseStr + stochStr + ".fits")
        else:
            libFile = (libDir + "/library" + fixedEkeV + "keV_PL" +
                       str(nsamples) + "_" + str(nSimPulsesLib) + "p" +
                       smprtStr + jitterStr + noiseStr + stochStr + ".fits")
    libFile = libFile.replace(".fits", libTmpl+".fits")

    root = ''.join([str(nSimPulses), 'p_SIRENA', str(nsamples), '_pL',
                    str(pulseLength), '_', mono1EkeV, 'keV_', mono2EkeV,
                    'keV_', TRIGG, "_", str(filterMeth), str(fdomain), '_',
                    str(labelLib), '_', str(reconMethod), str(filterLength),
                    smprtStr, jitterStr, noiseStr, stochStr])
    if mono2EkeV == "0":
        root = ''.join([str(nSimPulses), 'p_SIRENA', str(nsamples), '_pL',
                        str(pulseLength), '_', mono1EkeV, 'keV_', TRIGG, "_",
                        str(filterMeth), str(fdomain), '_', str(labelLib),
                        '_', str(reconMethod), str(filterLength), smprtStr,
                        jitterStr, noiseStr, stochStr])

    eresolFile = "eresol_" + root + ".json"
    eresolFile = eresolFile.replace(".json", libTmpl+".json")

    ofnoise = "NSD"
    if 'OPTFILTNM' in reconMethod:
        ofnoise = "WEIGHTM"
        reconMethod = "OPTFILT"

    # locate coefficients in table
    # -----------------------------
    alias = (detMethod + "_" + labelLib + "_" + reconMethod +
             str(pulseLength) + smprtStr + jitterStr + noiseStr + stochStr)
    if 'fixedlib' in labelLib and 'OF' not in labelLib:
        alias = (detMethod + "_" + labelLib + "OF_" + "_" + reconMethod +
                 str(pulseLength) + smprtStr + jitterStr + noiseStr + stochStr)
    print("Using alias=", alias)
    os.chdir(outDir)

    # ------------------------------------
    # --- Process input data files -------
    # ------------------------------------
    #  define/initialize structure to save all output data
    data = []  # list of dictionaries (one per pulses separation)
    for sep12 in sepsStr:
        dictSep = dict()
        dictBiasErecons = dict()
        dictBiasEreal = dict()
        dictFwhmErecons = dict()
        dictFwhmEreal = dict()
        dictFwhmErealErr = dict()

        sep = float(sep12)
        # calculate tstartPulse2 using separation:
        if tstartPulse1 > 0 and tstartPulse2 == -1:
            tstartPulse2 = tstartPulse1 + int(sep)

        inFile = (simDir + "/sep" + sep12 + "sam_" + str(nSimPulses) +
                  "p_" + mono1EkeV + "keV_" + mono2EkeV + "keV" +
                  smprtStr + jitterStr + noiseStr + stochStr + ".fits")
        if mono2EkeV == "0":  # for single Pulses
            inFile = (simDir + "/sep" + sep12 + "sam_" + str(nSimPulses) +
                      "p_" + mono1EkeV + "keV" +
                      smprtStr + jitterStr + noiseStr + stochStr + ".fits")

        evtFile = "events_sep" + sep12 + "sam_" + root + ".fits"
        evtFile = evtFile.replace(jitterStr + noiseStr + stochStr + ".fits",
                        libTmpl + jitterStr + noiseStr + stochStr + ".fits")
        print("=============================================")
        print("Using file: ", inFile)
        print("Using library: ", libFile)
        print("Using noisefile: ", noiseFile)
        print("Setting evtFile: ", evtFile)
        print("Setting eresolFile: ", eresolFile)
        print("=============================================")

        # -- SIRENA processing -----

        if os.path.isfile(evtFile):
            if os.stat(evtFile).st_size == 0:
                os.remove(evtFile)
                print("Event file", evtFile,
                      " already DOES exist but size 0: removing...")
            else:
                print("Event file", evtFile,
                      " already DOES exist: recalculating FWHMs...")

        if not os.path.isfile(evtFile):
            ndetpulses = 0
            print("\nRunning SIRENA for detection & reconstruction")
            print("filtEev=", filtEeV)
            print("XMLFile=", XMLfile)
            comm = ("tesreconstruction Recordfile=" + inFile +
                    " TesEventFile=" + evtFile +
                    " Rcmethod='SIRENA'" +
                    " PulseLength=" + str(pulseLength) +
                    " LibraryFile=" + libFile +
                    " scaleFactor=" + str(scaleFactor) +
                    " samplesUp=" + str(samplesUp) + " nSgms=" + str(nSgms) +
                    " samplesDown=" + str(samplesDown) +
                    " mode=1 NoiseFile=" + noiseFile +
                    " OFLib=" + OFLib + " FilterDomain=" + fdomain +
                    " detectionMode=" + detMethod +
                    " FilterMethod=" + filterMeth + " clobber=yes" +
                    " EventListSize=1000" + " EnergyMethod=" + reconMethod +
                    " LagsOrNot=1" + " tstartPulse1=" + str(tstartPulse1) +
                    " tstartPulse2=" + str(tstartPulse2) +
                    " OFNoise=" + ofnoise + " XMLFile=" + XMLfile +
                    " filtEeV=" + str(filtEeV)) + OFstrategy
            try:
                print(comm)
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except:
                print("Error running SIRENA for detection & reconstruction:")
                print(comm)
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise

            evtf = fits.open(evtFile)
            nTrigPulses = evtf[1].header["NETTOT"]
            ndetpulses = evtf[1].header["NAXIS2"]
            print("Checking DETECTION............................")
            assert ndetpulses > 0, "Empty evt file (%s): nrows=0 " % evtFile
            print("Reconstruction finished => sim/det: ",
                  nTrigPulses, "/", ndetpulses)
            evtf.close()

        # --------------------------------------------------------------
        # EVENT file processing to calculate FWHM (only for Hres pulses)
        # ---------------------------------------------------------------
        rootEvt = os.path.splitext(evtFile)[0]

        for aries in classAries:  # PRIMARIES / SECONDARIES / ALL
            print("Working with:", aries, "\n")
            if "aries" in aries:  # PRIMARIES // SECONDARIES
                evt = rootEvt + "_" + aries + ".fits"
                if os.path.isfile(evt):
                    os.remove(evt)
            elif aries == "all":
                # evt = evtFile
                evt = rootEvt + "_HR.fits"

            # use only rows with GRADE1==Hres && GRADE2>invalid (Hres pulses)
            try:
                monoEeV = mono1EeV
                if aries == "primaries":
                    comm = "fselect infile=" + evtFile + " outfile=" + evt + \
                           " expr='#ROW%2!=0 && GRADE1==" + str(Hres) +\
                           " && GRADE2>" + str(invalid) + "' clobber=yes"
                elif aries == "secondaries":
                    comm = "fselect infile=" + evtFile + " outfile=" + evt + \
                           " expr='#ROW%2==0 && GRADE1==" + str(Hres) +\
                           " && GRADE2>" + str(invalid) + "' clobber=yes"
                    monoEeV = mono2EeV
                else:
                    comm = "fselect infile=" + evtFile + " outfile=" + evt +\
                           " expr='GRADE1==" + str(Hres) + \
                           " && GRADE2>" + str(invalid) + "' clobber=yes"
                args = shlex.split(comm)
                check_call(args, stdout=open(os.devnull, 'wb'))
            except:
                print("Error selecting PRIM/SEC events in evtfile ", evt)
                print(comm)
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise
            f = fits.open(evt, memmap=True)
            nrows = f[1].header["NAXIS2"]
            assert nrows > 0, "Empty evt file (%s): nrows=0 " % evt

            # ---- if calibration parameters are provided,
            #              calculate also corrected energies
            # ---- if not, calculate only reconstructed energies FWMH and BIAS
            # FWHM & BIAS for Erecons
            # read Erecons (SIGNAL) column in numpy array (in keV)
            ftab = f[1].data
            SIGNALmean = ftab['SIGNAL'].mean()
            SIGNALsigma = ftab['SIGNAL'].std()

            fwhmEreal = 0
            fwhmEreal_err = 0
            biasEreal = 0

            fwhm = SIGNALsigma * 2.35 * 1000.
            # FWHM for reconstructed events in eV:
            fwhmErecons = '{0:0.5f}'.format(fwhm)
            # EBIAS = <Erecons> - monoEeV
            bias = SIGNALmean * 1000. - monoEeV
            biasErecons = '{0:0.5f}'.format(bias)
            # print("SIGNALmean=", SIGNALmean)
            # print("monoEeV=", monoEeV)
            # print("biasErecons=", biasErecons, "\n")

            # calculate corrected energies if polyfit coeffs are provided and
            # previous fwhm is not NaN (some NULL Eners)
            # FWHM only calculated for those pulses reconstructed with HRes
            # if any(a != 0 for a in coeffs) and not math.isnan(fwhm):

            if coeffsFile and not math.isnan(fwhm):
                EreconKeV = np.array(ftab['SIGNAL'])
                print("...Calculating corrected energies for ",
                      aries, " pulses")
                ErealKeV = auxpy.enerToCalEner(EreconKeV, coeffsFile, alias)

                ErealKeVsigma = ErealKeV.std()
                fwhm = ErealKeVsigma * 2.35 * 1000.
                fwhm_err = fwhm / math.sqrt(2 * nrows - 2)
                bias = ErealKeV.mean() * 1000. - monoEeV
                # FWHM for reconstructed PRIMARIES in eV:
                fwhmEreal = '{0:0.5f}'.format(fwhm)
                fwhmEreal_err = '{0:0.5f}'.format(fwhm_err)
                biasEreal = '{0:0.5f}'.format(bias)
                # print("bias=", bias, "\n")
                # print("...biasEreal=", biasEreal, "\n")

            f.close()
            del f[1].data
            dictFwhmErecons[aries] = fwhmErecons
            dictFwhmEreal[aries] = fwhmEreal
            dictFwhmErealErr[aries] = fwhmEreal_err
            dictBiasErecons[aries] = biasErecons
            dictBiasEreal[aries] = biasEreal
            # END of PRIM & SEC & ALL

        # write data for given separation to a dictionary
        dictSep["separation"] = sep12
        dictSep["fwhmErecons"] = dictFwhmErecons
        dictSep["biasErecons"] = dictBiasErecons
        dictSep["fwhmEreal"] = dictFwhmEreal
        dictSep["biasEreal"] = dictBiasEreal
        dictSep["fwhmEreal_err"] = dictFwhmErealErr

        # save dictionary to list of dictionaries "data"
        data.append(dictSep)

    # Dump total data to JSON file
    feresol = open(eresolFile, 'w')

    json.dump(obj=data, fp=feresol, sort_keys=True)
    feresol.close()
    os.chdir(cwd)
    shutil.rmtree(tmpDir)
#
# --------- Read input parameters and PROCESS simulations -----------------
#


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
            description='Get Energy resolution for single or pairs of pulses',
            prog='getEresolCurves_manySeps')

    parser.add_argument('--pixName', required=True,
                        help='Extension name in the FITS pixel definition file\
                        (SPA*, LPA1*, LPA2*, LPA3*)')
    parser.add_argument('--lib', required=True,
                        help='Label identifying the library (monolib,\
                        multilib, multilibOF, fixedlib, fixedlib2,...)')
    parser.add_argument('--libTmpl', choices=['SHORT', 'LONG'], default='LONG',
                        help='Label identifying the library for the templates:\
                        SHORT or LONG')
    parser.add_argument('--samprate', default="", choices=['', 'samprate2'],
                        help="baseline, half_baseline")
    parser.add_argument('--jitter', default="", choices=['', 'jitter'],
                        help="no jitter, jitter")
    parser.add_argument('--noise', default="", choices=['', 'nonoise'],
                        help="noisy, nonoise")
    parser.add_argument('--stoch', default="", choices=['', 'stoch'],
                        help="non-stochastic, stochastic")
    parser.add_argument('--monoEnergy1', help='Monochromatic energy (keV) of\
                        input primary simulated pulses', required=True)
    parser.add_argument('--monoEnergy2', help='Monochromatic energy (keV) of\
                        input secondary simulated pulse', default=0)
    parser.add_argument('--reconMethod', default='OPTFILT',
                        choices=['OPTFILT', 'OPTFILTNM', 'WEIGHT', 'WEIGHTN',
                                 'I2R', 'I2RALL', 'I2RNOL', 'I2RFITTED'],
                        help='Energy reconstruction Method (OPTFILT,\
                        OPTFILTNM, WEIGHT, WEIGHTN, I2R, I2RALL,\
                        I2RNOL, I2RFITTED)')
    parser.add_argument('--filter', choices=['F0', 'B0'], default='F0',
                        help='Optimal Filtering Method (F0, B0)\
                        [default %(default)s]')
    parser.add_argument('--filterLength', type=int,
                        help='filter length for reconstruction', default=0)
    parser.add_argument('--nsamples', type=int,
                        help='noise length samples', required=True)
    parser.add_argument('--pulseLength', type=int,
                        help='pulse length for reconstruction', required=True)
    parser.add_argument('--nSimPulses', type=int, required=True,
                        help='pulses in simulations (to locate filenames)')
    parser.add_argument('--nSimPulsesLib', type=int,
                        help='pulses used to build the library', required=True)
    parser.add_argument('--fdomain', choices=['F', 'T'], default='T',
                        help='Optimal Filtering domain (F(req) or T(ime))\
                        [default %(default)s]')
    parser.add_argument('--scaleFactor', type=float, default=0.0,
                        help='Param scaleFactor for detection\
                        [default %(default)s]')
    parser.add_argument('--samplesUp', type=int, default=0,
                        help='Param samplesUp for detection\
                        [default %(default)s]')
    parser.add_argument('--samplesDown', type=int, default=0,
                        help='Param samplesUp for detection\
                        [default %(default)s]')
    parser.add_argument('--nSgms', type=float, help='Param nSgms\
                        for detection [default %(default)s]', default=0)
    parser.add_argument('--tstartPulse1', type=int, default=0,
                        help='Start sample for first pulse\
                        [default %(default)s]')
    parser.add_argument('--tstartPulse2', type=int, default=0,
                        help='Start sample for second pulse\
                        [default %(default)s]')
    parser.add_argument('--coeffsFile', default='',
                        help='file with polynomial fit coefficients from\
                        R script polyfit2bias.R')
    parser.add_argument('--detMethod', default='AD', choices=['AD', 'STC'])
    parser.add_argument('--detectSP', default=1, type=int,
                        help='are Secondaries to be detected (1:yes; 0:no)')
    parser.add_argument('--resultsDir', default="",
                        help='directory for resulting evt and json files\
                        (from ..../PAIRS/eresol+pixName) tipically\
                        "nodetSP", "detSP" or "gainScale" or ""')
    parser.add_argument('--separations', nargs='+',
                        help='spaced list of separations between\
                        Prim & Sec pulses', default="")

    inargs = parser.parse_args()

    # print("array=",inargs.array)
    if (inargs.nSgms == 0) and (inargs.tstartPulse1 == 0
       and inargs.tstartPulse2 == 0):
        print("Start sample of pulses or detection parameters \
              must be provided")
        sys.exit()

    getEresolCurves(pixName=inargs.pixName, labelLib=inargs.lib,
                    samprate=inargs.samprate, jitter=inargs.jitter,
                    noise=inargs.noise, stoch=inargs.stoch,
                    mono1EkeV=inargs.monoEnergy1, mono2EkeV=inargs.monoEnergy2,
                    reconMethod=inargs.reconMethod, filterMeth=inargs.filter,
                    filterLength=inargs.filterLength,
                    nsamples=inargs.nsamples, pulseLength=inargs.pulseLength,
                    nSimPulses=inargs.nSimPulses, fdomain=inargs.fdomain,
                    scaleFactor=inargs.scaleFactor, samplesUp=inargs.samplesUp,
                    samplesDown=inargs.samplesDown, nSgms=inargs.nSgms,
                    detMethod=inargs.detMethod,
                    tstartPulse1=inargs.tstartPulse1,
                    tstartPulse2=inargs.tstartPulse2,
                    nSimPulsesLib=inargs.nSimPulsesLib,
                    coeffsFile=inargs.coeffsFile, libTmpl=inargs.libTmpl,
                    resultsDir=inargs.resultsDir, sepsStr=inargs.separations)
