"""
RESOLUTION CURVES for single/pairs of pulses

python getEresolCurves_manySeps.py

!!!!!!! Things to review evry time it is run (they can change):
        =======================================================
        ** Check resultsDir (nodetSP in case secondaries are not detected)
                            (detSP in case secondaries are also detected
                            gainScake in case results are for gain scale curves)
        ** Check library (LONG, SHORT) & filenames

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
from subprocess import check_call, check_output, STDOUT
from astropy.io import fits
from astropy.io import ascii


# ----GLOBAL VARIABLES -------------
PreBufferSize = 1000
# tstartPulse1 = int(PreBufferSize)
separation = '40000'  # if unique separation
cwd = os.getcwd()
tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"

XMLdir = os.environ["SIXTE"] + "/" + "share/sixte/instruments/athena/1469mm_xifu"
XMLfile = XMLdir + "/" + "xifu_detector_hex_baseline.xml"
#calibLibs = (0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)  # only for searching for libray intervals


def getEresolCurves(pixName, labelLib, samprate, mono1EkeV, mono2EkeV, reconMethod, filterMeth,
                    nsamples, pulseLength, nSimPulses, fdomain, scaleFactor, samplesUp,
                    samplesDown, nSgms, detMethod, tstartPulse1, tstartPulse2,
                    nSimPulsesLib, coeffsFile, libTmpl, resultsDir, sepsStr):
    """
    :param pixName: Extension name for FITS pixel definition file (SPA*, LPA1*, LPA2*, LPA3*)
    :param labelLib: Label identifying the library ( multilib, multilibOF, fixedlib1 )
    :param samprate: Samprate value with respect to baseline of 156250 Hz (base, half)
    :param mono1EkeV: Monochromatic energy (keV) of input primary simulated pulse
    :param mono2EkeV: Monochromatic energy (keV) of input secondary simulated pulse
    :param reconMethod: Energy reconstruction Method (OPTFILT, WEIGHT, WEIGHTN, I2R, I2RALL, I2RNOL)
    :param filterMeth: Optimal Filtering Method (F0 or B0)
    :param nsamples: noise length samples
    :param pulseLength: pulse length in samples
    :param nSimPulses: number of simulated pulses (to locate filenames to be analysed)
    :param nSimPulsesLib: number of simulated pulses for the library
    :param fdomain: Optimal Filtering domain (F(req) or T(ime))
    :param scaleFactor: Param scaleFactor for detection
    :param samplesUp: Initial Param samplesUp for detection
    :param samplesDown: Initial Param samplesDown for detection
    :param nSgms: Initial Param nSgms for detection
    :param detMethod: Detection method (AD or A1)
    :param tstartPulse1: Start sample for first pulse (PreBufferSize)
    :param tstartPulse2:
                    if /=0 Start sample for second pulse
                    if = -1 modified to (tstartPulse1 + separation)
                    if = 0 && tstartPulse1 /=0 => no secondary pulse
    :param nSimPulsesLib: number of simulated pulses for the library
    #:param coeffs: Polynomial coefficientes (a0 + a1*x + a2*x^2 + a3*x^3 +...) for fit of Erecons vs. Ecalib in keV
    :param coeffsFile: file with coefficients of polynomial fit to gain curves from polyfit2bias.R
    :param libTmpl: label to identify the library regarding the templates used (LONG or SHORT)
    :param resultsDir: directory for resulting evt and json files (from .../PAIRS/eresol+pixName ), tipically
            'nodetSP', 'detSP' or 'gainScale' or ''
    :param sepsStr: blank spaces separated list of pulses separations
    :return: file with energy resolutions for the input pairs of pulses
    """

    TRIGG = ""
    global XMLfile, PreBufferSize, separation
    # --- Define some initial values and conversions ----

    # samprate
    smprtStr = ""
    if samprate == 'samprate2':
        smprtStr = "_samprate2"
    XMLfile = XMLdir + "/" + "xifu_detector_hex_baseline" + smprtStr + ".xml"


    # pulses start
    if tstartPulse1 > 0:
        TRIGG = "_NTRIG"

    # pulses energies
    mono1EeV = float(mono1EkeV) * 1000.
    mono2EeV = float(mono2EkeV) * 1000.

    #pulses types
    tessim = "tessim" + pixName
    classAries = ["primaries", "secondaries", "all"]
    if tstartPulse2 == 0:
        classAries = ["all"]

    # separations
    if sepsStr == "":
        sepsStr = separation

    # data space and reconstruction method
    space = "ADC"
    if reconMethod in ("I2R", "I2RALL", "I2RNOL", "I2RFITTED"):
        space = reconMethod.lstrip('I2')

    # libraries
    OFLib = "no"
    if "OF" in labelLib:
        OFLib = "yes"

    # -- LIB & NOISE & SIMS & RESULTS dirs and files ----------
    simDir = cwd + "/PAIRS/" + tessim

    outDir = cwd + "/PAIRS/eresol" + pixName + "/" + resultsDir
    simSIXTEdir = "/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE"
    noiseDir = simSIXTEdir + "/NOISE/" + tessim
    noiseFile = noiseDir + "/noise" + str(nsamples) + "samples_" + tessim + "_B0_" + space + smprtStr + ".fits"
    libDirRoot = simSIXTEdir + "/LIBRARIES/" + tessim
    libDir = libDirRoot + "/GLOBAL/" + space + "/"
    if libTmpl == "SHORT":
        libTmpl = "_SHORT"
    else:
        libTmpl = ""

    if 'multilib' in labelLib:
        libFile = libDir + "/libraryMultiE_GLOBAL_PL" + str(nsamples) + "_" + str(nSimPulsesLib) + "p" + smprtStr + ".fits"
        filtEeV = 1000.  # eV
    elif 'fixedlib' in labelLib:  # fixedlib1,...
        fixedEkeV = labelLib.replace("OF", "")[8:]
        filtEeV = float(fixedEkeV) * 1E3  # eV
        if tstartPulse1 == 0: #detection to be performed -> require different models
            libFile = libDir + "/libraryMultiE_GLOBAL_PL" + str(nsamples) + "_" + str(nSimPulsesLib) + "p" + smprtStr \
                      + ".fits"
        else:
            libFile = libDir + "/library" + fixedEkeV + "keV_PL" + str(nsamples) + "_" + str(nSimPulsesLib) + "p" + \
                      smprtStr + ".fits"
    libFile = libFile.replace(".fits", libTmpl+".fits")

    root = ''.join([str(nSimPulses), 'p_SIRENA', str(nsamples), '_pL', str(pulseLength), '_',
                    mono1EkeV, 'keV_', mono2EkeV, 'keV_', detMethod, "_", str(filterMeth),
                    str(fdomain), '_', str(labelLib), '_', str(reconMethod), smprtStr])
    if mono2EkeV == 0:
        root = ''.join([str(nSimPulses), 'p_SIRENA', str(nsamples), '_pL', str(pulseLength), '_',
                        mono1EkeV, 'keV_', detMethod, "_", str(filterMeth), str(fdomain), '_',
                        str(labelLib), '_', str(reconMethod), smprtStr])
    if tstartPulse1 > 0:
        root += TRIGG
    eresolFile = "eresol_" + root + ".json"
    eresolFile = eresolFile.replace(".json", libTmpl+".json")

    ofnoise = "NSD"
    if 'OPTFILTNM' in reconMethod:
        ofnoise = "WEIGHTM"
        reconMethod = "OPTFILT"

    # locate coefficients in table
    # -----------------------------
    alias = labelLib + "_" + detMethod + "_" + reconMethod + smprtStr
    if 'fixedlib' in labelLib and 'OF' not in labelLib:
        alias = labelLib + "OF_" + detMethod + "_" + reconMethod + smprtStr
    print("Using alias=", alias)
    os.chdir(outDir)

    # -- Read gain scale coefficients if file is provided --- otherwise we are in the process of creating
    # the gain scale curve
    coeffsDict = dict()
    if coeffsFile:
        codata = ascii.read(coeffsFile, guess=False, format='basic')
        # codata[1] : row 2
        # codata[1][1]: row 2, col 2
        for i in range(0, len(codata)):
            #  METHOD   ALIAS  a0  a1  a2  a3  a4
            coeffsDict[codata[i][1]] = (codata[i][2], codata[i][3], codata[i][4], codata[i][5], codata[i][6])
    else:
        simDir += "/gainScale/"

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
        if tstartPulse1 > 0 and tstartPulse2 == -1:  # calculate tstartPulse2 using separation
            tstartPulse2 = tstartPulse1 + int(sep)
        inFile = simDir + "/sep" + sep12 + "sam_" + str(nSimPulses) + "p_" + mono1EkeV + "keV_" \
            + mono2EkeV + "keV" + smprtStr + "_jitter.fits"

        evtFile = "events_sep" + sep12 + "sam_" + root + "_jitter.fits"
        evtFile = evtFile.replace("_jitter.fits", libTmpl + "_jitter.fits")
        print("=============================================")
        print("Using file: ", inFile)
        print("Using library: ", libFile)
        print("Using noisefile: ", noiseFile)
        print("Setting evtFile: ", evtFile)
        print("Setting eresolFile: ", eresolFile)
        print("=============================================")


        # -- SIRENA processing -----
        if os.path.isfile(evtFile):
            # if events already exist get number nSgms used to get events
            print("Event file", evtFile, " already DOES exist: recalculating FWHMs...")
            try:
                comm = "fkeyprint infile=" + evtFile + "+0 keynam=HISTORY"
                args = shlex.split(comm)
                nSgmsStr = check_output(args)
                # print(nSgmsStr)
                # look for ocurrence of nSgms (before, nSgms, after) and take "after"
                nSgms = nSgmsStr.partition("nSgms")  # as in "HISTORY P19 nSgms = 0"
                # as "after" has a lot of lines, separate by newline and take first one
                nSgms = nSgms[2].splitlines()[0]  # = 0
                # remove "=" and get numerical value
                nSgms = nSgms.split()[1]
                # print(nSgms)
            except:
                print("Error reading pre-existing evtFile:", evtFile)
                print(comm)
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise
        else:
            ndetpulses = 0
            print("\nRunning SIRENA for detection & reconstruction")
            print("filtEev=", filtEeV)
            print("XMLFile=", XMLfile)
            comm = ("tesreconstruction Recordfile=" + inFile + " TesEventFile=" + evtFile +
                    " Rcmethod='SIRENA'" + " PulseLength=" + str(pulseLength) +
                    " LibraryFile=" + libFile + " scaleFactor=" + str(scaleFactor) +
                    " samplesUp=" + str(samplesUp) + " nSgms=" + str(nSgms) +
                    " samplesDown=" + str(samplesDown) + " mode=1 NoiseFile=" + noiseFile +
                    " OFLib=" + OFLib + " FilterDomain=" + fdomain + " detectionMode=" + detMethod +
                    " FilterMethod=" + filterMeth + " clobber=yes intermediate=0 " +
                    "EventListSize=1000" + " EnergyMethod=" + reconMethod +
                    " tstartPulse1=" + str(tstartPulse1) + " tstartPulse2=" +
                    str(tstartPulse2) + " OFNoise=" + ofnoise + " XMLFile=" + XMLfile +
                    " filtEeV=" + str(filtEeV))
                # sys.exit()
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
            nrows = evtf[1].header["NAXIS2"]
            nTrigPulses = evtf[1].header["NETTOT"]
            evtf.close()

            print("Checking DETECTION............................")
            assert nrows > 0, "Empty evt file (%s): nrows=0 " % evtFile

            try:
                comm = "fstatistic infile=" + evtFile + " colname='GRADE1' rows='-' minval=0"
                print("Running: ", comm)
                args = shlex.split(comm)
                check_call(args, stdout=open(os.devnull, 'wb'))  # >/dev/null" does not work
                comm = "pget fstatistic numb"
                args = shlex.split(comm)
                ndetpulses = int(check_output(args))
            except:
                print("Error checking number of detected pulses in evtfile:", evtFile)
                print(comm)
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise
            print('Reconstruction finished => sim/det: ", nTrigPulses, "/", ndetpulses')

        # -----------------------------------------
        # EVENT file processing to calculate FWHM
        # -----------------------------------------
        rootEvt = os.path.splitext(evtFile)[0]

        for aries in classAries:  # PRIMARIES / SECONDARIES / ALL
            print("Working with:", aries, "\n")
            if "aries" in aries:  # PRIMARIES // SECONDARIES
                evt = rootEvt + "_" + aries + "_jitter.fits"
                if os.path.isfile(evt):
                    os.remove(evt)
            elif aries == "all":
                evt = evtFile

            # use only rows with GRADE1>0 (non initially truncated)
            try:
                monoEeV = mono1EeV
                if aries == "primaries":
                    comm = "fselect infile=" + evtFile + " outfile=" + evt + " expr='#ROW%2!=0 && GRADE1>0' clobber=yes"
                elif aries == "secondaries":
                    comm = "fselect infile=" + evtFile + " outfile=" + evt + " expr='#ROW%2==0 && GRADE1>0' clobber=yes"
                    monoEeV = mono2EeV
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

            # ---- if calibration parameters are provided, calculate also corrected energies
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
            fwhmErecons = '{0:0.5f}'.format(fwhm)  # FWHM for reconstructed events in eV
            # EBIAS = <Erecons> - monoEeV
            bias = SIGNALmean * 1000. - monoEeV
            biasErecons = '{0:0.5f}'.format(bias)
            print("SIGNALmean=", SIGNALmean)
            print("monoEeV=", monoEeV)
            print("biasErecons=", biasErecons, "\n")

            # calculate corrected energies if polyfit coeffs are provided and previous fwhm is not NaN (some NULL Eners)
            # if any(a != 0 for a in coeffs) and not math.isnan(fwhm):
            if coeffsFile and not math.isnan(fwhm):
                print("...Calculating corrected energies for ", aries, " pulses")
                ErealKeV = np.zeros(ftab['SIGNAL'].size, dtype=float)
                ie = 0
                for SIGNALKeV in ftab['SIGNAL']:
                    # read fitting coeffs taken from polyfit2Bias.R (a0, a1, a2, a3)
                    #  as in y = a0 + a1*x + a2*x^2 + a3*x^3
                    # where y=E_reconstructed and x=Ecalibration (keV)

                    coeffs = coeffsDict[alias]
                    # print("SIGNAL=", SIGNALKeV, "coeffs=", coeffs)
                    npCoeffs = np.array(coeffs)
                    # print("npCoeffs=",npCoeffs)
                    npCoeffs[0] -= SIGNALKeV  # subtract "y" (0 = a0 + a1*x + a2*x^2 + ... - y)
                    npCoeffsRev = npCoeffs[::-1]  # reversed to say fit with poly1d definition
                    polyfit = np.poly1d(npCoeffsRev)
                    # get real root (value of Ereal for a given Erecons )
                    r = np.roots(polyfit)
                    # real && >0 roots
                    # print("r=", r)
                    rreal = r.real[abs(r.imag) < 1e-5]
                    # print("rreal=", rreal)
                    rrealpos = rreal[rreal > 0]
                    # print("rrealpos=", rrealpos)
                    # closest root
                    rclosest = min(enumerate(rrealpos), key=lambda x: abs(x[1]-SIGNALKeV))[1]  # (idx,value)
                    # print("For:", alias, " Recon energy=", rclosest, "npCoeffs=", npCoeffs)
                    ErealKeV[ie] = rclosest
                    ie += 1
                    # if aries == "secondaries":
                        #print("For:", alias, " Recon energy SEC=", SIGNALKeV, " Real energy SEC=", rclosest)
                        #print(rclosest)

                ErealKeVsigma = ErealKeV.std()
                fwhm = ErealKeVsigma * 2.35 * 1000.
                fwhm_err = fwhm / math.sqrt(2 * nrows - 2)
                bias = ErealKeV.mean() * 1000. - monoEeV
                fwhmEreal = '{0:0.5f}'.format(fwhm)  # FWHM for reconstructed PRIMARIES in eV
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

    parser = argparse.ArgumentParser(description='Get Energy resolution for single or pairs of pulses',
                                     prog='getEresolCurves_manySeps')

    parser.add_argument('--pixName', required=True,
                        help='Extension name in the FITS pixel definition file (SPA*, LPA1*, LPA2*, LPA3*)')
    parser.add_argument('--lib', required=True,
                        help='Label identifying the library (monolib, multilib, multilibOF, fixedlib, fixedlib2,...)')
    parser.add_argument('--libTmpl', required=True, choices=['SHORT','LONG'],
                        help='Label identifying the library for the templates: SHORT or LONG')
    parser.add_argument('--samprate', required=True, choices=['', 'samprate2'], help="baseline, half_baseline")
    parser.add_argument('--monoEnergy1', help='Monochromatic energy (keV) of input primary simulated pulses',
                        required=True)
    parser.add_argument('--monoEnergy2', help='Monochromatic energy (keV) of input secondary simulated pulse',
                        default=0)
    parser.add_argument('--reconMethod', default='OPTFILT',
                        choices=['OPTFILT', 'OPTFILTNM', 'WEIGHT', 'WEIGHTN', 'I2R', 'I2RALL', 'I2RNOL', 'I2RFITTED'],
                        help='Energy reconstruction Method (OPTFILT, OPTFILTNM, WEIGHT, WEIGHTN, I2R, I2RALL, '
                             'I2RNOL, I2RFITTED)')
    parser.add_argument('--filter', choices=['F0', 'B0'], default='F0',
                        help='Optimal Filtering Method (F0, B0) [default %(default)s]')
    parser.add_argument('--nsamples', type=int, help='noise length samples', required=True)
    parser.add_argument('--pulseLength', type=int, help='pulse length for reconstruction', required=True)
    parser.add_argument('--nSimPulses', type=int, help='pulses in simulations (to locate filenames)', required=True)
    parser.add_argument('--nSimPulsesLib', type=int, help='pulses used to build the library', required=True)
    parser.add_argument('--fdomain', choices=['F', 'T'], default='T',
                        help='Optimal Filtering domain (F(req) or T(ime)) [default %(default)s]')
    parser.add_argument('--scaleFactor', type=float, default=0.0,
                        help='Param scaleFactor for detection [default %(default)s]')
    parser.add_argument('--samplesUp', type=int, default=0,
                        help='Param samplesUp for detection [default %(default)s]')
    parser.add_argument('--samplesDown', type=int, default=0,
                        help='Param samplesUp for detection [default %(default)s]')
    parser.add_argument('--nSgms', type=float, help='Param nSgms for detection [default %(default)s]', default=0)
    parser.add_argument('--tstartPulse1', type=int, default=0,
                        help='Start sample for first pulse [default %(default)s]')
    parser.add_argument('--tstartPulse2', type=int, default=0,
                        help='Start sample for second pulse [default %(default)s]')
    parser.add_argument('--coeffsFile', help='file with polynomial fit coefficients from R script polyfit2bias.R',
                        default='')
    parser.add_argument('--detMethod', default='AD', choices=['AD', 'A1'])
    parser.add_argument('--detectSP', help='are Secondaries to be detected (1:yes; 0:no)', default=1, type=int)
    parser.add_argument('--resultsDir', help='directory for resulting evt and json files (from '
                        '.../PAIRS/eresol+pixName) tipically "nodetSP", "detSP" or "gainScale" or ""',
                        required=True)
    parser.add_argument('--separations', help='spaced list of separations between Prim & Sec pulses', default="",
                        nargs='+')

    inargs = parser.parse_args()

    # print("array=",inargs.array)
    if (inargs.nSgms == 0) and (inargs.tstartPulse1 == 0 and inargs.tstartPulse2 == 0):
        print("Start sample of pulses or detection parameters must be provided")
        sys.exit()

    getEresolCurves(pixName=inargs.pixName, labelLib=inargs.lib, samprate=inargs.samprate,
                    mono1EkeV=inargs.monoEnergy1, mono2EkeV=inargs.monoEnergy2,
                    reconMethod=inargs.reconMethod, filterMeth=inargs.filter,
                    nsamples=inargs.nsamples, pulseLength=inargs.pulseLength,
                    nSimPulses=inargs.nSimPulses, fdomain=inargs.fdomain,
                    scaleFactor=inargs.scaleFactor, samplesUp=inargs.samplesUp,
                    samplesDown=inargs.samplesDown, nSgms=inargs.nSgms,
                    detMethod=inargs.detMethod, tstartPulse1=inargs.tstartPulse1,
                    tstartPulse2=inargs.tstartPulse2, nSimPulsesLib=inargs.nSimPulsesLib,
                    coeffsFile=inargs.coeffsFile, libTmpl=inargs.libTmpl,
                    resultsDir=inargs.resultsDir, sepsStr=inargs.separations)
