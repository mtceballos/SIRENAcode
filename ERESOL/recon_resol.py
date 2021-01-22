"""
RECONSTRUCT multiples files (of pairs) of pulses (at different separations)
AND (optionally) calibrate energies AND (optionally) calculate Energy
resolutions (FWHM)

python recon_resol.py

1) reconstruct energies
2) (optionally) calibrate energies IF input gain scale coefficientes
    are provided
3) (optionally) produce json files with fwhm values pre/post calibration

!!!!!!! Things to review every time it is run (they can change):
        =======================================================
        ** Check resultsDir (nodetSP in case secondaries are not detected,
                           gainScale in case results are for gain scale curves)

        ** SINCE NEWSIMULATIONS (Feb2018): library is always the jitter one
        although applied to non-jitter simulations(bagplots)
"""

# ----IMPORT MODULES --------------
from __future__ import print_function
import os
import math
import shlex
import shutil
# import sys
import tempfile
import json
import sixtevars
from astropy.io import fits
from subprocess import check_call
from reconstruct import reconstruct
from convertEnergies import convertEnergies
import xml.etree.ElementTree as ET


# ----GLOBAL VARIABLES -------------
cwd = os.getcwd()
tmpDir = tempfile.mkdtemp()
print(tmpDir)
# input("stop here")
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"

#
# --------- Read input parameters and PROCESS simulations -----------------
#
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
            description='Reconstruct (and calibrate) (and get ener resol)',
            prog='recon_resol')

    parser.add_argument('--pixName', required=True,
                        help='Extension name in the FITS pixel definition file\
                        (SPA*, LPA1*, LPA2*, LPA3*)')
    parser.add_argument('--labelLib', required=True,
                        help='Label identifying the library (monolib,\
                        multilib, multilibOF, fixedlib, fixedlib2,...)')
    parser.add_argument('--libTmpl', default='LONG',
                        choices=['SHORT', 'LONG'],
                        help='Label identifying the library for the templates:\
                        SHORT or LONG')
    parser.add_argument('--samprate', default="",
                        choices=['', 'samprate2', 'samprate4'],
                        help="baseline, half_baseline, quarter baseline")
    parser.add_argument('--jitter', default="",
                        choices=['', 'jitter', 'jitter_M82_040'],
                        help="no jitter, jitter, M82_040")
    parser.add_argument('--noise', default="",
                        choices=['', 'nonoise'],
                        help="noisy, nonoise")
    parser.add_argument('--decimation', type=int, default=1,
                        help='xifusim decimation factor')
    parser.add_argument('--monoEnergy1', required=True,
                        help='Monochromatic energy (keV) of\
                        input primary simulated pulses')
    parser.add_argument('--monoEnergy2', default="0",
                        help='Monochromatic energy (keV) of\
                        input secondary simulated pulse')
    parser.add_argument('--reconMethod', default='OPTFILT',
                        help='Energy reconstruction Method (OPTFILT,\
                        OPTFILTNMnnnn, WEIGHT, WEIGHTN, I2R, I2RNMnnnn,\
                        I2RALL, I2RNOL, I2RFITTED)')
    parser.add_argument('--filterLength', type=int, default=0,
                        help='filter length for reconstruction')
    parser.add_argument('--nsamples', type=int,
                        help='noise length samples', required=True)
    parser.add_argument('--pulseLength', type=int, required=True,
                        help='pulse length for reconstruction')
    parser.add_argument('--nSimPulses', type=int, required=True,
                        help='pulses in simulations (to locate filenames)')
    parser.add_argument('--nSimPulsesLib', type=int, required=True,
                        help='pulses used to build the library')
    parser.add_argument('--fdomain', default='T',
                        choices=['F', 'T'],
                        help='Optimal Filtering domain (F(req) or T(ime))\
                        [default %(default)s]')
    parser.add_argument('--tstartPulse1', type=int, default=0,
                        help='Start sample for first pulse\
                        [default %(default)s]')
    parser.add_argument('--tstartPulse2', type=int, default=0,
                        help='Start sample for second pulse\
                        [default %(default)s]')
    parser.add_argument('--coeffsFile', default='',
                        help='file with polynomial fit coefficients from\
                        R script polyfit2bias.R')
    parser.add_argument('--detMethod', default='AD',
                        choices=['AD', 'STC'])
    parser.add_argument('--detSP', default=1, type=int,
                        help='are Secondaries to be detected (1:yes; 0:no)')
    parser.add_argument('--resultsDir', default="",
                        help='directory for resulting evt and json files\
                        (from ..../PAIRS/eresol+pixName) tipically\
                        "nodetSP", "detSP" or "gainScale" or ""')
    parser.add_argument('--separations', nargs='+',
                        help='spaced list of separations between\
                        Prim & Sec pulses', default="")
    parser.add_argument('--bbfb', default="",
                        choices=['', 'bbfb', 'bbfb_NewPar', 'bbfb_040',
                                 'bbfb_040_ct', 'fll', '8pix_tdm'],
                                help="dobbfb=n, dobbfb=y")
    parser.add_argument('--Lc', default="",
                        help="Inductance over critical value")
    parser.add_argument('--preBuffer', default=0, type=int,
                        help="preBuffer value for optimal filters")
    parser.add_argument('--Ifit', default=0, type=float,
                        help="Fitted constant for the I2RFITTED expression")
    parser.add_argument('--B0', default=0, type=int,
                        help="B0 (baseline subtraction; B0>0) of \
                        F0 (B0=0) optimal filter")
    parser.add_argument('--LbT', default="0",
                        help="Time (s) to average baseline")
    parser.add_argument('--Sum0Filt', default=0, type=int,
                        help="Optimal Filter SUM shoud be 0? (0=NO; 1=YES)")
    parser.add_argument('--lags', default=1, type=int,
                        help="Do parabola fit if lags=1")
    parser.add_argument('--ct', default="",
                        choices=['', 'ct', 'fit'],
                        help="Filters central part replaced by constant? or\
                        filters derived from largest 8192?")
    inargs = parser.parse_args()

    # print("array=",inargs.array)

    # rename input parameters
    pixName = inargs.pixName
    labelLib = inargs.labelLib
    samprate = inargs.samprate
    jitter = inargs.jitter
    noise = inargs.noise
    bbfb = inargs.bbfb
    Lc = inargs.Lc
    dcmt = inargs.decimation
    mono1EkeV = inargs.monoEnergy1
    mono2EkeV = inargs.monoEnergy2
    reconMethod = inargs.reconMethod
    filterLength = inargs.filterLength
    nsamples = inargs.nsamples
    pulseLength = inargs.pulseLength
    nSimPulses = inargs.nSimPulses
    fdomain = inargs.fdomain
    detMethod = inargs.detMethod
    tstartPulse1 = inargs.tstartPulse1
    tstartPulse2 = inargs.tstartPulse2
    nSimPulsesLib = inargs.nSimPulsesLib
    coeffsFile = inargs.coeffsFile
    libTmpl = inargs.libTmpl
    resultsDir = inargs.resultsDir
    sepsStr = inargs.separations
    detSP = inargs.detSP
    pB = inargs.preBuffer
    Ifit = inargs.Ifit
    LbT = inargs.LbT
    s0 = inargs.Sum0Filt
    lags = inargs.lags
    filterct = inargs.ct
    B0 = inargs.B0

    # general definitions
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA/"
    ERESOLdir = EURECAdir + "/ERESOL/"
    xifusim = "xifusim" + pixName
    simDir = ERESOLdir + "PAIRS/" + xifusim
    resultsDir = resultsDir.rstrip('\\')
    resultsDir = resultsDir.lstrip('\\')
    if ("gainScale" in resultsDir) or ("base" in resultsDir):
        simDir += "/" + resultsDir
    outDir = ERESOLdir + "PAIRS/eresol" + pixName + "/" + resultsDir
    classAries = ["primaries", "secondaries", "all"]
    if tstartPulse2 == 0:
        classAries = ["all"]

    idxsmp = sixtevars.sampids.index(samprate)
    separation = sixtevars.separations[idxsmp]
    if sepsStr == "":
        sepsStr = [separation]

    # read grading info
    XMLtree = ET.parse(sixtevars.XMLsixte)
    if bbfb == "fll":
        XMLtree = ET.parse(sixtevars.XMLfll)
    elif bbfb == "8pix_tdm":
        XMLtree = ET.parse(sixtevars.XMLtdm)
    XMLroot = XMLtree.getroot()
    for key in XMLroot.findall('grading'):
        num = key.get('num')
        if num == "1":
            invalid = key.get('pre')
            Hres = key.get('post')
            factor = sixtevars.sampfreqs[idxsmp]/sixtevars.sampfreqs[0]
            Hres = int(int(Hres) * factor)
            invalid = int(int(invalid) * factor)

    # to be able to use shorter filters or long filter shortened in FWHM curve
    Hres = min(filterLength, pulseLength)

    print("Hres=", Hres)
    print("filterLength=", filterLength)
    print("pulseLength=", pulseLength)

    # 1) Reconstruct energies
    # ------------------------
    smprtStr, jitterStr, noiseStr, bbfbStr, LcStr, pBStr, IfitStr, LbTStr, \
        s0Str, lagsStr, ctStr, B0str, evtFile, eresolFile = \
        reconstruct(pixName, labelLib, samprate, jitter, dcmt,
                          noise, bbfb, Lc, mono1EkeV, mono2EkeV, reconMethod,
                          filterLength, nsamples, pulseLength, nSimPulses,
                          fdomain, detMethod, tstartPulse1, tstartPulse2,
                          nSimPulsesLib, coeffsFile, libTmpl, simDir, outDir,
                          detSP, pB, Ifit, LbT, s0, lags, filterct, B0, sepsStr)

    # 2) Calibrate (AND/OR) extract Energy resolution info to .json files
    # --------------------------------------------------------------------
    # ---- if calibration parameters are provided:
    #              calculate corrected energies
    #              calculate FWHM, BIAS (calib)
    # ---- if not, calculate only FWMH and BIAS (recons)

    #  define/initialize structure to save all output data

    data = []  # list of dictionaries (one per pulses separation)
    for sep12 in sepsStr:
        dictSep = dict()
        dictBiasErecons = dict()
        dictBiasEreal = dict()
        dictFwhmErecons = dict()
        dictFwhmEreal = dict()
        dictFwhmErealErr = dict()

        print("=============================================")
        print("EXTRACTING FWHMs.........................")
        print("Working in:", outDir)
        print("Setting evtFile: ", evtFile)
        print("Setting eresolFile: ", eresolFile)
        print("=============================================")

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
                monoEeV = float(mono1EkeV) * 1000.
                if aries == "primaries":
                    comm = "fselect infile=" + evtFile + " outfile=" + evt + \
                           " expr='#ROW%2!=0 && GRADE1==" + str(Hres) +\
                           " && GRADE2>" + str(invalid) + "' clobber=yes"
                elif aries == "secondaries":
                    comm = "fselect infile=" + evtFile + " outfile=" + evt + \
                           " expr='#ROW%2==0 && GRADE1==" + str(Hres) +\
                           " && GRADE2>" + str(invalid) + "' clobber=yes"
                    monoEeV = float(mono2EkeV) * 1000.
                else:
                    # comm = "fselect infile=" + evtFile + " outfile=" + evt +\
                    #       " expr='GRADE1==" + str(Hres) + \
                    #       " && GRADE2>" + str(invalid) + "' clobber=yes"
                    comm = "fselect infile=" + evtFile + " outfile=" + evt +\
                           " expr='GRADE1>=" + str(Hres) + \
                           " && GRADE2>=" + str(Hres) + "' clobber=yes"
                print("Running:", comm)
                args = shlex.split(comm)
                check_call(args, stdout=open(os.devnull, 'wb'))

            except RuntimeError:
                print("Error selecting PRIM/SEC events in evtfile ", evt)
                print(comm)
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise

            # FWHM & BIAS for Erecons
            # read Erecons (SIGNAL) column in numpy array (in keV)
            fpre = fits.open(evt, memmap=True)
            nrows = fpre[1].header["NAXIS2"]
            assert nrows > 0, "Empty evt file (%s): nrows=0 " % evt
            ftab = fpre[1].data
            SIGNALmean = ftab['SIGNAL'].mean()
            SIGNALsigma = ftab['SIGNAL'].std()
            fpre.close()
            del fpre[1].data

            fwhmEreal = 0
            fwhmEreal_err = 0
            biasEreal = 0

            print("evt, SIGNALsigma=", evt, SIGNALsigma)

            fwhm = SIGNALsigma * 2.35 * 1000.
            print("evt, fwhm=", evt, fwhm)

            # FWHM for reconstructed events in eV:
            fwhmErecons = '{0:0.5f}'.format(fwhm)
            # EBIAS = <Erecons> - monoEeV
            bias = SIGNALmean * 1000. - monoEeV
            biasErecons = '{0:0.5f}'.format(bias)
            # print("SIGNALmean=", SIGNALmean)
            # print("monoEeV=", monoEeV)
            # print("biasErecons=", biasErecons, "\n")

            # calculate calibrated energies if polyfit coeffs are provided and
            # previous fwhm is not NaN (some NULL Eners)
            # FWHM only calculated for those pulses reconstructed with HRes
            # if any(a != 0 for a in coeffs) and not math.isnan(fwhm):

            if coeffsFile and not math.isnan(fwhm):
                print("=============================================")
                print("CALIBRATING ENERGIES.........................")
                print("=============================================")

                with open(coeffsFile, "rt") as f:
                    fileCont = f.read()   # JSON file
                    if 'surface' in fileCont:
                        ftype = "surface"
                    elif fileCont[0] == '{':
                        ftype = 'json'
                    else:
                        ftype = 'poly'

                if ftype == 'poly':
                    evtcalib = evt.replace(".fits", ".calib1D")
                    eresolFile = eresolFile.replace(".json", ".json1D")
                elif ftype == 'surface':
                    evtcalib = evt.replace(".fits", ".calib2D")
                    eresolFile = eresolFile.replace(".json", ".json2D")
                elif ftype == 'json':
                    evtcalib = evt.replace(".fits", ".calibItr")
                    eresolFile = eresolFile.replace(".json", ".jsonItr")
                else:
                    raise ValueError("Incorrect Coeffs file type")

                # locate coefficients in coeffs table
                # ------------------------------------
                alias = ("pL" + str(pulseLength) + "_" + detMethod + "_" +
                         fdomain + "_" + labelLib + "_" +
                         reconMethod + str(filterLength) + pBStr + IfitStr + LbTStr +
                         smprtStr + jitterStr + noiseStr + bbfbStr +
                         LcStr + s0Str + lagsStr + ctStr)
                # alias = ("pL" + str(pulseLength) + "_" + detMethod + "_" +
                #         fdomain + "_" + labelLib + "_" +
                #         reconMethod + str(filterLength) + pBStr + IfitStr +
                #         smprtStr + jitterStr + noiseStr + bbfbStr +
                #         LcStr + s0Str + "_base100")
                if "NM" in reconMethod:
                    alias = ("pL" + str(pulseLength) + "_" + detMethod + "_" +
                             fdomain + "_" + labelLib + "_" +
                             reconMethod.split("NM")[0] + str(filterLength) +
                             "NM" + reconMethod.split("NM")[1] + LbTStr +
                             smprtStr + jitterStr + noiseStr + bbfbStr + LcStr)
                # unchecked option...
                # if 'fixedlib' in labelLib and 'OF' not in labelLib:
                #    alias = (detMethod + "_" + labelLib + "OF_" + "_" +
                #             reconMethod + str(pulseLength) + smprtStr +
                #             jitterStr + noiseStr + bbfbStr)
                print("Using alias=", alias)
                convertEnergies(evt, evtcalib, coeffsFile, alias)

                fcal = fits.open(evtcalib, memmap=True)
                nrows = fcal[1].header["NAXIS2"]
                assert nrows > 0, "Empty file (%s): nrows=0 " % evtcalib

                fcaltab = fcal[1].data
                ErealKeVmean = fcaltab['SIGNAL'].mean()
                ErealKeVsigma = fcaltab['SIGNAL'].std()
                fcal.close()
                del fcal[1].data

                print("=============================================")
                print("EXTRACTING corrected FWHMs...................")
                print("=============================================")
                print("Output eresolFile: ", eresolFile, "\n")

                fwhm = ErealKeVsigma * 2.35 * 1000.
                fwhm_err = fwhm / math.sqrt(2 * nrows - 2)
                bias = ErealKeVmean * 1000. - monoEeV
                # FWHM for reconstructed PRIMARIES in eV:
                fwhmEreal = '{0:0.5f}'.format(fwhm)
                fwhmEreal_err = '{0:0.5f}'.format(fwhm_err)
                biasEreal = '{0:0.5f}'.format(bias)
                # print("bias=", bias, "\n")
                # print("...biasEreal=", biasEreal, "\n")

            # write FWHM, BIAS values to dictionary for PRIM/SEC/ALL
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
