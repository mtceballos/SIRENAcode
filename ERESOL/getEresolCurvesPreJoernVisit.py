"""
RESOLUTION CURVES for pairs of pulses

python getEresolCurves.py

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
from subprocess import Popen, check_call, check_output, STDOUT
from astropy.io import fits
from astropy.table import Table, Column

# ----GLOBAL VARIABLES -------------
PreBufferSize = 1000
tstartPulse1 = int(PreBufferSize - 1)
nSimPulses = 20000  # minimum number of simulated pulses in simulPairs.csh (approx.; to get filenames)
separation = '20000'  # if unique separation
cwd = os.getcwd()
tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"


def getEresolCurves(pixType, labelLib, monoEkeV, reconMethod, filterMeth, nsamples, fdomain, scaleFactor,
                    samplesUp, nSgms, ACDC, OFInterp, coeffs):
    """
    :param array: Array acronym (SPA, LPA1, LPA2, LPA3)
    :param labelLib: Label identifying the library (monolib, multilib, multilibOF, fixedlib, fixedlib3)
    :param monoEkeV: Monochromatic energy (keV) of input simulated pulses
    :param reconMethod: Energy reconstruction Method (OPTFILT, WEIGHT, WEIGHTN, I2R, I2RBIS)
    :param filterMeth: Optimal Filtering Method (F0 or B0)
    :param nsamples: Pulse length & noise length samples
    :param fdomain: Optimal Filtering domain (F(req) or T(ime))
    :param scaleFactor: Param scaleFactor for gennoise
    :param samplesUp: Param samplesUp for gennoise
    :param nSgms: Param nSgms for gennoise
    :param ACDC: Operating Current (AC or DC)
    :param OFInterp: OF interpolation method: MF (matched filter) or DAB (linear expansion of signal)
    :param coeffs: Polynomial coefficientes (a0 + a1*x + a2*x^2 + a3*x^3 +...) for fit of Erecons vs. Ecalib in keV
    :return: file with energy resolutions for the input pairs of pulses
    """

    # print("En rutina:")
    # print("array=", array)
    # print("labelLib=", labelLib)
    # print("monoEkeV=", monoEkeV)
    # print("reconMethod=", reconMethod)
    # print("filterMeth=", filterMeth)
    # print("nsamples=", nsamples)
    # print("fdomain=", fdomain)
    # print("scaleFactor=", scaleFactor)
    # print("samplesUp=", samplesUp)
    # print("nSgms=", nSgms)
    # print("ACDC=", ACDC)
    # print("OFInterp=", OFInterp)
    # print("coeffs=", coeffs)

    global sepsStr, tstartPulse1, tstartPulse2, nSimPulses
    if tstartPulse1 > 0:
        nSgms = 0
        samplesUp = 0
        scaleFactor = 0

    monoEeV = float(monoEkeV) * 1000.
    tessim = "tessim" + array
    Fil = ""
    OFLib = "no"
    if labelLib == "multilibOF":
        OFLib = "yes"
    if ('multi' in labelLib) and ('WEIGHT' not in reconMethod):
        labelLib += OFInterp
    ACDC = ACDC.lower()

    simDir = cwd + "/PAIRS/tessim" + array
    resultsDir = cwd + "/PAIRS/eresol" + array
    simSIXTEdir = "/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE"

    # -- LIB & NOISE dirs and files ----------
    noiseDir = simSIXTEdir + "/NOISE/" + tessim
    space = "ADC"
    if reconMethod in ("I2R", "I2RBIS"):
        space = reconMethod.lstrip('I2')
        space = space.rstrip('IS')

    noiseFile = noiseDir + "/noise" + str(nsamples) + "samples_" + tessim + "_B0_100s_pairscps_" + space + ".fits"
    libDirRoot = simSIXTEdir + "/LIBRARIES/" + tessim
    libDir = libDirRoot + "/GLOBAL/" + space
    if OFLib == "yes":
        libDir += "LIB"
    # OPTFILT (mono & multi & fixed lib) WEIGHT, WEIGHTN (multilib) I2R (multi, mono & fixed lib)
    if 'multi' in labelLib:
        libFile = libDir + "/libraryMultiE_GLOBAL_PL" + str(nsamples) + "_" + tessim + Fil + ".fits"
    elif "mono" in labelLib:
        libFile = libDir + "/libraryMultiE_PL" + str(nsamples) + "_" + monoEkeV + "keV_" + tessim + Fil + ".fits"
    elif 'fixedlib' in labelLib:
        libFile = libDir + "/libraryMultiE_PL" + str(nsamples) + "_1keV_" + tessim + Fil + ".fits"
    elif 'fixedlib3' in labelLib:
        libFile = libDir + "/libraryMultiE_PL" + str(nsamples) + "_3keV_" + tessim + Fil + ".fits"
    elif 'nonoise' in labelLib:
        libFile = libDir + "/libraryMultiE_PL" + str(nsamples) + "_" + monoEkeV + "keV_" + tessim + "Nonoise.fits"

    # IN case OLD tessim (no trigger must be used)
    if ACDC == 'none':
        nSimPulses = 3000  # minimum number of simulated pulses in simulPairs.csh (approx) -> old_tessim
        noiseFile = noiseDir + "/noise" + str(nsamples) + "samples_" + tessim + "_B0_100s_pairscps.fits"
        tstartPulse1 = int(PreBufferSize + 1)

    # --- get baseline from noise file (in current) ---
    #noisef = fits.open(noiseFile)
    #baseline = noisef[1].header["BASELINE"]
    #noisef.close()

    os.chdir(resultsDir)

    # -- Define pulses separations based on pixel type --------
    if array == "SPA":

        # sepsStr = ['00004', '00006', '00009', '00014', '00021', '00031', '00047', '00071', '00108', '00163', '00246',
        #           '00371', '00560', '00845', '01276', '01926', '02907', '04389']
        sepsStr = [separation]

    elif array in ("LPA1", "LPA2", "LPA3"):

        sepsStr = ['00004', '00006', '00009', '00014', '00021', '00031', '00047', '00071', '00108', '00163', '00246',
                   '00371', '00560', '00845', '01276', '01926', '02907', '04389', '06625', '10000']

    if array == "LPA3":
        scaleFactor = 0.02
        # if triggering is required, do filtering
        if tstartPulse1 == 0: Fil = "Fil"

    # -- Create output ERESOL file -----

    root = ''.join([str(nSimPulses), 'p_SIRENA', str(nsamples), '_', monoEkeV, 'keV_', str(filterMeth),
                    str(fdomain), '_', str(labelLib), '_', str(reconMethod)])
    if tstartPulse1 > 0: root += "_NTRIG"
    eresolFile = "eresol_" + root + ".dat"

    # write header in eresol file
    feresol = open(eresolFile, 'w')
    feresol.write(
        "samSep  FWHMP(Erecons)  FWHMS(Erecons)  EBIASP  EBIASS  FWHMP(Ereal)  FWHMS(Ereal)  FWHMP_ERR   FWHMS_ERR\n")

    # --- Process input data files -------
    for sep12 in sepsStr:
        sep = float(sep12)
        inFile = simDir + "/sep" + sep12 + "sam_" + str(nSimPulses) + "p_" + monoEkeV + "keV.fits"
        evtFile = "events_sep" + sep12 + "sam_" + root + ".fits"
        print("=============================================")
        print("Using file: ", inFile)
        print("Using library: ", libFile)
        print("Using noisefile: ", noiseFile)
        #print("Using Baseline: ", baseline)
        print("Setting evtFile: ", evtFile)
        print("Setting eresolFile: ", eresolFile)
        print("=============================================")

        # when run also in detection mode, run it iteratively to get the best possible combination
        # of sigmas/samples able to detect all pulses

        # tstartPulse1 = PreBufferSize + 1 # old tessim
        # tstartPulse2 = tstartPulse1 + sep

        # tstartPulse1 = int(PreBufferSize + 1)  # 1002 for fv, 1001 for GSL
        tstartPulse1 = int(PreBufferSize - 1 - 1)  # (there is currently a problem for primary pulses)
        # tstartPulse2 = int(tstartPulse1 + sep)  # 1000 + sep for fv; 1000+sep-1 for GSL
        tstartPulse2 = int(tstartPulse1 + 1 + sep)  # (there is currently a problem for primary pulses)
        tstartPulse3 = 0

        if reconMethod == "I2RBIS":  # In I2RBIS, pulse shape is a derivative and pulse template is incomplete otherwise
            tstartPulse1 -= 1
            tstartPulse2 -= 1

        if tstartPulse1 > 0:
            nSgms = 0
            samplesUp = 0
            scaleFactor = 0

        # -- SIRENA processing -----
        if os.path.isfile(evtFile):
            # if events already exist get number nSgms used to get events
            print("Event file", evtFile, " already DOES exist: recalculating FWHMs...")
            try:
                comm = "fkeyprint infile=" + evtFile + "+0 keynam=HISTORY"
                args = shlex.split(comm)
                nSgmsStr = check_output(args)
                #print(nSgmsStr)
                # look for ocurrence of nSgms (before, nSgms, after) and take "after"
                nSgms = nSgmsStr.partition("nSgms") # as in "HISTORY P19 nSgms = 0"
                # as "after" has a lot of lines, separate by newline and take first one
                nSgms = nSgms[2].splitlines()[0] # = 0
                # remove "=" and get numerical value
                nSgms = nSgms.split()[1]
                # print(nSgms)
            except:
                print("Error reading pre-exiting evtFile:", evtFile)
                print(comm)
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise
        else:
            sigmasMin = 0
            sigmasMax = 100
            if float(monoEkeV) == 1: sigmasMax = 50
            ite = 0
            nTrigPulses = nSimPulses # to initialize
            ndetpulses = 0
            print("\nRunning SIRENA for detection & reconstruction")
            while nTrigPulses != ndetpulses and ite < 50:
                ite += 1
                comm = ("tesreconstruction Recordfile=" + inFile + " TesEventFile=" + evtFile + " Rcmethod='SIRENA'" +
                        " PulseLength=" + str(nsamples) + " LibraryFile=" + libFile + " scaleFactor=" +
                        str(scaleFactor) + " samplesUp=" + str(samplesUp) + " nSgms=" + str(nSgms) +
                        " mode=1 NoiseFile=" + noiseFile + " OFLib=" + OFLib + " FilterDomain=" + fdomain +
                        " FilterMethod=" + filterMeth + " clobber=yes intermediate=0 " + " EnergyMethod=" +
                        reconMethod + " PixelType=" + array + " tstartPulse1=" + str(tstartPulse1) +
                        " tstartPulse2=" + str(tstartPulse2) + " tstartPulse3=" + str(tstartPulse3) +
                        " OFInterp=" + OFInterp)
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
                ndetpulses = 0
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

                print("   With nSgms=", nSgms, " => sim/det: ", nTrigPulses, "/", ndetpulses)
                if nTrigPulses != ndetpulses:
                    print("   Repeating detection process...")
                    if nTrigPulses < ndetpulses:
                        # increase sigmas
                        sigmasMin = nSgms
                    else:
                        # decrease sigmas
                        sigmasMax = nSgms
                    nSgms = sigmasMin + (sigmasMax - sigmasMin) / 2.
                    print("   Trying with nSgms=", nSgms, "(sigmasMin/Max:", sigmasMin, "/", sigmasMax)
                    if ite == 50:
                        os.remove(evtFile)
                        print("Reconstruction NOT finished with nSgms=", nSgms)
                    continue
                print('Reconstruction SUCCESSFULLY finished with nSgms={0:0.1f}'.format(nSgms))

                # -----------------------------------------
                # EVENT file processing to calculate FWHM
                # -----------------------------------------
        rootEvt = os.path.splitext(evtFile)[0]
        evt = [rootEvt + "_primaries.fits", rootEvt + "_secondaries.fits"]
        fwhmErecons = ["0.000", "0.000"]
        fwhmErecons_err = ["0.000", "0.000"]
        biasErecons = ["0.000", "0.000"]
        fwhmEreal = ["0.000", "0.000"]
        fwhmEreal_err = ["0.000", "0.000"]
        for ip in (0, 1):  # PRIMARIES & SECONDARIES
            if os.path.isfile(evt[ip]):
                os.remove(evt[ip])

            # use only rows with GRADE1>0 (non initially truncated)
            try:
                # PRIMARIES
                comm = "fselect infile=" + evtFile + " outfile=" + evt[ip] + " expr='#ROW%2!=0 && GRADE1>0' clobber=yes"
                if ip == 1: #SECONDARIES
                    comm = "fselect infile=" + evtFile + " outfile=" + evt[ip] + " expr='#ROW%2==0 && GRADE1>0' clobber=yes"
                args = shlex.split(comm)
                check_call(args, stdout=open(os.devnull, 'wb'))
            except:
                print("Error selecting PRIMARY events in evtfile ", evt[ip])
                print(comm)
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise
            f = fits.open(evt[ip], memmap=True)
            nrows = f[1].header["NAXIS2"]
            assert nrows > 0, "Empty evt file (%s): nrows=0 " % evt[ip]

            # ---- if calibration parameters are provided, calculate also corrected energies
            # ---- if not, calculate only reconstructed energies FWMH and BIAS
            # FWHM & BIAS for Erecons
            # read Erecons (SIGNAL) column in numpy array (in keV)
            ftab = f[1].data
            SIGNALmean = ftab['SIGNAL'].mean()
            SIGNALsigma = ftab['SIGNAL'].std()

            fwhm = 0
            fwhm_err = 0
            bias = 0

            fwhm = SIGNALsigma * 2.35 * 1000.
            fwhm_err = fwhm / math.sqrt(2 * nrows - 2)
            fwhmErecons[ip] = '{0:0.5f}'.format(fwhm)  # FWHM for reconstructed events in eV
            fwhmErecons_err[ip] = '{0:0.5f}'.format(fwhm_err)
            # EBIAS = <Erecons> - monoEeV
            bias = SIGNALmean * 1000. - monoEeV
            biasErecons[ip] = '{0:0.5f}'.format(bias)

            # calculate corrected energies if polyfit coeffs are provided and previous fwhm is not NaN (some NULL Eners)
            if any(a != 0 for a in coeffs) and not math.isnan(fwhm):
                i = 0
                fwhm = 0
                fwhm_err = 0
                bias = 0
                ErealKeV = []
                for SIGNALKeV in ftab['SIGNAL']:
                    # read fitting coeffs taken from polyfit2Bias.R (a0, a1, a2, a3, a4, a5, a6)
                    #  as in y = a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5
                    # where y=E_reconstructed and x=Ecalibration (keV)
                    npCoeffs = np.array(coeffs)
                    npCoeffs[0] -= SIGNALKeV  # subtract "y" (0 = a0 + a1*x + a2*x^2 + ... - y)
                    npCoeffsRev = npCoeffs[::-1]  # reversed to say fit with poly1d definition
                    polyfit = np.poly1d(npCoeffsRev)
                    # get real root (value of Ereal for a given Erecons )
                    r = np.roots(polyfit)
                    # real && >0 roots
                    rreal = r.real[abs(r.imag) < 1e-5]
                    rrealpos = rreal[rreal > 0]
                    # closest root
                    rclosest = min(enumerate(rrealpos), key=lambda x: abs(x[1]-SIGNALKeV))[1]  # (idx,value)
                    ErealKeV = np.append(ErealKeV, rclosest)
                ErealKeVsigma = ErealKeV.std()
                fwhm = ErealKeVsigma * 2.35 * 1000.
                fwhm_err = fwhm / math.sqrt(2 * nrows - 2)
                # print("FWHM=",fwhm)
                # print("FWHM_err=",fwhm_err)
                fwhmEreal[ip] = '{0:0.5f}'.format(fwhm)  # FWHM for reconstructed PRIMARIES in eV
                fwhmEreal_err[ip] = '{0:0.5f}'.format(fwhm_err)
            f.close()
            del f[1].data

    # write final data to eresol file
    feresol.write(' {0}     {1}            {2}       {3}  {4}     {5}         {6}        {7}       {8} \n'.format(sep12,
            fwhmErecons[0], fwhmErecons[1], biasErecons[0],  biasErecons[1], fwhmEreal[0],
            fwhmEreal[1], fwhmEreal_err[0], fwhmEreal_err[1]))


#
# --------- Read input parameters and PROCESS simulations -----------------
#
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Get Energy resolution for pairs of pulses', prog='getEresolCurves')

    parser.add_argument('--pixType', help='Extension name in the FITS pixel definition file (SPA*, LPA1*, LPA2*, LPA3*)')
    parser.add_argument('--lib', choices=['monolib', 'multilib', 'multilibOF', 'fixedlib', 'fixedlib3'],
                        default='multilib',
                        help='Label identifying the library (monolib, multilib, multilibOF, fixedlib, fixedlib3)')
    parser.add_argument('--monoEnergy', help='Monochromatic energy (keV) of input simulated pulses')
    parser.add_argument('--reconMethod', choices=['OPTFILT', 'WEIGHT', 'WEIGHTN', 'I2R', 'I2RBALL', 'I2RBNOL'],
                        default='OPTFILT',
                        help='Energy reconstruction Method (OPTFILT, WEIGHT, WEIGHTN, I2R, I2RBALL, I2RBNOL)')
    parser.add_argument('--filter', choices=['F0', 'B0'], default='F0',
                        help='Optimal Filtering Method (F0, B0) [default %(default)s]')
    parser.add_argument('--nsamples', type=int, help='pulse length & noise length samples')
    parser.add_argument('--fdomain', choices=['F', 'T'], default='T',
                        help='Optimal Filtering domain (F(req) or T(ime)) [default %(default)s]')
    parser.add_argument('--scaleFactor', type=float, default=0.005,
                        help='Param scaleFactor for gennoise [default %(default)s]')
    parser.add_argument('--samplesUp', type=float, default=2,
                        help='Param samplesUp for gennoise [default %(default)s]')
    parser.add_argument('--nSgms', type=float, help='Param nSgms for gennoise [default %(default)s]', default=20)
    parser.add_argument('--acbias', choices=['yes', 'no'],
                        help='Operating Current ("yes" for AC or "no" for DC)')
    parser.add_argument('--interp', default='DAB', choices=['MF', 'DAB'],
                        help=('OF interpolation method: MF (matched filter) or DAB (signal linear expansion)'
                              '[default %(default)s]'))
    parser.add_argument('--coeffs', type=float, nargs=6, default=(0, 0, 0, 0, 0, 0),
                        help='Polynomial coefficientes (a0+a1*x+a2*x2+a3*x3+a4*x4+a5*x5) for fit of Erecons vs. '
                             'Ecalib in keV [default %(default)s]')

    args = parser.parse_args()

    # print("array=",args.array)

    getEresolCurves(array=args.array, labelLib=args.lib, monoEkeV=args.monoEnergy, reconMethod=args.reconMethod,
                    filterMeth=args.filter, nsamples=args.nsamples, fdomain=args.fdomain, scaleFactor=args.scaleFactor,
                    samplesUp=args.samplesUp, nSgms=args.nSgms, ACDC=args.ACDC, OFInterp=args.interp,
                    coeffs=args.coeffs)
