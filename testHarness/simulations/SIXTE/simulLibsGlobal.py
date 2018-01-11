"""
CREATE GLOBAL LIBRARY for simulated pairs of pulses for calibration

python simulLibsGlobal.py

"""

# ----IMPORT MODULES --------------
from __future__ import print_function
import os
import sys
import shlex
import shutil
import tempfile
from subprocess import check_call, STDOUT
from astropy.io import fits
import xml.etree.ElementTree as ET


# ----GLOBAL VARIABLES -------------
# -- CALIB dirs & files --
simSIXTEdir = "/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE"
XMLdir = os.environ["SIXTE"] + "/" + "share/sixte/instruments/athena/1469mm_xifu"
# XMLfile = XMLdir + "/" + "xifu_baseline.xml"
XMLfile = XMLdir + "/" + "xifu_detector_hex_baseline.xml"
PIXIMPACTdir = simSIXTEdir + "/LIBRARIES/PIXIMPACT"
XMLtree = ET.parse(XMLfile)
XMLroot = XMLtree.getroot()
for samplefreq in XMLroot.findall('samplefreq'):
    samprate = samplefreq.get('value')
PreBufferSize = 1000
# separation = 20000  # LPA1shunt
separation = 40000    # LPA2shunt
dtaums = separation / float(samprate) * 1000.  # separation time (ms) between pulses
tstart = 0.5/float(samprate)  # added to solve floating point inaccuracies due to sampling rate (Christian's mail 31/03/2017)

# With triggerTH=20, 0.2 keV pulses trigger 1 sample late (1001 instead of 1000), but ALL of them
#                    0.5 keV : some pulses trigger 1 sample late and some pulses trigger ok (1000) => set to 50
#                    >= 1 keV:  ALL trigger OK

triggerTH = {'LPA1shunt': 50, 'LPA2shunt': 20}
pulsesPerRecord = {'LPA1shunt': 2, 'LPA2shunt': 1}
noise = ""
Fil = ""
pixel = 1

cwd = os.getcwd()
tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"


def rmLastAndFirst(simfile, ppr):

    """
    rm first (and LAST) record and update NETTOT (first starts high and last can be cut)
    (see Christian's email from 19 Jan 2017 @ EURECA)
    Also update number of pulses in NETTOT

    :type simfile: str
    :param simfile: simulated file where cleaning must be done
    :type ppr: int
    :param ppr: pulses per record in simulations

    """

    fsim = fits.open(simfile)
    nrows = fsim[1].header["NAXIS2"]
    fsim.close()
    assert nrows > 1, "Tessim failed for (%s): just one row present " % simfile

    try:
        comm = "fdelrow infile=" + simfile + "+1 firstrow=" + str(nrows) + " nrows=1 confirm=no proceed=yes"
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = "fdelrow infile=" + simfile + "+1 firstrow=1 nrows=1 confirm=no proceed=yes"
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        fsim = fits.open(simfile, mode='update')
        nrows2 = fsim[1].header['NAXIS2']
        assert nrows2 == nrows-2, "Failure removing initial & last rows in (%s): " % simfile
        nettot = nrows2 * ppr  # new number of pulses (=nofrecords in LPA2; ==2*nofrecords in LPA1)

        fsim[1].header['NETTOT'] = nettot
        fsim.close()

    except:
        print("Error running FTOOLS to remove initial & last rows in ", simfile)
        os.chdir(cwd)
        shutil.rmtree(tmpDir)
        raise


def simulGlobalLibs(pixName, space, pulseLength, libEnergies, largeFilter, nsamples, nSimPulses, acbias,
                    tstartPulse1All, tstartPulse2All, tstartPulse3All, createLib, noiseMat, weightMat):
    """
    :type pixName: str
    :param pixName: Extension name in FITS file pixel definition (SPA*, LPA1*, LPA2*, LPA3*)
    :type space: str
    :param space: Input Data Space (ADC, R, RALL, RNOL, RFITTED)
    :type pulseLength: int
    :param pulseLength: Pulse length
    :type nsamples: int
    :param nsamples: noise length samples
    :type libEnergies: str
    :param libEnergies: list of calibration energies (keV)
    :type largeFilter: int
    :param largeFilter: size of the extra-large filter to avoid record -length effects (samples)
    :type nSimPulses: int
    :param nSimPulses: number of pulses in simulated files (uses to create/choose files and to name library)
    :type acbias: str
    :param acbias: Operating Current: AC (acbias=yes) or DC (acbias=no)
    :type tstartPulse1All: int
    :param tstartPulse1All: list of initial sample for 1st pulses in each record. If empty, just simulate calib. files
                **  0 if detection is to be performed
    :type tstartPulse2All: int
    :param tstartPulse2All: list of initial sample for 1st pulses in each record. Required if PAIRS are simulated
                **  0 if detection is to be performed
                ** -1 if sample is to be calculated (tstartPulse1+ separation)
    :type tstartPulse3All: int
    :param tstartPulse3All: list of initial sample for 1st pulses in each record. Required if TRIOS are simulated
                **  0 if detection is to be performed
                ** -1 if sample is to be calculated (tstartPulse2+ separation)
    :type createLib: int
    :param createLib: if library should be created (1) or script should just create simulated files (0)
    :param noiseMat: should the Noise matrices HDU be created? (yes/no)
    :param weightMat: should the Weight matrices HDU be created? (yes/no)
    :return: simulated calibration pairs of pulses && Global library from them

    """

    global PreBufferSize, separation, Fil, cwd, simSIXTEdir, samprate, triggerTH, pulsesPerRecord
    tessim = "tessim" + pixName

    # Calibration energies and Tstarts of pulses
    tstartPulse1 = dict(zip(libEnergies, tstartPulse1All))
    tstartPulse2 = dict(zip(libEnergies, tstartPulse2All))
    tstartPulse3 = dict(zip(libEnergies, tstartPulse3All))

    # calculate simTime so that nSimPulses are simulated (1 pulses per record)
    simTime = nSimPulses * (PreBufferSize + separation) / float(samprate) # if tesgenimpacts
    simTime = '{0:0.0f}'.format(simTime)

    # Sigmas and Samples and scaleFactor for Detection
    samplesUpAll = [0 for i in range(0, len(libEnergies))]
    nSgmsAll = [0 for i in range(0, len(libEnergies))]
    if space == "ADC" or space == "R" or space == "RNOL" or space == "RFITTED":
        samplesUpAll = [3 for i in range(0, len(libEnergies))]
        nSgmsAll = [4 for i in range(0, len(libEnergies))]
    elif space == "RALL":
        samplesUpAll = [2 for i in range(0, len(libEnergies))]
        nSgmsAll = [2.5 for i in range(0, len(libEnergies))]
    scaleFactor = 0.  # no Filtering
    nSgms = dict(zip(libEnergies, nSgmsAll))
    samplesUp = dict(zip(libEnergies, samplesUpAll))
    maxFilterLength = pulseLength
    lF = ""
    if largeFilter > 0:
        lF = " largeFilter= " + str(largeFilter)
        maxFilterLength = largeFilter

    # Trigger sizes in tessim
    triggerSizeTS = PreBufferSize + maxFilterLength + 1000       # ___|\_______________ooo
    triggerTS3val = triggerSizeTS - PreBufferSize

    # -- CALIB dirs & files --
    SIMFILESdir = simSIXTEdir + "/LIBRARIES/" + tessim

    # -- LIB & NOISE dirs and files ----------
    noiseDir = simSIXTEdir + "/NOISE/" + tessim
    noiseFile = noiseDir + "/noise" + str(nsamples) + "samples_" + tessim + "_B0_" + space + ".fits"
    PixTypeFile = "'file:" + simSIXTEdir + "/newpixels.fits[" + pixName + "]'"

    #
    # Create simulated FILES for libraries and LIBRARY itself
    #
    libDirRoot = simSIXTEdir + "/LIBRARIES/" + tessim
    libDir = libDirRoot + "NEW/GLOBAL/" + space
    libFile = libDir + "/libraryMultiE_GLOBAL_PL" + str(pulseLength) + "_" + str(nSimPulses) + "p" + Fil + "_SHORT.fits"

    # libFile += "5-9"  # for library creation in steps

    # if os.path.isfile(libFile):
    #    os.remove(libFile)
    evttmpFile = tempfile.NamedTemporaryFile()
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
        if monoEkeV == "0.5":
            triggerTH["LPA2shunt"] = 50

        print("=============================================")
        print("Adding monochromatic energy", monoEkeV, "keV")
        print("=============================================")
        print("simTIme=", simTime, "\n")
        print("Temporary event file in: ", evttmpFile.name)
        # simulate SIXTE file with tesconstpileup + tessim
        root0 = "mono" + monoEkeV + "_sep" + str(separation) + "_pix" + str(pixel) + "_" + str(nSimPulses) + "p"
        root = root0 + "_" + str(pulseLength) + noise
        # pixFile = PIXIMPACTdir + "/" + root0 + "_trSz" + str(triggerSizeTC) + ".piximpact"
        pixFile = PIXIMPACTdir + "/" + root0 + ".piximpact"
        simFile = SIMFILESdir + "/" + root + ".fits"
	print(simFile)

        # -- TESGENIMPACTS: generate impacts for well separated single pulses --
        # -- TESSIM: simulate well separated pulses --
        if not os.path.isfile(simFile):
            print("Simulating & triggering pulses with TESSIM to ", simFile)
            comm = ("tesgenimpacts PixImpList=" + pixFile + " mode=const tstart=" + str(tstart) + " tstop=" + simTime +
                    " EConst=" + monoEkeV + " dtau=" + str(dtaums) + " clobber=yes")
            print("Generating PIXIMPACT ", pixFile, " running:\n", comm)
            try:
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except:
                print("Error running task for piximpact list generation")
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise

            comm = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFile + " Streamfile=" + simFile +
                    " tstart=0 tstop=" + str(simTime) + " triggerSize=" + str(triggerSizeTS) +
                    " preBuffer=" + str(PreBufferSize) + " acbias=" + acbias + " triggertype='diff:3:" +
                    str(triggerTH[pixName]) + ":" + str(triggerTS3val) + "'" + " sample_rate=" + samprate +
                    " PixType=" + PixTypeFile)

            print("Running tessim for simulation\n", comm)
            try:
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except:
                print("Error running TESSIM for simulation\n")
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise
            rmLastAndFirst(simFile, pulsesPerRecord[pixName])

        if not createLib:
            continue  # (to just simulate pulses files)

        # -- Create/update LIBRARY --
        comm = ("tesreconstruction Recordfile=" + simFile + " TesEventFile=" + evttmpFile.name + " Rcmethod=SIRENA" +
                " PulseLength=" + str(pulseLength) + " LibraryFile=" + libFile + " scaleFactor=" + str(scaleFactor) +
                " samplesUp=" + str(samplesUp[monoEkeV]) + " nSgms=" + str(nSgms[monoEkeV]) +
                " mode=0 clobber=yes intermediate=0" + " monoenergy=" + str(monoEeV) + " EventListSize=1000" +
                " tstartPulse1=" + str(tstartPulse1[monoEkeV]) + " tstartPulse2=" + str(tstartPulse2[monoEkeV]) +
                " tstartPulse3=" + str(tstartPulse3[monoEkeV]) + lF + " NoiseFile=" + noiseFile +
                " XMLFile=" + XMLfile + energyMethod + " hduPRECALWN=" + weightMat + " hduPRCLOFWM=" + noiseMat)
        args = shlex.split(comm)
        print("SIRENA reconstruction to add a line to the library, running command:\n", comm)
        try:
            check_call(args, stderr=STDOUT)
        except:
            print("Error running SIRENA to add new line to library")
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
            raise
        print("SIRENA reconstruction to add a line to the library: done\n", comm)
        # -- check that number of reconstructed pulses is the same than number of simulated pulses
        fsim = fits.open(simFile)
        nettot = fsim[1].header['NETTOT']  # number of simulated/triggered pulses
        fsim.close()
        fevt = fits.open(evttmpFile.name)
        ndetpulses = fevt[1].header['NAXIS2']  # number of detected pulses
        fevt.close()
        # assert nettot == ndetpulses, \
        #    "Detection failure: all pulses (%d) should be detected (%d) in %s: " % (nettot, ndetpulses, simFile)
        print("Simpulses=", nettot, "and detected pulses=", ndetpulses, "in ", simFile)

    os.chdir(cwd)
    shutil.rmtree(tmpDir)

#
# --------- Read input parameters and PROCESS simulations -----------------
#
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Create GLOBAL library from pairs of pulses', prog='simulLibsGlobal')

    parser.add_argument('--pixName', help='Extension name in FITS pixel definition file (SPA*, LPA1*, LPA2*, LPA3*)',
                        required=True)
    parser.add_argument('--space', choices=['ADC', 'R', 'RALL', 'RNOL', 'RFITTED'], required=True,
                        help='Input Data Space  (ADC, R, RALL, RNOL, RFITTED)')
    parser.add_argument('--pulseLength', type=int, help='pulse length samples', required=True)
    parser.add_argument('--nsamples', type=int, help='noise samples', required=True)
    parser.add_argument('--nSimPulses', type=int, help='Number of Pulses in simulated files', required=True)
    parser.add_argument('--calibEnergies', nargs='*', help="list of energies (keV) from calibration files",
                        required=True)
    parser.add_argument('--tstartPulse1All', nargs='*', type=int,
                        help="list of Tstarts for pulse 1 (if zeros, perform detection)")
    parser.add_argument('--tstartPulse2All', nargs='*', type=int,
                        help="list of Tstarts for pulse 2 (if zeros, perform detection)")
    parser.add_argument('--tstartPulse3All', nargs='*', type=int,
                        help="list of Tstarts for pulse 3 (if zeros, perform detection)")
    parser.add_argument('--largeFilter', type=int, help='Size of extra-large filter', default=0)
    parser.add_argument('--acbias', default='yes', choices=['yes', 'no'],
                        help='Operating Current (AC )(acbias=yes) or DC (acbias=no)) [default %(default)s]')
    parser.add_argument('--noiseMat', default='no', choices=['yes', 'no'],
                        help='Should the Noise Matrices HDU be created? [default %(default)s]')
    parser.add_argument('--weightMat', default='no', choices=['yes', 'no'],
                        help='Should the Weight Matrices HDU be created? [default %(default)s]')

    inargs = parser.parse_args()
    len1 = 0
    len2 = 0
    len3 = 0
    if inargs.tstartPulse1All is not None:
        len1 = len(inargs.tstartPulse1All)
    if inargs.tstartPulse2All is not None:
        len2 = len(inargs.tstartPulse2All)
    if inargs.tstartPulse3All is not None:
        len3 = len(inargs.tstartPulse3All)
    lenE = len(inargs.calibEnergies)
    if len1 == 0:  # just do file simulation
        cLib = 0
        inargs.tstartPulse1All = [0 for i in range(0, lenE)]
        inargs.tstartPulse2All = [0 for i in range(0, lenE)]
        inargs.tstartPulse3All = [0 for i in range(0, lenE)]
    else:
        cLib = 1
        if len2 == 0:
            inargs.tstartPulse2All = [0 for i in range(0, lenE)]
            inargs.tstartPulse3All = [0 for i in range(0, lenE)]
        if len3 == 0:
            inargs.tstartPulse3All = [0 for i in range(0, lenE)]

    if (len2 > 0 and (len2 != len1)) or (len3 > 0 and (len3 != len1 or len3 != len2)):
        print("Error: Inconsistent number of starting samples for 1st/2nd/3rd pulses")
        sys.exit()

    simulGlobalLibs(pixName=inargs.pixName, space=inargs.space, pulseLength=inargs.pulseLength,
                    largeFilter=inargs.largeFilter, libEnergies=inargs.calibEnergies,
                    tstartPulse1All=inargs.tstartPulse1All, tstartPulse2All=inargs.tstartPulse2All,
                    tstartPulse3All=inargs.tstartPulse3All, nsamples=inargs.nsamples,
                    nSimPulses=inargs.nSimPulses, acbias=inargs.acbias, createLib=cLib,
                    noiseMat=inargs.noiseMat, weightMat=inargs.weightMat)
