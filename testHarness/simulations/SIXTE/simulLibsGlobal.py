"""
CREATE GLOBAL LIBRARY for simulated pairs of pulses for calibration

python simulLibsGlobal.py

"""

# ----IMPORT MODULES --------------
from __future__ import print_function
import os
import shlex
import shutil
import tempfile
import auxpy
from time import gmtime, strftime
from subprocess import check_call, STDOUT
from astropy.io import fits
import xml.etree.ElementTree as ET


# ----GLOBAL VARIABLES -------------
PreBufferSize = 1000
# separation = 20000  # LPA1shunt
separation = 40000    # LPA2shunt

# With triggerTH=20, 0.2 keV pulses trigger 1 sample late (1001 instead
#                               of 1000), but ALL of them
#                    0.5 keV : some pulses trigger 1 sample late and some
#                               pulses trigger ok (1000) => set to 50
#                    >= 1 keV:  ALL trigger OK

triggerTH = {'LPA1shunt': 50, 'LPA2shunt': 20}
pixel = 1

cwd = os.getcwd()
tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"


def simulLibsGlobal(pixName, space, samprate, jitter, noise, stoch,
                    pulseLength, libEnergies, largeFilter, nsamples,
                    nSimPulses, acbias, tstartPulse1All,
                    createLib, noiseMat, weightMat):
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
    :type stoch: str
    :param stoch: ("stoch") stochastic or ("") non-stochastic simulations
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
    :type weigthMat: str
    :param weightMat: should the Weight matrices HDU be created? (yes/no)
    :return: simulated calibration pulses pulses && Global library from them

    """

    global PreBufferSize, separation, cwd, simSIXTEdir, triggerTH
    tessim = "tessim" + pixName

    # samprate
    smprtStr = ""
    if samprate == 'samprate2':
        smprtStr = "_samprate2"

    # jitter
    jitterStr = ""
    offset = ""
    if jitter == "jitter":
        jitterStr = "_jitter"
        offset = " offset=-1"

    # noise
    noiseStr = ""
    simnoise = " simnoise=y"
    if noise == "nonoise":
        noiseStr = "_nonoise"
        simnoise = " simnoise=n"

    # stochastic
    stochStr = ""
    stochastic = " stochastic_integrator=n"
    if stoch == "stoch":
        stochStr = "_stoch"
        stochastic = (" stochastic_integrator=y dobbfb=y" +
                      " carrier_frequency = 2e6 bbfb_delay = 40" +
                      " decimation_filter = y")

    # -- CALIB/SIM/NOISE/LIB/XML dirs & files --
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA"
    simSIXTEdir = EURECAdir + "/testHarness/simulations/SIXTE"
    XMLdir = EURECAdir + "/ERESOL"
    PIXIMPACTdir = simSIXTEdir + "/LIBRARIES/PIXIMPACT"
    SIMFILESdir = simSIXTEdir + "/LIBRARIES/" + tessim
    NOISEdir = simSIXTEdir + "/NOISE/" + tessim
    libDir = SIMFILESdir + "/GLOBAL/" + space

    # XMLfile = XMLdir + "/" + "xifu_baseline.xml"
    XMLfile = (XMLdir + "/" + "xifu_detector_hex_baselineNEWgrades" +
               smprtStr + ".xml")
    XMLtree = ET.parse(XMLfile)
    XMLroot = XMLtree.getroot()

    noiseFile = (NOISEdir + "/noise" + str(nsamples) + "samples_" + tessim +
                 "_B0_" + space + smprtStr + jitterStr + stochStr + ".fits")
    # PixTypeFile = "'file:"+ simSIXTEdir + "/newpixels.fits[" + pixName + "]'"
    PixTypeFile = ("'file:" + simSIXTEdir +
                   "/newpix_full.fits[" + pixName + "]'")
    libFile = (libDir + "/libraryMultiE_GLOBAL_PL" + str(pulseLength) +
               "_" + str(nSimPulses) + "p" + smprtStr + jitterStr +
               noiseStr + stochStr + ".fits")
    evttmpFile = tempfile.NamedTemporaryFile()

    for samplefreq in XMLroot.findall('samplefreq'):
        samplingrate = samplefreq.get('value')

    # added to solve floating point inaccuracies due to sampling rate
    # (Christian's mail 31/03/2017):
    tstart = 0.5/float(samplingrate)

    # Calibration energies and Tstarts of pulses
    tstartPulse1 = dict(zip(libEnergies, tstartPulse1All))

    # calculate simTime so that nSimPulses are simulated (1 pulses per record)
    simTime = nSimPulses * (PreBufferSize + separation) / float(samplingrate)
    tstop = tstart + simTime
    simTime = '{0:0.0f}'.format(simTime)
    tstop = '{0:0.0f}'.format(tstop)

    # Sigmas and Samples and scaleFactor for Detection
    samplesUp = 3
    samplesDown = 3
    nSgms = 5
    maxFilterLength = pulseLength
    lF = ""
    if largeFilter > 0:
        lF = " largeFilter= " + str(largeFilter)
        maxFilterLength = largeFilter

    # Trigger sizes in tesconstpileup
    triggerSizeTC = PreBufferSize + separation + separation + 1000
    # and tessim ___|\_______________ooo
    triggerSizeTS = PreBufferSize + maxFilterLength + 1000
    triggerTS3val = triggerSizeTS - PreBufferSize

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
        # print("Temporary event file in: ", evttmpFile.name)
        # simulate SIXTE file with tesconstpileup + tessim
        root0 = ("mono" + monoEkeV + "_sep" + str(separation) + "_pix" +
                 str(pixel) + "_" + str(nSimPulses) + "p")
        root = root0 + "_" + str(pulseLength)
        # pixFile = PIXIMPACTdir +
        #          "/" + root0 + "_trSz" + str(triggerSizeTC) + ".piximpact"
        pixFile = (PIXIMPACTdir + "/" + root0 + smprtStr + jitterStr +
                   ".piximpact")
        simFile = (SIMFILESdir + "/" + root + smprtStr + jitterStr + noiseStr +
                   stochStr + ".fits")

        # -- TESCONSTPILEUP: generate impacts for well separated single pulses
        # -- TESSIM: simulate well separated pulses --
        if not os.path.isfile(simFile):
            print("Simulating & triggering pulses with TESSIM to ", simFile)
            comm = ("tesconstpileup PixImpList=" + pixFile +
                    " XMLFile=" + XMLfile + " timezero=" + str(tstart) +
                    " tstop=" + str(tstop) + " energy=" + monoEkeV +
                    " pulseDistance=" + str(separation) + offset +
                    " TriggerSize=" + str(triggerSizeTC) + " clobber=yes")
            print("Generating PIXIMPACT ", pixFile, " running:\n", comm)
            try:
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except:
                print("Error running task for piximpact list generation")
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise

            commTessim = ("tessim PixID=" + str(pixel) +
                          " PixImpList=" + pixFile +
                          " Streamfile=" + simFile +
                          " tstart=0 tstop=" + str(simTime) +
                          " triggerSize=" + str(triggerSizeTS) +
                          " preBuffer=" + str(PreBufferSize) +
                          " acbias=" + acbias +
                          " triggertype='diff:3:" + str(triggerTH[pixName]) +
                          ":" + str(triggerTS3val) + "'" +
                          " sample_rate=" + samplingrate +
                          simnoise + stochastic + " PixType=" + PixTypeFile)

            print("Running tessim for simulation\n", commTessim)
            try:
                args = shlex.split(commTessim)
                check_call(args, stderr=STDOUT)
            except:
                print("Error running TESSIM for simulation\n")
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise

            auxpy.rmLastAndFirst(simFile, 1)

            # update HISTORY in header[0]
            dateTime = strftime("%Y-%m-%d %H:%M:%S", gmtime())
            history = ("Created & Updated by simulLibsGlobal.py on " +
                       dateTime + "with command: " + commTessim)
            auxpy.updateHISTORY(simFile, history)

        # print("Antes de evaluar createLib\n")
        if not createLib:
            continue  # (to just simulate pulses files)

        # print("Va a crear la libreria\n")
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
                " mode=0 clobber=yes intermediate=0" +
                " monoenergy=" + str(monoEeV) + " EventListSize=1000" +
                " XMLFile=" + XMLfile +
                " hduPRECALWN=" + weightMat + " hduPRCLOFWM=" + noiseMat)

        args = shlex.split(comm)
        print("SIRENA reconstruction to add a line to the library, "
              "running command:\n", comm)
        try:
            check_call(args, stderr=STDOUT)
        except:
            print("Error running SIRENA to add new line to library")
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
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

    os.chdir(cwd)
    shutil.rmtree(tmpDir)


#
# --------- Read input parameters and PROCESS simulations -----------------
#
if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
            description='Create GLOBAL library from pairs of pulses',
            prog='simulLibsGlobal')

    parser.add_argument('--pixName', help=('Extension name in FITS pixel \
                        definition file (SPA*, LPA1*, LPA2*, LPA3*)'),
                        required=True)
    parser.add_argument('--space', required=True,
                        choices=['ADC', 'R', 'RALL', 'RNOL', 'RFITTED'],
                        help='Input Data Space  (ADC, R, RALL, RNOL, RFITTED)')
    parser.add_argument('--samprate', default="", choices=['', 'samprate2'],
                        help="baseline, half_baseline")
    parser.add_argument('--jitter', default="", choices=['', 'jitter'],
                        help="no jitter, jitter")
    parser.add_argument('--noise', default="", choices=['', 'nonoise'],
                        help="noise, nonoise")
    parser.add_argument('--stoch', default="", choices=['', 'stoch'],
                        help="nonstochastic, stochastic")
    parser.add_argument('--pulseLength', type=int, required=True,
                        help='pulse length samples')
    parser.add_argument('--nsamples', type=int, required=True,
                        help='noise samples')
    parser.add_argument('--nSimPulses', type=int, required=True,
                        help='Number of Pulses in simulated files')
    parser.add_argument('--calibEnergies', nargs='*', required=True,
                        help="list of energies (keV) from calibration files")
    parser.add_argument('--tstartPulse1All', nargs='*', type=int,
                        help="list of Tstarts for pulse 1 \
                        (if zeros, perform detection)")
    parser.add_argument('--largeFilter', type=int, default=0,
                        help='Size of extra-large filter')
    parser.add_argument('--acbias', default='yes', choices=['yes', 'no'],
                        help='Operating Current (AC )(acbias=yes) or \
                        DC (acbias=no)) [default %(default)s]')
    parser.add_argument('--createLib', type=int, required=True,
                        help='Create library (1) or only the mono files (0)')
    parser.add_argument('--noiseMat', default='no', choices=['yes', 'no'],
                        help='Should the Noise Matrices HDU be created? \
                        [default %(default)s]')
    parser.add_argument('--weightMat', default='no', choices=['yes', 'no'],
                        help='Should the Weight Matrices HDU be created? \
                        [default %(default)s]')

    inargs = parser.parse_args()
    len1 = 0
    if inargs.tstartPulse1All is not None:
        len1 = len(inargs.tstartPulse1All)
    lenE = len(inargs.calibEnergies)
    cLib = 1
    if len1 == 0:
        inargs.tstartPulse1All = [0 for i in range(0, lenE)]

    simulLibsGlobal(pixName=inargs.pixName, space=inargs.space,
                    samprate=inargs.samprate, jitter=inargs.jitter,
                    noise=inargs.noise, stoch=inargs.stoch,
                    pulseLength=inargs.pulseLength,
                    largeFilter=inargs.largeFilter,
                    libEnergies=inargs.calibEnergies,
                    tstartPulse1All=inargs.tstartPulse1All,
                    nsamples=inargs.nsamples,
                    nSimPulses=inargs.nSimPulses,
                    acbias=inargs.acbias, createLib=inargs.createLib,
                    noiseMat=inargs.noiseMat, weightMat=inargs.weightMat)
