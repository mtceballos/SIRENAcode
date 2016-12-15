"""
CREATE GLOBAL LIBRARY for simulated pairs of pulses for calibration

python simulLibsSingle.py

"""

# ----IMPORT MODULES --------------
from __future__ import print_function
import os
# import sys
import shlex
import shutil
import tempfile
import numpy as np
from subprocess import check_call, STDOUT
from astropy.io import fits
import xml.etree.ElementTree as ET


# ----GLOBAL VARIABLES -------------
# -- CALIB dirs & files --
simSIXTEdir = "/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE"
XMLdir = os.environ["SIXTE"] + "/" + "share/sixte/instruments/athena/1469mm_xifu"
XMLfile = XMLdir + "/" + "xifu_detector_hex_baseline.xml"
PIXIMPACTdir = simSIXTEdir + "/LIBRARIES/PIXIMPACT"
XMLtree = ET.parse(XMLfile)
XMLroot = XMLtree.getroot()
for samplefreq in XMLroot.findall('samplefreq'):
    samprate = samplefreq.get('value')
PreBufferSize = 1000
separation = 100000  # separation between simulated pulses
maxPulseLength = 20000  # maximum length for reconstruction
libEnergiesAll = [0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9]

# perform detection only for 0.2 and 1 keV (bad starting samples in tessim)
tstartPulse1All = np.array([0, PreBufferSize, 0] + [PreBufferSize for i in range(3, len(libEnergiesAll))])
tstartPulse2All = np.array([0, PreBufferSize+separation, 0] +
                           [PreBufferSize+separation for i in range(3, len(libEnergiesAll))])
tstartPulse1 = dict(zip(libEnergiesAll, tstartPulse1All))
tstartPulse2 = dict(zip(libEnergiesAll, tstartPulse2All))
tstartPulse3 = 0

# use this to add individual entries to the library (comment out removal of lib!!!)
libEnergies = libEnergiesAll
# libEnergies = [5, 6, 7, 8, 9]

# length of simulation in "record" in tesconstpileup
triggerSizeTC = PreBufferSize + separation + separation + 1000

noise = ""
Fil = ""
pixel = 1

cwd = os.getcwd()
tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"


def simulGlobalLibs(pixName, space, pulseLength, nsamples, nSimPulses, acbias):
    """
    :type pixName: str
    :param pixName: Extension name in FITS file pixel definition (SPA*, LPA1*, LPA2*, LPA3*)
    :type space: str
    :param space: Input Data Space (ADC, R, RALL, RNOL, RFITTED)
    :type pulseLength: int
    :param pulseLength: Pulse length
    :type nsamples: int
    :param nsamples: noise length samples
    :type nSimPulses: int
    :param nSimPulses: number of pulses in simulated files (uses to create/choose files)
    :type acbias: str
    :param acbias: Operating Current: AC (acbias=yes) or DC (acbias=no)
    :return: simulated calibration pairs of pulses && Global library from them
    """

    global PreBufferSize, separation, Fil, cwd, simSIXTEdir, samprate, tstartPulse2, tstartPulse1, \
        tstartPulse3, triggerSizeTC, libEnergiesAll, maxPulseLength
    tessim = "tessim" + pixName

    # calculate simTime so that nSimPulses are simulated (2 pulses per simulated record)
    simTime = nSimPulses/2. * triggerSizeTC / float(samprate)
    simTime = '{0:0.0f}'.format(simTime)

    # do not use E=0.2 for RALL because of detection problems
    if space == "ADC" or space == "R" or space == "RNOL" or space == "RFITTED":
        samplesUpAll = np.array([3 for i in range(0, len(libEnergiesAll))])
        nSgmsAll = np.array([4 for i in range(0, len(libEnergiesAll))])
    elif space == "RALL":
        samplesUpAll = np.array([2 for i in range(0, len(libEnergiesAll))])
        nSgmsAll = np.array([2.5 for i in range(0, len(libEnergiesAll))])

    scaleFactor = 0. # no Filtering
    nSgms = dict()
    samplesUp = dict()
    nSgms = dict(zip(libEnergiesAll, nSgmsAll))
    samplesUp = dict(zip(libEnergiesAll, samplesUpAll))

    #pLforTS = max(pulseLength, 2048)  # used for triggerSizeTS
    triggerSizeTS = PreBufferSize + maxPulseLength + 1000  # ___|\...........__
    print("triggerSizeTS=", triggerSizeTS)
    triggerTS3val = triggerSizeTS - PreBufferSize
    assert triggerSizeTS < triggerSizeTC, "Record length for tessim is larger than simulated"

    # -- CALIB dirs & files --
    SIMFILESdir = simSIXTEdir + "/LIBRARIES/" + tessim

    # -- LIB & NOISE dirs and files ----------
    noiseDir = simSIXTEdir + "/NOISE/" + tessim
    noiseFile = noiseDir + "/noise" + str(nsamples) + "samples_" + tessim + "_B0_100s_pairscps_" + space + ".fits"
    PixTypeFile = "'file:" + simSIXTEdir + "/newpixels.fits[" + pixName + "]'"

    #
    # Create simulated FILES for libraries and LIBRARY itself
    #
    libDirRoot = simSIXTEdir + "/LIBRARIES/" + tessim
    libDir = libDirRoot + "/GLOBAL/" + space
    libFile = libDir + "/libraryMultiE_PL" + str(pulseLength) + "_" + str(nSimPulses) + "p" + Fil + ".fits"

    #libFile += ".5-9"  # for library creation in steps

    # if os.path.isfile(libFile):
    #    os.remove(libFile)
    evttmpFile = tempfile.NamedTemporaryFile()
    #
    # Select processing space
    #
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
        monoEeV = monoEkeV * 1.E3  # eV

        print("=============================================")
        print("Adding monochromatic energy", monoEkeV, "keV")
        print("=============================================")

        print("Temporary event file in: ", evttmpFile.name)
        # simulate SIXTE file with tesconstpileup + tessim
        root0 = "single" + str(monoEkeV) + "_sep" + str(separation) + "_pix" + str(pixel) + "_" + str(nSimPulses) + "p"
        root = root0 + "_" + str(pulseLength) + noise
        pixFile = PIXIMPACTdir + "/" + root0 + "_trSz" + str(triggerSizeTC) + ".piximpact"
        simFile = SIMFILESdir + "/" + root + ".fits"

        # -- TESCONSTPILEUP: generate impacts for well separated pairs of pulses --
        # (used her to simulate single pulses)
        # -- TESSIM: simulate well separated pairs of pulses --
        if not os.path.isfile(simFile):
            print("Simulating & triggering pulses with TESSIM to ", simFile)
            if not os.path.isfile(pixFile):
                comm = ("tesconstpileup PixImpList=" + pixFile + " XMLFile=" + XMLfile + " tstop=" + simTime +
                        " energy=" + str(monoEkeV) + " pulseDistance=" + str(separation) + " TriggerSize=" +
                        str(triggerSizeTC) + " clobber=yes")
                print("Generating PIXIMPACT ", pixFile, " running:\n", comm)
                try:
                    args = shlex.split(comm)
                    check_call(args, stderr=STDOUT)
                except:
                    print("Error running TESCONSTPILEUP for piximpact list generation")
                    os.chdir(cwd)
                    shutil.rmtree(tmpDir)
                    raise

            comm = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFile + " Streamfile=" + simFile +
                    " tstart=0" + " tstop=" + simTime + " triggerSize=" + str(triggerSizeTS) + " preBuffer=" +
                    str(PreBufferSize) + " acbias=" + acbias + " triggertype='diff:3:100:" + str(triggerTS3val) +
                    "'" + " sample_rate=" + samprate + " PixType=" + PixTypeFile)
            
            print("Running tessim for simulation\n", comm)
            try:
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except:
                print("Error running TESSIM for simulation\n")
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise
                
            # -- rm first (and LAST) record and update NETTOT (first starts high and last can be cut)
            fsim = fits.open(simFile)
            nrows = fsim[1].header["NAXIS2"]
            fsim.close()
            assert nrows > 1, "Tessim failed for (%s): just one row present " % simFile

            try:
                comm = "fdelrow infile=" + simFile + "+1 firstrow=" + str(nrows) + " nrows=1 confirm=no proceed=yes"
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
                comm = "fdelrow infile=" + simFile + "+1 firstrow=1 nrows=1 confirm=no proceed=yes"
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
                fsim = fits.open(simFile, mode='update')
                nrows2 = fsim[1].header['NAXIS2']
                assert nrows2 == nrows-2, "Failure removing initial & last rows in (%s): " % simFile
                nettot = nrows2 * 2  # new number of pulses
                fsim[1].header['NETTOT'] = nettot
                fsim.close()
            except:
                print("Error running FTOOLS to remove initial & last rows in ", simFile)
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise
        continue  # (to just simulate pulses files)
        # -- Create/update LIBRARY --

        comm = ("tesreconstruction Recordfile=" + simFile + " TesEventFile=" + evttmpFile.name + " Rcmethod=SIRENA" +
                " PulseLength=" + str(pulseLength) + " LibraryFile=" + libFile + " scaleFactor=" + str(scaleFactor) +
                " samplesUp=" + str(samplesUp[monoEkeV]) + " nSgms=" + str(nSgms[monoEkeV]) +
                " mode=0 clobber=yes intermediate=0" + " monoenergy=" + str(monoEeV) + " EventListSize=1000" +
                " tstartPulse1=" + str(tstartPulse1[monoEkeV]) + " tstartPulse2=" + str(tstartPulse2[monoEkeV]) +
                " tstartPulse3=" + str(tstartPulse3) + " NoiseFile=" + noiseFile + " XMLFile=" +
                XMLfile + energyMethod)
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

    parser.add_argument('--pixName', help='Extension name in FITS pixel definition file (SPA*, LPA1*, LPA2*, LPA3*)')
    parser.add_argument('--space', choices=['ADC', 'R', 'RALL', 'RNOL', 'RFITTED'],
                        help='Input Data Space  (ADC, R, RALL, RNOL, RFITTED)')
    parser.add_argument('--pulseLength', type=int, help='pulse length samples')
    parser.add_argument('--nsamples', type=int, help='noise samples')
    parser.add_argument('--nSimPulses', type=int, help='Number of Pulses in simulated files')
    parser.add_argument('--acbias', default='yes', choices=['yes', 'no'],
                        help='Operating Current (AC )(acbias=yes) or DC (acbias=no)) [default %(default)s]')

    args = parser.parse_args()

    simulGlobalLibs(pixName=args.pixName, space=args.space, pulseLength=args.pulseLength,
                    nsamples=args.nsamples, nSimPulses=args.nSimPulses, acbias=args.acbias)
