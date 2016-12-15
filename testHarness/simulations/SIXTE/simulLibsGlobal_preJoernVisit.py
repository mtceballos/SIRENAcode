"""
CREATE GLOBAL LIBRARY for simulated pairs of pulses for calibration

python simulGlobalLibs.py

"""

# ----IMPORT MODULES --------------
from __future__ import print_function
import os
# import math
import shlex
import shutil
import tempfile
# import numpy as np
from subprocess import check_call, check_output, STDOUT
from astropy.io import fits
# from astropy.table import Table, Column

# ----GLOBAL VARIABLES -------------
simTime = 3000
separation = 20000

PreBufferSize = 1000
# tstartPulse1 = int(PreBufferSize - 1)  # 1000 for fv, 999 for GSL
tstartPulse1 = int(PreBufferSize - 1 - 1)  # (there is currently a problem for primary pulses)
# tstartPulse2 = int(tstartPulse1 + separation)  # 1000 + sep for fv; 1000+sep-1 for GSL
tstartPulse2 = int(tstartPulse1 + 1 + separation)  # (there is currently a problem for primary pulses)
tstartPulse3 = 0
triggerSizeTC = PreBufferSize + separation + separation + PreBufferSize + 1000
triggerSizeTS = PreBufferSize + separation
noise = ""
Fil = ""
simnoise = "simnoise=y"
pixel = 1
samprate = 156250  # Hz
libEnergies = [0.1, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15.1, 16]
#libEnergies = [1]
lastEnergy = libEnergies[-1]

nSimPulses = 5000  # minimum number of simulated pulses in simulPairs.csh (approx.; to get filenames)
cwd = os.getcwd()
tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
if tstartPulse1 > 0:
    nSgms = 0
    samplesUp = 0
    scaleFactor = 0


def simulGlobalLibs(array, space, pulseLength, ACDC):
    """
    :type array: str
    :param array: Array acronym (SPA, LPA1, LPA2, LPA3)
    :type space: str
    :param space: Input Data Space (ADC, R, RB)
    :type pulseLength: int
    :param pulseLength: Pulse length & noise length samples
    :type ACDC: str
    :param ACDC: Operating Current (AC or DC)
    :return: simulated calibration pairs of pulses && Global library from them
    """

    global separation, simTime, Fil, nSgms, samplesUp, scaleFactor, cwd, tstartPulse2, tstartPulse1
    tessim = "tessim" + array
    ACDC = ACDC.lower()
    if space == "RB":  # In I2RBIS, pulse shape is a derivative and pulse template is incomplete otherwise
        tstartPulse1 -= 1
        tstartPulse2 -= 1

    # -- CALIB dirs & files --
    simSIXTEdir = "/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE"
    XMLrect = simSIXTEdir + "/CALIBRATION/newFormat/xifu_detector_300us_156kHz_rect.xml"
    PIXIMPACTdir = simSIXTEdir + "/LIBRARIES/PIXIMPACT"
    SIMFILESdir = simSIXTEdir + "/LIBRARIES/" + tessim
    PixType = "file:" + simSIXTEdir + "/newpixels.fits[" + array + ACDC + "]"

    # -- LIB & NOISE dirs and files ----------
    noiseDir = simSIXTEdir + "/NOISE/" + tessim
    noiseFile = noiseDir + "/noise" + str(pulseLength) + "samples_" + tessim + "_B0_100s_pairscps_" + space + ".fits"
    # --- get baseline from noise file (in current) --- NOT REQUIRED ANYMORE
    # noisef = fits.open(noiseFile)
    # baseline = noisef[1].header["BASELINE"]
    # noisef.close()

    #
    # Create LIBRARIES and FILES for libraries
    #
    libDirRoot = simSIXTEdir + "/LIBRARIES/" + tessim
    libDir = libDirRoot + "/GLOBAL/" + space
    libFile = libDir + "/libraryMultiE_GLOBAL_PL" + str(pulseLength) + "_" + tessim + Fil + ".fits"
    if os.path.isfile(libFile):
        os.remove(libFile)
    evttmpFile = libDir + "/evtcal.fits"

    #
    # Select processing space
    #
    if space == "R":
        energyMethod = "EnergyMethod=I2R"
    elif space == "RB":
        energyMethod = "EnergyMethod=I2RBIS"
    else:
        energyMethod = ""

    # --- Simulate and Process calibration data files -------
    for monoEkeV in libEnergies:  # keV
        monoEeV = monoEkeV * 1.E3  # eV
        # lastElibrary = 0
        # if monoEkeV == lastEnergy:
        #    lastElibrary = 1
        print("=============================================")
        print("Adding monochromatic energy", monoEkeV, "keV")
        print("=============================================")

        # modify if necessary also taking into account different energies (see .csh version)
        if (array == "SPA") or (array == "LPA1"):
            samplesUp = 3
            nSgms = 20
            scaleFactor = 0.005
        elif array == "LPA2":
            samplesUp = 2
            scaleFactor = 0.005
        elif array == "LPA3":
            samplesUp = 2
            scaleFactor = 0.02
            Fil = "Fil"

        if tstartPulse1 > 0:
            nSgms = 0
            samplesUp = 0
            scaleFactor = 0
            Fil = ""

        # simulate SIXTE file with tesconstpileup + tessim
        root0 = "mono" + str(monoEkeV) + "_sep" + str(separation) + "_pix" + str(pixel) + "_" + str(simTime) + "s"
        root = root0 + "_" + str(pulseLength) + noise
        pixFile = PIXIMPACTdir + "/" + root0 + "_trSz" + str(triggerSizeTC) + ".piximpact"
        simFile = SIMFILESdir + "/" + root + ".fits"
        # streamFile = SIMFILESdir + "/" + root + ".stream"

        # -- TESCONSTPILEUP: generate impacts for well separated pairs of pulses --

        # -- TESSIM: simulate well separated pairs of pulses --
        if not os.path.isfile(simFile):
            print("Simulating & triggering pulses with TESSIM to ", simFile)
            if monoEkeV == 15 and 0:
                tmplist = open('tmp.list', 'w')
                # do a trick (tessim stream files in 30s fractions + streamtotriggers) for 15 keV files
                simTimeNleft = simTime
                tstart = 0.
                tstop = min(30., simTimeNleft)
                i = 0
                # tmpFiles = ""
                while simTimeNleft > 0:
                    simTimeNleft -= 30.
                    i += 1
                    comm = ("tesconstpileup PixImpList=tmp.piximpact XMLFile=" + XMLrect + " tstop=" + str(tstop) +
                            " energy=" + str(monoEkeV) + " pulseDistance=" + str(separation) + " TriggerSize=" +
                            str(triggerSizeTC) + " clobber=yes")
                    print("Generating (<30s) PIXIMPACT for 15 keV tmp.piximpact running:\n", comm)
                    try:
                        args = shlex.split(comm)
                        check_call(args, stderr=STDOUT)
                    except:
                        print("Error running TESCONSTPILEUP for piximpact list generation (<30s)")
                        os.chdir(cwd)
                        shutil.rmtree(tmpDir)
                        raise

                    tmpFile = "tmp.fits" + str(i)
                    comm = ("tessim PixID=" + str(pixel) + " PixImpList=tmp.piximpact Streamfile=tmp.stream" +
                            " tstart=" + str(tstart) + " tstop=" + str(tstop) + " sample_rate=" + str(samprate) +
                            " acbias=yes triggertype=stream prebuffer=" + str(PreBufferSize) +
                            " clobber=yes PixType=" + PixType)
                    args = shlex.split(comm)
                    try:
                        check_call(args, stderr=STDOUT)
                    except:
                        print("Error running tool tessim for 15 keV files")
                        print(comm)
                        os.chdir(cwd)
                        shutil.rmtree(tmpDir)
                        raise
                    print("\nSTREAMTOTRIGGERS: Generating trigger file for tmp stream...")
                    comm = ("streamtotriggers PixImpList=tmp.piximpact XMLFile=" + XMLrect + " tstart=0. tstop=" +
                            str(tstop) + " Streamfile=tmp.stream" + " TesTriggerFile=" + tmpFile +
                            " TriggerSize=" + str(triggerSizeTC) + " PreBufferSize=" + str(PreBufferSize) + " pixels=" +
                            str(pixel) + " clobber=yes")
                    args = shlex.split(comm)
                    try:
                        check_call(args, stderr=STDOUT)
                    except:
                        print("Error running tool streamtotriggers (for tmp 15 keV):")
                        print(comm)
                        # print(sys.exc_info()[0])
                        os.chdir(cwd)
                        shutil.rmtree(tmpDir)
                        raise
                    print("STREAMTOTRIGGERS for 15keV: .................................END")
                    # tstart = tstop
                    # tstop = tstop + min(30, simTimeNleft)
                    tstop = min(30, simTimeNleft)
                    tmplist.write('{0}\n'.format(tmpFile))

                tmpFiles = tmpFiles + " " + tmpFile
                # fmerge or ftmerge do not work with var-length arrays...: stilts tcat + cphead to recover header
                # comm = ("fmerge infiles=@tmp.list outfile=" + simFile + " columns=- clobber=yes")
                try:
                    comm = "stilts tcat in=@tmp.list out=" + simFile
                    args = shlex.split(comm)
                    check_call(args, stderr=STDOUT)
                    comm = "cphead " + tmpFile + "+1 " + simFile + "+1"
                    args = shlex.split(comm)
                    check_call(args, stderr=STDOUT)
                except:
                    print("Error running tool fmerge (15 keV):")
                    print(comm)
                    os.chdir(cwd)
                    shutil.rmtree(tmpDir)
                    os.remove('tmp.list')
                    raise
                print("TESSIM + STREAMTOTRIGGERS for 15 keV: ................................END")
                # os.remove("tmp.stream")
                os.remove("tmp.piximpact")
            else:
                if not os.path.isfile(pixFile):
                    comm = ("tesconstpileup PixImpList=" + pixFile + " XMLFile=" + XMLrect + " tstop=" + str(simTime) +
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
                        " tstart=0. tstop=" + str(simTime) + " sample_rate=" + str(samprate) + " triggerSize=" +
                        str(triggerSizeTS) + " preBuffer=" + str(PreBufferSize) + " triggertype='movavg:5:1.1:0' " +
                        " PixType=" + PixType)
                print("Running tessim for simulation:\n", comm)
                try:
                    args = shlex.split(comm)
                    check_call(args, stderr=STDOUT)
                except:
                    print("Error running TESSIM for simulation")
                    os.chdir(cwd)
                    shutil.rmtree(tmpDir)
                    raise

            # -- rm first (and LAST) record and update NETTOT (first starts high and last can be cut)
            fsim = fits.open(simFile)
            nrows = fsim[1].header["NAXIS2"]
            fsim.close()

            # from http://docs.astropy.org/en/stable/io/unified.html#table-io (new AstroPy doc)
            # but it cannot operate over variable length array cols
            # tsim = Table.read(simFile)
            # tsim.remove_rows([0,nrows-1])
            # tsim.write(simFile, overwrite=True)

            try:
                comm = "fdelrow infile=" + simFile + "+1 firstrow=" + str(nrows) + " nrows=1 confirm=no proceed=yes"
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
                comm = "fdelrow infile=" + simFile + "+1 firstrow=1 nrows=1 confirm=no proceed=yes"
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
                fsim = fits.open(simFile)
                nrows2 = fsim[1].header['NAXIS2']
                assert nrows2 == nrows-2, "Failure removing initial & last rows in (%s): " % simFile
                nettot = nrows2 * 2 # new number of pulses
                fsim[1].header['NETTOT'] = nettot
                fsim.close()
            except:
                print("Error running FTOOLS to remove initial & last rows in ", simFile)
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise
        # continue
        # -- Create/update LIBRARY --
        if os.path.isfile(evttmpFile):
            os.remove(evttmpFile)
        # use only secondary pulses to avoid malformated primary ones in tessim
        # comm = ("tesreconstruction Recordfile=" + simFile + " TesEventFile=" + evttmpFile + " Rcmethod=SIRENA" +
        #         " PulseLength=" + str(pulseLength) + " LibraryFile=" + libFile + " scaleFactor=" + str(scaleFactor) +
        #         " samplesUp=" + str(samplesUp) + " nSgms=" + str(nSgms) + " mode=0 clobber=yes intermediate=0 " +
        #         " monoenergy=" + str(monoEeV) + " baseline=" + str(baseline) + " EventListSize=1000" +
        #         " lastElibrary=" + str(lastElibrary) + " tstartPulse1=" + str(tstartPulse1) + " tstartPulse2=" +
        #         str(tstartPulse2) + " tstartPulse3=" + str(tstartPulse3) + " NoiseFile=" + noiseFile + " OFInterp=DAB")
        comm = ("tesreconstruction Recordfile=" + simFile + " TesEventFile=" + evttmpFile + " Rcmethod=SIRENA" +
                " PulseLength=" + str(pulseLength) + " LibraryFile=" + libFile + " scaleFactor=" + str(scaleFactor) +
                " samplesUp=" + str(samplesUp) + " nSgms=" + str(nSgms) + " mode=0 clobber=yes intermediate=0" +
                " monoenergy=" + str(monoEeV) +  " EventListSize=1000" + " tstartPulse1=" + str(tstartPulse2) +
                " tstartPulse2=0 " + " tstartPulse3=" + str(tstartPulse3) + " NoiseFile=" + noiseFile +
                " OFInterp=DAB " + energyMethod)
        args = shlex.split(comm)
        print("SIRENA reconstruction to add a line to the library, running command:\n", comm)
        try:
            check_call(args, stderr=STDOUT)
        except:
            print("Error running SIRENA to add new line to library")
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
            raise

        # -- check that number of reconstructed pulses is the same than number of simulated pulses
        fsim = fits.open(simFile)
        nettot = fsim[1].header['NETTOT']  # number of simulated/triggered pulses
        fsim.close()
        fevt = fits.open(evttmpFile)
        ndetpulses = fevt[1].header['NAXIS2']  # number of detected pulses
        fevt.close()
        os.remove(evttmpFile)
        #assert nettot == ndetpulses, \
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

    parser.add_argument('--array', choices=['SPA', 'LPA1', 'LPA2', 'LPA3'], default='SPA',
                        help='Array acronym (SPA, LPA1, LPA2, LPA3)')
    parser.add_argument('--space', choices=['ADC', 'R', 'RB'], default='ADC',
                        help='Input Data Space  (ADC, R, RB)')
    parser.add_argument('--pulseLength', type=int, default=1024, help='pulse length & noise length samples [default %(default)s]')
    parser.add_argument('--ACDC', default='AC', choices=['AC', 'DC', 'None'],
                        help='Operating Current (AC or DC, or None for old tessim) [default %(default)s]')
    
    args = parser.parse_args()

    simulGlobalLibs(array=args.array, space=args.space, pulseLength=args.pulseLength, ACDC=args.ACDC)
