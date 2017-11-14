"""
# SINGLES simulation
#
#
#      |\       |\       |\       |\       |\     |\
#   ___| \______| \______| \______| \______| \____| \______
#    record   record   record
#   <-------><-------><------->
#
# python simulSingles.py
#
#  Input parameters:
#          array (SPA|LPA1|LPA2|LPA3)
#          monoEkeV: monochromatic energy (in keV)
#          ACDC: AC or DC
#
#   CAREFUL: use of "JITTER?"
"""

#
# --------- Read input parameters and RUN simulations -----------------
#
from __future__ import print_function
import os
import shlex
import shutil
import sys
import tempfile
from subprocess import check_call, check_output, STDOUT
from astropy.io import fits
import xml.etree.ElementTree as ET

cwd = os.getcwd()
tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"

nSimPulses = 20000
singleSeparation = "40000"
pulseLength = 4096  # only to calculate triggerSize

simSIXTEdir = "/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE"
PAIRSdir = "/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS"
XMLdir = os.environ["SIXTE"] + "/" + "share/sixte/instruments/athena/1469mm_xifu"
# XMLfile = XMLdir + "/" + "xifu_baseline.xml"
XMLfile = XMLdir + "/" + "xifu_detector_hex_baseline.xml"
pixel = 1
PreBufferSize = 1000

XMLtree = ET.parse(XMLfile)
XMLroot = XMLtree.getroot()
for samplefreq in XMLroot.findall('samplefreq'):
    samprate = samplefreq.get('value')

dtaums = int(singleSeparation) / float(samprate) * 1000.  # separation time (ms) between pulses
tstart = 0.5/float(samprate) # added to solve floating point inaccuracies due to sampling rate (Christian's mail 31/03/2017)

# With triggerTH=20, 0.2 keV pulses trigger 1 sample late (1001 instead of 1000), but ALL of them
#                    0.5 keV : some pulses trigger 1 sample late and some pulses trigger ok (1000)
#                    >= 1 keV:  ALL trigger OK
triggerTH = {'LPA1shunt': 50, 'LPA2shunt': 20}


def simulSingles(pixName, monoEkeV, acbias):
    """
    :param pixName: Extension name in the FITS pixel definition file (SPA*, LPA1*, LPA2*, LPA3*)
    :param monoEkeV: Monochromatic energy (keV) of input simulated pulses
    :param acbias: Operating Current (AC if acbias=yes or DC if acbias=no)
    :return: files with simulated PULSES
    """

    global cwd, nSimPulses, XMLfile, pixel, PreBufferSize, simSIXTEdir, samprate, triggerTH, tstart
    if monoEkeV == "0.5":
            triggerTH["LPA2shunt"] = 50

    tessim = "tessim" + pixName
    SIMFILESdir = PAIRSdir + "/" + tessim
    PixTypeFile = "file:" + simSIXTEdir + "/newpixels.fits[" + pixName + "]"

    triggerSizeTS = PreBufferSize + pulseLength + 1000  # ___|\_ooo
    triggerTS3val = triggerSizeTS - PreBufferSize

    # calculate sim time to have at least nSimPulses pulses:
    simTime = nSimPulses * int(singleSeparation) / float(samprate)  # if tesgenimpacts
    simTime = '{0:0.0f}'.format(simTime)

    root0 = "sep" + singleSeparation + "sam_" + simTime + "s_" + monoEkeV + "keV"  # for piximpact
    root = "sep" + singleSeparation + "sam_" + str(nSimPulses) + "p_" + monoEkeV + "keV"  # for fits
    pixFile = cwd + "/PIXIMPACT/" + root0 + ".piximpact"
    fitsFile = SIMFILESdir + "/" + root + "_jitter.fits"
    print("-------------------------------------------\n")
    print("Simulating ", fitsFile, "\n")
    print("-------------------------------------------\n")

    if not os.path.isfile(pixFile):
        comm = ("tesgenimpacts PixImpList=" + pixFile + " mode=const tstart=" + str(tstart) + " tstop=" + simTime +
                " EConst=" + monoEkeV + " dtau=" + str(dtaums) + " clobber=yes")
        print("\n##### Runing tesgenimpacts #########")
        print(comm, "\n")
        try:
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running tesgenimpacts for piximpact list generation")
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
            raise
        # continue  # to simulate only piximpact files

    if not os.path.isfile(fitsFile):
        comm = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFile + " Streamfile=" + fitsFile +
                " tstart=0 tstop=" + simTime + " triggerSize=" + str(triggerSizeTS) + " preBuffer=" +
                str(PreBufferSize) + " triggertype='diff:3:" + str(triggerTH[pixName]) + ":" +
                str(triggerTS3val) + "'" + " acbias=" + acbias + " sample_rate=" + samprate +
                " PixType=" + PixTypeFile)

        print("\n##### Runing tessim #########")
        print(comm, "\n")
        try:
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running TESSIM for data simulation")
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
            raise
        fsim = fits.open(fitsFile, mode='update')
        fsim[1].header["HISTORY"] = comm
        fsim.close()
        # continue
        # rm first (and LAST) record and update NETTOT
        fsim = fits.open(fitsFile)
        nrows = fsim[1].header["NAXIS2"]
        assert nrows > 1, "Tessim failed: just one huge row present!"
        fsim.close()
        # continue
        try:
            print("Removing first & last row, just in case, and updating NETTOS")
            comm = "fdelrow infile=" + fitsFile + "+1 firstrow=" + str(nrows) + " nrows=1 confirm=no proceed=yes"
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = "fdelrow infile=" + fitsFile + "+1 firstrow=1 nrows=1 confirm=no proceed=yes"
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            fsim = fits.open(fitsFile, mode='update')
            nrows2 = fsim[1].header['NAXIS2']
            assert nrows2 == nrows-2, "Failure removing initial & last rows in (%s): " % fitsFile
            nettot = nrows2  # new number of pulses
            fsim[1].header['NETTOT'] = nettot
            fsim.close()
        except:
            print("Error running FTOOLS to remove initial & last rows in ", fitsFile)
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
            raise

    # some cleaning before exiting the function
    os.chdir(cwd)
    shutil.rmtree(tmpDir)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Simulate pairs of pulses', prog='simulSingles')

    parser.add_argument('--pixName', help='Extension name in pixel definition FITS file (SPA*, LPA1*, LPA2*, LPA3*)')
    parser.add_argument('--monoEnergy', help='Monochromatic energy (keV) of input simulated pulses')
    parser.add_argument('--acbias', choices=['yes', 'no'],
                        help='Operating Current (acbias=yes for AC or acbias=no for DC)')

    inargs = parser.parse_args()
    simulSingles(pixName=inargs.pixName, monoEkeV=inargs.monoEnergy, acbias=inargs.acbias)

