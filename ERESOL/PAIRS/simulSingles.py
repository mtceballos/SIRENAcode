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
#          monoEkeV: monochromatic energy of pulses 1(in keV)
#          ACDC: AC or DC
#          samprate: "" or "samprate2"
#          jitter:  "" or "jitter"
#
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
import numpy
import auxpy
from time import gmtime, strftime
from subprocess import check_call, check_output, STDOUT
from astropy.io import fits
import xml.etree.ElementTree as ET

cwd = os.getcwd()
tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"

nSimPulses = 20000  # 10000 records = 10000 secondary pulses
singleSeparation = 40000  # separation from secondary --> primary for next record

simSIXTEdir = "/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE"
PAIRSdir = "/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS"
XMLdir = os.environ["SIXTE"] + "/" + "share/sixte/instruments/athena/1469mm_xifu"

pixel = 1
PreBufferSize = 1000
pulseLength = 8192  # only to calculate triggerSize
# with these thresholds, 0.2keV pulses are triggered at 999, 0.5keV @999 or @1000 and larger pulses @10000 
triggerTH = {'LPA1shunt': 50, 'LPA2shunt': 20} 


def simulSingles(pixName, monoEkeV, acbias, samprate, jitter):
    """
    :param pixName: Extension name in the FITS pixel definition file (SPA*, LPA1*, LPA2*, LPA3*)
    :param monoEkeV: Monochromatic energy (keV) of input simulated pulses
    :param acbias: Operating Current (AC if acbias=yes or DC if acbias=no)
    :param samprate: Samprate value with respect to baseline of 156250 Hz: "" (baseline), "samprate2" (half_baseline)
    :param jitter: jitter option ("" for no_jitter and "jitter" for jitter)
    :return: files with simulated SINGLES
    """

    global cwd, nSimPulses, XMLdir, pixel, PreBufferSize, simSIXTEdir, triggerTH, pulseLength
    if monoEkeV == "0.5":
            triggerTH["LPA2shunt"] = 50  # increase threshold so that they are all triggered @999 (if no jitter...)

    # samprate
    smprtStr = ""
    if samprate == 'samprate2':
        smprtStr = "_samprate2"
        triggerTH["LPA2shunt"] = 25  # increase threshold so that they are all triggered @999 (if no jitter...)
        if monoEkeV == "0.5":
            triggerTH["LPA2shunt"] = 60  # increase threshold so that they are all triggered @999 (if no jitter...)

    XMLfile = "/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/xifu_detector_hex_baselineNEWgrades" + smprtStr + ".xml"
    XMLtree = ET.parse(XMLfile)
    XMLroot = XMLtree.getroot()
    for samplefreq in XMLroot.findall('samplefreq'):
        samprate = samplefreq.get('value')
    # added to solve floating point inaccu. due to sampling rate (Christian's mail 31/03/2017)
    dtaums = int(singleSeparation) / float(samprate) * 1000.  # separation time (ms) between pulses
    tstart = 0.5 / float(samprate)

    # jitter
    jitterStr = ""
    offset = ""
    if jitter == "jitter":
        jitterStr = "_jitter"
        offset = " offset=-1"


    tessim = "tessim" + pixName
    SIMFILESdir = PAIRSdir + "/" + tessim
    # PixTypeFile = "file:" + simSIXTEdir + "/newpixels.fits[" + pixName + "]"  # pre-Feb 2018 telecon
    PixTypeFile = "file:" + simSIXTEdir + "/newpix_full.fits[" + pixName + "]"

    triggerSizeTC = PreBufferSize + singleSeparation + singleSeparation + 1000
    triggerSizeTS = PreBufferSize + pulseLength + 1000
    triggerTS3val = triggerSizeTS - PreBufferSize
    triggerSizeTC = int(triggerSizeTC)

    # calculate sim time to have at least nSimPulses pulses:
    #   simTime= (nSimPulses/2)recs * (triggerSize sam/rec)/(samprate sam/s)
    simTime = nSimPulses/2. * triggerSizeTC/float(samprate)
    simTime = '{0:0.0f}'.format(simTime)

    # for piximpact:
    root0 = "sep" + str(singleSeparation) + "sam_" + simTime + "s_" + monoEkeV + "keV" + smprtStr + jitterStr
    # for fits:
    root = "sep" + str(singleSeparation) + "sam_" + str(nSimPulses) + "p_" + monoEkeV + "keV" + smprtStr + jitterStr
    pixFile = cwd + "/PIXIMPACT/" + root0 + "_trSz" + str(triggerSizeTC) + ".piximpact"
    fitsFile = SIMFILESdir + "/" + root + ".fits"
    print("-------------------------------------------\n")
    print("Simulating ", fitsFile, "\n")
    print("-------------------------------------------\n")

    if not os.path.isfile(pixFile):
        if jitterStr == "":
            comm = ("tesgenimpacts PixImpList=" + pixFile + " mode=const tstart=" + str(tstart) + " tstop=" + simTime +
                    " EConst=" + monoEkeV + " dtau=" + str(dtaums) + " clobber=yes")
        else:
            comm = ("tesconstpileup PixImpList=" + pixFile + " XMLFile=" + XMLfile + " timezero=" + str(tstart) +
                    " tstop=" + str(simTime) + " energy=" + monoEkeV + offset +
                    " pulseDistance=" + str(singleSeparation) + " TriggerSize=" + str(triggerSizeTC) + " clobber=yes")

        print("\n##### Runing tesconstpileup/tesgenimpacts #########")
        print(comm, "\n")
        try:
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running tool for piximpact list generation")
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
            raise
        # continue  # to simulate only piximpact files

    if not os.path.isfile(fitsFile):
        commTessim = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFile + " Streamfile=" + fitsFile +
                      " tstart=0." + " tstop=" + simTime + " triggerSize=" + str(triggerSizeTS) + " preBuffer=" +
                      str(PreBufferSize) + " triggertype='diff:3:" + str(triggerTH[pixName]) + ":" +
                      str(triggerTS3val) + "'" + " acbias=" + acbias + " sample_rate=" + samprate +
                      " PixType=" + PixTypeFile)

        print("\n##### Runing tessim #########")
        print(commTessim, "\n")
        try:
            args = shlex.split(commTessim)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running TESSIM for data simulation")
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
            raise
        # continue
        # update HISTORY  for tesconstpileup run in header[0]
        auxpy.updateHISTORY(fitsFile, commTessim)

        # rm first (and LAST) record and update NETTOT
        fsim = fits.open(fitsFile)
        nrows = fsim[1].header["NAXIS2"]
        assert nrows > 1, "Tessim failed: just one huge row present!"
        fsim.close()
        # continue
        try:
            print("Removing first & last row, just in case, and updating NETTOT")
            auxpy.rmLastAndFirst(fitsFile, 2)
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
    parser.add_argument('--monoEnergy', help='Monochromatic energy (keV) of input simulated pulse')
    parser.add_argument('--acbias', choices=['yes', 'no'], default="yes",
                        help='Operating Current (acbias=yes for AC or acbias=no for DC)')
    parser.add_argument('--samprate', default="", choices=['', 'samprate2'], help="baseline, half_baseline")
    parser.add_argument('--jitter', default="", choices=['', 'jitter'], help="no jitter, jitter")

    inargs = parser.parse_args()
    simulSingles(pixName=inargs.pixName, monoEkeV=inargs.monoEnergy,
               acbias=inargs.acbias,
               samprate=inargs.samprate, jitter=inargs.jitter)
