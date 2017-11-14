"""
# PAIRS simulation
#
#      |\     |\                         |\     |\                         |\     |\
#   ___| \____| \________________________| \____| \________________________| \____| \____________
#           record                             record                             record
#   <----------------------->         <----------------------->         <----------------------->
#
# python simulPairs.py
#
#  Input parameters:
#          array (SPA|LPA1|LPA2|LPA3)
#          monoEkeV1: monochromatic energy of pulses 1(in keV)
#          monoEkeV2: monochromatic energy of pulses 2(in keV)
#          ACDC: AC or DC
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
from subprocess import check_call, check_output, STDOUT
from astropy.io import fits
import xml.etree.ElementTree as ET

cwd = os.getcwd()
tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"

nSimPulses = 20000  # 10000 records = 10000 secondary pulses
nSimPulses = 2000  # 1000 records = 1000 secondary pulses
recordSeparation = 40000  # separation from secondary --> primary for next record

#simSIXTEdir = "/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE"
#PAIRSdir = "/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS"
#XMLdir = os.environ["SIXTE"] + "/" + "share/sixte/instruments/athena/1469mm_xifu"
simSIXTEdir = "/dataj6/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE"
PAIRSdir = "/dataj6/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS"
XMLdir = os.environ["SIXTE"] + "/" + "share/sixte/instruments/athena/1469mm_xifu"
# XMLfile = XMLdir + "/" + "xifu_baseline.xml"
XMLfile = XMLdir + "/" + "xifu_detector_hex_baseline.xml"
pixel = 1
PreBufferSize = 1000

XMLtree = ET.parse(XMLfile)
XMLroot = XMLtree.getroot()
for samplefreq in XMLroot.findall('samplefreq'):
    samprate = samplefreq.get('value')

tstart = 0.5/float(samprate)  # added to solve floating point inaccu. due to sampling rate (Christian's mail 31/03/2017)

triggerTH = {'LPA1shunt': 50, 'LPA2shunt': 20}


def simulPairs(pixName, monoEkeV1, monoEkeV2, acbias):
    """
    :param pixName: Extension name in the FITS pixel definition file (SPA*, LPA1*, LPA2*, LPA3*)
    :param monoEkeV1: Monochromatic energy (keV) of input first simulated pulses
    :param monoEkeV2: Monochromatic energy (keV) of input second simulated pulses
    :param acbias: Operating Current (AC if acbias=yes or DC if acbias=no)
    :return: files with simulated PAIRS
    """

    global cwd, nSimPulses, XMLfile, pixel, PreBufferSize, simSIXTEdir, samprate, triggerTH, tstart
    if monoEkeV1 == "0.5" or monoEkeV2 == "0.5":
            triggerTH["LPA2shunt"] = 50

    tessim = "tessim" + pixName
    SIMFILESdir = PAIRSdir + "/" + tessim
    PixTypeFile = "file:" + simSIXTEdir + "/newpixels.fits[" + pixName + "]"
    sepsStr = ['0']
    pulseLength = 0
    if "SPA" in pixName:

        sepsStr = ['00004', '00005', '00007', '00010', '00013', '00017', '00023', '00031', '00042', '00056', '00075',
                   '00101', '00136', '00182', '00244', '00328', '00439', '00589', '00791', '01061', '01423', '01908']
        pulseLength = 1024  # only to calculate triggerSize

    elif "LPA1" in pixName:

        sepsStr = ['00004', '00005', '00007', '00010', '00013', '00017', '00023', '00031', '00042', '00056', '00075',
                   '00101', '00136', '00182', '00244', '00328', '00439', '00589', '00791', '01061', '01423', '01908'
                   ]
        pulseLength = 2048  # only to calculate triggerSize

    elif ("LPA2" in pixName) or ("LPA3" in pixName):

        #sepsStr = ['00004', '00006', '00008', '00011', '00044', '00060', '00090', '00100', '00120', '00200',
        #           '00250', '00300', '00400', '00500', '00600', '00800', '01000', '01600', '02000']

        sepsStr = ['{0:05d}'.format(member) for member
                   in list(map(int, numpy.ndarray.tolist(numpy.logspace(numpy.log10(4), numpy.log10(2000), num=20))))]

        pulseLength = 4096  # only to calculate triggerSize

    for sepA in sepsStr:
        sep12 = int(sepA)
        triggerSizeTC = PreBufferSize + sep12 + recordSeparation + 1000
        triggerSizeTS = PreBufferSize + sep12 + pulseLength + 1000
        triggerTS3val = triggerSizeTS - PreBufferSize
        triggerSizeTC = int(triggerSizeTC)

        # calculate sim time to have at least nSimPulses pulses:
        #   simTime= (nSimPulses/2)recs * (triggerSize sam/rec)/(samprate sam/s)
        simTime = nSimPulses/2. * triggerSizeTC/float(samprate)
        simTime = '{0:0.0f}'.format(simTime)

        root0 = "sep" + sepA + "sam_" + simTime + "s_" + monoEkeV1 + "keV_" + monoEkeV2 + "keV_jitter"  # for piximpact
        root = "sep" + sepA + "sam_" + str(nSimPulses) + "p_" + monoEkeV1 + "keV_" + monoEkeV2 + "keV_jitter"  # for fits
        pixFile = cwd + "/PIXIMPACT/" + root0 + "_trSz" + str(triggerSizeTC) + ".piximpact"
        fitsFile = SIMFILESdir + "/" + root + ".fits"
        print("-------------------------------------------\n")
        print("Simulating ", fitsFile, "\n")
        print("-------------------------------------------\n")

        if not os.path.isfile(pixFile):
            comm = ("tesconstpileup PixImpList=" + pixFile + " XMLFile=" + XMLfile + " timezero=" + str(tstart) +
                    " tstop=" + str(simTime) + " energy=" + monoEkeV1 + " energy2=" + monoEkeV2 +
                    " pulseDistance=" + str(sep12) + " TriggerSize=" + str(triggerSizeTC) + " clobber=yes")
            print("\n##### Runing tesconstpileup #########")
            print(comm, "\n")
            try:
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except:
                print("Error running TESCONSTPILEUP for piximpact list generation")
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise
            # continue  # to simulate only piximpact files

        if not os.path.isfile(fitsFile):
            comm = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFile + " Streamfile=" + fitsFile + " tstart=0." +
                    " tstop=" + simTime + " triggerSize=" + str(triggerSizeTS) + " preBuffer=" +
                    str(PreBufferSize) + " triggertype='diff:3:" + str(triggerTH[pixName]) + ":" + str(triggerTS3val) +
                    "'" + " acbias=" + acbias + " sample_rate=" + samprate + " PixType=" + PixTypeFile)

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
                nettot = nrows2 * 2  # new number of pulses
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

    parser = argparse.ArgumentParser(description='Simulate pairs of pulses', prog='simulPairs')

    parser.add_argument('--pixName', help='Extension name in pixel definition FITS file (SPA*, LPA1*, LPA2*, LPA3*)')
    parser.add_argument('--monoEnergy1', help='Monochromatic energy (keV) of input simulated first pulse')
    parser.add_argument('--monoEnergy2', help='Monochromatic energy (keV) of input simulated second pulse')
    parser.add_argument('--acbias', choices=['yes', 'no'],
                        help='Operating Current (acbias=yes for AC or acbias=no for DC)')

    inargs = parser.parse_args()
    simulPairs(pixName=inargs.pixName, monoEkeV1=inargs.monoEnergy1,
               monoEkeV2=inargs.monoEnergy2, acbias=inargs.acbias)
