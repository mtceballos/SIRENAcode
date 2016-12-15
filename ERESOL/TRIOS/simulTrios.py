"""
# TRIOS simulation
#
# python simulTrios.py
#
#  Input parameters:
#          pixType (SPA|LPA1*|LPA2*|LPA3*)
#          monoEkeV: monochromatic energy (in keV)
#          acbias: AC or DC
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
from subprocess import check_call, check_output, STDOUT
from astropy.io import fits
import xml.etree.ElementTree as ET

cwd = os.getcwd()
tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"

nSimPulses = 10000
simSIXTEdir = "/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE"
TRIOSdir = "/home/ceballos/INSTRUMEN/EURECA/ERESOL/TRIOS"
# XMLrect = simSIXTEdir + "/CALIBRATION/newFormat/xifu_detector_300us_156kHz_rect.xml"
XMLdir = os.environ["SIXTE"] + "/" + "share/sixte/instruments/athena/1469mm_xifu"
# XMLfile = XMLdir + "/" + "xifu_baseline.xml"
XMLfile = XMLdir + "/" + "xifu_detector_hex_baseline.xml"
pixel = 1
PreBufferSize = 1000

XMLtree = ET.parse(XMLfile)
XMLroot = XMLtree.getroot()
for samplefreq in XMLroot.findall('samplefreq'):
    samprate = samplefreq.get('value')


def simulTrios(pixType, monoEkeV, acbias):
    """
    :param pixType: Extension name in the FITS pixel definition file (SPA*, LPA1*, LPA2*, LPA3*)
    :param monoEkeV: Monochromatic energy (keV) of input simulated pulses
    :param acbias: Operating Current (AC if acbias=yes or DC if acbias=no)
    :return: files with simulated TRIPLETS
    """

    global cwd, nSimPulses, XMLfile, pixel, samprate, PreBufferSize, simSIXTEdir
    tessim = "tessim" + pixType
    SIMFILESdir = TRIOSdir + "/" + tessim
    # PixType = "file:" + simSIXTEdir + "/newpixels.fits[" + array + ACDC + "]"
    PixTypeFile = "file:" + simSIXTEdir + "/newpixels.fits[" + pixType + "]"

    thres = 100
    if (monoEkeV < 0.5):
        thres = 80

    if "SPA" in pixType:

        sepsStr = ['00004', '00005', '00007', '00010', '00013', '00017', '00023', '00031', '00042', '00056', '00075',
                   '00101', '00136', '00182', '00244', '00328', '00439', '00589', '00791', '01061', '01423', '01908']
        pulseLength = 1024  # only to calculate triggerSize

    elif "LPA1" in pixType:

        sepsStr = ['00004', '00005', '00007', '00010', '00013', '00017', '00023', '00031', '00042', '00056', '00075',
                   '00101', '00136', '00182', '00244', '00328', '00439', '00589', '00791', '01061', '01423', '01908',
                   '02560', '03433', '04605', '06178', '08287']
        pulseLength = 2048  # only to calculate triggerSize

    elif ("LPA2" in pixType) or ("LPA3" in pixType):

        sepsStr = ['00004', '00005', '00007', '00010', '00013', '00017', '00023', '00031', '00042', '00056', '00075',
                   '00101', '00136', '00182', '00244', '00328', '00439', '00589', '00791', '01061', '01423', '01908',
                   '02560', '03433', '04605', '06178', '08287', '11115', '14910', '20000']

        pulseLength = 2048  # only to calculate triggerSize

    for sepA in sepsStr:
        sep12 = int(sepA)

        for sepB in sepsStr:
            sep23 = int(sepB)
            triggerSizeTC = PreBufferSize + sep12 + sep23 + 20000 + PreBufferSize + 1000
            triggerSizeTS = PreBufferSize + sep12 + sep23 + pulseLength + 1000
            triggerTS3val = triggerSizeTS - PreBufferSize
            triggerSizeTC = int(triggerSizeTC)

            # calculate sim time to have at least nSimPulses pulses:
            #   simTime= (nSimPulses/3)recs * (triggerSize sam/rec)/(samprate sam/s)
            simTime = nSimPulses/3. * triggerSizeTC/float(samprate)
            simTime = '{0:0.0f}'.format(simTime)

            root0 = "sep" + sepA + "sep" + sepB + "sam_" + simTime + "s_" + str(monoEkeV) + "keV"  # for piximpact
            root = "sep" + sepA + "sep" + sepB + "sam_" + str(nSimPulses) + "p_" + str(monoEkeV) + "keV"  # for fits
            pixFile = cwd + "/PIXIMPACT/" + root0 + "_trSz" + str(triggerSizeTC) + ".piximpact"
            fitsFile = SIMFILESdir + "/" + root + ".fits"
            print("-------------------------------------------\n")
            print("Simulating ", fitsFile, "\n")
            print("-------------------------------------------\n")

            if not os.path.isfile(pixFile):
                comm = ("tesconstpileup PixImpList=" + pixFile + " XMLFile=" + XMLfile + " tstop=" + str(simTime) +
                        " energy=" + str(monoEkeV) + " pulseDistance=" + str(sep12) + " pulseDistance2=" + str(sep23) +
                        " energy2=" + str(monoEkeV) + " energy3=" + str(monoEkeV) + " TriggerSize=" +
                        str(triggerSizeTC) + " clobber=yes")
                print("\n##### Runing tesconstpileup #########")
                print(comm,"\n")
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
                # comm = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFile + " Streamfile=" + fitsFile +
                #        " tstart=0. tstop=" + simTime + " sample_rate=" + str(samprate) + " triggerSize=" +
                #        str(triggerSizeTS) + " preBuffer=" + str(PreBufferSize) + " triggertype='movavg:5:1.1:0' " +
                #        " PixType=" + PixType)
                comm = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFile + " Streamfile=" + fitsFile +
                        " tstart=0. tstop=" + simTime + " triggerSize=" + str(triggerSizeTS) + " preBuffer=" +
                        str(PreBufferSize) + " triggertype='diff:3:" + str(thres) + ":" + str(triggerTS3val) +
                        "' acbias=" + acbias + " sample_rate=" + samprate + " PixType=" + PixTypeFile)
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
                    comm = "fdelrow infile=" + fitsFile + "+1 firstrow=" + str(nrows) + " nrows=1 confirm=no proceed=yes"
                    args = shlex.split(comm)
                    check_call(args, stderr=STDOUT)
                    comm = "fdelrow infile=" + fitsFile + "+1 firstrow=1 nrows=1 confirm=no proceed=yes"
                    args = shlex.split(comm)
                    check_call(args, stderr=STDOUT)
                    fsim = fits.open(fitsFile, mode='update')
                    nrows2 = fsim[1].header['NAXIS2']
                    assert nrows2 == nrows-2, "Failure removing initial & last rows in (%s): " % fitsFile
                    nettot = nrows2 * 3  # new number of pulses
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

    parser = argparse.ArgumentParser(description='Simulate triplets of pulses', prog='simulTrios')

    parser.add_argument('--pixType', help='Extension name in pixel definition FITS file (SPA*, LPA1*, LPA2*, LPA3*)')
    parser.add_argument('--monoEnergy', default=1,  help='Monochromatic energy (keV) of input simulated pulses')
    parser.add_argument('--acbias', choices=['yes', 'no'],
                        help='Operating Current (acbias=yes for AC or acbias=no for DC)')

    args = parser.parse_args()
    simulTrios(pixType=args.pixType, monoEkeV=args.monoEnergy, acbias=args.acbias)

