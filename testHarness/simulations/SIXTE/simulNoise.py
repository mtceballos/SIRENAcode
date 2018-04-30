"""
# NOISE spectrum simulation
#
# python simulNoise.py
#
#  Input parameters:
#          pixName (SPA*|LPA1*|LPA2*|LPA3*)
#          pulseLength: length of pulses (in samples)
*          nsamples: samples for the noise
#          space: ADC (current space) or R or RALL or RNOL or RFITTED
#                                        (resistance space)
#          simTimeN: simulation time (s)
#
#
#
# 1) Simulate 10s stream with no events to calculate Baseline
#                             (pixdetillum + tessim) --> not requiered anymore
# 2) Simulate 100s stream with no events as input for gennoisespec
#                             (tesconstpileup + tessim)
# 3) Obtain noise spectrum with gennoisespec
#
"""

# ----IMPORT MODULES --------------
from __future__ import print_function
import os, auxpy
import shlex
import shutil
import sys
import tempfile
import xml.etree.ElementTree as ET
from subprocess import check_call, STDOUT
from astropy.io import fits


cpsN = "pairs"  # counts per second or pairs of pulses
preBuffer = 1000
Ifit = 45.3E-6
tmpDir = tempfile.mkdtemp()
tmpFile = tempfile.TemporaryFile()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]


# Different R transformation are requiered due to the different format in
# stream/trigger files (derivatives)

# ----MAIN routine definition ------

def simulNoise(pixName, samprate, jitter, stoch, pulseLength, space, acbias,
               scaleFactor, samplesUp, nSgms, nintervals, simTimeN, pixel):
    """ simulate data in input parameter space, calculate the data baseline
            and create Noise file to be ingested in SIRENA processing tasks

          :param pixName : name of extension in pixel definitin file
                          SPA*, LPA1*, LPA2*, LPA3*
          :param samprate: Samprate value with respect to baseline
              of 156250 Hz: "" (baseline), "samprate2" (half_baseline)
          :param jitter: jitter option ("" for no_jitter and
                                        "jitter" for jitter)
          :param stoch: stochastic option ("" for no_stoch and
                                        "stoch" for stochastic)
          :param pulseLength : pulse length in samples to select sampling
                  in noise spectrum
          :param space :  Data space: ADC for current, R or RALL or RNOL
                  or RFITTED for resistance
          :param acbias : AC (acbias=yes) or DC (acbias=no)
          :param scaleFactor : Param scaleFactor for gennoise tool
          :param samplesUp : Param samplesUp for gennoise tool
          :param nSgms : Param nSgms for gennoise tool
          :param nintervals : Number of intervals for gennoisespec
                  noise calculation
          :param simTimeN : Simulation time (s) for noise spectra calculation
          :param pixel : Pixel number
          :return fits file with Noise spectra in specified space and sampling

    """

    global samplingrate, preBuffer

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

    # stochastic
    stochStr = ""
    stochastic = " stochastic_integrator=n"
    if stoch == "stoch":
        stochStr = "_stoch"
        stochastic = (" stochastic_integrator=y dobbfb=y" +
                      " carrier_frequency=2e6 bbfb_delay=40" +
                      " decimation_filter=y")

    # ----GLOBAL VARIABLES -------------
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA"
    XMLdir = EURECAdir + "/ERESOL"
    XMLfile = (XMLdir + "/" + "xifu_detector_hex_baselineNEWgrades" +
               smprtStr + ".xml")
    XMLtree = ET.parse(XMLfile)
    XMLroot = XMLtree.getroot()
    for samplefreq in XMLroot.findall('samplefreq'):
        samplingrate = samplefreq.get('value')

    simSIXTEdir = EURECAdir + "/testHarness/simulations/SIXTE"
    cwd = os.getcwd()

    # set input params dependent variables
    # -------------------------------------
    triggerSize = max(10000, pulseLength+preBuffer)
    # PixTypeFile = "'file:" + simSIXTEdir + "/tespixels.fits[" + array + "]'"
    PixTypeFile = ("'file:" + simSIXTEdir +
                   "/newpix_full.fits[" + pixName + "]'")
    tessim = "tessim" + pixName
    wdir = "NOISE/" + tessim
    os.chdir(wdir)
    # define files
    rootN = ("forNoise" + str(pulseLength) + "samples_" + tessim +
             "_" + str(simTimeN) + "s_" + cpsN + "cps_" + space)
    noiseFile = ("noise" + str(pulseLength) + "samples_" + tessim + "_B0_" +
                 space + smprtStr + jitterStr + stochStr + ".fits")
    pixFileN = rootN + smprtStr + jitterStr + ".piximpact"
    fitsFileN = rootN + smprtStr + jitterStr + stochStr + ".fits"

    # -------------------------------------------------------------------------
    # get NOISE file: process empty trigger file
    # -------------------------------------------------------------------------
    print("\nTESCONSPILEUP: Generating no-events file for NOISE...")
    comm = ("tesconstpileup PixImpList=" + pixFileN + " XMLFile=" + XMLfile
            + " tstop=" + str(simTimeN) + offset +
            " energy=0 pulseDistance=1 TriggerSize=" +
            str(triggerSize) + " clobber=yes")
    args = shlex.split(comm)
    try:
        check_call(args, stderr=STDOUT)
    except:
        print("Error running tool tesconstpileup (NOISE):")
        print(comm)
        os.chdir(cwd)
        shutil.rmtree(tmpDir)
        raise
    print("TESCONSPILEUP: ...................................END")

    print("\nTESSIM: Generating stream file for NOISE...")
    tstart = 0.
    tstop = simTimeN
    comm = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFileN +
            " Streamfile=" + fitsFileN + " sample_rate=" + samplingrate +
            " tstart=" + str(tstart) + " tstop=" + tstop +
            " triggertype=noise " + " triggersize=" + str(triggerSize) +
            " prebuffer=" + str(preBuffer) + " acbias=" + acbias +
            stochastic + " clobber=yes PixType=" + PixTypeFile)

    args = shlex.split(comm)
    try:
        print("Running tool tessim (NOISE):", comm)
        check_call(args, stderr=STDOUT)
    except:
        print("Error running tool tessim (NOISE):")
        print(comm)
        os.chdir(cwd)
        shutil.rmtree(tmpDir)
        raise
    print("TESSIM: ................................END")
    # sys.exit()

    try:
        # rm first record in fits file (noisy)
        print("rm first record in fits file (noisy)")
        auxpy.rmLastAndFirst(fitsFileN, 0)
    except:
        print("Error making up trigger file  (NOISE):")
        os.chdir(cwd)
        shutil.rmtree(tmpDir)
        raise
    
    # Modify EXTNAME if not EXTNAME=RECORDS in empty trigger file
    fN = fits.open(fitsFileN, mode='update')
    fNhdr = fN[1].header
    oldExtName = fN[1].header["EXTNAME"]
    if not oldExtName == "RECORDS":
        print("Changing name...")
        fNhdr['EXTNAME'] = 'RECORDS'
    fN.close()

    if space in ("R", "RFITTED", "RALL", "RNOL"):
        print("RfromItrigger: Adding Resistance column for NOISE...")
        try:
            auxpy.RfromItrigger(infile=fitsFileN, Icol="ADC",
                                current="ADC", Rcol="ADCR",
                                Rmethod=space, samplingrate=samplingrate)
            # delete ADC col
            comm = ("fdelcol infile=" + fitsFileN +
                    "+1 colname=ADC confirm=no proceed=YES")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            # copy ADCR in ADC new col (required by gennoise)
            comm = ("ftcalc infile=" + fitsFileN + " outfile=" + fitsFileN +
                    " clobber=yes column=ADC expr='ADCR'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error adding resistance column in ADCR -> ADC (NOISE):")
            print(comm)
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
            raise
        print("FCALC: ..................................END")

    print("\nGENNOISESPEC: Generating NOISE spectrum file "
          "in (", space, " space)")
    comm = ("gennoisespec --inFile=" + fitsFileN + " --outFile=" + noiseFile +
            " --intervalMinSamples=" + str(pulseLength) +
            " --nintervals=" + str(nintervals) +
            " --scaleFactor=" + str(scaleFactor) +
            " --samplesUp=" + str(samplesUp) +
            " --nSgms=" + str(nSgms) + " --pulse_length=" + str(pulseLength) +
            " --clobber=yes verbosity=0 --weightMS=no ")
    args = shlex.split(comm)
    try:
        check_call(args, stderr=STDOUT)
    except:
        print("Error running tool gennoise (NOISE):")
        print(comm)
        os.chdir(cwd)
        shutil.rmtree(tmpDir)
        raise
    print("GENNOISESPEC: .......................................END")
    print("Noise file ", noiseFile, " has been successfully created\n")
    os.chdir(cwd)
    shutil.rmtree(tmpDir)


#
# --------- Read input parameters and PROCESS simulations -----------------
#
if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
            description='Create NOISE spectrum', prog='simulNoise')

    parser.add_argument('--pixName', required=True,
                        help='Extension name in FITS pixel definition file:\
                        SPA*, LPA1*, LPA2*, LPA3*')
    parser.add_argument('--samprate', default="", choices=['', 'samprate2'],
                        help='baseline, half_baseline')
    parser.add_argument('--jitter',  default="", choices=['', 'jitter'],
                        help='no_jitter, jitter')
    parser.add_argument('--stoch', default="", choices=['', 'stoch'],
                        help="nonstochastic, stochastic")
    parser.add_argument('--pulseLength', required=True, type=int,
                        help='pulse length in samples')
    parser.add_argument('--space', required=True,
                        choices=['ADC', 'R', 'RALL', 'RNOL', 'RFITTED'],
                        help='Data space: ADC for current, R or RALL or\
                        RNOL or RFITTED for resistance')
    parser.add_argument('--acbias', default="yes",
                        help='AC (acbias=yes) or DC (acbias=no)\
                        [default %default]')
    parser.add_argument('--scaleFactor', default=0.0, type=float,
                        help='Param for gennoise[default %default]')
    parser.add_argument('--samplesUp', default=2, type=int,
                        help='Param samplesUp for gennoise [default %default]')
    parser.add_argument('--nSgms', default=5.,
                        help='Param nSgms for gennoise [default %default]')
    parser.add_argument('--simTimeN', default=100,
                        help='Simulation time (s) for noise spectra\
                        calculation [default %default]')
    parser.add_argument('--nintervals', default=1000, type=int,
                        help='Number of intervals in gennoisespec for spectra\
                        calculation [default %default]')
    parser.add_argument('--pixel', default=1,
                        help='Pixel Number [default %default]')

    inargs = parser.parse_args()

    simulNoise(pixName=inargs.pixName, samprate=inargs.samprate,
               jitter=inargs.jitter, stoch=inargs.stoch,
               pulseLength=inargs.pulseLength,
               space=inargs.space, acbias=inargs.acbias,
               scaleFactor=inargs.scaleFactor, samplesUp=inargs.samplesUp,
               nSgms=inargs.nSgms, nintervals=inargs.nintervals,
               simTimeN=inargs.simTimeN, pixel=inargs.pixel)
