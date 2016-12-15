"""
# NOISE spectrum simulation
#
# python simulNoise.py 
#
#  Input parameters: 
#          array (SPA|LPA1|LPA2|LPA3)
#          pulseLength: length of pulses (in samples)
#          space: ADC (current space) or R (resistance space)
#          [baseline]: baseline value
#
#
#
# 1) Simulate 10s stream with no events to calculate Baseline (pixdetillum + tessim)
# 2) Simulate 100s stream with no events as input for gennoisespec (tesconstpileup + tessim + streamtotriggers)
# 3) Obtain noise spectrum with gennoisespec
#
"""

# ----IMPORT MODULES --------------
from __future__ import print_function
import os
import shlex
import shutil
import sys
import tempfile
from subprocess import Popen, check_call, check_output, STDOUT
from astropy.io import fits
from astropy.table import Table, Column

# ----GLOBAL VARIABLES -------------
XMLrect = "/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/CALIBRATION/newFormat/xifu_detector_300us_156kHz_rect.xml"
XMLdir = os.environ["SIXTE"] + "/" + "share/sixte/instruments/athena/1469mm_xifu"
XMLfile = XMLdir + "/" + "xifu_baseline.xml"
XMLfile = XMLdir + "/" + "xifu_detector_hex_baseline.xml"
# XMLfile = XMLrect

cpsN = "pairs"  # counts per second or pairs of pulses
triggerSize = 10000
preBuffer = 1000
tmpDir = tempfile.mkdtemp()
# tmpDir = "/tmp/simulNoise/pfiles"
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]


# ----Error Class definitions ------
# class Error(Exception):
#    """Base class for exceptions in this module."""
#    pass
#
# class SixteError(Error):
#    """Exception raised for errors running SIXTE tools.
#
#    Attributes:
#        tool -- SIXTE tool in which the error occurred
#        msg  -- explanation of the error
#    """
#
#    def __init__(self, tool, msg):
#        self.tool = tool
#        self.msg = msg

def subprocess_cmd(command):
    """ Run concatenated commands in shell (and get output)
    """
    process = Popen(command, shell=True)
    proc_stdout = process.communicate()[0].strip()
    return proc_stdout
    # to also run processes waiting for them to end, use:
    # command = "bla bla bla"
    # args = shlex.split(command)
    # proc = call(args)


def getBaselines(simTimeB, pixel, samprate, PixType, space, pulseLength, array):
    """
    :param simTimeB: simulation time (s) for Baseline
    :param pixel: pixel number
    :param samprate: sampling rate (Hz)
    :param PixType: fits file with pixel definitions
    :param space: "ADC" or "R" or "RB" for noise calculations
    :param pulseLength: pulse length for noise samples (1024?)
    :param array: SPA, LPA1, LPA2, LPA3
    :return: (baseline in current, stddev in current, baseline in R space-if requested)
    """

    global tmpDir
    tmpDir = tempfile.mkdtemp()
    os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
    baselineI = 0.
    baselineR = 0.
    sigmaI = 0.
    sigmaR = 0.
    pixelStr = "%05d" % pixel
    ADCcol = "PXL" + pixelStr
    AMPcol = "PULSE" + pixelStr
    # dataDump = "TIME," + ADCcol
    tessim = "tessim" + array
    root1 = "forBaseline" + str(pulseLength) + "samples_" + tessim + "_" + str(simTimeB) + "s_" + space
    pixFile1 = root1 + ".piximpact"
    streamFile1 = root1 + ".stream"
    # txtFile1 = root1 + ".txt"

    print("\nPIXDETILLUM: Generating file with 0 events for BASELINE...")
    comm = ("pixdetillum PixImpList=" + pixFile1 + " XMLFile=" + XMLfile + " tstart=0. Tstop=" +
            str(simTimeB) + " pixels=" + str(pixel) + " rate=0 energy=0 clobber=yes seed=-1")
    args = shlex.split(comm)
    try:
        check_call(args, stderr=STDOUT)
    except:
        print("Error running tool pixdetillum (BASELINE):")
        print(comm)
        # print(sys.exc_info()[0])
        # os.chdir(cwd)
        shutil.rmtree(tmpDir)
        raise
    print("PIXDETILLUM: .....................................END")

    print("\nTESSIM: Generating stream file for BASELINE...")
    comm = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFile1 + " Streamfile=" + streamFile1 +
            " tstart=0. tstop=" + str(simTimeB) + " sample_rate=" + str(samprate) + " acbias=yes" +
            " triggertype=stream prebuffer=" + str(preBuffer) + " clobber=yes PixType=" + PixType)
    # comm = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFile1 +
    #         " Streamfile=" + streamFile1 + " tstart=0. tstop=" +
    #         str(simTimeB) + " sample_rate=" + str(samprate) +
    #         " clobber=yes PixType=" + PixType)
    args = shlex.split(comm)
    try:
        check_call(args, stderr=STDOUT)
    except:
        print("Error running tool tessim (BASELINE):")
        print(comm)
        shutil.rmtree(tmpDir)
        raise
    print("TESSIM: ...................................END")

    # calculate baseline(s) with a fit to a constant: mean + stddev
    fstr = fits.open(streamFile1)
    tbdata = fstr[1].data
    # -- for ADC: interested in ADCcol
    baselineI = tbdata[ADCcol].mean()
    sigmaI = tbdata[ADCcol].std()
    fstr.close()
    # -- for R/RB: interested in Rcol
    if space in ("R", "RB"):
        print("\nFCALC: Adding Resistance column for BASELINE...")
        Rcol = "AMP2" + space
        RfromIstream(infile=streamFile1, Icol=AMPcol, current="AMP", Rcol=Rcol, Rmethod=space, samprate=samprate)
        # dataDump += "," + Rcol
        fstr = fits.open(streamFile1)
        tbdata = fstr[1].data
        baselineR = tbdata[Rcol].mean()
        sigmaR = tbdata[Rcol].std()
        print("FCALC: .....................................END")

    """ Baseline calculation with R fit to a line
    print("\nFDUMP: Dumping stream file for BASELINE...")
    comm = ("fdump infile=" + streamFile1 + " outfile=" + txtFile1 +
            " prhead=no showrow=no showunit=no clobber=yes showcol=no" +
            " columns='" + dataDump + "' rows='-'")
    args = shlex.split(comm)
    try:
        check_call(args, stderr=STDOUT)
    except:
        print("Error running tool fdump (BASELINE):")
        print(comm)
        raise
    print("FDUMP: ................................END")

    # run R script to get baseline for ADC (data in column 2)
    print("\nGETBASELINE: Calculating BASELINE in Current...")
    Rdir = "/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/NOISE/"
    comm = Rdir + "getBaseline.R " + txtFile1 + " 2"
    print("Comm=", comm)
    args = shlex.split(comm)
    baselineI = float(check_output(args, stderr=STDOUT))
    if baselineI == 0:
        print("Error running R tool getBaseline in ADC:")
        print(comm)
        os.chdir(cwd)
        shutil.rmtree(tmpDir)
        sys.exit()

    if space in ("R", "RB"):
        # run R script to get baseline for AMP2R(B) (data in column 3)
        print("\nGETBASELINE: Calculating BASELINE in Resistance...")
        comm = Rdir + "getBaseline.R " + txtFile1 + " 3"
        print("Comm=", comm)
        args = shlex.split(comm)
        baselineR = float(check_output(args, stderr=STDOUT))
        if baselineR == 0:
            print("Error running R tool getBaseline in AMP2R(B):")
            print(comm)
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
            sys.exit()
    """

    # do some cleaning
    # os.remove(streamFile1)
    # os.remove(txtFile1)
    baselines = (baselineI, sigmaI, baselineR, sigmaR)
    return baselines


def RfromAMP(infile, AMPcol, Rcol):
    """
    :param infile: fits input file with PULSE????? column
    :param AMPcol: existing AMP column name
    :return: infile transformed with a new column AMP2R
    Rcol = R = R0 -R0*((abs(AMPcol -I0)/I0)/(1 + abs(AMPcol-I0)/I0) # Bandler?
    """
    try:
        comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=" + Rcol +
                " expr='#R0 - #R0*((abs(" + AMPcol + "-#I0_START)/#I0_START)" +
                "/(1 + abs(" + AMPcol + "-#I0_START)/#I0_START))'")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except:
        print("Error running ftool for R calculation:")
        print(comm)
        raise
    return


def RBfromAMP(infile, AMPcol, RBcol, samprate):
    """
    :param infile: fits input file with PULSE????? column
    :param AMPcol: existing AMP column name
    :param RBcol: output (new) RB column name
    :param samprate: sampling rate (Hz)
    :return: infile transformed with a new column RBcol
    RBcol = RBIS = (V0 - L*dI/dt)/I (see PP email from 16/11/2015
    V0 = I0 * (R0+RL)
    L = 2*LFILTER/(TTR*TTR)
    """
    fstr = fits.open(infile)
    TRIGGSZ = fstr[1].header["TRIGGSZ"]
    fstr.close()

    try:
        comm = (
        "fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=I expr='#I0_START-" + AMPcol + "'")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DER expr='seqdiff(I)'")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=UNOS expr=1")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fcollen fitsfile=" + infile + " colname=UNOS collen=" + TRIGGSZ)
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DER1COL expr=DER")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fcollen fitsfile=" + infile + " colname=DER1COL collen=1")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fcollen fitsfile=" + infile + " colname=DER1COL collen=" + TRIGGSZ)
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DER0COL1 expr='DER-DER1COL'")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DER2 expr='UNOS*DER[2]'")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DERFINAL expr='DER2+DER0COL1'")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DERIT expr='DERFINAL*" +
                str(samprate) + "'")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fdelcol infile=" + infile + "+1 colname=UNOS confirm=no proceed=yes")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fdelcol infile=" + infile + "+1 colname=DER confirm=no proceed=yes")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fdelcol infile=" + infile + "+1 colname=DER1COL confirm=no proceed=yes")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fdelcol infile=" + infile + "+1 colname=DER0COL1 confirm=no proceed=yes")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = ("fdelcol infile=" + infile + "+1 colname=DER2 confirm=no proceed=yes")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        # R = (V0 - L*dI/dt)/I
        comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=" + RBcol +
                " expr='((#I0_START*(#R0+#RL)) - (2*#LFILTER/(#TTR*#TTR))*DERIT)/I'")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except:
        print("Error running ftool for RB calculation:")
        print(comm)
        raise
    return


def RfromItrigger(infile, Icol, current, Rcol, Rmethod, samprate):
    """
    :param infile: fits input file with current column Icol
    :param Icol: existing current column name
    :param current: input current: ADC or AMP (PULSE0000)
    :param Rcol: new output R column name
    :param Rmethod: resistance method calculation (R or RB)
    :param samprate: sampling rate (Hz)
    :return: infile transformed with a new column in resistance space (Rcol)

    R = R0 -R0*((abs(AMP -I0)/I0)/(1 + abs(AMP-I0)/I0) # Bandler?

    RB = RBIS = (V0 - L*dI/dt)/I (see PP email from 16/11/2015), I=I0-AMP
    V0 = I0 * (R0+RL)
    L = 2*LFILTER/(TTR*TTR)

    """
    if current == "ADC":
        # create col AMP in fits file from ADC column (fcalc does not work)
        try:
            comm = ("ftcalc infile=" + infile + " outfile=" + infile + " clobber=yes column=AMP expr='" + Icol +
                    "*#ADUCNV+#Imin'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftcalc for AMP column calculation:")
            print(comm)
            raise
        Icol = "AMP"

    if Rmethod == "R":
        try:
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=" + Rcol +
                    " expr='#R0 - #R0*((abs(" + Icol + "-#I0_START)/#I0_START)" +
                    "/(1 + abs(" + Icol + "-#I0_START)/#I0_START))'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for R calculation:")
            print(comm)
            raise
    elif Rmethod == "RB":
        fstr = fits.open(infile)
        TRIGGSZ = fstr[1].header["TRIGGSZ"]
        fstr.close()
        try:
            comm = (
            "fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=I expr='#I0_START-" + Icol + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DER expr='seqdiff(I)'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=UNOS expr=1")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcollen infile=" + infile + " colname=UNOS collen=" + str(TRIGGSZ))
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DER1COL expr=DER")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcollen infile=" + infile + " colname=DER1COL collen=1")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcollen infile=" + infile + " colname=DER1COL collen=" + str(TRIGGSZ))
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DER0COL1 expr='DER-DER1COL'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DER2 expr='UNOS*DER[2]'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = (
            "fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DERFINAL expr='DER2+DER0COL1'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DERIT expr='DERFINAL*" +
                    str(samprate) + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fdelcol infile=" + infile + "+1 colname=UNOS confirm=no proceed=yes")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fdelcol infile=" + infile + "+1 colname=DER confirm=no proceed=yes")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fdelcol infile=" + infile + "+1 colname=DER1COL confirm=no proceed=yes")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fdelcol infile=" + infile + "+1 colname=DER0COL1 confirm=no proceed=yes")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fdelcol infile=" + infile + "+1 colname=DER2 confirm=no proceed=yes")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            # R = (V0 - L*dI/dt)/I
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=" + Rcol +
                    " expr='((#I0_START*(#R0+#RL)) - (2*#LFILTER/(#TTR*#TTR))*DERIT)/I'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for RB calculation:")
            print(comm)
            raise
    return


def RfromIstream(infile, Icol, current, Rcol, Rmethod, samprate):
    """
    :param infile: stream fits input file with current column Icol
    :param Icol: existing current column name
    :param current: input current: ADC or AMP (PULSE0000)
    :param Rcol: new output R column name
    :param Rmethod: resistance method calculation (R or RB)
    :param samprate: sampling rate (Hz)
    :return: infile transformed with a new column in resistance space (Rcol)

    R = R0 -R0*((abs(AMP -I0)/I0)/(1 + abs(AMP-I0)/I0) # Bandler?

    RB = RBIS = (V0 - L*dI/dt)/I (see PP email from 16/11/2015), I=I0-AMP
    V0 = I0 * (R0+RL)
    L = 2*LFILTER/(TTR*TTR)

    """
    if current == "ADC":
        # create col AMP in fits file from ADC column (fcalc does not work)
        try:
            comm = ("ftcalc infile=" + infile + " outfile=" + infile + " clobber=yes column=AMP expr='" + Icol +
                    "*#ADUCNV+#Imin'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftcalc for AMP column calculation:")
            print(comm)
            raise
        Icol = "AMP"

    if Rmethod == "R":
        try:
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=" + Rcol +
                    " expr='#R0 - #R0*((abs(" + Icol + "-#I0_START)/#I0_START)" +
                    "/(1 + abs(" + Icol + "-#I0_START)/#I0_START))'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for R calculation:")
            print(comm)
            raise
    elif Rmethod == "RB":
        try:
            comm = (
            "fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=I expr='#I0_START-" + Icol + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DER expr='seqdiff(I)'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DERIT expr='DER*" +
                    str(samprate) + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            # R = (V0 - L*dI/dt)/I
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=" + Rcol +
                    " expr='((#I0_START*(#R0+#RL)) - (2*#LFILTER/(#TTR*#TTR))*DERIT)/I'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for RB calculation:")
            print(comm)
            raise
    return


# ----MAIN routine definition ------

def simulNoise(array, pulseLength, space, scaleFactor, samplesUp, nSgms, simTimeB, simTimeN, samprate, pixel):
    """ simulate data in input parameter space, calculate the data baseline and 
          create Noise file to be ingested in SIRENA processing tasks
    """
    cwd = os.getcwd()

    # set input params dependent variables
    # -------------------------------------

    # PixType = "'file:/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/tespixels.fits[" + array + "]'"
    PixType = "'file:/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/newpixels.fits[" + array + "ac]'"
    tessim = "tessim" + array
    wdir = "NOISE/" + tessim
    os.chdir(wdir)
    # define files
    root2 = "forNoise" + str(pulseLength) + "samples_" + tessim + "_" + str(simTimeN) + "s_" + cpsN + "cps_" + space
    noiseFile = "noise" + str(pulseLength) + "samples_" + tessim + "_B0_" + str(simTimeN) + "s_" + cpsN + "cps_" + \
                space + ".fits"
    pixFile2 = root2 + ".piximpact"
    streamFile2 = root2 + ".stream"
    fitsFile2 = root2 + ".fits"

    # -------------------------------------------------------------------------
    # get baseline
    # -------------------------------------------------------------------------
    (baselineI, sigmaI, baselineR, sigmaR) = getBaselines(simTimeB, pixel, samprate, PixType, space, pulseLength, array)
    # baselineI = 567.36
    # baselineR = 0.0005497636

    print("########################################\n")
    print("    BaselineI/R=", baselineI, baselineR)
    print("    SigmaI/R=", sigmaI,sigmaR)
    print("\n######################################\n")
    # sys.exit()

    # -------------------------------------------------------------------------
    # get NOISE file: process empty trigger file
    # -------------------------------------------------------------------------
    print("\nTESCONSPILEUP: Generating no-events file for NOISE...")
    comm = ("tesconstpileup PixImpList=" + pixFile2 + " XMLFile=" + XMLfile + " tstop=" + str(simTimeN) +
            " energy=0 pulseDistance=1 TriggerSize=" + str(triggerSize) + " clobber=yes")
    # print("\nPIXDETILLUM: Generating no events file for NOISE...not usable if required streamtotriggers")
    # comm = ("pixdetillum PixImpList=" + pixFile2 + " XMLFile=" + XMLrect + " tstart=0. Tstop=" +
    #         str(simTimeN) + " pixels=" + str(pixel) + " rate=0 energy=0 clobber=yes seed=-1")
    args = shlex.split(comm)
    try:
        check_call(args, stderr=STDOUT)
    except:
        print("Error running tool tesconstpileup (NOISE):")
        print(comm)
        # print(sys.exc_info()[0])
        os.chdir(cwd)
        shutil.rmtree(tmpDir)
        raise
    print("TESCONSPILEUP: ...................................END")

    print("\nTESSIM: Generating stream file for NOISE...")
    simTimeNleft = simTimeN
    tstart = 0.
    tstop = min(30., simTimeNleft)
    tstop = simTimeN
    i = 0
    streamFiles = ""
    while simTimeNleft > 0:
        simTimeNleft -= 30.
        i += 1
        streamFile = "noise.stream" + str(i)
        #comm = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFile2 + " Streamfile=" + streamFile +
                " tstart=" + str(tstart) + " tstop=" + str(tstop) + " sample_rate=" + str(samprate) +
                " acbias=yes triggertype=stream " + " triggersize=" + str(triggerSize) + " prebuffer=" +
                str(preBuffer) + " clobber=yes PixType=" + PixType)
        comm = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFile2 + " Streamfile=" + streamFile +
                " tstart=" + str(tstart) + " tstop=" + str(tstop) + " sample_rate=" + str(samprate) +
                " acbias=yes triggertype=stream " + " triggersize=" + str(triggerSize) + " prebuffer=" +
                str(preBuffer) + " clobber=yes PixType=" + PixType)
        args = shlex.split(comm)
        try:
            check_call(args, stderr=STDOUT)
        except:
            print("Error running tool tessim (NOISE):")
            print(comm)
            # print(sys.exc_info()[0])
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
            raise
        tstart = tstop
        tstop = tstop + min(30, simTimeNleft)
        streamFiles = streamFiles + " " + streamFile
    comm = ("fmerge infiles='" + streamFiles + "' outfile=" + streamFile2 + " columns=- clobber=yes")
    args = shlex.split(comm)
    try:
        check_call(args, stderr=STDOUT)
    except:
        print("Error running tool fmerge (NOISE):")
        print(comm)
        os.chdir(cwd)
        shutil.rmtree(tmpDir)
        raise
    print("TESSIM: ................................END")

    print("\nSTREAMTOTRIGGERS: Generating trigger file for NOISE...")
    comm = ("streamtotriggers PixImpList=" + pixFile2 + " XMLFile=" +
            XMLfile + " tstart=0. tstop=" + str(simTimeN) + " Streamfile=" +
            streamFile2 + " TesTriggerFile=" + fitsFile2 + " TriggerSize=10000 " +
            "PreBufferSize=1000 pixels=" + str(pixel) + " clobber=yes")
    args = shlex.split(comm)
    try:
        check_call(args, stderr=STDOUT)
    except:
        print("Error running tool streamtotriggers (NOISE):")
        print(comm)
        # print(sys.exc_info()[0])
        os.chdir(cwd)
        shutil.rmtree(tmpDir)
        raise
    print("STREAMTOTRIGGERS: .................................END")
    # sys.exit()

    print("\nMAKINGUP (update header + remove tricky records + fixed-lenghtify) trigger file for NOISE...")
    try:
        # save fitsFile2 header (modified by stilts)
        shutil.copyfile(fitsFile2, 'pp.fits')
        # astropy & ftools & stilts does not work with variable-length arrays, so a previous
        # conversion to fixed-format is required)
        comm = "stilts tpipe cmd='addcol TIME2 \"(TIME*2)\"' cmd='delcols TIME2' in=" + fitsFile2 + " out=" + fitsFile2
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        # Add stream & old-fits header to new-fits file header
        comm = "cphead " + streamFile2 + "+1 " + fitsFile2 + "+1"
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = "cphead pp.fits+1 " + fitsFile2 + "+1"
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        # rm first record in fits file (noisy)
        comm = ("fdelrow infile=" + fitsFile2 + "+1 firstrow=1 nrows=1 confirm=no proceed=yes")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except:
        print("Error making up trigger file  (NOISE):")
        print(comm)
        # print(sys.exc_info()[0])
        os.chdir(cwd)
        shutil.rmtree(tmpDir)
        raise
    os.remove('pp.fits')
    # os.remove(streamFile2)
    print("MAKINGUP: ...............................................................END")
    # sys.exit()

    # Modify EXTNAME if not EXTNAME=RECORDS in empty trigger file
    f2 = fits.open(fitsFile2, mode='update')
    f2hdr = f2[1].header
    oldExtName = f2[1].header["EXTNAME"]
    if not oldExtName == "RECORDS":
        print("Changing name...")
        f2hdr['EXTNAME'] = 'RECORDS'
    f2.close()

    if space in ("R", "RB"):

        # calculate AMP & ADCR columns with astropy.fits (does not work with varlenght cols)
        # f2tb    = f2[1].data
        # f2cols  = f2[1].columns
        # ADUCNV  = f2[1].header["ADUCNV"]
        # Imin    = f2[1].header["Imin"]
        # Ibias   = f2[1].header["I0_START"]
        # R0      = f2[1].header["R0"]
        # colAMP  = fits.Column(name='AMP', format=f2cols[f2cols.names.index('ADC')].format,
        #                       unit='Amp', array=(f2tb['ADC']*ADUCNV+Imin))
        # f2cols.add_col(colAMP)
        # colADCR = fits.Column(name='ADCR', format=f2cols[f2cols.names.index('ADC')].format,
        #                       unit='Ohm', array=(R0 - R0*((abs(f2tb['ADC']*ADUCNV+Imin-Ibias)/Ibias) /
        #                                       (1 + abs(f2tb['ADC']*ADUCNV+Imin-Ibias)/Ibias))))
        # f2cols.add_col(colADCR)
        # f2tbnew = fits.BinTableHDU.from_columns(columns=f2cols, header=f2hdr)
        # f2tbnew.writeto(fitsFile2, clobber=True)
        # f2.close()
        print("RfromItrigger: Adding Resistance column for NOISE...")
        try:
            # # create col AMP in fits file (fcalc does not work ok)
            # comm = ("ftcalc infile=" + fitsFile2 + " outfile=" + fitsFile2 +
            #         " clobber=yes column=AMP expr='ADC*#ADUCNV+#Imin'")
            # args = shlex.split(comm)
            # print(comm)
            # check_call(args, stderr=STDOUT)
            # # sys.exit()
            # # convert ADC to R --> ADCR
            # comm = ("ftcalc infile=" + fitsFile2 + " outfile=" +
            #         fitsFile2 + " clobber=yes column=ADCR " +
            #         "expr='#R0 - #R0*((abs(AMP-#I0_START)/#I0_START)" +
            #         "/(1 + abs(AMP-#I0_START)/#I0_START))'")
            # args = shlex.split(comm)
            # check_call(args, stderr=STDOUT)
            RfromItrigger(infile=fitsFile2, Icol="ADC", current="ADC", Rcol="ADCR", Rmethod=space, samprate=samprate)
            # delete ADC col
            comm = ("fdelcol infile=" + fitsFile2 + "+1 colname=ADC confirm=no proceed=YES")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            # copy ADCR in ADC new col (required by gennoise)
            comm = ("ftcalc infile=" + fitsFile2 + " outfile=" + fitsFile2 +
                    " clobber=yes column=ADC expr='ADCR'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)

        except:
            print("Error adding resistance column in ADCR -> ADC (NOISE):")
            print(comm)
            # print(sys.exc_info()[0])
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
            raise
        print("FCALC: ..................................END")

    # sys.exit()
    baseline = baselineI
    if space in ("R", "RB"):
        baseline = baselineR

    print("\nGENNOISESPEC: Generating NOISE spectrum file in (", space, " space)")
    comm = ("gennoisespec --inFile=" + fitsFile2 + " --outFile=" + noiseFile
            + " --intervalMinSamples=" + str(pulseLength) +
            " --nintervals=1000 " + " --scaleFactor=" + str(scaleFactor) +
            " --samplesUp=" + str(samplesUp) + " --nSgms=" + str(nSgms) +
            " --pulse_length=" + str(pulseLength) + " --baseline=" +
            str(baseline) + " --ntausPF=0 --clobber=yes verbosity=0 ")
    args = shlex.split(comm)
    try:
        check_call(args, stderr=STDOUT)
    except:
        print("Error running tool gennoise (NOISE):")
        print(comm)
        # print(sys.exc_info()[0])
        os.chdir(cwd)
        shutil.rmtree(tmpDir)
        raise
    print("GENNOISESPEC: .......................................END")
    # sys.exit()

    # add keywords BASELINI & BASELINR to noise file header
    hdulist = fits.open(noiseFile, mode='update')
    hdulist[1].header["BASELINE"] = baselineI
    hdulist[1].header["BASELINR"] = baselineR
    hdulist[1].header["NOISESTD"] = sigmaI
    hdulist[1].header["NOISESTR"] = sigmaR
    hdulist.close()
    print("Noise file ", noiseFile, " has been successfully created\n")
    os.chdir(cwd)
    shutil.rmtree(tmpDir)


#
# --------- Read input parameters and PROCESS simulations -----------------
#
if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option('--array', '-a', action='store', dest='array', type='string',
                      help='Array acronym: SPA, LPA1, LPA2, LPA3')
    parser.add_option('--pulseLength', '-p', action='store', dest='pulseLength', type='int',
                      help='pulse length in samples')
    parser.add_option('--space', '-s', action='store', dest='space', type='string',
                      help='Data space: ADC for current, R or RB for resistance')
    parser.add_option('--scaleFactor', '-f', action='store', dest='scaleFactor', type='float',
                      help='Param scaleFactor for gennoise [default %default]', default=0.005)
    parser.add_option('--samplesUp', '-U', action='store', dest='samplesUp', type='float',
                      help='Param samplesUp for gennoise [default %default]', default=2)
    parser.add_option('--nSgms', '-g', action='store', dest='nSgms', type='float',
                      help='Param nSgms for gennoise [default %default]', default=20)
    parser.add_option('--samprate', '-r', action='store', dest='samprate', type='float',
                      help='Param samprate (Hz) for tessim [default %default]', default=156250)
    parser.add_option('--simTimeB', '-B', action='store', dest='simTimeB', type='float',
                      help='Simulation time (s) for baseline calculation [default %default]', default=10)
    parser.add_option('--simTimeN', '-N', action='store', dest='simTimeN', type='float',
                      help='Simulation time (s) for noise spectra calculation [default %default]', default=100)
    parser.add_option('--pixel', '-x', action='store', dest='pixel', type='int', help='Pixel [default %default]',
                      default=1)

    (opts, args) = parser.parse_args()

    if not opts.array:

        message = "ERROR:  Data array must be provided (SPA|LPA1|LPA2|LPA3)."
        sys.exit(message)

    elif not opts.pulseLength:

        message = "ERROR:  Pulse length (in samples) must be provided."
        sys.exit(message)

    elif not opts.space:

        message = "ERROR:  Data space (ADC or R or RB) must be provided."
        sys.exit(message)

    else:

        simulNoise(array=opts.array, pulseLength=opts.pulseLength,
                   space=opts.space, scaleFactor=opts.scaleFactor,
                   samplesUp=opts.samplesUp, nSgms=opts.nSgms,
                   simTimeB=opts.simTimeB, simTimeN=opts.simTimeN,
                   samprate=opts.samprate, pixel=opts.pixel)
