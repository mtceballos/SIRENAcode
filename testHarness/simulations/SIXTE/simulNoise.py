"""
# NOISE spectrum simulation
#
# python simulNoise.py 
#
#  Input parameters: 
#          pixName (SPA*|LPA1*|LPA2*|LPA3*)
#          pulseLength: length of pulses (in samples)
*          nsamples: samples for the noise
#          space: ADC (current space) or R or RALL or RNOL or RFITTED(resistance space)
#
#
#
# 1) Simulate 10s stream with no events to calculate Baseline (pixdetillum + tessim) --> not requiered anymore
# 2) Simulate 100s stream with no events as input for gennoisespec (tesconstpileup + tessim)
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
import xml.etree.ElementTree as ET
from subprocess import Popen, check_call, STDOUT
from astropy.io import fits

# ----GLOBAL VARIABLES -------------
# XMLrect = "$HOME/INSTRUMEN/EURECA/testHarness/simulations/SIXTE/CALIBRATION/newFormat/xifu_detector_300us_156kHz_rect.xml"
XMLdir = os.environ["SIXTE"] + "/" + "share/sixte/instruments/athena/1469mm_xifu"
# XMLfile = XMLdir + "/" + "xifu_baseline.xml"
XMLfile = XMLdir + "/" + "xifu_detector_hex_baseline.xml"
# XMLfile = XMLrect

XMLtree = ET.parse(XMLfile)
XMLroot = XMLtree.getroot()
for samplefreq in XMLroot.findall('samplefreq'):
    samprate = samplefreq.get('value')

cpsN = "pairs"  # counts per second or pairs of pulses
preBuffer = 1000
Ifit = 45.3E-6
tmpDir = tempfile.mkdtemp()
tmpFile = tempfile.TemporaryFile()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]


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


# Baseline calculation is no longer necessary as it is done by default in gennoisespec...
def getBaselines(simTimeB, pixel, PixTypeFile, space, pulseLength, pixName, acbias):
    """
    :param simTimeB: simulation time (s) for Baseline
    :param pixel: pixel number
    :param PixTypeFile: fits file with pixel definitions
    :param space: "ADC" or "R" or "RALL" or "RNOL" for noise calculations
    :param pulseLength: pulse length for noise samples (1024?)
    :param pixName: SPA*, LPA1*, LPA2*, LPA3*
    :param acbias: AC (acbias=yes) or DC (acbias=no)
    :return: (baseline in current, stddev in current, baseline in R space-if requested)
    """

    global tmpDir, samprate
    baselineI = 0.
    baselineR = 0.
    sigmaI = 0.
    sigmaR = 0.
    pixelStr = "%05d" % pixel
    ADCcol = "PXL" + pixelStr
    AMPcol = "PULSE" + pixelStr
    tessim = "tessim" + pixName
    rootB = "forBaseline" + str(pulseLength) + "samples_" + tessim + "_" + str(simTimeB) + "s_" + space
    pixFileB = rootB + ".piximpact"
    streamFileB = rootB + ".stream"

    # generate PIXIMPACT for noise stream (pixdetillum  or tesconstpileup can be used)
    print("\nPIXDETILLUM: Generating file with 0 events for BASELINE...")
    comm = ("pixdetillum PixImpList=" + pixFileB + " XMLFile=" + XMLfile + " tstart=0. Tstop=" +
            str(simTimeB) + " pixels=" + str(pixel) + " rate=0 energy=0 clobber=yes seed=-1")
    args = shlex.split(comm)
    try:
        check_call(args, stderr=STDOUT)
    except:
        print("Error running tool pixdetillum (BASELINE):")
        print(comm)
        shutil.rmtree(tmpDir)
        raise
    print("PIXDETILLUM: .....................................END")

    print("\nTESSIM: Generating stream file for BASELINE...")
    comm = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFileB + " Streamfile=" + streamFileB +
            " tstart=0. tstop=" + str(simTimeB) + " acbias=" + acbias + " sample_rate=" + samprate +
            " triggertype=stream prebuffer=" + str(preBuffer) + " clobber=yes PixType=" + PixTypeFile)
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
    fstrB = fits.open(streamFileB)
    tbdataB = fstrB[1].data
    # -- for ADC: interested in ADCcol
    baselineI = tbdataB[ADCcol].mean()
    sigmaI = tbdataB[ADCcol].std()
    fstrB.close()
    # -- for R/RB: interested in Rcol
    if "R" in space:
        print("\nFCALC: Adding Resistance column for BASELINE...")
        Rcol = "AMP2" + space
        RfromIstream(infile=streamFileB, Icol=AMPcol, current="AMP", Rcol=Rcol, Rmethod=space)
        fstrB = fits.open(streamFileB)
        tbdataB = fstrB[1].data
        baselineR = tbdataB[Rcol].mean()
        sigmaR = tbdataB[Rcol].std()
        fstrB.close()
        print("FCALC: .....................................END")

    baselines = (baselineI, sigmaI, baselineR, sigmaR)
    return baselines

# Different R transformation are requiered due to the different format in stream/trigger files (derivatives)

def RfromItrigger(infile, Icol, current, Rcol, Rmethod):
    """
    :param infile: fits input file with current column Icol
    :param Icol: existing current column name
    :param current: input current: ADC or AMP (PULSE0000)
    :param Rcol: new output R column name
    :param Rmethod: resistance method calculation (R or RALL or RNOL or RFITTED)
    :return: infile transformed with a new column in resistance space (Rcol)

    R = R0 -R0*((abs(AMP -I0)/I0)/(1 + abs(AMP-I0)/I0) # Bandler?

    RALL = (V0 - I*RL -L*dI/dt)/I (see PP SPIE paper)
    RNOL = (V0 - I*RL)/I (see PP SPIE paper)
    RFITTED = V0/(Ifit + I) (see PP SPIE paper)
    I=I0-AMP
    V0 = I0 * (R0+RL)
    L = LFILTER/(TTR*TTR)
    RL = RPARA/(TTR*TTR)
    Ifit = 45.3E-6 AMP

    """
    global Ifit
    comm = ""
    fstr = fits.open(infile)
    triggsz = fstr[1].header['TRIGGSZ']
    RPARA = fstr[1].header['RPARA']
    TTR = fstr[1].header['TTR']
    LFILTER = fstr[1].header['LFILTER']
    RL = RPARA/(TTR*TTR)
    L = LFILTER/(TTR*TTR)
    fstr.close()

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
    elif Rmethod == "RALL":
        try:
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=I expr='#I0_START-" +
                    Icol + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DER expr='seqdiff(I)'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=UNOS expr=1")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcollen infile=" + infile + " colname=UNOS collen=" + str(triggsz))
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DER1COL expr=DER")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcollen infile=" + infile + " colname=DER1COL collen=1")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcollen infile=" + infile + " colname=DER1COL collen=" + str(triggsz))
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DER0COL1 expr='DER-DER1COL'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DER2 expr='UNOS*DER[2]'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=DERFINAL expr='DER2+DER0COL1'")
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
            # R = (V0 - I*RL - L*dI/dt)/I
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=" + Rcol +
                    " expr='((#I0_START*(#R0+" + str(RL) + ")) - (I*" + str(RL) + ") - (" + str(L) + "*DERIT))/I'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for RALL calculation:")
            print(comm)
            raise
    elif Rmethod == "RNOL":

        try:
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=I expr='#I0_START-" +
                    Icol + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)

            # R = (V0 - I*RL)/I
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=" + Rcol +
                    " expr='((#I0_START*(#R0+" + str(RL) + ")) - (I*" + str(RL) + "))/I'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for RNOL calculation:")
            print(comm)
            raise
    elif Rmethod == "RFITTED":

        try:
            # IIF = I + Ifit = (I0 - AMP) + Ifit
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=IIF expr='#I0_START-" +
                    Icol + "+" + str(Ifit) + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)

            # R = V0/IIF
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=" + Rcol +
                    " expr='#I0_START*(#R0+" + str(RL) + ")/IIF'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for RFITTED calculation:")
            print(comm)
            raise

    return


def RfromIstream(infile, Icol, current, Rcol, Rmethod):
    """
    :param infile: stream fits input file with current column Icol
    :param Icol: existing current column name
    :param current: input current: ADC or AMP (PULSE0000)
    :param Rcol: new output R column name
    :param Rmethod: resistance method calculation (R or RALL or RNOL)
    :return: infile transformed with a new column in resistance space (Rcol)

    R = R0 -R0*((abs(AMP -I0)/I0)/(1 + abs(AMP-I0)/I0) # Bandler?

    RALL = (V0 - I*RL -L*dI/dt)/I (see PP SPIE paper)
    RNOL = (V0 - I*RL)/I (see PP SPIE paper)
    RFITTED = V0/(Ifit + I) (see PP SPIE paper)
    I=I0-AMP
    V0 = I0 * (R0+RL)
    L = LFILTER/(TTR*TTR)
    RL = RPARA/(TTR*TTR)
    Ifit = 45.3E-6 AMP

    """
    global samprate, Ifit
    comm = ""

    triggsz = fstr[1].header['TRIGGSZ']
    RPARA = fstr[1].header['RPARA']
    TTR = fstr[1].header['TTR']
    LFILTER = fstr[1].header['LFILTER']
    RL = RPARA / (TTR * TTR)
    L = LFILTER / (TTR * TTR)

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
    elif Rmethod == "RALL":
        try:
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=I expr='#I0_START-" +
                    Icol + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DER expr='seqdiff(I)'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=DERIT expr='DER*" +
                    str(samprate) + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            # R = (V0 - I*RL - L*dI/dt)/I
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=" + Rcol +
                    " expr='((#I0_START*(#R0+" + str(RL) + ")) - (I*" + str(RL) + ") - (" + str(L) + "*DERIT))/I'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for RALL calculation:")
            print(comm)
            raise
    elif Rmethod == "RNOL":
        try:
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=I expr='#I0_START-" +
                    Icol + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            # R = (V0 - I*RL)/I
            comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=" + Rcol +
                    " expr='((#I0_START*(#R0+" + str(RL) + ")) - (I*" + str(RL) + "))/I'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for RNOL calculation:")
            print(comm)
            raise
    elif Rmethod == "RFITTED":
            try:
                # IIF = I + Ifit = (I0 - AMP) + Ifit
                comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=IIF expr='#I0_START-" +
                        Icol + "+" + str(Ifit) + "'")
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)

                # R = V0/IIF
                comm = ("fcalc infile=" + infile + " outfile=" + infile + " clobber=yes clname=" + Rcol +
                        " expr='#I0_START*(#R0+" + str(RL) + ")/IIF'")
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except:
                print("Error running ftool for RFITTED calculation:")
                print(comm)
                raise

    return


# ----MAIN routine definition ------

def simulNoise(pixName, pulseLength, space, acbias, scaleFactor, samplesUp, nSgms, simTimeN, pixel):
    """ simulate data in input parameter space, calculate the data baseline and 
          create Noise file to be ingested in SIRENA processing tasks

          :param pixName : name of extension in pixel definitin file SPA*, LPA1*, LPA2*, LPA3*
          :param pulseLength : pulse length in samples to select sampling in noise spectrum
          :param space :  Data space: ADC for current, R or RALL or RNOL or RFITTED for resistance
          :param acbias : AC (acbias=yes) or DC (acbias=no)
          :param scaleFactor : Param scaleFactor for gennoise tool
          :param samplesUp : Param samplesUp for gennoise tool
          :param nSgms : Param nSgms for gennoise tool
          :param simTimeN : Simulation time (s) for noise spectra calculation
          :param pixel : Pixel number
          :return fits file with Noise spectra in specified space and sampling

    """

    global samprate, preBuffer

    cwd = os.getcwd()

    # set input params dependent variables
    # -------------------------------------
    triggerSize = max(10000, pulseLength+preBuffer)
    simSIXTEdir = "/home/ceballos/INSTRUMEN/EURECA/testHarness/simulations/SIXTE"
    # PixTypeFile = "'file:" + simSIXTEdir + "/tespixels.fits[" + array + "]'"
    PixTypeFile = "'file:" + simSIXTEdir + "/newpixels.fits[" + pixName + "]'"
    tessim = "tessim" + pixName
    wdir = "NOISE/" + tessim
    os.chdir(wdir)
    # define files
    rootN = "forNoise" + str(pulseLength) + "samples_" + tessim + "_" + str(simTimeN) + "s_" + cpsN + "cps_" + space
    #noiseFile = "noise" + str(pulseLength) + "samples_" + tessim + "_B0_" + str(simTimeN) + "s_" + cpsN + "cps_" +
    # space + ".fits"
    noiseFile = "noise" + str(pulseLength) + "samples_" + tessim + "_B0_" + space + ".fits"
    pixFileN = rootN + ".piximpact"
    fitsFileN = rootN + ".fits"

    # -------------------------------------------------------------------------
    # get baseline
    # -------------------------------------------------------------------------
    # (baselineI, sigmaI, baselineR, sigmaR) = \
    #    getBaselines(simTimeB, pixel, PixTypeFile, space, pulseLength, pixName, acbias)
    # # baselineI = 567.36
    # # baselineR = 0.0005497636

    # print("########################################\n")
    # print("    BaselineI/R=", baselineI, baselineR)
    # print("    SigmaI/R=", sigmaI,sigmaR)
    # print("\n######################################\n")
    # sys.exit()

    # -------------------------------------------------------------------------
    # get NOISE file: process empty trigger file
    # -------------------------------------------------------------------------
    print("\nTESCONSPILEUP: Generating no-events file for NOISE...")
    # Also possible with pixdetillum if no streamtotriggers must be done...
    comm = ("tesconstpileup PixImpList=" + pixFileN + " XMLFile=" + XMLfile + " tstop=" + str(simTimeN) +
            " energy=0 pulseDistance=1 TriggerSize=" + str(triggerSize) + " clobber=yes")
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
    comm = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFileN + " Streamfile=" + fitsFileN + " sample_rate=" +
            samprate + " tstart=" + str(tstart) + " tstop=" + str(tstop) + " triggertype=noise " + " triggersize=" +
            str(triggerSize) + " prebuffer=" + str(preBuffer) + " acbias=" + acbias + " clobber=yes PixType=" +
            PixTypeFile)
    args = shlex.split(comm)
    try:
        check_call(args, stderr=STDOUT)
    except:
        print("Error running tool tessim (NOISE):")
        print(comm)
        os.chdir(cwd)
        shutil.rmtree(tmpDir)
        raise
    print("TESSIM: ................................END")

    print("\nMAKINGUP (remove tricky first records + fixed-lenghtify) trigger file for NOISE...")
    try:
        # save fitsFileN header (will be modified by stilts)
        shutil.copyfile(fitsFileN, tmpFile.name)
        # astropy & ftools & stilts does not work with variable-length arrays, so a previous
        # conversion to fixed-format is required)
        comm = "stilts tpipe cmd='addcol TIME2 \"(TIME*2)\"' cmd='delcols TIME2' in=" + fitsFileN + " out=" + fitsFileN
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        comm = "cphead " + tmpFile.name + "+1 " + fitsFileN + "+1"
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        # rm first record in fits file (noisy)
        comm = ("fdelrow infile=" + fitsFileN + "+1 firstrow=1 nrows=1 confirm=no proceed=yes")
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except:
        print("Error making up trigger file  (NOISE):")
        print(comm)
        os.chdir(cwd)
        shutil.rmtree(tmpDir)
        raise
    print("MAKINGUP: ...............................................................END")
    # sys.exit()

    # Modify EXTNAME if not EXTNAME=RECORDS in empty trigger file
    # And ADD MONOEN=0 keyword
    fN = fits.open(fitsFileN, mode='update')
    fNhdr = fN[1].header
    fNhdr['MONOEN'] = 0.
    oldExtName = fN[1].header["EXTNAME"]
    if not oldExtName == "RECORDS":
        print("Changing name...")
        fNhdr['EXTNAME'] = 'RECORDS'
    fN.close()

    if space in ("R", "RFITTED", "RALL", "RNOL"):
        print("RfromItrigger: Adding Resistance column for NOISE...")
        try:
            RfromItrigger(infile=fitsFileN, Icol="ADC", current="ADC", Rcol="ADCR", Rmethod=space)
            # delete ADC col
            comm = ("fdelcol infile=" + fitsFileN + "+1 colname=ADC confirm=no proceed=YES")
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

    # sys.exit()
    # baseline = baselineI
    # if space in ("R", "RALL", "RNOL"):
    #    baseline = baselineR

    print("\nGENNOISESPEC: Generating NOISE spectrum file in (", space, " space)")
    comm = ("gennoisespec --inFile=" + fitsFileN + " --outFile=" + noiseFile +
            " --intervalMinSamples=" + str(pulseLength) +
            " --nintervals=1000 " + " --scaleFactor=" + str(scaleFactor) +
            " --samplesUp=" + str(samplesUp) + " --nSgms=" + str(nSgms) +
            " --pulse_length=" + str(pulseLength) + " --clobber=yes verbosity=0 ")
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
    # sys.exit()

    # add keywords BASELINI & BASELINR to noise file header
    # hdulist = fits.open(noiseFile, mode='update')
    # hdulist[1].header["BASELINE"] = baselineI
    # hdulist[1].header["BASELINR"] = baselineR
    # hdulist[1].header["NOISESTD"] = sigmaI
    # hdulist[1].header["NOISESTR"] = sigmaR
    # hdulist.close()
    print("Noise file ", noiseFile, " has been successfully created\n")
    os.chdir(cwd)
    shutil.rmtree(tmpDir)


#
# --------- Read input parameters and PROCESS simulations -----------------
#
if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()

    parser.add_option('--pixName', action='store', dest='pixName', type='string',
                      help='Extension name in FITS pixel definition file: SPA*, LPA1*, LPA2*, LPA3*')
    parser.add_option('--pulseLength', action='store', dest='pulseLength', type='int',
                      help='pulse length in samples')
    parser.add_option('--space', action='store', dest='space', type='string',
                      help='Data space: ADC for current, R or RALL or RNOL or RFITTED for resistance')
    parser.add_option('--acbias', action='store', dest='acbias', type='string',
                      help='AC (acbias=yes) or DC (acbias=no)  [default %default]', default='yes')
    parser.add_option('--scaleFactor', action='store', dest='scaleFactor', type='float',
                      help='Param scaleFactor for gennoise [default %default]', default=0.005)
    parser.add_option('--samplesUp', action='store', dest='samplesUp', type='float',
                      help='Param samplesUp for gennoise [default %default]', default=2)
    parser.add_option('--nSgms', action='store', dest='nSgms', type='float',
                      help='Param nSgms for gennoise [default %default]', default=20)
    parser.add_option('--simTimeN',  action='store', dest='simTimeN', type='float',
                      help='Simulation time (s) for noise spectra calculation [default %default]', default=100)
    parser.add_option('--pixel', action='store', dest='pixel', type='int', help='Pixel [default %default]',
                      default=1)

    (opts, args) = parser.parse_args()

    if not opts.pixName:

        message = "ERROR:  pixName must be provided (SPA*|LPA1*|LPA2*|LPA3*)."
        sys.exit(message)

    elif not opts.pulseLength:

        message = "ERROR:  Pulse length (in samples) must be provided."
        sys.exit(message)

    elif not opts.space:

        message = "ERROR:  Data space (ADC or R or RALL or RNOL or RFITTED) must be provided."
        sys.exit(message)

    elif not opts.acbias:

        message = "ERROR:  AC (acbias=yes) or DC (acbias=no) must be provided."
        sys.exit(message)
    else:

        simulNoise(pixName=opts.pixName, pulseLength=opts.pulseLength,
                   space=opts.space, acbias=opts.acbias,
                   scaleFactor=opts.scaleFactor,
                   samplesUp=opts.samplesUp, nSgms=opts.nSgms,
                   simTimeN=opts.simTimeN, pixel=opts.pixel)
