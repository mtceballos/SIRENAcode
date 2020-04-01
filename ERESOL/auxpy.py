from __future__ import print_function
import os
# from sys import exit
import shlex
from shutil import copy, rmtree
import tempfile
from time import gmtime, strftime
import numpy as np
import numpy.polynomial.polynomial as poly
import math
from subprocess import check_call, STDOUT
from astropy.io import fits, ascii
# from astropy.table import Table
import json
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import seaborn as sns
from scipy.signal import argrelextrema
import scipy.integrate as integrate



# import xml.etree.ElementTree as ET

tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"
SIXTEinst = os.environ["SIXTE"] + "/share/sixte/instruments/athena-xifu/"
XIFUSIMinst = os.environ["SIXTE"] + "/share/xifusim/instruments/"
XMLsixte = (SIXTEinst +
            "/xifu_detector_lpa_75um_AR0.5_pixoffset_mux40_pitch275um.xml")
# samplerate-dependent quantities
sampfreqs = (156250., 78125, 39062.5)  # Hz
sampids = ("", "samprate2", "samprate4")
sampStrs = ("", "_samprate2", "_samprate4")
separations = ('40000', '20000', '10000')
samplesUps = (3, 2, 2)
samplesDowns = (4, 3, 3)
nSigmss = (3.5, 4.5, 4)
preBufferPulses = 1000
scaleFactor = 0


def concatrows(inFile, ext, outFile, ncon):
    """
    Concatenate nrows contiguous rows in a fitsfile to get larger records

    :type inFile: str
    :param inFile: 1-col input FITS file name (with short records)

    :type outFile: str
    :param outFile: output FITS file name (with larger records)

    :type ext: int
    :param ext: extension number

    :type ncon: int
    :param ncon: number of rows to be concatenated

    """

    f = fits.open(inFile)
    data = f[ext].data
    colname = f[ext].header["TTYPE1"]
    coldim0 = f[ext].header["TFORM1"]
    coldim0 = coldim0.rstrip('E')
    nrows0 = f[ext].header["NAXIS2"]
    coldata = data[colname]
    f.close()

    nrows1 = nrows0//ncon
    coldim1 = int(coldim0) * ncon
    data2 = np.zeros((nrows1, coldim1))
    coldim1 = str(coldim1) + 'E'

    j = 0
    for i in range(0, nrows0, ncon):
        if j > nrows1-1:
            break
        # data2[j] = np.concatenate((coldata[i], coldata[i+1], coldata[i+2]))
        pp = coldata[i]
        for k in range(1, ncon):
            pp = np.concatenate((pp, coldata[i+k]))
        data2[j] = pp
        j = j + 1

    newcol1 = fits.Column(name='ADC', format=coldim1, array=data2)
    t = fits.BinTableHDU.from_columns([newcol1, ])
    t.writeto(outFile)


def addcolumn(inFile, ext, colname, datacol):
    """
    Add column to existing file
    :type inFile: str
    :param inFile: input FITS file name

    :type colname: str
    :param colname: column name

    :type outFile: str
    :param outFile: output FITS file name (with larger records)

    :type ext: int
    :param ext: extension number

    :type datacol: np array
    :param datacol: data for the column (same length as existing columns)


    """

    f = fits.open(inFile)
    nrows0 = f[ext].header["NAXIS2"]
    f.close()

    coldim1 = str(nrows0) + 'E'
    newcol1 = fits.Column(name=colname, format=coldim1, array=datacol)
    t = fits.BinTableHDU.from_columns([newcol1, ])
    t.writeto(inFile)


def addkeys(fitsfile, ext, keynames, keyvals):
    """
    Add key to FITS file name in given extension
    :type fitsfile: str
    :param fitsfile: FITS file name

    :type ext: int
    :param ext: extension number

    :type keynames: list
    :param keynames: keyword name

    :type keyvals: list
    :param keyvals: keyword value

    """

    f = fits.open(fitsfile, mode='update')
    for i in range(len(keynames)):
        print("Open file ", fitsfile, "to add", keynames[i])
        f[ext].header[keynames[i]] = keyvals[i]
    f.close()
    print("Closed")


def tsfn(i, nproc, comm, simTimeN):
    """
    Funcion to run a single process of noise generation
    """
    if i == 0:
        tstart = 0.
    else:
        tstart = i * float(simTimeN)/10. + 0.01
    tend = float(simTimeN)/10. * (1 + i)
    outfile = "tmpbbfb" + str(i) + ".fits"
    comm += (" Streamfile=" + outfile +
             " tstart=" + str(tstart) + " tstop=" + str(tend))

    try:
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error creating ForNoise i=", i)
        raise
    updateHISTORY(outfile, comm)
    rmLastAndFirst(outfile, 0)


def simulNoise(pixName, samprate, jitter, bbfb, pulseLength,
               space, acbias, scaleFactor, samplesUp, nSgms,
               nintervals, dcmt):
    """ simulate data in input parameter space, calculate the data baseline
            and create Noise file to be ingested in SIRENA processing tasks

          :param pixName : name of extension in pixel definition file
                          SPA*, LPA1*, LPA2*, LPA3*
          :param samprate: Samprate value with respect to baseline
              of 156250 Hz: "" (baseline), "samprate2" (half_baseline)
          :param jitter: jitter option ("" for no_jitter and
                                        "jitter" for jitter)
          :param bbfb: ("") for dobbfb=n or ("bbfb") for dobbfb=y
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
          :param dcmt: decimation factor for xifusim jitter simulation
          :return fits file with Noise spectra in specified space and sampling

    """
    # ----GLOBAL VARIABLES -------------
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA"
    simSIXTEdir = EURECAdir + "/testHarness/simulations/SIXTE"
    cwd = os.getcwd()

    preBuffer = 0
    diffth = 0

    # ---- param dependent variables
    triggerSize = pulseLength + 10

    # samprate
    smprtStr = ""
    idxsmp = sampids.index(samprate)
    smprtStr = sampStrs[idxsmp]
    samplingrate = sampfreqs[idxsmp]

    # set simulation time based on number of intervals required
    simTimeN = math.ceil(nintervals * triggerSize / samplingrate)

    # jitter
    jitterStr = ""
    offset = ""
    if jitter == "jitter":
        jitterStr = "_jitter"
        jitterStrPix = "_jitter"
        offset = " offset=-1"
    if dcmt > 1:
        jitterStr = "_jitter_dcmt" + str(dcmt)

    # bbfb
    bbfbStr = ""
    XMLxifusim = XIFUSIMinst + "1pix_noise_nobbfb.xml"
    assert not(jitter == "" and bbfb == "bbfb"),\
        "Error: BBFB will produce jitter"

    if bbfb == "bbfb":
        XMLxifusim = XIFUSIMinst + "1pix_noise_bbfb.xml"
        bbfbStr = "_bbfb"

    # energyMethod
    spacePar = "NONE"
    if space == "R":
        spacePar = "I2R=I2R"
    elif space == "RALL":
        spacePar = "I2R=RALL"
    elif space == "RNOL":
        spacePar = "I2R=RNOL"
    elif space == "RFITTED":
        spacePar = "I2R=RFITTED"
    elif space == "ADC":
        spacePar = "I2R=I"

    xifusim = "xifusim" + pixName
    wdir = simSIXTEdir + "/NOISE/" + xifusim
    os.chdir(wdir)

    # define files
    rootN = ("forNoise" + str(pulseLength) + "samples_" + xifusim +
             "_" + str(simTimeN) + "s")
    noiseFile = ("noise" + str(pulseLength) + "samples_" + xifusim + "_B0_" +
                 space + smprtStr + jitterStr + bbfbStr + ".fits")
    pixFileN = rootN + smprtStr + jitterStrPix + ".piximpact"
    fitsFileN = rootN + smprtStr + jitterStr + bbfbStr + ".fits"

    # -------------------------------------------------------------------------
    # get NOISE file: process empty trigger file
    # -------------------------------------------------------------------------
    if not os.path.isfile(pixFileN):
        print("\nTESCONSPILEUP: Generating no-events file for NOISE...")
        comm = ("tesconstpileup PixImpList=" + pixFileN +
                " sample_freq=" + str(samplingrate) +
                " XMLFile=" + XMLsixte +
                " tstop=" + str(simTimeN) + offset +
                " energy=0 pulseDistance=1 TriggerSize=" +
                str(triggerSize) + " clobber=yes")
        args = shlex.split(comm)
        try:
            check_call(args, stderr=STDOUT)
        except RuntimeError:
            print("Error running tool tesconstpileup (NOISE):")
            print(comm)
            os.chdir(cwd)
            rmtree(tmpDir)
            raise

        print("TESCONSPILEUP: ...................................END")

    if not os.path.isfile(fitsFileN):
        print("\nxifusim: Generating file for NOISE:", fitsFileN)
        tstart = 0.
        tstop = float(simTimeN)

        comm = ("xifusim PixImpList=" + pixFileN +
                " Streamfile=" + fitsFileN +
                " sample_rate=" + str(samplingrate) +
                " tstart=" + str(tstart) + " tstop=" + str(tstop) +
                " trig_reclength=" + str(triggerSize) +
                " trig_n_pre=" + str(preBuffer) +
                " trig_thresh=" + str(diffth) +
                " trig_n_suppress=" + str(pulseLength) +
                " acbias=" + acbias +
                " decimate_factor=" + str(dcmt) +
                " clobber=yes XMLfilename=" + XMLxifusim)

        args = shlex.split(comm)
        try:
            print("Running tool xifusim (NOISE):", comm)
            check_call(args, stderr=STDOUT)
        except RuntimeError:
            print("Error running tool xifusim (NOISE):")
            print(comm)
            os.chdir(cwd)
            rmtree(tmpDir)
            raise
        print("xifusim: ................................END")
        # sys.exit()

        # it seems (27/03/2019) that there are no problems now
#        try:
#            # rm first record in fits file (noisy)
#            print("rm first record in fits file (noisy)")
#            rmLastAndFirst(fitsFileN, 'TESRECORDS', 0)
#            updateHISTORY(fitsFileN, comm)
#        except RuntimeError:
#            print("Error making up trigger file  (NOISE):")
#            os.chdir(cwd)
#            rmtree(tmpDir)
#            raise

        # and Add TRIGGSZ to xifusim simulated file
        fN = fits.open(fitsFileN, mode='update')
        # fNhdr = fN[1].header
        fN[0].header["TRIGGSZ"] = triggerSize
        fN[0].header["DELTAT"] = dcmt * fN['TESRECORDS'].header["DELTA_T"]
        fN['TESRECORDS'].header["TRIGGSZ"] = triggerSize
        fN['TESRECORDS'].header["DELTAT"] = \
            dcmt * fN['TESRECORDS'].header["DELTA_T"]
        fN.close()

    print("\nGENNOISESPEC: Generating NOISE spectrum file "
          "in (", space, " space)")
    comm = ("gennoisespec --inFile=" + fitsFileN +
            " --" + spacePar +
            " --outFile=" + noiseFile +
            " --intervalMinSamples=" + str(pulseLength) +
            " --nintervals=" + str(nintervals) +
            " --scaleFactor=" + str(scaleFactor) +
            " --samplesUp=" + str(samplesUp) +
            " --nSgms=" + str(nSgms) +
            " --pulse_length=" + str(pulseLength) +
            " --clobber=yes verbosity=0")
    print("             ", comm)
    args = shlex.split(comm)
    try:
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running tool gennoise (NOISE):")
        print(comm)
        os.chdir(cwd)
        rmtree(tmpDir)
        raise
    print("GENNOISESPEC: .......................................END")
    print("Noise file ", noiseFile, " has been successfully created\n")
    os.chdir(cwd)


def simulSingles(pixName, monoEkeV, acbias, samprate, jitter, noise,
                 bbfb, nSimPulses, pulseLength, dcmt):
    """
    :param pixName: Extension name in the FITS pixel definition file
                    (SPA*, LPA1*, LPA2*, LPA3*)
    :param monoEkeV: Monochromatic energy (keV) of input simulated pulses
    :param acbias: Operating Current (AC if acbias=yes or DC if acbias=no)
    :param samprate: Samprate value with respect to baseline of 156250 Hz:
                    "" (baseline), "samprate2" (half_baseline)
    :param jitter: jitter option ("" for no_jitter and "jitter" for jitter)
    :param noise: noise option ("" for noise and "nonoise" for nonoise)
    :param bbfb: ("") for dobbfb=n or ("bbfb") for dobbfb=y
    :param nSimPulses: number of pulses to be simulated
    :param pulseLength: length of pulses to calculate trigger size
    :param dcmt: decimation factor for xifusim jitter simulations
    :return: files with simulated SINGLES
    """

    cwd = os.getcwd()
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA"
    ERESOLdir = EURECAdir + "/ERESOL"
    PAIRSdir = ERESOLdir + "/PAIRS"
    preBufferSize = preBufferPulses

    # samprate
    smprtStr = ""
    idxsmp = sampids.index(samprate)
    smprtStr = sampStrs[idxsmp]
    samplingrate = sampfreqs[idxsmp]

    # deal with separations from definitions above
    idxsmp = sampids.index(samprate)
    singleSeparation = separations[idxsmp]

    XMLxifusim = XIFUSIMinst + "1pix_nobbfb.xml"

    # added to solve floating point inaccu. due to sampling rate
    # (Christian's mail 31/03/2017)
    tstart = 0.5 / samplingrate

    # jitter
    jitterStr = ""
    jitterStrPix = ""
    offset = ""
    if jitter == "jitter":
        jitterStr = "_jitter_dcmt" + str(dcmt)
        jitterStrPix = "_jitter"
        offset = " offset=-1"

    # noise
    noiseStr = ""
    simnoise = " simnoise=y"
    if noise == "nonoise":
        noiseStr = "_nonoise"
        simnoise = " simnoise=n"

    # bbfb
    bbfbStr = ""
    dobbfb = " dobbfb=n"
    if jitter == "jitter":
        dobbfb = ""
    if bbfb == "bbfb":
        bbfbStr = "_bbfb"
        XMLxifusim = XIFUSIMinst + "1pix_bbfb.xml"
        dobbfb = (" dobbfb=y decimation_filter=y")

    xifusim = "xifusim" + pixName
    SIMFILESdir = PAIRSdir + "/" + xifusim

    triggerSizeTC = preBufferSize + singleSeparation + singleSeparation + 1000
    triggerSizeXF = preBufferSize + pulseLength + 1000
    triggerSuppXF = triggerSizeXF - preBufferSize
    triggerSizeTC = int(triggerSizeTC)
    diffth = 60

    # calculate sim time to have at least nSimPulses pulses:
    #   simTime= (nSimPulses/2)recs * (triggerSize sam/rec)/(samprate sam/s)
    simTime = nSimPulses / 2. * triggerSizeTC / samplingrate
    simTime = '{0:0.0f}'.format(simTime)

    # for piximpact:
    root0 = "sep" + str(singleSeparation) + "sam_" + simTime + "s_" + \
        str(monoEkeV) + "keV" + smprtStr + jitterStrPix
    # for fits:
    root = ("sep" + str(singleSeparation) + "sam_" + str(nSimPulses) + "p_" +
            str(monoEkeV) + "keV" + smprtStr + jitterStr + noiseStr +
            bbfbStr)
    pixFile = (PAIRSdir + "/PIXIMPACT/" + root0 + "_trSz" +
               str(triggerSizeTC) + ".piximpact")
    fitsFile = SIMFILESdir + "/" + root + ".fits"
    print("-------------------------------------------\n")
    print("Simulating ", fitsFile, "\n")
    print("-------------------------------------------\n")

    if not os.path.isfile(pixFile):
        comm = ("tesconstpileup PixImpList=" + pixFile +
                " XMLFile=" + XMLsixte + " timezero=" + str(tstart) +
                " tstop=" + str(simTime) + " energy=" + str(monoEkeV) +
                offset + " pulseDistance=" + str(singleSeparation) +
                " TriggerSize=" + str(triggerSizeTC) + " clobber=yes" +
                " sample_freq=" + str(samplingrate))

        print("\n##### Runing tesconstpileup #########")
        print(comm, "\n")
        try:
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except RuntimeError:
            print("Error running tool for piximpact list generation")
            os.chdir(cwd)
            rmtree(tmpDir)
            raise
        # continue  # to simulate only piximpact files

    if os.path.isfile(fitsFile):
        # verify existing file
        errs, warns = fitsVerify(fitsFile)
        print("        Num errors/warnings=", errs, warns)
        if errs > 0:
            print("numerrs = ", errs, " for ", fitsFile,
                  ": repeating simulation")
            os.remove(fitsFile)

    if not os.path.isfile(fitsFile):
        commxifusim = ("xifusim PixImpList=" + pixFile +
                       " Streamfile=" + fitsFile +
                       " tstart=0." + " tstop=" + simTime +
                       " trig_reclength=" + triggerSizeXF +
                       " trig_n_pre=" + preBufferSize +
                       " trig_thresh=" + str(diffth) +
                       " trig_n_suppress=" + triggerSuppXF +
                       " acbias=" + acbias + " sample_rate=" +
                       str(samplingrate) + simnoise + dobbfb +
                       " decimate_factor=" + str(dcmt) +
                       " XMLfilename=" + XMLxifusim)

        print("\n##### Runing xifusim #########")
        print(commxifusim, "\n")
        try:
            args = shlex.split(commxifusim)
            check_call(args, stderr=STDOUT)
        except RuntimeError:
            print("Error running xifusim for data simulation")
            os.chdir(cwd)
            rmtree(tmpDir)
            raise
        # continue
        # update HISTORY  for xifusim run in header[0]
        updateHISTORY(fitsFile, commxifusim)

        # rm first (and LAST) record and update NETTOT
        fsim = fits.open(fitsFile)
        nrows = fsim['TESRECORDS'].header["NAXIS2"]
        assert nrows > 1, "xifusim failed: just one huge row present!"
        fsim.close()
        # continue
        try:
            print("Removing first & last row, just in case, ",
                  "and updating NETTOT")
            rmLastAndFirst(fitsFile, 1)
        except RuntimeError:
            print("Error running FTOOLS to remove initial & last rows in ",
                  fitsFile)
            os.chdir(cwd)
            rmtree(tmpDir)
            raise

        # Modify EXTNAME if not EXTNAME=RECORDS in empty trigger file
        # and Add TRIGGSZ to xifusim simulated file
        fsim = fits.open(fitsFile, mode='update')
        fsim[1].header["TRIGGSZ"] = triggerSizeXF
        fsim[1].header["DELTAT"] = dcmt * fsim[1].header["DELTA_T"]
        fsim[0].header["TRIGGSZ"] = triggerSizeXF
        fsim[0].header["DELTAT"] = dcmt * fsim[1].header["DELTA_T"]
        oldExtName = fsim[1].header["EXTNAME"]
        if not oldExtName == "RECORDS":
            print("Changing name...")
            fsim[1].header['EXTNAME'] = 'RECORDS'
            fsim.close()

    os.chdir(cwd)


def simulPairs(pixName, monoEkeV1, monoEkeV2, acbias, samprate, jitter,
               noise, bbfb, nSimPulses, pulseLength, dcmt, sepsStr):
    """
    :param pixName: Extension name in the FITS pixel definition file
                     (SPA*, LPA1*, LPA2*, LPA3*)
    :param monoEkeV1: Monochromatic energy (keV) of input first
                    simulated pulses
    :param monoEkeV2: Monochromatic energy (keV) of input second
                    simulated pulses
    :param acbias: Operating Current (AC if acbias=yes or DC if acbias=no)
    :param samprate: Samprate value with respect to baseline of 156250 Hz:
                    "" (baseline), "samprate2" (half_baseline)
    :param jitter: jitter option ("" for no_jitter and "jitter" for jitter)
    :param noise: noise option ("" for noisy and "nonoise" for no_noise)
    :param bbfb: ("") for dobbfb=n or ("bbfb") for dobbfb=y
    :param nSimPulses: number of pulses to be simulated
    :param pulseLength: length of pulses to calculate trigger size
    :param sepsStr: separations between pulses
    :return: files with simulated PAIRS
    """
    cwd = os.getcwd()
    print("cwd=", cwd)
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA"
    ERESOLdir = EURECAdir + "/ERESOL"
    PAIRSdir = ERESOLdir + "/PAIRS"
    preBufferSize = preBufferPulses

    # samprate
    smprtStr = ""
    idxsmp = sampids.index(samprate)
    smprtStr = sampStrs[idxsmp]
    samplingrate = sampfreqs[idxsmp]

    XMLxifusim = XIFUSIMinst + "1pix_nobbfb.xml"

    # deal with record separations from definitions above
    recordSeparation = separations[idxsmp]

    # added to solve floating point inaccu. due to sampling rate
    # (Christian's mail 31/03/2017)
    tstart = 0.5 / samplingrate

    # jitter
    jitterStr = ""
    jitterStrPix = ""
    offset = ""
    if jitter == "jitter":
        jitterStr = "_jitter_dcmt" + str(dcmt)
        jitterStrPix = "_jitter"
        offset = " offset=-1"

    # noise
    noiseStr = ""
    simnoise = " simnoise=y"
    if noise == "nonoise":
        noiseStr = "_nonoise"
        simnoise = " simnoise=n"

    # bbfb
    bbfbStr = ""
    assert not(jitter == "" and bbfb == "bbfb"),\
        "Error: BBFB will produce jitter"

    if bbfb == "bbfb":
        bbfbStr = "_bbfb"
        XMLxifusim = XIFUSIMinst + "1pix_bbfb.xml"

    xifusim = "xifusim" + pixName
    SIMFILESdir = PAIRSdir + "/" + xifusim

    for sepA in sepsStr:
        sep12 = int(sepA)
        triggerSizeTC = preBufferSize + sep12 + recordSeparation + 1000
        triggerSizeXF = preBufferSize + sep12 + pulseLength + 1000
        triggerSuppXF = triggerSizeXF - preBufferSize
        triggerSizeTC = int(triggerSizeTC)
        diffth = 60

        # calculate sim time to have at least nSimPulses pulses:
        #  simTime= (nSimPulses/2)recs * (triggerSize sam/rec)/(samprate sam/s)
        simTime = nSimPulses / 2. * triggerSizeTC / samplingrate
        simTime = '{0:0.0f}'.format(simTime)

        # for piximpact:
        root0 = ("sep" + sepA + "sam_" + simTime + "s_" + monoEkeV1 +
                 "keV_" + monoEkeV2 + "keV" + smprtStr + jitterStrPix)
        # for fits:
        root = ("sep" + sepA + "sam_" + str(nSimPulses) + "p_" +
                monoEkeV1 + "keV_" + monoEkeV2 + "keV" +
                smprtStr + jitterStr + noiseStr + bbfbStr)
        pixFile = (PAIRSdir + "/PIXIMPACT/" + root0 + "_trSz" +
                   str(triggerSizeTC) + ".piximpact")
        fitsFile = SIMFILESdir + "/" + root + ".fits"
        print("-------------------------------------------\n")
        print("Simulating ", fitsFile, "\n")
        print("  with ", pixFile, "\n")
        print("-------------------------------------------\n")

        if not os.path.isfile(pixFile):
            comm = ("tesconstpileup PixImpList=" + pixFile +
                    " XMLFile=" + XMLsixte + " timezero=" + str(tstart) +
                    " tstop=" + str(simTime) + " energy=" + monoEkeV1 +
                    " energy2=" + monoEkeV2 + offset + " clobber=yes" +
                    " pulseDistance=" + str(sep12) + " sample_freq=" +
                    str(samplingrate) + " TriggerSize=" +
                    str(triggerSizeTC))
            print("\n##### Runing tesconstpileup #########")
            print(comm, "\n")
            try:
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except RuntimeError:
                print("Error running TESCONSTPILEUP for piximpact list",
                      " generation with:", comm)
                os.chdir(cwd)
                rmtree(tmpDir)
                raise
            # continue  # to simulate only piximpact files

        if os.path.isfile(fitsFile):
            # verify existing file
            numerrs, numwrns = fitsVerify(fitsFile)
            if numerrs > 0:
                print("numerrs = ", numerrs, " for ", fitsFile,
                      ": repeating simulation")
                os.remove(fitsFile)

        if not os.path.isfile(fitsFile):
            commxifusim = ("xifusim PixImpList=" + pixFile +
                           " Streamfile=" + fitsFile +
                           " tstart=0." + " tstop=" + simTime +
                           " trig_reclength=" + triggerSizeXF +
                           " trig_n_pre=" + preBufferSize +
                           " trig_thresh=" + str(diffth) +
                           " trig_n_suppress=" + triggerSuppXF +
                           " acbias=" + acbias +
                           " sample_rate=" + str(samplingrate) + simnoise +
                           " decimate_factor=" + str(dcmt) +
                           " XMLfilename=" + XMLxifusim)

            print("\n##### Runing xifusim #########")
            print(commxifusim, "\n")
            try:
                args = shlex.split(commxifusim)
                check_call(args, stderr=STDOUT)
            except RuntimeError:
                print("Error running xifusim for data simulation")
                os.chdir(cwd)
                rmtree(tmpDir)
                raise
            # continue
            # update HISTORY  for tesconstpileup run in header[0]
            updateHISTORY(fitsFile, commxifusim)

            # rm first (and LAST) record and update NETTOT
            fsim = fits.open(fitsFile)
            nrows = fsim[1].header["NAXIS2"]
            assert nrows > 1, "xifusim failed: just one huge row present!"
            fsim.close()
            # continue
        try:
            print("Removing first & last row, just in case, ",
                  "and updating NETTOT")
            rmLastAndFirst(fitsFile, 1)
        except RuntimeError:
            print("Error running FTOOLS to remove initial & last rows in ",
                  fitsFile)
            os.chdir(cwd)
            rmtree(tmpDir)
            raise
        # Modify EXTNAME if not EXTNAME=RECORDS in empty trigger file
        # and Add TRIGGSZ to xifusim simulated file
        fsim = fits.open(fitsFile, mode='update')
        fsim[0].header["TRIGGSZ"] = triggerSizeXF
        fsim[0].header["DELTAT"] = dcmt * fsim[1].header["DELTA_T"]
        fsim[1].header["TRIGGSZ"] = triggerSizeXF
        fsim[1].header["DELTAT"] = dcmt * fsim[1].header["DELTA_T"]
        oldExtName = fsim[1].header["EXTNAME"]
        if not oldExtName == "RECORDS":
            print("Changing name...")
            fsim[1].header['EXTNAME'] = 'RECORDS'
            fsim.close()

    # some cleaning before exiting the function
    os.chdir(cwd)


def simulLibsGlobal(pixName, space, samprate, jitter, noise, bbfb,
                    pulseLength, libEnergies, largeFilter, nsamples,
                    nSimPulses, acbias, tstartPulse1All,
                    createLib, noiseMat, nintervals, weightMat, dcmt):
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
    :type bbfb: str
    :param bbfb: ("") for dobbfb=n or ("bbfb") for dobbfb=y
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
    :type nintervals: int
    :param nintervals: number of intervals if WEIGHTM noise used (inst.of NSD)
    :type weigthMat: str
    :param weightMat: should the Weight matrices HDU be created? (yes/no)
    :param dcmt: decimation factor for jitter in xifusim
    :return: simulated calibration pulses pulses && Global library from them

    """

    cwd = os.getcwd()
    preBufferSize = preBufferPulses
    xifusim = "xifusim" + pixName

    # samprate
    smprtStr = ""
    idxsmp = sampids.index(samprate)
    smprtStr = sampStrs[idxsmp]
    samplingrate = sampfreqs[idxsmp]

    # deal with separations from definitions above
    separation = int(separations[idxsmp])

    # Sigmas and Samples and scaleFactor for Detection
    maxFilterLength = pulseLength
    if largeFilter > 0:
        maxFilterLength = largeFilter

    samplesUp = samplesUps[idxsmp]
    samplesDown = samplesDowns[idxsmp]
    nSgms = nSigmss[idxsmp]

    XMLxifusim = XIFUSIMinst + "1pix_nobbfb.xml"

    # jitter
    jitterStr = ""
    jitterStrPix = ""
    offset = ""
    if jitter == "jitter":
        jitterStr = "_jitter"
        jitterStrPix = "_jitter"
        offset = " offset=-1"
    if dcmt > 1:
        jitterStr = "_jitter_dcmt" + str(dcmt)

    # noise
    noiseStr = ""
    simnoise = " simnoise=y"
    if noise == "nonoise":
        noiseStr = "_nonoise"
        simnoise = " simnoise=n"

    # noiseMat
    if noiseMat == "yes" and nintervals > 0:
        noiseMatStr = "_" + str(nintervals) + "int"
    # bbfb
    bbfbStr = ""
    assert not(jitter == "" and bbfb == "bbfb"),\
        "Error: BBFB will produce jitter"

    if bbfb == "bbfb":
        bbfbStr = "_bbfb"
        XMLxifusim = XIFUSIMinst + "1pix_bbfb.xml"

    # -- CALIB/SIM/NOISE/LIB/XML dirs & files --
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA"
    simSIXTEdir = EURECAdir + "/testHarness/simulations/SIXTE"
    PIXIMPACTdir = simSIXTEdir + "/LIBRARIES/PIXIMPACT"
    SIMFILESdir = simSIXTEdir + "/LIBRARIES/" + xifusim
    NOISEdir = simSIXTEdir + "/NOISE/" + xifusim
    libDir = SIMFILESdir + "/GLOBAL/" + space

    noiseFile = (NOISEdir + "/noise" + str(nsamples) + "samples_" + xifusim +
                 "_B0_" + space + smprtStr + jitterStr +
                 bbfbStr + noiseMatStr + ".fits")
    libFile = (libDir + "/libraryMultiE_GLOBAL_PL" + str(pulseLength) +
               "_" + str(nSimPulses) + "p" + smprtStr + jitterStr +
               noiseStr + bbfbStr + noiseMatStr + ".fits")
    evttmpFile = tempfile.NamedTemporaryFile()

    # added to solve floating point inaccuracies due to sampling rate
    # (Christian's mail 31/03/2017):
    tstart = 0.5/samplingrate

    # Calibration energies and Tstarts of pulses
    tstartPulse1 = dict(zip(libEnergies, tstartPulse1All))

    # calculate simTime so that nSimPulses are simulated (1 pulses per record)
    simTime = nSimPulses * (preBufferSize + separation) / float(samplingrate)
    tstop = tstart + simTime
    simTime = '{0:0.0f}'.format(simTime)
    tstop = '{0:0.0f}'.format(tstop)

    # Trigger sizes in tesconstpileup
    triggerSizeTC = preBufferSize + separation + separation + 1000
    # and xifusim ___|\_______________ooo
    triggerSizeXF = preBufferSize + maxFilterLength + 1000
    triggerSuppXF = triggerSizeXF - preBufferSize
    diffth = 60
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

        print("=============================================")
        print("Adding monochromatic energy", monoEkeV, "keV")
        print("=============================================")
        print("simTime=", simTime, "\n")
        # print("Temporary event file in: ", evttmpFile.name)
        # simulate SIXTE file with tesconstpileup + xifusim
        root0 = ("mono" + str(monoEkeV) + "_sep" + str(separation) + "_pix1" +
                 "_" + str(nSimPulses) + "p")
        root = root0 + "_" + str(pulseLength)
        # pixFile = PIXIMPACTdir +
        #          "/" + root0 + "_trSz" + str(triggerSizeTC) + ".piximpact"
        pixFile = (PIXIMPACTdir + "/" + root0 + smprtStr + jitterStrPix +
                   ".piximpact")
        simFile = (SIMFILESdir + "/" + root + smprtStr + jitterStr + noiseStr +
                   bbfbStr + ".fits")

        # -- TESCONSTPILEUP: generate impacts for well separated single pulses
        # -- xifusim: simulate well separated pulses --
        if not os.path.isfile(simFile):
            print("Simulating & triggering pulses with xifusim to ", simFile)
            comm = ("tesconstpileup PixImpList=" + pixFile +
                    " XMLFile=" + XMLsixte + " timezero=" + str(tstart) +
                    " tstop=" + str(tstop) + " energy=" + str(monoEkeV) +
                    " pulseDistance=" + str(separation) + offset +
                    " TriggerSize=" + str(triggerSizeTC) + " clobber=yes" +
                    " sample_freq= " + str(samplingrate))
            print("Generating PIXIMPACT ", pixFile, " running:\n", comm)
            try:
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except RuntimeError:
                print("Error running task for piximpact list generation")
                os.chdir(cwd)
                rmtree(tmpDir)
                raise

            # continue # to only generate piximpact
            commxifusim = ("xifusim PixImpList=" + pixFile +
                           " Streamfile=" + simFile +
                           " tstart=0 tstop=" + str(simTime) +
                           " trig_reclength=" + triggerSizeXF +
                           " trig_n_pre=" + preBufferSize +
                           " trig_thresh=" + str(diffth) +
                           " trig_n_suppress=" + triggerSuppXF +
                           " acbias=" + acbias +
                           " sample_rate=" + str(samplingrate) +
                           " decimate_factor=" + str(dcmt) + simnoise +
                           " XMLfilename=" + XMLxifusim)

            print("Running xifusim for simulation\n", commxifusim)
            try:
                args = shlex.split(commxifusim)
                check_call(args, stderr=STDOUT)
            except RuntimeError:
                print("Error running xifusim for simulation\n")
                os.chdir(cwd)
                rmtree(tmpDir)
                raise

            rmLastAndFirst(simFile, 1)

            # update HISTORY in header[0]
            dateTime = strftime("%Y-%m-%d %H:%M:%S", gmtime())
            history = ("Created & Updated by simulLibsGlobal.py on " +
                       dateTime + "with command: " + commxifusim)
            updateHISTORY(simFile, history)

        # print("Antes de evaluar createLib\n")
        if not createLib:
            print("NO Library creation\n")
            continue  # (to just simulate pulses files)

        print("Creating library\n")
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
                " opmode=0 clobber=yes intermediate=0" +
                " monoenergy=" + str(monoEeV) + " EventListSize=1000" +
                " XMLFile=" + XMLsixte +
                " hduPRECALWN=" + weightMat + " hduPRCLOFWM=" + noiseMat)

        args = shlex.split(comm)
        print("SIRENA reconstruction to add a line to the library, "
              "running command:\n", comm)
        try:
            check_call(args, stderr=STDOUT)
        except RuntimeError:
            print("Error running SIRENA to add new line to library")
            os.chdir(cwd)
            rmtree(tmpDir)
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


def reconstruct(pixName, labelLib, samprate, jitter, dcmt, noise, bbfb, Lc,
                mono1EkeV, mono2EkeV, reconMethod, filterLength,
                nsamples, pulseLength, nSimPulses, fdomain, detMethod,
                tstartPulse1, tstartPulse2, nSimPulsesLib, coeffsFile,
                libTmpl, simDir, outDir, detSP, pB, LbT, s0, lags,
                filterct, B0, sepsStr):

    """
    :param pixName: Extension name for FITS pixel definition file
                    (SPA*, LPA1*, LPA2*, LPA3*)
    :param labelLib: Label identifying the library
                    ( multilib, multilibOF, fixedlib1 )
    :param samprate: Samprate value with respect to baseline of 156250 Hz:
                    "" (baseline), "samprate2" (half_baseline),
                    samprate4 (quarter baseline)
    :param jitter: jitter option ("" for no_jitter and "jitter" for jitter or
                                  "jitter_M82" for M82 special case)
    :param dcmt: decimation factor for xifusim jitter
    :param bbfb: ("") for dobbfb=n or ("bbfb") for dobbfb=y
    :param Lc: ("" for L=Lcrit; otherwise, L/Lcrit)
    :param noise: simulations done with ("") or without ("nonoise") noise
    :param mono1EkeV: Monochromatic energy (keV)
                        of input primary simulated pulse
    :param mono2EkeV: Monochromatic energy (keV)
                        of input secondary simulated pulse
    :param reconMethod: Energy reconstruction Method
                        (OPTFILT, WEIGHT, WEIGHTN, I2R, I2RALL, I2RNOL,...)
    :param filterLength: if set to 0 (default) filters are taken
                        from library in BASE2, otherwise filterStartegy=FIXED
    :param nsamples: noise length samples
    :param pulseLength: pulse length in samples
    :param nSimPulses: number of simulated pulses
                        (to locate filenames to be analysed)
    :param nSimPulsesLib: number of simulated pulses for the library
    :param fdomain: Optimal Filtering domain (F(req) or T(ime))
    :param detMethod: Detection method (AD or STC)
    :param tstartPulse1: Start sample for first pulse (PreBufferSize)
    :param tstartPulse2:
                    if /=0 Start sample for second pulse
                    if = -1 modified to (tstartPulse1 + separation)
                    if = 0 && tstartPulse1 /=0 => no secondary pulse
    :param nSimPulsesLib: number of simulated pulses for the library
    :param libTmpl: label to identify the library regarding
                    the templates used (LONG or SHORT)
    :param simDir: directory for simulated files
    :param outDir: directory for resulting evt and json files
    :param detSP: 1 secondary pulses will be detected (default), 0 otherwise
    :param pB: preBuffer value for optimal filters
    :param LbT: Time (s) to average baseline to be subtracted
    :param s0: Optimal Filters' SUM should be '0'?: s0=0 (NO); s0=1 (YES)
    :param lags: Do parabola fit? lags=1 (YES), lags=0 (NO)
    :param filterct: Filters central part have been replaced by ct value?
                    or use derived filters from largest 8192 filter (fit)))
    :param B0: if B0>0 processing will be done by B0 (baseline subtraction)
                Otherwise, F0 (bin 0 freq.)
    :param sepsStr: blank spaces separated list of pulses separations
    :return: file with energy resolutions for the input pairs of pulses
    """

    # --- Define some initial values and conversions ----
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA/"
    simSIXTEdir = EURECAdir + "testHarness/simulations/SIXTE"

    # samprate
    smprtStr = ""
    idxsmp = sampids.index(samprate)
    smprtStr = sampStrs[idxsmp]

    # separations
    separation = separations[idxsmp]
    if sepsStr == "":
        sepsStr = [separation]

    # detection related params
    samplesUp = samplesUps[idxsmp]
    samplesDown = samplesDowns[idxsmp]
    nSgms = nSigmss[idxsmp]

    # jitter
    jitterStr = ""
    if jitter == "jitter":
        jitterStr = "_jitter"
    if jitter == "jitter" and dcmt > 1:
        jitterStr = "_jitter_dcmt" + str(dcmt)
    if jitter == "jitter_M82":
        jitterStr = "_jitter_M82"

    # noise
    noiseStr = ""
    if noise == "nonoise":
        noiseStr = "_nonoise"

    # bbfb
    bbfbStr = ""
    if bbfb == "bbfb":
        bbfbStr = "_bbfb"
    if "NewPar" in bbfb or "040" in bbfb:
        bbfbStr = "_" + bbfb
        bbfb = "bbfb"

    # optimal filters' preBuffer
    pBStr = ""
    if pB > 0:
        pBStr = "_pB" + str(pB)
    # baseline average samples
    LbTStr = ""
    if float(LbT) > 0:
        LbTStr = "_LbT" + str(LbT)
    # optimal filter's SUM
    s0Str = ""
    s0Param = ""
    if s0 == 1:
        s0Str = "_Sum0Filt"
        s0Param = " Sum0Filt=1"
    # lags
    lagsStr = ""
    if lags == 0:
        lagsStr = "_nolags"

    # optimal filters' (constant central part or derived-from-8192 filters)
    ctStr = ""
    if filterct:
        ctStr = "_" + filterct
    # optimal filters B0 (baseline subtraction) or F0 (bin f=0) processing
    B0str = ""
    FMparam = " FilterMethod=F0"
    if B0 > 0:
        B0str = "_B0-" + str(B0)
        FMparam = " FilterMethod=B0"

    # Lcrit
    LcStr = ""
    if not Lc == "":
        LcStr = "_" + str(Lc) + "Lc"

    # pulses start
    TRIGG = ""
    if tstartPulse1 > 0:
        TRIGG = "NTRIG"
    else:
        TRIGG = detMethod

    # pulses types
    xifusim = "xifusim" + pixName

    # data space and reconstruction method
    space = "ADC"
    if "I2R" in reconMethod:
        space = reconMethod.lstrip('I2')
        space = space.split('NM')[0]

    ofnoise = "NSD"
    reconMethod2 = ""
    if "NM" in reconMethod:
        ofnoise = "WEIGHTM"
        reconMethod2 = "NM" + reconMethod.split('NM')[1]
        reconMethod = reconMethod.split('NM')[0]

    # if "WEIGHT" in reconMethod:
    #    lags = 0
    # libraries
    OFLib = "no"
    OFstrategy = ""
    # noiseDir = simSIXTEdir + "/NOISE/" + xifusim
    # noiseFile = (noiseDir + "/noise" + str(nsamples) + "samples_" +
    #             xifusim + "_B0_" + space + smprtStr +
    #             jitterStr + bbfbStr + LcStr + ".fits")
    # noiseParam = " NoiseFile=" + noiseFile
    if "OF" in labelLib:
        OFLib = "yes"
        OFstrategy = " OFStrategy=FIXED OFLength=" + str(filterLength)
        noiseParam = ""
        # if "WEIGHT" in reconMethod:
        #    noiseParam = " NoiseFile=" + noiseFile

    # -- LIB & NOISE & SIMS & RESULTS dirs and files ----------
    # simDir = ERESOLdir + "PAIRS/" + xifusim
    # resultsDir = resultsDir.rstrip('\\')
    # resultsDir = resultsDir.lstrip('\\')
    # if ("gainScale" in resultsDir) or ("base" in resultsDir):
    #    simDir += "/" + resultsDir + "/"

    # outDir = ERESOLdir + "PAIRS/eresol" + pixName + "/" + resultsDir
    libDirRoot = simSIXTEdir + "/LIBRARIES/" + xifusim
    libDir = libDirRoot + "/GLOBAL/" + space + "/"
    if libTmpl == "SHORT":
        libTmpl = "_SHORT"
    else:
        libTmpl = ""

    if 'multilib' in labelLib:
        libFile = (libDir + "/libraryMultiE_GLOBAL_PL" + str(nsamples) +
                   "_" + str(nSimPulsesLib) + "p" + smprtStr +
                   jitterStr + noiseStr + bbfbStr + LcStr + ".fits")
        filtEeV = 6000.  # eV
    elif 'fixedlib' in labelLib:  # fixedlib1,...
        fixedEkeV = labelLib.replace("OF", "")[8:]
        filtEeV = float(fixedEkeV) * 1E3  # eV
        # detection to be performed -> require different models:
        if tstartPulse1 == 0 and detMethod == "AD":
            libFile = (libDir + "/libraryMultiE_GLOBAL_PL" + str(nsamples) +
                       "_" + str(nSimPulsesLib) + "p" + smprtStr + B0str +
                       jitterStr + noiseStr + bbfbStr + LcStr + pBStr +
                       ctStr + ".fits")
        else:
            libFile = (libDir + "/library" + fixedEkeV + "keV_PL" +
                       str(nsamples) + "_" + str(nSimPulsesLib) + "p" +
                       smprtStr + B0str + jitterStr + noiseStr + bbfbStr +
                       LcStr + pBStr + ctStr + ".fits")
    libFile = libFile.replace(".fits", libTmpl + ".fits")

    root = ''.join([str(nSimPulses), 'p_SIRENA', str(nsamples), '_pL',
                    str(pulseLength), '_', mono1EkeV, 'keV_', mono2EkeV,
                    'keV_', TRIGG, "_", str(fdomain), '_',
                    str(labelLib), '_', str(reconMethod), str(filterLength),
                    reconMethod2, B0str, pBStr, LbTStr, smprtStr, jitterStr,
                    noiseStr, bbfbStr, LcStr, s0Str, lagsStr, ctStr])
    if mono2EkeV == "0":
        root = ''.join([str(nSimPulses), 'p_SIRENA', str(nsamples), '_pL',
                        str(pulseLength), '_', mono1EkeV, 'keV_', TRIGG, "_",
                        str(fdomain), '_', str(labelLib),
                        '_', str(reconMethod), str(filterLength),
                        reconMethod2, B0str, pBStr, LbTStr, smprtStr,
                        jitterStr, noiseStr, bbfbStr, LcStr, s0Str,
                        lagsStr, ctStr])

    eresolFile = "eresol_" + root + ".json"
    eresolFile = eresolFile.replace(".json", libTmpl+".json")

    os.chdir(outDir)

    # ------------------------------------
    # --- Process input data files -------
    # ------------------------------------
    #
    for sep12 in sepsStr:
        sep = float(sep12)
        # calculate tstartPulse2 using separation:
        if tstartPulse1 > 0 and tstartPulse2 == -1:
            tstartPulse2 = tstartPulse1 + int(sep)

        inFile = (simDir + "/sep" + sep12 + "sam_" + str(nSimPulses) +
                  "p_" + str(mono1EkeV) + "keV_" + str(mono2EkeV) + "keV" +
                  smprtStr + jitterStr + noiseStr + bbfbStr + LcStr +
                  ".fits")
        if mono2EkeV == "0":  # for single Pulses
            inFile = (simDir + "/sep" + sep12 + "sam_" + str(nSimPulses) +
                      "p_" + str(mono1EkeV) + "keV" +
                      smprtStr + jitterStr + noiseStr + bbfbStr + LcStr +
                      ".fits")

        evtFile = "events_sep" + sep12 + "sam_" + root + ".fits"
        evtFile = evtFile.replace(
                jitterStr + noiseStr + bbfbStr + LcStr + lagsStr + ".fits",
                libTmpl + jitterStr + noiseStr + bbfbStr + LcStr + lagsStr +
                ".fits")

        print("=============================================")
        print("RECONSTRUCTING ENERGIES.....................")
        print("Working in:", outDir)
        print("Using file: ", inFile)
        print("Using library: ", libFile)
        # print("Using noisefile: ", noiseFile)
        print("Setting evtFile: ", evtFile)
        print("=============================================")

        # -- SIRENA processing -----

        if os.path.isfile(evtFile):
            if os.stat(evtFile).st_size == 0:
                os.remove(evtFile)
                print("Event file", evtFile,
                      " already DOES exist but size 0: removing...")
            else:
                print("Event file", evtFile,
                      " already DOES exist: recalculating FWHMs...")

        if not os.path.isfile(evtFile):
            ndetpulses = 0
            print("\nRunning SIRENA for detection & reconstruction")
            print("filtEev=", filtEeV)
            print("XMLFile=", XMLsixte)
            comm = ("tesreconstruction Recordfile=" + inFile +
                    " TesEventFile=" + evtFile +
                    " Rcmethod='SIRENA'" +
                    " PulseLength=" + str(pulseLength) +
                    " LibraryFile=" + libFile +
                    " scaleFactor=" + str(scaleFactor) +
                    " samplesUp=" + str(samplesUp) + " nSgms=" + str(nSgms) +
                    " samplesDown=" + str(samplesDown) +
                    " opmode=1 " + noiseParam +
                    " OFLib=" + OFLib + " FilterDomain=" + fdomain +
                    " detectionMode=" + detMethod + " detectSP=" + str(detSP) +
                    " clobber=yes" +
                    " EventListSize=1000" + " EnergyMethod=" + reconMethod +
                    " LagsOrNot=" + str(lags) +
                    " tstartPulse1=" + str(tstartPulse1) +
                    " tstartPulse2=" + str(tstartPulse2) +
                    " OFNoise=" + ofnoise + " XMLFile=" + XMLsixte +
                    " filtEeV=" + str(filtEeV) + OFstrategy +
                    " preBuffer=" + str(pB) +
                    " LbT=" + str(LbT) + s0Param + FMparam)
            try:
                print(comm)
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except RuntimeError:
                print("Error running SIRENA for detection & reconstruction:")
                print(comm)
                rmtree(tmpDir)
                raise

            # Check that detection went well
            evtf = fits.open(evtFile)
            nTrigPulses = evtf[1].header["NETTOT"]
            ndetpulses = evtf[1].header["NAXIS2"]
            print("Checking DETECTION............................")
            assert ndetpulses > 0, "Empty evt file (%s): nrows=0 " % evtFile
            print("Reconstruction finished => sim/det: ",
                  nTrigPulses, "/", ndetpulses)
            evtf.close()

    return smprtStr, jitterStr, noiseStr, bbfbStr, LcStr, pBStr, LbTStr,\
        s0Str, lagsStr, ctStr, B0str, evtFile, eresolFile


def convertEnergies(inFile, outFile, coeffsFile, alias):
    """
    :param inFile: absolute path to file with reconstructed
                    (non-calibrated) energies
    :param coeffsFile: file with coefficients of polynomial
                    fit to gain curves from polyfit2bias.R
                    or JSON file with data points for spline fitting
    :param alias: string to select reconstruction type in
                    the coefficients table
    :param outFile: file with reconstructed/calibrated energies
                    for the input pulses
    """

    # ------------------------------------
    # --- Process input data file  -------
    # ------------------------------------

    f = fits.open(inFile, memmap=True)
    nrows = f[1].header["NAXIS2"]
    assert nrows > 0, "Empty evt file (%s): nrows=0 " % inFile

    # read Erecons (SIGNAL) column in numpy array (in keV)
    ftab = f[1].data
    EreconKeV = np.array(ftab['SIGNAL'])
    # reconPhase = np.array(ftab['PHI'])
    reconPhase = np.array(ftab['PHI']) + np.array(ftab['LAGS'])

    # calculate corrected energies with polyfit coeffs
    # 1) check if pulses are of different length (different gain scale)
    print("...Calculating corrected energies for pulses in ", inFile)
    EcorrKeV = enerToCalEner(EreconKeV, reconPhase, coeffsFile, alias)

    # close input file
    f.close()
    del f[1].data

    # ------------------------------------
    # --- Create and populate output file
    # ------------------------------------
    copy(inFile, outFile)
    f = fits.open(outFile, memmap=True, mode='update')
    hdr = f[0].header
    hdr['HISTORY'] = ("Updated by convertEnergies.py to correct "
                      "energies by gainscale using", coeffsFile)
    ftab = f[1].data
    ftab['SIGNAL'] = EcorrKeV
    f.close()
    del f[1].data


def find_nearest_idx(array, values):
    """
    finds the index where a numpy array has the element with the closest value
    and corrects for borders.
    Value is on the left of the index
    Designed for events simulation crossmatch
    :param array: array of values
    :param values: values to locate in array
    :return: indices of array elements with the closest values
    """

    # if 1 value => returns 1 numpy integer:
    idxs = np.searchsorted(array, values, side="left")

    # convert to numpy array if idxs is an integer
    if idxs.size == 1:
        idxs = np.array([idxs])
    if values.size == 1:
        values = np.array([values])

    for i in range(0, idxs.size):
        index = idxs[i]
        index1 = index-1
        if (index > 0 and (index == len(array) or
                           math.fabs(values[i] - array[index1]) <
                           math.fabs(values[i] - array[index]))):
            idxs[i] -= 1

    return idxs


def arrivalPhase(impF, evtF, samprate, arrmethod):
    """
        Calculate the arrival phase in samples of the simulated photons:
            distance from true nonjitter simulated time and detected time
        :param impF: FITS file with impact times(from simput or tesconstpileup)
        :param evtF: FITS file with SIRENA-detected events
        :param samprate: sampling rate (Hz)
        :param arrmethod: integer for arrival phase method calculation
             0: jitter impact time (imp) - no_jitter impact time(imp0)
             1: (evtTimes - impTimes + (impTimes-imp0Times))
             2: (evtTimes - imp0Times)
             3: round to the closest sample
             4: round always to previous sample
             5: phi (distance to trigger point) It can be >1 or <-1
             6: phi (distance to trigger point) + n (additional lags). Ex.:
                *    | *      *   => phi=-0.3, n=0
                *      * |    *   => phi= 0.2, n=0
                *      *      *    | *  => phi=0.8, n=1
        :return arrPhase: list
    """
    arrPhase = [None] * 7
    # PIXIMPACT
    imp = fits.open(impF, memmap=True)
    impTab = imp[1].data
    impTimes = impTab['TIME']
    imp.close()

    # SIRENA reconstruction
    evt = fits.open(evtF, memmap=True)
    evtTab = evt[1].data
    evtTimes = evtTab['TIME']
    evtPHI = evtTab['PHI']
    evtN = evtTab['LAGS']
    evt.close()

    clsidx = find_nearest_idx(impTimes, evtTimes)

    # Offsets with no jitter times (known from simulations)
    # =====================================================
    # offset with imp times
    off1 = (evtTimes - impTimes[clsidx]) * samprate
    # offset of imp times with nonjitter times: dec part of imptime
    off2 = np.modf(impTimes[clsidx] * samprate)[0]
    # ArrPhase1: (evtTimes - impTimes + (impTimes-imp0Times))
    arrPhase[1] = off1 + off2

    # PIXIMPACT NO JITTER (if it exists):
    imp0F = impF.replace("_jitter_", "_")
    # print("impF=", impF)
    # print("imp0F=", imp0F)
    # print("evtF=", evtF)
    if os.path.isfile(imp0F):
        imp0 = fits.open(imp0F, memmap=True)
        imp0Tab = imp0[1].data
        imp0Times = imp0Tab['TIME']
        imp0.close()
        arrPhase[0] = (impTimes - imp0Times)*samprate
        clsidx0 = find_nearest_idx(imp0Times, evtTimes)
        # arrPhase 2 should be the same that arrPhase1
        arrPhase[2] = (evtTimes - imp0Times[clsidx0]) * samprate
        # print("max(Arr1-Arr2)=", max(abs(arrPhase[1]-arrPhase[2])))
        # print("Arr1=", arrPhase[1])
        # print("Arr2=", arrPhase[2])
        # assert np.allclose(arrPhase[1], arrPhase[2], atol=5E-6), \
        #    "Incompatible arrival Phases calculations"
    else:
        if arrmethod in (0, 2):
            raise ValueError('Arrival method cannot be 0 or 2 if',
                             'piximpact0 does not exist')

    # Offset: round to the closest sample
    # ====================================
    arrPhase[3] = evtTimes*samprate - np.round(evtTimes*samprate, 0)
    # Offset: round always to previous sample
    # ========================================
    arrPhase[4] = evtTimes*samprate - np.floor(evtTimes*samprate)
    # arrPhase4 = np.modf(evtTimes*samprate)[0]  # dec part of evtTime
    # return arrPhase0[0:19]
    arrPhase[6] = evtPHI
    arrPhase[5] = evtPHI + evtN
    arrival = arrPhase[arrmethod]
    return arrival


def fitsVerify(fitsfile):
    """
    :param fitsfile: input file to be verified
    :return: number of errors+ warnings in verification process
    """

    def getNumerrs(file):
        infile = open(file, "r")
        nume = 0
        numw = 0
        for line in infile:
            if line.find("Verification found") >= 0:
                words = line.split()
                nume = words[6]
                numw = words[3]
        infile.close()
        return int(nume), int(numw)

    tmpFile = tmpDir + "/pp.stdout"
    commVerify = "fverify infile=" + fitsfile + " outfile=" + \
        tmpFile + " clobber=yes"
    try:
        args = shlex.split(commVerify)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running fverify command")
        raise
    numerrs, numwrns = getNumerrs(tmpFile)
    return numerrs, numwrns

    # does not work properly if prstat=no in command line
    # getNumerrs = "pget fverify numerrs"
    # getNumwrns = "pget fverify numwrns"
    # numerrs = 0
    # numwrns = 0
    # try:
    #     args = shlex.split(getNumerrs)
    #     numerrsStr = check_output(args, stderr=STDOUT)
    #     print("numerrsStr=", numerrs)
    # except:
    #     print("Error running fverify(2) to get numerrs")
    #     rmtree(tmpDir)
    #     raise
    # try:
    #     args = shlex.split(getNumwrns)
    #     numwrnsStr = check_output(args, stderr=STDOUT)
    #     print("numwrnsStr=", numwrnsStr)
    #     numwrns = int(numwrnsStr)
    #     print("numwrns=", numwrns)
    # except:
    #     print("Error running fverify(2) to get numwrns")
    #     rmtree(tmpDir)
    #     raise
    # if os.path.isfile("pp.stdout"):
    #    os.remove("pp.stdout")
    # if fverify detects errors/warnings, return numerrs+numwrns
    # return numerrs+numwrns


def rmLastAndFirst(simfile, ext, ppr):
    """
    rm first (and LAST) record of SIXTE/xifusim simulated file and
    update NETTOT (first starts high and last can be cut)
    (see Christian's email from 19 Jan 2017 @ EURECA).
    Also update number of pulses in NETTOT keyword

    :type simfile: str
    :param simfile: simulated file where cleaning must be done
    :type ext: str
    :param ext: extension name
    :type ppr: int
    :param ppr: pulses per record in simulations
    """

    fsim = fits.open(simfile, mode='update')
    # nrows = fsim[1].header["NAXIS2"]
    nrows = fsim[ext].data.shape[0]

    # remove first and last row with astropy: it does not read ADC properly
    # binTable = fsim[1].data
    # binTable = binTable[2:nrows-1]
    # fsim.flush()
    fsim.close()

    assert nrows > 1, "Xifusim failed for %s: just one row present " % simfile
    tmpFile = tmpDir + "/pp.fits"
    copy(simfile, tmpFile)

    # try to remove last row
    comm1 = "fdelrow infile=" + simfile + "['" + ext + "'] firstrow=" + \
        str(nrows) + " nrows=1 confirm=no proceed=yes"
    try:
        args = shlex.split(comm1)
        check_call(args, stderr=STDOUT)
        print("     Final row removed in ", simfile)
    except RuntimeError:
        print("Error running ", comm1, " to remove final row in ", simfile)
        raise
    errs1, warns1 = fitsVerify(simfile)
    print("        Num errors/warnings=", errs1, warns1)
    updateHISTORY(simfile, comm1)
    # update NETTOT
    fsim = fits.open(simfile, mode='update')
    nrows2 = fsim[ext].header['NAXIS2']
    nettot = nrows2 * ppr  # new number of pulses (==2*nofrecords)
    fsim[ext].header['NETTOT'] = nettot
    fsim.close()

    # remove first row
    if errs1 > 0:
        copy(tmpFile, simfile)
        print("Errors running FTOOLS to remove final row (abandon) in ",
              simfile, " (returned to original file)")
    elif warns1 > 0:
        print("Warnings running FTOOLS to remove final row (abandon) in ",
              simfile, " (keep modified file)")
    else:
        comm2 = ("fdelrow infile=" + simfile + "['" + ext +
                 "'] firstrow=1 nrows=1 confirm=no proceed=yes")
        try:
            args = shlex.split(comm2)
            check_call(args, stderr=STDOUT)
            print("     Initial row removed in ", simfile)
        except RuntimeError:
            print("Error running ", comm2, " to remove initial row in ",
                  simfile)
            raise
        errs2, warns2 = fitsVerify(simfile)
        print("        Num errors/warnings=", errs2, warns2)
        updateHISTORY(simfile, comm2)
        # update NETTOT
        fsim = fits.open(simfile, mode='update')
        nrows2 = fsim[ext].header['NAXIS2']
        nettot = nrows2 * ppr  # new number of pulses (==2*nofrecords)
        fsim[ext].header['NETTOT'] = nettot
        fsim.close()

        if sum(fitsVerify(simfile)) > 0:
            copy(tmpFile, simfile)
            print("Errors/warnings after FTOOLS removal of initial row in ",
                  simfile)


def updateHISTORY(file, history):
    """

    :param file: File for which the keyword HISTORY in the HEADER[0]
                will be updated
    :param history: String for the HISTORY keyword (if large it will be
                splitted in several lines
    :return:

    """

    dateTime = strftime("%Y-%m-%d %H:%M:%S", gmtime())
    f = fits.open(file, mode='update')
    f[0].header['HISTORY'] = ("Created/Updated on " + dateTime +
                              "with commands: ")
    # split command line in blocks of 60 chars for HISTORY keyword
    histSplit = [history[i:i + 60] for i in range(0, len(history), 60)]
    for i in range(len(histSplit)):
        f[0].header['HISTORY'] = histSplit[i]
    f.close()
    print(histSplit)


def enerToCalEner(inEner, inPhase, coeffsFile, alias):
    """
    :param inEner: numpy array with input uncorrected energies
    :param inPhase: numpy array with input phases (jitter)
    :param coeffsFile: file with coefficients of polynomial fit to gain curves
                        from polyfit2bias.R or surface 2D polynomial or
                        JSON file with data points for spline fit
    :param alias: string to select reconstruction type in the
                        coefficients table (if gainScale curve)
    :return calEner: numpy vector with calibrated energies
    """

    # locate coefficients in calibration table
    # ----------------------------------------
    coeffsDict = dict()
    with open(coeffsFile, "rt") as f:
        fileCont = f.read()   # JSON file
        if 'surface' in fileCont:
            ftype = "surface"
        elif fileCont[0] == '{':
            ftype = 'json'
        else:
            ftype = 'poly'

    if ftype == 'poly':
        codata = ascii.read(coeffsFile, guess=False, format='basic')
        # codata[1] : row 2
        # codata[1][1]: row 2, col 2
        print("Reading curve coefficients from", coeffsFile, "\n")
        for i in range(0, len(codata)):
            #  METHOD   ALIAS  a0  a1  a2  a3  a4
            coeffsDict[codata[i][1]] = (codata[i][2], codata[i][3],
                                        codata[i][4], codata[i][5],
                                        codata[i][6])
        npCoeffs = np.array(coeffsDict[alias])
    elif ftype == 'surface':
        print("Reading surface coefficients from", coeffsFile, "\n")
        npCoeffs = np.loadtxt(coeffsFile, comments="#")
    elif ftype == 'json':
        with open(coeffsFile, 'r') as f:
            jsonDict = json.load(f)
        nalias = len(jsonDict["ALIAS"])
        aliasidx = jsonDict["ALIAS"].index(alias)
        if len(jsonDict["xdata"]) % nalias == 0:
            nEnerCal = len(jsonDict["xdata"])//nalias
        else:
            raise ValueError("Length of xdata is not a multiple of number"
                             " of calibration energies")
        stridx = nEnerCal*aliasidx
        endidx = stridx + nEnerCal
        xdata = jsonDict["xdata"][stridx:endidx]
        ydata = jsonDict["ydata"][stridx:endidx]
        funinterp = interp1d(xdata, ydata, kind="linear",
                             fill_value="extrapolate")
    else:
        raise ValueError("Incorrect Coeffs file type")

    calEner = np.zeros(inEner.size, dtype=float)
    ie = 0

    # print("npCoeffs=",npCoeffs)
    #
    # convert energies
    #
    if ftype == 'surface':
        # print("Using surface to calibrate reconstructed energies\n")
        calEner = np.polynomial.polynomial.polyval2d(inEner, inPhase, npCoeffs)
    elif ftype == 'poly':
        for ie in range(0, inEner.size):
            # read fitting coeffs taken from polyfit2Bias.R (a0, a1, a2, a3)
            #  as in y = a0 + a1*x + a2*x^2 + a3*x^3
            # where y=E_reconstructed and x=Ecalibration (keV)

            # print("Using curve to calibrate reconstructed energies\n")
            # subtract "y" (0 = a0 + a1*x + a2*x^2 + ... - y) :
            npCoeffs = np.array(coeffsDict[alias])
            npCoeffs[0] -= inEner[ie]
            # reversed to say fit with poly1d definition:
            npCoeffsRev = npCoeffs[::-1]
            polyfit = np.poly1d(npCoeffsRev)
            # get real root (value of Ereal for a given Erecons )
            r = np.roots(polyfit)
            # real && >0 roots
            # print("r=", r)
            rreal = r.real[abs(r.imag) < 1e-5]
            # print("rreal=", rreal)
            rrealpos = rreal[rreal > 0]
            # print("rrealpos=", rrealpos)
            # closest root
            # print(inEner[ie])
            rclosest = min(enumerate(rrealpos),
                           key=lambda x: abs(x[1]-inEner[ie]))[1] # (idx,value)
            # print("For:", alias, " Recon energy=", rclosest)

            calEner[ie] = rclosest

    elif ftype == 'json':
        calEner = inEner/funinterp(inEner)

    # return calibrated energies
    # ------------------------------------
    return calEner


def VLtoFL(inputFile, extnum, outputFile):
    """
    Transform variable length column to fixed length using stilts
    Astropy & ftools & stilts does not work with variable-length arrays,
    so a previous conversion to fixed-format is required
    """

    print("\n Converting variable length column to fixed length column")
    print("Using:", inputFile, extnum, outputFile)
    try:
        # save inputFile into outFile
        copy(inputFile, outputFile)
        # comm = ("java -jar /home/sw/stilts.jar tpipe cmd='addcol " +
        #        "TIME2 \"(TIME*2)\"' cmd='delcols TIME2' in=" + outputFile +
        #        " out=" + outputFile)
        comm = ("stilts tpipe cmd='addcol " +
                "TIME2 \"(TIME*2)\"' cmd='delcols TIME2' ofmt=fits in=" +
                outputFile + "#" + str(extnum) + " out=" + outputFile)
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
        # copy header (modified by stilts) from one FITS header to another
        # existing FITS file header
        comm = ("cphead " + inputFile + "+" + str(extnum) + " " + outputFile +
                "+" + str(extnum))
        args = shlex.split(comm)
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error Coverting columns: ")
        print(comm)
        raise


class RxLines:
    def __init__(self, complabel, ilabels, energies_eV,
                 fwhms_eV, rel_amplitudes):
        self.complabel = complabel
        self.ilabels = ilabels
        self.energies_eV = energies_eV
        self.fwhms_eV = fwhms_eV
        self.rel_amplitudes = rel_amplitudes

    def getNumber(self):
        '''Get number of lines'''
        return len(self.ilabels)

    def plotLorentz(self):
        '''Ploting Lorentzian profiles for every line'''
        fig = plt.figure(figsize=(9, 6))
        ax1 = fig.add_subplot(1, 1, 1)
        minx = min(self.energies_eV)-20
        maxx = max(self.energies_eV)+20
        x_interval = np.linspace(minx, maxx, 1000)
        LmodSum = models.Lorentz1D(x_0=0, amplitude=0., fwhm=0.)
        nlines = self.getNumber()
        for i in range(nlines):
            Lmod = models.Lorentz1D(x_0=self.energies_eV[i], amplitude=self.rel_amplitudes[i],
                                    fwhm=self.fwhms_eV[i])
            LmodSum += Lmod
            ax1.plot(x_interval, Lmod(x_interval), label=self.ilabels[i])

        title = self.complabel + ' line complex'
        ax1.plot(x_interval, LmodSum(x_interval), marker='.', label=(self.complabel + ' complex'))
        ax1.set_xlabel("Energy (eV)")
        ax1.set_ylabel("Relative Intensity")
        ax1.set_title(title)
        ax1.legend()

    def broadGauss(self, sigma):
        '''Broad Lorentz lines profile with a Gaussian of given std_dev and plot results

        Parameters:
        sigma: standard deviation of the Gaussian to broaden the line profile

        '''
        fwhm = sigma*2*np.sqrt(2*np.log(2))
        fig = plt.figure(figsize=(9, 6))
        ax1 = fig.add_subplot(1, 1, 1)
        minx = min(self.energies_eV)-20
        maxx = max(self.energies_eV)+20
        x_interval = np.linspace(minx, maxx, 1000)
        LmodSum = models.Voigt1D(x_0=0, amplitude_L=0., fwhm_L=1., fwhm_G=1.)
        nlines = self.getNumber()
        for i in range(nlines):
            Lmod = models.Voigt1D(x_0=self.energies_eV[i], amplitude_L=self.rel_amplitudes[i],
                                  fwhm_L=self.fwhms_eV[i], fwhm_G=fwhm)
            LmodSum += Lmod
            ax1.plot(x_interval, Lmod(x_interval), label=self.ilabels[i])

        title = (self.complabel + ' Gauss-widened (sigma=' + str(sigma) + 'eV) line complex')
        ax1.plot(x_interval, LmodSum(x_interval), marker='.', label=(self.complabel + ' complex'))
        ax1.set_xlabel("Energy (eV)")
        ax1.set_ylabel("Relative Intensity")
        ax1.set_title(title)
        ax1.legend()


def fit2GaussAndRatio(data=None, a1=50, a2=90, mean1=5800, mean2=5900, sig1=5, sig2=5, nbins1=200, ratio=None,
                      xlab=None, xlim=(0, 0), ylim=(0, 0)):

    """"

    Fit 2 Gaussians (Ka1, Ka2) to Kas histogram
    Histograms are created and plotted with matplotlib.pyplot.hist in Density
    Gaussians functions from astropy.fitting module (fitting byLevMarLSQFitter)

    data1: (array) data for 1st histogram
    a1: (float) initial amplitude for Gaussian1
    a2: (float)initial amplitude for Gaussian2
    mean1: (float)initial mean for Gaussian1
    mean2: (float)initial mean for Gaussian2
    sig1: (float)std dev for Gaussian1
    sig2: (float)std dev for Gaussian2
    nbins1: number of bins for first (Kas) histogram
    ratio: GaussKa2/GaussKa1 ratio to select Ka2 photons
    xlab: xlabel of histogram plot
    xlim: (xmin,xmax) limits of X axis
    xlim: (ymin,ymax) limits of Y axis

    returns:
        PHmin,PHmax: x limiting values for Ka2 complex
    """

    fig = plt.figure(figsize=(16, 6))
    ax1 = fig.add_subplot(1, 2, 1)

    # create histogram
    bin_heights, bin_borders, _ = ax1.hist(data, bins=nbins1, density=True, label="Histogram", alpha=0.5)
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
    x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)

    # fit two gaussians to density-histogram (also "curve_fit" ?)
    gg_init = (models.Gaussian1D(amplitude=a1, mean=mean1, stddev=sig1) +
               models.Gaussian1D(amplitude=a2, mean=mean2, stddev=sig2))
    # fitter = fitting.SLSQPLSQFitter()
    fitter = fitting.LevMarLSQFitter()
    gg_fit = fitter(gg_init, bin_centers, bin_heights, maxiter=300)
    print("Message (Kas)=", fitter.fit_info['message'])

    # C1 = gg_fit.param_sets[0][0]
    # mean1 = gg_fit.param_sets[1][0] #u.a.
    # sigma1 = gg_fit.param_sets[2][0] #u.a.
    # C2 = gg_fit.param_sets[3][0]
    # mean2 = gg_fit.param_sets[4][0] #u.a.
    # sigma2 = gg_fit.param_sets[5][0] #u.a.

    C1 = gg_fit.amplitude_0[0]
    mean1 = gg_fit.mean_0[0]
    sigma1 = gg_fit.stddev_0[0]
    # fwhm1 = sigma1 * 2 * np.sqrt(2*np.log(2))
    C2 = gg_fit.amplitude_1[0]
    mean2 = gg_fit.mean_1[0]
    sigma2 = gg_fit.stddev_1[0]
    # fwhm2 = sigma2 * 2 * np.sqrt(2*np.log(2))

    # gg_fit.fit_info['residuals']
    g1 = models.Gaussian1D(amplitude=C1, mean=mean1, stddev=sigma1)
    g2 = models.Gaussian1D(amplitude=C2, mean=mean2, stddev=sigma2)
    ratioGG = g2(x_interval_for_fit)/g1(x_interval_for_fit)
    # print("minGG=",min(ratioGG), "maxGG=",max(ratioGG))

    # plot histogram and Gaussians fit
    ax1.plot(x_interval_for_fit, gg_fit(x_interval_for_fit), label='Gauss fit')
    ax1.plot(x_interval_for_fit, g1(x_interval_for_fit), label="Gauss Ka1")
    ax1.plot(x_interval_for_fit, g2(x_interval_for_fit), label="Gauss Ka2")
    ax1.plot(x_interval_for_fit, ratioGG, label="ratio G2/G1")
    ax1.set_xlabel(xlab)
    ax1.set_ylabel("Density")
    if (xlim[0] > 0 or xlim[1] > 0):
        ax1.set_xlim(xlim)
    if (ylim[0] > 0 or ylim[1] > 0):
        ax1.set_ylim(ylim)
    PHmin = x_interval_for_fit[ratioGG >= ratio][0]
    PHmax = min(x_interval_for_fit[ratioGG >= ratio][-1], (mean2+10*sigma2))
    ax1.axvline(PHmin, linestyle='--', color='tab:purple', label="Region for ratio")
    ax1.axvline(PHmax, linestyle='--', color='tab:purple')
    plt.legend()
    plt.show()

    print("Ka2 PHs in [" + '{:0.3f}'.format(PHmin) + "," + '{:0.3f}'.format(PHmax) + "] a.u.")
    return((PHmin, PHmax))


def fit3gauss2hist(data1=None, data2=None, a1=0.05, a2=0.09, a3=0.05, mean1=5800, mean2=5900, mean3=6500,
                   sig1=5, sig2=5, sig3=5, nbins1=200, nbins2=200, xlab=None, plot=True):
    """"

    Fit 2 Gaussians (Ka1, Ka2) to Kas histogram and 1 Gaussian to Kb histogram
    Histograms are created and plotted with matplotlib.pyplot.hist in Density
    Gaussians functions from astropy.fitting module (fitting byLevMarLSQFitter)

    data1: (array) data for 1st histogram
    data2: (array) data for 2st histogram
    a1: (float) initial amplitude for Gaussian1
    a2: (float)initial amplitude for Gaussian2
    a3: (float)initial amplitude for Gaussian3
    mean1: (float)initial mean for Gaussian1
    mean2: (float)initial mean for Gaussian2
    mean3: (float)initial mean for Gaussian3
    sig1: (float)std dev for Gaussian1
    sig2: (float)std dev for Gaussian2
    sig3: (float)std dev for Gaussian3
    nbins1: number of bins for first (Kas) histogram
    nbins2: number of bins for second (Kb) histogram
    xlab: xlabel of histogram plot
    plot: (bool) should histogram and fit be plotted?

    returns:
        (mean1, mean2, mean3): tuple of gaussians centres
    """
    fig = plt.figure(figsize=(20, 6))
    ax1 = fig.add_subplot(1, 2, 1)

    # create histogram
    bin_heights, bin_borders, _ = ax1.hist(data1, bins=nbins1, density=True, label="Histogram", alpha=0.5)
    if not plot:
        plt.clf()

    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
    x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)

    # fit two gaussians to density-histogram (also "curve_fit" ?)
    gg_init = (models.Gaussian1D(amplitude=a1, mean=mean1, stddev=sig1) +
               models.Gaussian1D(amplitude=a2, mean=mean2, stddev=sig2))
    # fitter = fitting.SLSQPLSQFitter()
    # fitter = fitting.LinearLSQFitter()
    #     -> model is not linear in parameters: linear fit should not be used
    fitter = fitting.LevMarLSQFitter()
    gg_fit = fitter(gg_init, bin_centers, bin_heights, maxiter=300)
    # check there are no errors in fitting
    if not fitter.fit_info['ierr'] in [1,2,3,4]:
        print("Solution not found for fit")
        print("Message (Kas) =", fitter.fit_info['message'])

    C1 = gg_fit.amplitude_0[0]
    mean1 = gg_fit.mean_0[0]
    sigma1 = gg_fit.stddev_0[0]
    fwhm1 = sigma1 * 2 * np.sqrt(2*np.log(2))
    C2 = gg_fit.amplitude_1[0]
    mean2 = gg_fit.mean_1[0]
    sigma2 = gg_fit.stddev_1[0]
    fwhm2 = sigma2 * 2 * np.sqrt(2*np.log(2))
    g1 = models.Gaussian1D(amplitude=C1, mean=mean1, stddev=sigma1)
    g2 = models.Gaussian1D(amplitude=C2, mean=mean2, stddev=sigma2)

    if plot:
        # plot 1st histogram and 2 Gaussians fit
        ax1.plot(x_interval_for_fit, gg_fit(x_interval_for_fit), label='Gauss fit')
        ax1.plot(x_interval_for_fit, g1(x_interval_for_fit), label="Gauss Ka1")
        ax1.plot(x_interval_for_fit, g2(x_interval_for_fit), label="Gauss Ka2")
        maxy = max(bin_heights)
        xtxt = min(data1)
        ax1.text(xtxt, maxy-0.03, "Double Gaussian", color='tab:orange')
        ax1.text(xtxt, maxy-0.04, ("Mean(Ka1)=" + '{:0.3f}'.format(mean1) + "a.u"), color='tab:green')
        ax1.text(xtxt, maxy-0.05, ("FWHM(Ka1)=" + '{:0.3f}'.format(fwhm1) + "ma.u"), color='tab:green')
        ax1.text(xtxt, maxy-0.06, ("Mean(Ka2)=" + '{:0.3f}'.format(mean2) + "a.u"), color='tab:red')
        ax1.text(xtxt, maxy-0.07, ("FWHM(Ka2)=" + '{:0.3f}'.format(fwhm2) + "ma.u"), color='tab:red')
        ax1.set_xlabel(xlab)
        ax1.set_ylabel("Density")
        ax1.legend()

    ax2 = fig.add_subplot(1, 2, 2)

    # create 2nd histogram
    bin_heights, bin_borders, _ = ax2.hist(data2, bins=nbins2, density=True, label="Histogram", alpha=0.5)
    if not plot:
        plt.clf()
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
    x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)
    # fit 1 gaussian to density-histogram (also "curve_fit" ?)
    g_init = models.Gaussian1D(amplitude=a3, mean=mean3, stddev=sig3)
    fitter = fitting.LevMarLSQFitter()
    g_fit = fitter(g_init, bin_centers, bin_heights)

    # check there are no errors in fitting
    if not fitter.fit_info['ierr'] in [1,2,3,4]:
        print("Solution not found for fit")
        print("Message (Kb) =", fitter.fit_info['message'])

    # C3 = g_fit.amplitude[0]
    mean3 = g_fit.mean[0]
    sigma3 = g_fit.stddev[0]
    fwhm3 = sigma3 * 2 * np.sqrt(2*np.log(2))

    if plot:
        # plot histogram and Gaussians fit
        ax2.plot(x_interval_for_fit, g_fit(x_interval_for_fit), label='Gauss fit')
        maxy = max(bin_heights)
        xtxt = min(data2)
        ax2.text(xtxt, maxy-0.03, ("Mean(Kb)=" + '{:0.3f}'.format(mean3) + "a.u"), color='tab:orange')
        ax2.text(xtxt, maxy-0.04, ("FWHM(Kb)=" + '{:0.3f}'.format(fwhm3) + "ma.u"), color='tab:orange')
        ax2.set_xlabel(xlab)
        ax2.set_ylabel("Density")
        ax2.legend()

    return((mean1, mean2, mean3))


def checkLine(lineComplex, lineComplex_fit):
    """ check consistency of theoretical and fitted line
        In particular:
            * ratio of intensities
            * ratio of line centres
            * widths of Lorentzians (fixed)
        Returns:
            1 : if problems found
            0 : if successful
    """
    # check intensity ratios
    nlines = lineComplex.getNumber()
    nlines_fit = lineComplex_fit.getNumber()
    if nlines != nlines_fit:
        print("Inconsistent number of lines:", nlines, nlines_fit)
        return 1
    for i in range(1, nlines):
        amp_0 = lineComplex.rel_amplitudes[0]
        amp_i = lineComplex.rel_amplitudes[i]
        ratioAmp = amp_i/amp_0
        amp_fit_0 = lineComplex_fit.rel_amplitudes[0]
        amp_fit_i = lineComplex_fit.rel_amplitudes[i]
        ratioAmp_fit = amp_fit_i/amp_fit_0
        if abs(ratioAmp-ratioAmp_fit)/ratioAmp > 1e-3:
            print("Inconsistent ratio of amplitudes for", lineComplex.ilabel[i], ":", ratioAmp, ratioAmp_fit)
            return 1

        x0 = lineComplex.energies_eV[0]
        xi = lineComplex.energies_eV[i]
        ratio_x0 = xi/x0
        x0_fit = lineComplex_fit.energies_eV[0]
        xi_fit = lineComplex_fit.energies_eV[i]
        ratio_x0_fit = xi_fit/x0_fit
        # print("Ratio of line centres for", lineComplex.ilabels[i], ":", ratio_x0, ratio_x0_fit)
        # print("Line:", lineComplex.energies_eV[i],lineComplex.energies_eV[0])
        # print("Line fit:", lineComplex_fit.energies_eV[i], lineComplex_fit.energies_eV[0])

        if abs(ratio_x0-ratio_x0_fit)/ratio_x0 > 1e-5:
            print("Inconsistent ratio of line centres for", lineComplex.ilabels[i], ":", ratio_x0, ratio_x0_fit)
            return 1

        if lineComplex.fwhms_eV[i] != lineComplex_fit.fwhms_eV[i]:
            print("Inconsistent fwhm_L for line", lineComplex.ilabels[i], ":",
                  lineComplex.fwhms_eV[i], lineComplex_fit.fwhms_eV[i])
            return 1
    return 0


def ecdf(data):
    """Compute the Empirical Cumulative Distribution Function
         or from statsmodels.distributions.empirical_distribution import ECDF
    """
    x = np.sort(data)
    n = x.size
    y = np.arange(1, n + 1) / n
    return x, y

def cdf_interp(model, minx, maxx, nsteps):
    """ Compute approx. (by cumsum) CDF for a model function over equispaced-data
        between minx and maxx over nsteps

        xdata: (array) x data to evaluate model
        model: (func or astropy model)
        returns:
            model_cdf_int: (func) interpolated & computed CDF

    """
    xdata = np.linspace(minx, maxx, nsteps)
    # calculate CDF of model fit
    bsize = xdata[1] - xdata[0]
    model_xdata = model(xdata)
    model05_xdata = model(xdata-0.5*bsize)
    # remove triangles in excess of the approx. CDF
    model_cdf  = np.cumsum(model_xdata)*bsize - model_xdata *bsize/2. - (model_xdata-model05_xdata)*bsize/8.
    model_cdf_int = interp1d(xdata, model_cdf)
    return model_cdf_int

def fitVoigt2hist(data=None, lines=None, nbins=None, outfig=None, ax=None):
    """
    Fit Voigt profiles (according to description in lines to data histogram) or
    compute curve of residuals vs. number of bins
    Histogram is created and plotted with matplotlib.pyplot.hist in Density
    Voigt functions from astropy.fitting module (fitting by LevMarLSQFitter)

    data: (array) data for 1st histogram
    lines: (RxLines) lines complex for histo
    nbins: (array) number of bins for histogram.
        If single value, perform fit and plot for this number of bins
        If array: do fit for each binning, calculate CDF for ech fit and
                  residuals. Plot comparative of CDF and residuals plot
    ax: (object) ax object from figure to plot results
    outfig: file to save final figure

    returns:
        fwhm_G: fwhm of Gaussian broadening (if run in single model (scalar nbins))

    Voigt relative intensities keep tied
    Lorentzian FWHMs are kept fixed
    Gaussian broadening are the same for all lines
    """

    # define functions and function-builders to tie parameters
    # Tie FWHM of Gaussians:
    def tie_gauss(model):
        return model.fwhm_G_0

    # Tie amplitudes of Voigt
    def tie_amplitude_builder(line, i):
        def tie_ampl(model):
            rel0 = line.rel_amplitudes[0]
            reli = line.rel_amplitudes[i]
            return model.amplitude_L_0 / rel0 * reli
        return tie_ampl

    # Tie line centres
    def tie_x0_builder(line, i):
        def tie_x0(model):
            return model.x_0_0 / line.energies_eV[0] * line.energies_eV[i]
        return tie_x0


    def fitVoigtLines(nbins=0, plotFit=False):
        """Fit Voigt profiles (according to description in lines to data histogram
        Histogram is created and plotted with matplotlib.pyplot.hist in Density
        Voigt functions from astropy.fitting module (fitting by LevMarLSQFitter)

        data: (array) data for 1st histogram
        lines: (RxLines) lines complex for histo
        nbins: number of bins for histogram.
        plotFit: if fit is to be plot

        returns:
            fitted model (astropy.fitting Voigt profiles object model)

        Voigt relative intensities keep tied
        Lorentzian FWHMs are kept fixed
        Gaussian broadening are the same for all lines
        """

        nlines = lines.getNumber()
        # create dictionary of functions for L_amplitude ties
        dict_tie_ampl = {}
        for i in range(0, nlines):
            func_name = "tie_ampl_" + str(i)
            dict_tie_ampl[func_name] = tie_amplitude_builder(lines, i)

        # create dict of functions for centres ties
        dict_tie_x0 = {}
        for i in range(0, nlines):
            func_name = "tie_x0_" + str(i)
            dict_tie_x0[func_name] = tie_x0_builder(lines, i)

        # define models (initial values, fixed params and ties)
        Vmods = list()
        Vmods.append(models.Voigt1D(x_0=lines.energies_eV[0], amplitude_L=lines.rel_amplitudes[0],
                                    fwhm_L=lines.fwhms_eV[0], fwhm_G=2., fixed={'fwhm_L': True}))
        #Vmods[0].fwhm_G.min = 1.e-6
        sumVoigt = Vmods[0]
        # print("Model=", sumVoigt)
        for i in range(1, nlines):
            Vmods.append(models.Voigt1D(x_0=lines.energies_eV[i], amplitude_L=lines.rel_amplitudes[i],
                                        fwhm_L=lines.fwhms_eV[i], fwhm_G=2., fixed={'fwhm_L': True}))
            Vmods[i].fwhm_G.tied = tie_gauss
            Vmods[i].amplitude_L.tied = dict_tie_ampl["tie_ampl_" + str(i)]
            Vmods[i].x_0.tied = dict_tie_x0["tie_x0_" + str(i)]

            sumVoigt += Vmods[i]
            # print("Adding Model=", Vmods[i])

        # create histogram
        bin_heights, bin_borders = np.histogram(data, bins=nbins,  density=True)
        if plotFit:
            ax.hist(data, bins=bin_borders, density=True, alpha=0.5, label="Histogram")

        bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
        x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)

        # fit nlines1 Voigts to density-histogram
        fitter = fitting.LevMarLSQFitter()
        vv_fit = fitter(sumVoigt, bin_centers, bin_heights, maxiter=300)

        # check there are no errors in fitting
        if not fitter.fit_info['ierr'] in [1,2,3,4]:
            print("Solution not found for nbins=", nbins)
            print("Message (", lines.complabel, ")=", fitter.fit_info['message'])

        # check fitting (ratios ,etc)
        ilabels_fit = [s + "_fit" for s in lines.ilabels]
        energies_eV_fit = np.zeros(nlines, dtype=np.float64)
        fwhms_eV_fit = np.zeros(nlines, dtype=np.float64)
        rel_amplitudes_fit = np.zeros(nlines, dtype=np.float64)

        # plot individual line and global model
        for i in range(nlines):
            # print("Line:", lines1.ilabels[i])
            energies_eV_fit[i] = vv_fit.param_sets[i*4][0]
            rel_amplitudes_fit[i] = vv_fit.param_sets[i*4+1][0]
            fwhms_eV_fit[i] = vv_fit.param_sets[i*4+2][0]
            fwhm_G = vv_fit.param_sets[i*4+3][0]
            # print("x0=", energies_eV_fit[i])
            # print("Ampl=", rel_amplitudes_fit[i])
            # print("FWHM_L=", fwhms_eV_fit[i])
            # print("FWHM_G=", fwhm_G)
            v = models.Voigt1D(x_0=energies_eV_fit[i], amplitude_L=rel_amplitudes_fit[i],
                           fwhm_L=fwhms_eV_fit[i], fwhm_G=fwhm_G)
            if plotFit:
                ax.plot(x_interval_for_fit, v(x_interval_for_fit), label="Voigt " + lines.ilabels[i])


        err_fwhm_G = 0.
        param_cov = fitter.fit_info['param_cov']
        if (param_cov is None):
            print("Matrix for ", lines.complabel, " is singular for nbins=", nbins)
        else:
            npars = param_cov.shape[0]
            err_fwhm_G = np.sqrt(param_cov[(npars-1), (npars-1)])

        # check fitting consistency
        lines_fit = RxLines(complabel=lines.complabel + "_fit", ilabels=ilabels_fit, energies_eV=energies_eV_fit,
                         fwhms_eV=fwhms_eV_fit, rel_amplitudes=rel_amplitudes_fit)
        status = checkLine(lines, lines_fit)
        if status:
            raise RuntimeError
            print("Line consistency check status:", status, "(0 = OK)")

        if plotFit:
            ax.plot(x_interval_for_fit, vv_fit(x_interval_for_fit), label='Voigt fit', color='black')
            maxy = max(bin_heights)
            xtxt = min(data)
            ax.text(xtxt, maxy-0.01, "FWHM_G=" + '{:0.2f}'.format(fwhm_G) + "+/-" + '{:0.2f}'.format(err_fwhm_G) + "eV")
            ax.set_xlabel("Energy (eV)")
            ax.set_ylabel("Density")
            ax.legend()

        return (fwhm_G, err_fwhm_G, vv_fit)


    if isinstance(nbins, int):  # if just a single value for nbins

        fwhm_G, err_fwhm_G, vv_model_fit = fitVoigtLines(nbins=nbins, plotFit=True)
        return (fwhm_G, err_fwhm_G)

    else:  # calculate residuals in diff of CDFs for different number of bins
        # calculate ECDF of data
        x_ecdf, y_ecdf = ecdf(data)
        residuals = np.zeros(len(nbins))

        for ib in range(0, len(nbins)):
            # calculate CDF of model fit
            fwhm_G, err_fwhm_G, vv_model_fit = fitVoigtLines(nbins=nbins[ib], plotFit=False)
            vv_model_cdf_int = cdf_interp(vv_model_fit, min(data)-30, max(data)+30, 10000)
            model_cdf = vv_model_cdf_int(x_ecdf)
            #ax.plot(x_ecdf, y_ecdf)
            #ax.plot(x_ecdf, model_cdf)
            residuals[ib] = np.sum((y_ecdf-model_cdf)**2)
        ax.plot(nbins, residuals)
        ax.set_xlabel("nbins")
        ax.set_ylabel("residuals: SUM(data_ECDF-model_CDF)²")


def gainScaleLinearFit(xData=None, yData=None, ylab="Lines energies (eV)"):
    """
    Fit a linear model to the input data (lines) and return slope and intercept

    xData: (array) x data
    yData: (array) y data
    return: slope, intercept
    """

    """
    # fit polynomial
    coefs = poly.polyfit(x=recon_lines, y=lines, deg=2)
    ffit = poly.polyval(x_interval_for_fit, coefs)
    """

    # define a model for a line; initialize a linear model
    line_init = models.Linear1D()
    # initialize a linear fitter
    fitter = fitting.LinearLSQFitter()
    # fit the data with the fitter
    fitted_line = fitter(line_init, xData, yData)
    print(fitted_line)
    print("Residuals=", fitter.fit_info['residuals'])
    print("Params=", fitter.fit_info['params'])
    slope = fitted_line.slope[0]
    inter = fitted_line.intercept[0]
    # check R²
    absError = fitted_line(xData) - yData
    SE = np.square(absError)  # squared errors
    MSE = np.mean(SE)  # mean squared errors
    RMSE = np.sqrt(MSE)  # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (np.var(absError) / np.var(yData))
    print('RMSE:', RMSE)
    print('R-squared:', Rsquared)

    # plot the model
    fig = plt.figure(figsize=(16, 6))
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.plot(xData, yData, 'ko', label='Data')
    ax1.plot(xData, fitted_line(xData), 'k-', label='Fitted Model', color="red")
    ax1.set_xlabel('reconstructed lines (a.u.)')
    ax1.set_ylabel(ylab)
    ax1.set_xlim(xData[0]-5, xData[1]+5)
    ax1.set_ylim(yData[0]-5, yData[1]+5)
    ax1.legend()
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.plot(xData, yData, 'ko', label='Data')
    ax2.plot(xData, fitted_line(xData), 'k-', label='Fitted Model', color="red")
    ax2.set_xlabel('reconstructed lines (a.u.)')
    ax2.set_ylabel(ylab)
    ax2.legend()

    return (slope, inter)


def jitterCorr(reconPH=None, phase=None):

    """
    Calculate jitter correction: dependance of recons PH vs Phase
    Fit a polynomial (deg 2) and subtract effect
    Return: reconstructed PH with jitter removed

    reconPH: (array) reconstructed PH
    phase1: (array) phase information

    """

    # plot phases
    fig = plt.figure(figsize=(20,6))
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.scatter(phase, reconPH, marker="o")
    # fit polynomial
    x_interval_for_fit = np.linspace(min(phase), max(phase), 10000)
    ##poly1 = np.poly1d(np.polyfit(phaseKas_HR, enerKas_HR, 2))
    # poly.polyfit recommended instead of np.polyfit + np.poly1d
    coefs = poly.polyfit(x=phase, y=reconPH, deg=2)
    ffit = poly.polyval(x_interval_for_fit, coefs)
    ax1.plot(x_interval_for_fit, ffit,'-', color="red")
    ax1.set_xlabel("Phase (samples)")
    ax1.set_ylabel("PH Kas (a.u.)")
    print("Fit Kas=",'{:0.3f}'.format(coefs[0]) + "+ (" + '{:0.3f}'.format(coefs[1]) + ")*x" +
          "+(" + '{:0.3f}'.format(coefs[2]) + ")*x²" )
    # subtract polynomial (flat jitter effect)
    ax2 = fig.add_subplot(1, 2, 2)
    reconPH_jitter = reconPH - coefs[1]*phase - coefs[2]*phase**2
    ax2.scatter(phase, reconPH_jitter, marker="o")
    ax2.set_xlabel("Phase (samples)")
    ax2.set_ylabel("Corrected PH Kas (a.u.)")
    coefsJ = poly.polyfit(x=phase, y=reconPH_jitter, deg=2)
    print("Fit Kas (corrected)=",'{:0.3f}'.format(coefsJ[0]) + "+ (" + '{:0.3f}'.format(coefsJ[1]) +
          ")*x" + "+(" + '{:0.3f}'.format(coefsJ[2]) + ")*x²" )
    ffitJ = poly.polyval(x_interval_for_fit, coefsJ)
    ax2.plot(x_interval_for_fit, ffitJ,'-', color="white")
    return(reconPH_jitter)


def getMaximaDensity(data=None, nmaxima=0, intervalmin=0., intervalmax=0.):
    """ Get x coordinate(s) for maxima of density plot in histogram
        Kernel Density Curve obtained by seaborn.distplot

    data: (array) data in histogram
    nmaxima: number of maxima to locate
    intervalmin.intervalmax:
        (float) limits of range in x where maxima must be enclosed
    return
        centres: (numpy array) with maxima location
    """
    # For a complete plot:
    # fig = plt.figure(figsize=(20, 6))
    # ax1 = fig.add_subplot(1, 2, 1)
    # sns.distplot(1e3*dataKas_HR.SIGNAL, hist=True, kde=True,
    #              bins=nbinsKas, color = 'darkblue',
    #              hist_kws={'edgecolor':'black'},
    #              kde_kws={'linewidth': 3, "label": "KDE"})
    x,y = sns.distplot(data).get_lines()[0].get_data()
    plt.clf() # clear all Axes
    xlines = x[(x>intervalmin) & (x<intervalmax)]
    ylines = y[(x>intervalmin) & (x<intervalmax)]
    ymaxima = argrelextrema(ylines, np.greater)
    centres = xlines[ymaxima]
    # [ax1.axvline(c, color="gray", linestyle='--') for c in centres]
    if len(centres) > nmaxima:
        print("Too many maxima locations:", centres)
        raise RuntimeError
    return centres
