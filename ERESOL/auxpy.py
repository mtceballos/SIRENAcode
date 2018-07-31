from __future__ import print_function
import os
import sys
import shlex
import shutil
import tempfile
from time import gmtime, strftime
import numpy as np
import math
from subprocess import check_call, check_output,STDOUT
from astropy.io import fits, ascii
from multiprocessing import Pool
from itertools import repeat
from glob import glob
import xml.etree.ElementTree as ET


tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"


def tsfn(i, nproc, comm, simTimeN):
    """
    Funcion to run a single process of noise generation
    """
    if i == 0:
        tstart = 0.
    else:
        tstart = i * simTimeN/10. + 0.01
    tend = simTimeN/10. * (1 + i)
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


def simulNoise(pixName, samprate, jitter, stoch, bbfb, pulseLength,
               space, acbias, scaleFactor, samplesUp, nSgms,
               nintervals, simTimeN, pixel, preBufferSize):
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
          :param simTimeN : Simulation time (s) for noise spectra calculation
          :param pixel : Pixel number
          :return fits file with Noise spectra in specified space and sampling

    """
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
    if jitter == "jitter":
        stochastic = ""
    if stoch == "stoch":
        stochStr = "_stoch"
        stochastic = (" stochastic_integrator=y")
    # bbfb
    bbfbStr = ""
    dobbfb = " dobbfb=n"
    if jitter == "jitter":
        dobbfb = ""
    if bbfb == "bbfb":
        bbfbStr = "_bbfb"
        dobbfb = (" dobbfb=y carrier_frequency=2e6 bbfb_delay=40" +
                  " decimation_filter=y")

    # ----GLOBAL VARIABLES -------------
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA"
    simSIXTEdir = EURECAdir + "/testHarness/simulations/SIXTE"
    XMLdir = EURECAdir + "/ERESOL"
    XMLfile = (XMLdir + "/" + "xifu_detector_hex_baselineNEWgrades" +
               smprtStr + ".xml")
    XMLtree = ET.parse(XMLfile)
    XMLroot = XMLtree.getroot()
    for samplefreq in XMLroot.findall('samplefreq'):
        samplingrate = samplefreq.get('value')
    cwd = os.getcwd()

    # set input params dependent variables
    # -------------------------------------
    triggerSize = max(10000, pulseLength + preBufferSize)
    # PixTypeFile = "'file:" + simSIXTEdir + "/tespixels.fits[" + array + "]'"
    PixTypeFile = ("'file:" + simSIXTEdir +
                   "/newpix_full.fits[" + pixName + "]'")
    tessim = "tessim" + pixName
    wdir = simSIXTEdir + "/NOISE/" + tessim
    os.chdir(wdir)
    # define files
    rootN = ("forNoise" + str(pulseLength) + "samples_" + tessim +
             "_" + str(simTimeN) + "s_pairscps_" + space)
    noiseFile = ("noise" + str(pulseLength) + "samples_" + tessim + "_B0_" +
                 space + smprtStr + jitterStr + stochStr + bbfbStr + ".fits")
    pixFileN = rootN + smprtStr + jitterStr + ".piximpact"
    fitsFileN = rootN + smprtStr + jitterStr + stochStr + bbfbStr + ".fits"

    # -------------------------------------------------------------------------
    # get NOISE file: process empty trigger file
    # -------------------------------------------------------------------------
    if not os.path.isfile(pixFileN):    
        print("\nTESCONSPILEUP: Generating no-events file for NOISE...")
        comm = ("tesconstpileup PixImpList=" + pixFileN + " XMLFile=" + XMLfile
                + " tstop=" + str(simTimeN) + offset +
                " energy=0 pulseDistance=1 TriggerSize=" +
                str(triggerSize) + " clobber=yes")
        args = shlex.split(comm)
        try:
            check_call(args, stderr=STDOUT)
        except RuntimeError:
            print("Error running tool tesconstpileup (NOISE):")
            print(comm)
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
            raise
        print("TESCONSPILEUP: ...................................END")

    if not os.path.isfile(fitsFileN):        
        print("\nTESSIM: Generating stream file for NOISE...")
        tstart = 0.
        tstop = simTimeN

        if bbfb == "bbfb":
            # in parallel if bbfb
            print("Running in parallel for BBFB...")
            comm = ("tessim PixID=1 PixImpList=" + pixFileN +
                    " sample_rate=" + samplingrate +
                    " triggertype=noise triggersize=10000 " +
                    " prebuffer=1000 acbias=yes" + stochastic + dobbfb +
                    " clobber=yes PixType='file:" +
                    simSIXTEdir + "/newpix_full.fits[LPA2shunt]'")
            nproc = 10

            with Pool(processes=nproc) as pool:
                pool.starmap(tsfn, zip(range(nproc), repeat(nproc),
                                       repeat(comm), repeat(simTimeN)))
                pool.close()
                pool.join()

                # merge files
                formerge = "fmerge '"
                for i in range(nproc):
                    tmpFile = "tmpbbfb" + str(i) + ".fits "
                    formerge += tmpFile
                    formerge += "'"
                    comm = (formerge + " outfile=" + fitsFileN +
                            " columns='-' clobber=yes")
                try:
                    args = shlex.split(comm)
                    check_call(args, stderr=STDOUT)
                except RuntimeError:
                    print("Error merging tessim noises")
                    os.chdir(cwd)
                    shutil.rmtree(tmpDir)
                    raise
                updateHISTORY(fitsFileN, comm)
                for sf in glob("tmpbbfb*.fits"):
                    os.remove(sf)
        else:
            comm = ("tessim PixID=" + str(pixel) + " PixImpList=" + pixFileN +
                    " Streamfile=" + fitsFileN + " sample_rate=" + samplingrate +
                    " tstart=" + str(tstart) + " tstop=" + str(tstop) +
                    " triggertype=noise " + " triggersize=" + str(triggerSize) +
                    " prebuffer=" + str(preBufferSize) + " acbias=" + acbias +
                    stochastic + dobbfb + " clobber=yes PixType=" + PixTypeFile)

            args = shlex.split(comm)
            try:
                print("Running tool tessim (NOISE):", comm)
                check_call(args, stderr=STDOUT)
            except RuntimeError:
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
                rmLastAndFirst(fitsFileN, 0)
                updateHISTORY(fitsFileN, comm)
            except RuntimeError:
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
                    RfromItrigger(infile=fitsFileN, Icol="ADC",
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
                except RuntimeError:
                    print("Error adding resistance column in ADCR -> ADC (NOISE):")
                    print(comm)
                    os.chdir(cwd)
                    shutil.rmtree(tmpDir)
                    raise
                print("FCALC: ..................................END")

    print("\nGENNOISESPEC: Generating NOISE spectrum file "
          "in (", space, " space)")
    comm = ("gennoisespec --inFile=" + fitsFileN +
            " --outFile=" + noiseFile +
            " --intervalMinSamples=" + str(pulseLength) +
            " --nintervals=" + str(nintervals) +
            " --scaleFactor=" + str(scaleFactor) +
            " --samplesUp=" + str(samplesUp) +
            " --nSgms=" + str(nSgms) +
            " --pulse_length=" + str(pulseLength) +
            " --clobber=yes verbosity=0 --weightMS=no ")
    args = shlex.split(comm)
    try:
        check_call(args, stderr=STDOUT)
    except RuntimeError:
        print("Error running tool gennoise (NOISE):")
        print(comm)
        os.chdir(cwd)
        shutil.rmtree(tmpDir)
        raise
    print("GENNOISESPEC: .......................................END")
    print("Noise file ", noiseFile, " has been successfully created\n")
    os.chdir(cwd)


def simulSingles(pixName, monoEkeV, acbias, samprate, jitter, noise, stoch,
                 bbfb, nSimPulses, singleSeparation, pixel, preBufferSize,
                 pulseLength):
    """
    :param pixName: Extension name in the FITS pixel definition file
                    (SPA*, LPA1*, LPA2*, LPA3*)
    :param monoEkeV: Monochromatic energy (keV) of input simulated pulses
    :param acbias: Operating Current (AC if acbias=yes or DC if acbias=no)
    :param samprate: Samprate value with respect to baseline of 156250 Hz:
                    "" (baseline), "samprate2" (half_baseline)
    :param jitter: jitter option ("" for no_jitter and "jitter" for jitter)
    :param noise: noise option ("" for noise and "nonoise" for nonoise)
    :param stoch: stochastic option ("" for non-stochastic and
                                     "stoch" for stochastic)
    :param bbfb: ("") for dobbfb=n or ("bbfb") for dobbfb=y
    :param nSimPulses: number of pulses to be simulated
    :param singleSeparation: separation from secondary->primary for next record
    :param pixel: pixel number
    :param preBufferSize: pre-buffersize
    :param pulseLength: length of pulses to calculate trigger size
    :return: files with simulated SINGLES
    """

    cwd = os.getcwd()
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA"
    simSIXTEdir = EURECAdir + "/testHarness/simulations/SIXTE"
    ERESOLdir = EURECAdir + "/ERESOL"
    PAIRSdir = ERESOLdir + "/PAIRS"

    # increase threshold so that they are all triggered @999 (if no jitter...)
    triggerTH = {'LPA1shunt': 50, 'LPA2shunt': 20}
    if monoEkeV == "0.5":
        triggerTH["LPA2shunt"] = 50

    # samprate
    smprtStr = ""
    if samprate == 'samprate2':
        smprtStr = "_samprate2"
        # increase threshold so that they are all triggered @999 (if !jitter..)
        triggerTH["LPA2shunt"] = 25
        # increase threshold so that they are all triggered @999 (if !jitter..)
        if monoEkeV == "0.5":
            triggerTH["LPA2shunt"] = 60

    XMLfile = ERESOLdir + "/xifu_detector_hex_baselineNEWgrades" + \
        smprtStr + ".xml"
    XMLtree = ET.parse(XMLfile)
    XMLroot = XMLtree.getroot()
    for samplefreq in XMLroot.findall('samplefreq'):
        samprate = samplefreq.get('value')
    # added to solve floating point inaccu. due to sampling rate
    # (Christian's mail 31/03/2017)
    tstart = 0.5 / float(samprate)

    # jitter
    jitterStr = ""
    offset = ""
    if jitter == "jitter":
        jitterStr = "_jitter"
        offset = " offset=-1"

    # noise
    noiseStr = ""
    simnoise = " simnoise=y"
    if noise == "nonoise":
        noiseStr = "_nonoise"
        simnoise = " simnoise=n"

    # stochastic
    stochStr = ""
    stochastic = " stochastic_integrator=n"
    if jitter == "jitter":
        stochastic = ""
    if stoch == "stoch":
        stochStr = "_stoch"
        stochastic = (" stochastic_integrator=y")
    # bbfb
    bbfbStr = ""
    dobbfb = " dobbfb=n"
    if jitter == "jitter":
        dobbfb = ""
    if bbfb == "bbfb":
        bbfbStr = "_bbfb"
        dobbfb = (" dobbfb=y carrier_frequency=2e6 bbfb_delay=40" +
                  " decimation_filter=y")

    tessim = "tessim" + pixName
    SIMFILESdir = PAIRSdir + "/" + tessim
    # pre-Feb 2018 telecon:
    # PixTypeFile = "file:" + simSIXTEdir + "/newpixels.fits[" + pixName + "]"
    PixTypeFile = "file:" + simSIXTEdir + "/newpix_full.fits[" + pixName + "]"

    triggerSizeTC = preBufferSize + singleSeparation + singleSeparation + 1000
    triggerSizeTS = preBufferSize + pulseLength + 1000
    triggerTS3val = triggerSizeTS - preBufferSize
    triggerSizeTC = int(triggerSizeTC)

    # calculate sim time to have at least nSimPulses pulses:
    #   simTime= (nSimPulses/2)recs * (triggerSize sam/rec)/(samprate sam/s)
    simTime = nSimPulses/2. * triggerSizeTC/float(samprate)
    simTime = '{0:0.0f}'.format(simTime)

    # for piximpact:
    root0 = "sep" + str(singleSeparation) + "sam_" + simTime + "s_" + \
        str(monoEkeV) + "keV" + smprtStr + jitterStr
    # for fits:
    root = ("sep" + str(singleSeparation) + "sam_" + str(nSimPulses) + "p_" +
            str(monoEkeV) + "keV" + smprtStr + jitterStr + noiseStr +
            stochStr + bbfbStr)
    pixFile = (PAIRSdir + "/PIXIMPACT/" + root0 + "_trSz" +
               str(triggerSizeTC) + ".piximpact")
    fitsFile = SIMFILESdir + "/" + root + ".fits"
    print("-------------------------------------------\n")
    print("Simulating ", fitsFile, "\n")
    print("-------------------------------------------\n")

    if not os.path.isfile(pixFile):
        comm = ("tesconstpileup PixImpList=" + pixFile +
                " XMLFile=" + XMLfile + " timezero=" + str(tstart) +
                " tstop=" + str(simTime) + " energy=" + str(monoEkeV) +
                offset + " pulseDistance=" + str(singleSeparation) +
                " TriggerSize=" + str(triggerSizeTC) + " clobber=yes")

        print("\n##### Runing tesconstpileup #########")
        print(comm, "\n")
        try:
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except RuntimeError:
            print("Error running tool for piximpact list generation")
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
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
        commTessim = ("tessim PixID=" + str(pixel) +
                      " PixImpList=" + pixFile +
                      " Streamfile=" + fitsFile +
                      " tstart=0." + " tstop=" + simTime +
                      " triggerSize=" + str(triggerSizeTS) +
                      " preBuffer=" + str(preBufferSize) +
                      " triggertype='diff:3:" + str(triggerTH[pixName]) +
                      ":" + str(triggerTS3val) + "'" +
                      " acbias=" + acbias + " sample_rate=" + samprate +
                      simnoise + stochastic + dobbfb +
                      " PixType=" + PixTypeFile)

        print("\n##### Runing tessim #########")
        print(commTessim, "\n")
        try:
            args = shlex.split(commTessim)
            check_call(args, stderr=STDOUT)
        except RuntimeError:
            print("Error running TESSIM for data simulation")
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
            raise
        # continue
        # update HISTORY  for tessim run in header[0]
        updateHISTORY(fitsFile, commTessim)

        # rm first (and LAST) record and update NETTOT
        fsim = fits.open(fitsFile)
        nrows = fsim[1].header["NAXIS2"]
        assert nrows > 1, "Tessim failed: just one huge row present!"
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
            shutil.rmtree(tmpDir)
            raise
    os.chdir(cwd)


def simulLibsGlobal(pixName, space, samprate, jitter, noise, stoch, bbfb,
                    pulseLength, libEnergies, largeFilter, nsamples,
                    nSimPulses, acbias, tstartPulse1All,
                    createLib, noiseMat, weightMat, preBufferSize,
                    separation, pixel):
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
    :type stoch: str
    :param stoch: ("stoch") stochastic or ("") non-stochastic simulations
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
    :type weigthMat: str
    :param weightMat: should the Weight matrices HDU be created? (yes/no)
    :param separation: separation from secondary->primary for next record
    :param pixel: pixel number
    :param preBufferSize: pre-buffersize
    :return: simulated calibration pulses pulses && Global library from them

    """

    cwd = os.getcwd()
    tessim = "tessim" + pixName
    triggerTH = {'LPA1shunt': 50, 'LPA2shunt': 20}

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

    # noise
    noiseStr = ""
    simnoise = " simnoise=y"
    if noise == "nonoise":
        noiseStr = "_nonoise"
        simnoise = " simnoise=n"

    # stochastic
    stochStr = ""
    stochastic = " stochastic_integrator=n"
    if jitter == "jitter":
        stochastic = ""
    if stoch == "stoch":
        stochStr = "_stoch"
        stochastic = (" stochastic_integrator=y")
    # bbfb
    bbfbStr = ""
    dobbfb = " dobbfb=n"
    if jitter == "jitter":
        dobbfb = ""
    if bbfb == "bbfb":
        bbfbStr = "_bbfb"
        dobbfb = (" dobbfb=y carrier_frequency=2e6 bbfb_delay=40" +
                  " decimation_filter=y")

    # -- CALIB/SIM/NOISE/LIB/XML dirs & files --
    EURECAdir = "/dataj6/ceballos/INSTRUMEN/EURECA"
    simSIXTEdir = EURECAdir + "/testHarness/simulations/SIXTE"
    XMLdir = EURECAdir + "/ERESOL"
    PIXIMPACTdir = simSIXTEdir + "/LIBRARIES/PIXIMPACT"
    SIMFILESdir = simSIXTEdir + "/LIBRARIES/" + tessim
    NOISEdir = simSIXTEdir + "/NOISE/" + tessim
    libDir = SIMFILESdir + "/GLOBAL/" + space

    # XMLfile = XMLdir + "/" + "xifu_baseline.xml"
    XMLfile = (XMLdir + "/" + "xifu_detector_hex_baselineNEWgrades" +
               smprtStr + ".xml")
    XMLtree = ET.parse(XMLfile)
    XMLroot = XMLtree.getroot()

    noiseFile = (NOISEdir + "/noise" + str(nsamples) + "samples_" + tessim +
                 "_B0_" + space + smprtStr + jitterStr + stochStr + bbfbStr +
                 ".fits")
    # PixTypeFile = "'file:"+ simSIXTEdir + "/newpixels.fits[" + pixName + "]'"
    PixTypeFile = ("'file:" + simSIXTEdir +
                   "/newpix_full.fits[" + pixName + "]'")
    libFile = (libDir + "/libraryMultiE_GLOBAL_PL" + str(pulseLength) +
               "_" + str(nSimPulses) + "p" + smprtStr + jitterStr +
               noiseStr + stochStr + bbfbStr + ".fits")
    evttmpFile = tempfile.NamedTemporaryFile()

    for samplefreq in XMLroot.findall('samplefreq'):
        samplingrate = samplefreq.get('value')

    # added to solve floating point inaccuracies due to sampling rate
    # (Christian's mail 31/03/2017):
    tstart = 0.5/float(samplingrate)

    # Calibration energies and Tstarts of pulses
    tstartPulse1 = dict(zip(libEnergies, tstartPulse1All))

    # calculate simTime so that nSimPulses are simulated (1 pulses per record)
    simTime = nSimPulses * (preBufferSize + separation) / float(samplingrate)
    tstop = tstart + simTime
    simTime = '{0:0.0f}'.format(simTime)
    tstop = '{0:0.0f}'.format(tstop)

    # Sigmas and Samples and scaleFactor for Detection
    samplesUp = 3
    samplesDown = 3
    nSgms = 5
    maxFilterLength = pulseLength
    lF = ""
    if largeFilter > 0:
        lF = " largeFilter= " + str(largeFilter)
        maxFilterLength = largeFilter

    # Trigger sizes in tesconstpileup
    triggerSizeTC = preBufferSize + separation + separation + 1000
    # and tessim ___|\_______________ooo
    triggerSizeTS = preBufferSize + maxFilterLength + 1000
    triggerTS3val = triggerSizeTS - preBufferSize

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
        if monoEkeV == "0.5":
            triggerTH["LPA2shunt"] = 50

        print("=============================================")
        print("Adding monochromatic energy", monoEkeV, "keV")
        print("=============================================")
        print("simTime=", simTime, "\n")
        # print("Temporary event file in: ", evttmpFile.name)
        # simulate SIXTE file with tesconstpileup + tessim
        root0 = ("mono" + str(monoEkeV) + "_sep" + str(separation) + "_pix" +
                 str(pixel) + "_" + str(nSimPulses) + "p")
        root = root0 + "_" + str(pulseLength)
        # pixFile = PIXIMPACTdir +
        #          "/" + root0 + "_trSz" + str(triggerSizeTC) + ".piximpact"
        pixFile = (PIXIMPACTdir + "/" + root0 + smprtStr + jitterStr +
                   ".piximpact")
        simFile = (SIMFILESdir + "/" + root + smprtStr + jitterStr + noiseStr +
                   stochStr + bbfbStr + ".fits")

        # -- TESCONSTPILEUP: generate impacts for well separated single pulses
        # -- TESSIM: simulate well separated pulses --
        if not os.path.isfile(simFile):
            print("Simulating & triggering pulses with TESSIM to ", simFile)
            comm = ("tesconstpileup PixImpList=" + pixFile +
                    " XMLFile=" + XMLfile + " timezero=" + str(tstart) +
                    " tstop=" + str(tstop) + " energy=" + str(monoEkeV) +
                    " pulseDistance=" + str(separation) + offset +
                    " TriggerSize=" + str(triggerSizeTC) + " clobber=yes")
            print("Generating PIXIMPACT ", pixFile, " running:\n", comm)
            try:
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except RuntimeError:
                print("Error running task for piximpact list generation")
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise
                
            # continue # to only generate piximpact
            commTessim = ("tessim PixID=" + str(pixel) +
                          " PixImpList=" + pixFile +
                          " Streamfile=" + simFile +
                          " tstart=0 tstop=" + str(simTime) +
                          " triggerSize=" + str(triggerSizeTS) +
                          " preBuffer=" + str(preBufferSize) +
                          " acbias=" + acbias +
                          " triggertype='diff:3:" + str(triggerTH[pixName]) +
                          ":" + str(triggerTS3val) + "'" +
                          " sample_rate=" + samplingrate + dobbfb +
                          simnoise + stochastic + " PixType=" + PixTypeFile)

            print("Running tessim for simulation\n", commTessim)
            try:
                args = shlex.split(commTessim)
                check_call(args, stderr=STDOUT)
            except RuntimeError:
                print("Error running TESSIM for simulation\n")
                os.chdir(cwd)
                shutil.rmtree(tmpDir)
                raise

            rmLastAndFirst(simFile, 1)

            # update HISTORY in header[0]
            dateTime = strftime("%Y-%m-%d %H:%M:%S", gmtime())
            history = ("Created & Updated by simulLibsGlobal.py on " +
                       dateTime + "with command: " + commTessim)
            updateHISTORY(simFile, history)

        # print("Antes de evaluar createLib\n")
        if not createLib:
            print("NO Va a crear la libreria\n")
            continue  # (to just simulate pulses files)

        print("Va a crear la libreria\n")
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
                " mode=0 clobber=yes intermediate=0" +
                " monoenergy=" + str(monoEeV) + " EventListSize=1000" +
                " XMLFile=" + XMLfile +
                " hduPRECALWN=" + weightMat + " hduPRCLOFWM=" + noiseMat)

        args = shlex.split(comm)
        print("SIRENA reconstruction to add a line to the library, "
              "running command:\n", comm)
        try:
            check_call(args, stderr=STDOUT)
        except RuntimeError:
            print("Error running SIRENA to add new line to library")
            os.chdir(cwd)
            shutil.rmtree(tmpDir)
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
             5: phi (distance to central parabola)
             6: phi (distance to central parabola) + n (additional lags). Ex.:
                *    | *      *   => phi=-0.3, n=0
                *      * |    *   => phi= 0.2, n=0
                *      *      *    | *  => phi=0.8, n=1
        :return arrPhase:
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
    #print("impF=", impF)
    #print("imp0F=", imp0F)
    #print("evtF=", evtF)
    if os.path.isfile(imp0F):
        imp0 = fits.open(imp0F, memmap=True)
        imp0Tab = imp0[1].data
        imp0Times = imp0Tab['TIME']
        imp0.close()
        arrPhase[0] = (impTimes - imp0Times)*samprate
        clsidx0 = find_nearest_idx(imp0Times, evtTimes)
        # arrPhase 2 should be the same that arrPhase1
        arrPhase[2] = (evtTimes - imp0Times[clsidx0]) * samprate
        #print("max(Arr1-Arr2)=", max(abs(arrPhase[1]-arrPhase[2])))
        #print("Arr1=", arrPhase[1])
        #print("Arr2=", arrPhase[2])
        assert np.allclose(arrPhase[1], arrPhase[2], atol=5E-6), \
            "Incompatible arrival Phases calculations"
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
    arrPhase[5] = evtPHI
    arrPhase[6] = evtPHI + evtN
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
    #     shutil.rmtree(tmpDir)
    #     raise
    # try:
    #     args = shlex.split(getNumwrns)
    #     numwrnsStr = check_output(args, stderr=STDOUT)
    #     print("numwrnsStr=", numwrnsStr)
    #     numwrns = int(numwrnsStr)
    #     print("numwrns=", numwrns)
    # except:
    #     print("Error running fverify(2) to get numwrns")
    #     shutil.rmtree(tmpDir)
    #     raise
    # if os.path.isfile("pp.stdout"):
    #    os.remove("pp.stdout")
    # if fverify detects errors/warnings, return numerrs+numwrns
    # return numerrs+numwrns


def rmLastAndFirst(simfile, ppr):
    """
    rm first (and LAST) record of SIXTE/tessim simulated file and update NETTOT
    (first starts high and last can be cut) (see Christian's email from
    19 Jan 2017 @ EURECA). Also update number of pulses in NETTOT keyword

    :type simfile: str
    :param simfile: simulated file where cleaning must be done
    :type ppr: int
    :param ppr: pulses per record in simulations
    """

    fsim = fits.open(simfile, mode='update')
    # nrows = fsim[1].header["NAXIS2"]
    nrows = fsim[1].data.shape[0]

    # remove first and last row with astropy: it does not read ADC properly
    # binTable = fsim[1].data
    # binTable = binTable[2:nrows-1]
    # fsim.flush()
    fsim.close()

    assert nrows > 1, "Tessim failed for (%s): just one row present " % simfile
    tmpFile = tmpDir + "/pp.fits"
    shutil.copy(simfile, tmpFile)

    # try to remove last row
    comm1 = "fdelrow infile=" + simfile + "+1 firstrow=" + str(nrows) +\
        " nrows=1 confirm=no proceed=yes"
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
    nrows2 = fsim[1].header['NAXIS2']
    nettot = nrows2 * ppr  # new number of pulses (==2*nofrecords)
    fsim[1].header['NETTOT'] = nettot
    fsim.close()

    # remove first row
    if errs1 > 0:
        shutil.copy(tmpFile, simfile)
        print("Errors running FTOOLS to remove final row (abandon) in ",
              simfile, " (returned to original file)")
    elif warns1 > 0:
        print("Warnings running FTOOLS to remove final row (abandon) in ",
              simfile, " (keep modified file)")
    else:
        comm2 = ("fdelrow infile=" + simfile +
                 "+1 firstrow=1 nrows=1 confirm=no proceed=yes")
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
        nrows2 = fsim[1].header['NAXIS2']
        nettot = nrows2 * ppr  # new number of pulses (==2*nofrecords)
        fsim[1].header['NETTOT'] = nettot
        fsim.close()

        if sum(fitsVerify(simfile)) > 0:
            shutil.copy(tmpFile, simfile)
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


def enerToCalEner(inEner, coeffsFile, alias):
    """
    :param inEner: numpy array with input uncorrected energies
    :param coeffsFile: file with coefficients of polynomial fit to gain curves
                        from polyfit2bias.R
    :param alias: string to select reconstruction type in the
                        coefficients table
    :return calEner: numpy vector with calibrated energies
    """

    # locate coefficients in calibration table
    # ----------------------------------------
    coeffsDict = dict()
    if coeffsFile:
        codata = ascii.read(coeffsFile, guess=False, format='basic')
        # codata[1] : row 2
        # codata[1][1]: row 2, col 2
        for i in range(0, len(codata)):
            #  METHOD   ALIAS  a0  a1  a2  a3  a4
            coeffsDict[codata[i][1]] = (codata[i][2], codata[i][3],
                                        codata[i][4], codata[i][5],
                                        codata[i][6])

    calEner = np.zeros(inEner.size, dtype=float)
    ie = 0

    #
    # convert energies
    #
    for ie in range(0, inEner.size):
        # read fitting coeffs taken from polyfit2Bias.R (a0, a1, a2, a3)
        #  as in y = a0 + a1*x + a2*x^2 + a3*x^3
        # where y=E_reconstructed and x=Ecalibration (keV)

        coeffs = coeffsDict[alias]
        # print("SIGNAL=", inEner[ie], "coeffs=", coeffs)
        npCoeffs = np.array(coeffs)
        # print("npCoeffs=",npCoeffs)

        # subtract "y" (0 = a0 + a1*x + a2*x^2 + ... - y) :
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
                       key=lambda x: abs(x[1]-inEner[ie]))[1]  # (idx,value)
        # print("For:", alias, " Recon energy=",rclosest)

        calEner[ie] = rclosest

    # return calibrated energies
    # ------------------------------------
    return calEner


def RfromItrigger(infile, Icol, current, Rcol, Rmethod, samplingrate):
    """
    :param infile: fits input file with current column Icol
    :param Icol: existing current column name
    :param current: input current: ADC or AMP (PULSE0000)
    :param Rcol: new output R column name
    :param Rmethod: resistance method calculation(R or RALL or RNOL or RFITTED)
    :param samplingrate: sampling rate (Hz-1)
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
    Ifit = 45.3E-6
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
            comm = ("ftcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes column=AMP expr='" + Icol +
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
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=" + Rcol +
                    " expr='#R0 - #R0*((abs(" + Icol +
                    "-#I0_START)/#I0_START)" +
                    "/(1 + abs(" + Icol + "-#I0_START)/#I0_START))'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for R calculation:")
            print(comm)
            raise
    elif Rmethod == "RALL":
        try:
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=I expr='#I0_START-" +
                    Icol + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=DER expr='seqdiff(I)'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=UNOS expr=1")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcollen infile=" + infile + " colname=UNOS collen="
                    + str(triggsz))
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=DER1COL expr=DER")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcollen infile=" + infile + " colname=DER1COL collen=1")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcollen infile=" + infile + " colname=DER1COL collen="
                    + str(triggsz))
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=DER0COL1 expr='DER-DER1COL'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=DER2 expr='UNOS*DER[2]'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=DERFINAL expr='DER2+DER0COL1'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=DERIT expr='DERFINAL*" +
                    str(samplingrate) + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fdelcol infile=" + infile +
                    "+1 colname=UNOS confirm=no proceed=yes")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fdelcol infile=" + infile +
                    "+1 colname=DER confirm=no proceed=yes")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fdelcol infile=" + infile +
                    "+1 colname=DER1COL confirm=no proceed=yes")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fdelcol infile=" + infile +
                    "+1 colname=DER0COL1 confirm=no proceed=yes")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fdelcol infile=" + infile +
                    "+1 colname=DER2 confirm=no proceed=yes")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            # R = (V0 - I*RL - L*dI/dt)/I
            comm = ("fcalc infile=" + infile + " outfile=" +
                    infile + " clobber=yes clname=" + Rcol +
                    " expr='((#I0_START*(#R0+" +
                    str(RL) + ")) - (I*" + str(RL) + ") - (" +
                    str(L) + "*DERIT))/I'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for RALL calculation:")
            print(comm)
            raise
    elif Rmethod == "RNOL":

        try:
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=I expr='#I0_START-" +
                    Icol + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)

            # R = (V0 - I*RL)/I
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=" + Rcol +
                    " expr='((#I0_START*(#R0+" + str(RL) + ")) - (I*" +
                    str(RL) + "))/I'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for RNOL calculation:")
            print(comm)
            raise
    elif Rmethod == "RFITTED":

        try:
            # IIF = I + Ifit = (I0 - AMP) + Ifit
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=IIF expr='#I0_START-" +
                    Icol + "+" + str(Ifit) + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)

            # R = V0/IIF
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=" + Rcol +
                    " expr='#I0_START*(#R0+" + str(RL) + ")/IIF'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for RFITTED calculation:")
            print(comm)
            raise

    return


def RfromIstream(infile, Icol, current, Rcol, Rmethod, samplingrate):
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
     AMP

    """
    Ifit = 45.3E-6
    comm = ""

    fstr = fits.open(infile)
    RPARA = fstr[1].header['RPARA']
    TTR = fstr[1].header['TTR']
    LFILTER = fstr[1].header['LFILTER']
    RL = RPARA / (TTR * TTR)
    L = LFILTER / (TTR * TTR)

    if current == "ADC":
        # create col AMP in fits file from ADC column (fcalc does not work)
        try:
            comm = ("ftcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes column=AMP expr='" + Icol +
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
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=" + Rcol +
                    " expr='#R0 - #R0*((abs(" + Icol +
                    "-#I0_START)/#I0_START)" +
                    "/(1 + abs(" + Icol + "-#I0_START)/#I0_START))'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for R calculation:")
            print(comm)
            raise
    elif Rmethod == "RALL":
        try:
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=I expr='#I0_START-" +
                    Icol + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=DER expr='seqdiff(I)'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=DERIT expr='DER*" +
                    str(samplingrate) + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            # R = (V0 - I*RL - L*dI/dt)/I
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=" + Rcol +
                    " expr='((#I0_START*(#R0+" + str(RL) + ")) - (I*" +
                    str(RL) + ") - (" + str(L) + "*DERIT))/I'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for RALL calculation:")
            print(comm)
            raise
    elif Rmethod == "RNOL":
        try:
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=I expr='#I0_START-" +
                    Icol + "'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
            # R = (V0 - I*RL)/I
            comm = ("fcalc infile=" + infile + " outfile=" + infile +
                    " clobber=yes clname=" + Rcol +
                    " expr='((#I0_START*(#R0+" + str(RL) + ")) - (I*" +
                    str(RL) + "))/I'")
            args = shlex.split(comm)
            check_call(args, stderr=STDOUT)
        except:
            print("Error running ftool for RNOL calculation:")
            print(comm)
            raise
    elif Rmethod == "RFITTED":
            try:
                # IIF = I + Ifit = (I0 - AMP) + Ifit
                comm = ("fcalc infile=" + infile + " outfile=" + infile +
                        " clobber=yes clname=IIF expr='#I0_START-" +
                        Icol + "+" + str(Ifit) + "'")
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)

                # R = V0/IIF
                comm = ("fcalc infile=" + infile + " outfile=" + infile +
                        " clobber=yes clname=" + Rcol +
                        " expr='#I0_START*(#R0+" + str(RL) + ")/IIF'")
                args = shlex.split(comm)
                check_call(args, stderr=STDOUT)
            except:
                print("Error running ftool for RFITTED calculation:")
                print(comm)
                raise

    return
