from __future__ import print_function
import os
import shlex
import shutil
import tempfile
from time import gmtime, strftime
import numpy as np
import math
from subprocess import check_call, check_output,STDOUT
from astropy.io import fits
from astropy.io import ascii


tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"


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


def arrivalPhase(impF, evtF, samprate):
    """
        Calculate the arrival phase of the simulated photons:
            distance from true nonjitter simulated time and detected time
        :param impF: FITS file with impact times(from simput or tesconstpileup)
        :param evtF: FITS file with SIRENA-detected events
        :param samprate: sampling rate (Hz)
        :return arrPhase:
    """
    # PIXIMPACT
    imp = fits.open(impF, memmap=True)
    impTab = imp[1].data
    impTimes = impTab['TIME']

    # SIRENA reconstruction
    evt = fits.open(evtF, memmap=True)
    evtTab = evt[1].data
    evtTimes = evtTab['TIME']

    clsidx = find_nearest_idx(impTimes, evtTimes)

    # offset with imp times
    off1 = (evtTimes - impTimes[clsidx]) * samprate
    # offset of imp times with nonjitter times: dec part of imptime
    off2 = np.modf(impTimes[clsidx] * samprate)[0]
    arrPhase1 = off1 + off2

    # PIXIMPACT NO JITTER (if it exists):
    imp0F = impF.replace("keV_jitter_", "keV_")
    if os.path.isfile(imp0F):
        imp0 = fits.open(imp0F, memmap=True)
        imp0Tab = imp0[1].data
        imp0Times = imp0Tab['TIME']
        clsidx0 = find_nearest_idx(imp0Times, evtTimes)
        # arrPhase 2 should be the same that arrPhase1
        arrPhase2 = (evtTimes - imp0Times[clsidx0]) * samprate
        assert np.allclose(arrPhase1, arrPhase2), \
            "Incompatible arrival Phases calculations"

    arrPhase3 = evtTimes*samprate - np.round(evtTimes*samprate, 0)
    # arrPhase3 = np.modf(evtTimes*samprate)[0]  # dec part of evtTime
    return arrPhase1
    # return arrPhase3


def fitsVerify(fitsfile):
    """
    :param fitsfile: input file to be verified
    :return: number of errors+ warnings in verification process
    """

    def getNumerrs(file):
        infile = open(file, "r")
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
    except:
        print("Error running fverify command")
        shutil.rmtree(tmpDir)
        raise
    numerrs, numwrns = getNumerrs(tmpFile)
    os.remove(tmpFile)
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
    except:
        print("Error running ", comm1, " to remove final row in ", simfile)
        shutil.rmtree(tmpDir)
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
        except:
            print("Error running ", comm2, " to remove initial row in ",
                  simfile)
            shutil.rmtree(tmpDir)
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

    os.remove(tmpFile)


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
        # print("For:", alias, " Recon energy=",rclosest,"npCoeffs=",npCoeffs)
        calEner[ie] = rclosest

    # return calibrated energies
    # ------------------------------------
    return calEner
