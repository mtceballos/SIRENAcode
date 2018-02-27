from __future__ import print_function
import os
import shlex
import shutil
import tempfile
from time import gmtime, strftime
from subprocess import check_call, check_output,STDOUT
from astropy.io import fits

tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"


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
    commVerify = "fverify infile=" + fitsfile + " outfile=" + tmpFile + " clobber=yes"
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
    #getNumerrs = "pget fverify numerrs"
    #getNumwrns = "pget fverify numwrns"
    #numerrs = 0
    #numwrns = 0
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
    #if os.path.isfile("pp.stdout"):
    #    os.remove("pp.stdout")
    # if fverify detects errors/warnings, return numerrs+numwrns
    #return numerrs+numwrns


def rmLastAndFirst(simfile, ppr):
    """
    rm first (and LAST) record of SIXTE/tessim simulated file and update NETTOT
    (first starts high and last can be cut) (see Christian's email from 19 Jan 2017 @ EURECA)
    Also update number of pulses in NETTOT keyword

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
    comm1 = "fdelrow infile=" + simfile + "+1 firstrow=" + str(nrows) + " nrows=1 confirm=no proceed=yes"
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
        print("Errors running FTOOLS to remove final row (abandon) in ", simfile, " (returned to original file)")
    elif warns1 > 0:
        print("Warnings running FTOOLS to remove final row (abandon) in ", simfile, " (keep modified file)")
    else:
        comm2 = "fdelrow infile=" + simfile + "+1 firstrow=1 nrows=1 confirm=no proceed=yes"
        try:
            args = shlex.split(comm2)
            check_call(args, stderr=STDOUT)
            print("     Initial row removed in ", simfile)
        except:
            print("Error running ", comm2, " to remove initial row in ", simfile)
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
            print("Errors/warnings after FTOOLS removal of initial row in ", simfile)

    os.remove(tmpFile)

def updateHISTORY(file, history):
    """
    
    :param file: File for which the keyword HISTORY in the HEADER[0] will be updated 
    :param history: String for the HISTORY keyword (if large it will be splitted in several lines
    :return:
     
    """

    dateTime = strftime("%Y-%m-%d %H:%M:%S", gmtime())
    f = fits.open(file, mode='update')
    f[0].header['HISTORY'] = "Created/Updated on " + dateTime + "with commands: "
    # split command line in blocks of 60 chars for HISTORY keyword
    histSplit = [history[i:i + 60] for i in range(0, len(history), 60)]
    for i in range(len(histSplit)):
        f[0].header['HISTORY'] = histSplit[i]
    f.close()
