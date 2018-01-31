from __future__ import print_function
import os
import shlex
import shutil
import tempfile
from time import gmtime, strftime
from subprocess import check_call, STDOUT
from astropy.io import fits

tmpDir = tempfile.mkdtemp()
os.environ["PFILES"] = tmpDir + ":" + os.environ["PFILES"]
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null/"


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

    fsim = fits.open(simfile)
    nrows = fsim[1].header["NAXIS2"]
    fsim.close()
    assert nrows > 1, "Tessim failed for (%s): just one row present " % simfile

    try:
        comm1 = "fdelrow infile=" + simfile + "+1 firstrow=" + str(nrows) + " nrows=1 confirm=no proceed=yes"
        args = shlex.split(comm1)
        check_call(args, stderr=STDOUT)
        comm2 = "fdelrow infile=" + simfile + "+1 firstrow=1 nrows=1 confirm=no proceed=yes"
        args = shlex.split(comm2)
        check_call(args, stderr=STDOUT)
        fsim = fits.open(simfile, mode='update')
        nrows2 = fsim[1].header['NAXIS2']
        assert nrows2 == nrows-2, "Failure removing initial & last rows in (%s): " % simfile
        nettot = nrows2 * ppr  # new number of pulses (=nofrecords in LPA2; ==2*nofrecords in LPA1)
        fsim[1].header['NETTOT'] = nettot
        # update HISTORY in header[0]
        dateTime = strftime("%Y-%m-%d %H:%M:%S", gmtime())
        fsim[0].header['HISTORY'] = "Created & Updated by rmLastAndFirst.py on " + dateTime + "with commands: "
        # split command line in blocks of 60 chars for HISTORY keyword
        commSplit1 = [comm1[i:i + 60] for i in range(0, len(comm1), 60)]
        for i in range(len(commSplit1)):
            fsim[0].header['HISTORY'] = commSplit1[i]
        # split command line in blocks of 60 chars for HISTORY keyword
        commSplit2 = [comm1[i:i + 60] for i in range(0, len(comm1), 60)]
        for i in range(len(commSplit2)):
            fsim[0].header['HISTORY'] = commSplit2[i]
        fsim.close()

    except:
        print("Error running FTOOLS to remove initial & last rows in ", simfile)
        shutil.rmtree(tmpDir)
        raise

