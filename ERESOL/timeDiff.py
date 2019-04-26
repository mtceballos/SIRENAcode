"""

Calculate pulse TIME difference (in s and samples) with reference timeline
in a given sample
"""
# ----IMPORT MODULES --------------
from __future__ import print_function
from __future__ import division
from astropy.io import fits


def timeDiff(file, sample, plsnum, samprate):

    imp = fits.open(file, memmap=True)
    impTab = imp[1].data
    impTimes = impTab['TIME']

    # Calculate TIME for sample (TIME=0.0 for sample=1)
    timeForSample = (sample-1) / samprate  # in s
    timeDiffsec = impTimes[plsnum-1] - timeForSample
    timeDiffsmp = timeDiffsec * samprate
    # print("Timediff (sec,sam)=", timeDiffsec, timeDiffsmp)
    print("FileTime(sec)= %.8f" % (impTimes[plsnum-1]))
    print("RefTime(sec)= %.8f" % (timeForSample))
    print("Timediff(sam)= %.3f" % (timeDiffsmp))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
            description='Time differences',
            prog='timeDiff.py')
    parser.add_argument('--file', help='file1',
                        required=True)
    parser.add_argument('--sample', type=int, help='reference sample',
                        required=True)
    parser.add_argument('--plsnum', type=int, help='pulse number',
                        required=True)
    parser.add_argument('--samprate', type=float, help='sampling rate (Hz)',
                        default=156250.)

    inargs = parser.parse_args()

    timeDiff(file=inargs.file, sample=inargs.sample, plsnum=inargs.plsnum,
             samprate=inargs.samprate)
