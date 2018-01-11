"""

Identify PHotons detected/missed by SIRENA using time impacts in piximact files and start times in event files
The used criterium  is that time difference <= accuracy

"""
# ----IMPORT MODULES --------------
from __future__ import print_function
from __future__ import division
import math
import sys
import numpy as np
from astropy.io import fits
from collections import Counter


def find_nearest_idx(array,value):

    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        idx -= 1
    return idx


def identifyPH(impF, evtF, tssF, acc):

    """
        :param impF: FITS file with impact times (from simput or tesconstpileup)
        :param evtF: FITS file with SIRENA-detected events
        :param tssF: FITS file with tessim-simulated events
        :param acc: accuracy for time comparison
    """

    imp = fits.open(impF, memmap=True)
    impTab = imp[1].data
    impTimes = impTab['TIME']
    impIDs = impTab['PH_ID']

    evt = fits.open(evtF, memmap=True, mode='update')
    evtTab = evt[1].data
    evtTimes = evtTab['TIME']
    # evtIDs = np.zeros(evtTab.size, dtype=int)
    evtIDs = evtTab['PH_ID']

    for i in range(0, evtTab.size):
        clsidx = find_nearest_idx(impTimes, evtTimes[i])
        if math.fabs(impTimes[clsidx] - evtTimes[i]) > acc:
            print("No Photon in IMP File close enough to Photon in row ", i+1, " of EVT file")
            print("Closest row in impact=", clsidx+1)
            print("imp=", impTimes[clsidx], "evt=", evtTimes[i])
            print("dist=",impTimes[clsidx] - evtTimes[i])
            #sys.exit()
        evtIDs[i] = impIDs[clsidx]

    # check that identifications are unique
    # only python 3!
    # vals, inverse, count = np.unique(evtIDs, return_inverse=True, return_counts=True) # python3
    # idx_phs_repeated = np.where(count > 1)[0]
    # phs_repeated = vals[idx_phs_repeated]
    phs_repeated = [item for item, count in Counter(evtIDs).iteritems() if count > 1]

    if len(phs_repeated) > 0:
        print("These photons have multiple identifications:", phs_repeated)
    else:
        # check missed photons (in piximpact but not in evtlist)
        phMissed = np.setdiff1d(impIDs, evtIDs, assume_unique=True)
        print("Piximpact photons missed in evtlist):\n", phMissed)

        # check photons unidentified (fake pulses)
        # (python3)
        # idx_phs_fake = np.where(vals == 0)[0]
        # nFake = len(idx_phs_fake)
        phs_fake = [count for item, count in Counter(evtIDs).iteritems() if item == 0]
        nFake = len(phs_fake)
        print("Number of Fake pulses:", nFake)

        # Number of detected pulses
        print("Number of 'detected' pulses:", len(evtIDs))

        if tssF:
            tss = fits.open(tssF, memmap=True)
            tssPhs = tss[1].header["NETTOT"]
            tss.close()
            print("Number of tessim-simulated pulses:", tssPhs)
            print("Detected/tessim-simulated=", '{:6.2f}'.format(100*len(evtIDs)/tssPhs),"%")

        print("Detected/piximpact=", '{:6.2f}'.format(100*len(evtIDs)/len(impIDs)),"%")

    imp.close()
    evt.close()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Identify PHotons detected/missed by SIRENA',
                                     prog='identifyPH')
    parser.add_argument('--impactFile', help='file with impact times', required=True)
    parser.add_argument('--eventFile', help='file with events arrival times', required=True)
    parser.add_argument('--tessimFile', help='file with tessim simulated events')
    parser.add_argument('--accuracy', type=float, help='accuracy for time comparison', default=9.7E-6)

inargs = parser.parse_args()

identifyPH(impF=inargs.impactFile, evtF=inargs.eventFile, tssF=inargs.tessimFile, acc=inargs.accuracy)
