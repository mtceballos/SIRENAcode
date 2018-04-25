"""

Identify PHotons detected/missed by SIRENA using time impacts in piximact files
and start times in event files. The used criterium  is that
time difference <= accuracy

"""
# ----IMPORT MODULES --------------
from __future__ import print_function
from __future__ import division
import math
import numpy as np
import auxpy
from astropy.io import fits
from collections import Counter


def identifyPH(impF, evtF, tssF, acc, samprate, ppr):

    """
        :param impF: FITS file with impact times(from simput or tesconstpileup)
        :param evtF: FITS file with SIRENA-detected events
        :param tssF: FITS file with tessim-simulated events
        :param acc: accuracy for time comparison
        :param samprate: sampling rate (Hz)
    """

    # Impacts
    imp = fits.open(impF, memmap=True)
    impTab = imp[1].data
    impTimes = impTab['TIME']
    impIDs = impTab['PH_ID']
    print("Number of piximpact photons: ", len(impTimes))

    # tessim-simulated pulses
    tss = fits.open(tssF, memmap=True)
    tssTab = tss[1].data
    # tssRecords = tssTab['TIME']
    nettot = tss[1].header["NETTOT"]
    # to account for PAIRS, SINGLES... :
    # tssPhs = tssTab.size * len(re.findall("keV", tssF))
    tssPhs = tssTab.size * int(ppr)
    print("Number of simulated photons:", tssPhs)
    if nettot != tssPhs:
        print("     Inconsitent NETTOT value (", nettot, ") for ",
              tssPhs, "simulated photons")

    # Detected events
    evt = fits.open(evtF, memmap=True)
    evtTab = evt[1].data
    evtTimes = evtTab['TIME']
    # evtIDs = np.zeros(evtTab.size, dtype=int)
    evtIDs = evtTab['PH_ID']
    print("Number of detected photons:", len(evtTimes))
    # print("Detected/tessim-simulated=",
    #                              '{:6.2f}'.format(100*len(evtIDs)/tssPhs),"%")
    # print("Detected/piximpact=",
    #                         '{:6.2f}'.format(100*len(evtIDs)/len(impIDs)),"%")

    # IDENTIFICATION
    # ===============
    clsidx = auxpy.find_nearest_idx(impTimes, evtTimes)
    evtIDs = impIDs[clsidx]

    # print(evtIDs)
    # print(impIDs)

    # ARRIVAL PHASES
    # ===============
    arrPhases = auxpy.arrivalPhase(impF, evtF, samprate)
    print("Arrival Phases=", arrPhases)

    for i in range(0, evtTab.size):
        # print("Arrival Phase: ",arrPhases[i])
        if math.fabs(impTimes[clsidx[i]] - evtTimes[i]) > acc:
            print("No Photon in IMP File close enough to Photon in row ", i+1,
                  " of EVT file")
            print("Closest row in impact=", clsidx[i]+1)
            print("imp=", impTimes[clsidx[i]], "evt=", evtTimes[i])
            print("dist=", impTimes[clsidx[i]] - evtTimes[i])

    # check that identifications are unique
    # =======================================
    # only python 3!
    # vals, inverse, count =
    #               np.unique(evtIDs, return_inverse=True, return_counts=True)
    # idx_phs_repeated = np.where(count > 1)[0]
    # phs_repeated = vals[idx_phs_repeated]
    phs_repeated = [item for item, count in Counter(evtIDs).iteritems()
                    if count > 1]
    if len(phs_repeated) > 0:
        print(
                "These photons have multiple identifications (",
                len(phs_repeated), "):", phs_repeated)

    # check missed photons (in IMPlist but not in EVTlist)
    # ======================================================
    phMissed = np.setdiff1d(impIDs, evtIDs, assume_unique=True)
    # print(phMissed)
    # print(phMissed.size, len(impTimes), tssPhs)
    # print(phMissed.size - (len(impTimes) - tssPhs),"/", len(phs_repeated))
    # print("Piximpact photons missed in evtlist(", len(phMissed) ,"):\n", phMissed)
    print("Number of missed photons: ", phMissed.size-(len(impTimes)-tssPhs))

    if ppr == "2":
        phMissed_Evens = []
        phMissed_Odds = []
        index_Odds = 0
        index_Evens = 0

        for i in range(0, phMissed.size):
            if phMissed[i] % 2 == 0:
                # print(phMissed[i])
                phMissed_Evens.insert(index_Evens, phMissed[i])
                index_Evens = index_Evens + 1
            else:
                # print(phMissed[i])
                phMissed_Odds.insert(index_Odds, phMissed[i])
                index_Odds = index_Odds + 1

        # print("Piximpact P photons missed in evtlist (",index_Odds,"):\n", phMissed_Odds)
        # print("Piximpact S photons missed in evtlist (",index_Evens,"):\n", phMissed_Evens)
        # print("Piximpact P photons missed in evtlist = ",index_Odds,"\n")
        # print("Piximpact S photons missed in evtlist = ", index_Evens, "\n")
        print("Missed Primaries/Secondaries/repeated")
        print(index_Odds - 2 - (len(impTimes) - tssPhs - 4) / 2, "/",
              index_Evens - 2 - (len(impTimes) - tssPhs - 4) / 2,
              "/", len(phs_repeated))

    # check photons unidentified (fake pulses)
    # ========================================
    # (python3)
    # idx_phs_fake = np.where(vals == 0)[0]
    # nFake = len(idx_phs_fake)
    # (end python3)
    phs_fake = [count for item, count in Counter(evtIDs).iteritems()
                if item == 0]
    nFake = len(phs_fake)
    print("Number of Fake pulses:", nFake)
    
    imp.close()
    evt.close()
    tss.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
            description='Identify PHotons detected/missed by SIRENA',
            prog='identifyPH')
    parser.add_argument('--impactFile', help='file with impact times',
                        required=True)
    parser.add_argument('--eventFile', help='file with events arrival times',
                        required=True)
    parser.add_argument('--tessimFile', help='file with tessim simul events',
                        required=True)
    parser.add_argument('--accuracy', type=float,
                        help='accuracy for time comparison', default=9.7E-6)
    parser.add_argument('--samprate', type=float, help='sampling rate (Hz)',
                        default=156250.)
    parser.add_argument('--ppr', help='pulses per record', choices=['1', '2'],
                        default=2)

    inargs = parser.parse_args()

    identifyPH(impF=inargs.impactFile, evtF=inargs.eventFile,
               tssF=inargs.tessimFile, acc=inargs.accuracy,
               samprate=inargs.samprate, ppr=inargs.ppr)
