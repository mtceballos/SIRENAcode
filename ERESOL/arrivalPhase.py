#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 12:51:29 2020

@author: ceballos
"""

from astropy.io import fits
from os import path
import numpy as np
from find_nearest_idx import find_nearest_idx

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
    if path.isfile(imp0F):
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
