#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 11:27:06 2019

@author: ceballos
"""

from __future__ import print_function
from os import rename

goodList = [1010, 1060, 1110, 115, 1160, 118, 1210, 122, 125, 1260, 128,
            130, 1310, 1360, 1410, 1460, 15, 150, 1510, 1560, 1610, 1660,
            170, 1710, 172, 1760, 178, 1810, 182, 1860, 188, 190, 1910,
            192, 1960, 198, 2010, 202, 2060, 208, 210, 2110, 212, 2160,
            218, 2210, 222, 2260, 228, 230, 2310, 232, 2360, 2410, 2460,
            250,2510, 255, 2560, 258, 260, 2610, 2660, 2710, 2760, 2810,
            2860, 2910, 2960, 3010, 3060, 310, 3110, 3160, 3210, 3260,
            3310, 3360, 3410, 3460, 35, 3510, 3560, 360, 3610, 3660, 3710,
            3760, 3810, 3860, 3910, 3950, 3960, 3975, 4000, 4010, 4025, 4050,
            4060, 4075, 4095, 4098, 410, 4100, 460, 510, 55, 560, 610, 660,
            710, 75, 760, 810, 860, 910, 95, 960]

newList = ['{:05d}'.format(item) for item in goodList]


oldStr = ("00015 00035 00055 00075 00095 00115 00118 00122 00125 00128 00130" +
          " 00150 00170 00172 00178 00182 00188 00190 00192 00198 00202 00208" +
          " 00210 00212 00218 00222 00228 00230 00232 00250 00255 00258 00260" +
          " 00310 00360 00410 00460 00510 00560 00610 00660 00710 00760 00810" +
          " 00860 00910 00960 01010 01060 01110 01160 01210 01260 01310 01360" +
          " 01410 01460 01510 01560 01610 01660 01710 01760 01810 01860 01910" +
          " 01960 02010 02060 02110 02160 02210 02260 02310 02360 02410 02460" +
          " 02510 02560 02610 02660 02710 02760 02810 02860 02910 02960 03010" +
          " 03060 03110 03160 03210 03260 03310 03360 03410 03460 03510 03560" +
          " 03610 03660 03710 03760 03810 03860 03910 03950 03960 03975 04000" +
          " 04010 04025 04050 04060 04075 04095 04098 04100")
oldList = oldStr.split()

if not len(newList) == len(oldList):
    raise Exception("Different list sizes")

Eprims = [0.2, 0.5, 1, 4, 6, 8]
Esecs = [0.2, 0.5, 0.8, 1, 1.3, 2, 2.3, 3, 3.3,
         4, 4.3, 5, 5.3, 6, 6.3, 7, 7.3, 8]
evtdir = "/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/eresol/nodetSP/"
simdir = "/home/ceballos/INSTRUMEN/EURECA/ERESOL/PAIRS/xifusimLPA75um/"

# ============================================================
# Renaming EVT files
for io in range(len(oldList)):
    for pL in (4096, 256, 128):
        for E1 in Eprims:
            for E2 in Esecs:
                oldFile = (evtdir + "events_sep" + oldList[io] +
                           "sam_400p_SIRENA4096_pL" + str(pL) + "_" +
                           str(E1) + "keV_" + str(E2) +
                           "keV_STC_T_fixedlib6OF_" +
                           "OPTFILT4096_samprate2_jitter_dcmt200.fits")
                tmpFile = (evtdir + "tmp_sep" + newList[io] +
                           "sam_400p_SIRENA4096_pL" + str(pL) + "_" +
                           str(E1) + "keV_" + str(E2) +
                           "keV_STC_T_fixedlib6OF_" +
                           "OPTFILT4096_samprate2_jitter_dcmt200.fits")
                rename(oldFile, tmpFile)

for io in range(len(oldList)):
    for pL in (4096, 256, 128):
        for E1 in Eprims:
            for E2 in Esecs:
                tmpFile = (evtdir + "tmp_sep" + newList[io] +
                           "sam_400p_SIRENA4096_pL" + str(pL) + "_" +
                           str(E1) + "keV_" + str(E2) +
                           "keV_STC_T_fixedlib6OF_" +
                           "OPTFILT4096_samprate2_jitter_dcmt200.fits")
                newFile = (evtdir + "events_sep" + newList[io] +
                           "sam_400p_SIRENA4096_pL" + str(pL) + "_" +
                           str(E1) + "keV_" + str(E2) +
                           "keV_STC_T_fixedlib6OF_" +
                           "OPTFILT4096_samprate2_jitter_dcmt200.fits")
                rename(tmpFile, newFile)
# ============================================================
#  Renaming SIM files
#for io in range(len(oldList)):
#    for E1 in Eprims:
#        for E2 in Esecs:
#
#            oldFile = (simdir + "sep" + oldList[io] + "sam_400p_" +
#                       str(E1) + "keV_" + str(E2) + "keV_" +
#                       "samprate2_jitter_dcmt200.fits")
#            tmpFile = (simdir + "tmp" + newList[io] + "sam_400p_" +
#                       str(E1) + "keV_" + str(E2) + "keV_" +
#                       "samprate2_jitter_dcmt200.fits")
#            rename(oldFile, tmpFile)
#for io in range(len(oldList)):
#    for E1 in Eprims:
#        for E2 in Esecs:
#            tmpFile = (simdir + "tmp" + newList[io] + "sam_400p_" +
#                       str(E1) + "keV_" + str(E2) + "keV_" +
#                       "samprate2_jitter_dcmt200.fits")
#            newFile = (simdir + "sep" + newList[io] + "sam_400p_" +
#                       str(E1) + "keV_" + str(E2) + "keV_" +
#                       "samprate2_jitter_dcmt200.fits")
#            print("Renaming:", oldFile, " to", newFile)
#            rename(tmpFile, newFile)
