#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 16:44:30 2020

@author: ceballos

# Calibration lines

Source: CALDB20161122 + https://xdb.lbl.gov/Section1/Table_1-2.pdf

"""



import numpy as np
from RxLines import RxLines

# MnKa
ilabels = ["Ka11","Ka12", "Ka13", "Ka14", "Ka15", "Ka16",  "Ka21", "Ka22"]
energies_eV = np.array([5898.882, 5897.898, 5894.864, 5896.566, 5899.444, 5902.712, 5887.772, 5886.528], dtype=np.float64)
fwhms_eV = np.array([1.7145, 2.0442, 4.4985, 2.6616, 0.97669, 1.5528, 2.3604, 4.2168], dtype=np.float64)
rel_amplitudes = np.array([0.784, 0.263, 0.067, 0.095, 0.071, 0.011, 0.369, 0.1], dtype=np.float64)
areas = np.array([0.3523,0.1409,0.07892,0.06624,0.01818,0.004475,0.2283,0.1106], dtype=np.float64)
MnKas = RxLines(complabel="MnKa", ilabels=ilabels, energies_eV=energies_eV, fwhms_eV=fwhms_eV, rel_amplitudes=rel_amplitudes)
MnKas_cmass = '{:0.2f}'.format(np.sum(energies_eV*areas)/np.sum(areas))
MnKa1_cmass = '{:0.2f}'.format(np.sum(energies_eV[0:6]*areas[0:6])/np.sum(areas[0:6]))
MnKa2_cmass = '{:0.2f}'.format(np.sum(energies_eV[6:]*areas[6:])/np.sum(areas[6:]))
#Ka2_cmass = 5887.65 # eV wikipedia nominal
#Ka1_cmass = 5898.75 # eV wikipedia nominal
#Ka2_cmass = 5887.772 # eV most intense
#Ka1_cmass = 5898.882 # eV most intense

# MnKb
ilabels = ['Kb1', 'Kb2', 'Kb3', 'Kb4', 'Kb5']
energies_eV = np.array([6490.89, 6486.31, 6477.73, 6490.06, 6488.83], dtype=np.float64)
fwhms_eV = np.array([1.83, 9.4, 13.22, 1.81, 2.81], dtype=np.float64)
rel_amplitudes = np.array([0.608, 0.109, 0.077, 0.397, 0.176], dtype=np.float64)
areas = np.array([0.254,0.234,0.234,0.164,0.114], dtype=np.float64)
MnKb = RxLines(complabel="MnKb", ilabels=ilabels, energies_eV=energies_eV, fwhms_eV=fwhms_eV, rel_amplitudes=rel_amplitudes)
MnKb_cmass = '{:0.2f}'.format(np.sum(energies_eV*areas)/np.sum(areas))
#MnKb_cmass = 6486.38 # eV
#MnKb_cmass = 6490.89 # eV "most intense"

# Cr-Ka
ilabels = ["Ka11","Ka12", "Ka13", "Ka14", "Ka15", "Ka21", "Ka22"]
energies_eV = np.array([5414.874,5414.099,5412.745,5410.583,5418.304,5405.551,5403.986], dtype=np.float64)
fwhms_eV = np.array([1.457,1.760,3.138,5.149,1.988,2.2224,4.740], dtype=np.float64)
rel_amplitudes = np.array([0.822,0.237,0.085,0.045,0.015,0.386,0.036], dtype=np.float64)
areas = np.array([0.378,0.132,0.084,0.073,0.009,0.271,0.054], dtype=np.float64)
CrKas = RxLines(complabel="CrKa", ilabels=ilabels, energies_eV=energies_eV, fwhms_eV=fwhms_eV, rel_amplitudes=rel_amplitudes)
CrKas_cmass = '{:0.2f}'.format(np.sum(energies_eV*areas)/np.sum(areas))
CrKa1_cmass = '{:0.2f}'.format(np.sum(energies_eV[0:5]*areas[0:5])/np.sum(areas[0:5]))
CrKa2_cmass = '{:0.2f}'.format(np.sum(energies_eV[5:]*areas[5:])/np.sum(areas[5:]))
#print(CrKas_cmass, CrKa1_cmass, CrKa2_cmass)

#Cr-Kb
ilabels = ['Kb1', 'Kb2', 'Kb3', 'Kb4', 'Kb5']
energies_eV = np.array([5947.00,5935.31,5946.24,5942.04,5944.93], dtype=np.float64)
fwhms_eV = np.array([1.70,15.98,1.90,6.69,3.37], dtype=np.float64)
rel_amplitudes = np.array([0.670,0.055,0.337,0.082,0.151], dtype=np.float64)
areas = np.array([0.307,0.236,0.172,0.148,0.137], dtype=np.float64)
CrKb = RxLines(complabel="CrKb", ilabels=ilabels, energies_eV=energies_eV, fwhms_eV=fwhms_eV, rel_amplitudes=rel_amplitudes)
CrKb_cmass = '{:0.2f}'.format(np.sum(energies_eV*areas)/np.sum(areas))

# Cu-Ka
ilabels = ["Ka11","Ka12", "Ka21", "Ka22"]
energies_eV = np.array([8047.837,8045.367,8027.993,8026.504], dtype=np.float64)
fwhms_eV = np.array([2.285,3.358,2.666,3.571], dtype=np.float64)
rel_amplitudes = np.array([0.957,0.090,0.334,0.111], dtype=np.float64)
areas = np.array([0.579,0.080,0.236,0.105], dtype=np.float64)
CuKas = RxLines(complabel="CuKa", ilabels=ilabels, energies_eV=energies_eV, fwhms_eV=fwhms_eV, rel_amplitudes=rel_amplitudes)
CuKas_cmass = '{:0.2f}'.format(np.sum(energies_eV*areas)/np.sum(areas))
CuKa1_cmass = '{:0.2f}'.format(np.sum(energies_eV[0:2]*areas[0:2])/np.sum(areas[0:2]))
CuKa2_cmass = '{:0.2f}'.format(np.sum(energies_eV[2:]*areas[2:])/np.sum(areas[2:]))
#print(CuKas_cmass)

#Cu-Kb
ilabels = ['Kb1', 'Kb2', 'Kb3', 'Kb4', 'Kb5']
energies_eV = np.array([8905.532,8903.109,8908.462,8897.387,8911.393], dtype=np.float64)
fwhms_eV = np.array([3.52,3.52,3.55,8.08,5.31], dtype=np.float64)
rel_amplitudes = np.array([0.757,0.388,0.171,0.068,0.055], dtype=np.float64)
areas = np.array([0.485,0.248,0.110,0.100,0.055], dtype=np.float64)
CuKb = RxLines(complabel="CuKb", ilabels=ilabels, energies_eV=energies_eV, fwhms_eV=fwhms_eV, rel_amplitudes=rel_amplitudes)
CuKb_cmass = '{:0.2f}'.format(np.sum(energies_eV*areas)/np.sum(areas))
#print(CuKb_cmass)

# other lines
ScKa1_cmass = 4090.6
ScKa2_cmass = 4086.1
ScKb_cmass = 4460.5

TiKa1_cmass = 4510.84
TiKa2_cmass = 4504.86
TiKb_cmass = 4931.81

VKa1_cmass = 4952.20
VKa2_cmass = 4944.64
VKb_cmass = 5427.29

FeKa1_cmass = 6403.84
FeKa2_cmass = 6390.84
FeKb_cmass = 7057.98

CoKa1_cmass = 6930.32
CoKa2_cmass = 6915.30
CoKb_cmass = 7649.43

NiKa1_cmass = 7478.15
NiKa2_cmass = 7460.89
NiKb_cmass = 8264.66

ZnKa1_cmass = 8638.86
ZnKa2_cmass = 8615.78
ZnKb_cmass = 9572.0

GeKa1_cmass = 9886.42
GeKa2_cmass = 9855.32
GeKb_cmass = 10982.1

BrKa1_cmass = 11924.2
BrKa2_cmass = 11877.6
BrKb_cmass = 13291.4
