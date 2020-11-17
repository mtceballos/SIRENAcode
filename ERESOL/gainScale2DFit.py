#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 17:16:32 2020

@author: ceballos + Peille


"""

import numpy as np

def fit_2D_gain_scale(ph_values,phase_values,energy_list,deg=6,show=False):
    '''Fits a 2D gain scale giving the energy as a function of Pulse Heigth (reconstructed energy)
    and phase (distance in samples from triggering point to maximu-energy-point)

    Arguments:
        - ph_values: 2D array of the measured PH values
        - phase_values: 2D array of the measured phase values
        - energy_list: 1D array of the simulated energies
        - deg: degree of the polynom to fit (in both X and Y)
        - show: option to plot the results
    '''
    ph_values = np.array(ph_values)
    phase_values = np.array(phase_values)
    energies = np.array([energy_list for _ in range(len(ph_values[0]))]).transpose()
    poly_coeffs = polyfit2d(ph_values, phase_values, energies, [deg,deg])
    return poly_coeffs

def polyfit2d(x, y, f, deg):
    x = np.asarray(x)
    y = np.asarray(y)
    f = np.asarray(f)
    deg = np.asarray(deg)
    vander = np.polynomial.polynomial.polyvander2d(x, y, deg)
    #print(vander.size)
    #print(vander)
    vander = vander.reshape((-1,vander.shape[-1]))
    f = f.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, f,rcond=None)[0]
    return c.reshape(deg+1)

