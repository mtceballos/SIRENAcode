#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 13:22:39 2020

@author: ceballos
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models


class RxLines:
    def __init__(self, complabel, ilabels, energies_eV,
                 fwhms_eV, rel_amplitudes):
        self.complabel = complabel
        self.ilabels = ilabels
        self.energies_eV = energies_eV
        self.fwhms_eV = fwhms_eV
        self.rel_amplitudes = rel_amplitudes

    def getNumber(self):
        '''Get number of lines'''
        return len(self.ilabels)

    def plotLorentz(self):
        '''Ploting Lorentzian profiles for every line'''
        fig = plt.figure(figsize=(9, 6))
        ax1 = fig.add_subplot(1, 1, 1)
        minx = min(self.energies_eV)-20
        maxx = max(self.energies_eV)+20
        x_interval = np.linspace(minx, maxx, 1000)
        LmodSum = models.Lorentz1D(x_0=0, amplitude=0., fwhm=0.)
        nlines = self.getNumber()
        for i in range(nlines):
            Lmod = models.Lorentz1D(x_0=self.energies_eV[i], amplitude=self.rel_amplitudes[i],
                                    fwhm=self.fwhms_eV[i])
            LmodSum += Lmod
            ax1.plot(x_interval, Lmod(x_interval), label=self.ilabels[i])

        title = self.complabel + ' line complex'
        ax1.plot(x_interval, LmodSum(x_interval), marker='.', label=(self.complabel + ' complex'))
        ax1.set_xlabel("Energy (eV)")
        ax1.set_ylabel("Relative Intensity")
        ax1.set_title(title)
        ax1.legend()

    def broadGauss(self, sigma):
        '''Broad Lorentz lines profile with a Gaussian of given std_dev and plot results

        Parameters:
        sigma: standard deviation of the Gaussian to broaden the line profile

        '''
        fwhm = sigma*2*np.sqrt(2*np.log(2))
        fig = plt.figure(figsize=(9, 6))
        ax1 = fig.add_subplot(1, 1, 1)
        minx = min(self.energies_eV)-20
        maxx = max(self.energies_eV)+20
        x_interval = np.linspace(minx, maxx, 1000)
        LmodSum = models.Voigt1D(x_0=0, amplitude_L=0., fwhm_L=1., fwhm_G=1.)
        nlines = self.getNumber()
        for i in range(nlines):
            Lmod = models.Voigt1D(x_0=self.energies_eV[i], amplitude_L=self.rel_amplitudes[i],
                                  fwhm_L=self.fwhms_eV[i], fwhm_G=fwhm)
            LmodSum += Lmod
            ax1.plot(x_interval, Lmod(x_interval), label=self.ilabels[i])

        title = (self.complabel + ' Gauss-widened (sigma=' + str(sigma) + 'eV) line complex')
        ax1.plot(x_interval, LmodSum(x_interval), marker='.', label=(self.complabel + ' complex'))
        ax1.set_xlabel("Energy (eV)")
        ax1.set_ylabel("Relative Intensity")
        ax1.set_title(title)
        ax1.legend()
