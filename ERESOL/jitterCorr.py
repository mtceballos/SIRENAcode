#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 16:27:00 2020

@author: ceballos
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly


def jitterCorr(reconPH=None, phase=None):

    """
    Calculate jitter correction: dependance of recons PH vs Phase
    Fit a polynomial (deg 2) and subtract effect
    Return: reconstructed PH with jitter removed

    reconPH: (array) reconstructed PH
    phase1: (array) phase information

    """

    # plot phases
    fig = plt.figure(figsize=(20,6))
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.scatter(phase, reconPH, marker="o")
    # fit polynomial
    x_interval_for_fit = np.linspace(min(phase), max(phase), 10000)
    ##poly1 = np.poly1d(np.polyfit(phaseKas_HR, enerKas_HR, 2))
    # poly.polyfit recommended instead of np.polyfit + np.poly1d
    coefs = poly.polyfit(x=phase, y=reconPH, deg=2)
    ffit = poly.polyval(x_interval_for_fit, coefs)
    ax1.plot(x_interval_for_fit, ffit,'-', color="red")
    ax1.set_xlabel("Phase (samples)")
    ax1.set_ylabel("PH Kas (a.u.)")
    print("Fit Kas=",'{:0.3f}'.format(coefs[0]) + "+ (" + '{:0.3f}'.format(coefs[1]) + ")*x" +
          "+(" + '{:0.3f}'.format(coefs[2]) + ")*x²" )
    # subtract polynomial (flat jitter effect)
    ax2 = fig.add_subplot(1, 2, 2)
    reconPH_jitter = reconPH - coefs[1]*phase - coefs[2]*phase**2
    ax2.scatter(phase, reconPH_jitter, marker="o")
    ax2.set_xlabel("Phase (samples)")
    ax2.set_ylabel("Corrected PH Kas (a.u.)")
    coefsJ = poly.polyfit(x=phase, y=reconPH_jitter, deg=2)
    print("Fit Kas (corrected)=",'{:0.3f}'.format(coefsJ[0]) + "+ (" + '{:0.3f}'.format(coefsJ[1]) +
          ")*x" + "+(" + '{:0.3f}'.format(coefsJ[2]) + ")*x²" )
    ffitJ = poly.polyval(x_interval_for_fit, coefsJ)
    ax2.plot(x_interval_for_fit, ffitJ,'-', color="white")
    return(reconPH_jitter)
