#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 16:27:00 2020

@author: ceballos
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly


def jitterCorr(reconPH=None, phase=None, deg=2, xsize=20, ysize=6, alpha=0.5):

    """
    Calculate jitter correction: dependance of recons PH vs Phase
    Fit a polynomial (deg 2) and subtract effect
    Return: reconstructed PH with jitter removed

    reconPH: (array) reconstructed PH
    phase1: (array) phase information
    xsize: horizontal size of plotting area
    ysize: vertical size of plotting area
    alpha: transparency in scatter plot

    """

    # plot phases
    fig = plt.figure(figsize=(xsize,ysize))
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.scatter(phase, reconPH, marker="o", alpha=alpha)
    # fit polynomial
    x_interval_for_fit = np.linspace(min(phase), max(phase), 10000)
    ##poly1 = np.poly1d(np.polyfit(phaseKas_HR, enerKas_HR, 2))
    # poly.polyfit recommended instead of np.polyfit + np.poly1d
    coefs = poly.polyfit(x=phase, y=reconPH, deg=deg)
    ffit = poly.polyval(x_interval_for_fit, coefs)
    ax1.plot(x_interval_for_fit, ffit,'-', color="red")
    ax1.set_xlabel("Phase (samples)")
    ax1.set_ylabel("PH (a.u.)")
    # subtract polynomial (flat jitter effect)
    ax2 = fig.add_subplot(1, 2, 2)

    fit_txt = "Fit PH=" + '{:0.3f}'.format(coefs[0])
    reconPH_jitter = np.copy(reconPH)
    for i in range(1,len(coefs)):
        fit_txt +=  "+ (" + '{:0.3f}'.format(coefs[i]) + ")*x**" + str(i)
        reconPH_jitter -= coefs[i]*phase**i
    if min(reconPH_jitter) < min(reconPH):
        reconPH_jitter += min(reconPH)-min(reconPH_jitter)
    print(fit_txt)

    ax2.scatter(phase, reconPH_jitter, marker="o", alpha=alpha)
    ax2.set_xlabel("Phase (samples)")
    ax2.set_ylabel("Corrected PH (a.u.)")
    coefsJ = poly.polyfit(x=phase, y=reconPH_jitter, deg=deg)

    fit_txt = "Fit PH corrected=" + '{:0.3f}'.format(coefsJ[0])
    for i in range(1,len(coefsJ)):
        fit_txt +=  "+ (" + '{:0.3f}'.format(coefsJ[i]) + ")*x**" + str(i)
    print(fit_txt)
    ffitJ = poly.polyval(x_interval_for_fit, coefsJ)
    ax2.plot(x_interval_for_fit, ffitJ,'-', color="white")
    fig.tight_layout()

    return(reconPH_jitter)
