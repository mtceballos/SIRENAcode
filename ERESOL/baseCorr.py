#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 16:27:00 2020

@author: ceballos
"""

import numpy as np
#import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly


#def baseCorr(reconPH=None, base=None, deg=2, xsize=20, ysize=6, alpha=0.5):
def baseCorr(reconPH=None, base=None, deg=2, ax0=None, ax1=None, alpha=0.5):

    """
    Calculate baseline correction: dependance of recons PH vs baseline
    Fit a polynomial (deg 2) and subtract effect
    Return: reconstructed PH with baseline dependance removed

    reconPH: (array) reconstructed PH
    base: (array) phase information
    xsize: horizontal size of plotting area
    ysize: vertical size of plotting area
    alpha: transparency in scatter plot

    """

    # plot phases
    #fig = plt.figure(figsize=(xsize,ysize))
    #ax1 = fig.add_subplot(1, 2, 1)
    ax0.scatter(base, reconPH, marker="o", alpha=alpha)
    # fit polynomial
    x_interval_for_fit = np.linspace(min(base), max(base), 10000)

    coefs = poly.polyfit(x=base, y=reconPH, deg=deg)
    ffit = poly.polyval(x_interval_for_fit, coefs)
    ax0.plot(x_interval_for_fit, ffit,'-', color="red")
    ax0.set_xlabel("Renormalized Baseline (ADC a.u.)")
    ax0.set_ylabel("PH of events (a.u.)")

    # subtract polynomial (flat baseline effect)
    #ax2 = fig.add_subplot(1, 2, 2)
    fit_txt = "Fit PH=" + '{:0.3f}'.format(coefs[0])
    reconPH_base = np.copy(reconPH)
    for i in range(1,len(coefs)):
        fit_txt +=  "+ (" + '{:0.3f}'.format(coefs[i]) + ")*x**" + str(i)
        reconPH_base -= coefs[i]*base**i
    if min(reconPH_base) < min(reconPH):
        reconPH_base += min(reconPH)-min(reconPH_base)
    print(fit_txt)

    ax1.scatter(base, reconPH_base, marker="o", alpha=alpha)
    ax1.set_xlabel("Renormalized Baseline (ADC a.u.)")
    ax1.set_ylabel("Corrected PH of events (a.u.)")
    coefsJ = poly.polyfit(x=base, y=reconPH_base, deg=deg)
    fit_txt = "Fit PH corrected=" + '{:0.3f}'.format(coefsJ[0])
    for i in range(1,len(coefsJ)):
        fit_txt +=  "+ (" + '{:0.3f}'.format(coefsJ[i]) + ")*x**" + str(i)
    print(fit_txt)
    ffitJ = poly.polyval(x_interval_for_fit, coefsJ)
    ax1.plot(x_interval_for_fit, ffitJ,'-', color="white")
    #fig.tight_layout()


    return(reconPH_base)
