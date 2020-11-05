#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 13:24:42 2020

@author: ceballos
"""

import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import numpy as np

def fit2GaussAndRatio(data=None, a1=50, a2=90, mean1=5800, mean2=5900, sig1=5, sig2=5, nbins1=200, ratio=None,
                      xlab=None, xlim=(0, 0), ylim=(0, 0), xsize=10, ysize=4):

    """"

    Fit 2 Gaussians (Ka1, Ka2) to Kas histogram
    Histograms are created and plotted with matplotlib.pyplot.hist in Density
    Gaussians functions from astropy.fitting module (fitting byLevMarLSQFitter)

    data1: (array) data for 1st histogram
    a1: (float) initial amplitude for Gaussian1
    a2: (float)initial amplitude for Gaussian2
    mean1: (float)initial mean for Gaussian1
    mean2: (float)initial mean for Gaussian2
    sig1: (float)std dev for Gaussian1
    sig2: (float)std dev for Gaussian2
    nbins1: number of bins for first (Kas) histogram
    ratio: GaussKa2/GaussKa1 ratio to select Ka2 photons
    xlab: xlabel of histogram plot
    xlim: (xmin,xmax) limits of X axis
    xlim: (ymin,ymax) limits of Y axis
    xsize: horizontal size of plotting area
    ysize: vertical size of plotting area

    returns:
        PHmin,PHmax: x limiting values for Ka2 complex
    """

    fig = plt.figure(figsize=(xsize, ysize))
    ax1 = fig.add_subplot(1, 2, 1)

    # create histogram
    bin_heights, bin_borders, _ = ax1.hist(data, bins=nbins1, density=True, label="Histogram", alpha=0.5)
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
    x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)

    # fit two gaussians to density-histogram (also "curve_fit" ?)
    gg_init = (models.Gaussian1D(amplitude=a1, mean=mean1, stddev=sig1) +
               models.Gaussian1D(amplitude=a2, mean=mean2, stddev=sig2))
    # fitter = fitting.SLSQPLSQFitter()
    fitter = fitting.LevMarLSQFitter()
    gg_fit = fitter(gg_init, bin_centers, bin_heights, maxiter=300)
    print("Message (Kas)=", fitter.fit_info['message'])

    # C1 = gg_fit.param_sets[0][0]
    # mean1 = gg_fit.param_sets[1][0] #u.a.
    # sigma1 = gg_fit.param_sets[2][0] #u.a.
    # C2 = gg_fit.param_sets[3][0]
    # mean2 = gg_fit.param_sets[4][0] #u.a.
    # sigma2 = gg_fit.param_sets[5][0] #u.a.

    C1 = gg_fit.amplitude_0[0]
    mean1 = gg_fit.mean_0[0]
    sigma1 = gg_fit.stddev_0[0]
    # fwhm1 = sigma1 * 2 * np.sqrt(2*np.log(2))
    C2 = gg_fit.amplitude_1[0]
    mean2 = gg_fit.mean_1[0]
    sigma2 = gg_fit.stddev_1[0]
    # fwhm2 = sigma2 * 2 * np.sqrt(2*np.log(2))

    # gg_fit.fit_info['residuals']
    g1 = models.Gaussian1D(amplitude=C1, mean=mean1, stddev=sigma1)
    g2 = models.Gaussian1D(amplitude=C2, mean=mean2, stddev=sigma2)
    ratioGG = g2(x_interval_for_fit)/g1(x_interval_for_fit)
    # print("minGG=",min(ratioGG), "maxGG=",max(ratioGG))

    # plot histogram and Gaussians fit
    ax1.plot(x_interval_for_fit, gg_fit(x_interval_for_fit), label='Gauss fit')
    ax1.plot(x_interval_for_fit, g1(x_interval_for_fit), label="Gauss Ka2")
    ax1.plot(x_interval_for_fit, g2(x_interval_for_fit), label="Gauss Ka1")
    ax1.plot(x_interval_for_fit, ratioGG, label="ratio G1/G2")
    ax1.set_xlabel(xlab)
    ax1.set_ylabel("Density")
    if (xlim[0] > 0 or xlim[1] > 0):
        ax1.set_xlim(xlim)
    if (ylim[0] > 0 or ylim[1] > 0):
        ax1.set_ylim(ylim)
    PHmin = x_interval_for_fit[ratioGG >= ratio][0]
    PHmax = min(x_interval_for_fit[ratioGG >= ratio][-1], (mean2+10*sigma2))
    ax1.axvline(PHmin, linestyle='--', color='tab:purple', label="Region for ratio")
    ax1.axvline(PHmax, linestyle='--', color='tab:purple')
    plt.legend()
    plt.show()
    fig.tight_layout()
    print("Ka1 PHs in [" + '{:0.3f}'.format(PHmin) + "," + '{:0.3f}'.format(PHmax) + "] a.u.")
    return((PHmin, PHmax))
