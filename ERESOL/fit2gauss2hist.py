#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 13:26:24 2020

@author: ceballos
"""

import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import numpy as np

def fit2gauss2hist(data=None, a1=0.05, a2=0.09, mean1=5800, mean2=5900,
                   sig1=5, sig2=5, nbins=200, xlab=None, plot=True, xsize=20, ysize=6):
    """"

    Fit 2 Gaussians (Ka1, Ka2) to Kas histogram
    Histogram is created and plotted with matplotlib.pyplot.hist in Density
    Gaussians functions from astropy.fitting module (fitting byLevMarLSQFitter)

    data: (array) data for histogram
    a1: (float) initial amplitude for Gaussian1
    a2: (float)initial amplitude for Gaussian2
    mean1: (float)initial mean for Gaussian1
    mean2: (float)initial mean for Gaussian2
    sig1: (float)std dev for Gaussian1
    sig2: (float)std dev for Gaussian2
    nbins1: number of bins for first (Kas) histogram
    xlab: xlabel of histogram plot
    plot: (bool) should histogram and fit be plotted?
    xsize: horizontal size of plotting area
    ysize: vertical size of plotting area

    returns:
        (mean1, mean2): tuple of gaussians centres
    """
    fig = plt.figure(figsize=(xsize, ysize))
    ax1 = fig.add_subplot(1, 1, 1)

    # create histogram
    bin_heights, bin_borders, _ = ax1.hist(data, bins=nbins, density=True, label="Histogram", alpha=0.5)
    if not plot:
        plt.clf()

    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
    x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 10000)

    # fit two gaussians to density-histogram
    gg_init = (models.Gaussian1D(amplitude=a1, mean=mean1, stddev=sig1) +
               models.Gaussian1D(amplitude=a2, mean=mean2, stddev=sig2))
    # fitter = fitting.SLSQPLSQFitter()
    # fitter = fitting.LinearLSQFitter()
    #     -> model is not linear in parameters: linear fit should not be used
    fitter = fitting.LevMarLSQFitter()
    gg_fit = fitter(gg_init, bin_centers, bin_heights, maxiter=300)
    # check there are no errors in fitting
    if not fitter.fit_info['ierr'] in [1,2,3,4]:
        print("Solution not found for fit")
        print("Message (Kas) =", fitter.fit_info['message'])

    C1 = gg_fit.amplitude_0[0]
    mean1 = gg_fit.mean_0[0]
    sigma1 = gg_fit.stddev_0[0]
    fwhm1 = sigma1 * 2 * np.sqrt(2*np.log(2))
    C2 = gg_fit.amplitude_1[0]
    mean2 = gg_fit.mean_1[0]
    sigma2 = gg_fit.stddev_1[0]
    fwhm2 = sigma2 * 2 * np.sqrt(2*np.log(2))
    g1 = models.Gaussian1D(amplitude=C1, mean=mean1, stddev=sigma1)
    g2 = models.Gaussian1D(amplitude=C2, mean=mean2, stddev=sigma2)

    if plot:
        # plot 1st histogram and 2 Gaussians fit
        ax1.plot(x_interval_for_fit, gg_fit(x_interval_for_fit), label='Gauss fit')
        ax1.plot(x_interval_for_fit, g1(x_interval_for_fit), label="Gauss Ka2")
        ax1.plot(x_interval_for_fit, g2(x_interval_for_fit), label="Gauss Ka1")
        maxy = max(bin_heights)
        xtxt = min(data)
        ax1.text(xtxt, maxy-0.03, "Double Gaussian", color='tab:orange')
        ax1.text(xtxt, maxy-0.04, ("Mean(Ka2)=" + '{:0.3f}'.format(mean1) + "a.u"), color='tab:green')
        ax1.text(xtxt, maxy-0.05, ("FWHM(Ka2)=" + '{:0.3f}'.format(fwhm1) + "ma.u"), color='tab:green')
        ax1.text(xtxt, maxy-0.06, ("Mean(Ka1)=" + '{:0.3f}'.format(mean2) + "a.u"), color='tab:red')
        ax1.text(xtxt, maxy-0.07, ("FWHM(Ka1)=" + '{:0.3f}'.format(fwhm2) + "ma.u"), color='tab:red')
        npulses = len(data)
        txt ="Total number of pulses:" + str(npulses)
        ax1.text(xtxt, maxy-0.02, txt, color='black')
        ax1.set_xlabel(xlab)
        ax1.set_ylabel("Density")
        ax1.legend()
        fig.tight_layout()

    return((mean1, mean2))
