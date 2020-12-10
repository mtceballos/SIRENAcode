#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 16:27:00 2020

@author: ceballos
"""

import numpy as np
import numpy.polynomial.polynomial as poly


class applyCorr:
    """
    Calculate correction: dependance of ydata vs xdata (i.e. jitter correction, drift correction, etc.)
    Fit a polynomial (deg 2) and subtract effect
    Return: ydata with dependance removed

    ydata: (1D array) Y data
    xdata: (1D array) X data
    ninter: (integer) fitting is perfomed dividing the data en 'ninter' equidistant intervals
    fit: (string) type of fitting:
        "poly": polynomial fit (numpy.polynomial.polynomial)
        "cubic_spline": natural cubic smoothing spline
    plot: (bool) plotting?
    ax0: (axis) first plot axis for plotting area
    ax1: (axis) second plot axis for plotting area
    alpha: (float) transparency in scatter plot
    size: (float) size of plotting area
    verbose:(integer) verbosity level (>0 => chatty)
    """

    def __init__(self, xdata, ydata, deg, verbose):
        self.xdata = xdata
        self.ydata = ydata
        self.deg = 2
        self.verbose = 0

        #
        # calculate coefficients of polynomial fit
        #
        self.coefs = poly.polyfit(x=self.xdata, y=self.ydata, deg=self.deg)

        #
        # flatten data with polynomial correction
        #
        self.ydata_corr = np.copy(self.ydata)
        # subtract polynomial (flat xdataline effect)
        self.ydata_corr -= poly.polyval(self.xdata, self.coefs)
        #if min(ydata_inter_corr) < min(ydata_inter):
        #    ydata_inter_corr += min(ydata_inter)-min(ydata_inter_corr)
        self.ydata_corr += min(self.ydata)-min(self.ydata_corr) # rescale to same values
        if self.verbose:
            fit_txt = "Fit Y=" + '{:0.3e}'.format(self.coefs[0])
            print(fit_txt)

        #
        # get fit curve to plot
        #
        xmin = min(self.xdata)
        xmax = max(self.xdata)
        x_interval_for_fit = np.linspace(xmin, xmax, 1000)
        self.ffit = poly.polyval(x_interval_for_fit, self.coefs)

        #
        # Fit polynomial to flattened data
        #
        self.coefsCorr = poly.polyfit(x=self.xdata, y=self.ydata_corr, deg=self.deg)
        fit_txt = "Fit Y corrected=" + '{:0.3e}'.format(self.coefsCorr[0])
        for ii in range(1,len(self.coefsCorr)):
            fit_txt +=  "+ (" + '{:0.3e}'.format(self.coefsCorr[ii]) + ")*x**" + str(ii)
        if verbose:
            print(fit_txt)
        self.ffitCorr = poly.polyval(x_interval_for_fit, self.coefsCorr)


    def plotDataFit(self, ax0, alpha=0.5, size=1):
        """
        Plot initital data + polynomial fit

        Parameters
        ----------
        ax0 : TYPE
            DESCRIPTION.
        alpha : TYPE, optional
            DESCRIPTION. The default is 0.5.
        size : TYPE, optional
            DESCRIPTION. The default is 1.

        Returns
        -------
        None.

        """
        xmin = min(self.xdata)
        xmax = max(self.xdata)
        x_interval_for_fit = np.linspace(xmin, xmax, 1000)
        ax0.scatter(self.xdata, self.ydata, marker=".", alpha=alpha, s=size)
        ax0.plot(x_interval_for_fit, self.ffit,'-', color="red")

    def plotDataCorrFit(self, ax0, alpha=0.5, size=1):
        """
        Plot flattened data + polynomial fit

        Parameters
        ----------
        ax0 : TYPE
            DESCRIPTION.
        alpha : TYPE, optional
            DESCRIPTION. The default is 0.5.
        size : TYPE, optional
            DESCRIPTION. The default is 1.

        Returns
        -------
        None.

        """
        xmin = min(self.xdata)
        xmax = max(self.xdata)
        x_interval_for_fit = np.linspace(xmin, xmax, 1000)
        ax0.scatter(self.xdata, self.ydata_corr, marker=".", alpha=alpha, s=size)
        ax0.plot(x_interval_for_fit, self.ffitCorr,'-', color="red")

def applyCorr2(xdata=None, ydata=None, deg=2,
              plot=True, ax0=None, ax1=None, alpha=0.5, size=1, verbose=0):

    """
    Calculate correction: dependance of ydata vs xdata (i.e. jitter correction, drift correction, etc.)
    Fit a polynomial (deg 2) and subtract effect
    Return: ydata with dependance removed

    ydata: (1D array) Y data
    xdata: (1D array) X data
    ninter: (integer) fitting is perfomed dividing the data en 'ninter' equidistant intervals
    fit: (string) type of fitting:
        "poly": polynomial fit (numpy.polynomial.polynomial)
        "cubic_spline": natural cubic smoothing spline
    plot: (bool) plotting?
    ax0: (axis) first plot axis for plotting area
    ax1: (axis) second plot axis for plotting area
    alpha: (float) transparency in scatter plot
    size: (float) size of plotting area
    verbose:(integer) verbosity level (>0 => chatty)

    """

    interval_size = max(xdata) - min(xdata)


    # fit polynomial to interval
    xmin = min(xdata) + interval_size
    xmax = xmin + interval_size
    xdata_inter = xdata[(xdata>=xmin) & (xdata<=xmax)]
    ydata_inter = ydata[(xdata>=xmin) & (xdata<=xmax)]
    x_interval_for_fit = np.linspace(xmin, xmax, 1000)
    ydata_inter_corr = np.copy(ydata_inter)

    coefs = poly.polyfit(x=xdata_inter, y=ydata_inter, deg=deg)
    ffit = poly.polyval(x_interval_for_fit, coefs)
    # subtract polynomial (flat xdataline effect)
    fit_txt = "Fit Y=" + '{:0.3e}'.format(coefs[0])
    ydata_inter_corr -= poly.polyval(xdata_inter, coefs)


    #if min(ydata_inter_corr) < min(ydata_inter):
    #    ydata_inter_corr += min(ydata_inter)-min(ydata_inter_corr)
    ydata_inter_corr += min(ydata_inter)-min(ydata_inter_corr) # rescale to same values
    if verbose:
        print(fit_txt)

    coefsCorr = poly.polyfit(x=xdata_inter, y=ydata_inter_corr, deg=deg)
    fit_txt = "Fit Y corrected=" + '{:0.3e}'.format(coefsCorr[0])
    for ii in range(1,len(coefsCorr)):
        fit_txt +=  "+ (" + '{:0.3e}'.format(coefsCorr[ii]) + ")*x**" + str(ii)
    if verbose:
        print(fit_txt)
    ffitCorr = poly.polyval(x_interval_for_fit, coefsCorr)

    if plot:
        # plot
        ax0.scatter(xdata_inter, ydata_inter, marker=".", alpha=alpha, s=size)
        ax0.plot(x_interval_for_fit, ffit,'-', color="red")
        ax1.scatter(xdata_inter, ydata_inter_corr, marker=".", alpha=alpha, s=size, color="green")
        ax1.plot(x_interval_for_fit, ffitCorr,'-', color="orange")


    ydata_corr = np.copy(ydata_inter_corr)

    return(ydata_corr)
