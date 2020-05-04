#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 16:25:27 2020

@author: ceballos
"""

from astropy.modeling import models, fitting
import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly

def gainScalePolyFit(xData=None, yData=None, deg=2, ylab="Lines energies (eV)"):
    """
    Fit a polynomial model to the input data (lines) and return slope and intercept

    xData: (array) x data
    yData: (array) y data
    deg: (integer) degree of polynomial
    return: coefficients of polynomial fit
    """

    coefs = poly.polyfit(x=xData, y=yData, deg=deg)
    x_interval_for_fit = np.linspace(xData[0], xData[-1], 10000)
    ffit = poly.polyval(x_interval_for_fit, coefs)

    # check R²
    absError = poly.polyval(xData, coefs) - yData
    SE = np.square(absError)  # squared errors
    MSE = np.mean(SE)  # mean squared errors
    RMSE = np.sqrt(MSE)  # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (np.var(absError) / np.var(yData))
    print('RMSE:', RMSE)
    print('R-squared:', Rsquared)

    # plot the model
    fig = plt.figure(figsize=(16, 6))
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.plot(xData, yData, 'ko', label='Data')
    ax1.plot(x_interval_for_fit, ffit, 'k-', label='Fitted Model', color="red")
    ax1.set_xlabel('reconstructed lines (a.u.)')
    ax1.set_ylabel(ylab)
    ax1.set_xlim(xData[0]-5, xData[1]+5)
    ax1.set_ylim(yData[0]-5, yData[1]+5)
    ax1.legend()
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.plot(xData, yData, 'ko', label='Data')
    ax2.plot(x_interval_for_fit, ffit, 'k-', label='Fitted Model', color="red")
    ax2.set_xlabel('reconstructed lines (a.u.)')
    ax2.set_ylabel(ylab)
    ax2.legend()

    return (coefs)


def gainScaleLinearFit(xData=None, yData=None, ylab="Lines energies (eV)"):
    """
    Fit a linear model to the input data (lines) and return slope and intercept

    xData: (array) x data
    yData: (array) y data
    return: slope, intercept
    """

    # define a model for a line; initialize a linear model
    line_init = models.Linear1D()
    # initialize a linear fitter
    fitter = fitting.LinearLSQFitter()
    # fit the data with the fitter
    fitted_line = fitter(line_init, xData, yData)
    print(fitted_line)
    print("Residuals=", fitter.fit_info['residuals'])
    print("Params=", fitter.fit_info['params'])
    slope = fitted_line.slope[0]
    inter = fitted_line.intercept[0]
    # check R²
    absError = fitted_line(xData) - yData
    SE = np.square(absError)  # squared errors
    MSE = np.mean(SE)  # mean squared errors
    RMSE = np.sqrt(MSE)  # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (np.var(absError) / np.var(yData))
    print('RMSE:', RMSE)
    print('R-squared:', Rsquared)

    # plot the model
    fig = plt.figure(figsize=(16, 6))
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.plot(xData, yData, 'ko', label='Data')
    ax1.plot(xData, fitted_line(xData), 'k-', label='Fitted Model', color="red")
    ax1.set_xlabel('reconstructed lines (a.u.)')
    ax1.set_ylabel(ylab)
    ax1.set_xlim(xData[0]-5, xData[1]+5)
    ax1.set_ylim(yData[0]-5, yData[1]+5)
    ax1.legend()
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.plot(xData, yData, 'ko', label='Data')
    ax2.plot(xData, fitted_line(xData), 'k-', label='Fitted Model', color="red")
    ax2.set_xlabel('reconstructed lines (a.u.)')
    ax2.set_ylabel(ylab)
    ax2.legend()

    return (slope, inter)
