#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 16:28:13 2020

@author: ceballos
"""

import seaborn as sns
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
import numpy as np

def getMaximaDensity(data=None, nmaxima=0, intervalmin=0., intervalmax=0.):
    """ Get x coordinate(s) for maxima of density plot in histogram
        Kernel Density Curve obtained by seaborn.distplot

    data: (array) data in histogram
    nmaxima: number of maxima to locate
    intervalmin.intervalmax:
        (float) limits of range in x where maxima must be enclosed
    return
        centres: (numpy array) with maxima location
    """
    # For a complete plot:
    # fig = plt.figure(figsize=(20, 6))
    # ax1 = fig.add_subplot(1, 2, 1)
    # sns.distplot(1e3*dataKas_HR.SIGNAL, hist=True, kde=True,
    #              bins=nbinsKas, color = 'darkblue',
    #              hist_kws={'edgecolor':'black'},
    #              kde_kws={'linewidth': 3, "label": "KDE"})
    x,y = sns.distplot(data).get_lines()[0].get_data()
    plt.clf() # clear all Axes
    xlines = x[(x>intervalmin) & (x<intervalmax)]
    ylines = y[(x>intervalmin) & (x<intervalmax)]
    ymaxima = argrelextrema(ylines, np.greater)
    centres = xlines[ymaxima]
    # [ax1.axvline(c, color="gray", linestyle='--') for c in centres]
    if len(centres) > nmaxima:
        print("Too many maxima locations:", centres)
        raise RuntimeError
    return centres
