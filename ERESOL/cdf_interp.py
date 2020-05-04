#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 16:20:46 2020

@author: ceballos
"""

import numpy as np 
from scipy.interpolate import interp1d


def cdf_interp(model, minx, maxx, nsteps):
    """ Compute approx. (by cumsum) CDF for a model function over equispaced-data
        between minx and maxx over nsteps

        xdata: (array) x data to evaluate model
        model: (func or astropy model)
        returns:
            model_cdf_int: (func) interpolated & computed CDF

    """
    xdata = np.linspace(minx, maxx, nsteps)
    # calculate CDF of model fit
    bsize = xdata[1] - xdata[0]
    model_xdata = model(xdata)
    model05_xdata = model(xdata-0.5*bsize)
    # remove triangles in excess of the approx. CDF
    model_cdf  = np.cumsum(model_xdata)*bsize - model_xdata *bsize/2. - (model_xdata-model05_xdata)*bsize/8.
    model_cdf_int = interp1d(xdata, model_cdf)
    return model_cdf_int
